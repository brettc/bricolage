from bricolage.analysis_ext import MIAnalyzer
#     MutualInfoAnalyzer, AverageControlAnalyzer, CausalFlowAnalyzer,
#     OutputControlAnalyzer, RelevantControlAnalyzer, Information,
#     _set_max_category_size, _get_max_category_size)
from bricolage.core import InterventionState
import numpy
numpy.set_printoptions(linewidth=120)


def _max_entropy(cats):
    # Assuming 1 and 0's, equiproble
    arr = numpy.asarray(cats)
    tot = float(len(cats))
    num_1s = float((arr == 1).sum())

    p_on = num_1s / tot
    p_off = 1.0 - p_on
    
    probs = numpy.array([p_on, p_off])
    return -(numpy.log2(probs) * probs).sum()


def _mutual_info(joint):
    assert numpy.isclose(joint.sum(), 1.0)

    info = 0.0
    rownum, colnum = joint.shape
    colsum = joint.sum(axis=0)
    rowsum = joint.sum(axis=1)
    for row in range(rownum):
        for col in range(colnum):
            p_xy = joint[row, col]
            p_x = rowsum[row]
            p_y = colsum[col]
            if p_xy != 0:
                info += p_xy * numpy.log2(p_xy / (p_x * p_y))
    return info


def calc_mutual_info(n, categories):
    w = n.factory.world
    assert len(categories) == len(w.environments)

    # Features should be consecutive numbers from 0
    all_feat = set(categories)
    assert all_feat == set(range(len(all_feat)))

    # Okay, we want to generate mutual information for every category
    # separately. We RE-categorize everything into just two categories, our
    # focal category and anything else.
    for focal in range(len(all_feat)):
        bin_cats = [int(x is focal) for x in categories]
        print "entropy", _max_entropy(bin_cats)
        print calc_mutual_info_binary(n, bin_cats)


def calc_mutual_info_binary(n, binary_cats):
    w = n.factory.world

    reg_base, reg_to = w.reg_range
    channel_dim = reg_to - reg_base
    feat_dim = 2
    state_dim = 2  # on or off
    env_dim = len(w.environments)

    # Now we have the dimensions of our array
    probs = numpy.zeros((channel_dim, feat_dim, state_dim))

    for i, channel in enumerate(range(*w.reg_range)):
        for attrs, feat in zip(n.attractors, binary_cats):
            p_event = 1.0 / float(env_dim)

            # If there is more than one state in the attractor, distribute
            # over all (equiprobable)
            p_event /= len(attrs)

            for a in attrs:
                val = a.test(channel)
                probs[i, feat, val] += p_event

    info = numpy.zeros(channel_dim)
    for i, channel in enumerate(range(*w.reg_range)):
        info[i] = _mutual_info(probs[i])

    return info


def calc_natural(env_probs, net, w):
    # Calculate the "natural" probabilities of each regulatory signals
    # without intervention
    pdist_regs = numpy.zeros(w.reg_channels)
    reg_base, reg_to = w.reg_range
    for env_i, attr in enumerate(net.attractors):
        # Reduce probability by number of attractor states
        p_state = env_probs[env_i] / float(len(attr))
        for st in attr:
            # For all attractors states, and each regulatory channel
            for i in range(w.reg_channels):
                # If it is on, then add the probability
                if st.test(reg_base + i):
                    pdist_regs[i] += p_state

    return pdist_regs


class RateCategorizer(object):
    def __init__(self, targets):
        # Uniqize
        self.targets = list(set(targets))
        self.num_cats = len(self.targets)

    def categorize(self, index, current):
        assert 0 <= index < self.num_cats
        # Return is either
        #  1: same as the index
        #  0: One of the other targets
        # -1: NONE of the targets
        for i, t in enumerate(self.targets):
            if current == t:
                if i == index:
                    return 1
                return 0

        # Not found
        return -1


def get_relevant_control(net, targets):
    """Calculate the causal specificity for the entire output in each
    environment"""
    wrld = net.factory.world

    # For now, states are equiprobable
    env_count = len(wrld.environments)
    env_probs = numpy.ones(env_count) / env_count

    pdist_regs = calc_natural(env_probs, net, wrld)

    rate_cat = RateCategorizer(targets)

    info = numpy.zeros((wrld.reg_channels, env_count, rate_cat.num_cats))
    cat_for_off = numpy.zeros(info.shape, dtype=int)
    cat_for_off[:, :, :] = -1

    for i in range(wrld.reg_channels):
        gene = net.genes[i]
        p_on = pdist_regs[i]
        p_off = 1.0 - p_on
        # ----- GENE IS MANIPULATED OFF
        # NOTE: This automatically updates everything (changing the rates!).
        gene.intervene = InterventionState.INTERVENE_OFF
        for j, rate in enumerate(net.rates):
            if not numpy.isclose(p_off, 0.0):
                for k in range(rate_cat.num_cats):
                    cat = rate_cat.categorize(k, tuple(rate))
                    cat_for_off[i, j, k] = cat
                    if cat != -1:
                        # Record the info, assuming this is the only entry in
                        # the row / column of the joint dist. We'll correct
                        # this below if we were wrong.
                        info[i, j, k] += p_off * numpy.log2(1.0 / p_off)

        # ----- GENE IS MANIPULATED ON
        gene.intervene = InterventionState.INTERVENE_ON
        for j, rate in enumerate(net.rates):
            if not numpy.isclose(p_on, 0.0):
                for k in range(rate_cat.num_cats):
                    cat = rate_cat.categorize(k, tuple(rate))
                    # If the off category was the same as the on category,
                    # then information is zero (the manipulation makes no
                    # difference).
                    if cat != -1:
                        if cat_for_off[i, j, k] == cat:
                            info[i, j, k] = 0.0
                        else:
                            info[i, j, k] += p_on * numpy.log2(1.0 / p_on)

        # Reset this gene
        gene.intervene = InterventionState.INTERVENE_NONE

    # Now just average over the different environments
    return info.mean(axis=1)


# def test_1(three_database):
#     from bricolage.dot_layout import save_network_as_fullgraph
#     from bricolage.graph_maker import GraphType
#     p = three_database.population
#     t = three_database.targets[0]
#     cats = t.calc_categories()
#     n = p.get_best()[0]
#     print n.fitness
#     print n.identifier
#     save_network_as_fullgraph(n, graph_type=GraphType.GENE)
#     print cats
#     print calc_mutual_info(n, cats)
#
#     cmi = MIAnalyzer(n.factory.world, cats)
#     print cmi.numpy_info_from_network(n)

    # targs = t.calc_distinct_outputs() 
    # rc = get_relevant_control(n, targs)
    # print rc

# def test_3(three_database):
#     p = three_database.population
#     n = p.get_best()[1]
#     print n.fitness
#     print n.attractor_robustness
#     print n.rates
#     print n.attractors[-1]
#     # n.genes[1].intervene = InterventionState.INTERVENE_ON
#     n.genes[5].intervene = InterventionState.INTERVENE_OFF
#     print n.rates
#     print n.attractors[-1]


# def test_3(three_database):
#     p = three_database.population
#     t = three_database.targets[0]
#     n = p.get_best()[0]
#     cats = t.calc_categories()
#     calc_mutual_info(n, cats)
#
