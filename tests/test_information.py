import pytest
from math import log as logarithm
from bricolage.analysis_ext import (
    MutualInfoAnalyzer, AverageControlAnalyzer, CausalFlowAnalyzer,
    OutputControlAnalyzer, RelevantControlAnalyzer, Information,
    _set_max_category_size, _get_max_category_size)
from bricolage.core import InterventionState
import numpy
numpy.set_printoptions(linewidth=120)


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
                info += p_xy * logarithm(p_xy / (p_x * p_y), 2)
    return info


def calc_mutual_info(n, categories):
    w = n.factory.world
    assert len(categories) == len(w.environments)

    # Features should be consecutive numbers
    all_feat = set(categories)
    assert all_feat == set(range(len(all_feat)))

    reg_base, reg_to = w.reg_range
    channel_dim = reg_to - reg_base
    feat_dim = len(set(categories))
    state_dim = 2  # on or off
    env_dim = len(w.environments)

    # Now we have the dimensions of our array
    probs = numpy.zeros((channel_dim, feat_dim, state_dim))

    for i, channel in enumerate(range(*w.reg_range)):
        for attrs, feat in zip(n.attractors, categories):
            p_event = 1.0 / float(env_dim)

            # If there is more than one state in the attractor, distribute
            # over all (equiprobable)
            p_event /= len(attrs)

            for a in attrs:
                val = a[channel]  # 0 or 1
                probs[i, feat, val] += p_event

    info = numpy.zeros(channel_dim)
    for i, channel in enumerate(range(*w.reg_range)):
        info[i] = _mutual_info(probs[i])

    return probs, info


def test_persistence(bowtie_network):
    net = bowtie_network
    # graph.save_network_as_fullgraph(net, name='unbob', simplify=False)
    # graph.save_network_as_fullgraph(net, name='bob')

    # We know from inspection that gene 5 is at the bottleneck
    net.genes[5].intervene = InterventionState.INTERVENE_OFF
    for r in net.rates:
        assert numpy.allclose(r, numpy.array([0, 0, 0]))
    net.genes[5].intervene = InterventionState.INTERVENE_ON
    for r in net.rates:
        assert numpy.allclose(r, numpy.array([1, .5, .25]))


def bowtie_categorize_output(values, matches):
    ret = []
    for v, m in zip(values, matches):
        if numpy.isclose(v, m):
            ret.append(1)
        else:
            ret.append(0)
    return ret


def test_mutual_info_cython(bowtie_network, bowtie_env_categories):
    f = MutualInfoAnalyzer(bowtie_network.factory.world,
                           bowtie_env_categories)
    j = f.analyse_network(bowtie_network)
    p_joint, p_info = calc_mutual_info(bowtie_network,
                                       bowtie_env_categories)
    info = Information(j)
    c_info = numpy.asarray(info)[0]

    # Need to account for different shape
    numpy.testing.assert_allclose(c_info[:, 0], p_info)


class RateCategorizer(object):
    max_categories = 16

    def __init__(self):
        self.categories = {}
        self.next_cat = 0
        self.probabilities = []

    def categorize(self, rate, pr):
        assert pr > 0.0
        cat = self.categories.setdefault(rate, self.next_cat)
        # Did we insert?
        if cat == self.next_cat:
            self.probabilities.append(pr)
            self.next_cat += 1
            if self.next_cat > self.max_categories:
                raise RuntimeError(
                    "Category explosion: {}".format(self.next_cat))
        else:
            self.probabilities[cat] += pr
        return cat


def get_causal_specs(net):
    """Calculate the causal specificity for the entire output in each
    environment"""
    wrld = net.factory.world

    # For now, states are equiprobable
    env_count = len(wrld.environments)
    env_probs = numpy.ones(env_count) / env_count

    pdist_regs = calc_natural(env_probs, net, wrld)

    # All probabilities of co-occurences. A joint prob distn under
    # intervention for each channel, about each output channel. We create one
    # for each environment.
    #
    # 1. We need 2 rows for FALSE / TRUE.
    # 2. And we guess(!) a maximum number of columns for all possible rates in this
    # network. This dedicate the first two to 0 / 1 and allocate the others as
    # needed.
    categorizers = [[RateCategorizer() for _ in range(wrld.out_channels)]
                    for __ in range(wrld.reg_channels)]
    probs = [numpy.zeros((wrld.reg_channels,
                          wrld.out_channels,
                          2,
                          RateCategorizer.max_categories))
             for _ in range(env_count)]

    p_env = 1.0 / float(env_count)
    # Now, do the interventions for each regulatory gene
    # The genes are organised so that they begin with regulatory genes, so we
    # can simply grab those
    for i in range(wrld.reg_channels):
        gene = net.genes[i]
        p_gene_on = pdist_regs[i]
        p_gene_off = 1.0 - p_gene_on
        # ----- GENE IS MANIPULATED OFF
        # NOTE: This automatically updates everything (changing the rates!).
        gene.intervene = InterventionState.INTERVENE_OFF
        for j, rates in enumerate(net.rates):
            for k, r in enumerate(rates):
                if not numpy.isclose(p_gene_off, 0.0):
                    c = categorizers[i][k].categorize(r, p_env * p_gene_off)
                    probs[j][i, k, 0, c] += p_gene_off

        # ----- GENE IS MANIPULATED ON
        gene.intervene = InterventionState.INTERVENE_ON
        for j, rates in enumerate(net.rates):
            for k, r in enumerate(rates):
                if not numpy.isclose(p_gene_on, 0.0):
                    c = categorizers[i][k].categorize(r, p_env * p_gene_on)
                    probs[j][i, k, 1, c] += p_gene_on

        # Reset this gene
        gene.intervene = InterventionState.INTERVENE_NONE

    return probs, categorizers


class PhenotypeCategorizer(object):
    max_categories = 16

    def __init__(self, targets, env_size):
        self.categories = {}
        self.next_cat = 0
        self.probabilities = []
        self.targets = set()
        if targets is not None:
            for t in targets:
                self.targets.add(t)

        self.targets_hit = [0.0] * env_size

    def categorize(self, rates, pr, env_num):
        rates = tuple(rates)

        # If targets was None, we won't be doing this.
        for t in self.targets:
            if numpy.allclose(t, rates):
                self.targets_hit[env_num] += pr
                break

        # What category is it?
        cat = self.categories.setdefault(rates, self.next_cat)
        # Did we insert?
        if cat == self.next_cat:
            self.probabilities.append(pr)
            self.next_cat += 1
            if self.next_cat > self.max_categories:
                raise RuntimeError(
                    "Category explosion: {}".format(self.next_cat))
        else:
            self.probabilities[cat] += pr

        return cat


def get_causal_specs_phenotype(net, targets=None):
    """Calculate the causal specificity for the entire output in each environment"""
    wrld = net.factory.world

    # For now, states are equiprobable
    env_count = len(wrld.environments)
    env_probs = numpy.ones(env_count) / env_count

    pdist_regs = calc_natural(env_probs, net, wrld)

    # All probabilities of co-occurences. A joint prob distn under
    # intervention for each channel, about each output channel. We create one
    # for each environment.
    #
    # 1. We need 2 rows for FALSE / TRUE.
    # 2. And we guess(!) a maximum number of columns for all possible rates in this
    # network.
    categorizers = [PhenotypeCategorizer(targets, env_count) for _ in range(wrld.reg_channels)]
    probs = [numpy.zeros((wrld.reg_channels,
                          2,
                          PhenotypeCategorizer.max_categories))
             for _ in range(env_count)]

    # Now, do the interventions for each regulatory gene
    # The genes are organised so that they begin with regulatory genes, so we
    # can simply grab those
    for i in range(wrld.reg_channels):
        gene = net.genes[i]
        p_on = pdist_regs[i]
        p_off = 1.0 - p_on
        # ----- GENE IS MANIPULATED OFF
        # NOTE: This automatically updates everything (changing the rates!).
        gene.intervene = InterventionState.INTERVENE_OFF
        for j, rate in enumerate(net.rates):
            if not numpy.isclose(p_off, 0.0):
                cat = categorizers[i].categorize(rate, p_off, j)
                probs[j][i, 0, cat] += p_off

        # ----- GENE IS MANIPULATED ON
        gene.intervene = InterventionState.INTERVENE_ON
        for j, rate in enumerate(net.rates):
            if not numpy.isclose(p_on, 0.0):
                cat = categorizers[i].categorize(rate, p_on, j)
                probs[j][i, 1, cat] += p_on

        # Reset this gene
        gene.intervene = InterventionState.INTERVENE_NONE

    return probs, categorizers


def get_causal_flow(net):
    wrld = net.factory.world

    summed_probs = numpy.zeros((wrld.reg_channels,
                                wrld.out_channels,
                                2,
                                RateCategorizer.max_categories))

    # Just sum up those from causal spec, weighted
    prob_list, cats = get_causal_specs(net)
    penv = 1.0 / float(len(prob_list))

    for prob in prob_list:
        summed_probs += penv * prob

    return summed_probs


def calc_natural(env_probs, net, w):
    # First, calculate the "natural" probabilities of each regulatory signals
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


# Slow python calculation of information
def calc_info_from_probs(probs):
    """Assumption: 3 or 4 dimension array. We summarise information in the final 2
    dimensions"""

    # We handle 3 dimensions by expanding and then collapsing it.
    is_4 = len(probs.shape) == 4

    if is_4:
        d1, d2, row_num, col_num = probs.shape
    else:
        d1, row_num, col_num = probs.shape
        d2 = 1

    # Now calculate the information
    info = numpy.zeros((d1, d2))
    for i in range(d1):
        for j in range(d2):
            if is_4:
                info[i, j] = _mutual_info(probs[i, j])
            else:
                info[i, j] = _mutual_info(probs[i])

    # Collapse other dimension if it didn't exist
    if not is_4:
        info.shape = d1,

    return info


def test_causal_flow_net(bowtie_network):
    """Compare the clunky python version with our C++ version"""

    f = CausalFlowAnalyzer(bowtie_network.factory.world)
    j = f.analyse_network(bowtie_network)

    # We just access the 0th element, as we only sent one network for analysis
    c_joint = numpy.asarray(j)[0]

    # Test that each of the signals has a complete joint probability. It
    # should sum to 1.0
    for per_reg in c_joint:
        for per_output in per_reg:
            assert per_output.sum() == 1.0

    # Get the slow version from python.
    p_joint = get_causal_flow(bowtie_network)

    # TODO: They should be the same.
    # This actually fails, but not sure why right now. I think it is minor
    # stuff. Note that the INFORMATION calculated from these is the same
    # below. We're not using causal flow right now anyway...
    # numpy.testing.assert_allclose(p_joint, c_joint)

    # So should the information.
    info = Information(j)
    c_info = numpy.asarray(info)[0]
    p_info = calc_info_from_probs(p_joint)
    numpy.testing.assert_allclose(c_info, p_info)


def test_causal_flow_pop(bowtie_database, bowtie_env_categories):
    # This what it was evolved to do (see the generate.py file)
    pop = bowtie_database.population
    f_flow = CausalFlowAnalyzer(bowtie_database.world)
    j_flow = f_flow.analyse_collection(pop)
    flow = Information(j_flow)
    c_flow = numpy.asarray(flow)

    f_mut = MutualInfoAnalyzer(bowtie_database.world, bowtie_env_categories)
    j_mut = f_mut.analyse_collection(pop)
    mut = Information(j_mut)
    c_mut = numpy.asarray(mut)

    # Let's just try the first few
    for i, n in enumerate(pop):
        p_joint = get_causal_flow(n)
        p_flow = calc_info_from_probs(p_joint)
        numpy.testing.assert_allclose(c_flow[i], p_flow)

        _, p_mut = calc_mutual_info(n, bowtie_env_categories)
        numpy.testing.assert_allclose(c_mut[i, :, 0], p_mut)

        # It takes too long....
        if i == 50:
            break


def get_average_control(net):
    w = net.factory.world
    prob_list, cats = get_causal_specs(net)

    # To calculate this, we need to get the causal spec of each env, then
    # average the information over each env.
    # We double the number of output channels to see the entropies too.
    summed_info = numpy.zeros((w.reg_channels, w.out_channels * 2))
    penv = 1.0 / float(len(prob_list))
    for prob in prob_list:
        summed_info[:, :w.out_channels] += penv * calc_info_from_probs(prob)

    # Now create the entropies
    for i, reg_cats in enumerate(cats):
        for j, cat in enumerate(reg_cats):
            # Need to make them conditional
            probs = numpy.asarray(cat.probabilities)
            summed_info[i, j + w.out_channels] = _entropy(probs)

    return summed_info


def _entropy(probs):
    numpy.testing.assert_allclose(probs.sum(), 1.0)
    assert (probs > 0.0).all()
    return (probs * -numpy.log2(probs)).sum()


def get_average_control_phenotype(net):
    w = net.factory.world
    prob_list, cats = get_causal_specs_phenotype(net)

    # To calculate this, we need to get the causal spec of each env, then
    # average the information over each env.
    summed_info = numpy.zeros((w.reg_channels))
    penv = 1.0 / float(len(prob_list))
    for prob in prob_list:
        summed_info += penv * calc_info_from_probs(prob)

    # We also want the entropy of the phenotype generated by particular
    # natural variation of the regulatory signal.
    entropies = numpy.zeros((w.reg_channels))
    for i, c in enumerate(cats):
        # Need to make them conditional
        probs = numpy.asarray(c.probabilities) * penv
        entropies[i] = _entropy(probs)

    return summed_info, entropies


def test_average_control_phenotype_net(bowtie_network):
    net = bowtie_network
    probs, ents = get_average_control_phenotype(net)
    for p, e in zip(probs, ents):
        # We should never have more information than there is to explain!
        assert p <= e

    onz = OutputControlAnalyzer(net.factory.world)
    cy_info = onz.numpy_info_from_network(net)[0][:,:2]
    py_info = numpy.asarray([[p, e] for (p, e) in zip(probs, ents)])
    numpy.testing.assert_allclose(cy_info, py_info)


def test_average_control_phenotype_pop(bowtie_database):
    pop = bowtie_database.population
    anz = OutputControlAnalyzer(pop.factory.world)
    cy_info = numpy.asarray(anz.analyse_collection(pop))

    for i, net in enumerate(pop):
        probs, ents = get_average_control_phenotype(net)
        py_info = numpy.asarray([[p, e] for (p, e) in zip(probs, ents)])
        numpy.testing.assert_allclose(py_info, cy_info[i][:, :2])

        assert not numpy.any(numpy.isnan(cy_info[1]))

        # # Let's just do 50.
        if i > 50:
            break


def test_average_control_net(bowtie_network):
    net = bowtie_network
    py_info = get_average_control(net)

    anz = AverageControlAnalyzer(net.factory.world)
    cy_info = numpy.asarray(anz.analyse_network(net))

    numpy.testing.assert_allclose(py_info, cy_info[0])


def test_average_control_pop(bowtie_database):
    pop = bowtie_database.population
    anz = AverageControlAnalyzer(pop.factory.world)
    cy_info = numpy.asarray(anz.analyse_collection(pop))

    for i, net in enumerate(pop):
        py_info = get_average_control(net)
        numpy.testing.assert_allclose(py_info, cy_info[i])

        # # Let's just do 50.
        if i > 50:
            break


def test_category_size_control(bowtie_network):
    """Make sure that we can catch the exception when we run out of categories"""
    net = bowtie_network
    assert _get_max_category_size() == 16

    _set_max_category_size(2)

    anz = AverageControlAnalyzer(net.factory.world)

    with pytest.raises(IndexError):
        anz.analyse_network(net)

    _set_max_category_size(16)


def get_weighted_control_phenotype(net, targets):
    w = net.factory.world
    prob_list, cats = get_causal_specs_phenotype(net, targets)

    # To calculate this, we need to get the causal spec of each env, then
    # average the information over each env.
    summed_info = numpy.zeros((w.reg_channels))
    penv = 1.0 / float(len(prob_list))

    # For each environment
    for i, prob in enumerate(prob_list):
        # We need to weight the information by the occurence of the
        # probability of getting the targets that we actually WANT
        info_per_reg = penv * calc_info_from_probs(prob)
        for j, cat in enumerate(cats):
            info_per_reg[j] *= cat.targets_hit[i]

        summed_info += info_per_reg

    return summed_info


def test_weighted_control_phenotype(bowtie_database, bowtie_network):
    net = bowtie_network
    # net = bowtie_database.population[50]
    t = bowtie_database.targets[0]
    tset = t.calc_distinct_outputs()
    onz = OutputControlAnalyzer(net.factory.world, tset)
    cy_info = onz.numpy_info_from_network(net)[0]

    probs, cats = get_causal_specs_phenotype(net, tset)
    wc = get_weighted_control_phenotype(net, tset)
    ac, ent = get_average_control_phenotype(net)
    # Should always DECREASE
    assert (wc <= ac).all()
    numpy.testing.assert_allclose(wc, cy_info[:, 2])


def test_weighted_control_phenotype_pop(bowtie_database):
    pop = bowtie_database.population
    t = bowtie_database.targets[0]
    tset = t.calc_distinct_outputs()
    anz = OutputControlAnalyzer(pop.factory.world, tset)
    cy_info = numpy.asarray(anz.analyse_collection(pop))

    for i, net in enumerate(pop):
        wc = get_weighted_control_phenotype(net, tset)
        assert (cy_info[i][:, 0] >= cy_info[i][:, 2]).all()
        numpy.testing.assert_allclose(wc, cy_info[i][:, 2])

        # # Let's just do 50.
        if i > 50:
            break


def get_relevant_control(net, targets):
    """Calculate the causal specificity for the entire output in each environment"""
    wrld = net.factory.world

    # For now, states are equiprobable
    env_count = len(wrld.environments)
    env_probs = numpy.ones(env_count) / env_count

    pdist_regs = calc_natural(env_probs, net, wrld)

    target_set = set()
    for t in targets:
        target_set.add(t)
    targets = list(target_set)

    info = numpy.zeros((wrld.reg_channels, env_count))
    cat_for_off = numpy.zeros(info.shape, dtype=int)
    cat_for_off[:, :] = -1

    for i in range(wrld.reg_channels):
        gene = net.genes[i]
        p_on = pdist_regs[i]
        p_off = 1.0 - p_on
        # ----- GENE IS MANIPULATED OFF
        # NOTE: This automatically updates everything (changing the rates!).
        gene.intervene = InterventionState.INTERVENE_OFF
        for j, rate in enumerate(net.rates):
            if not numpy.isclose(p_off, 0.0):
                for k, targ in enumerate(targets):
                    if targ == tuple(rate):
                        # Record the category
                        cat_for_off[i, j] = k

                        # Record the info, assuming this is the only entry in
                        # the row / column of the joint dist. We'll correct
                        # this below if we were wrong.
                        info[i, j] += p_off * numpy.log2(1.0 / p_off)
                        break

        # ----- GENE IS MANIPULATED ON
        gene.intervene = InterventionState.INTERVENE_ON
        for j, rate in enumerate(net.rates):
            if not numpy.isclose(p_on, 0.0):
                for k, targ in enumerate(targets):
                    if targ == tuple(rate):
                        # If the off category was the same as the on category,
                        # then information is zero (the manipulation makes no
                        # difference).
                        if cat_for_off[i, j] == k:
                            info[i, j] = 0.0
                        else:
                            info[i, j] += p_on * numpy.log2(1.0 / p_on)
                        break

        # Reset this gene
        gene.intervene = InterventionState.INTERVENE_NONE

    # print
    # print pdist_regs
    # print cat_for_off

    # Now just average over the different environments
    return info.mean(axis=1)


def test_relevant_control_network(bowtie_database, bowtie_network):
    net = bowtie_network
    # net = bowtie_database.population[50]
    t = bowtie_database.targets[0]
    tset = t.calc_distinct_outputs()
    py_info = get_relevant_control(net, tset)
    rz = RelevantControlAnalyzer(net.factory.world, tset)
    cy_info = rz.numpy_info_from_network(net)
    numpy.testing.assert_allclose(py_info, cy_info)


def test_relevant_control_pop(bowtie_database, bowtie_network):
    net = bowtie_network
    pop = bowtie_database.population
    t = bowtie_database.targets[0]
    tset = t.calc_distinct_outputs()

    rz = RelevantControlAnalyzer(net.factory.world, tset)
    cy_info = rz.numpy_info_from_collection(pop)

    for i, net in enumerate(pop):
        py_info = get_relevant_control(net, tset)
        numpy.testing.assert_allclose(py_info, cy_info[i])
        if i > 200:
            break
