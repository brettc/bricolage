import pytest
from math import log as logarithm
from bricolage.analysis_ext import (
    MutualInfoAnalyzer, AverageControlAnalyzer, CausalFlowAnalyzer,
    Information, _set_max_category_size, _get_max_category_size)
from bricolage.core import InterventionState
import numpy


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
        self.next_cat = 2

    def categorize(self, rates):
        cats = []
        for r in rates:
            if r == 0.0 or r == 1.0:
                cat = int(r)
            else:
                # r is our rate. What category is it?
                cat = self.categories.setdefault(r, self.next_cat)
                # Did we insert?
                if cat == self.next_cat:
                    self.next_cat += 1
                    if self.next_cat > self.max_categories:
                        raise RuntimeError(
                            "Category explosion: {}".format(self.next_cat))

            cats.append(cat)

        return cats


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
    categorizer = RateCategorizer()
    probs = [numpy.zeros((wrld.reg_channels, 
                          wrld.out_channels, 
                          2,
                          categorizer.max_categories))
             for _ in range(env_count)]

    # Now, do the interventions for each regulatory gene
    # The genes are organised so that they begin with regulatory genes, so we
    # can simply grab those
    for i in range(wrld.reg_channels):
        gene = net.genes[i]
        # ----- GENE IS MANIPULATED OFF
        # NOTE: This automatically updates everything (changing the rates!).
        gene.intervene = InterventionState.INTERVENE_OFF
        for j, rates in enumerate(net.rates):
            cats = categorizer.categorize(rates)
            for k, c in enumerate(cats):
                p_state = (1.0 - pdist_regs[i])
                probs[j][i, k, 0, c] += p_state

        # ----- GENE IS MANIPULATED ON
        gene.intervene = InterventionState.INTERVENE_ON
        for j, rates in enumerate(net.rates):
            cats = categorizer.categorize(rates)
            for k, c in enumerate(cats):
                p_state = pdist_regs[i]
                probs[j][i, k, 1, c] += p_state

        # Reset this gene
        gene.intervene = InterventionState.INTERVENE_NONE

    return probs


class PhenotypeCategorizer(object):
    max_categories = 64

    def __init__(self):
        self.categories = {}
        self.next_cat = 0
        self.probabilities = []

    def categorize(self, rates, pr):
        if float(pr) == 0.0:
            print "FUCK"

        rates = tuple(rates)
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


def get_causal_specs_phenotype(net):
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
    categorizers = [PhenotypeCategorizer() for _ in range(wrld.reg_channels)]
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
        for j, rates in enumerate(net.rates):
            if p_off != 0.0:
                cat = categorizers[i].categorize(rates, p_off)
                probs[j][i, 0, cat] += p_off

        # ----- GENE IS MANIPULATED ON
        gene.intervene = InterventionState.INTERVENE_ON
        for j, rates in enumerate(net.rates):
            if p_on != 0.0:
                cat = categorizers[i].categorize(rates, p_on)
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
    prob_list = get_causal_specs(net)
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

    # They should be the same.
    numpy.testing.assert_allclose(p_joint, c_joint)

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
    prob_list = get_causal_specs(net)

    # To calculate this, we need to get the causal spec of each env, then
    # average the information over each env.
    summed_info = numpy.zeros((w.reg_channels, w.out_channels))
    penv = 1.0 / float(len(prob_list))
    for prob in prob_list:
        summed_info += penv * calc_info_from_probs(prob)
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
        # We should never have more information than there is!!
        assert p <= e
        print e - p
    

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
        


