import pytest
from math import log as logarithm
from bricolage.analysis_ext import (MutualInfoAnalyzer,
                                    CausalFlowAnalyzer, Information)

# from bricolage import threshold3, lineage, graph
from bricolage.core import InterventionState

from generate import get_database
import numpy

@pytest.yield_fixture
def bowtie_database():
    db = get_database('bowtie')
    yield db
    db.close()

@pytest.fixture
def bowtie_env_categories(bowtie_database):
    """Categorise the targets"""
    targ = bowtie_database.targets[0]
    cat_dict = {}
    cats = []
    cat_n = 0
    for et in targ.as_array():
        et = tuple(et)
        if et in cat_dict:
            cats.append(cat_dict[et])
        else:
            cat_dict[et] = cat_n
            cats.append(cat_n)
            cat_n += 1
    return cats

@pytest.fixture
def bowtie_network(bowtie_database):
    return bowtie_database.population.get_best()[0]

def calc_mutual_info(n, categories):
    w = n.constructor.world
    assert len(categories) == len(w.environments)

    # Features should be consecutive numbers
    all_feat = set(categories)
    assert all_feat == set(range(len(all_feat)))

    reg_base, reg_to = w.reg_range
    channel_dim = reg_to - reg_base
    feat_dim = len(set(categories))
    state_dim = 2 # on or off
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
                val = a[channel] # 0 or 1
                probs[i, feat, val] += p_event

    info = numpy.zeros(channel_dim)
    for i, channel in enumerate(range(*w.reg_range)):
        c_prob = probs[i]
        I = 0.0
        colsum = c_prob.sum(axis=0)
        rowsum = c_prob.sum(axis=1)
        for row in range(feat_dim):
            for col in 0, 1:
                p_xy = c_prob[row, col]
                p_x = rowsum[row]
                p_y = colsum[col]
                if p_xy != 0:
                    I += p_xy * logarithm(p_xy / (p_x * p_y), 2)
        info[i] = I

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
    f = MutualInfoAnalyzer(bowtie_network.constructor.world, 
                           bowtie_env_categories)
    j = f.analyse_network(bowtie_network)
    p_joint, p_info = calc_mutual_info(bowtie_network, 
                                       bowtie_env_categories)
    info = Information(j)
    c_info = numpy.asarray(info)[0]

    # Need to account for different shape
    numpy.testing.assert_allclose(c_info[:, 0], p_info)
        
def get_causal_flow(net, matches):
    # graph.save_network_as_fullgraph(net, name='bob')
    w = net.constructor.world

    # For now, states are equiprobable
    env_count = len(w.environments)
    env_probs = numpy.ones(env_count) / env_count

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

    # All probabilities of co-occurences. A joint prob distn under
    # intervention for each channel, about each output channel.
    probs = numpy.zeros((w.reg_channels, w.out_channels, 2, 2))

    # Now, do the interventions for each regulatory gene
    # The genes are organised so that they begin with regulatory genes, so we
    # can simply grab those
    for i in range(w.reg_channels):
        gene = net.genes[i]
        # ----- GENE IS MANIPULATED OFF
        # (this automatically updates everything)
        gene.intervene = InterventionState.INTERVENE_OFF
        # Record everything the state of this genes output against each of the
        # outputs. Go through each environments, comparing output...
        for j, r in enumerate(net.rates):
            cats = bowtie_categorize_output(r, matches)
            for k, c in enumerate(cats):
                # What is the probability here? ith Gene is off in jth environment
                p_state = env_probs[j] * (1.0 - pdist_regs[i])
                # Gene is OFF. Output is category c
                probs[i, k, 0, c] += p_state

        # ----- GENE IS MANIPULATED ON
        # Manipulate the gene ON
        gene.intervene = InterventionState.INTERVENE_ON
        for j, r in enumerate(net.rates):
            cats = bowtie_categorize_output(r, matches)
            for k, c in enumerate(cats):
                # What is the probability here? ith Gene is ON in jth environment
                p_state = env_probs[j] * pdist_regs[i]
                # Gene is ON. Output is category c
                probs[i, k, 1, c] += p_state

        # Reset this gene
        gene.intervene = InterventionState.INTERVENE_NONE
    return probs

# Slow python calculation of information 
def calc_info_from_causal_flow(net, probs):
    w = net.constructor.world

    # Now calculate the information
    info = numpy.zeros((w.reg_channels, w.out_channels))
    for i in range(w.reg_channels):
        for j in range(w.out_channels):
            c_prob = probs[i, j]
            assert numpy.isclose(c_prob.sum(), 1.0)

            I = 0.0
            colsum = c_prob.sum(axis=0)
            rowsum = c_prob.sum(axis=1)
            for row in 0, 1:
                for col in 0, 1:
                    p_xy = c_prob[row, col]
                    p_x = rowsum[row]
                    p_y = colsum[col]
                    if p_xy != 0:
                        I += p_xy * logarithm(p_xy / (p_x * p_y), 2)
            info[i, j] = I

    return info

def test_causal_flow_cython(bowtie_network):
    """Compare the clunky python version with our C++ version"""

    # This what it was evolved to do (see the generate.py file)
    matches = [1, .5, .25]

    f = CausalFlowAnalyzer(bowtie_network.constructor.world, matches)
    j = f.analyse_network(bowtie_network)

    # We just access the 0th element, as we only sent one network for analysis
    c_joint = numpy.asarray(j)[0]

    # Test that each of the signals has a complete joint probability. It
    # should sum to 1.0
    for per_reg in c_joint:
        for per_output in per_reg:
            assert per_output.shape == (2, 2)
            assert per_output.sum() == 1.0

    # Get the slow version from python.
    p_joint = get_causal_flow(bowtie_network, matches)

    # They should be the same.
    numpy.testing.assert_allclose(p_joint, c_joint)

    info = Information(j)
    c_info = numpy.asarray(info)[0]
    p_info = calc_info_from_causal_flow(bowtie_network, p_joint)
    numpy.testing.assert_allclose(c_info, p_info)

def test_causal_flow_pop(bowtie_database, bowtie_env_categories):
    # This what it was evolved to do (see the generate.py file)
    matches = [1, .5, .25]
    pop = bowtie_database.population
    f_flow = CausalFlowAnalyzer(bowtie_database.world, matches)
    j_flow = f_flow.analyse_collection(pop)
    flow = Information(j_flow)
    c_flow = numpy.asarray(flow)

    f_mut = MutualInfoAnalyzer(bowtie_database.world, bowtie_env_categories)
    j_mut = f_mut.analyse_collection(pop)
    mut = Information(j_mut)
    c_mut = numpy.asarray(mut)

    # Let's just try the first few
    for i, n in enumerate(pop):
        p_joint = get_causal_flow(n, matches)
        p_flow = calc_info_from_causal_flow(n, p_joint)
        numpy.testing.assert_allclose(c_flow[i], p_flow)

        _, p_mut = calc_mutual_info(n, bowtie_env_categories)
        numpy.testing.assert_allclose(c_mut[i, :, 0], p_mut)

        # It takes too long....
        if i == 10:
            break
