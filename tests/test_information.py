import pytest
import pathlib
from math import log as logarithm
from bricolage.core_ext import InfoE
from bricolage import threshold3, lineage, graph
from bricolage.core import InterventionState
import numpy

def feature1(a, b):
    """compute an xor function"""
    if (a or b) and not (a and b):
        return 0
    return 1

# Output rates for depending on environment
rates = [0, 0], [0.5, 1]

def xor_target(a, b):
    """compute an xor function"""
    return rates[feature1(a, b)]

def bowtie_target(a, b, c):
    if (a and not c) or (b and c):
        return [1, .5, .25]
    return [0, 0, 0]

def bowtie_categorize(a, b, c):
    ca = cb = cc = 0
    if numpy.isclose(a, 1.0):
        ca = 1
    if numpy.isclose(b, .5):
        cb = 1
    if numpy.isclose(c, .25):
        cc = 1
    return [ca, cb, cc]

@pytest.fixture
def p_2x2():
    return threshold3.Parameters(
        seed=3, 
        cis_count=2, 
        reg_channels=4,
        out_channels=2, 
        cue_channels=2, 
        population_size=1000,
        mutation_rate=.001,
    )

@pytest.fixture
def lineage_2x2(p_2x2, tmpdir):
    # base = pathlib.Path(str(tmpdir))
    base = pathlib.Path('.')
    path = base / 'env_info.db'
    if path.exists():
        a = lineage.FullLineage(path)
    else:
        a = lineage.FullLineage(path, params=p_2x2)
        a.add_target(xor_target)
        while 1:
            a.next_generation()
            w, b = a.population.worst_and_best()
            if b == 1.0:
                break
    return a

@pytest.fixture
def complex_network1(lineage_2x2):
    b = lineage_2x2.population.best_indexes()[0]
    return lineage_2x2.population[b]

@pytest.fixture
def p_3x3():
    return threshold3.Parameters(
        seed=8, 
        cis_count=2, 
        reg_channels=8,
        out_channels=3, 
        cue_channels=3, 
        population_size=1000,
        mutation_rate=.002,
    )

@pytest.fixture
def lineage_3x3(p_3x3, tmpdir):
    # base = pathlib.Path(str(tmpdir))
    base = pathlib.Path('.')
    path = base / '3x3_info8.db'
    cnt = 0
    if path.exists():
        a = lineage.FullLineage(path)
    else:
        a = lineage.FullLineage(path, params=p_3x3)
        a.add_target(bowtie_target)
        while 1:
            a.next_generation()
            w, b = a.population.worst_and_best()
            if b == 1.0:
                cnt += 1
            if cnt == 1000:
                break
    return a

@pytest.fixture
def complex_bowtie(lineage_3x3):
    b = lineage_3x3.population.best_indexes()[0]
    return lineage_3x3.population[b]

def calc_info_for_network(n, func):
    w = n.constructor.world
    cr = range(*w.cue_range)
    features = []

    for e in w.environments:
        inputs = [e[i] for i in cr]
        features.append(func(*inputs))
    assert len(features) == len(w.environments)

    # Features should be consecutive numbers
    all_feat = set(features)
    assert all_feat == set(range(len(all_feat)))

    reg_base, reg_to = w.reg_range
    channel_dim = reg_to - reg_base
    feat_dim = len(set(features))
    state_dim = 2 # on or off
    env_dim = len(w.environments)

    # Now we have the dimensions of our array
    probs = numpy.zeros((channel_dim, feat_dim, state_dim))

    for i, channel in enumerate(range(*w.reg_range)):
        for attrs, feat in zip(n.attractors, features):
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

def test_cython_pop_info(lineage_2x2):
    pop = lineage_2x2.population
    w = lineage_2x2.world

    e = InfoE(w, feature1)

    # Calculate it via cython
    c_probs = e.collection_probs(pop)
    c_info = e.collection_info(pop)

    # Calculate it via cython
    c_probs_net = [] 
    for n in pop:
        c_probs_net.append(e.network_probs(n))

    # Calculate it using the test function
    p_probs = []
    p_info = []
    for n in pop:
        p, i = calc_info_for_network(n, feature1)
        p_probs.append(p)
        p_info.append(i)

    for c, cn, pn in zip(c_probs, c_probs_net, p_probs):
        numpy.testing.assert_allclose(c, cn)
        numpy.testing.assert_allclose(c, pn)

    for c, pn in zip(c_info, p_info):
        numpy.testing.assert_allclose(c, pn)


def test_env_info(complex_network1):
    n = complex_network1
    graph.save_network_as_fullgraph(n, name='bob')
    w = n.constructor.world
    for e, r in zip(w.environments, n.rates):
        inputs = [e[i] for i in range(*w.cue_range)]
        print e, xor_target(*inputs), r
        # for i in range(5)

# def test_target_info(complex_network1, lineage_2x2):
def test_target_info(complex_bowtie, lineage_3x3):
    net = complex_bowtie
    lin = lineage_3x3
    graph.save_network_as_fullgraph(net, name='bob')
    w = net.constructor.world

    # For now, states are equiprobable
    env_count = len(w.environments)
    env_probs = numpy.ones(env_count) / env_count

    # First, calculate the "natural" probabilites of each regulatory signals
    # without intervention
    pdist_regs = numpy.zeros(w.reg_channels)
    reg_base, reg_to = w.reg_range
    for env_i, attr in enumerate(net.attractors):
        p_state = env_probs[env_i] / float(len(attr))
        # For each attractor ...
        for st in attr:
            # For all attractors states, and each regulatory channel
            for i in range(w.reg_channels):
                # If it is on, then add the probability
                if st.test(reg_base + i):
                    pdist_regs[i] += p_state

    # Get the target
    target = lin.targets[0].as_array()
    print target

    # All probabilities of co-occurences. A joint prob distn under
    # intervention for each channel, about each output channel.
    # TODO: second to last channel should be categorise MORE THAN 2
    probs = numpy.zeros((w.reg_channels, w.out_channels, 2, 2))

    print 'dist', pdist_regs

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
        for r in net.rates:
            cats = bowtie_categorize(*r)
            for j, c in enumerate(cats):
                # What is the probability here? ith Gene is off in jth environment
                p_state = env_probs[j] * (1.0 - pdist_regs[i])
                # Gene is OFF. Output is category c
                probs[i, j, 0, c] += p_state

        # ----- GENE IS MANIPULATED ON
        # Manipulate the gene ON
        gene.intervene = InterventionState.INTERVENE_ON
        for r in net.rates:
            cats = bowtie_categorize(*r)
            for j, c in enumerate(cats):
                # What is the probability here? ith Gene is ON in jth environment
                p_state = env_probs[j] * pdist_regs[i]
                # Gene is ON. Output is category c
                probs[i, j, 1, c] += p_state

        # Reset this gene
        gene.intervene = InterventionState.INTERVENE_NONE

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

    print info

    # target = lineage_2x2.targets[0]
    # target_outputs = target.as_array()
    
def test_manip(complex_bowtie):
    net = complex_bowtie
    graph.save_network_as_fullgraph(net, name='unbob', simplify=False)
    print net.rates
    net.genes[7].intervene = InterventionState.INTERVENE_OFF
    print net.rates
    net.genes[7].intervene = InterventionState.INTERVENE_ON
    print net.rates

