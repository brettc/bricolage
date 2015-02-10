import pytest
import pathlib
from math import log as logarithm
from bricolage.core_ext import EnvironmentI
from bricolage import threshold3, lineage, graph
import numpy

def feature1(a, b):
    """compute an xor function"""
    if (a or b) and not (a and b):
        return 0
    return 1

rates = [0, 0], [0.5, 1]

def target1(a, b):
    """compute an xor function"""
    return rates[feature1(a, b)]

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
        a.add_target(target1)
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

def numpy_compare(a, b):
    assert max(a - b) < 0e-10

def test_cython_pop_info(lineage_2x2):
    pop = lineage_2x2.population
    w = lineage_2x2.world

    e = EnvironmentI(w, feature1)

    # Calculate it via cython
    c_probs = e.calc_collection(pop)
    c_info = e.calc_info(pop)

    # Calculate it via cython
    c_probs_net = [] 
    for n in pop:
        c_probs_net.append(e.calc_network(n))

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
        print e, target1(*inputs), r
        # for i in range(5)
    
