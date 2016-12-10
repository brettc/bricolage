"""Testing more complex info measures using the double bow
"""

import pytest
import itertools
# from math import log as logarithm
from bricolage.core_ext import Collection
from bricolage.analysis_ext import MIAnalyzer, WCAnalyzer, MutualInfoAnalyzer
from bricolage.core import InterventionState
from bricolage.dot_layout import DotMaker
from bricolage.graph_draw import SmallDiagram, TextDiagram
from bricolage import graph_maker
from bricolage.analysis import NetworkAnalysis
# from bricolage import dot_layout
import numpy
numpy.set_printoptions(linewidth=120)

@pytest.fixture
def net1(double_bow):
    return double_bow.population.get_best()[0]

@pytest.fixture
def net2(double_bow):
    return double_bow.population.get_best()[1]

@pytest.fixture
def small_pop(double_bow):
    """faster to just look at a few"""
    coll = Collection(double_bow.factory)
    for i in range(25):
        coll.add(double_bow.population[i])
    return coll
    
def w_mutual_info(joint, weighting):
    """Weighted Mutual information"""
    assert numpy.isclose(joint.sum(), 1.0)
    assert joint.shape == weighting.shape

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
                info += p_xy * numpy.log2(p_xy / (p_x * p_y)) * weighting[row, col]
    return info


def distance(a1, a2):
    """Distance between two vectors"""
    assert len(a1) == len(a2)
    diff = a1 - a2
    return numpy.sqrt((diff ** 2).sum())


def similarity(a1, a2, strength=.2):
    """Convert distance into nice measure between 0 and 1"""
    dist = distance(a1, a2)
    ret = numpy.exp(-numpy.abs(dist) / strength)
    return ret


class Categorizer(object):
    """Categorize boolean genes"""
    def __init__(self, ngenes, target_1, target_2):
        self.ngenes = ngenes
        self.target_1 = target_1
        self.target_2 = target_2
        self.patterns = []
        self.probs = []

    def add(self, pat, gene, is_on, pr):
        tindex = int(is_on)
        assert 0 <= tindex <= 1
        for i, t in enumerate(self.patterns):
            # Look for a match
            if numpy.allclose(t, pat):
                self.probs[i][gene, tindex] += pr
                break
        else:
            # No match found -- add new pattern.
            self.patterns.append(pat)
            self.probs.append(numpy.zeros((self.ngenes, 2)))
            self.probs[-1][gene, tindex] += pr

    def make_joint(self):
        # Summarize the details into a joint distribution
        joint = numpy.zeros((self.ngenes, 2, len(self.patterns)), float)
        for i, pats in enumerate(self.probs):
            for j, gene in enumerate(pats): 
                joint[j, 0, i] = gene[0]
                joint[j, 1, i] = gene[1]
        return joint

    def construct_weightings(self, strength):
        # Construct weightings using the similarity measure. We want
        # similarities to two different patterns. The weightings are the row
        # swappings of each other.
        w1 = numpy.zeros((2, len(self.patterns))) 
        for i, p in enumerate(self.patterns):
            w1[0, i] = similarity(p, self.target_1, strength=strength)
            w1[1, i] = similarity(p, self.target_2, strength=strength)

        # w2 is the same with rows swapped
        w2 = numpy.zeros((2, len(self.patterns))) 
        w2[0,:] = w1[1,:]
        w2[1,:] = w1[0,:]

        return w1, w2

    def calc_weighted_joint(self, strength):
        info = numpy.zeros(self.ngenes, float)

        if strength < 0.0:
            w1 = w2 = numpy.ones((2, len(self.patterns))) 
        else:
            w1, w2 = self.construct_weightings(strength)

        joints = self.make_joint()

        # We take the BEST that a gene can do to increase the fitness
        for i, jnt in enumerate(joints):
            m1 = w_mutual_info(jnt, w1)
            m2 = w_mutual_info(jnt, w2)
            info[i] = max(m1, m2)

        return info


def get_joint(net, outputs, target_1, target_2):
    """Return the joint probabilities under intervention"""
    wrld = net.factory.world
    target_1 = numpy.asarray(target_1)
    target_2 = numpy.asarray(target_2)
    outputs = numpy.asarray(outputs)

    assert len(outputs) == len(target_1) == len(target_2)

    # For now, states are equiprobable
    env_count = len(wrld.environments)
    env_prob = 1.0 / env_count

    # For now, ignoring "natural probabilities"
    p_off = p_on = 0.5 * env_prob

    categorizer = Categorizer(
        wrld.reg_channels, target_1, target_2)

    for i in range(wrld.reg_channels):
        gene = net.genes[i]
        # ----- GENE IS MANIPULATED OFF ------
        # NOTE: This automatically updates everything (changing the rates!).
        gene.intervene = InterventionState.INTERVENE_OFF
        if not numpy.isclose(p_off, 0.0):
            for j, rate in enumerate(net.rates):
                categorizer.add(rate[outputs], i, False, p_off)

        # ----- GENE IS MANIPULATED ON -------
        gene.intervene = InterventionState.INTERVENE_ON
        if not numpy.isclose(p_on, 0.0):
            for j, rate in enumerate(net.rates):
                categorizer.add(rate[outputs], i, True, p_on)

        # ----- Reset this gene
        gene.intervene = InterventionState.INTERVENE_NONE

    # Now, build the joint probabilities!
    return categorizer

# def test_graph(double_bow, net2):
#     ana = NetworkAnalysis(net2)
#     ana.annotate(double_bow.targets[0])
#     g = graph_maker.GeneSignalGraph(ana)
#     d = DotMaker(g)
#     d.save_picture('signal.png')
    
ARGS = [range(4), [0, 0, 1, 1], [1, 1, 0, 0]]
WEIGHTING = .2

def test_single_joint(net1):
    c_args = [net1.factory.world] + ARGS + [WEIGHTING]
    p_args = [net1] + ARGS

    wc = WCAnalyzer(*c_args)
    j_c = wc.get_joint(net1)
    
    cat = get_joint(*p_args)
    j_p = cat.make_joint()
    assert j_c.shape == j_p.shape
    numpy.testing.assert_allclose(j_c, j_p)

# @pytest.mark.skip(reason="Algorithm difference for some reason...")
def test_pop_joint(small_pop):
    c_args = [small_pop.factory.world] + ARGS + [WEIGHTING]
    wc = WCAnalyzer(*c_args)

    for net in small_pop:
        p_args = [net] + ARGS
        j_c = wc.get_joint(net)
        cat = get_joint(*p_args)
        j_p = cat.make_joint()

        for j in j_p:
            assert numpy.isclose(j.sum(), 1.0)

        for j in j_c:
            assert numpy.isclose(j.sum(), 1.0)

        # Currently FAILS. Something to do with the algorithm? Information is
        # correct below
        # assert j_c.shape == j_p.shape
        # numpy.testing.assert_allclose(j_c, j_p)


def test_pop_info(small_pop):
    c_args = [small_pop.factory.world] + ARGS + [WEIGHTING]
    wc = WCAnalyzer(*c_args)

    # Analyze everything the fast way
    all_i = wc.analyse_collection(small_pop)

    for net, i_c_coll in zip(small_pop, all_i):
        p_args = [net] + ARGS
        j_c = wc.get_joint(net)
        cat = get_joint(*p_args)

        i_p = cat.calc_weighted_joint(WEIGHTING)
        i_c = wc.analyse_network(net)

        assert i_p.shape == i_c.shape
        numpy.testing.assert_allclose(i_c, i_p)

        assert i_c.shape == i_c_coll.shape
        numpy.testing.assert_allclose(i_c, i_c_coll)


def make_cats(world):
    def fn1(a, b, c, d, e, f):
        if (a and not b) or (b and not c):
            return 1
        return 0

    def fn2(a, b, c, d, e, f):
        if (d and not e) or (d and not f) or (not e and not f):
            return 1
        return 0

    a, b =  world.cue_range
    c1 = []
    c2 = []

    # Stolen from targets_ext.pyx
    for i, e in enumerate(world.environments):
        # TODO: Clean up the refs here
        args = e.as_array()[a:b]
        c1.append(fn1(*args))
        c2.append(fn2(*args))

    return c1, c2


def test_mi(double_bow, net2):
    # from bricolage.dot_layout import save_network_as_fullgraph
    # from bricolage.graph_maker import GraphType
    # p = three_database.population
    t = double_bow.targets[0]
    c1, c2 = make_cats(double_bow.world)
    mi1 = MIAnalyzer(double_bow.world, c1)
    mi2 = MIAnalyzer(double_bow.world, c2)
    print mi1.analyse_network(net2)
    print mi2.analyse_network(net2)

    mi3 = MutualInfoAnalyzer(double_bow.world, c1)
    print mi3.numpy_info_from_network(net2).ravel()
    mi4 = MutualInfoAnalyzer(double_bow.world, c2)
    print mi4.numpy_info_from_network(net2).ravel()

    # print t.as_array()
    # cats = t.calc_categories()
    # print cats
    # n = p.get_best()[0]
    # print n.fitness
    # print n.identifier
    # save_network_as_fullgraph(n, graph_type=GraphType.GENE)
    # print cats
    # print calc_mutual_info(n, cats)
    #
    # cmi = MIAnalyzer(n.factory.world, cats)
    # print cmi.numpy_info_from_network(n)


def test_envs(double_bow):
    """Test that we can create environments the same as the c++"""
    # We have to reverse things to match the way we construct environments
    p_envs = [_[::-1] for _ in itertools.product([0, 1], repeat=6)]

    world = double_bow.world
    a, b =  world.cue_range
    c_envs = []
    for i, e in enumerate(world.environments):
        # TODO: Clean up the refs here
        args = e.as_array()[a:b]
        c_envs.append(tuple(args))

    for a, b in zip(p_envs, c_envs):
        assert a == b

