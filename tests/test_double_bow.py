"""Testing more complex info measures using the double bow
"""

import pytest
# from math import log as logarithm
from bricolage.core_ext import Collection
from bricolage.analysis_ext import MIAnalyzer, WCAnalyzer
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


def test_graph(double_bow, net2):
    ana = NetworkAnalysis(net2)
    ana.annotate(double_bow.targets[0])
    g = graph_maker.GeneSignalGraph(ana)
    d = DotMaker(g)
    d.save_picture('signal.png')
    # print n.fitness
    # tset = t.calc_distinct_outputs()
    # py_info = get_relevant_control(net, tset)
    # rz = RelevantControlAnalyzer(net.factory.world, tset)
    # cy_info = rz.numpy_info_from_network(net)
    # numpy.testing.assert_allclose(py_info, cy_info)
    

def test_basic(double_bow, net2):
    cat = get_joint(net2, range(4), [0, 0, 1, 1], [1, 1, 0, 0])
    print cat.calc_weighted_joint(.25)
    cat = get_joint(net2, range(4, 8), [0, 0, 1, 1], [1, 1, 0, 0])
    print cat.calc_weighted_joint(.25)
    # print cat.calc_weighted_joint(-.5)

def test_cpp(net1):
    wc = WCAnalyzer(net1.factory.world, range(4), [0, 0, 1, 1], [1, 1, 0, 0], .2)
    # j_c = wc.get_joint(net1)
    #
    cat = get_joint(net1, range(4), [0, 0, 1, 1], [1, 1, 0, 0])
    # j_p = cat.make_joint()
    # assert j_c.shape == j_p.shape
    # numpy.testing.assert_allclose(j_c, j_p)

    j_i = wc.analyse_network(net1)
    print j_i
    print cat.calc_weighted_joint(.2)

def test_array(double_bow):
    wc = WCAnalyzer(double_bow.world, range(4), [0, 0, 1, 1], [1, 1, 0, 0], .2)

    # net = double_bow.population[0]
    # j_c = wc.analyse_network(net)
    # print j_c
    coll = Collection(double_bow.factory)
    for i in range(20):
        coll.add(double_bow.population[i])
    #
    j_c = wc.analyse_collection(coll)
    print j_c[0]
    print j_c[1]

