"""Testing the final simple measures
"""

import pytest
from bricolage.core_ext import Collection
from bricolage.analysis_ext import FastCAnalyzer
from bricolage.core import InterventionState
import numpy

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
    
def w_mutual_info(joint, rowsum):
    """Mutual information"""
    # assert numpy.isclose(joint.sum(), 1.0)

    info = 0.0
    rownum, colnum = joint.shape
    colsum = joint.sum(axis=0)
    for row in range(rownum):
        for col in range(colnum):
            p_xy = joint[row, col]
            p_x = rowsum[row]
            p_y = colsum[col]
            if p_xy != 0:
                info += p_xy * numpy.log2(p_xy / (p_x * p_y)) 
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

    joint = numpy.zeros((wrld.reg_channels, 2, 2))

    for i in range(wrld.reg_channels):
        gene = net.genes[i]
        # ----- GENE IS MANIPULATED OFF ------
        # NOTE: This automatically updates everything (changing the rates!).
        gene.intervene = InterventionState.INTERVENE_OFF
        for j, rate in enumerate(net.rates):
            if numpy.allclose(rate[outputs], target_1):
                joint[i, 0, 0] += p_off
            elif numpy.allclose(rate[outputs], target_2):
                joint[i, 0, 1] += p_off


        # ----- GENE IS MANIPULATED ON -------
        gene.intervene = InterventionState.INTERVENE_ON
        for j, rate in enumerate(net.rates):
            if numpy.allclose(rate[outputs], target_1):
                joint[i, 1, 0] += p_on
            elif numpy.allclose(rate[outputs], target_2):
                joint[i, 1, 1] += p_on

        # ----- Reset this gene
        gene.intervene = InterventionState.INTERVENE_NONE

    # Now, build the joint probabilities!
    # Pass in the row-sum (our intervention probabilities), so that the mutual
    # info score is weighted
    rowsum = numpy.asarray([.5, .5])
    info = numpy.zeros(wrld.reg_channels)
    for i, j in enumerate(joint):
        info[i] = w_mutual_info(j, rowsum)

    return joint, info


ARGS = [range(4), [0, 0, 1, 1], [1, 1, 0, 0]]

def test_single_joint(net1):
    c_args = [net1.factory.world] + ARGS
    p_args = [net1] + ARGS

    fc = FastCAnalyzer(*c_args)
    j_p, p_info = get_joint(*p_args)

    c_info = fc.analyse_network(net1)
    numpy.testing.assert_allclose(p_info, c_info)


def test_pop_info(small_pop):
    c_args = [small_pop.factory.world] + ARGS 
    fc = FastCAnalyzer(*c_args)

    # Analyze everything the fast way
    all_info = fc.analyse_collection(small_pop)

    for net, c_pop_info in zip(small_pop, all_info):
        p_args = [net] + ARGS
        j, p_info = get_joint(*p_args)

        c_info = fc.analyse_network(net)

        numpy.testing.assert_allclose(p_info, c_info)
        numpy.testing.assert_allclose(c_pop_info, p_info)


