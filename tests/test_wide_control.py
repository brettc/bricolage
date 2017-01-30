"""Testing the final simple measures
"""

import pytest
from bricolage.core_ext import Collection
from bricolage.analysis_ext import FastCandBAnalyzer
from bricolage.core import InterventionState
import numpy
numpy.set_printoptions(linewidth=150)

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
    

def mutual_info(joint, rowsum):
    """Mutual information"""
    # assert numpy.isclose(joint.sum(), 1.0)
    assert rowsum.shape[0] == joint.shape[0]
    assert rowsum.sum() == 1.0

    info = 0.0
    rownum, colnum = joint.shape
    colsum = joint.sum(axis=0)
    # rowsum = joint.sum(axis=1)
    for row in range(rownum):
        for col in range(colnum):
            p_xy = joint[row, col]
            p_x = rowsum[row]
            p_y = colsum[col]
            if p_xy != 0:
                info += p_xy * numpy.log2(p_xy / (p_x * p_y)) 
    return info


def get_fast_info(net, outputs, target_1, target_2):
    """Return the joint probabilities under intervention"""
    wrld = net.factory.world
    target_1 = numpy.asarray(target_1)
    target_2 = numpy.asarray(target_2)
    outputs = numpy.asarray(outputs)

    non_outputs = [i for i in range(wrld.out_channels) if i not in outputs]
    non_outputs = numpy.asarray(non_outputs)

    assert len(outputs) == len(target_1) == len(target_2)
    assert len(non_outputs) == wrld.out_channels - len(outputs)

    # For now, states are equiprobable
    env_count = len(wrld.environments)
    env_prob = 1.0 / env_count

    # For now, ignoring "natural probabilities"
    p_off = p_on = 0.5 * env_prob

    joint = numpy.zeros((wrld.reg_channels, 2, 2))
    all_rates = numpy.zeros((2, env_count, wrld.out_channels))

    for i in range(wrld.reg_channels):
        gene = net.genes[i]
        # ----- GENE IS MANIPULATED OFF ------
        gene.intervene = InterventionState.INTERVENE_OFF
        all_rates[0] = net.rates
        # ----- GENE IS MANIPULATED ON -------
        gene.intervene = InterventionState.INTERVENE_ON
        all_rates[1] = net.rates
        # ----- Reset this gene
        gene.intervene = InterventionState.INTERVENE_NONE

        for j, (off_rate, on_rate) in enumerate(zip(all_rates[0], all_rates[1])):
            # No change on the non-target rates?
            if numpy.allclose(off_rate[non_outputs], on_rate[non_outputs]):
                if numpy.allclose(off_rate[outputs], target_1):
                    joint[i, 0, 0] += p_off
                elif numpy.allclose(off_rate[outputs], target_2):
                    joint[i, 0, 1] += p_off

                if numpy.allclose(on_rate[outputs], target_1):
                    joint[i, 1, 0] += p_on
                elif numpy.allclose(on_rate[outputs], target_2):
                    joint[i, 1, 1] += p_on

        # ----- Reset this gene
        gene.intervene = InterventionState.INTERVENE_NONE

    # Now, build the joint probabilities!
    # Pass in the row-sum (our intervention probabilities), so that the mutual
    # info score is weighted
    rowsum = numpy.asarray([.5, .5])
    info = numpy.zeros(wrld.reg_channels)
    for i, target_j in enumerate(joint):
        info[i] = mutual_info(target_j, rowsum)

    return joint, info


class Background(object):
    def __init__(self, pat):
        self.pat = pat
        self.row_sum = numpy.zeros(1)
        self.joint = numpy.zeros((1, 2))


def get_target_and_background_info(net, outputs, target_1, target_2):
    """Return the joint probabilities under intervention"""
    wrld = net.factory.world
    target_1 = numpy.asarray(target_1)
    target_2 = numpy.asarray(target_2)
    outputs = numpy.asarray(outputs)

    non_outputs = [i for i in range(wrld.out_channels) if i not in outputs]
    non_outputs = numpy.asarray(non_outputs)

    assert len(outputs) == len(target_1) == len(target_2)
    assert len(non_outputs) == wrld.out_channels - len(outputs)

    # For now, states are equiprobable
    env_count = len(wrld.environments)
    env_prob = 1.0 / env_count

    # For now, ignoring "natural probabilities"
    p_off = p_on = 0.5 * env_prob

    joint = numpy.zeros((wrld.reg_channels, 2, 2))
    observe_bg = [dict() for _ in range(wrld.reg_channels)]

    all_rates = numpy.zeros((2, env_count, wrld.out_channels))


    for i in range(wrld.reg_channels):
        gene = net.genes[i]
        # ----- GENE IS MANIPULATED OFF ------
        gene.intervene = InterventionState.INTERVENE_OFF
        all_rates[0] = net.rates
        # ----- GENE IS MANIPULATED ON -------
        gene.intervene = InterventionState.INTERVENE_ON
        all_rates[1] = net.rates
        # ----- Reset this gene
        gene.intervene = InterventionState.INTERVENE_NONE

        for j, (off_rate, on_rate) in enumerate(zip(all_rates[0], all_rates[1])):
            target_rates = off_rate[outputs]
            bg_key = tuple(off_rate[non_outputs])

            bg = observe_bg[i].setdefault(bg_key, Background(bg_key))
            bg.row_sum[0] += p_off

            # No change on the non-target rates?
            if numpy.allclose(target_rates, target_1):
                joint[i, 0, 0] += p_off
                bg.joint[0, 1] += p_off
            elif numpy.allclose(target_rates, target_2):
                joint[i, 0, 1] += p_off
                bg.joint[0, 1] += p_off

            target_rates = on_rate[outputs]
            bg_key = tuple(on_rate[non_outputs])

            bg = observe_bg[i].setdefault(bg_key, Background(bg_key))
            bg.row_sum[0] += p_on

            if numpy.allclose(target_rates, target_1):
                joint[i, 1, 0] += p_on
                bg.joint[0, 0] += p_on
            elif numpy.allclose(target_rates, target_2):
                joint[i, 1, 1] += p_on
                bg.joint[0, 1] += p_on

        # ----- Reset this gene
        gene.intervene = InterventionState.INTERVENE_NONE

    # Now, build the joint probabilities!
    # Pass in the row-sum (our intervention probabilities), so that the mutual
    # info score is weighted
    rowsum = numpy.asarray([.5, .5])
    target_info = numpy.zeros(wrld.reg_channels)
    bg_info = numpy.zeros(wrld.reg_channels)
    for i, target_j in enumerate(joint):
        target_info[i] = mutual_info(target_j, rowsum)
        bg = observe_bg[i]
        bg_rowsum = numpy.concatenate([b.row_sum for b in bg.values()])
        bg_joint = numpy.concatenate([b.joint for b in bg.values()])
        bg_info[i] = mutual_info(bg_joint, bg_rowsum)

        print i
        print target_j
        print target_info[i]
        print [b.pat for b in bg.values()]
        print bg_rowsum
        print bg_joint
        print bg_info[i]
        print '--------'


    return target_info, bg_info

    
ARGS = [range(4), [0, 0, 1, 1], [1, 1, 0, 0]]

def test_single_joint(net1):
    c_args = [net1.factory.world] + ARGS
    p_args = [net1] + ARGS

    fc = FastCandBAnalyzer(*c_args)
    j_p, p_info = get_fast_info(*p_args)
    c_info = fc.analyse_network(net1)
    print
    print p_info
    print c_info
    numpy.testing.assert_allclose(p_info, c_info)


def test_pop_info(small_pop):
    c_args = [small_pop.factory.world] + ARGS 
    fc = FastCandBAnalyzer(*c_args)

    # Analyze everything the fast way
    all_info = fc.analyse_collection(small_pop)

    for net, c_pop_info in zip(small_pop, all_info):
        p_args = [net] + ARGS
        j, p_info = get_fast_info(*p_args)
        c_info = fc.analyse_network(net)
        print c_info

        numpy.testing.assert_allclose(p_info, c_info)
        numpy.testing.assert_allclose(c_pop_info, p_info)


