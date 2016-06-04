import pytest
from bricolage.threshold3 import Parameters
from bricolage.core import InputType
from bricolage.dot_layout import save_network_as_fullgraph
import numpy.random as rand
import numpy

@pytest.fixture
def p_3x1():
    return Parameters(
        seed=1,
        cis_count=2,
        reg_channels=4,
        out_channels=1,
        cue_channels=2,
        population_size=1000,
        input_type=InputType.PULSE,
        mutation_rate=.001,
        replicates=10,
    )

def target1(a, b):
    if (a or b) and not (a and b):
        return 1
    return 0

def test_robustness(bowtie_network):
    net = bowtie_network
    save_network_as_fullgraph(net)
    w = net.factory.world
    low, high = w.reg_range

    print w.reg_range
    for e in w.environments:
        orig = net.stabilise(e)[2]
        i = rand.randint(low, high)
        e.flip(i)
        pert = net.stabilise(e)[2] 
        print orig, pert

def test_caching(bowtie_network):
    # TODO: check that all of the transient and attractors stuff is in the
    # cache
    net = bowtie_network
    # print len(net.rates_cache)


    # w = net.factory.world
    # for e in w.environments:
    #     a, t, r = net.stabilise(e)
    #     assert (r == net.get_rates(e)).all()
    #
    print len(net.rates_cache)
    ks = net.rates_cache.keys()
    ks.sort()
    for k in ks:
        print k, net.rates_cache[k]
    
    # print net.attractors
    #
def get_robustness_attr(net):
    w = net.factory.world
    low, high = w.cue_range[0], w.reg_range[1]
    per_channel = 1.0 / float((high - low) * len(w.environments))
    rob = 0.0

    # For each env
    for env in w.environments:
        attr, trans, orig_rates = net.stabilise(env)

        # These better be the same
        other_orig_rates = net.get_rates(env)
        numpy.testing.assert_allclose(orig_rates, other_orig_rates)

        per_test = per_channel / float(len(attr))

        # For each state in the attractor
        for attr_state in attr:

            # For each possible regulatory channel
            for i in range(low, high):
                ec = attr_state.copy()
                ec.flip(i)
                new_rates = net.get_rates(ec)
                nattr, ntrans, alt_new_rates = net.stabilise(ec)
                numpy.testing.assert_allclose(alt_new_rates, new_rates)

                if numpy.allclose(orig_rates, new_rates):
                    rob += per_test

    return rob

def test_robustness_attr(bowtie_network):
    net = bowtie_network
    py_rob = get_robustness_attr(net)
    cy_rob = net.attractor_robustness
    assert py_rob == cy_rob


def test_robustness_pop(bowtie_database):
    pop = bowtie_database.population

    x = pop.robustness()
    k = numpy.isclose(x, 1.0)
    print numpy.where(k)[0]

# def test1(tmpdir, p_3x1):
#     base = Path(str(tmpdir))
#     path = base / 'test1.db'
#     with FullLineage(path, p_3x1) as a:
#         a.add_target(target1)
#         for i in range(100):
#             a.next_generation()
#
#     net = a.population.get_best()[0]
#     print net.fitness
#     print net.attractors
#     for i in range(5):
#         net.calc_perturbation(False)
#         print net.pert_attractors
#
# def test2(bowtie_database):
#     db = bowtie_database
#     net = db.population.get_best()[0]
#     print net.fitness
#     # print net.attractors[0]
#     print net.rates[0], net.rates[1]
#     print '---'
#     for i in range(150):
#         net.calc_perturbation(False)
#     for i in range(5):
#         net.calc_perturbation(False)
#         print net.pert_rates[0], net.pert_rates[1]
#         # print net.pert_attractors[0]
        
def perturb_target(a, b, c):
    if (a and not c) or (b and c):
        return [0, 1, 0]
    return [1, 0, 1]
        
# def test3(perturb_database):
#     db = perturb_database
#     t1 = DefaultTarget(db.world, perturb_target)
#     t2 = NoisyTarget(db.world, perturb_target)
#     print t1.assess_collection(db.population)[:10]
#     print t2.assess_collection(db.population)[:10]
