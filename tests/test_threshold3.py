import cPickle as pickle
from bricolage.threshold3 import (
    World, Parameters, Factory, Population,
    SelectionModel,
)

from bricolage.core import InputType, ScoringMethod
from bricolage.targets_ext import NoisyTarget, DefaultTarget


def xor_func(a, b):
    return (a or b) and not (a and b)


def fitness_func1(a, b):
    return xor_func(a, b)


def test_attractors():
    params = Parameters(seed=4, cis_count=2, reg_channels=5, out_channels=2,
                        cue_channels=3,)
    world = World(params)
    factory = Factory(world)
    net = factory.create_network()
    alens = [len(a) for a in net.attractors]
    print alens
    print net.attractors_size


def test_network_pickle():
    params = Parameters(seed=4, cis_count=2, reg_channels=5, out_channels=2,
                        cue_channels=3,)
    world = World(params)
    factory = Factory(world)
    n1 = factory.create_network()
    out = pickle.dumps(n1, -1)
    n2 = pickle.loads(out)
    for g1, g2 in zip(n1.genes, n2.genes):
        assert g1.pub == g2.pub
        for m1, m2 in zip(g1.modules, g2.modules):
            assert m1.same_as(m2)


def test_population():
    params = Parameters(seed=4, cis_count=2, reg_channels=5, out_channels=2,
                        cue_channels=3,)
    world = World(params)
    factory = Factory(world)
    popul = Population(factory, 1000)
    print popul[0].attractors[1]


def test_xor():
    """We should be able to find xor!"""
    params = Parameters(seed=2, cis_count=2, cue_channels=2, reg_channels=1,
                        out_channels=1)
    world = World(params)
    factory = Factory(world)
    target = DefaultTarget(world, fitness_func1)
    select = SelectionModel(world)
    pop = Population(factory, 10000)

    while 1:
        pop.assess(target)
        pop.select(select)
        fits = pop.fitnesses
        max_fit = max(fits)
        if max_fit == 1.0:
            break
        pop.mutate(.05)

    # net = pop.get_best()[0]
    # dot_layout.save_network_as_fullgraph(net, name='xor')



def fitness_func2(a, b):
    return xor_func(a, b)


def test_stabilise():
    # params = Parameters(seed=3, input_type=InputType.PULSE, cis_count=3,
    params = Parameters(seed=5, input_type=InputType.PULSE, cis_count=3,
                        cue_channels=2, reg_channels=4, out_channels=1)
    world = World(params)
    factory = Factory(world)
    target = DefaultTarget(world, fitness_func2,
                           scoring_method=ScoringMethod.EXPONENTIAL_VEC,
                           strength=.2)
    ntarget = NoisyTarget(world, fitness_func2,
                           scoring_method=ScoringMethod.EXPONENTIAL_VEC,
                           strength=.2, perturb_count=1, perturb_prop=.2)
    select = SelectionModel(world)
    pop = Population(factory, 5000)
    while 1:
        pop.assess(ntarget)
        w, b = pop.worst_and_best()
        print b
        if b == 1.0:
            break
        pop.select(select)
        pop.mutate(.002)

    # Get a network!
    n = pop.get_best()[0]

    print '-->', target.assess(n)
    print '-->', ntarget.assess(n)

    # Get an attractor state
    orig = n.attractors[2][0]
    perturbed = orig.copy()

    # Flip it out of that state
    # See if it goes back in to the attractor
    perturbed.flip(5)
    trans, attr, rates = n.stabilise(perturbed)
    assert orig != attr[0]

    # Push it HARDER
    perturbed.flip(6)
    trans, attr, rates = n.stabilise(perturbed)
    assert orig != attr[0]
