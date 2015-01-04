from bricolage.threshold3 import (
    World, Parameters, Target, Constructor, Network, Population)

def xor_func(a, b):
    return (a or b) and not (a and b)

def fitness_func1(a, b):
    return xor_func(a, b)

def test_network():
    params = Parameters(seed=4, cis_count=2, reg_channels=5, out_channels=2,
                        cue_channels=3,)
    world = World(params)
    const = Constructor(world, params)
    net = Network(const)
    print net.attractors
    print world.pub_signals
    # net.genes[0].pub = 3
    #

def test_population():
    params = Parameters(seed=4, cis_count=2, reg_channels=5, out_channels=2,
                        cue_channels=3,)
    world = World(params)
    const = Constructor(world, params)
    popul = Population(const, 1000)
    print popul[0].attractors[1]

def test_xor():
    params = Parameters(seed=2, cis_count=2, cue_channels=2, reg_channels=0,
                        out_channels=1)
    world = World(params)
    print world.environments
    const = Constructor(world, params)
    target = Target(world, fitness_func1)
    pop = Population(const, 10000)

    while 1:
        pop.select(target)
        fits = [n.fitness for n in pop]
        max_fit = max(fits)
        if max_fit == 1.0:
            break
        pop.mutate(.05)

    # print w.environments
    for n in pop:
        if n.fitness == 1.0:
            print n
            for g in n.genes:
                for m in g.modules:
                    print m.channels
                    print m.bindings
                    print m.operation

