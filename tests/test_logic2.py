import pytest
import numpy
from bricolage import logic2 as T
from bricolage import operand

@pytest.fixture
def basic_params():
    return T.Parameters()

@pytest.fixture
def p_3x2():
    o = T.Operand
    ops = o.NOT_A_AND_B, o.A_AND_NOT_B, o.NOR, o.AND
    return T.Parameters(seed=1, operands=ops, cis_count=2, reg_channels=5,
                        out_channels=2, cue_channels=3)
@pytest.fixture
def c_3x2(p_3x2):
    world = T.World(p_3x2)
    return T.Constructor(world, p_3x2)

def test_constructor(p_3x2):
    world = T.World(p_3x2)
    const = T.Constructor(world, p_3x2)
    assert set(const.operands) == set(p_3x2.operands)

def test_network_ids(c_3x2):
    for i in range(10):
        n = T.Network(c_3x2)
        assert n.identifier == i

    pop = T.NetworkCollection(c_3x2, 10)
    assert pop[9].identifier == 19

def test_network_construction(c_3x2):
    w = c_3x2.world
    pop = T.NetworkCollection(c_3x2, 1000)
    for n in pop:
        assert len(n.genes) == c_3x2.gene_count
        for g in n.genes:
            assert len(g.modules) == c_3x2.module_count
            assert g.pub in range(*w.pub_range)
            for m in g.modules:
                assert T.Operand(m.op) in c_3x2.operands
                for i in range(m.site_count()):
                    assert w.sub_range[0] <= m.get_site(i) < w.sub_range[1]

def test_referencing(c_3x2):
    original = T.Network(c_3x2)
    nc = T.NetworkCollection(c_3x2, 0)

    # Get the "pointer" value
    id_original = id(original)

    nc.add(original)
    del original

    # These should be the same, as we re-use python references if they exist
    a = nc[0]

    # This should be the same python object that we deleted above, as the
    # reference was maintained internally 
    assert id(a) == id_original

    # Getting it again should work too.
    b = nc[0]
    assert a is b

def test_bad_access(c_3x2):
    nc = T.NetworkCollection(c_3x2, 0)
    with pytest.raises(IndexError):
        a = nc[0]

    nc.add(T.Network(c_3x2))
    a = nc[0]
    b = nc[0]

    assert a is b

# ---------- HERE 
def network_cycle(network, curstate):
    """A Python version of what the C++ cycle does."""
    nextstate = network.factory.create_state()
    for g in network.genes:
        for m in g.modules:
            if operand.calculate(m.op, curstate[m.sub1], curstate[m.sub2]):
                nextstate.set(g.pub)
                break
    return nextstate

def construct_attractor(net, env):
    """Make an attractor by cycling repeatedly"""
    cur = env.copy()
    path = []
    path.append(cur)
    while 1:
        cur = network_cycle(net, cur)
        cur.merge(env)
        for i, prev in enumerate(path):
            if cur == prev:
                return tuple(path[i:])
        path.append(cur)

def test_attractors(p_3x2, basic_params):
    f = T.Factory(basic_params)
    # f = T.Factory(p_3x2)
    nc = f.create_collection(1)
    for net in nc:
        pattractors = [construct_attractor(net, env) for env in f.environments]
        # Make sure the attractors created in C++ are the same
        assert tuple(pattractors) == net.attractors

        # It should be readonly
        with pytest.raises(ValueError):
            net.rates[0, 0] = 10.

def test_collection(p_3x2):
    f = T.Factory(p_3x2)

    nets = f.create_collection(1000)

    old_nets = [_ for _ in nets]
    mutated = nets.mutate(.01)
    #
    # nmutations = 0
    # for i, n in enumerate(nets):
    #     if i in mutated:
    #         nmutations += 1
    #         assert n is not old_nets[i]
    #     else:
    #         assert n is old_nets[i]
    #
    # assert nmutations > 0

def test_mutation(p_3x2):
    f = T.Factory(p_3x2)
    orig = f.create_network()
    mute = orig.mutated(1)
    print orig, mute
    for g1, g2 in zip(orig.genes, mute.genes):
        for m1, m2 in zip(g1.modules, g2.modules):
            if m1 != m2:
                print g1, m1
                print g2, m2

    for i, (a1, a2) in enumerate(zip(orig.attractors, mute.attractors)):
        if a1 != a2:
            print [str(x) for x in a1]
            print [str(x) for x in a2]

def test_rates(p_3x2):
    f = T.Factory(p_3x2)
    orig = f.create_network()
    rates = orig.rates

    # It should be readonly
    with pytest.raises(ValueError):
        rates[0, 0] = 10.

    # rates only include the output channels
    outc = p_3x2.out_channels
    nc = f.create_collection(100)
    for net in nc:
        for attr, rate in zip(net.attractors, net.rates):
            test_rate = numpy.zeros(outc)
            for state in attr: 
                test_rate += state.as_array()[-outc:]
            test_rate /= len(attr)
            assert (test_rate == rate).all()

@pytest.fixture
def target_3x2():
    """Return a function for initialising a target that has 3 inputs and 2
    outputs"""
    def make_target(x):
        a, b, c = x
        f1 = .5 if a and b or not c else 1.0
        f2 = 1 if ((a or c) and not (a and b)) and b else 0
        return f1, f2
    return make_target

@pytest.fixture
def target_3x3():
    """Return a function for initialising a target that has 3 inputs and 2
    outputs"""
    def make_target(x):
        a, b, c = x
        res = (a and b) or (b and c)
        if res:
            return 0, 1, .5
        return 1, .5, 0
    return make_target

def test_targets(p_3x2, target_3x2):
    f = T.Factory(p_3x2)
    targ = T.Target(f, target_3x2)
    nc = f.create_collection(100)
    for net in nc:
        diffs = abs(targ.as_array() - net.rates)
        scores = 1.0 - diffs

        # These should be the same (approx)
        assert abs(scores.mean() - targ.assess(net)) < 1e-12

        # summed = scores.sum(axis=1) * ch.fitness_contribution
        # summed = scores.sum(axis=1)

def test_random_engine(p_3x2, target_3x2):
    f = T.Factory(p_3x2)
    f.seed_random_engine(1)
    first_time = [f.get_random_double(0, 1) for _ in range(20)]
    first_time += [f.get_random_int(0, 100) for _ in range(20)]
    
    f.seed_random_engine(1)
    second_time = [f.get_random_double(0, 1) for _ in range(20)]
    second_time += [f.get_random_int(0, 100) for _ in range(20)]

    assert first_time == second_time
        
def test_selection(p_3x2, target_3x2):
    # TODO: Finish testing selection by python
    factory = T.Factory(p_3x2)
    target = T.Target(factory, target_3x2)
    population = factory.create_collection(1000)

    # What does c++ do?
    factory.seed_random_engine(1)
    c_selection_indexes = population.selection_indexes(target)

    # Now, do it in python
    factory.seed_random_engine(1)
    p_selection_indexes = []
    cum_scores = []
    cum_score = 0.0
    for net in population:
        cum_scores.append(cum_score)
        cum_score += target.assess(net)

    indexes = []

@pytest.fixture
def params3x3():
    return T.Parameters(
        seed=1,
        operands = [
            T.Operand.NOT_A_AND_B,
            T.Operand.A_AND_NOT_B,
            T.Operand.NOR,
            T.Operand.AND,
        ],
        cis_count=2,
        reg_channels=3,
        out_channels=3,
        cue_channels=3,
    )


def not_test_analysis(params3x3, target_3x3):
    factory = T.Factory(params3x3)
    target = T.Target(factory, target_3x3)
    net = factory.create_network()
    ana = T.NetworkAnalysis(net)
    edges = ana.get_edges()
    ko = ana.get_knockouts()
    # TODO: Extend this to a large selection 
    for e in ko:
        assert e in edges

def nottest_play(target_3x3):
    p = T.Parameters(
        seed=4,
        operands = [
            T.Operand.NOT_A_AND_B,
            T.Operand.A_AND_NOT_B,
            T.Operand.AND,
            T.Operand.NOR,
            T.Operand.A_AND_NOT_B,
            T.Operand.NOT_A_AND_B,
        ],
        cis_count=1,
        reg_channels=8,
        out_channels=3,
        cue_channels=3,
    )

    factory = T.Factory(p)
    target = T.Target(factory, target_3x3)
    pop = factory.create_collection(1000)
    while 1:
        pop.select(target)
        maxf = max([n.fitness for n in pop])
        if maxf == 1.0:
            break
        pop.mutate()

    winners = []
    for net in pop:
        if net.fitness == 1.0:
            winners.append(net)

    for net in winners:
        for g in net.genes:
            print g
            for m in g.modules:
                print '   ', m
        ana = T.NetworkAnalysis(net)
        for k in ana.knockouts:
            print k

        print '--edges'
        for e in ana.get_edges():
            print e

        break

        print '--'

