import pytest
import numpy
from organismal import pubsub2 as T
from organismal import operand

@pytest.fixture
def basic_params():
    return T.Parameters()

@pytest.fixture
def params3x2():
    return T.Parameters(
        seed=4,
        operands = [
            T.Operand.NOT_A_AND_B,
            T.Operand.A_AND_NOT_B,
            T.Operand.NOR,
            T.Operand.AND,
        ],
        cis_count=2,
        reg_channels=5,
        out_channels=2,
        cue_channels=3,
    )

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

def test_factory(params3x2):
    f = T.Factory(params3x2)

    # Make sure we can stuff from the Factory
    assert params3x2 == f.params
    e = f.environments
    for i in e:
        print repr(i)

def test_network_ids(params3x2):
    f = T.Factory(params3x2)
    for i in range(10):
        n = f.create_network()
        assert n.identifier == i

    pop = f.create_collection(10)
    assert pop[9].identifier == 19

def test_network_construction(params3x2):
    p = params3x2
    f = T.Factory(p)
    pop = f.create_collection(1000)
    for n in pop:
        assert len(n.genes) == p.gene_count
        for g in n.genes:
            print g.modules
            assert len(g.modules) == p.cis_count
            assert g.pub in p.pub_signals
            for c in g.modules:
                assert T.Operand(c.op) in p.operands
                for i in range(c.site_count()):
                    assert p.sub_range[0] <= c.get_site(i) < p.sub_range[1]

def test_referencing(basic_params):
    f = T.Factory(basic_params)
    nc = T.NetworkCollection(f)
    net = f.create_network()

    nc.add(net)
    del net

    # These should be the same, as we re-use python references if they exist
    a = nc[0]
    b = nc[0]

    assert a is b
    assert a.factory is b.factory
    assert a.factory.params is b.factory.params

def test_bad_access(params3x2):
    f = T.Factory(params3x2)
    nc = T.NetworkCollection(f)
    with pytest.raises(IndexError):
        a = nc[0]

    nc.add(f.create_network())
    a = nc[0]
    b = nc[0]

    assert a is b

def test_channelstate(params3x2):
    f = T.Factory(params3x2)
    e2 = f.environments[-1]
    e2_again = f.environments[-1]

    # We should get the same channels states out.
    assert e2 == e2_again
    assert e2 is e2_again

    # When we copy, they should be the same, but not identical.
    copy_e2 = e2.copy()
    assert e2 == copy_e2
    assert e2 is not copy_e2

    # Modify the state -- testing still work
    copy_e2.flip(0)
    assert e2 != copy_e2
    copy_e2.flip(0)
    assert e2 == copy_e2
    
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

def test_attractors(params3x2, basic_params):
    f = T.Factory(basic_params)
    # f = T.Factory(params3x2)
    nc = f.create_collection(1)
    for net in nc:
        pattractors = [construct_attractor(net, env) for env in f.environments]
        # Make sure the attractors created in C++ are the same
        assert tuple(pattractors) == net.attractors

        # It should be readonly
        with pytest.raises(ValueError):
            net.rates[0, 0] = 10.

def test_collection(params3x2):
    f = T.Factory(params3x2)

    nets = f.create_collection(1000)

    old_nets = [_ for _ in nets]
    mutated = nets.mutate()

    nmutations = 0
    for i, n in enumerate(nets):
        if i in mutated:
            nmutations += 1
            assert n is not old_nets[i]
        else:
            assert n is old_nets[i]

    assert nmutations > 0

def test_mutation(params3x2):
    f = T.Factory(params3x2)
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

def test_rates(params3x2):
    f = T.Factory(params3x2)
    orig = f.create_network()
    rates = orig.rates

    # It should be readonly
    with pytest.raises(ValueError):
        rates[0, 0] = 10.

    # rates only include the output channels
    outc = params3x2.out_channels
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

def test_targets(params3x2, target_3x2):
    f = T.Factory(params3x2)
    targ = T.Target(f, target_3x2)
    nc = f.create_collection(100)
    for net in nc:
        diffs = abs(targ.as_array() - net.rates)
        scores = 1.0 - diffs

        # These should be the same (approx)
        assert abs(scores.mean() - targ.assess(net)) < 1e-12

        # summed = scores.sum(axis=1) * ch.fitness_contribution
        # summed = scores.sum(axis=1)

def test_random_engine(params3x2, target_3x2):
    f = T.Factory(params3x2)
    f.seed_random_engine(1)
    first_time = [f.get_random_double(0, 1) for _ in range(20)]
    first_time += [f.get_random_int(0, 100) for _ in range(20)]
    
    f.seed_random_engine(1)
    second_time = [f.get_random_double(0, 1) for _ in range(20)]
    second_time += [f.get_random_int(0, 100) for _ in range(20)]

    assert first_time == second_time
        
def test_selection(params3x2, target_3x2):
    # TODO: Test selection
    factory = T.Factory(params3x2)
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

def test_play(target_3x3):
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
    for i in range(1000):
        pop.select(target)
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
        break

        print '--'

