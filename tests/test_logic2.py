import cPickle as pickle
import pytest
import numpy
from bricolage import logic2 as T
from bricolage import operand

@pytest.fixture
def c_one_unit():
    o = T.Operand
    ops = o.NOT_A_AND_B, o.A_AND_NOT_B, o.NOR, o.AND
    p = T.Parameters(seed=1, operands=ops, cis_count=1, reg_channels=0,
                        out_channels=1, cue_channels=2)
    w = T.World(p)
    return T.Constructor(w)

@pytest.fixture
def p_3x2():
    o = T.Operand
    ops = o.NOT_A_AND_B, o.A_AND_NOT_B, o.NOR, o.AND
    return T.Parameters(seed=1, operands=ops, cis_count=2, reg_channels=5,
                        out_channels=2, cue_channels=3)


@pytest.fixture
def c_3x2(p_3x2):
    world = T.World(p_3x2)
    return T.Constructor(world)

def test_constructor(p_3x2):
    world = T.World(p_3x2)
    const = T.Constructor(world)
    assert set(const.operands) == set(p_3x2.operands)

def test_network_ids(c_3x2):
    for i in range(10):
        n = c_3x2.create_network() #T.Network(c_3x2)
        assert n.identifier == i

    pop = T.Population(c_3x2, 10)
    assert pop[9].identifier == 19

def test_network_construction(c_3x2):
    w = c_3x2.world
    pop = T.Population(c_3x2, 1000)
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
    original = c_3x2.create_network()
    nc = T.Population(c_3x2, 0)

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
    nc = T.Population(c_3x2, 0)
    with pytest.raises(IndexError):
        a = nc[0]

    nc.add(c_3x2.create_network())
    a = nc[0]
    b = nc[0]
    assert a is b

def test_network_pickle():
    params = T.Parameters(seed=4, cis_count=2, reg_channels=5, out_channels=2,
                        cue_channels=3,)
    world = T.World(params)
    const = T.Constructor(world)
    n1 = const.create_network()
    out = pickle.dumps(n1, -1)
    n2 = pickle.loads(out)
    for g1, g2 in zip(n1.genes, n2.genes):
        assert g1.pub == g2.pub
        for m1, m2 in zip(g1.modules, g2.modules):
            assert m1.same_as(m2)

def test_population_mutation(c_3x2):
    psize = 1000
    pop = T.Population(c_3x2, psize)
    assert pop.size == psize
    assert pop.mutated == []

    nm = pop.mutate(.01)
    assert len(pop.mutated) == nm
    print pop.mutated
    for i, n in enumerate(pop):
        if n.parent_identifier != -1:
            assert i in pop.mutated


def network_cycle(network, curstate):
    """A Python version of what the C++ cycle does."""
    nextstate = network.constructor.world.create_state()
    for g in network.genes:
        for m in g.modules:
            a, b = m.channels
            if operand.calculate(m.op, curstate[a], curstate[b]):
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

def test_attractors(c_3x2):
    nc = T.Population(c_3x2, 100)
    for net in nc:
        pattractors = [construct_attractor(net, env) 
                       for env in c_3x2.world.environments]
        # Make sure the attractors created in C++ are the same
        assert tuple(pattractors) == net.attractors

def test_rates(c_3x2):
    network = c_3x2.create_network()
    rates = network.rates

    # It should be readonly
    with pytest.raises(ValueError):
        rates[0, 0] = 10.

    nc = T.Population(c_3x2, 100)
    # rates only include the output channels
    outc = c_3x2.world.out_channels
    for net in nc:
        for attr, rate in zip(net.attractors, net.rates):
            test_rate = numpy.zeros(outc)
            for state in attr: 
                test_rate += state.as_array()[-outc:]
            test_rate /= len(attr)
            assert (test_rate == rate).all()

# def test_mutation(c_3x2):
#     orig = f.create_network()
#     mute = orig.mutated(1)
#     print orig, mute
#     for g1, g2 in zip(orig.genes, mute.genes):
#         for m1, m2 in zip(g1.modules, g2.modules):
#             if m1 != m2:
#                 print g1, m1
#                 print g2, m2
#
#     for i, (a1, a2) in enumerate(zip(orig.attractors, mute.attractors)):
#         if a1 != a2:
#             print [str(x) for x in a1]
#             print [str(x) for x in a2]

def test_cis_manipulation(c_one_unit):
    net = c_one_unit.create_network()
    cis = net.genes[0].modules[0]
    for i in range(16):
        op = T.Operand(i)
        cis.op = op
        for a, b in ((0, 0), (0, 1), (1, 0), (1, 1)):
            assert operand.calculate(op, a, b) == cis.test(a, b)

def test_cis_mutation(c_one_unit):
    net = c_one_unit.create_network()
    cis = net.genes[0].modules[0]

    # Get originals
    oo = cis.op
    oa, ob = cis.channels

    # ROUGH: But everything should change within 20 cycles...
    for i in range(20):
        cis.mutate()
        o = cis.op
        a, b = cis.channels
        if o != oo and oa != a and ob != b:
            break
    else:
        assert False

@pytest.fixture
def target_3x2():
    """Return a function for initialising a target that has 3 inputs and 2
    outputs"""
    def make_target(a, b, c):
        f1 = .5 if a and b or not c else 1.0
        f2 = 1 if ((a or c) and not (a and b)) and b else 0
        return f1, f2
    return make_target

def test_targets(c_3x2, target_3x2):
    targ = T.Target(c_3x2.world, target_3x2)
    nc = T.Population(c_3x2, 1000)
    for net in nc:
        diffs = abs(targ.as_array() - net.rates)
        scores = 1.0 - diffs

        # These should be the same (approx)
        assert abs(scores.mean() - targ.assess(net)) < 1e-12

        # summed = scores.sum(axis=1) * ch.fitness_contribution
        # summed = scores.sum(axis=1)
        
# @pytest.fixture
# def target_3x3():
#     """Return a function for initialising a target that has 3 inputs and 2
#     outputs"""
#     def make_target(a, b, c):
#         res = (a and b) or (b and c)
#         if res:
#             return 0, 1, .5
#         return 1, .5, 0
#     return make_target
#
# def test_selection(p_3x2, target_3x2):
#     # TODO: Finish testing selection by python
#     factory = T.Factory(p_3x2)
#     target = T.Target(factory, target_3x2)
#     population = factory.create_collection(1000)
#
#     # What does c++ do?
#     factory.seed_random_engine(1)
#     c_selection_indexes = population.selection_indexes(target)
#
#     # Now, do it in python
#     factory.seed_random_engine(1)
#     p_selection_indexes = []
#     cum_scores = []
#     cum_score = 0.0
#     for net in population:
#         cum_scores.append(cum_score)
#         cum_score += target.assess(net)
#
#     indexes = []
