import cPickle as pickle
import pytest
import numpy
from bricolage import logic2 as T
from bricolage import operand
from bricolage.core import Channels
from bricolage.core_ext import Collection

@pytest.fixture
def c_one_unit():
    o = T.Operand
    ops = o.NOT_A_AND_B, o.A_AND_NOT_B, o.NOR, o.AND
    p = T.Parameters(seed=1, operands=ops, cis_count=1, reg_channels=0,
                        out_channels=1, cue_channels=2)
    w = T.World(p)
    return T.Factory(w)

@pytest.fixture
def p_3x2():
    o = T.Operand
    ops = o.NOT_A_AND_B, o.A_AND_NOT_B, o.NOR, o.AND
    return T.Parameters(seed=1, operands=ops, cis_count=2, reg_channels=5,
                        out_channels=2, cue_channels=3)


@pytest.fixture
def c_3x2(p_3x2):
    world = T.World(p_3x2)
    return T.Factory(world)

def test_factory(p_3x2):
    world = T.World(p_3x2)
    const = T.Factory(world)
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
    const = T.Factory(world)
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
    nextstate = Channels(network.factory.world)
    for g in network.genes:
        for m in g.modules:
            a, b = m.channels
            if operand.calculate(m.op, curstate.test(a), curstate.test(b)):
                nextstate.set(g.pub)
                break
    return nextstate

def construct_attractor(net, env):
    """Make an attractor by cycling repeatedly"""
    cur = env.copy()
    cur.set(1)
    path = []
    path.append(cur)
    while 1:
        cur = network_cycle(net, cur)
        cur.set(1)
        # cur.merge(env)
        for i, prev in enumerate(path):
            if cur == prev:
                return path[i:]
        path.append(cur)

def test_attractors(c_3x2):
    nc = T.Population(c_3x2, 100)
    for net in nc:
        pattractors = [construct_attractor(net, env) 
                       for env in c_3x2.world.environments]
        # Make sure the attractors created in C++ are the same
        assert pattractors == net.attractors

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
        print op
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
    targ = T.DefaultTarget(c_3x2.world, target_3x2)
    nc = T.Population(c_3x2, 1000)
    for net in nc:
        diffs = abs(targ.as_array() - net.rates)
        scores = 1.0 - diffs

        # These should be the same (approx)
        assert abs(scores.mean() - targ.assess(net)) < 1e-12
        # summed = scores.sum(axis=1) * ch.fitness_contribution summed =
        # scores.sum(axis=1)

def test_network_dup(c_3x2):
    n = c_3x2.create_network() #T.Network(c_3x2)
    n.duplicate(1)
    pubs = [g.pub for g in n.genes]
    assert len(set(pubs)) < len(pubs)
    
def test_pop_dup(c_3x2):
    psize = 1000
    pop = T.Population(c_3x2, psize)
    assert pop.size == psize
    assert pop.mutated == []

    nm = pop.mutate(.01, .1, 0)
    assert len(pop.mutated) == nm
    dupped = 0
    for i in pop.mutated:
        n = pop[i]
        pubs = [g.pub for g in n.genes]
        if len(set(pubs)) < len(pubs):
            dupped += 1

    # There should be lots
    assert dupped > 200

def test_network_trans_mutation(c_3x2):
    n = c_3x2.create_network() #T.Network(c_3x2)
    n.mutate(0, 1)
    pubs = [g.pub for g in n.genes]
    assert len(set(pubs)) < len(pubs)

def test_reg_gene_size_big():
    params = T.Parameters(seed=4, cis_count=2, reg_channels=5, out_channels=2,
                        cue_channels=3, reg_gene_count=6)
    world = T.World(params)
    const = T.Factory(world)
    n1 = const.create_network()
    print world.reg_range
    print world.out_range
    print [g.pub for g in n1.genes]

def test_reg_gene_size_small():
    params = T.Parameters(seed=4, cis_count=2, reg_channels=20, out_channels=2,
                        cue_channels=3, reg_gene_count=5)
    world = T.World(params)
    const = T.Factory(world)
    n1 = const.create_network()
    print world.reg_range
    print world.out_range
    print [g.pub for g in n1.genes]
    n1.mutate(0, 1)
    # n1.calc_attractors()
    print [g.pub for g in n1.genes]
    n1.mutate(0, 1)
    # n1.calc_attractors()

def test_trans_pop():
    params = T.Parameters(seed=4, cis_count=2, reg_channels=30, out_channels=2,
                        cue_channels=3, reg_gene_count=6)
    world = T.World(params)
    fact = T.Factory(world)
    pop = T.Population(fact, 5)
    for n in pop:
        print [g.pub for g in n.genes]
    pop.mutate(0, .3, 0)
    for n in pop:
        print [g.pub for g in n.genes]
    # pop.mutate(0, .01, 0)
    # pop.mutate(0, .01, 0)
    # n1 = const.create_network()
    # print world.reg_range
    # print world.out_range
    # print [g.pub for g in n1.genes]
    #

def test_diff():
    params = T.Parameters(seed=1, cis_count=3, reg_channels=4, out_channels=2,
                        cue_channels=2)
    world = T.World(params)
    fact = T.Factory(world)
    n = fact.create_network()
    coll = Collection(fact)
    coll.add(n)
    np = fact.to_numpy(coll)
    coll2 = Collection(fact)
    coll2.fill_with_mutations(n, numpy.asarray(range(20)))

    q = fact.to_numpy(coll2)

    for i in range(20):
        # print (q['sub'][i] != np['sub'][0]).sum(axis=2)
        # print (q['op'][i] != np['op'][0])
        # print q['op'][i] != np['op'][0]
        # print (q['sub'][i] != np['sub'][0]).sum()
        # print (q['op'][i] != np['op'][0]).sum()

        # print T.count_diff_networks(n, coll2[i])
        print T.modules_changed(n, coll2[i])
