# import sys
# from pathlib import Path
# sys.path.append(str(Path('.')))
import pytest
from organismal import pubsub2 as T
from organismal import operand

@pytest.fixture
def defaults():
    p = T.Parameters()
    f = T.Factory(p)
    return p, f

def test_factory(defaults):
    p = T.Parameters(cue_channels=3, reg_channels=5)
    f = T.Factory(p)

    # Make sure we can stuff from the Factory
    assert p == f.params
    e = f.environments
    print
    for i in e:
        print repr(i)

def test_network_ids(defaults):
    p, f = defaults
    for i in range(10):
        n = f.create_network()
        assert n.identifier == i

def test_network_construction():
    p = T.Parameters(
        seed=21, 
        cis_count=4, 
        cue_channels=3, 
        reg_channels=10, 
        out_channels=3
    )

    f = T.Factory(p)
    n = f.create_network()
    assert n.ready
    assert len(n.genes) == p.gene_count
    for g in n.genes:
        assert len(g.modules) == p.cis_count
        assert g.pub in p.pub_signals
        for c in g.modules:
            assert T.Operand(c.op) in p.operands
            assert p.sub_range[0] <= c.sub1 < p.sub_range[1]

def test_referencing(defaults):
    p, f = defaults
    nc = T.NetworkCollection(f)
    net = f.create_network()

    nc.add(net)
    del net

    # These should be the same -- as we re-use python references if they
    # exist
    a = nc[0]
    b = nc[0]

    assert a is b
    assert a.factory is b.factory
    assert a.factory.params is b.factory.params

def test_bad_access(defaults):
    p, f = defaults
    nc = T.NetworkCollection(f)
    with pytest.raises(IndexError):
        a = nc[0]

    nc.add(f.create_network())
    a = nc[0]
    b = nc[0]

    assert a is b

def test_channelstate(defaults):
    p, f = defaults
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
    nextstate = network.factory.create_state()
    for g in network.genes:
        for m in g.modules:
            if operand.calculate(m.op, curstate[m.sub1], curstate[m.sub2]):
                nextstate.set(g.pub)
                break
    return nextstate

def construct_attractor(net, env):
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


def test_attractors():
    # TODO: make a load of these
    p = T.Parameters(
        seed=4,
        operands = [
            T.Operand.NOT_A_AND_B,
            T.Operand.A_AND_NOT_B,
            T.Operand.NOR,
            T.Operand.AND,
            T.Operand.FALSE,

        ],
        cis_count=2,
        reg_channels=5,
        out_channels=2,
        cue_channels=3,
    )
    f = T.Factory(p)
    net = f.create_network()
    pattractors = [construct_attractor(net, env) for env in f.environments]
    assert tuple(pattractors) == net.attractors
    # for patt, catt in zip(pattractors, net.attractors):
    #     assert patt == catt

def test_mutator(defaults):
    p, f = defaults
    net = f.create_network()
    print
    print net.genes[0].modules[0]
    f.mutate_network(net)
    print net.genes[0].modules[0]

def test_collection(defaults):
    p = T.Parameters(
        gene_mutation_rate=.1,
        seed=4,
        operands = [
            T.Operand.NOT_A_AND_B,
            T.Operand.A_AND_NOT_B,
            T.Operand.NOR,
            T.Operand.AND,
            T.Operand.FALSE,

        ],
        cis_count=2,
        reg_channels=5,
        out_channels=3,
        cue_channels=3,
    )
    f = T.Factory(p)
    nc = f.create_collection(10)
    orig = nc[0]
    nc.mutate()
    mute = nc[0]
    print
    assert orig is not mute
    for g1, g2 in zip(orig.genes, mute.genes):
        for m1, m2 in zip(g1.modules, g2.modules):
            if m1 != m2:
                print g1, m1
                print g2, m2

    for i, (a1, a2) in enumerate(zip(orig.attractors, mute.attractors)):
        if a1 != a2:
            print [str(x) for x in a1]
            print [str(x) for x in a2]



