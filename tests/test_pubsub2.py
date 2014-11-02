# import sys
# from pathlib import Path
# sys.path.append(str(Path('.')))
import pytest
from organismal import pubsub2 as T

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
    print a
    b = nc[-1]

    assert a is b


def test_network():
    p = T.Parameters(
        seed=3,
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
    envs = f.environments
    print 
    for g in net.genes:
        print g,
        for m in g.modules:
            print m,
        print
    for e in envs:
        cur = e.__copy__()
        for i in range(4):
            print cur
            net.cycle(cur)
            cur.merge(e)
        print '--'








