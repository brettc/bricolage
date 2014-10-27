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
    p, f = defaults

    # Make sure we can stuff from the Factory
    assert p == f.params

def test_network_construction(defaults):
    p, f = defaults
    n = f.create_network()
    assert n.ready
    assert len(n.genes) == p.gene_count

def test_cismodules(defaults):
    p, f = defaults
    n = f.create_network()
    print 
    g = n[0]
    print g
    print "modules", g.module_count
    cm = g[1]
    print cm.gene
    # print n.export_genes()

def test_referencing(defaults):
    p, f = defaults
    pop = f.create_population()

    # These should be the same
    a = pop.get(0)
    b = pop.get(0)

    assert a is b
    assert a.factory is b.factory
    assert a.factory.params is b.factory.params

def test_5(defaults):
    p, f = defaults
    n = f.create_network()
    for i, g in enumerate(n):
        g.pub = i
    print n.export_genes()

def test_genes(defaults):
    p, f = defaults
    n = f.create_network()
    g = n[1]
    g.pub = 2
    g.sub1 = 1
    g.sub2 = 2
    g.op = 8 + 4 + 2
    print
    p = T.Products(5)
    p.set(4)
    p.set(3)
    print g.active(p)

def test_products():
    p = T.Products(5)
    p.set(1)
    print str(p)



