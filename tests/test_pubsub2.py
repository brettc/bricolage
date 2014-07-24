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

def test_1(defaults):
    p, f = defaults
    n = T.Network(f)
    assert n.ready == False

    n = f.create_network()
    assert n.ready == True

def test_3():
    p = T.Parameters()
    q = T.Factory(p)
    nok = q.create_network()
    x = nok.export_genes()
    print x
    x[0,0] = 10
    nok.import_genes(x)
    x = nok.export_genes()
    print x
    # print q.cparams
    # nok.test()
    # print nok.export_genes()

def test_35():
    p = T.Parameters()
    f = T.Factory(p)
    n = f.create_network()
    n[4].pub = 4
    n[2].pub = 2

    print n.export_genes()

def test_4():
    p = T.Parameters()
    f = T.Factory(p)
    nets = [f.create_network() for _ in range(10)]
    pop = f.create_population()
    for n in nets:
        print n, n.identifier
        pop.add(n)

    n2 = pop.get(2)
    print n2.export_genes()

    # knok = pop.get(0)
    # print knok
    # print knok.export_genes()
    #
    # del newnok
    # del knok
    # # nok.test()
    # # print nok.export_genes()

def test_5(defaults):
    p, f = defaults
    n = f.create_network()
    for i, g in enumerate(n):
        g.pub = i
    print n.export_genes()
