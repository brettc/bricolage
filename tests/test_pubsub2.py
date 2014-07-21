# import sys
# from pathlib import Path
# sys.path.append(str(Path('.')))
from organismal import pubsub2 as T

def test_1():
    n = T.Network(1)
    print n
    print n.params
    print n.ready

def test_2():
    p = T.Parameters()
    q = T.Factory(p)
    nok = q.create_network()
    assert nok.ready == True

    nbad = T.Network(p)
    assert nbad.ready == False

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

