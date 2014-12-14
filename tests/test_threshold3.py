
# import pytest
from organismal.core_ext import Target
from organismal.threshold3 import Parameters, Factory

def not_fun(x):
    return [not x[0]]

def test_binding():
    p = Parameters(cis_count=1, cue_channels=1, reg_channels=0, out_channels=1)
    f = Factory(p)
    target = Target(f, not_fun)
    pop = f.create_collection(1000)
    while 1:
        pop.select(target)
        fits = [n.fitness for n in pop]
        mfit = max(fits)
        print mfit
        if mfit == 1.0:
            break
        pop.mutate()
        break

    print f.environments
    for n in pop:
        for g in n.genes:
            for m in g.modules:
                print m.channels
                print m.bindings
                print m.operation

