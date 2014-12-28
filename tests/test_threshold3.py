from bricolage.threshold3 import *

def xor_func(a, b):
    return (a or b) and not (a and b)

def not_fun(x):
    a, b = x
    return [xor(a, b)]

def test_binding():
    p = Parameters(seed=2, cis_count=2, cue_channels=2, reg_channels=0, out_channels=1)
    w = World(p)
c = Constructor(w, pgtjjck
    target = Target(w, xor_func)
    pop = NetworkCollection(w, 1000)
    while 1:
        pop.select(target)
        fits = [n.fitness for n in pop]
        mfit = max(fits)
        if mfit == 1.0:
            break
        pop.mutate(.05)
        print mfit

    # print w.environments
    for n in pop:
        if n.fitness == 1.0:
            print n
            for g in n.genes:
                for m in g.modules:
                    print m.channels
                    print m.bindings
                    print m.operation

