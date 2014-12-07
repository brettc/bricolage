
from organismal.logic2 import Parameters
from organismal.threshold3_ext import Factory
from organismal.core_ext import NetworkAnalysis, Target
from organismal.graph import *

p = Parameters(
    cis_count=1,
    reg_channels=5,
    out_channels=3,
    cue_channels=3,
)

spat = [1, 0, 0], [0, 1, 1]
cpat = [0, 0, .5], [1, .5, 0]
pat = spat

def make_target(x):
    a, b, c = x
    # res = (a and b) or (b and c)
    res = ((a and b) or (not b and not c))
    if res:
        return pat[0]
    return pat[1]

factory = Factory(p)
target = Target(factory, make_target)
factory.seed_random_engine(5)
pop = factory.create_collection(1000)

generation = 0
drift = 0
while 1:
    pop.select(target)
    fits = [n.fitness for n in pop]
    mfit = max(fits)
    if generation % 100 == 0:
        print generation, mfit
    if mfit == 1.0:
        drift += 1 
    else:
        drift = 0
    if drift == 1000:
        break
    pop.mutate()
    generation += 1
    if generation > 3000:
        break

maxf = max([n.fitness for n in pop])
winners = []
for net in pop:
    if net.fitness == maxf:
        winners.append(net)

for i, net in enumerate(winners):
    # for g in net.genes:
    #     print g
    #     for m in g.modules:
    #         print '   ', m
    ana = NetworkAnalysis(net)
    g = FullGraph(ana)
    d = DotMaker(g)
    print 'output', i
    d.draw('./test-{:02d}.png'.format(i))

    if i > 20:
        break
    # ko = ana.get_knockouts()
    # edges = ana.get_edges()
    # x = [e for e in edges]
    # x.sort()
    # for e in x:
    #     print e

    # print '--edges'
    # for e in edges:
    #     if e in ko:
    #         print "NOT!", 
    #     print e

    # grp = FullGraph(ana)
    # dm = DotMaker(grp)
    # dm.draw('./test01.png')
    #
    # grp = FullGraph(ana, knockouts=False)
    # dm = DotMaker(grp)
    # dm.save('./test02.dot')
    # dm.draw('./test02.png')
    #
    # gcg = GCGraph(ana)
    # dm = DotMaker(gcg)
    # dm.draw('./test03.png')
    # gcg = GCGraph(ana, knockouts=False)
    # dm = DotMaker(gcg)
    # dm.draw('./test04.png')
    # break
