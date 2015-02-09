import pytest
import pathlib
import cPickle as pickle
from bricolage import threshold3, graph
from bricolage.core import InterventionState
from bricolage.core_ext import NetworkLockedError

@pytest.fixture
def net_1():
    here = pathlib.Path(__file__).parent
    with open(str(here / 'data' / 'problem_1.network')) as f:
        net = pickle.load(f)
    return net

# def network_cycle(network, curstate):
#     """A Python version of what the C++ cycle does."""
#     nextstate = network.constructor.world.create_state()
#     for g in network.genes:
#         for m in g.modules:
#             a, b = m.channels
#             if operand.calculate(m.op, curstate[a], curstate[b]):
#                 nextstate.set(g.pub)
#                 break
#     return nextstate
#
# def construct_attractor(net, env):
#     """Make an attractor by cycling repeatedly"""
#     cur = env.copy()
#     path = []
#     path.append(cur)
#     while 1:
#         cur = network_cycle(net, cur)
#         cur.merge(env)
#         for i, prev in enumerate(path):
#             if cur == prev:
#                 return tuple(path[i:])
#         path.append(cur)

def test_1(net_1):
    print "BEGIN"
    # ana_1 = threshold3.NetworkAnalysis(net_1)
    # ana_1
    assert net_1.locked
    with pytest.raises(NetworkLockedError):
        net_1.genes[0].intervene = InterventionState.INTERVENE_OFF

    ana = threshold3.NetworkAnalysis(net_1)
    # print ana.get_active_edges()
    gph = graph.FullGraph(ana, knockouts=True)
    # gph = graph.GeneSignalGraph(ana, knockouts=False)
    # gph = graph.GeneGraph(ana, knockouts=False)
    dot = graph.DotMaker(gph)
    dot.save_picture('problem_1.png')

    # Get an unlocked copy
    n = net_1.unlocked_copy()
    # edges = []

    # This is what analysis does...
    print net_1.attractors
    # for i, g in enumerate(n.genes):
    for i in range(n.gene_count-1, -1, -1):
        g  = n.genes[i]
        # print 'testing', g
        g.intervene = InterventionState.INTERVENE_OFF
        # if g.sequence == 0:
        #     print net_1.rates 
        #     print n.rates
            # print x
            # print x.all()

        if (n.rates == net_1.rates).all():
            print i, "doesn't matter"
            continue
        # Reset
        g.intervene = InterventionState.INTERVENE_NONE

        # for j, m in enumerate(g.modules):
        #     # print 'testing', g, m
        #     m.intervene = InterventionState.INTERVENE_OFF
        #     # print 'testing ', g, m
        #     if (n.rates == net_1.rates).all():
        #         # print i, j, "doesn't matter"
        #         continue
        #     # if g.sequence == 3:
        #     #     print n.rates == net_1.rates
        #     m.intervene = InterventionState.INTERVENE_NONE


    # Round 2
    print 'round 2'
    for i, g in enumerate(n.genes):
        # print 'testing', g
        g.intervene = InterventionState.INTERVENE_OFF
        # if g.sequence == 0:
        #     print net_1.rates 
        #     print n.rates
            # print x
            # print x.all()

        if (n.rates == net_1.rates).all():
            print i, "doesn't matter"
            continue
        # Reset
        g.intervene = InterventionState.INTERVENE_NONE

def test_2(net_1):
    print 
    n = net_1.unlocked_copy()
    for i in 2, 5, 7:
        n.genes[i].intervene = InterventionState.INTERVENE_OFF

    # n.genes[0].intervene = InterventionState.INTERVENE_OFF
    # print net_1.rates
    # print n.rates
    print n.rates == net_1.rates

    print net_1.attractors[1]
    print n.attractors[1]
    # env = n.constructor.world.environments[1].copy()
    # for i in range(5):
    #     n.cycle_with_intervention(env)
    #     print env



