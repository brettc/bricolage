from bricolage import threshold3
from bricolage.core import InputType
# from bricolage.lineage import SnapshotLineage
from bricolage.targets_ext import MultiTarget
from bricolage.core import Channels

def test_create():
    p = threshold3.Parameters(
        seed=1, cis_count=2, reg_channels=4, out_channels=1, cue_channels=2,
        population_size=1000, mutation_rate=.002, input_type=InputType.PULSE)
        
    w = threshold3.World(p)
    factory = threshold3.Factory(w)
    # w 
    # m = MultiTarget(w)
    # print m
    #
    c = Channels(w)
    # c.set(1)
    # print c
    n = factory.create_network()
    print n
    print n.stabilise(c)


def recurse_attr(start, inputs, net):
    start.merge(inputs[0])
    more = inputs[1:]
    attr, trans, rates = net.stabilise(start)
    if not more:
        return rates

    count_rates = len(rates)
    mult = 1.0/float(count_rates)
    averaged_rates = [0.0 for _ in range(count_rates)]
    for a in attr:
        rates = recurse_attr(a, more, net)
        for i, r in enumerate(rates):
            averaged_rates[i] = r * mult

    return averaged_rates

def test_network(bowtie_network):
    net = bowtie_network
    print net
    f = net.factory
    w = f.world
    pulses = [Channels(w) for _ in range(2)]
    pulses[1].set(2)
    pulses[0].set(3)
    pulses[1].set(4)

    # for i, p in enumerate(pulses):
    #     p.set(2 + i)

    for e in w.environments:
        starts = [e.copy().filter(p) for p in pulses]
        rate = recurse_attr(Channels(w), starts, net)
        print rate

        # for p in pulses:
        #     c = e.copy()
        #     c.filter(p)
        #     attr, trans, rates = net.stabilise(c)
        #     print attr
    


