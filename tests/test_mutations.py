from bricolage.neighbourhood import Collection
from bricolage.threshold3 import World, Parameters, Factory, MutateType
import numpy as np
np.set_printoptions(linewidth=150)

def get_binding_values(fact, net):
    # Create a collection just so we can get the mutations in numpy format
    base = Collection(fact)
    base.add(net)
    np = fact.to_numpy(base)
    return np['binding'].ravel()

def get_cis_values(fact, net):
    # Create a collection just so we can get the mutations in numpy format
    base = Collection(fact)
    base.add(net)
    np = fact.to_numpy(base)
    return np['sub'].ravel()

def test_progressive():
    params = Parameters(seed=4, cis_count=2, reg_channels=5, out_channels=2,
                        cue_channels=3, mutate_type=MutateType.PROGRESSIVE)
    world = World(params)
    fact = Factory(world)
    net = fact.create_network()
    binding1 = get_binding_values(fact, net)
    for i in range(1000):
        net.mutate(1)
        binding2 = get_binding_values(fact, net)
        diffs = (binding1 != binding2).sum()
        assert diffs == 1
        binding1 = binding2

def test_variable():
    params = Parameters(seed=4, cis_count=2, reg_channels=5, out_channels=2,
                        cue_channels=3, mutate_type=MutateType.PROGRESSIVE, add_zeros=20)
    world = World(params)
    fact = Factory(world)
    print fact.draw_from_subs
    net = fact.create_network()
    subs = get_cis_values(fact, net)
    print subs


def make_with_zeros(addz):
    params = Parameters(cis_count=2, reg_channels=5, out_channels=2,
                        cue_channels=3, add_zeros=addz)
    world = World(params)
    fact = Factory(world)
    return fact, world, params

def test_add_zeros():
    f, w, p = make_with_zeros(0)
    assert f.draw_from_subs == range(*w.sub_range)

    f, w, p = make_with_zeros(5)
    assert f.draw_from_subs == range(*w.sub_range) + [0] * 5

    env_reg = range(*w.sub_range)[2:]
    assert len(env_reg) == w.reg_channels + w.cue_channels
    assert [2 + i for i in range(w.reg_channels + w.cue_channels)] == env_reg

    new_subs = env_reg + [0] * len(env_reg)
    f, w, p = make_with_zeros(new_subs)
    assert f.draw_from_subs == new_subs
    print new_subs

