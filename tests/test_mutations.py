from bricolage.neighbourhood import NeighbourhoodSample, Collection
from bricolage.threshold3 import World, Parameters, Constructor, MutateType
import numpy as np
np.set_printoptions(linewidth=150)

def get_binding_values(construct, net):
    # Create a collection just so we can get the mutations in numpy format
    base = Collection(construct)
    base.add(net)
    np = construct.to_numpy(base)
    return np['binding'].ravel()

def get_cis_values(construct, net):
    # Create a collection just so we can get the mutations in numpy format
    base = Collection(construct)
    base.add(net)
    np = construct.to_numpy(base)
    return np['sub'].ravel()

def test_progressive():
    params = Parameters(seed=4, cis_count=2, reg_channels=5, out_channels=2,
                        cue_channels=3, mutate_type=MutateType.PROGRESSIVE)
    world = World(params)
    const = Constructor(world)
    net = const.create_network()
    binding1 = get_binding_values(const, net)
    for i in range(1000):
        net.mutate(1)
        binding2 = get_binding_values(const, net)
        diffs = (binding1 != binding2).sum()
        assert diffs == 1
        binding1 = binding2

def test_variable():
    params = Parameters(seed=4, cis_count=2, reg_channels=5, out_channels=2,
                        cue_channels=3, mutate_type=MutateType.PROGRESSIVE, add_zeros=20)
    world = World(params)
    const = Constructor(world)
    print const.draw_from_subs
    net = const.create_network()
    subs = get_cis_values(const, net)
    print subs

def test_neighbourhood():
    # Note we use the PROGRESSIVE mutation as this guarantees a change!
    params = Parameters(seed=4, cis_count=2, reg_channels=5, out_channels=2,
                        cue_channels=3, mutate_type=MutateType.PROGRESSIVE)
    world = World(params)
    const = Constructor(world)
    net = const.create_network()

    base_bnd = get_binding_values(const, net)

    # This should just take 1 step neighbours
    nayb = NeighbourhoodSample(net, 1000)
    npformat = const.to_numpy(nayb.neighbours)
    # Get all of the bindings out
    bnd = npformat['binding']
    for b in bnd:
        # There should be just one difference!
        assert (base_bnd != b.ravel()).sum() == 1
        
        


