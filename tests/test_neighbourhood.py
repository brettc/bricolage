import pytest
from bricolage.neighbourhood import NetworkNeighbourhood, Collection, \
    PopulationNeighbourhood
from bricolage import threshold3, logic2
import numpy as np
# from bricolage.frames import get_population_neighbourhood_fitness
from generate import get_database

np.set_printoptions(linewidth=150)


@pytest.yield_fixture
def bowtie_database():
    db = get_database('bowtie', readonly=True)
    yield db
    db.close()


def get_binding_values(construct, net):
    # Create a collection just so we can get the mutations in numpy format
    base = Collection(construct)
    base.add(net)
    np = construct.to_numpy(base)
    return np['binding'].ravel()


def test_neighbourhood_thresh():
    # Note we use the PROGRESSIVE mutation as this guarantees a change!
    params = threshold3.Parameters(
        seed=4, cis_count=2, reg_channels=5, out_channels=2,
        cue_channels=3, 
        mutate_type=threshold3.MutateType.PROGRESSIVE)
    world = threshold3.World(params)
    const = threshold3.Factory(world)
    net = const.create_network()

    base_bnd = get_binding_values(const, net)

    # This should just take 1 step neighbours
    nayb = NetworkNeighbourhood(net, 1000)
    npformat = const.to_numpy(nayb.neighbours)
    # Get all of the bindings out
    bnd = npformat['binding']
    for b in bnd:
        # There should be just one difference!
        assert (base_bnd != b.ravel()).sum() == 1


def test_population(bowtie_database):
    targ = bowtie_database.targets[0]
    popul = bowtie_database.population
    nayb = PopulationNeighbourhood(popul, 10, .5)
    print popul.size, nayb.neighbours.size
    # mean = get_population_neighbourhood_fitness(popul, target,
    # sample_per_network=100)
    # fits = np.zeros(popul.size)
    f1 = targ.assess_collection(popul)
    print f1.mean(), sum(f1 == 1.0)

    f2 = targ.assess_collection(nayb.neighbours)
    print f2.mean(), sum(f2 == 1.0)

def test_logic():
    # Note we use the PROGRESSIVE mutation as this guarantees a change!
    params = logic2.Parameters(seed=4, cis_count=2, reg_channels=5, out_channels=2,
                        cue_channels=3)
    world = logic2.World(params)
    const = logic2.Factory(world)
    net = const.create_network()

    # This should just take 1 step neighbours
    nayb = NetworkNeighbourhood(net, 1000)
    # Get all of the bindings out
    nc = 0
    for n in nayb.neighbours:
        # There should be just one difference!
        ch = logic2.modules_changed(net, n)
        if len(ch) == 0:
            nc += 1
    print nc / float(1000)
        # assert (base_bnd != b.ravel()).sum() == 1

