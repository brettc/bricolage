import pytest
from bricolage.neighbourhood import NeighbourhoodSample, Collection
from bricolage.threshold3 import World, Parameters, Constructor, MutateType
import numpy as np
from bricolage.analysis import get_population_neighbourhood_fitness
from generate import get_database

np.set_printoptions(linewidth=150)

@pytest.yield_fixture
def bowtie_database():
    db = get_database('bowtie')
    yield db
    db.close()

def get_binding_values(construct, net):
    # Create a collection just so we can get the mutations in numpy format
    base = Collection(construct)
    base.add(net)
    np = construct.to_numpy(base)
    return np['binding'].ravel()

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
        
def test_population(bowtie_database):
    target = bowtie_database.targets[0]
    popul = bowtie_database.population
    mean = get_population_neighbourhood_fitness(popul, target, sample_per_network=100)
    fits = np.zeros(popul.size)
    popul.get_fitnesses(fits)
    print fits.mean(), mean




