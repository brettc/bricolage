import pytest
from bricolage.neighbourhood import NetworkNeighbourhood, Collection, \
    PopulationNeighbourhood
from bricolage.threshold3 import World, Parameters, Factory, MutateType
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


def test_neighbourhood():
    # Note we use the PROGRESSIVE mutation as this guarantees a change!
    params = Parameters(seed=4, cis_count=2, reg_channels=5, out_channels=2,
                        cue_channels=3, mutate_type=MutateType.PROGRESSIVE)
    world = World(params)
    const = Factory(world)
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
    print popul.fitnesses.mean(), sum(popul.fitnesses == 1.0)

    targ.assess_collection(nayb.neighbours)
    print nayb.neighbours.fitnesses.mean(), sum(
        nayb.neighbours.fitnesses == 1.0)
    # print fits.mean(), mean
