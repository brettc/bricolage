import pytest
from generate import get_database
from bricolage.experiment import Experiment, Treatment
from bricolage.threshold3 import Parameters
import pathlib


@pytest.yield_fixture
def bowtie_database():
    db = get_database('bowtie', readonly=True)
    yield db
    db.close()


@pytest.fixture
def bowtie_env_categories(bowtie_database):
    """Categorise the targets"""
    targ = bowtie_database.targets[0]
    return targ.calc_categories()


@pytest.fixture
def bowtie_network(bowtie_database):
    return bowtie_database.population.get_best()[0]


@pytest.yield_fixture
def xor_database():
    db = get_database('xor', readonly=True)
    yield db
    db.close()


@pytest.fixture
def xor_network(xor_database):
    return xor_database.population.get_best()[0]


@pytest.yield_fixture
def perturb_database():
    db = get_database('perturb', readonly=True)
    yield db
    db.close()

@pytest.yield_fixture
def three_database():
    db = get_database('three', readonly=True)
    yield db
    db.close()

class TestTreatment(Treatment):
    def __init__(self, name, params, count, target):
        super(TestTreatment, self).__init__(name, params, count)
        self.target = target

    def run_replicate(self, replicate, lineage):
        if len(lineage.targets) == 0:
            lineage.add_target(self.target)
        while lineage.generation < 100:
            lineage.next_generation()


# A Makeshift class until we clean up the experiment class
class Args(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


@pytest.fixture(scope='session')
def simple_experiment(tmpdir_factory):
    tmpdir = tmpdir_factory.mktemp('experiments').strpath

    p = Parameters(
        cis_count=2,
        reg_channels=4,
        out_channels=2,
        cue_channels=2,
        population_size=100,
        mutation_rate=.001,
    )

    def target_and_or(a, b):
        return [a and b, a or b]

    def target_or_not(a, b):
        return [a or b, a and not b]

    tmp_path = pathlib.Path(tmpdir)

    treats = [
        TestTreatment('and_', p, 3, target_and_or),
        TestTreatment('or_', p, 3, target_or_not),
    ]
    e = Experiment(tmp_path, treats, seed=1)
    e.run(Args(overwrite=True, dry=False))

    return e
