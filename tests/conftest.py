import pytest
from generate import get_database


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
