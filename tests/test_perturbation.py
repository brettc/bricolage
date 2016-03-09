import pytest
from pathlib import Path
from bricolage.threshold3 import Parameters
from bricolage.lineage import FullLineage
from bricolage.core import InputType

@pytest.fixture
def p_3x1():
    return Parameters(
        seed=1,
        cis_count=2,
        reg_channels=4,
        out_channels=1,
        cue_channels=2,
        population_size=1000,
        input_type=InputType.PULSE,
        mutation_rate=.001,
        replicates=10,
    )

def target1(a, b):
    if (a or b) and not (a and b):
        return 1
    return 0

def test1(tmpdir, p_3x1):
    base = Path(str(tmpdir))
    path = base / 'test1.db'
    with FullLineage(path, p_3x1) as a:
        a.add_target(target1)
        for i in range(10000):
            a.next_generation()

    net = a.population.get_best()[0]
    print net.fitness
    print net.attractors
    for i in range(5):
        net.calc_perturbation()
        print net.pert_attractors

def test2(perturb_database):
    db = perturb_database
    net = db.population.get_best()[0]
    print net.fitness
    # print net.attractors[0]
    print net.rates[0], net.rates[1]
    print '---'
    for i in range(150):
        net.calc_perturbation()
    for i in range(5):
        net.calc_perturbation()
        print net.pert_rates[0], net.pert_rates[1]
        # print net.pert_attractors[0]
