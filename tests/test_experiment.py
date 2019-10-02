from bricolage.experiment import Experiment, Treatment
from pathlib import Path
from bricolage.threshold3 import Parameters, Population, DefaultTarget


def target1(a, b):
    if a and b:
        return 1.0
    return 0.0


class MyTreatment(Treatment):
    def run_replicate(self, replicate, lineage):
        if len(lineage.targets) == 0:
            lineage.add_target(DefaultTarget(lineage.world, target1))
        while lineage.generation < 100:
            lineage.next_generation()


_params = Parameters(
    cis_count=2,
    reg_channels=1,
    out_channels=1,
    cue_channels=2,
    population_size=100,
    mutation_rate=0.001,
)


def test_exp1(tmpdir):
    tmpdir = Path(str(tmpdir))
    treats = [MyTreatment("bob", _params, 10)]
    e = Experiment(tmpdir, treats, seed=1)
    e.run()


class CloningTreatment(Treatment):
    def run_replicate(self, replicate, lineage):
        if len(lineage.targets) == 0:
            lineage.add_target(DefaultTarget(lineage.world, target1))
        while lineage.generation < 100:
            lineage.next_generation()

    def make_initial_population(self, replicate, factory, size):
        p = Population(factory)
        n = factory.create_network()

        # Note: we need to do this because the identifiers of the networks
        # must match their index in the database (and we just created one
        # above).
        factory.world.next_network_id = 0
        p.fill(n, size)
        return p


def test_exp2(tmpdir):
    tmpdir = Path(str(tmpdir))
    treats = [CloningTreatment("bob", _params, 10)]
    e = Experiment(tmpdir, treats, seed=1)
    e.run()
