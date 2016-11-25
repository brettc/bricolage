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
    mutation_rate=.001,
)


def test_exp1(tmpdir):
    tmpdir = Path(str(tmpdir))
    # pth = Path('.')
    treats = [MyTreatment('bob', _params, 10)]
    e = Experiment(tmpdir, treats, seed=1, analysis_path="/Users/brett/Desktop")
    e.run()

    # with e.treatments[0].replicates[5].get_lineage(True) as l:
    #     print l.generation
    #     print l.population.worst_and_best()


class CloningTreatment(Treatment):

    def run_replicate(self, replicate, lineage):
        if len(lineage.targets) == 0:
            lineage.add_target(DefaultTarget(lineage.world, target1))
        while lineage.generation < 100:
            lineage.next_generation()

    def make_initial_population(self, replicate, factory, size):
        p = Population(factory)
        n = factory.create_network()
        p.fill(n, size)
        return p


def test_exp2(tmpdir):
    tmpdir = Path(str(tmpdir))
    treats = [CloningTreatment('bob', _params, 10)]
    e = Experiment(tmpdir, treats, seed=1)
    e.run()

# def test_treatment(tmpdir, p_3x2):
#     tmpdir = pathlib.Path(str(tmpdir))
#     name = 'treat'
#     path = tmpdir / name
#     apath = tmpdir / (name + '_analysis')
#     treat = L.Treatment(path, p_3x2, analysis_path=apath
#     max_gen = 100
#
#     for rep in treat.iter_replicates():
#         with rep.get_lineage() as lin:
#             lin.add_target(target1)
#             while lin.generation < max_gen:
#                 lin.next_generation()
#
#         assert rep.analysis_path is not rep.path
#         assert rep.analysis_path.exists()
#
#     with treat.replicates[5].get_lineage() as l5:
#         assert l5.population
#
#     # Kill this
#     del treat
#
#     treat = L.Treatment(path, p_3x2)
#     for rep in treat.iter_replicates():
#         with rep.get_lineage() as lin:
#             assert lin.generation == max_gen
#
#
