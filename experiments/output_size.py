import sys

# Make bricolage available...
from pathlib import Path
sys.path.append(str(Path(__file__).absolute().parent.parent))

from bricolage.experiment import Experiment, Treatment, add_argument
from bricolage.threshold3 import Parameters
from boolean_functions import get_functions

USE_FUNCTION = 5

def target(size):
    desc, func = get_functions()[USE_FUNCTION]
    def new_f(a, b, c):
        if func(a, b, c):
            return [0, 1, 1, 0][:size]
        return [1, 0, 0, 1][:size]
    return new_f


params = [Parameters(
    cis_count=2, 
    reg_channels=8,
    out_channels=N, 
    cue_channels=3 ,
    population_size=1000,
    mutation_rate=.001) for N in range(1, 5)]


class MyTreatment(Treatment):
    def run_replicate(self, replicate, lineage):
        if len(lineage.targets) == 0:
            lineage.add_target(target(lineage.params.out_channels))
        size = lineage.params.out_channels
        lineage.extra.flow = [0, 1, 1, 0][:size]
        while lineage.generation < 50000:
            if lineage.generation % 100 == 0:
                print replicate.treatment.seq, replicate.seq, lineage.generation,\
                    lineage.population.worst_and_best()
            lineage.next_generation()


treats = [MyTreatment("output_{}".format(N+1), params[N], 20) for N in range(4)]
the_exp = Experiment('/Users/Brett/Desktop', 'output_size', seed=5).add_all(*treats)


if __name__ == '__main__':
    the_exp.run_from_commandline()

