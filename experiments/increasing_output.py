import sys

# Make bricolage available...
from pathlib import Path
sys.path.append(str(Path(__file__).absolute().parent.parent))

from bisect import bisect

from bricolage.experiment import Experiment, Treatment
from bricolage.threshold3 import Parameters, MutateType
from boolean_functions import get_functions

USE_FUNCTION = 5

def target():
    desc, func = get_functions()[USE_FUNCTION]
    def new_f(a, b, c):
        if func(a, b, c):
            return [0, 1, 1, 0]
        return [1, 0, 0, 1]
    return new_f


params = Parameters(
    cis_count=3, 
    reg_channels=8,
    out_channels=4, 
    cue_channels=3 ,
    population_size=1000,
    mutation_rate=.001)

params2 = Parameters(
    cis_count=3, 
    reg_channels=8,
    out_channels=4, 
    cue_channels=3 ,
    population_size=1000,
    mutation_rate=.005,
    mutate_type=MutateType.PROGRESSIVE,
    add_zeros=5)


class MyTreatment(Treatment):
    def run_replicate(self, replicate, lineage):
        if len(lineage.targets) == 0:
            w = [0, 0, 0, 0]
            for i in range(4):
                w[i] = 1
                lineage.add_target(target(), weighting=w)

        # The same as the first target
        lineage.extra.flow = [0, 1, 1, 0]
        lineage.set_target(0)

        while lineage.generation < 50000:
            t = bisect([3000, 10000, 20000], lineage.generation)
            if t != lineage.target_index:
                lineage.set_target(t)
                
            # Complexify the targets
            if lineage.generation % 5000 == 0:
                replicate.draw_winners(lineage)

            if lineage.generation % 100 == 0:
                print replicate.treatment.seq, replicate.seq, lineage.generation,\
                    lineage.population.worst_and_best(), 'target:', lineage.target_index
            lineage.next_generation()

        replicate.draw_winners(lineage)


# t = MyTreatment("x", params, 4)
t = MyTreatment("y", params2, 20)
the_exp = Experiment(
    '/Users/Brett/Desktop', 
    'increasing_size', 
    analysis_path='/Users/Brett/Dropbox/SimulationOutput/',
    seed=1
)
the_exp.add(t)


if __name__ == '__main__':
    the_exp.run_from_commandline()

