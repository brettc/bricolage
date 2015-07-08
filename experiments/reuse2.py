import sys

# Make bricolage available...
from pathlib import Path
sys.path.append(str(Path(__file__).absolute().parent.parent))

from bisect import bisect

from bricolage.experiment import Experiment, Treatment
from bricolage.threshold3 import Parameters, MutateType
from boolean_functions import get_functions

F1 = 5
F2 = 5

def target():
    _, f1 = get_functions()[F1]
    _, f2 = get_functions()[F2]
    def new_f(a, b, c, d, e, f):
        ret = []
        if f1(a, b, c):
            ret.extend([0, 1, 1])
        else:
            ret.extend([1, 0, 0])
        if f2(d, e, f):
            ret.extend([0, 1, 1])
        else:
            ret.extend([1, 0, 0])
        return ret
    return new_f

params = Parameters(
    cis_count=3, 
    reg_channels=16,
    out_channels=6, 
    cue_channels=6 ,
    population_size=1000,
    # mutation_rate=.001)
    mutation_rate=.005,
    mutate_type=MutateType.PROGRESSIVE,
    add_zeros=4)

class MyTreatment(Treatment):
    def run_replicate(self, replicate, lineage):
        if len(lineage.targets) == 0:
            lineage.add_target(target(), weighting=[1,1,1,0,0,0])
            lineage.add_target(target(), weighting=[1,1,1,1,1,1])

        # The same as the first target
        # lineage.extra.flow = [0, 1, 1, 0]
        lineage.set_target(0)

        while lineage.generation < 100000:
            g = lineage.generation
            if g >= 30000:
                lineage.set_target(1)

            # Complexify the targets
            if g > 0 and g % 10000 == 0:
                replicate.draw_winners(lineage)

            if g % 100 == 0:
                print replicate.treatment.seq, replicate.seq, lineage.generation,\
                    lineage.population.worst_and_best(), 'target:', lineage.target_index
            lineage.next_generation()

        replicate.draw_winners(lineage)


# t = MyTreatment("x", params, 4)
t = MyTreatment("double", params, 20)
the_exp = Experiment(
    '/Users/Brett/Desktop', 
    'reuse2', 
    analysis_path='/Users/Brett/Dropbox/SimulationOutput/',
    seed=1
)
the_exp.add(t)


if __name__ == '__main__':
    the_exp.run_from_commandline()

