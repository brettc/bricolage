import sys

# Make bricolage available...
from pathlib import Path
sys.path.append(str(Path(__file__).absolute().parent.parent))
from bisect import bisect

from bricolage.experiment import Experiment, Treatment
from bricolage.threshold3 import Parameters, MutateType
from boolean_functions import get_functions

USE_FUNCTION = 10

THIRD = 1.0/3.0

def target():
    desc, func = get_functions()[USE_FUNCTION]
    def new_f(a, b, c):
        if func(a, b, c):
            return 1, 1, 1
        return 0, 0, 0
    return new_f

regs = range(6, 10, 2)

param_list = [("regs_{}".format(N), 
           Parameters(
               cis_count=2, 
               reg_channels=N,
               out_channels=3, 
               cue_channels=3 ,
               population_size=5000,
               mutation_rate=.001,
               # mutation_rate=.005,
               # mutate_type=MutateType.PROGRESSIVE,
               # add_zeros=5,
           ))
          for N in regs
          ]

class MyTreatment(Treatment):
    def run_replicate(self, replicate, lineage):
        if len(lineage.targets) == 0:
            # lineage.add_target(target(), weighting=[1, 0, 0])
            # lineage.add_target(target(), weighting=[1, 1, 0])
            lineage.add_target(target(), weighting=[1, 1, 1])
            print lineage.targets[0].as_array()
        lineage.set_target(0)
        while lineage.generation < 60000:
            g = lineage.generation
            # t = bisect([2000, 5000], g)
            # if t != lineage.target_index:
            #     lineage.set_target(t)

            if g % 100 == 0:
                print replicate.treatment.seq, replicate.seq, lineage.generation,\
                    lineage.population.worst_and_best()
            if g > 0 and g % 5000 == 0:
                replicate.draw_winners(lineage)
            lineage.next_generation()
        replicate.draw_winners(lineage)


treats = [MyTreatment(n, p, 20) for n, p in param_list]
the_exp = Experiment('/Users/Brett/Desktop', name='output_first', seed=5).add_all(*treats)


if __name__ == '__main__':
    the_exp.run_from_commandline()

