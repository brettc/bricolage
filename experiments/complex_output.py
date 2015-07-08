import sys

# Make bricolage available...
from pathlib import Path
sys.path.append(str(Path(__file__).absolute().parent.parent))

from bricolage.experiment import Experiment, Treatment, add_argument
from bricolage.threshold3 import Parameters
from boolean_functions import get_functions

USE_FUNCTION = 5

def target():
    desc, func = get_functions()[USE_FUNCTION]
    def new_f(a, b, c):
        if func(a, b, c):
            return 1, .5, .25
        return .5, 0, 0
    return new_f

regs = range(8, 10, 2)

param_list = [("regs_{}".format(N), 
           Parameters(
               cis_count=2, 
               reg_channels=N,
               out_channels=3, 
               cue_channels=3 ,
               population_size=1000,
               mutation_rate=.001))
          for N in regs
          ]

class MyTreatment(Treatment):
    def run_replicate(self, replicate, lineage):
        if len(lineage.targets) == 0:
            lineage.add_target(target())
            print lineage.targets[0].as_array()
        lineage.extra.flow = [1, .5, .25]
        while lineage.generation < 40000:
            g = lineage.generation
            if g % 100 == 0:
                print replicate.treatment.seq, replicate.seq, lineage.generation,\
                    lineage.population.worst_and_best()
            if g > 0 and g % 5000 == 0:
                replicate.draw_winners(lineage)
            lineage.next_generation()
        replicate.draw_winners(lineage)


treats = [MyTreatment(n, p, 20) for n, p in param_list]
the_exp = Experiment('/Users/Brett/Desktop', seed=5).add_all(*treats)


if __name__ == '__main__':
    the_exp.run_from_commandline()

