import sys

# Make bricolage available...
from pathlib import Path
sys.path.append(str(Path(__file__).absolute().parent.parent))
from functools import partial

from bricolage.experiment import Experiment, Treatment
from bricolage.threshold3 import Parameters
from bricolage.boolean_functions import get_functions

USE_FUNCTION = 4

def target():
    desc, func = get_functions()[USE_FUNCTION]
    def new_f(a, b, c):
        if func(a, b, c):
            return 1, 0, 1
        return 0, 1, 0
    return new_f

# Set defaults 
make_params = partial(
    Parameters, 
    cis_count=2, 
    out_channels=3, 
    cue_channels=3 ,
    population_size=5000,
    mutation_rate=.001,
)

class MyTreatment(Treatment):
    def run_replicate(self, replicate, lineage):
        if len(lineage.targets) == 0:
            lineage.add_target(target(), weighting=[1, 1, 1])
            print lineage.targets[0].as_array()
        lineage.set_target(0)
        lineage.extra.flow = [1, 0, 1]
        while lineage.generation < 60000:
            g = lineage.generation
            if g % 100 == 0:
                print replicate.treatment.seq, replicate.seq, lineage.generation,\
                    lineage.population.worst_and_best()
            if g > 0 and g % 5000 == 0:
                replicate.draw_winners(lineage)
            lineage.next_generation()
        replicate.draw_winners(lineage)


the_exp = Experiment('/Users/Brett/Desktop')
# the_exp.add(MyTreatment('reg_4', make_params(reg_channels=4), 20))
# the_exp.add(MyTreatment('reg_6', make_params(reg_channels=6), 20))
the_exp.add(MyTreatment('reg_8', make_params(reg_channels=8), 2))

if __name__ == '__main__':
    the_exp.run_from_commandline()

