import sys
from pathlib import Path
sys.path.append(str(Path(__file__).absolute().parent.parent))

from bricolage.logtools import init_logging, get_logger
init_logging()
log = get_logger()

from bricolage.experiment import Experiment, Treatment, add_argument
from bricolage.threshold3 import Parameters
from bricolage.analysis import InfoSummarizer

def target1(a, b):
    if a and b:
        return 1.0
    return 0.0

class TestTreatment(Treatment):
    def run_replicate(self, replicate, lineage):
        if len(lineage.targets) == 0:
            lineage.add_target(target1)
        lineage.extra.flow = [1.0]
        while lineage.generation < 10000:
            lineage.next_generation()

p = Parameters(
    cis_count=2, 
    reg_channels=2,
    out_channels=1, 
    cue_channels=2, 
    population_size=100,
    mutation_rate=.001,
)


class MyExperiment(Experiment):
    @add_argument('--verbose', action='store_true')
    @add_argument('--every', type=int, default=10)
    def test_info(self, args):
        self.database.create(args.verbose)
        for rep, lin in self.iter_lineages():
            targ = lin.targets[-1]
            flow = lin.extra.flow
            infosum = InfoSummarizer(lin, targ, flow)
            rep.clean_stats(infosum.get_names())
            for n, gen in lin.all_generations(every=args.every):
                rep.write_stats(n, infosum.get_values(gen))
            self.database.session.commit()


the_exp = MyExperiment('/Users/brett/Desktop', 'bob', seed=5).add_all(
    TestTreatment('one', p, 3),
    TestTreatment('two', p, 3),
)

if __name__ == '__main__':
    the_exp.run_from_commandline()

