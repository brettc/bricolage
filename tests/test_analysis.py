from bricolage.experiment import Experiment
from bricolage.analysis import (StatsAverageControl, StatsFitness,
                                StatsVisitor, StatsMutualInformation)


def test_database(simple_experiment):
    assert isinstance(simple_experiment, Experiment)

    exp = simple_experiment
    exp.database.create(True)

    visitor = StatsVisitor(exp, [StatsMutualInformation])
    exp.visit_generations(visitor, every=10)

    # This reloads everything, so it will exclude what has been saved.
    visitor = StatsVisitor(exp, [StatsAverageControl, StatsFitness])
    exp.visit_generations(visitor, every=10)
