from bricolage.experiment import Experiment
from bricolage.stats import (StatsAverageControl, StatsFitness,
                             StatsVisitor, StatsMutualInformation,
                             StatsOutputControl)


# def test_database(simple_experiment):
#     assert isinstance(simple_experiment, Experiment)
#
#     exp = simple_experiment
#     exp.database.create(True)
#
#     visitor = StatsVisitor(exp, [StatsMutualInformation])
#     exp.visit_generations(visitor, every=10)

    # This reloads everything, so it will exclude what has been saved.
    # visitor = StatsVisitor(exp, [StatsAverageControl, StatsFitness])
    # visitor = StatsVisitor(exp, [StatsOutputControl])
    # exp.visit_generations(visitor, every=10)
