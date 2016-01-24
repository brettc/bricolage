#!env python
import click
from logtools import set_logging, get_logger
from bricolage.stats import (
    StatsFitness, StatsVisitor, StatsMutualInformation, StatsOutputControl,
    StatsAverageControl, StatsLag)
from .analysis_ext import OutputControlAnalyzer, MutualInfoAnalyzer
from experiment import ExperimentError
from experimentdb import StatsReplicateRecord
import numpy

log = get_logger()


class NS(object):
    """Just a namespace"""
    experiment = None


def run_from_commandline(exp):
    """Assign the global experiment, and then run the commands using click"""
    NS.experiment = exp
    bricolage()


@click.group(chain=True)
def bricolage():
    pass


verbose_ = click.option('--verbose', is_flag=True, default=False,
                        help="Show debug output.")
every_ = click.option('--every', default=1000,
                      help="Only do every N generations")
treatment_ = click.option('--treatment', default="",
                          help="Filter treatments by name.")
replicate_ = click.option('--replicate', default=-1,
                          help="Filter replicates by number.")


@bricolage.command()
@verbose_
@click.option('--overwrite', is_flag=True, default=False,
              help="Trash the experiment and start again.")
def run(overwrite, verbose):
    """Run the simulation.

    This will create a new simulation or complete an existing one (if
    unfinished).
    """
    set_logging(verbose)
    NS.experiment.run(overwrite=overwrite)


@bricolage.command()
@every_
@treatment_
@replicate_
@verbose_
def stats(every, treatment, replicate, verbose):
    """Gather statistics about the simulation"""
    set_logging(verbose)

    try:
        the_t, the_rep = NS.experiment.find_matching(treatment, replicate)
    except ExperimentError as e:
        raise click.BadParameter(e.message)

    visitor = StatsVisitor(NS.experiment,
                           [StatsOutputControl, StatsFitness,
                            StatsMutualInformation, StatsAverageControl])
    NS.experiment.visit_generations(visitor,
                                    only_treatment=the_t,
                                    only_replicate=the_rep,
                                    every=every)


@bricolage.command()
def trash():
    """Trash the experiment.

    If you have move2trash installed you can recover it from the trash.
    """
    if click.confirm("Are you sure?"):
        NS.experiment._remove_paths()


class DrawVisitor(object):
    def wants_generation(self, gen_num):
        return True

    def visit_lineage(self, rep, lin):
        self.replicate = rep
        self.lineage = lin

    def visit_generation(self, gen_num, pop):
        winners = [(n.fitness, n.identifier, n) for n in pop]
        winners.sort(reverse=True)
        for i, (fit, ident, net) in enumerate(winners):
            if i == 1:
                break
            self.replicate.draw_net('best', net, gen_num,
                                    knockouts=True,
                                    target=self.lineage.targets[0])


@bricolage.command()
@every_
@treatment_
@replicate_
@verbose_
def draw(every, treatment, replicate, verbose):
    """Draw graphs of the best networks."""
    try:
        the_t, the_rep = NS.experiment.find_matching(treatment, replicate)
    except ExperimentError as e:
        raise click.BadParameter(e.message)

    NS.experiment.visit_generations(DrawVisitor(),
                                    only_treatment=the_t,
                                    only_replicate=the_rep,
                                    every=every)


@bricolage.command()
@verbose_
@treatment_
@replicate_
def calc_lag(verbose, treatment, replicate):
    set_logging(verbose)

    # try:
    #     the_t, the_rep = NS.experiment.find_matching(treatment, replicate)
    # except ExperimentError as e:
    #     raise click.BadParameter(e.message)

    # First, let's find the best fitness
    NS.experiment.visit_lineages(StatsLag(NS.experiment))
                                 # only_treatment=the_t,
                                 # only_replicate=the_rep)


# def status(verbose):
#     """Show status of experiment"""
#     exp = NS.experiment
#     if not exp.path.exists() or not exp.database.path.exists():
#         raise click.UsageError(
#             "Experiment does not exist yet. You need to run it")
#
#     exp.database.create(verbose)
#     trecords = exp.database.session.query(TreatmentRecord).all()
#     for t in trecords:
#         print t
#         for r in t.replicates:
#             print r
#
# @add_command()
# def draw_top(self, args):
#     """Draw top graphs from final populations."""
#     for rep, lin in self.iter_lineages():
#         rep.draw_winners(lin)
#
#
# @verbose
# @add_argument('--every', type=int, default=25)
# def calc_stats(self, args):
#     self.database.create(args.verbose)
#     for rep, lin in self.iter_lineages():
#         targ = lin.targets[-1]
#         infosum = InfoSummarizer(lin, targ)
#         rep.clean_stats(infosum.get_names())
#         for n, gen in lin.all_generations(every=args.every):
#             rep.write_stats(n, infosum.get_values(gen))
#             print rep.treatment.seq, rep.seq, n
#         self.database.session.commit()
#
# @verbose
# def calc_neighbours(self, args):
#     self.database.create(args.verbose)
#     for rep, lin in self.iter_lineages():
#         targ = lin.targets[-1]
#         naysum = NeighbourhoodSummarizer(lin, targ)
#         rep.clean_stats(naysum.get_names())
#         vals = naysum.get_values(lin.population)
#         rep.write_stats(lin.generation, vals)
#         print rep.treatment.seq, rep.seq, vals
#         self.database.session.commit()
#
# @verbose
# @add_argument('--every', type=int, default=25)
# def calc_attractor_len(self, args):
#     self.database.create(args.verbose)
#
#     maxlen = 0
#     for rep, lin in self.iter_lineages():
#         print rep, lin
#         for gnum, gen in lin.all_generations(every=args.every):
#             for net in gen:
#                 l = net.attractors_size.max()
#                 if l > maxlen:
#                     print "Max", l
#                     for a in net.attractors:
#                         print a
#                     print net.rates
#                     maxlen = l
