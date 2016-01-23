#!env python
import click
from logtools import set_logging, get_logger
from bricolage.stats import (
    StatsFitness, StatsVisitor, StatsMutualInformation, StatsOutputControl,
    StatsBindings)
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
                            StatsMutualInformation, StatsBindings])
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


class FindFirstFitVisitor(object):
    def __init__(self, experiment):
        self.experiment = experiment
        self.session = experiment.database.session
        self.replicate = None
        self.lineage = None
        self.oc_analyzer = None
        self.mi_analyzer = None
        self.done = {}

        self.first_best = 'FIRST_BEST'
        self.first_control = 'FIRST_CONTROL'

        # TODO: should only load relevant ones
        for srep in self.session.query(StatsReplicateRecord).all():
            self.done[(srep.treatment_id, srep.replicate_id, srep.kind)] = srep

    def is_done(self, rep, kind):
        return (rep.treatment.seq, rep.seq, kind) in self.done

    def visit_lineage(self, rep, lin):
        log.info("{}".format(rep)).push().add()

        # Setup
        self.replicate = rep
        self.lineage = lin
        self.oc_analyzer = OutputControlAnalyzer(lin.world)
        self.mi_analyzer = MutualInfoAnalyzer(lin.world, lin.targets[0].calc_categories())

        # Analysis
        if not self.is_done(rep, self.first_best):
            fgen = self.find_first_winner()
            self.session.add(StatsReplicateRecord(rep, self.first_best, fgen))

        if not self.is_done(rep, self.first_control):
            cgen = self.find_first_control()
            self.session.add(StatsReplicateRecord(rep, self.first_control, cgen))

        self.session.commit()

        log.pop()

    def find_first_winner(self):
        first_win = None
        for g in self.lineage._generations.where("best == 1.0"):
            first_win = g['generation']
            break

        log.info("Found first winner at generation {}".format(first_win))
        return first_win

    def find_first_control(self):
        # Load last generation and check if we got a master gene
        ct = self.get_controlled_indexes(self.lineage.population)
        if ct is None:
            log.debug("No control evolved")
            return None
        log.debug("We got a master gene at the last generation.")

        # Grab the first one and load the ancestry
        generations = []
        for i, net_index in enumerate(ct):
            net = self.lineage.population[net_index]
            first = self.get_ancestry_control(net)
            generations.append(first)
            if i == 0:
                # Just one for now. It looks like the coalesce!
                break

        # Just one
        log.info("First control is at generation {}.".format(generations[0]))
        return generations[0]

    def get_controlled_indexes(self, collection):
        ai = self.oc_analyzer.numpy_info_from_collection(collection)
        mi = self.mi_analyzer.numpy_info_from_collection(collection)
        mi.shape = mi.shape[:-1]
        left_to_explain = ai[:, :, 1] - ai[:, :, 0]
        control = mi - left_to_explain
        controlled = numpy.isclose(control, 1.0)
        if not controlled.any():
            return None

        # Just grab the network indexes (we don't are which gene it was)
        where = numpy.where(controlled)[0]
        return where

    def get_ancestry_control(self, network):
        log.debug("loading ancestry of winner {}".format(network.identifier))
        anc = self.lineage.get_ancestry(network.identifier)

        log.debug("Calculating ancestry control of {} networks.".format(anc.size))
        ct = self.get_controlled_indexes(anc)

        # We know that at least the last one should be good!
        assert ct is not None

        first_ancestor = anc[ct[0]]

        # Make sure the network has a fitnes
        self.lineage.targets[0].assess(first_ancestor)
        log.debug("First ancestor with master gene is at {}.".format(
            first_ancestor.generation))
        self.replicate.draw_net('first-control', first_ancestor,
                                first_ancestor.generation, signals=False)
        return first_ancestor.generation


@bricolage.command()
@verbose_
def calc_lag(verbose):
    set_logging(verbose)

    # First, let's find the best fitness
    NS.experiment.visit_lineages(FindFirstFitVisitor(NS.experiment))


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
