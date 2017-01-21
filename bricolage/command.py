#!env python
import click
import cPickle as pickle
from logtools import set_logging, get_logger
from bricolage.stats import (
    StatsFitness, StatsVisitor, StatsMutualInformation, StatsRelevantControl,
    StatsLag, StatsRobustness, StatsGenerations, StatsBindings,
    StatsEnvironmental, StatsFirstWinner)
from experiment import ExperimentError
from bricolage.dot_layout import DotMaker
from bricolage.graph_maker import get_graph_by_type, GraphType
from bricolage.analysis_ext import NetworkAnalysis
import pathlib

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
only_ = click.option('--only', help="Only do this generation")
treatment_ = click.option('--treatment', default=-1,
                          help="Filter treatments by name.")
replicate_ = click.option('--replicate', default=-1,
                          help="Filter replicates by number.")

def composed(*decs):
    def deco(f):
        for dec in reversed(decs):
            f = dec(f)
        return f

stats_ = composed(verbose_, treatment_, replicate_, every_, only_)

@bricolage.command()
@verbose_
@treatment_
@replicate_
@click.option('--overwrite', is_flag=True, default=False,
              help="Trash the experiment and start again.")
def run(overwrite, verbose, treatment, replicate):
    """Run the simulation.

    This will create a new simulation or complete an existing one (if
    unfinished).
    """
    set_logging(verbose)
    try:
        the_t, the_rep = NS.experiment.find_matching(treatment, replicate)
    except ExperimentError as e:
        raise click.BadParameter(e.message)

    NS.experiment.run(overwrite=overwrite,
                      only_treatment=the_t,
                      only_replicate=the_rep)


@bricolage.command()
def list():
    """List Treatments and Replicates.
    """

    for t in NS.experiment.treatments:
        info = "{:02d} : Seed {}, {} replicates, name {}".format(
            t.seq,
            t.seed,
            len(t.replicates),
            t.name)
        click.echo(info)


# @bricolage.command()
# @every_
# @treatment_
# @replicate_
# @verbose_
# def stats(every, treatment, replicate, verbose):
#     """Gather statistics about the simulation"""
#     set_logging(verbose)
#
#     try:
#         the_t, the_rep = NS.experiment.find_matching(treatment, replicate)
#     except ExperimentError as e:
#         raise click.BadParameter(e.message)
#
#     visitor = StatsVisitor(NS.experiment,
#                            [StatsRelevantControl, StatsFitness,
#                             StatsMutualInformation, StatsBindings])
#     NS.experiment.visit_generations(visitor,
#                                     only_treatment=the_t,
#                                     only_replicate=the_rep,
#                                     every=every)


@bricolage.command()
def trash():
    """Trash the experiment.

    If you have move2trash installed you can recover it from the trash.
    """
    if click.confirm("Are you sure?"):
        NS.experiment._remove_paths()


class DrawVisitor(object):

    def __init__(self):
        super(DrawVisitor, self).__init__()

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
            self.replicate.draw_net('g_best', net,
                                    target=self.lineage.targets[0],
                                    graph_type=GraphType.GENE)


@bricolage.command()
@every_
@treatment_
@replicate_
@verbose_
@only_
def draw(every, treatment, replicate, verbose, only):
    """Draw graphs of the best networks."""
    try:
        the_t, the_rep = NS.experiment.find_matching(treatment, replicate)
    except ExperimentError as e:
        raise click.BadParameter(e.message)

    NS.experiment.visit_generations(DrawVisitor(),
                                    only_treatment=the_t,
                                    only_replicate=the_rep,
                                    every=every,
                                    only=only)


@bricolage.command()
@verbose_
@treatment_
@replicate_
def calc_lag(verbose, treatment, replicate):
    set_logging(verbose)

    try:
        the_t, the_rep = NS.experiment.find_matching(treatment, replicate)
    except ExperimentError as e:
        raise click.BadParameter(e.message)

    # First, let's find the best fitness
    NS.experiment.visit_lineages(StatsLag(NS.experiment),
                                 only_treatment=the_t,
                                 only_replicate=the_rep)

@bricolage.command()
@verbose_
@treatment_
@replicate_
def calc_gens(verbose, treatment, replicate):
    set_logging(verbose)

    try:
        the_t, the_rep = NS.experiment.find_matching(treatment, replicate)
    except ExperimentError as e:
        raise click.BadParameter(e.message)

    # First, let's find the best fitness
    NS.experiment.visit_lineages(StatsGenerations(NS.experiment),
                                 only_treatment=the_t,
                                 only_replicate=the_rep)


@bricolage.command()
@verbose_
@treatment_
@replicate_
@every_
@only_
def pop_robustness(verbose, treatment, replicate, every, only):
    set_logging(verbose)

    try:
        the_t, the_rep = NS.experiment.find_matching(treatment, replicate)
    except ExperimentError as e:
        raise click.BadParameter(e.message)

    # First, let's find the best fitness
    visitor = StatsVisitor(NS.experiment,
                           [StatsRobustness])
    NS.experiment.visit_generations(visitor,
                                    only_treatment=the_t,
                                    only_replicate=the_rep,
                                    every=every,
                                    only=only)

@bricolage.command()
@verbose_
@treatment_
@replicate_
@every_
@only_
def env_robustness(verbose, treatment, replicate, every, only):
    set_logging(verbose)

    try:
        the_t, the_rep = NS.experiment.find_matching(treatment, replicate)
    except ExperimentError as e:
        raise click.BadParameter(e.message)

    # First, let's find the best fitness
    visitor = StatsVisitor(NS.experiment,
                           [StatsEnvironmental])
    NS.experiment.visit_generations(visitor,
                                    only_treatment=the_t,
                                    only_replicate=the_rep,
                                    every=every,
                                    only=only)

@bricolage.command()
@verbose_
@treatment_
@replicate_
def first_winner(verbose, treatment, replicate):
    set_logging(verbose)

    try:
        the_t, the_rep = NS.experiment.find_matching(treatment, replicate)
    except ExperimentError as e:
        raise click.BadParameter(e.message)

    # First, let's find the first fitness
    NS.experiment.visit_lineages(StatsFirstWinner(NS.experiment),
                                 only_treatment=the_t,
                                 only_replicate=the_rep)

@bricolage.command()
@click.argument('treatment', type=int)
@click.argument('replicate', type=int)
@click.argument('network', type=int)
def export_network(treatment, replicate, network):
    """Export an network."""
    try:
        the_t, the_rep = NS.experiment.find_matching(treatment, replicate)
    except ExperimentError as e:
        raise click.BadParameter(e.message)

    with the_rep.get_lineage() as lin:
        net = lin.get_network(network)
        name = "network-{}".format(net.identifier)
        with open('{}.pickle'.format(name), 'wb') as f:
            pickle.dump(net, f, -1)


@bricolage.command()
@click.argument('treatment', type=int)
@click.argument('replicate', type=int)
@click.argument('network', type=int)
def export_dot(treatment, replicate, network):
    """Export an network."""
    try:
        the_t, the_rep = NS.experiment.find_matching(treatment, replicate)
    except ExperimentError as e:
        raise click.BadParameter(e.message)

    with the_rep.get_lineage() as lin:
        net = lin.get_network(network)
        ana = NetworkAnalysis(net)
        curpath = pathlib.Path('.')
        name = curpath / "T{}-R{}-network-{}".format(treatment, replicate, net.identifier)
        g = get_graph_by_type(GraphType.GENE, ana, simple_labels=True)
        d = DotMaker(g)
        d.save_dot(name.with_suffix('.dot').as_posix())
        g = get_graph_by_type(GraphType.GENE, ana)
        d = DotMaker(g)
        d.save_picture(name.with_suffix('.png').as_posix())

@bricolage.command()
@click.argument('treatment', type=int)
@click.argument('replicate', type=int)
def best_diff(treatment, replicate):
    """Export an network."""
    try:
        the_t, the_rep = NS.experiment.find_matching(treatment, replicate)
    except ExperimentError as e:
        raise click.BadParameter(e.message)

    with the_rep.get_lineage() as lin:
        gfirst = lin.first_winning_generation()
        gen = lin.get_generation(gfirst)
        log.info("First generation is {}".format(gen))
        cur_net = gen.get_best(1)[0]

        anc = lin.get_ancestry(cur_net.identifier)
        first_net = anc[0]
        log.info("Loaded ancestry")

        cur_ana = NetworkAnalysis(cur_net)
        first_ana = NetworkAnalysis(first_net)
        cur_graph = get_graph_by_type(GraphType.GENE, cur_ana)
        first_graph = get_graph_by_type(GraphType.GENE, first_ana)

        curpath = pathlib.Path('.')
        name = curpath / "T{}-R{}-network-{}".format(treatment, replicate, cur_net.identifier)
        d = DotMaker(cur_graph)
        d.save_diff_dot(name.with_suffix('.dot').as_posix(), first_graph)
        d.save_diff_picture(name.with_suffix('.png').as_posix(), first_graph)
        # g = get_graph_by_type(GraphType.GENE, ana)
        # d = DotMaker(g)
        # d.save_picture(name.with_suffix('.png').as_posix())

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
