#!env python
import click
import logging
from bricolage.stats import (StatsFitness, StatsVisitor,
                             StatsMutualInformation, StatsOutputControl)


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

verbose = click.option('--verbose', is_flag=True, default=False,
                       help="Show debug output.")


@bricolage.command()
@verbose
@click.option('--overwrite', is_flag=True, default=False,
              help="Trash the experiment and start again.")
def run(overwrite, verbose):
    """Run the simulation.

    This will create a new simulation or complete an existing one (if unfinished).
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

        # We need to do this separately. This is equivalent to setting "echo"
        # on the database creation.
        logging.getLogger('sqlalchemy.engine').setLevel(logging.INFO)

    NS.experiment.run(overwrite=overwrite)


@bricolage.command()
@click.option('--every', default=100)
@click.option('--treatment', default="", help="Filter treatments by name.")
@click.option('--replicate', default=-1, help="Filter replicates by number.")
@verbose
def stats(every, treatment, replicate, verbose):
    """Gather statistics about the simulation"""
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

        # We need to do this separately. This is equivalent to setting "echo"
        # on the database creation.
        logging.getLogger('sqlalchemy.engine').setLevel(logging.INFO)

    # Find the treatment.
    the_t = None
    if treatment != "":
        the_t = NS.experiment.find_matching_treatment(treatment)
        if the_t is None:
                raise click.BadParameter(
                    "Can't find unique treatment name starting: '{}".format(treatment))
    if replicate < 1:
        replicate = None

    visitor = StatsVisitor(NS.experiment,
                           [StatsOutputControl, StatsFitness, StatsMutualInformation])
    NS.experiment.visit_generations(visitor,
                                    only_treatment=the_t,
                                    only_replicate=replicate,
                                    every=every)


@bricolage.command()
def trash():
    """Trash the experiment.

    If you have move2trash installed you can recover it from the trash.
    """
    if click.confirm("Are you sure?"):
        NS.experiment._remove_paths()

# @bricolage.command()
# @click.option('--verbose', is_flag=True)
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
