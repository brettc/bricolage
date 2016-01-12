#!env python
import click
from bricolage.analysis import (StatsAverageControl, StatsFitness,
                                StatsVisitor, StatsMutualInformation)


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


@bricolage.command()
@click.option('--overwrite', is_flag=True, default=False,
              help="Trash the experiment and start again.")
def run(overwrite):
    """Run the simulation.

    This will create a new simulation or complete an existing one (if unfinished).
    """
    NS.experiment.run(overwrite=overwrite)


@bricolage.command()
@click.option('--every', default=100)
def stats(every):
    """Gather statistics about the simulation"""
    visitor = StatsVisitor(NS.experiment,
                           [StatsAverageControl, StatsFitness, StatsMutualInformation])
    NS.experiment.visit_generations(visitor, every=every)

# VERBOSITY
# logging.getLogger('sqlalchemy.engine').setLevel(logging.INFO)

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
# @add_command()
# def trash(self, args):
#     """Trash the experiment.
#     NOTE: If you have move2trash installed you can recover it from the
#     trash.
#     """
#     self._remove_paths()
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
