import logtools
log = logtools.get_logger()

import copy
import shutil
import random
import argparse

from .analysis_ext import NetworkAnalysis
import graph

try:
    from send2trash import send2trash
except ImportError:
    send2trash = None


from pathlib import Path
from lineage import FullLineage, SnapshotLineage
from experimentdb import Database, TreatmentRecord, ReplicateRecord

class ExperimentError(RuntimeError):
    pass


def _remove_folder(path):
    if not path.exists():
        log.warning("Trying to remove non-existent path: {}".format(str(path)))
        return
    
    if send2trash is not None:
        log.info("Trashing {} (you can find it in the trash)".format(str(path)))
        send2trash(str(path))
    else:
        log.info("Deleting {}".format(str(path)))
        shutil.rmtree(str(path))


def _construct_path(path, name):
    if not isinstance(path, Path):
        path = Path(path)
    path = path / name
    return path.absolute()


class Replicate(object):
    """Wraps a lineage into a sub-folder"""

    FILENAME = 'lineage.h5'

    def __init__(self, treatment, seq):
        self.treatment = treatment
        self.seq = seq
        self.seed = treatment.rng.randint(0, 1 << 16)
        self.fresh = True
        self.params = copy.deepcopy(self.treatment.params)
        self.params.seed = self.seed
        self.generations = -1

    @property
    def session(self):
        return self.treatment.experiment.database.session

    def run(self):
        # TODO: Maybe use this, instead of merging
        # record = self.session.query(ReplicateRecord)\
        #     .filter(ReplicateRecord.treatment_id == self.treatment.seq)\
        #     .filter(ReplicateRecord.replicate_id == self.seq).one()
        if not self.path.exists():
            self.path.mkdir()
        with self.get_lineage() as lin:
            try:
                self.treatment.run_replicate(self, lin)
            except (KeyboardInterrupt, SystemExit):
                # Try and exit gracefully with a keyboard interrupt
                self.treatment.experiment.user_interrupt = True

            self.generations = lin.generation
            self.session.merge(ReplicateRecord(self))
            self.session.commit()

    @property
    def dirname(self):
        return "{:03d}".format(self.seq)

    @property
    def path(self):
        return self.treatment.path / self.dirname

    @property
    def analysis_path(self):
        return self.treatment.analysis_path / self.dirname

    @property
    def lineage_name(self):
        return self.path / self.FILENAME

    def has_started(self):
        return self.lineage_name.exists()

    def get_lineage(self, readonly=False):
        p = self.path
        if not p.exists():
            if readonly:
                log.error("Replicate doesn't exist, {}".format(self.name))
                raise ExperimentError

            p.mkdir()

        # We don't use this, but users of the class might depend on it
        if not self.analysis_path.exists():
            self.analysis_path.mkdir()

        db_path = self.lineage_name

        # Load the actual lineage
        if db_path.exists():
            lin = self.treatment.experiment.lineage_class(
                db_path, readonly=readonly)

            # If we reloaded then the seed should be the same as what was
            # generated by the original parameters
            if not lin.params.same_as(self.params):
                lin.close()
                log.error("Parameters have changed!")
                raise ExperimentError

            log.info("Loaded lineage {}".format(lin))

            # Kill this flag
            self.fresh = False

        if self.fresh:
            lin = self.treatment.experiment.lineage_class(db_path, self.params)
            log.info("Created lineage {}".format(lin))

        return lin


class Treatment(object):
    """Replicate a set of experiments into different folders"""

    def __init__(self, name, params, count):
        self.name = name
        self.params = params
        self.count = count
        self.experiment = None
        self.seq = None
        self.seed = None
        self.rng = None
        self.replicates = None

    def _bind(self, experiment):
        assert self.experiment is None
        self.experiment = experiment
        self.seq = experiment.get_next_treatment_seq()
        self.seed = experiment.rng.randint(0, 1 << 16)
        self.rng = random.Random(self.seed)
        self.replicates = [Replicate(self, i) for i in range(1, self.count + 1)]

    @property
    def dirname(self):
        return "{0.seq:02d}_{0.name}".format(self)

    @property
    def path(self):
        return self.experiment.path / self.dirname

    @property
    def analysis_path(self):
        return self.experiment.analysis_path / self.dirname

    def run(self):
        # Check paths...
        if not self.path.exists():
            self.path.mkdir()
        if not self.analysis_path.exists():
            self.analysis_path.mkdir()

        for r in self.replicates:
            if not self.experiment.user_interrupt:
                r.run()


class add_command(object):
    subparsers = {}
    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    def __call__(self, func):
        if 'help' not in self.kwargs:
            self.kwargs['help'] = func.__doc__
        self.subparsers[func] = self.args, self.kwargs
        return func


class add_argument(object):
    options = {}
    def __init__(self, name, *args, **kwargs):
        self.name = name
        self.args = args
        self.kwargs = kwargs

    def __call__(self, func):
        ops = self.options.setdefault(func, [])
        ops.append((self.name, self.args, self.kwargs))

        # Make a command if it isn't there
        if func not in add_command.subparsers:
            add_command()(func)

        return func


# TODO: Externalise the root
class Experiment(object):
    def __init__(self, path, name, seed, analysis_path=None, full=True):
        logtools.init_logging()

        self.path = _construct_path(path, name)
        if analysis_path:
            self.analysis_path = _construct_path(analysis_path, name)
        else:
            self.analysis_path = self.path

        self.next_treatment_seq = 0
        self.treatments = []

        if full:
            self.lineage_class = FullLineage
        else:
            self.lineage_class = SnapshotLineage

        self.seed = seed
        self.rng = random.Random(seed)

        # We create the database on the analysis path --- oooh! Clever.
        self.database = Database(self.analysis_path)
        self.user_interrupt = False

    def add(self, treatment):
        self.treatments.append(treatment)
        treatment._bind(self)
        return self

    def add_all(self, *treats):
        for t in treats:
            self.add(t)
        return self

    def get_next_treatment_seq(self):
        self.next_treatment_seq += 1
        return self.next_treatment_seq

    def _create_path(self, path):
        if not path.parent.exists():
            log.error("No parent path: {}".format(str(path.parent)))
            raise ExperimentError

        log.info("Creating path {}".format(str(path)))
        path.mkdir()

    def _remove_paths(self):
        if self.path.exists():
            _remove_folder(self.path)
        if self.path is not self.analysis_path:
            if self.analysis_path.exists():
                _remove_folder(self.analysis_path)
            self.analysis_path.mkdir()

    def _init_db(self):
        for t in self.treatments:
            self.database.session.merge(TreatmentRecord(t))
            for r in t.replicates:
                self.database.session.merge(ReplicateRecord(r))
        self.database.session.commit()

    def iter_lineages(self):
        for t in self.treatments:
            for r in t.replicates:
                with r.get_lineage() as lin:
                    yield r, lin

    def _make_parser(self):
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers(title="Available Commands")
        for func, (args, kwargs) in add_command.subparsers.items():
            subp = subparsers.add_parser(func.__name__, *args, **kwargs)
            subp.set_defaults(func=func)
            if func in add_argument.options:
                for name, optargs, optkwargs in add_argument.options[func]:
                    subp.add_argument(name, *optargs, **optkwargs)
        return parser

    def run_from_commandline(self):
        parser = self._make_parser()
        args = parser.parse_args()
        try:
            args.func(self, args)
        except ExperimentError:
            exit(1)
        exit(0)

    # --------------------
    # Commands
    @add_argument('--overwrite', action="store_true", 
                  help="Trash the entire experiment and start again.")
    @add_argument('--dry', action="store_true", 
                  help="Just create the database, don't run the simulations.")
    @add_argument('--verbose', action="store_true", help="Say lots.")
    def run(self, args):
        """Run the experiment."""
        if args.overwrite:
            self._remove_paths()

        if not self.path.exists():
            self._create_path(self.path)
            if self.path != self.analysis_path:
                self._create_path(self.analysis_path)

        self.database.create()
        self._init_db()

        if args.dry:
            log.info("Dry run -- exiting without simulating")
            return

        # Now run
        log.info("Beginning experiment in {}".format(str(self.path)))
        for treat in self.treatments:
            if not self.user_interrupt:
                treat.run()

        if self.user_interrupt:
            log.info("User interrupted --- quitting")

    @add_argument('--verbose', action="store_true", help="Say lots.")
    def status(self, args):
        """Show status of experiment"""
        if not self.path.exists() or not self.database.path.exists():
            log.error("Experiment does not exist yet. You need to run it")
            raise ExperimentError

        self.database.create(args.verbose)
        trecords = self.database.session.query(TreatmentRecord).all()
        for t in trecords:
            print t
            for r in t.replicates:
                print r

    @add_command()
    def draw_top(self, args):
        """Draw top graphs from final populations."""
        for rep, lin in self.iter_lineages():
            winners = [(n.fitness, n) for n in lin.population]
            winners.sort(reverse=True)
            for i, (fit, net) in enumerate(winners):
                if i == 3:# args.maxn:
                    break
                draw_net(rep, 'winner', net, lin.generation)


    @add_command()
    def trash(self, args):
        """Trash the experiment.
        NOTE: If you have move2trash installed you can recover it from the trash.
        """
        self._remove_paths()


def draw_net(rep, prefix, net, gen):
        ana = NetworkAnalysis(net)
        g = graph.FullGraph(ana)#, knockouts=False)
        d = graph.DotMaker(g)
        p = rep.analysis_path / '{}-G{:07d}-N{:02d}-F{}.png'.format(
            prefix, gen, net.identifier, net.fitness)
        log.info("writing {}".format(str(p)))
        d.save_picture(str(p))
        d.save_dot(str(p.with_suffix('.dot')))
