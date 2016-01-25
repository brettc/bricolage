from bricolage.graph_layout import DotMaker
from bricolage.graph_maker import get_graph_by_type, GraphType
import logtools

import copy
import shutil
import random
import inspect

from analysis import NetworkAnalysis

try:
    from send2trash import send2trash
except ImportError:
    send2trash = None

from pathlib import Path
from lineage import FullLineage, SnapshotLineage
from experimentdb import (Database, TreatmentRecord, ReplicateRecord,
                          StatsRecord)
from sqlalchemy.sql import and_, delete

log = logtools.get_logger()


class ExperimentError(RuntimeError):
    pass


def _make_name_from_program():
    # Get caller of the caller ...
    stackinfo = inspect.stack()[2]
    fname = stackinfo[1]
    # should be a python path
    prog = Path(fname)
    return prog.name[:prog.name.find('.')]


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

    def __repr__(self):
        return "<Replicate: {}-{}>".format(
            self.treatment.name,
            self.seq,
        )

    @property
    def session(self):
        if not hasattr(self, '_session'):
            self._session = self.treatment.experiment.database.session
        return self._session

    def run(self):
        # TODO: Maybe use this, instead of merging
        # record = self.session.query(ReplicateRecord)\
        #     .filter(ReplicateRecord.treatment_id == self.treatment.seq)\
        #     .filter(ReplicateRecord.replicate_id == self.seq).one()

        log.info("Running replicate {}".format(self.path))
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
                raise ExperimentError("Replicate doesn't exist, {}".format(self.path.name))

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
                raise ExperimentError("Parameters have changed!")

            log.info("Loaded lineage {}".format(lin))

            # Kill this flag
            self.fresh = False

        if self.fresh:
            lin = self.treatment.experiment.lineage_class(db_path, self.params)
            log.info("Created lineage {}".format(lin))

        return lin

    def clean_stats(self, tags=None):
        eng = self.treatment.experiment.database.engine
        se = delete(StatsRecord).where(and_(
            # StatsRecord.kind.in_(names),
            StatsRecord.replicate_id == self.seq,
            StatsRecord.treatment_id == self.treatment.seq))
        eng.execute(se)

        ge = delete(StatsRecord).where(and_(
            # StatsRecord.kind.in_(names),
            StatsRecord.replicate_id == self.seq,
            StatsRecord.treatment_id == self.treatment.seq))
        eng.execute(ge)

    def write_stats(self, gen, namevalues):
        stats = [StatsRecord(self, gen, n, v) for n, v in namevalues]
        self.session.add_all(stats)
        # log.info("Writing stats for generation {}".format(gen))
        # self.session.bulk_save_objects(stats)

    def draw_winners(self, lin, maxw=3):
        winners = [(n.fitness, n.identifier, n) for n in lin.population]
        winners.sort(reverse=True)
        for i, (fit, ident, net) in enumerate(winners):
            if i == maxw:
                break
            self.draw_net('winner', net, target=lin.targets[0])

    def draw_net(self, prefix, net,
                 graph_type=GraphType.GENE_SIGNAL,
                 knockouts=True,
                 target=None):
        ana = NetworkAnalysis(net)
        ana.calc_output_control()
        if target:
            ana.calc_mutual_info(target)
        g = get_graph_by_type(graph_type, ana)
        d = DotMaker(g)
        p = self.analysis_path / '{}-G{:07d}-N{:02d}-F{}.png'.format(
            prefix, net.generation, net.identifier, net.fitness)
        log.info("Writing {}".format(str(p)))
        d.save_picture(str(p))
        # d.save_dot(str(p.with_suffix('.dot')))


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
        self.replicates = [Replicate(self, i) for i in
                           range(1, self.count + 1)]

    @property
    def dirname(self):
        return "{0.seq:02d}_{0.name}".format(self)

    @property
    def path(self):
        return self.experiment.path / self.dirname

    @property
    def analysis_path(self):
        return self.experiment.analysis_path / self.dirname

    def with_replicate_id(self, rep_id):
        if rep_id > len(self.replicates):
            raise ExperimentError("{} has no replicate with id {}".format(self, rep_id))

        rep = self.replicates[rep_id - 1]
        assert rep.seq == rep_id
        return rep

    def run(self, only_replicate=None):
        # Check paths...
        log.debug("Running Treatment {}".format(self))
        if not self.path.exists():
            self.path.mkdir()
        if not self.analysis_path.exists():
            self.analysis_path.mkdir()

        for r in self.replicates:
            if only_replicate is not None and r is not only_replicate:
                continue
            if not self.experiment.user_interrupt:
                r.run()


class Experiment(object):
    def __init__(self, path, treatments, name=None, seed=1, analysis_path=None, full=True):
        logtools.set_logging()

        if name is None:
            name = _make_name_from_program()

        self.path = _construct_path(path, name)
        if analysis_path:
            self.analysis_path = _construct_path(analysis_path, name)
        else:
            self.analysis_path = self.path

        self.next_treatment_seq = 0
        if full:
            self.lineage_class = FullLineage
        else:
            self.lineage_class = SnapshotLineage

        self.seed = seed
        self.rng = random.Random(seed)

        self.user_interrupt = False

        self.treatments = []
        for t in treatments:
            assert isinstance(t, Treatment)
            self.treatments.append(t)
            t._bind(self)

        self._construct()

    def _construct(self):
        # We create the database on the analysis path --- oooh! Clever.
        if not self.path.exists():
            self._create_path(self.path)
            if self.path != self.analysis_path:
                if not self.analysis_path.exists():
                    self._create_path(self.analysis_path)
        self.database = Database(self.analysis_path)
        self.database.create()
        self._init_db()

    def get_next_treatment_seq(self):
        self.next_treatment_seq += 1
        return self.next_treatment_seq

    def _create_path(self, path):
        if not path.parent.exists():
            raise ExperimentError("No parent path: {}".format(str(path.parent)))

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
        # Use merge to allow for nice things to happen
        # TODO: We don't deal with the deletion of old records if the
        # replicate number are changed
        for t in self.treatments:
            self.database.session.merge(TreatmentRecord(t))
            for r in t.replicates:
                self.database.session.merge(ReplicateRecord(r))

        # TODO: we should really update the status of this stuff too?
        self.database.session.commit()

    def iter_lineages(self, readonly=False, skip_errors=True):
        for t in self.treatments:
            for r in t.replicates:
                with r.get_lineage(readonly=readonly) as lin:
                    yield r, lin

    def run(self, overwrite=False, verbose=False, dry=False,
            only_treatment=None, only_replicate=None):
        """Run the experiment."""
        if overwrite:
            self._remove_paths()
            self._construct()

        if dry:
            log.info("Dry run -- exiting without simulating")
            return

        # Now run
        log.info("Beginning experiment in {}".format(str(self.path)))
        for treat in self.treatments:
            if only_treatment is not None and treat != only_treatment:
                continue
            if not self.user_interrupt:
                treat.run(only_replicate)

        if self.user_interrupt:
            log.info("User interrupted --- quitting")

    def visit_lineages(self, visitor,
                       only_treatment=None,
                       only_replicate=None):

        for treat in self.treatments:
            if only_treatment is not None and treat != only_treatment:
                continue

            for rep in treat.replicates:
                if only_replicate is not None and rep != only_replicate:
                    continue

                with rep.get_lineage() as lin:
                    visitor.visit_lineage(rep, lin)

    def visit_generations(self, visitor, every=1,
                          only_treatment=None,
                          only_replicate=None):
        """Try and visit everything with the least amount of loading"""

        for treat in self.treatments:
            if only_treatment is not None and treat != only_treatment:
                continue

            for rep in treat.replicates:
                if only_replicate is not None and rep != only_replicate:
                    continue

                with rep.get_lineage() as lin:
                    visitor.visit_lineage(rep, lin)

                    # Now iterate through the generations
                    gen_num = 0
                    while gen_num <= lin.generation:
                        # Check if the visitor wants this generation (as loading is
                        # expensive)
                        if visitor.wants_generation(gen_num):
                            pop = lin.get_generation(gen_num)
                            visitor.visit_generation(gen_num, pop)
                        gen_num += every

    def find_matching(self, text, repnum):
        # Default to ALL
        matching_treatment = None
        matches = []

        if text == "":
            # If they've selected a replicate too, then we need a treatment.
            if len(self.treatments) > 1:
                raise ExperimentError("No unique Treatment found.")
            else:
                matches.append(self.treatments[0])
        else:
            look = text.lower()
            len_look = len(text)
            for t in self.treatments:
                current = t.name.lower()
                if len_look <= len(current):
                    if look == current[:len_look]:
                        matches.append(t)

        if not matches:
            raise ExperimentError("No match for {}.".format(text))

        if len(matches) > 1:
            raise Experiment("More than one match for {}.".format(text))

        matching_treatment = matches[0]

        if repnum >= 1:
            matching_replicate = matching_treatment.with_replicate_id(repnum)
        else:
            matching_replicate = None

        return matching_treatment, matching_replicate

    # def visit_replicates(self, visitor,
    #                      filter_treatments=None,
    #                      filter_replicates=None):
    #     try:
    #         for rep, lin in self.iter_lineages():
    #             visitor.visit_lineage(rep, lin)
    #     except (KeyboardInterrupt, SystemExit):
    #         log.info("User interrupted --- quitting")
    #         self.user_interrupt = True
