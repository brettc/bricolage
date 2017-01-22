from bricolage.dot_layout import DotMaker
from bricolage.graph_maker import get_graph_by_type, GraphType
from bricolage.core_ext import Population
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
        return "<Replicate: {}-{} Seed-{}>".format(
            self.treatment.name,
            self.seq,
            self.seed,
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

        log.info("Running replicate {}".format(self))
        log.info("In Path {}".format(self.path))
        if not self.path.exists():
            self.path.mkdir(parents=True)
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

            p.mkdir(parents=True)

        # We don't use this, but users of the class might depend on it
        if not self.analysis_path.exists():
            self.analysis_path.mkdir(parents=True)

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
            # We need to create a new lineage. Let's see if the treatment
            # class has some special construction
            if hasattr(self.treatment, 'make_initial_population'):
                def init_pop_fun(factory, size):
                    return self.treatment.make_initial_population(self, factory, size)
            else:
                init_pop_fun = None
            lin = self.treatment.experiment.lineage_class(db_path, self.params,
                                                          init_pop_fun=init_pop_fun)
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

    def draw_winners(self, lin, maxw=1, graph_type=GraphType.GENE, 
                     with_dot=False, knockouts=True):
        winners = [(n.fitness, n.identifier, n) for n in lin.population]
        winners.sort(reverse=True)
        for i, (fit, ident, net) in enumerate(winners):
            if i == maxw:
                break
            self.draw_net('winner', net, target=lin.targets[0], graph_type=graph_type,
                          knockouts=knockouts, with_dot=with_dot)

    def draw_net(self, prefix, net,
                 graph_type=GraphType.GENE_SIGNAL,
                 knockouts=True,
                 target=None,
                 with_dot=False):
        ana = NetworkAnalysis(net)
        ana.annotate(target)
        g = get_graph_by_type(graph_type, ana, knockouts=knockouts)
        d = DotMaker(g)
        p = self.analysis_path / '{}-G{:07d}-N{:02d}-F{}.png'.format(
            prefix, net.generation, net.identifier, net.fitness)
        log.info("Writing {}".format(str(p)))
        d.save_picture(str(p))
        if with_dot:
            d.save_dot(str(p.with_suffix('.dot')))



class Treatment(object):
    """Replicate a set of experiments into different folders"""

    def __init__(self, name, params, count, seed=None):
        self.name = name
        self.params = params
        self.count = count
        self.experiment = None
        self.seq = None
        self.seed = seed
        self.rng = None
        self.replicates = None

    def __repr__(self):
        return "<Treatment: {} R-{} Seed-{}>".format(
            self.name,
            self.count,
            self.seed,
        )

    def _bind(self, experiment):
        assert self.experiment is None
        self.experiment = experiment
        self.seq = experiment.get_next_treatment_seq()
        if self.seed is None:
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
        log.info("Running Treatment {}".format(self))
        if not self.path.exists():
            self.path.mkdir(parents=True)
        if not self.analysis_path.exists():
            self.analysis_path.mkdir(parents=True)

        for r in self.replicates:
            if only_replicate is not None and r is not only_replicate:
                continue
            if not self.experiment.user_interrupt:
                r.run()

    def run_variable(self, replicate, lineage, 
                     targfun, 
                     streak_len,
                     overrun_len,
                     max_gens,
                     log_every=1000,
                     ):

        # assert isinstance(lineage, lineage.FullLineage)
        if len(lineage.targets) == 0:
            lineage.add_target(targfun())

        streak_try = None
        streak_begin = lineage.first_winning_streak(streak_len)

        log.info("First streak found at {}".format(streak_begin))

        while True:
            g = lineage.generation

            if g >= max_gens:
                log.info("Reached max gens")
                if g > max_gens:
                    raise RuntimeError("Generations are > {}".format(max_gens))
                break

            w, b = lineage.population.worst_and_best()

            # Have we gone as far as we want?
            if streak_begin is None:
                # Look for a streak
                if b == 1.0:
                    if streak_try is None:
                        streak_try = g
                else:
                    streak_try = None

                if streak_try and ((g - streak_try + 1) == streak_len):
                    streak_begin = streak_try
                    log.info("Found a streak at {} in generation".format(streak_begin, g))

            # Now see if we can die yet
            if streak_begin and g == (streak_begin + overrun_len):
                log.info("Breaking as streak_begin is {} and generation is {}".format(streak_begin, g))
                break

            if g % log_every == 0:
                replicate.draw_winners(lineage)

            lineage.next_generation()

        replicate.draw_winners(lineage)

class DependentTreatment(Treatment):
    def __init__(self, name, params, count, exp, tnum, repnum, generation=-1, **kw):
        super(DependentTreatment, self).__init__(name, params, count, **kw)
        assert isinstance(exp, Experiment)
        self.starting_rep = exp.get_replicate(tnum, repnum)
        self.starting_generation = generation

    def make_initial_population(self, replicate, factory, size):
        # This is a very ugly way to load a population into a new factory!
        with self.starting_rep.get_lineage() as lin:
            if self.starting_generation > 0:
                other_pop = lin.get_generation(self.starting_generation)
            else:
                other_pop = lin.population

            pop = Population(factory, size)
            arr_old = lin.factory.pop_to_numpy(other_pop)
            arr_new = factory.pop_to_numpy(pop)

            arr_old['generation'] = arr_new['generation']
            arr_old['parent'] = arr_new['parent']
            arr_old['id'] = arr_new['id']

            actual_pop = Population(factory)
            factory.from_numpy(arr_old, actual_pop)

        return actual_pop

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

    def with_treatment_id(self, t_id):
        if t_id > len(self.treatments):
            raise ExperimentError("{} has no treatment with id {}".format(self, t_id))

        tr = self.treatments[t_id - 1]
        assert tr.seq == t_id
        return tr

    def get_replicate(self, t_id, r_id):
        t_idx = t_id - 1
        r_idx = r_id - 1
        if t_id < 1 or t_id > len(self.treatments):
            raise ExperimentError("No treatment {}".format(str(t_id)))

        tr = self.treatments[t_idx]
        assert isinstance(tr, Treatment)
        assert tr.seq == t_id

        if r_id < 1 or r_id > len(tr.replicates):
            raise ExperimentError("No Replicate {}".format(str(r_id)))

        rep = tr.replicates[r_idx]
        assert isinstance(rep, Replicate)
        assert rep.seq == r_id

        return rep

    def get_next_treatment_seq(self):
        self.next_treatment_seq += 1
        return self.next_treatment_seq

    def _create_path(self, path):
        if not path.parent.exists():
            raise ExperimentError("No parent path: {}".format(str(path.parent)))

        log.info("Creating path {}".format(str(path)))
        path.mkdir(parents=True)

    def _remove_paths(self):
        if self.path.exists():
            _remove_folder(self.path)
        if self.path is not self.analysis_path:
            if self.analysis_path.exists():
                _remove_folder(self.analysis_path)
            self.analysis_path.mkdir(parents=True)

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
            raise RuntimeError("Overwrite not currently supported")
            # self._remove_paths()
            # self._construct()

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
                    # Check for list of visitors
                    try:
                        for v in visitor:
                            v.visit_lineage(rep, lin)
                    except TypeError:
                        visitor.visit_lineage(rep, lin)

    def visit_generations(self, visitor, every=1, only=None,
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

                    if only is not None:
                        only = int(only)
                        if only == -1:
                            only = lin.generation

                        if visitor.wants_generation(only):
                            pop = lin.get_generation(only)
                            visitor.visit_generation(only, pop)
                    else:
                        # Now iterate through the generations
                        gen_num = 0
                        last_gen = False
                        while gen_num <= lin.generation:
                            # Check if the visitor wants this generation (as loading is
                            # expensive)
                            if visitor.wants_generation(gen_num):
                                pop = lin.get_generation(gen_num)
                                visitor.visit_generation(gen_num, pop)

                            # Make sure we always do the last generation
                            if gen_num == lin.generation:
                                last_gen = True

                            gen_num += every

                            if gen_num > lin.generation and not last_gen:
                                # We overshot, let's go back
                                gen_num = lin.generation

                    if hasattr(visitor, 'leave_lineage'):
                        visitor.leave_lineage(rep, lin)

    def find_matching(self, treatnum, repnum):
        treatnum = int(treatnum)
        repnum = int(repnum)
        # Default to ALL
        matching_treatment = None

        if treatnum < 1:
            # If they've selected a replicate too, then we need a treatment.
            if repnum >= 1:
                if len(self.treatments) > 1:
                    raise ExperimentError("No treatment supplied and more than one is possible.")
                    # If they've supplied a replicate...
                    
            matching_treatment = self.with_treatment_id(1)
        else:
            matching_treatment = self.with_treatment_id(treatnum)

        if repnum >= 1:
            matching_replicate = matching_treatment.with_replicate_id(repnum)
        else:
            matching_replicate = None

        return matching_treatment, matching_replicate
