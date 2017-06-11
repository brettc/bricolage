import logtools
import numpy
import pathlib
import tables
import pandas as pd
import numpy as np
from . import core_ext

log = logtools.get_logger()


class LineageError(Exception):
    pass


# Utility class for storing attributes for pickling
class Attributes(object):
    def __init__(self, **kw):
        self.__dict__.update(kw)


class BaseLineage(object):
    def __init__(self, path, params=None, readonly=False):
        if not isinstance(path, pathlib.Path):
            path = pathlib.Path(path)
        self.path = path

        if readonly:
            assert params is None

        self.params = params
        self.readonly = readonly
        self.world = None
        self.factory = None
        self.population = None
        self.targets = []
        self.generation = 0
        self.extra = Attributes()

        # You need to set a target before doing anything
        self.current_target = None
        self.target_index = -1

    def __enter__(self):
        log.debug("Opening lineage {} in __enter__".format(self))
        return self

    def __exit__(self, type, value, traceback):
        log.debug("Closing lineage {} in __exit__".format(self))
        self.close()

    def add_target(self, targ):
        if self.readonly:
            raise LineageError("Network is readonly")

        self.targets.append(targ)

        # By default, set the target to the latest
        self.set_target(-1)
        self._save_header()

    def set_target(self, index=-1):
        if index == -1:
            index = len(self.targets) - 1

        if index == self.target_index:
            return

        assert 0 <= index < len(self.targets)

        self.target_index = index
        self.current_target = self.targets[self.target_index]

        # Make sure our population is assessed by the current target
        self.population.assess(self.current_target)

        # Update the generations record
        self.post_set_target(index)

    def post_set_target(self, index):
        pass

    def next_generation(self):
        """Make a new generation"""
        assert self.current_target is not None
        if self.readonly:
            raise LineageError("Network is readonly")

        # At this stage, the fitnesses MUST have already been calculated when
        # the initial target was set, or from the previous call to
        # next_generation
        self.generation += 1
        self.population.select(self.selection_model)
        n = self.population.mutate(self.params.mutation_rate, 
                                   self.params.trans_mutation_rate,
                                   self.params.duplication_rate,
                                   self.generation)
        log.debug("{} Mutations at generation {}".format(n, self.generation))

        # Now we re-assess the population to ensure that each of them has
        # fitnesses, ready for the next round of selection. Plus, we always
        # want the population to be in a state that everything has a fitness
        self.population.assess(self.current_target)

    def _generation_dtype(self):
        return numpy.dtype([
            ('generation', int),
            ('target', int),
            ('best', float),
            ('indexes', int, self._size),
        ])

    def _commit_generation(self, networks, indexes):
        # Add a generations row
        w, b = self.population.worst_and_best()
        row = numpy.zeros(1, self._generation_dtype())
        row['best'] = b
        row['generation'] = self.generation
        row['target'] = self.target_index
        row['indexes'] = indexes

        # Now commit everything
        try:
            pass
        finally:
            # Force this to always happen together, even if somebady raised a
            # KeyboardException or SystemExit
            #
            # Save the networks, and the generations, and update the header
            self._networks.append(networks)
            self._generations.append(row)
            self._save_header()
            self._h5.flush()

    def _new_database(self):
        # Create H5 table with high compression using blosc
        filters = tables.Filters(complib='blosc', complevel=5)
        h5 = tables.open_file(str(self.path), 'w', filters=filters)
        try:
            attrs = h5.root._v_attrs

            # This is fixed
            attrs.storage_class = self.__class__.__name__

            # Save the networks
            z = numpy.zeros
            n = h5.create_table('/', 'networks',
                                z(0, dtype=self.factory.dtype()))
            g = h5.create_table('/', 'generations',
                                z(0, dtype=self._generation_dtype()))
        except:
            # If things fail, then close the h5 file
            h5.close()
            raise

        self._h5 = h5
        self._networks = n
        self._generations = g
        self._attrs = attrs

        # Save the current set of attributes
        self._save_header()

    def _open_database(self):
        test = tables.is_pytables_file(str(self.path))
        if test <= 0:
            raise LineageError("File is not a pytables")

        if self.readonly:
            mode = 'r'
        else:
            mode = 'r+'

        log.debug("Opening {} in mode {}".format(self.path.name, mode))
        h5 = tables.open_file(str(self.path), mode=mode)
        attrs = h5.root._v_attrs

        assert attrs.storage_class == self.__class__.__name__

        # Assign the arrays ------------
        self._h5 = h5
        self._attrs = attrs
        self._networks = h5.root.networks
        self._generations = h5.root.generations

    def _create(self, init_pop_fun=None):
        assert self.params is not None
        self.world = core_ext.World(self.params)
        self.factory = self.params.factory_class(self.world)

        self._size = self.params.population_size

        if init_pop_fun is None:
            self.population = core_ext.Population(self.factory, self._size)
        else:
            self.population = init_pop_fun(self.factory, self._size)
            assert self.population.size == self._size
            assert self.factory == self.population.factory

        if not hasattr(self.params, 'relative_selection'):
            self.params.relative_selection = False
        self.selection_model = self.params.selection_class(self.world, 
                                                           self.params.relative_selection)

    def _load(self):
        self._load_header()
        if not hasattr(self.params, 'relative_selection'):
            self.params.relative_selection = False

        # Don't bother creating anything, we're about to fill it out
        self.population = core_ext.Population(self.factory, 0)
        self.selection_model = self.params.selection_class(self.world,
                                                           self.params.relative_selection)

        # Load the last generation
        rec = self._generations[-1]
        self.generation = rec['generation']
        indexes = rec['indexes']

        arr = self._networks.read_coordinates(indexes)
        self.factory.from_numpy(arr, self.population)

        target_index = rec['target']
        if target_index >= 0:
            self.set_target(target_index)

        if hasattr(self._attrs, 'extra'):
            self.extra = self._attrs.extra

        log.info("---- Current generation is {}".format(self.generation))

    def _load_header(self):
        # Recover the python objects that we need to reconstruct everything
        # TODO: Fix this so that we look at each of these separately
        data = self._attrs.data
        self.params = data.params
        self.world = data.world
        self.factory = data.factory
        self.targets = data.targets
        self._size = self.params.population_size

    def _save_header(self):
        """Update the header information"""

        # Currently poorly written, as only the world and targets ever get
        # updated
        if self.readonly:
            return

        # Only things that can be pickled go in here
        data = Attributes(
            params=self.params,
            world=self.world,
            factory=self.factory,
            targets=self.targets,
        )
        # This pickles it all
        self._attrs.data = data
        self._attrs.extra = self.extra
        self._h5.flush()

    def close(self):
        raise NotImplementedError


class SnapshotLineage(BaseLineage):
    def __init__(self, path, params=None, readonly=False, init_pop_fun=None):
        BaseLineage.__init__(self, path, params=params, readonly=readonly)
        if params is None:
            self._open_database()
            self._load()
        else:
            self._create(init_pop_fun)
            self._new_database()

    def __repr__(self):
        return "<SnapShotLineage: {}".format(str(self.path.as_posix()))

    def save_snapshot(self):
        # Check how big networks is now...
        start = len(self._networks)

        # Add all networks in the population
        arr = self.factory.pop_to_numpy(self.population)

        # Calculate the indexes for the population
        indexes = numpy.arange(start, start + self._size, dtype=int)
        self._commit_generation(arr, indexes)

    def _load_generation(self, rec):
        gen = rec['generation']
        t_index = rec['target']
        indexes = rec['indexes']
        arr = self._networks.read_coordinates(indexes)
        gen_pop = core_ext.Population(self.factory, 0)
        self.factory.from_numpy(arr, gen_pop)

        # Assess the population by the same target
        gen_pop.assess(self.targets[t_index])
        return gen, gen_pop

    def get_generation(self, wanted_g):
        results = self._generations.read_where(
            'generation == {}'.format(wanted_g))
        if not results:
            return None
        rec = results[0]
        g, gen = self._load_generation(rec)
        assert g == wanted_g
        return gen

    def all_generations(self):
        # Just iterate over what is there
        for rec in self._generations:
            yield self._load_generation(rec)

    def close(self):
        # If the latest generation isn't saved, the save it automatically.
        if not self.readonly:
            g = -1
            if len(self._generations) != 0:
                rec = self._generations[-1]
                g = rec['generation']
            if g < self.generation:
                self.save_snapshot()

            self._save_header()
        self._h5.close()


class FullLineage(BaseLineage):
    def __init__(self, path, params=None, readonly=False, init_pop_fun=None):
        BaseLineage.__init__(self, path, params=params, readonly=readonly)
        if params is None:
            self._open_database()
            self._load()
            self._check_integrity()
        else:
            self._create(init_pop_fun)
            self._new_database()
            self.save_generation(initial=True)

    def __repr__(self):
        return "<FullLineage: {}>".format(str(self.path.as_posix()))

    def post_set_target(self, index):
        if not self.readonly:
            w, b = self.population.worst_and_best()
            self._h5.root.generations.cols.best[self.generation] = b
            self._h5.root.generations.cols.target[self.generation] = index
            self._check_integrity()

    def save_generation(self, initial=False):
        # If it is not the initial save, then just save the mutations only
        arr = self.factory.pop_to_numpy(self.population, not initial)
        # Add a generations row
        idents = numpy.zeros(self._size, int)
        self.population._get_identifiers(idents)
        self._commit_generation(arr, idents)

    def next_generation(self):
        """Make a new generation"""
        BaseLineage.next_generation(self)
        self.save_generation()

    def get_generation(self, wanted_g):
        if wanted_g < 0:
            # Translate negative indices
            actual_g = self.generation + wanted_g + 1
        else:
            actual_g = wanted_g

        assert actual_g <= self.generation
        rec = self._generations[actual_g]
        gen = rec['generation']
        assert gen == actual_g
        t_index = rec['target']
        indexes = rec['indexes']

        arr = self._networks.read_coordinates(indexes)
        gen_pop = core_ext.Population(self.factory, 0)
        self.factory.from_numpy(arr, gen_pop)

        # Assess the population by the same target
        gen_pop.assess(self.targets[t_index])
        return gen_pop

    def all_generations(self, every=1):
        g = 0
        while g <= self.generation:
            yield g, self.get_generation(g)
            g += every

    def get_ancestry(self, ident):
        nets = []
        findident = ident
        while findident != -1:
            # Work our way backwards find the parents
            found = self._networks[findident]
            nets.append(found)
            findident = found[1]

        # Create a numpy array (we could not do that above, as the number was
        # not certain. Note, we do it in reverse order...
        arr = numpy.zeros(len(nets), dtype=self.factory.dtype())
        for i, values in enumerate(nets[::-1]):
            arr[i] = values

        anc = core_ext.Ancestry(self.factory)
        self.factory.from_numpy(arr, anc)
        return anc

    def get_network(self, ident):
        data = self._networks[ident]
        arr = numpy.zeros(1, dtype=self.factory.dtype())
        arr[0] = data
        temp = core_ext.Collection(self.factory)
        self.factory.from_numpy(arr, temp)
        return temp[0]

    def first_winning_generation(self):
        for g in self._h5.root.generations.where("best == 1.0"):
            return g['generation']
        return None

    def first_winning_streak(self, length):
        bests = pd.Series(self._h5.root.generations.cols.best[:])

        # The sum of the window where everything is 1.0
        streaks = np.where((bests.rolling(length, win_type='boxcar').sum()) == length)[0]
        if len(streaks) == 0:
            return None

        # Get the beginning of the streak
        gen = streaks[0] - (length - 1)
        assert bests.values[gen] == 1.0
        if gen != 0:
            assert bests.values[gen - 1] != 1.0
        return gen

    def final_streak_length(self, sample_size):
        best_col = self._h5.root.generations.cols.best

        # Send a -1 if we can't get our sample size
        if sample_size > len(best_col):
            sample_size = len(best_col)

        # Load the last sample_size of best scores
        sample = best_col[-sample_size:]

        # Find the last one that is not equal to 1.0
        found, = np.where(sample != 1.0)
        if len(found) == 0:
            return sample_size

        last_not_1 = found[-1]

        # Get the actual generation of this.
        last_actual_gen = self.generation - sample_size + last_not_1 + 1

        # Check we've found the switch from < 1.0 to 1.0
        assert best_col[last_actual_gen] != 1.0
        # if sample_size > 1:
        #     assert best_col[last_actual_gen + 1] == 1.0

        # Everything after this is 1.0, so the length is...
        return self.generation - last_actual_gen

    def _check_integrity(self):
        num_networks = len(self._networks)

        if len(self._generations) == 0:
            assert num_networks == 0
            assert self.world.next_target_id == 0
        else:
            rec = self._generations[-1]
            max_index = max(rec['indexes'])
            assert max_index < num_networks
            assert num_networks == self.world.next_network_id

    def close(self):
        if not self.readonly:
            self._check_integrity()
            self._h5.flush()

        self._h5.close()


def new_lineage_from_old(db_path, old_lineage, copy_targets=False,
                         kls=SnapshotLineage):

    def init_pop_fun(factory, size):
        pop_old = old_lineage.population
        size_old = pop_old.size
        assert size == size_old

        pop_new = core_ext.Population(factory, size)

        # This does the conversion safely
        arr_old = pop_old.factory.pop_to_numpy(pop_old)
        arr_new = factory.pop_to_numpy(pop_new)

        arr_old['generation'] = arr_new['generation']
        arr_old['parent'] = arr_new['parent']
        arr_old['id'] = arr_new['id']

        # Overwrite the current pop
        factory.from_numpy(arr_old, pop_new)
        return pop_new

    new_lineage = kls(db_path, old_lineage.params, init_pop_fun=init_pop_fun)

    if copy_targets:
        for targ in old_lineage.targets:
            new_lineage.targets.append(targ)
        new_lineage.set_target(-1)

    return new_lineage
