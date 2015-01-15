import logging
log = logging.getLogger('history')

import numpy
import pathlib
import tables
import copy
import cPickle as pickle
import shutil
import random
from . import core_ext

class BaseLineage(object):
    def __init__(self, path, params=None):
        if not isinstance(path, pathlib.Path):
            path = pathlib.Path(path)
        self.path = path
        self.params = params
        self.world = None
        self.factory = None
        self.population = None
        self.targets = []
        self.generation = 0

        # Just created or just loaded? This means that the networks have no
        # fitness assigned to them yet (as we don't yet have a target)
        self.fresh = True

    # TODO: target adding
    # def add_target(self, func, target_class=None):
    #     if target_class is None:
    #         target_class = self.params.target_class
    #
    #     # Add locally and save
    #     t = target_class(self.world, func)
    #     self._targets.append(t)
    #     self.targets.append(t)

    def next_generation(self, target):
        """Make a new generation"""
        assert target.world is self.world

        # TODO: THIS IS STILL NOT COOL. We could change targets across loading
        # and this could screw things up. Really need to keep targets in the
        # lineage...?
        #
        # If we're completely fresh, then we need to regenerate the fitnesses,
        # as these are not stored.
        if self.fresh:
            self.population.assess(target)
            self.fresh = False

        self.generation += 1
        self.population.select(self.selection_model)
        self.population.mutate(self.params.mutation_rate, self.generation)

        # Now we re-assess the population to ensure that each of them has
        # fitnesses, ready for the next round of selection, plus the
        # population can also be assessed for fitness characteristics
        self.population.assess(target)

    def _create(self, loading=False):
        assert self.params is not None
        self.world = core_ext.World(self.params)
        self.factory = self.params.factory_class(self.world, self.params)
        self.selection_model = self.params.selection_class(self.world)
        self._size = self.params.population_size

        if loading:
            # Don't bother creating anything, we're about to fill it out
            self.population = core_ext.Population(self.factory, 0)
        else:
            self.population = core_ext.Population(self.factory, self._size)

    def _save_mutable(self):
        self._attrs.random_state = self.world.get_random_state()
        self._attrs.next_network_id = self.world.next_network_id

    def _load_mutable(self):
        self.world.set_random_state(self._attrs.random_state)
        self.world.next_network_id = self._attrs.next_network_id

    def _generation_dtype(self):
        return numpy.dtype([
            ('generation', int),
            ('indexes', int, self._size),
        ])

    def _new_database(self):
        # Create H5 table with high compression using blosc
        filters = tables.Filters(complib='blosc', complevel=5)
        h5 = tables.open_file(str(self.path), 'w', filters=filters)
        attrs = h5.root._v_attrs

        attrs.storage_class = self.__class__.__name__
        attrs.params = self.params
        attrs.generation = self.generation

        # Save the networks
        z = numpy.zeros
        n = h5.create_table('/', 'networks', z(0, dtype=self.factory.dtype()))
        g = h5.create_table('/', 'generations', z(0, dtype=self._generation_dtype()))
        t = h5.create_vlarray('/', 'targets', tables.ObjectAtom())

        self._h5 = h5
        self._networks = n
        self._targets = t
        self._generations = g
        self._attrs = attrs

        self._save_mutable()

    def _open_database(self):
        h5 = tables.open_file(str(self.path), mode='r+')
        attrs = h5.root._v_attrs

        assert attrs.storage_class == self.__class__.__name__

        # Recover the python objects that we need to reconstruct everything
        self.params = attrs.params
        self.generation = attrs.generation

        self._h5 = h5
        self._attrs = attrs 
        self._networks = h5.root.networks
        self._generations = h5.root.generations
        self._targets = h5.root.targets

    def _load(self):
        self._open_database()
        self._create(loading=True)
        # Load the last generation
        g, indexes = self._generations[-1]
        self.generation = g
        arr = self._networks.read_coordinates(indexes)
        self.factory.from_numpy(arr, self.population)
        self._load_mutable()
        # for i in range(len(self._targets)):
        #     self.targets.append(self._targets[i])

    def close(self):
        # Avoid messages from tables
        self._h5.close()

    def __del__(self):
        self.close()


class SnapshotLineage(BaseLineage):
    def __init__(self, path, params=None):
        BaseLineage.__init__(self, path, params)
        if params is None:
            self._load()
        else:
            self._create()
            self._new_database()

    def __repr__(self):
        return "<SnapShotLineage: '{}', {}S, {}N>".format(
            str(self.path.name),
            len(self._generations),
            len(self._networks),
        )

    def save_snapshot(self):
        start = len(self._networks)
        
        # Add all networks in the population
        arr = self.factory.to_numpy(self.population)
        self._networks.append(arr)

        # Add a generations row
        row = self._generations.row
        row['generation'] = self.generation
        row['indexes'] = numpy.arange(start, start + self._size, dtype=int)
        row.append()

        self._save_mutable()
        self._h5.flush()

    def get_generation(self, g):
        grecord = self._generations.read_where('generation == {}'.format(g))
        if not grecord:
            return None
        g, indexes = grecord[0]
        arr = self._networks.read_coordinates(indexes)
        gen_pop = core_ext.Population(self.factory, 0)
        self.factory.from_numpy(arr, gen_pop)
        return gen_pop

    def close(self):
        # If the latest generation isn't saved, the save it automatically.
        g = -1 
        if len(self._generations) != 0:
            g, indexes = self._generations[-1]
        if g < self.generation:
            self.save_snapshot()

        BaseLineage.close(self)

class FullLineage(BaseLineage):
    def __init__(self, path, params=None):
        BaseLineage.__init__(self, path, params)
        if params is None:
            self._load()
        else:
            # if not overwrite:
            #     assert not path.exists()
            self._create()
            self._new_database()
            self.save_generation(initial=True)

    def __repr__(self):
        return "<FullLineage: '{}', {}/{}G, {}N>".format(
            str(self.path.name),
            len(self._generations),
            self._size,
            len(self._networks),
        )

    def save_generation(self, initial=False):
        # If it is not the initial save, then just save the mutations only
        arr = self.factory.to_numpy(self.population, not initial)
        self._networks.append(arr)

        # TODO: make this more efficient?
        idents = numpy.zeros(self._size, int)
        self.population._get_identifiers(idents)

        # Add a generations row
        row = self._generations.row
        row['generation'] = self.generation
        row['indexes'] = idents
        row.append()

        self._h5.flush()

    def next_generation(self, target):
        """Make a new generation"""
        BaseLineage.next_generation(self, target)
        self.save_generation()

    def get_generation(self, wanted_g):
        assert wanted_g <= self.generation
        actual_g, identifiers = self._generations[wanted_g]
        assert actual_g == wanted_g
        arr = self._networks.read_coordinates(identifiers)
        gen_pop = core_ext.Population(self.factory, 0)
        self.factory.from_numpy(arr, gen_pop)
        return gen_pop

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

    def close(self):
        self._save_mutable()
        self._h5.close()

    def __del__(self):
        self.close()


class Replicate(object):
    """Wraps a lineage into a sub-folder"""

    LINEAGE_FILENAME = 'lineage.h5'

    def __init__(self, treatment, seq):
        self.treatment = treatment
        self.sequence = seq
        self.seed = treatment.rand.randint(0, 1 << 16)

    @property
    def name(self):
        return "{:03d}".format(self.sequence)

    @property
    def path(self):
        return self.treatment.path / self.name

    def get_lineage(self):
        p = self.path
        if not p.exists():
            p.mkdir()

        db_path = p / self.LINEAGE_FILENAME

        # Load the lineage
        if db_path.exists():
            lin = FullLineage(db_path)
            # If we reloaded then the seed should be the same as what was
            # generated by the original parameters
            assert lin.params.seed == self.seed
        else:
            # Make a copy of the parameters, but update the seed. This is the
            # only thing that changes in a replicate! Take a deepcopy so that
            # we don't modify anything from others.
            new_params = copy.deepcopy(self.treatment.params)
            new_params.seed = self.seed
            lin = FullLineage(db_path, new_params)

        return lin


class Treatment(object):
    """Replicate a set of experiments into different folders"""

    TREATMENT_FILENAME = 'treatment.pickle'

    def __init__(self, path, params=None, overwrite=False):
        if not isinstance(path, pathlib.Path):
            path = pathlib.Path(path)
        self.path = path
        self.params = params
        if overwrite or not self.path.exists():
            self._new()
        else:
            self._load()
        self._create()

    @property
    def filename(self):
        return str(self.path / self.TREATMENT_FILENAME)

    def run(self, callback):
        for r in self.replicates:
            callback(r)

    def _create(self):
        self.rand = random.Random(self.params.seed)
        self.replicates = [Replicate(self, i) for i in range(self.params.replicates)]

    def _new(self):
        if self.path.exists():
            shutil.rmtree(str(self.path))
        self.path.mkdir()
        with open(self.filename, 'wb') as f:
            pickle.dump(self.params, f, protocol=pickle.HIGHEST_PROTOCOL)

    def _load(self):
        with open(self.filename, 'rb') as f:
            self.params = pickle.load(f)



