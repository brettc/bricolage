import logging
log = logging.getLogger('history')

import numpy
import pathlib
import tables
from . import core_ext

class BaseLineage(object):
    def __init__(self, path, params=None, factory_class=None):
        if not isinstance(path, pathlib.Path):
            path = pathlib.Path(path)
        self.path = path
        self.params = params
        self.factory_class = factory_class
        self.world = None
        self.factory = None
        self.population = None
        self.generation = 0

    def assess(self, target):
        assert target.world is self.world
        self.population.assess(target)

    def next_generation(self, mutation_rate, selection_model):
        """Make a new generation"""
        assert selection_model.world is self.world
        self.generation += 1
        self.population.select(selection_model)
        self.population.mutate(mutation_rate, self.generation)

    def _create(self, loading=False):
        assert self.params is not None
        assert self.factory_class is not None
        self.world = core_ext.World(self.params)
        self.factory = self.factory_class(self.world, self.params)
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
        attrs.factory_class = self.factory_class
        attrs.generation = self.generation

        # Save the networks
        z = numpy.zeros
        n = h5.create_table('/', 'networks', z(0, dtype=self.factory.dtype()))
        g = h5.create_table('/', 'generations', z(0, dtype=self._generation_dtype()))

        self._h5 = h5
        self._networks = n
        self._generations = g
        self._attrs = attrs

        self._save_mutable()

    def _open_database(self):
        h5 = tables.open_file(str(self.path), mode='r+')
        attrs = h5.root._v_attrs

        assert attrs.storage_class == self.__class__.__name__

        # Recover the python objects that we need to reconstruct everything
        self.params = attrs.params
        self.factory_class = attrs.factory_class
        self.generation = attrs.generation

        self._h5 = h5
        self._attrs = attrs 
        self._networks = h5.root.networks
        self._generations = h5.root.generations

    def _load(self):
        self._open_database()
        self._create(loading=True)
        # Load the last generation
        g, indexes = self._generations[-1]
        self.generation = g
        arr = self._networks.read_coordinates(indexes)
        self.factory.from_numpy(arr, self.population)
        self._load_mutable()

    def close(self):
        # Avoid messages from tables
        self._h5.close()

    def __del__(self):
        self.close()


class SnapshotLineage(BaseLineage):
    def __init__(self, path, params=None, factory_class=None):
        BaseLineage.__init__(self, path, params, factory_class)
        if params is None:
            self._load()
        else:
            assert factory_class is not None
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
    def __init__(self, path, params=None, factory_class=None, overwrite=False):
        BaseLineage.__init__(self, path, params, factory_class)
        if params is None:
            self._load()
        else:
            assert factory_class is not None
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

    def next_generation(self, mutation_rate, selection_model):
        """Make a new generation"""
        BaseLineage.next_generation(self, mutation_rate, selection_model)
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

