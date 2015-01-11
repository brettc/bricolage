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

        # In addition, if we're loading, we need to set a few mutable thing
        # from the world
        if loading:
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

        # Recover the python objects that we need to reconstruct everything
        self.params = attrs.params
        self.factory_class = attrs.factory_class
        self.generation = attrs.generation

        self._h5 = h5
        self._attrs = attrs 
        self._networks = h5.root.networks
        self._generations = h5.root.generations

    def __del__(self):
        # Avoid messages from tables
        del self._attrs
        self._h5.close()


class SnapshotLineage(BaseLineage):
    def __init__(self, path, params=None, factory_class=None):
        BaseLineage.__init__(self, path, params, factory_class)
        if params is None:
            self._load()
        else:
            assert factory_class is not None
            self._create()
            self._new_database()

    def _load(self):
        self._open_database()
        self._create(loading=True)
        # Load the last generation
        g, indexes = self._generations[-1]
        self.generation = g
        arr = self._networks.read_coordinates(indexes)
        self.factory.from_numpy(arr, self.population)
        self._load_mutable()

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

class Lineage(object):
    def __init__(self, path, params=None, factory_class=None, overwrite=False):
        if not isinstance(path, pathlib.Path):
            path = pathlib.Path(path)

        # TODO: make some of this "private"
        self.path = path
        self.params = params
        self.factory_class = factory_class
        self.world = None
        self.factory = None
        self.population = None
        self.current_gen = 0

        if params is None:
            # We must be loading
            self.load()
        else:
            # We're making a new one
            assert hasattr(params, 'population_size')
            if not overwrite:
                assert not path.exists()
            assert factory_class is not None
            self._create()
            self.create()

    def __repr__(self):
        return "<Lineage: {}, {}N>".format(
            str(self.path.name),
            # self.factory_class.__name__,
            len(self.networks),
        )

    def _create(self):
        self.world = core_ext.World(self.params)
        self.factory = self.factory_class(self.world, self.params)
        self._size = self.params.population_size
        self.population = core_ext.Population(self.factory, self._size)

    def generation_dtype(self):
        return numpy.dtype([
            ('generation', int),
            ('networks', int, self._size),
        ])

    def create(self):
        # Create H5 table with high compression using blosc
        f = tables.Filters(complib='blosc', complevel=5)
        self.h5 = tables.open_file(str(self.path), 'w', filters=f)

        # All the individuals we go through
        initial = numpy.zeros(0, dtype=self.factory.dtype())
        self.networks = self.h5.create_table('/', 'networks', initial)
        gens = numpy.zeros(0, dtype=self.generation_dtype())
        self.generations = self.h5.create_table('/', 'generations', gens)

        # Create an array that can store pickled items. We keep these VERY
        # simple, only storing the initial parameters and the factory class.
        # We rely on reconstructing everything.
        self._registry = self.h5.create_vlarray('/', 'registry', tables.ObjectAtom())
        self._registry.append(self.params)
        self._registry.append(self.factory_class)
        self._registry.append(self.current_gen)
        self._gen_record = numpy.zeros(1, dtype=self.generation_dtype())

        self.save_generation(initial=True)

    def load(self):
        # Load H5 table, r+ mode means that it is writeable, but assumes
        # something is already there (unlike 'a')
        self.h5 = tables.open_file(str(self.path), mode='r+')
        self.networks = self.h5.root.networks
        self.generations = self.h5.root.generations
        self._registry = self.h5.root.store

        # Recover the python objects that we need to reconstruct everything
        self.params = self._registry[0]
        self.factory_class = self._registry[1]
        self._gen_record = self.generations[-1:]
        self.current_gen = self._gen_record['generation'][0]
        self._create()

        # Load the last population, using the coordinates that were written
        # into the last generation 
        arr = self.networks.read_coordinates(self._gen_record['networks'][0])
        self.factory.from_numpy(arr, self.population)

    def save_generation(self, initial=False):
        # assert population.parameters is self.parameters
        # If it is not the initial save, then just save the mutations only
        arr = self.factory.to_numpy(self.population, not initial)
        self.networks.append(arr)
        self.generations.append(self._gen_record)

        self.h5.flush()

    def assess(self, target):
        self.population.assess(target)

    def next_generation(self, mutation_rate, selection_model):
        """Make a new generation"""

        # Update the population
        self.current_gen += 1
        self.population.select(selection_model)
        self.population.mutate(mutation_rate, self.current_gen)

        # Update the generation record, let the population fill out the
        # identifiers
        self._gen_record['generation'][0] = self.current_gen
        idents = self._gen_record['networks'][0]
        self.population._get_identifiers(idents)

        self.save_generation()

    def get_generation(self, wanted_g):
        assert wanted_g <= self.current_gen
        actual_g, identifiers = self.generations[wanted_g]
        assert actual_g == wanted_g
        arr = self.networks.read_coordinates(identifiers)
        gen_pop = core_ext.Population(self.factory, 0)
        self.factory.from_numpy(arr, gen_pop)
        return gen_pop

    def get_ancestry(self, ident):
        nets = []
        findident = ident
        while findident != -1:
            # Work our way backwards find the parents
            found = self.networks[findident]
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

