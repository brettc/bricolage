import logging
log = logging.getLogger('history')

import numpy
import pathlib
import tables
from . import core_ext

class BaseLineage(object):
    def __init__(self, params=None, factory_class=None):
        self.params = params
        self.factory_class = factory_class
        self.world = None
        self.factory = None
        self.population = None
        self.current_gen = 0

    def _create(self):
        assert self.params is not None
        assert self.factory_class is not None
        self.world = core_ext.World(self.params)
        self.factory = self.factory_class(self.world, self.params)
        self._size = self.params.population_size
        self.population = core_ext.Population(self.factory, self._size)

    def _new_database(self, path):
        # Create H5 table with high compression using blosc
        filters = tables.Filters(complib='blosc', complevel=5)
        h5 = tables.open_file(str(path), 'w', filters=filters)

        # Create the internal registry for pickled python objects
        reg = h5.create_vlarray('/', 'registry', tables.ObjectAtom())
        reg.append(self.params)
        reg.append(self.factory_class)
        reg.append(self.current_gen)

        # Save the networks
        h5.create_table('/', 'networks', numpy.zeros(0, dtype=self.factory.dtype()))
        return h5

    def _load_database(self, path, mode):
        h5 = tables.open_file(str(path), mode=mode)
        reg = h5.root.registry

        # Recover the python objects that we need to reconstruct everything
        self.params = reg[0]
        self.factory_class = reg[1]
        self.current_gen = reg[2]

        return h5


class SnapshotLineage(BaseLineage):
    def __init__(self, params=None, factory_class=None, path=None):
        BaseLineage.__init__(self, params, factory_class)
        if path is None:
            self._create()
        else:
            self._load(path)

    def _load(self, path):
        h5 = self._load_database(path, 'r')
        self._create()
        # Load everything
        self.factory.from_numpy(h5.root.networks[:], self.population)
        h5.close()

    def assess(self, target):
        assert target.world is self.world
        self.population.assess(target)

    def next_generation(self, mutation_rate, selection_model):
        """Make a new generation"""
        assert selection_model.world is self.world
        self.current_gen += 1
        self.population.select(selection_model)
        self.population.mutate(mutation_rate, self.current_gen)

    def save_snapshot(self, path):
        h5 = self._new_database(path)
        arr = self.factory.to_numpy(self.population)
        h5.root.networks.append(arr)
        h5.close()


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
                assert not path.exists
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

    def __del__(self):
        # Avoid messages from tables
        self.h5.close()
