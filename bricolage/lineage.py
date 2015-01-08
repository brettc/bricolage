import logging
log = logging.getLogger('history')

import numpy
import pathlib
import tables
from . import core_ext

class Lineage(object):
    """Wrap a population into a saveable lineage, using pytables"""
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
            assert factory_class is not None
            self.init_factory()
            self.create()

    def __repr__(self):
        return "<Lineage: {}, {}N>".format(
            str(self.path.name),
            # self.factory_class.__name__,
            len(self.networks),
        )

    def init_factory(self):
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
        self._store = self.h5.create_vlarray('/', 'store', tables.ObjectAtom())
        self._store.append(self.params)
        self._store.append(self.factory_class)
        self._store.append(self.current_gen)
        self._gen_record = numpy.zeros(1, dtype=self.generation_dtype())

        self.save_generation(initial=True)

    def load(self):
        # Load H5 table, r+ mode means that it is writeable, but assumes
        # something is already there (unlike 'a')
        self.h5 = tables.open_file(str(self.path), mode='r+')
        self.networks = self.h5.root.networks
        self.generations = self.h5.root.generations
        self._store = self.h5.root.store

        # Recover the python objects that we need to reconstruct everything
        self.params = self._store[0]
        self.factory_class = self._store[1]
        self._gen_record = self.generations[-1:]
        self.current_gen = self._gen_record['generation'][0]
        self.init_factory()

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

    def next_generation(self, mutation_rate, target, selection_model):
        """Make a new generation"""

        # Update the population
        self.current_gen += 1
        self.population.select(target, selection_model)
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

    def get_lineage(self, ident):
        nets = []
        findident = ident
        while findident != -1:
            print findident
            # Do it in reverse order
            found = self.networks[findident]
            nets.append(found)
            findident = found[1]
        return nets[::-1]

    def __del__(self):
        # Avoid messages from tables
        self.h5.close()
