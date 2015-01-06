import logging
log = logging.getLogger('history')

import numpy
import pathlib
import tables
from . import core_ext

class Lineage(object):
    def __init__(self, path, params=None, factory_class=None, overwrite=False):
        if not isinstance(path, pathlib.Path):
            path = pathlib.Path(path)
        self.path = path
        self.params = params
        self.factory_class = factory_class
        self.world = None
        self.factory = None
        self.population = None

        if params is None:
            # We must be loading
            self.load()
        else:
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
        self.size = self.params.population_size
        self.population = core_ext.Population(self.factory, self.size)

    def create(self):
        # Create H5 table with high compression using blosc
        f = tables.Filters(complib='blosc', complevel=5)
        self.h5 = tables.openFile(str(self.path), 'w', filters=f)

        # All the individuals we go through
        initial = numpy.zeros(0, dtype=self.factory.dtype())
        self.networks = self.h5.createTable(
            '/', 'networks', initial)

        self._store = self.h5.createVLArray('/', 'store', tables.ObjectAtom())
        self._store.append(self.params)
        self._store.append(self.factory_class)

        # initial = numpy.zeros(0, dtype=self.parameters.generation_dtype())
        # self.generations = self.h5.createTable(
        #     '/', 'generations', initial)
        
        self.save_generation()

    def load(self):
        # Load H5 table
        self.h5 = tables.openFile(str(self.path), mode='r+')
        self.networks = self.h5.root.networks
        self._store = self.h5.root.store
        self.params = self._store[0]
        self.factory_class = self._store[1]
        self.init_factory()

        # Load the last population
        arr = self.networks[:]
        self.factory.from_numpy(arr, self.population)

    def save_generation(self):
        # assert population.parameters is self.parameters
        arr = self.factory.to_numpy(self.population)
        self.networks.append(arr)

        # grecord = numpy.zeros(1, self.parameters.generation_dtype())
        # grecord['challenge'] = population.challenge_identifier
        # grecord['generation'] = population.generation
        # self.generations.append(grecord)
        self.h5.flush()

    # def create_indexes(self):
    #     log.info('Creating Indexes...')
    #     self.individuals.cols.id.createIndex(kind='light')
    #     self.individuals.cols.parent.createIndex(kind='light')
    #     self.individuals.cols.generation.createIndex(kind='light')
    #     # log.info("Finished Creating Indexes ... %d seconds", self.report_time())
    
    def __del__(self):
        # Avoid messages from tables
        self.h5.close()

    # def close(self):
        # I think it is better to create the indexes at the end of a run ---
        # the simulations seem to slow down significantly if the indexes are
        # there during the run. And we only need them at the end. 
        # if self.mode == 'w':
        #     self.create_indexes()
        #
        # Close the h5 history file
        # self.h5.flush()
        # self.h5.close()
    #
    #
    # # These two functions might be be called multiple times during an
    # # analysis, so let's make it cheap to get the same thing.
    # @lru_cache(maxsize=100)
    # def get_lineage(self, ident):
    #     lastgen = self.individuals[-1]['generation']
    #     netsize = lastgen + 1
    #     nets = numpy.zeros(netsize, self.parameters.network_dtype())
    #     current = -1
    #     findident = ident
    #     while findident != -1:
    #         # Do it in reverse order
    #         nets[current] = self.individuals[findident]
    #         findident = nets[current]['parent']
    #         current -= 1
    #     return nets
    #
    # @lru_cache(maxsize=100)
    # def get_generation(self, g):
    #     if g == -1:
    #         g = self.individuals[-1]['generation']
    #
    #     n = self.parameters.pop_size
    #     f, t = g * n, (g + 1) * n
    #
    #     nets = self.individuals[f:t]
    #     return nets
    #
    # def get_where(self, query):
    #     nets = self.individuals.readWhere(query)
    #     log.debug("Loaded history where '%s', %d results", query, len(nets))
    #     return nets

