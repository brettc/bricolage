import logging

log = logging.getLogger('lineage')

import numpy
import pathlib
import tables
from . import core_ext


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

        # You need to set a target before doing anything
        self.current_target = None
        self.target_index = -1

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def add_target(self, func, name="", target_class=None):
        if self.readonly:
            raise LineageError("Network is readonly")
        if target_class is None:
            target_class = self.params.target_class

        # Make up a name if none is given
        if name == "":
            name = str(len(self.targets))

        # Add locally and save
        t = target_class(self.world, func, name=name)
        self.targets.append(t)

        # By default, set the target to the latest
        self.set_target(len(self.targets) - 1)
        self._save()

    # TODO: add by name?
    def set_target(self, index):
        assert 0 <= index < len(self.targets)
        self.target_index = index
        self.current_target = self.targets[self.target_index]

        # Make sure our population is assessed by the current target
        self.population.assess(self.current_target)

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
        self.population.mutate(self.params.mutation_rate, self.generation)

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

    def _add_generation(self, indexes):
        # Add a generations row
        row = self._generations.row

        w, b = self.population.worst_and_best()
        row['best'] = b
        row['generation'] = self.generation
        row['target'] = self.target_index
        row['indexes'] = indexes
        row.append()

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
        self._save()

    def _open_database(self):
        test = tables.is_pytables_file(str(self.path))
        if test <= 0:
            raise LineageError("File is not a pytables")

        if self.readonly:
            mode = 'r'
        else:
            mode = 'r+'

        h5 = tables.open_file(str(self.path), mode=mode)
        attrs = h5.root._v_attrs

        assert attrs.storage_class == self.__class__.__name__

        # Assign the arrays
        self._h5 = h5
        self._attrs = attrs
        self._networks = h5.root.networks
        self._generations = h5.root.generations

    def _create(self):
        assert self.params is not None
        self.world = core_ext.World(self.params)
        self.factory = self.params.factory_class(self.world)

        self._size = self.params.population_size
        self.population = core_ext.Population(self.factory, self._size)
        self.selection_model = self.params.selection_class(self.world)

    def _load(self):
        # Recover the python objects that we need to reconstruct everything
        try:
            data = self._attrs.data

            self.params = data.params
            self.world = data.world
            self.factory = data.factory
            self.targets = data.targets
            self._size = self.params.population_size

            # Don't bother creating anything, we're about to fill it out
            self.population = core_ext.Population(self.factory, 0)
            self.selection_model = self.params.selection_class(self.world)

            # Load the last generation
            rec = self._generations[-1]
            self.generation = rec['generation']
            indexes = rec['indexes']

            arr = self._networks.read_coordinates(indexes)
            self.factory.from_numpy(arr, self.population)
            target_index = rec['target']
            if target_index >= 0:
                self.set_target(target_index)
        except:
            self._h5.close()
            raise

    def _save(self):
        if self.readonly:
            return
            # raise LineageError("Network is readonly")

        # Only things that can be pickled go in here
        data = Attributes(
            params=self.params,
            world=self.world,
            factory=self.factory,
            targets=self.targets,
        )
        # This pickles it all
        self._attrs.data = data
        self._h5.flush()

    def close(self):
        raise NotImplementedError


class SnapshotLineage(BaseLineage):
    def __init__(self, path, params=None, readonly=False):
        BaseLineage.__init__(self, path, params=params, readonly=readonly)
        if params is None:
            self._open_database()
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
        # TODO: is this needed
        self._save()

        # Check how big networks is now...
        start = len(self._networks)

        # Add all networks in the population
        arr = self.factory.pop_to_numpy(self.population)
        self._networks.append(arr)

        # Add a generations row
        indexes = numpy.arange(start, start + self._size, dtype=int)
        self._add_generation(indexes)

        self._h5.flush()

    def get_generation(self, wanted_g):
        results = self._generations.read_where(
            'generation == {}'.format(wanted_g))
        if not results:
            return None
        rec = results[0]
        gen = rec['generation']
        assert gen == wanted_g
        t_index = rec['target']
        indexes = rec['indexes']
        arr = self._networks.read_coordinates(indexes)
        gen_pop = core_ext.Population(self.factory, 0)
        self.factory.from_numpy(arr, gen_pop)

        # Assess the population by the same target 
        gen_pop.assess(self.targets[t_index])
        return gen_pop

    def close(self):
        # If the latest generation isn't saved, the save it automatically.
        if not self.readonly:
            g = -1
            if len(self._generations) != 0:
                rec = self._generations[-1]
                g = rec['generation']
            if g < self.generation:
                self.save_snapshot()

            self._save()
        self._h5.close()


class FullLineage(BaseLineage):
    def __init__(self, path, params=None, readonly=False):
        BaseLineage.__init__(self, path, params=params, readonly=readonly)
        if params is None:
            self._open_database()
            self._load()
        else:
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
        arr = self.factory.pop_to_numpy(self.population, not initial)
        self._networks.append(arr)

        # TODO: make this more efficient?
        idents = numpy.zeros(self._size, int)
        self.population._get_identifiers(idents)

        # Add a generations row
        self._add_generation(idents)
        self._h5.flush()

    def next_generation(self):
        """Make a new generation"""
        BaseLineage.next_generation(self)
        self.save_generation()

    def get_generation(self, wanted_g):
        assert wanted_g <= self.generation
        rec = self._generations[wanted_g]
        gen = rec['generation']
        assert gen == wanted_g
        t_index = rec['target']
        indexes = rec['indexes']
        arr = self._networks.read_coordinates(indexes)
        gen_pop = core_ext.Population(self.factory, 0)
        self.factory.from_numpy(arr, gen_pop)

        # Assess the population by the same target 
        gen_pop.assess(self.targets[t_index])
        return gen_pop

    def all_generations(self):
        g = 0
        while g <= self.generation:
            yield g, self.get_generation(g)
            g += 1

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
        # Avoid messages from tables
        if not self.readonly:
            self._save()
        self._h5.close()


