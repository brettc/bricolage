import logging

log = logging.getLogger('lineage')

import numpy
import pathlib
import tables
import copy
import cPickle as pickle
import shutil
import random
import yaml
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
            self.set_target(rec['target'])
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


class Replicate(object):
    """Wraps a lineage into a sub-folder"""

    LINEAGE_FILENAME = 'lineage.h5'

    def __init__(self, treatment, seq):
        self.treatment = treatment
        self.sequence = seq
        self.seed = treatment.rand.randint(0, 1 << 16)
        self.fresh = True

    @property
    def name(self):
        return "{:03d}".format(self.sequence)

    def __repr__(self):
        return "<Replicate:{}>".format(self.name)

    @property
    def path(self):
        return self.treatment.path / self.name

    @property
    def analysis_path(self):
        return self.treatment.analysis_path / self.name

    def get_lineage(self, readonly=False):
        p = self.path
        if not p.exists():
            if readonly:
                log.error("Replicate doesn't exist, {}".format(self.name))
                raise LineageError

            p.mkdir()

        # We don't use this, but users of the class might depend on it
        if not self.analysis_path.exists():
            self.analysis_path.mkdir()

        db_path = p / self.LINEAGE_FILENAME

        # Load the lineage
        if db_path.exists():
            try:
                lin = self.treatment.lineage_class(db_path, readonly=readonly)
                # If we reloaded then the seed should be the same as what was
                # generated by the original parameters
                assert lin.params.seed == self.seed

                # Kill this flag
                self.fresh = False
            except ValueError:
                log.error("Error opening existing lineage."
                          "I'm going to create a new one.")

        if self.fresh:
            # Make a copy of the parameters, but update the seed. This is the
            # only thing that changes in a replicate! Take a deepcopy so that
            # we don't modify anything from others.
            new_params = copy.deepcopy(self.treatment.params)
            new_params.seed = self.seed
            lin = self.treatment.lineage_class(db_path, new_params)

        return lin


class StopReplicates(Exception):
    pass


class Treatment(object):
    """Replicate a set of experiments into different folders"""

    TREATMENT_FILENAME = 'treatment.pickle'
    DESCRIPTION_FILENAME = 'description.txt'

    def __init__(self, path, params=None, analysis_path=None, overwrite=False,
                 full=True, readonly=False):
        # Sort out class to use
        if full:
            self.lineage_class = FullLineage
        else:
            self.lineage_class = SnapshotLineage

        # Set up paths
        if not isinstance(path, pathlib.Path):
            path = pathlib.Path(path)
        assert path.parent.exists()
        self.path = path

        if analysis_path is None:
            self.analysis_path = self.path
        else:
            if not isinstance(analysis_path, pathlib.Path):
                analysis_path = pathlib.Path(analysis_path)
            self.analysis_path = analysis_path
        assert self.analysis_path.parent.exists()

        # Now establish how we load
        self.params = params
        if not readonly and (overwrite or not self.path.exists()):
            self._new()
        else:
            if not self.exists():
                log.error("No such path {}".format(self.filename))
                raise LineageError

            self._load()

            if self.params is None:
                log.error("Failed to load {}".format(self.filename))
                raise LineageError

        self._create()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def close(self):
        # TODO: Check that everything has closed?
        pass

    @property
    def filepath(self):
        return self.path / self.TREATMENT_FILENAME

    @property
    def filename(self):
        return str(self.filepath)

    def exists(self):
        return (self.path / self.TREATMENT_FILENAME).exists()

    @property
    def name(self):
        try:
            return self.params.name
        except AttributeError:
            return self.path.parts[-1]

    def __repr__(self):
        return "<Treatment:{}>".format(self.name)

    @property
    def description_filename(self):
        return str(self.path / self.DESCRIPTION_FILENAME)

    def run(self, callback, **kwargs):
        for r in self.replicates:
            try:
                callback(r, **kwargs)
            except StopReplicates:
                break

    def iter_replicates(self):
        for r in self.replicates:
            yield r

    def _create(self):
        self.rand = random.Random(self.params.seed)
        self.replicates = [Replicate(self, i) for i in
                           range(self.params.replicates)]

    def _new(self):
        if self.path.exists():
            shutil.rmtree(str(self.path))
        self.path.mkdir()
        if self.path is not self.analysis_path:
            if self.analysis_path.exists():
                shutil.rmtree(str(self.analysis_path))
            self.analysis_path.mkdir()
        with open(self.filename, 'wb') as f:
            pickle.dump(self.params, f, protocol=pickle.HIGHEST_PROTOCOL)

        with open(self.description_filename, 'wb') as f:
            f.write(yaml.dump(self.params.__dict__, default_flow_style=False))

    def _load(self):
        try:
            f = open(self.filename, 'rb')
            self.params = pickle.load(f)
        except:
            log.error("Error loading treatment file {}".format(self.filename))
            raise LineageError
