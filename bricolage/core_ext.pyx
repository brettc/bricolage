# distutils: include_dirs = NUMPY_PATH
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# distutils: define_macros=NPY_NO_DEPRECATED_API

import cython
import numpy
import copy
import sys
# from operand import Operand

# cimports
cimport numpy as np
import numpy
from cython.operator import dereference as deref, preincrement as preinc

import random

reserved_channels = 2
off_channel = 0
on_channel = 1

cdef class ChannelStateFrozen:
    def __cinit__(self):
        pass

    cdef init(self, cWorld_ptr &w, cChannelState &p):
        self.world = w
        self._this = p

    def test(self, size_t i):
        return self._this.test(i)

    def __getitem__(self, size_t i):
        return self._this.test(i)

    def __cmp__(self, ChannelStateFrozen other):
        return bitset_cmp(self._this, other._this)

    property size:
        def __get__(self):
            return self._this.size()

    def __str__(self):
        # TODO: Clean this function up---it is really crappy
        cdef:
            string cstr
            cWorld *w = self.world.get()

        if w == NULL:
            # TODO: Some proper exceptions would be good...
            raise RuntimeError

        cdef size_t cuereg = w.cue_channels + w.reg_channels

        to_string(self._this, cstr)

        # I think it is much easier to understand if we reverse it
        # Also, clip the reserved channels 
        # TODO: fix this nasty hack
        s = cstr[::-1][2:]
        env = s[:w.cue_channels]
        reg = s[w.cue_channels:cuereg]
        out = s[cuereg:]

        return env + '|' + reg + '|' + out

    def __copy__(self):
        other = ChannelState()
        other.init(self.world, self._this)
        return other

    def copy(self):
        return self.__copy__()

    def as_array(self):
        vals = numpy.zeros(self.size, dtype=numpy.int32)
        cdef: 
            np.npy_int32[:] v = vals
            size_t i

        for i in range(self._this.size()):
            v[i] = self._this.test(i)

        return vals

    def __repr__(self):
        return "<ChannelsRO: {}>".format(self.__str__())


cdef class ChannelState(ChannelStateFrozen):

    def __cinit__(self):
        pass

    def set(self, size_t i):
        self._this.set(i)

    def reset(self, size_t i):
        self._this.reset(i)

    def flip(self, size_t i):
        self._this.flip(i)

    def __setitem__(self, size_t i, bint b):
        if b:
            self._this.set(i)
        else:
            self._this.reset(i)

    def merge(self, ChannelStateFrozen other):
        self._this |= other._this

    def __repr__(self):
        return "<Channels: {}>".format(self.__str__())


def _construct_world(params, net_id, target_id, r_state):
    w = World(params)
    w._this.next_network_identifier = net_id
    w._this.next_target_identifier = target_id
    w._this.set_random_state(r_state)
    return w


cdef class World:
    def __cinit__(self, params):
        self._params = copy.deepcopy(params)
        self._shared = cWorld_ptr(new cWorld(
            params.seed, 
            params.cue_channels, 
            params.reg_channels, 
            params.out_channels,
        ))
        self._this = self._shared.get()

        if hasattr(params, 'input_type'):
            self._this.input_type = params.input_type

        self.reserved_signals = set([on_channel, off_channel])
        self.cue_signals = set(range(*self._this.cue_range))
        self.reg_signals = set(range(*self._this.reg_range))
        self.out_signals = set(range(*self._this.out_range))
        self.sub_signals = set(range(*self._this.sub_range))
        self.pub_signals = set(range(*self._this.pub_range))

    def __reduce__(self):
        return _construct_world, (self._params, 
                                  self._this.next_network_identifier,
                                  self._this.next_target_identifier,
                                  self._this.get_random_state())
    property params:
        def __get__(self):
            return copy.deepcopy(self._params)

    def create_state(self):
        c = ChannelState()
        c._this.resize(self._this.channel_count)
        c.world = self._shared
        return c

    property environments:
        def __get__(self):
            if self._environments is not None:
                return self._environments

            cdef vector[cChannelState].iterator i =  self._this.environments.begin()
            envs = []
            while i != self._this.environments.end():
                p = ChannelStateFrozen()
                p.init(self._shared, deref(i))
                envs.append(p)
                preinc(i)

            self._environments = envs
            return envs

    # Some routines just for testing, mainly for testing (NOT fast)
    def seed_random_engine(self, int s):
        self._this.rand.seed(s)

    def get_random_state(self):
        return self._this.get_random_state()

    # TODO: Remove me
    def set_random_state(self, string s):
        self._this.set_random_state(s)

    def get_random_double(self, double low, double high):
        return self._this.get_random_double(low, high)

    def get_random_int(self, int low, int high):
        return self._this.get_random_int(low, high)

    property next_network_id:
        def __get__(self):
            return self._this.next_network_identifier
        # def __set__(self, sequence_t i):
        #     self._this.next_network_identifier = i

    property next_target_id:
        def __get__(self):
            return self._this.next_target_identifier
        # def __set__(self, sequence_t i):
        #     self._this.next_network_identifier = i

    property cue_channels:
        def __get__(self):
            return self._this.cue_channels
    property reg_channels:
        def __get__(self):
            return self._this.reg_channels
    property out_channels:
        def __get__(self):
            return self._this.out_channels
    property channel_count:
        def __get__(self):
            return self._this.channel_count
    # property off_channel:
    #     def __get__(self):
    #         return off_channel
    # property on_channel:
    #     def __get__(self):
    #         return off_channel

    property cue_range:
        def __get__(self):
            return self._this.cue_range
    property reg_range:
        def __get__(self):
            return self._this.reg_range
    property out_range:
        def __get__(self):
            return self._this.out_range
    property sub_range:
        def __get__(self):
            return self._this.sub_range
    property pub_range:
        def __get__(self):
            return self._this.pub_range
    
    def name_for_channel(self, c):
        sz = 1
        if c == off_channel:
            return "OFF"
        if c == on_channel:
            return "ON"
        # We use "mathematical" number, starting at 1.0
        if c in self.out_signals:
            return "P{0:0{1:}d}".format(
                c + 1 - self._this.out_range.first, sz)
        if c in self.reg_signals:
            return "T{0:0{1:}d}".format(
                c + 1 - self._this.reg_range.first, sz)
        return "E{0:0{1:}d}".format(
            c + 1 - self._this.cue_range.first, sz)

    property input_type:
        def __get__(self):
            return self._this.input_type


cdef class Factory:
    def __cinit__(self, World w):
        self.world = w
        self.gene_class = Gene
        self.module_class = CisModule
        self.network_class = Network
        self._secret_key = 0x8008008

    def create_network(self):
        cdef Network n = self.network_class(self, self._secret_key)
        n.bind_to(self._this.construct(True))
        return n


def _construct_network(Factory factory, array):
    c = Collection(factory)
    factory.from_numpy(array, c)
    return c[0]


cdef class Network:
    def __cinit__(self, Factory c, int key=0):
        self.factory = c
        self._this = NULL
        if key != self.factory._secret_key:
            raise RuntimeError("You cannot create a Network directly."
                               " Use Factory.create_network")

    def __reduce__(self):
        c = Collection(self.factory)
        c.add(self)
        return _construct_network, (self.factory, 
                                    self.factory.to_numpy(c))

    cdef bind_to(self, cNetwork_ptr ptr):
        self._shared = ptr
        self._this = self._shared.get()
        self._this.pyobject = <void *>(self)

    def __dealloc__(self):
        if self._this != NULL:
            self._this.pyobject = <void *>0

    property identifier:
        def __get__(self):
            return self._this.identifier

    property parent_identifier:
        def __get__(self):
            return self._this.parent_identifier

    property generation:
        def __get__(self):
            return self._this.generation

    property gene_count:
        def __get__(self):
            return self._this.gene_count()

    property genes:
        def __get__(self):
            if self._genes is None:
                self._genes = tuple(self.factory.gene_class(self, i) for i in range(self.gene_count))
            return self._genes

    property attractors_size:
        def __get__(self):
            r = numpy.zeros(self._this.attractors.size())
            cdef:
                vector[cChannelStateVector].iterator cattr_iter
                np.npy_double[:] c_r = r
                size_t i = 0

            cattr_iter = self._this.attractors.begin()
            while cattr_iter != self._this.attractors.end():
                c_r[i] = deref(cattr_iter).size()
                i += 1
                preinc(cattr_iter)

            return r

    def cycle(self, ChannelState c):
        self._this.cycle(c._this)

    def cycle_with_intervention(self, ChannelState c):
        self._this.cycle_with_intervention(c._this)

    def calc_perturbation(self):
        self._this.calc_perturbation()

    def recalculate(self, with_intervention=False):
        if with_intervention:
            self._this.calc_attractors_with_intervention()
        else:
            self._this.calc_attractors()

        self._attractors = None
        self._rates = None

    def mutate(self, size_t nmutations):
        """Mutate the network. 

        This invalidates lots of assumptions required for selection to work
        correctly. Use only if you understand what you are doing.
        """
        self._this.mutate(nmutations)
        self.recalculate()

    cdef _make_python_attractors(self, cAttractors &attrs):
        cdef:
            vector[cChannelStateVector].iterator cattr_iter
            vector[cChannelState].iterator cstate_iter

        py_attrs = []
        cattr_iter = attrs.begin()
        while cattr_iter != attrs.end():
            attr = []
            cstate_iter = deref(cattr_iter).begin()
            while cstate_iter != deref(cattr_iter).end():
                c = ChannelStateFrozen()
                c.init(self.factory._this.world, deref(cstate_iter))
                attr.append(c)
                preinc(cstate_iter)

            # Make everything a tuple -- you shouldn't be able to mess
            # with it!
            py_attrs.append(tuple(attr))
            preinc(cattr_iter)

        return tuple(py_attrs)

    cdef _make_python_rates(self, cRatesVector &rates):
        cdef cWorld *world = self.factory._this.world.get()

        # Construct the numpy array via python
        r = numpy.zeros((world.environments.size(), world.out_channels))

        cdef:
            np.npy_double[:,:] c_r = r
            size_t i, j

        # Copy in the goods using a memory array
        for i in range(world.environments.size()):
            for j in range(world.out_channels):
                c_r[i, j] = rates[i][j]

        # Don't mess with it!
        r.flags.writeable = False
        return r

    property attractors:
        """A tuple containing the attractors for each environment"""
        def __get__(self):
            # Have if we're detached, we have to recalc every time
            if self._attractors is not None:
                return self._attractors

            self._attractors = self._make_python_attractors(self._this.attractors)
            return self._attractors

    property pert_attractors:
        def __get__(self):
            return self._make_python_attractors(self._this.pert_attractors)

    property rates:
        """Return a readonly numpy array of the rates"""
        def __get__(self):
            # Lazy evaluation
            if self._rates is not None:
                return self._rates
            self._rates = self._make_python_rates(self._this.rates)
            return self._rates

    property pert_rates:
        def __get__(self):
            return self._make_python_rates(self._this.pert_rates)

    property fitness:
        def __get__(self):
            return self._this.fitness

    property target:
        def __get__(self):
            return self._this.target

    def __repr__(self):
        return "<Network id:{} pt:{}>".format(self.identifier, self.parent_identifier)


cdef class Gene:
    """A proxy to a gene.
    """
    def __cinit__(self, Network n, size_t g):
        """This is not meant to be called publicly"""
        self.network = n
        self.gene_number = g
        self._this = self.network._this.get_gene(g)

    property sequence:
        def __get__(self):
            return self._this.sequence

    property pub:
        def __get__(self):
            return self._this.pub

    property intervene:
        def __get__(self):
            return self._this.intervene
        def __set__(self, InterventionState i):
            self._this.intervene = i
            self.network.recalculate(with_intervention=True)

    property module_count:
        def __get__(self):
            return self._this.module_count()

    property modules:
        def __get__(self):
            # Lazy construction
            if self._modules is None:
                self._modules = tuple(self.network.factory.module_class(self, i) 
                    for i in range(self.module_count))
            return self._modules

    def __repr__(self):
        w = self.network.factory.world
        return "<Gene[{}]: {}>".format(self.sequence, w.name_for_channel(self.pub))


cdef class CisModule:
    """A proxy to a CisModule.
    """
    def __cinit__(self, Gene g, size_t i):
        self.gene = g
        assert i < g._this.module_count()
        self._this = g._this.get_module(i)

    def site_count(self):
        return self._this.site_count()

    def get_site(self, size_t i):
        assert i < self._this.site_count()
        return self._this.get_site(i)

    def set_site(self, size_t i, size_t c):
        assert c in self.gene.network.world.sub_signals
        old = self._this.set_site(i, c)
        self.gene.network.recalculate()
        return old

    property intervene:
        def __get__(self):
            return self._this.intervene
        # def __set__(self):
        #     return self._this.intervene

    property channels:
        def __get__(self):
            return tuple(self._this.channels[i]
                         for i in range(self._this.site_count()))
        def __set__(self, t):
            assert len(t) == self._this.site_count()
            valid = self.gene.network.factory.world.sub_signals
            cdef size_t i, c
            for i in range(self._this.site_count()):
                assert t[i] in valid
                self._this.channels[i] = t[i]
            self.gene.network.recalculate()

    property channel_names:
        def __get__(self):
            w = self.gene.network.factory.world
            return [w.name_for_channel(c) for c in self.channels]


cdef class SelectionModel:
    def __cinit__(self, World w):
        self.world = w
        self._this = new cSelectionModel(w._shared)

    def __dealloc__(self):
        del self._this


cdef class CollectionBase:
    def __cinit__(self, Factory c, size_t size=0):
        self.factory = c

    cdef object get_at(self, size_t i):
        if i >= self._collection.size():
            raise IndexError

        cdef:
            cNetwork_ptr ptr = deref(self._collection)[i]
            cNetwork *net = ptr.get()

        # Is there an existing python object?
        # Note: cython automatically increments the reference count when we do
        # this (which is what we want). This move just means we get object
        # identity, at least while one python reference continues to exist.
        if net.pyobject:
            return <object>(net.pyobject)

        # Nope. We need to create a new python wrapper object
        cdef Network n = self.factory.network_class(self.factory, 
                                                        self.factory._secret_key)
        n.bind_to(ptr)
        return n

    def add(self, Network n):
        assert n.factory is self.factory
        self._collection.push_back(n._shared)

    property size:
        def __get__(self):
            return self._collection.size()

    def __getitem__(self, size_t i):
        return self.get_at(i)

    def __iter__(self):
        for i in range(self._collection.size()):
            yield self.get_at(i)

    def assess(self, BaseTarget target):
        target._base.assess_networks(deref(self._collection))

    def get_fitnesses(self, np.double_t[:] fits):
        assert fits.shape[0] == self.size
        cdef size_t i
        for i in range(self._collection.size()):
            fits[i] = deref(self._collection)[i].get().fitness

    def extend(self, CollectionBase other):
        for i in range(other._collection.size()):
            self._collection.push_back(deref(other._collection)[i])

    def fill(self, Network n, size_t size):
        cdef size_t i
        for i in range(size):
            self._collection.push_back(n._this.clone())
            self._collection.back().get().calc_attractors()

    def fill_with_mutations(self, Network n, np.int_t[:] mutations):
        cdef:
            size_t i
            size_t sz = mutations.shape[0]
            cFactory *con = n._this.factory.get()

        for i in range(sz):
            self._collection.push_back(
                con.clone_and_mutate_network(n._shared, mutations[i], 1))


    property fitnesses:
        def __get__(self):
            fits = numpy.zeros(self.size)
            self.get_fitnesses(fits)
            return fits

    property generations:
        def __get__(self):
            gens = numpy.zeros(self.size, dtype=int)
            cdef:
                size_t i
                np.int_t[:] c_gens = gens
            for i in range(self._collection.size()):
                c_gens[i] = deref(self._collection)[i].get().generation
            return gens

    property active_bindings:
        def __get__(self):
            bindings = numpy.zeros(self.size, dtype=int)
            cdef:
                size_t i
                np.int_t[:] c_bindings = bindings
                cNetworkAnalysis *analysis

            for i in range(self._collection.size()):
                # TODO: make this more sensible 
                analysis = new cNetworkAnalysis(deref(self._collection)[i])
                c_bindings[i] = analysis.calc_active_bindings()
                del analysis

            return bindings


cdef class Collection(CollectionBase):
    def __cinit__(self, Factory c, size_t size=0):
        self._this = new cNetworkVector()
        self._collection = self._this

    def __dealloc__(self):
        del self._this


cdef class Ancestry(CollectionBase):
    def __cinit__(self, Factory c, size_t size=0):
        self._this = new cNetworkVector()
        self._collection = self._this

    def __repr__(self):
        if self.size == 0:
            return "<Ancestry: EMPTY>"

        return "<Ancestry: {}N, G{}-G{}>".format(
            self.size,
            self._this.front().get().generation,
            self._this.back().get().generation,
        )

    def network_at_generation(self, size_t gen):
        # Should really use binary search ...
        cdef size_t i
        for i in range(self._this.size() - 1):
            if deref(self._this)[i].get().generation > gen:
                break
        else:
            # If we hit the end actually return the last one
            i += 1
        return self.get_at(i - 1)
        

    def __dealloc__(self):
        del self._this


cdef class Population(CollectionBase):
    def __cinit__(self, Factory c, size_t size=0):
        self._this = new cPopulation(c._shared, size)
        self._collection = &(self._this.networks)

    def __dealloc__(self):
        del self._this

    def worst_and_best(self):
        return self._this.worst_and_best()

    def best_indexes(self):
        cdef cIndexes best;
        self._this.best_indexes(best)
        return best

    def get_best(self, maxn=sys.maxint):
        best = []
        for i, ndx in enumerate(self.best_indexes()):
            if i == maxn:
                break
            best.append(self[ndx])
        return best

    property site_count:
        def __get__(self):
            return self._this.factory.get().site_count(self._this.networks)

    def _get_identifiers(self, np.int_t[:] output):
        """A list of identifiers for the current population"""
        # TODO: maybe this should be pre-allocated?
        assert output.shape[0] == self._this.networks.size()
        cdef size_t i
        for i in range(self._this.networks.size()):
            output[i] = self._this.networks[i].get().identifier

    def __repr__(self):
        return "<Population: {}>".format(self.size)

    def mutate(self, double site_rate, int generation=0):
        return self._this.mutate(site_rate, generation)

    def select(self, SelectionModel sm, size=None):
        cdef size_t s
        if size is None:
            s = self._this.networks.size()
        else:
            s = size
        return self._this.select(deref(sm._this), s)

    property mutated:
        def __get__(self):
            return self._this.mutated

    property selected:
        def __get__(self):
            return self._this.selected


cdef class BaseTarget:
    def __cinit__(self):
        self._base = NULL

    def __dealloc__(self):
        if self._base != NULL:
            del self._base

    def assess(self, Network net):
        # assert net.factory.world is self.world
        return self._base.assess(deref(net._this));

    def assess_collection(self, CollectionBase coll):
        self._base.assess_networks(deref(coll._collection));

    property weighting:
        def __get__(self):
            return self._base.weighting

        def __set__(self, vector[double] wghts):
            self._base.set_weighting(wghts)

    property identifier:
        def __get__(self):
            return self._base.identifier

    property name:
        def __get__(self):
            return self._base.name


# NOTE: allow weighting to be None for backward compatibility
def _default_target(World w, name, ident, rates, weighting=None,
                      scoring_method=None, strength=None):
    t = DefaultTarget(w, None, name, ident=ident)
    # Manually construct these
    t._this.optimal_rates = rates
    if weighting is not None:
        t.weighting = weighting
    if scoring_method is not None:
        t.scoring_method = scoring_method
    if strength is not None:
        t.strength = strength

    return t

cdef class DefaultTarget(BaseTarget):
    def __cinit__(self, World w, init_func=None, name="", ident=-1,
                  scoring_method=0, strength=0.0):
        self.world = w
        self._this = new cDefaultTarget(w._shared, name, ident, 
                                        scoring_method, strength)
        self._base = self._this
        if init_func:
            self._construct_from_function(init_func)

    def _construct_from_function(self, init_func):
        a, b = self.world._this.cue_range

        # Slow and cumbersome, but it doesn't matter
        for i, e in enumerate(self.world.environments):
            # TODO: Clean up the refs here
            outputs = init_func(*e.as_array()[a:b])
            try:
                s = len(outputs)
            except TypeError:
                # Must be a single value...
                outputs = [outputs]
                s = 1

            if len(outputs) != self.world.out_channels:
                raise RuntimeError(
                    "return value of Target function must be length %s" \
                    % self.world.out_channels)

            for j, val in enumerate(outputs):
                self._this.optimal_rates[i][j] = float(val)


    def __reduce__(self):
        return _default_target, (self.world, 
                                   self._this.name, 
                                   self._this.identifier,
                                   self._this.optimal_rates,
                                   self._this.weighting,
                                   self._this.scoring_method,
                                   self._this.strength)

    def as_array(self):
        return numpy.array(self._this.optimal_rates)

    property scoring_method:
        def __get__(self):
            return self._this.scoring_method
        def __set__(self, ScoringMethod method):
            self._this.scoring_method = method

    property strength:
        def __get__(self):
            return self._this.strength
        def __set__(self, double s):
            self._this.strength = s

    def calc_categories(self):
        """Categorise the targets"""
        cat_dict = {}
        cats = []
        cat_n = 0
        for et in self.as_array():
            et = tuple(et)
            if et in cat_dict:
                cats.append(cat_dict[et])
            else:
                cat_dict[et] = cat_n
                cats.append(cat_n)
                cat_n += 1
        return cats
            
    def calc_distinct_outputs(self):
        out = set()
        for et in self.as_array():
            out.add(tuple(et))
        return out
            
