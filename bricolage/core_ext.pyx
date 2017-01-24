# distutils: include_dirs = NUMPY_PATH
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# distutils: define_macros=NPY_NO_DEPRECATED_API

from analysis cimport *
from targets_ext cimport *

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

cdef class Channels:
    def __cinit__(self, World w, bits_t init=0):
        self.world = w
        self._this.bits = init
        self.size = w._this.channel_count

    # cdef _assign(self, bits_t bits):
    #     self._this.bits = bits

    def __hash__(self):
        return self._this.bits

    def max(self):
        return self._this.max()

    def set(self, index_t i):
        self._this.set(i, self.size)

    def clear(self, index_t i):
        self._this.clear(i, self.size)

    def flip(self, index_t i):
        self._this.flip(i, self.size)

    def test(self, index_t i):
        return self._this.test(i, self.size)

    def __str__(self):
        return self._this.to_string(self.size)

    def __repr__(self):
        bit_string = self._this.to_string(self.size)
        return "<Channels: {}>".format(bit_string)

    def __cmp__(self, Channels other):
        return bitset_cmp(self._this, other._this)

    # def __str__(self):
    #     # TODO: Clean this function up---it is really crappy
    #     cdef:
    #         string cstr
    #         cWorld *w = self.world.get()
    #
    #     if w == NULL:
    #         # TODO: Some proper exceptions would be good...
    #         raise RuntimeError
    #
    #     cdef size_t cuereg = w.cue_channels + w.reg_channels
    #
    #     to_string(self._this, cstr)
    #
    #     # I think it is much easier to understand if we reverse it
    #     # Also, clip the reserved channels 
    #     # TODO: fix this nasty hack
    #     s = cstr[::-1]
    #     res = s[:2]
    #     env = s[2:2+w.cue_channels]
    #     reg = s[2+w.cue_channels:2+cuereg]
    #     out = s[2+cuereg:]
    #
    #     return "|".join((res, env, reg, out))

    def __copy__(self):
        other = Channels(self.world, self._this.bits)
        return other

    def copy(self):
        return self.__copy__()

    def as_array(self):
        vals = numpy.zeros(self.size, dtype=numpy.int32)
        cdef: 
            np.npy_int32[:] v = vals
            size_t i

        for i in range(self.size):
            v[i] = self._this.test(i, self.size)

        return vals

    def merge(self, Channels other):
        assert self.size == other.size
        self._this.unchecked_union(other._this)
        return self

    def filter(self, Channels other):
        assert self.size == other.size
        self._this.unchecked_intersection(other._this)
        return self


def _construct_world(params, net_id, target_id, r_state):
    w = World(params)
    w._this.next_network_identifier = net_id
    w._this.next_target_identifier = target_id
    w._this.set_random_state(r_state)
    return w


cdef class World:
    def __cinit__(self, params):
        self._params = copy.deepcopy(params)
        # Hack for old stuff
        if not hasattr(self._params, 'reg_gene_count'):
            self._params.reg_gene_count = self._params.reg_channels

        self._shared = cWorld_ptr(new cWorld(
            self._params.seed, 
            self._params.cue_channels, 
            self._params.reg_channels, 
            self._params.out_channels,
            self._params.reg_gene_count,
        ))
        self._this = self._shared.get()

        if hasattr(self._params, 'input_type'):
            self._this.input_type = params.input_type

        if hasattr(self._params, 'pulse_for'):
            self._this.pulse_for = params.pulse_for

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

    property environments:
        def __get__(self):
            intvec = deref(to_cAttractorBits(&self._this.environments))
            return [Channels(self, x) for x in intvec]

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

    property draw_from_subs:
        def __get__(self):
            return self._this.draw_from_subs
        def __set__(self, cIndexes ch):
            self._this.set_draw_from_subs(ch)

    property draw_from_regs:
        def __get__(self):
            return self.this.draw_from_regs
        def __set__(self, cIndexes ch):
            self._this.set_draw_from_regs(ch)

    property gene_count:
        def __get__(self):
            return self._this.gene_count
    property module_count:
        def __get__(self):
            return self._this.module_count

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
                vector[cAttractor].iterator cattr_iter
                np.npy_double[:] c_r = r
                size_t i = 0

            cattr_iter = self._this.attractors.begin()
            while cattr_iter != self._this.attractors.end():
                c_r[i] = deref(cattr_iter).size()
                i += 1
                preinc(cattr_iter)

            return r

    def cycle(self, Channels c):
        self._this.cycle(c._this)

    def cycle_with_intervention(self, Channels c):
        self._this.cycle_with_intervention(c._this)

    # def calc_perturbation(self, bint env_only):
    #     self._this.calc_perturbation(env_only)

    def recalculate(self, with_intervention=False):
        if with_intervention:
            self._this.calc_attractors_with_intervention()
        else:
            self._this.calc_attractors()

        self._attractors = None
        self._rates = None

    def mutate(self, size_t n_cis, size_t n_trans=0):
        """Mutate the network. 

        This invalidates lots of assumptions required for selection to work
        correctly. Use only if you understand what you are doing.
        """
        self._this.mutate(n_cis, n_trans)
        self.recalculate()

    def duplicate(self, size_t nmutations):
        """Mutate the network. 

        This invalidates lots of assumptions required for selection to work
        correctly. Use only if you understand what you are doing.
        """
        self._this.duplicate(nmutations)
        self.recalculate()

    cdef _make_python_attractor(self, cAttractor &c_attr):
        w = self.factory.world
        intvec = deref(to_cAttractorBits(&c_attr))
        return [Channels(w, x) for x in intvec]

    cdef _make_python_attractors(self, cAttractorSet &attrs):
        w = self.factory.world
        intvecvec = deref(to_cAttractorSetBits(&attrs))
        return [[Channels(w, x) for x in intvec] for intvec in intvecvec]

    cdef _make_python_rate(self, cRates &rates):
        cdef cWorld *world = self.factory._this.world.get()

        # Construct the numpy array via python
        r = numpy.zeros((world.out_channels))

        cdef:
            np.npy_double[:] c_r = r
            size_t i

        # Copy in the goods using a memory array
        for i in range(world.out_channels):
            c_r[i] = rates[i]

        # Don't mess with it!
        r.flags.writeable = False
        return r

    cdef _make_python_rates(self, cRatesVector &rates_vec):
        cdef cWorld *world = self.factory._this.world.get()

        # Construct the numpy array via python
        r = numpy.zeros((world.environments.size(), world.out_channels))

        cdef:
            np.npy_double[:,:] c_r = r
            size_t i, j

        # Copy in the goods using a memory array
        for i in range(world.environments.size()):
            for j in range(world.out_channels):
                c_r[i, j] = rates_vec[i][j]

        # Don't mess with it!
        r.flags.writeable = False
        return r

    def stabilise(self, external, intervention=False):
        cdef:
            cAttractor c_attr, c_trans, c_external
            cRates c_rates

        if isinstance(external, Channels):
            external = [external]

        for channels in external:
            c_channels = <Channels?>(channels)
            c_external.push_back(c_channels._this)

        self._this.stabilise(c_external, intervention, c_attr, c_trans, c_rates)
        attr = self._make_python_attractor(c_attr)
        trans = self._make_python_attractor(c_trans)
        rates = self._make_python_rate(c_rates)
        return attr, trans, rates

    def get_rates(self, Channels c, bint use_cache=True):
        cdef:
            cRates c_rates
        self._this.get_rates(c._this, c_rates, use_cache)
        return c_rates

    property attractors:
        """A tuple containing the attractors for each environment"""
        def __get__(self):
            # Have if we're detached, we have to recalc every time
            if self._attractors is not None:
                return self._attractors

            self._attractors = self._make_python_attractors(self._this.attractors)
            return self._attractors

    property rates:
        """Return a readonly numpy array of the rates"""
        def __get__(self):
            # Lazy evaluation
            if self._rates is not None:
                return self._rates
            self._rates = self._make_python_rates(self._this.rates)
            return self._rates

    property rates_cache:
        def __get__(self):
            themap = deref(to_cChannelRatesMapBits(&self._this.cached_mappings))
            return dict([(Channels(self.factory.world, k), v) 
                         for k, v in themap.items()])

    def clear_rate_cache(self):
        self._this.clear_rate_cache()

    property fitness:
        def __get__(self):
            return self._this.fitness

    property target:
        def __get__(self):
            return self._this.target

    property attractor_robustness:
        def __get__(self):
            return self._this.attractor_robustness()

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
        assert c in self.gene.network.factory.world.sub_signals
        old = self._this.set_site(i, c)
        self.gene.network.recalculate(with_intervention=True)
        return old

    property intervene:
        def __get__(self):
            return self._this.intervene
        def __set__(self, InterventionState i):
            self._this.intervene = i
            self.gene.network.recalculate(with_intervention=True)

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
        cdef Network n = self.factory.network_class(
                self.factory, 
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
                con.clone_and_mutate_network(n._shared, mutations[i], 0, 0, 1))

    property generations:
        def __get__(self):
            gens = numpy.zeros(self.size, dtype=int)
            cdef:
                size_t i
                np.int_t[:] c_gens = gens
            for i in range(self._collection.size()):
                c_gens[i] = deref(self._collection)[i].get().generation
            return gens

    # TODO: turn this into a function (it is deceptive)
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

    def robustness(self):
        robustness = numpy.zeros(self.size, dtype=float)
        cdef:
            size_t i
            np.double_t[:] c_robustness = robustness

        for i in range(self._collection.size()):
            # TODO: make this more sensible 
            c_robustness[i] = deref(self._collection)[i].get().attractor_robustness()

        return robustness

    def active_cis(self):

        cdef:
            cFactory *fact = self.factory._this

        mod_count = fact.gene_count * fact.module_count
        bindings = numpy.zeros((mod_count), dtype=float)

        cdef:
            size_t i, offset #, gene_id, cis_id, module_id
            np.npy_double[:] c_bindings = bindings
            Edge_t edge 
            cEdgeList edges
            std_set[Edge_t].iterator edges_iter
            cNetworkAnalysis *analysis

        for i in range(self._collection.size()):
            # TODO: make this more sensible 
            analysis = new cNetworkAnalysis(deref(self._collection)[i])
            analysis.make_active_edges(edges)
            edges_iter = edges.begin()
            while edges_iter != edges.end():
                edge = deref(edges_iter)
                if edge.second.first == NT_MODULE:
                    module_id = edge.second.second
                    gene_id = (0xff00 & module_id) >> 8
                    cis_id = 0xff & module_id
                    offset = gene_id * fact.module_count + cis_id
                    c_bindings[offset] += 1.0
                preinc(edges_iter)


            del analysis

        return bindings / self._collection.size()


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

    def assess(self, BaseTarget target):
        self._this.assess(deref(target._base))

    def select(self, SelectionModel sm, size=None):
        cdef size_t s
        if size is None:
            s = self._this.networks.size()
        else:
            s = size
        return self._this.select(deref(sm._this), s)

    def mutate(self, double cis_rate, double trans_rate=0.0, double dup_rate=0.0, int generation=0):
        return self._this.mutate(cis_rate, trans_rate, dup_rate, generation)

    property mutated:
        def __get__(self):
            return self._this.mutated

    property selected:
        def __get__(self):
            return self._this.selected

    property fitnesses:
        def __get__(self):
            ret = numpy.empty(self._this.fitnesses.size())
            cdef:
                np.npy_double[:] np_ret = ret
                size_t i

            for i in range(self._this.fitnesses.size()):
                np_ret[i] = self._this.fitnesses[i]
            return ret


cdef class SelectionModel:
    def __cinit__(self, World w, bint relative=False):
        self.world = w
        self._this = new cSelectionModel(w._shared, relative)

    def __dealloc__(self):
        del self._this



