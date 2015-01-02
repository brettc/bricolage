# distutils: include_dirs = NUMPY_PATH
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

import cython
import numpy
# from operand import Operand

# cimports
cimport numpy as np
from cython.operator import dereference as deref, preincrement as preinc

import random

reserved_channels = 2
off_channel = 0
on_channel = 1

cdef class ChannelStateFrozen:
    def __cinit__(self):
        pass

    cdef init(self, cWorld_ptr &w, cChannelState &p):
        self.cworld_ptr = w
        self.cchannel_state = p

    def test(self, size_t i):
        return self.cchannel_state.test(i)

    def __getitem__(self, size_t i):
        return self.cchannel_state.test(i)

    def __cmp__(self, ChannelStateFrozen other):
        return bitset_cmp(self.cchannel_state, other.cchannel_state)

    property size:
        def __get__(self):
            return self.cchannel_state.size()

    def __str__(self):
        # TODO: Clean this function up---it is really crappy
        cdef:
            string cstr
            cWorld *w = self.cworld_ptr.get()

        if w == NULL:
            # TODO: Some proper exceptions would be good...
            raise RuntimeError

        cdef size_t cuereg = w.cue_channels + w.reg_channels

        to_string(self.cchannel_state, cstr)

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
        other.init(self.cworld_ptr, self.cchannel_state)
        return other

    def copy(self):
        return self.__copy__()

    def as_array(self):
        vals = numpy.zeros(self.size, dtype=numpy.int32)
        cdef: 
            np.npy_int32[:] v = vals
            size_t i

        for i in range(self.cchannel_state.size()):
            v[i] = self.cchannel_state.test(i)

        return vals

    def __repr__(self):
        return "<ChannelsRO: {}>".format(self.__str__())


cdef class ChannelState(ChannelStateFrozen):

    def __cinit__(self):
        pass

    def set(self, size_t i):
        self.cchannel_state.set(i)

    def reset(self, size_t i):
        self.cchannel_state.reset(i)

    def flip(self, size_t i):
        self.cchannel_state.flip(i)

    def __setitem__(self, size_t i, bint b):
        if b:
            self.cchannel_state.set(i)
        else:
            self.cchannel_state.reset(i)

    def merge(self, ChannelStateFrozen other):
        self.cchannel_state |= other.cchannel_state

    def __repr__(self):
        return "<Channels: {}>".format(self.__str__())


cdef class World:
    def __cinit__(self, params):
        self.cworld_ptr = cWorld_ptr(new cWorld(
            params.seed, 
            params.cue_channels, 
            params.reg_channels, 
            params.out_channels,
        ))
        self.cworld = self.cworld_ptr.get()

        self.reserved_signals = set([on_channel, off_channel])
        self.cue_signals = set(range(*self.cworld.cue_range))
        self.reg_signals = set(range(*self.cworld.reg_range))
        self.out_signals = set(range(*self.cworld.out_range))
        self.sub_signals = set(range(*self.cworld.sub_range))
        self.pub_signals = set(range(*self.cworld.pub_range))


    def create_state(self):
        c = ChannelState()
        c.cchannel_state.resize(self.cworld.channel_count)
        c.cworld_ptr = self.cworld_ptr
        return c

    property environments:
        def __get__(self):
            if self._environments is not None:
                return self._environments

            cdef vector[cChannelState].iterator i =  self.cworld.environments.begin()
            envs = []
            while i != self.cworld.environments.end():
                p = ChannelStateFrozen()
                p.init(self.cworld_ptr, deref(i))
                envs.append(p)
                preinc(i)

            self._environments = envs
            return envs

    # Some routines just for testing, mainly for testing (NOT fast)
    def seed_random_engine(self, int s):
        self.cworld.rand.seed(s)

    def get_random_double(self, double low, double high):
        return self.cworld.get_random_double(low, high)

    def get_random_int(self, int low, int high):
        return self.cworld.get_random_int(low, high)

    property cue_channels:
        def __get__(self):
            return self.cworld.cue_channels
    property reg_channels:
        def __get__(self):
            return self.cworld.reg_channels
    property out_channels:
        def __get__(self):
            return self.cworld.out_channels
    property channel_count:
        def __get__(self):
            return self.cworld.channel_count
    # property off_channel:
    #     def __get__(self):
    #         return off_channel
    # property on_channel:
    #     def __get__(self):
    #         return off_channel

    property cue_range:
        def __get__(self):
            return self.cworld.cue_range
    property reg_range:
        def __get__(self):
            return self.cworld.reg_range
    property out_range:
        def __get__(self):
            return self.cworld.out_range
    property sub_range:
        def __get__(self):
            return self.cworld.sub_range
    property pub_range:
        def __get__(self):
            return self.cworld.pub_range
    
    def name_for_channel(self, c):
        sz = 1
        if c == off_channel:
            return "OFF"
        if c == on_channel:
            return "ON"
        # We use "mathematical" number, starting at 1.0
        if c in self.out_signals:
            return "P{0:0{1:}d}".format(
                c + 1 - self.cworld.out_range.first, sz)
        if c in self.reg_signals:
            return "T{0:0{1:}d}".format(
                c + 1 - self.cworld.reg_range.first, sz)
        return "E{0:0{1:}d}".format(
            c + 1 - self.cworld.cue_range.first, sz)

cdef class Constructor:
    def __cinit__(self, World w, params):
        self.world = w
        self.gene_class = Gene
        self.module_class = CisModule

cdef class Network:
    def __cinit__(self, Constructor c):
        self.constructor = c
        cdef cNetwork_ptr ptr = c._this.construct()
        self.bind_to(ptr)

    cdef bind_to(self, cNetwork_ptr &ptr):
        self.ptr = ptr
        self.cnetwork = self.ptr.get()
        self.cnetwork.pyobject = <void *>(self)

    def __dealloc__(self):
        self.cnetwork.pyobject = <void *>0

    property identifier:
        def __get__(self):
            return self.cnetwork.identifier

    property parent_identifier:
        def __get__(self):
            return self.cnetwork.parent_identifier

    property gene_count:
        def __get__(self):
            return self.cnetwork.gene_count()

    property genes:
        def __get__(self):
            if self._genes is None:
                self._genes = tuple(self.constructor.gene_class(self, i) for i in range(self.gene_count))
            return self._genes

    def cycle(self, ChannelState c):
        self.cnetwork.cycle(c.cchannel_state)

    def cycle_with_intervention(self, ChannelState c):
        self.cnetwork.cycle_with_intervention(c.cchannel_state)

    def _evil_mutate(self, size_t nmutations):
        """Mutate the network. 

        This invalidates lots of assumptions required for selection to work
        correctly. Use only if you understand what you are doing.
        """
        self.cnetwork.mutate(nmutations)
        self._invalidate_cached()

    def _invalidate_cached(self):
        self.cnetwork.calc_attractors()
        self._attractors = None
        self._rates = None

    property attractors:
        """A tuple containing the attractors for each environment"""
        def __get__(self):
            # Have if we're detached, we have to recalc every time
            if self._attractors is not None:
                return self._attractors

            cdef:
                vector[cChannelStateVector].iterator cattr_iter
                vector[cChannelState].iterator cstate_iter

            attrs = []
            cattr_iter = self.cnetwork.attractors.begin()
            while cattr_iter != self.cnetwork.attractors.end():
                attr = []
                cstate_iter = deref(cattr_iter).begin()
                while cstate_iter != deref(cattr_iter).end():
                    c = ChannelStateFrozen()
                    c.init(self.constructor._this.world, deref(cstate_iter))
                    attr.append(c)
                    preinc(cstate_iter)

                # Make everything a tuple -- you shouldn't be able to mess
                # with it!
                attrs.append(tuple(attr))
                preinc(cattr_iter)

            # Again: tuples, so you can't mess with it.
            self._attractors = tuple(attrs)
            return self._attractors

    property rates:
        """Return a readonly numpy array of the rates"""
        def __get__(self):
            # Lazy evaluation
            if self._rates is not None:
                return self._rates

            cdef cWorld *world = self.constructor._this.world.get()

            # Construct the numpy array via python
            r = numpy.zeros((world.environments.size(), world.out_channels))

            cdef:
                np.npy_double[:,:] c_r = r
                size_t i, j

            # Copy in the goods using a memory array
            for i in range(world.environments.size()):
                for j in range(world.out_channels):
                    c_r[i, j] = self.cnetwork.rates[i][j]

            # Don't mess with it!
            r.flags.writeable = False
            self._rates = r
            return self._rates

    property fitness:
        def __get__(self):
            return self.cnetwork.fitness

    property target:
        def __get__(self):
            return self.cnetwork.target

    def __repr__(self):
        return "<Network id:{} pt:{}>".format(self.identifier, self.parent_identifier)


cdef class Gene:
    """A proxy to a gene.
    """
    def __cinit__(self, Network n, size_t g):
        """This is not meant to be called publicly"""
        self.network = n
        self.gene_number = g
        self.cgene = self.network.cnetwork.get_gene(g)

    property sequence:
        def __get__(self):
            return self.cgene.sequence

    property pub:
        def __get__(self):
            return self.cgene.pub

    property intervene:
        def __get__(self):
            return self.cgene.intervene

    # def _evil_set_pub(self, size_t p):
    #     assert p in self.network.world.pub_signals
    #     # self.cgene.pub = p
    #     self.network._invalidate_cached()

    property module_count:
        def __get__(self):
            return self.cgene.module_count()

    property modules:
        def __get__(self):
            # Lazy construction
            if self._modules is None:
                self._modules = tuple(self.network.constructor.module_class(self, i) 
                    for i in range(self.module_count))
            return self._modules

    def __repr__(self):
        w = self.network.constructor.world
        return "<Gene[{}]: {}>".format(self.sequence, w.name_for_channel(self.pub))


cdef class CisModule:
    """A proxy to a CisModule.
    """
    def __cinit__(self, Gene g, size_t i):
        self.gene = g
        assert i < g.cgene.module_count()
        self.ccismodule = g.cgene.get_module(i)

    def site_count(self):
        return self.ccismodule.site_count()

    def get_site(self, size_t i):
        assert i < self.ccismodule.site_count()
        return self.ccismodule.get_site(i)

    def set_site(self, size_t i, size_t c):
        assert c in self.gene.network.world.sub_signals
        old = self.ccismodule.set_site(i, c)
        self.gene.network._invalidate_cached()
        return old

    property intervene:
        def __get__(self):
            return self.ccismodule.intervene
        # def __set__(self):
        #     return self.ccismodule.intervene

    property channels:
        def __get__(self):
            return tuple(self.ccismodule.channels[i]
                         for i in range(self.ccismodule.site_count()))
        def __set__(self, t):
            assert len(t) == self.ccismodule.site_count()
            valid = self.gene.network.constructor.world.sub_signals
            cdef size_t i, c
            for i in range(self.ccismodule.site_count()):
                assert t[i] in valid
                self.ccismodule.channels[i] = t[i]
            self.gene.network._invalidate_cached()

    property channel_names:
        def __get__(self):
            w = self.gene.network.constructor.world
            return [w.name_for_channel(c) for c in self.channels]


cdef class NetworkCollection:
    def __cinit__(self, Constructor c, size_t size):
        self.constructor = c
        self.cnetworks.reserve(size)
        cdef size_t i
        for i in range(size):
            self.cnetworks.push_back(c._this.construct())

    def add(self, Network n):
        # TODO: Mix different network types?
        assert n.constructor is self.constructor
        self.cnetworks.push_back(n.ptr)

    cdef object get_at(self, size_t i):
        if i >= self.cnetworks.size():
            raise IndexError

        cdef:
            cNetwork_ptr ptr = self.cnetworks[i]
            cNetwork *net = ptr.get()

        # Is there an existing python object?
        # Note: cython automatically increments the reference count when we do
        # this (which is what we want). This move just means we get object
        # identity, at least while one python reference continues to exist.
        if net.pyobject:
            return <object>(net.pyobject)

        # Nope. We need to create a new python wrapper object
        n = Network(self.constructor)
        n.bind_to(ptr)
        return n

    property size:
        def __get__(self):
            return self.cnetworks.size()

    property site_count:
        def __get__(self):
            return self.constructor._this.site_count(self.cnetworks)

    def __getitem__(self, size_t i):
        return self.get_at(i)

    def __iter__(self):
        for i in range(self.size):
            yield self.get_at(i)

    def __repr__(self):
        return "<NetworkCollection: {}>".format(self.size)

    def mutate(self, double site_rate):
        cdef cIndexes mutated
        self.constructor._this.mutate_collection(
            self.cnetworks, mutated, site_rate)
        # Return indexes of the mutated networks. Automatic conversion (thank
        # you Cython)
        return mutated

    def select(self, Target target):
        """Inplace selection of networks, replacing networks"""
        cdef:
            cSelectionModel *sm
            bint sel
            cIndexes indexes
            cNetworkVector new_networks

        sm = new cSelectionModel(self.constructor._this.world)
        sel = sm.select(self.cnetworks, deref(target.ctarget), 
                     self.cnetworks.size(), indexes)
        if sel:
            sm.copy_using_indexes(self.cnetworks, new_networks, indexes)

            # Replace everything -- this is fast
            self.cnetworks.swap(new_networks)

        # Get rid of this
        del sm

        return sel

    def selection_indexes(self, Target target):
        """Just return the indices where selection would happen.

        ***This is used for testing***
        """
        cdef:
            cSelectionModel *sm
            bint sel
            cIndexes indexes
            cNetworkVector new_networks

        sm = new cSelectionModel(self.constructor._this.world)

        sel = sm.select(self.cnetworks, deref(target.ctarget), 
                     self.cnetworks.size(), indexes)

        return indexes

cdef class Target:
    def __cinit__(self, World w, init_func):
        self.world = w
        self.ctarget = new cTarget(w.cworld_ptr)
        a, b = w.cworld.cue_range

        # Slow and cumbersome, but it doesn't matter
        for i, e in enumerate(w.environments):
            # TODO: Clean up the refs here
            outputs = init_func(*e.as_array()[a:b])
            try:
                s = len(outputs)
            except TypeError:
                # Must be a single value...
                outputs = [outputs]
                s = 1

            if len(outputs) != w.out_channels:
                raise RuntimeError(
                    "return value of Target function must be length %s" % w.out_channels)

            for j, val in enumerate(outputs):
                self.ctarget.optimal_rates[i][j] = float(val)

    def __dealloc__(self):
        del self.ctarget

    def assess(self, Network net):
        assert net.constructor.world is self.world
        return self.ctarget.assess(deref(net.cnetwork));

    def as_array(self):
        return numpy.array(self.ctarget.optimal_rates)


cdef class NetworkAnalysis:
    def __cinit__(self, Network net):
        self.network = net
        self.canalysis = new cNetworkAnalysis(net.ptr)

    def __dealloc__(self):
        del self.canalysis

    def get_edges(self):
        cdef:
            cEdgeList edges
        self.canalysis.make_edges(edges)
        return edges

    def get_active_edges(self):
        cdef:
            cEdgeList edges
        self.canalysis.make_active_edges(edges)
        return edges
