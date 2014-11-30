# distutils: include_dirs = NUMPY_PATH
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

import cython
import numpy
from operand import Operand

# cimports
cimport numpy as np
from cython.operator import dereference as deref, preincrement as preinc
from _cpp cimport *
from pubsub2_ext cimport *

import random
cdef int magic_number = random.randint(0, 100000)

cdef class ChannelStateFrozen:
    cdef:
        cChannelState cchannel_state
        cFactory_ptr cfactory_ptr

    def __cinit__(self):
        pass

    cdef init(self, cFactory_ptr &f, cChannelState &p):
        self.cfactory_ptr = f
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
        # TODO: Clean this function up---it is a bit rough.
        cdef:
            string cstr
            cFactory *f = self.cfactory_ptr.get()

        if f == NULL:
            # TODO: Some proper exceptions would be good...
            raise RuntimeError

        cdef size_t cuereg = f.cue_channels + f.reg_channels

        to_string(self.cchannel_state, cstr)

        # I think it is much easier to understand if we reverse it
        # Also, clip the "silencing" channel 0
        # TODO: fix this nasty hack
        s = cstr[::-1][1:]
        env = s[:f.cue_channels]
        reg = s[f.cue_channels:cuereg]
        out = s[cuereg:]

        return env + '|' + reg + '|' + out

    def __copy__(self):
        other = ChannelState()
        other.init(self.cfactory_ptr, self.cchannel_state)
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
        return "<ChannelStateFrozen: {}>".format(self.__str__())


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
        return "<ChannelState: {}>".format(self.__str__())


cdef class Factory:
    cdef:
        cFactory_ptr cfactory_ptr
        cFactory *cfactory
        cGeneFactory *cgenefactory
        readonly:
            object params

        object _environments

    def __cinit__(self, params):
        self.params = params
        self.cfactory_ptr = cFactory_ptr(new cFactory(params.seed))
        self.cfactory = self.cfactory_ptr.get()

        # Now translate everything to cpp -- this needs to be done first.
        # There is automatic conversion for most of this (stl::vector etc)
        self.cfactory.gene_count = params.gene_count
        self.cfactory.cis_count = params.cis_count
        self.cfactory.operands = params.operands
        self.cfactory.sub_range = params.sub_range
        self.cfactory.pub_range = params.pub_range
        self.cfactory.out_range = params.out_range
        self.cfactory.cue_channels = params.cue_channels
        self.cfactory.reg_channels = params.reg_channels
        self.cfactory.out_channels = params.out_channels

        self.cfactory.init_environments()

        # Create a mutator
        self.cgenefactory = new cGeneFactoryLogic2(self.cfactory, 
                                         params.gene_mutation_rate)

    def __dealloc__(self):
        del self.cgenefactory

    def create_state(self):
        c = ChannelState()
        c.cchannel_state.resize(self.cfactory.total_channels)
        c.cfactory_ptr = self.cfactory_ptr
        return c

    def create_network(self):
        cdef cNetwork_ptr ptr = cNetwork_ptr(new cNetwork(self.cfactory_ptr))
        self.cgenefactory.construct_network(deref(ptr.get()))
        n = Network(self)
        n.bind_to(ptr)
        return n

    def create_collection(self, size_t size):
        cdef:
            cNetwork_ptr ptr
            size_t i

        nc = NetworkCollection(self)
        for i in range(size):
            ptr = cNetwork_ptr(new cNetwork(self.cfactory_ptr))
            self.cgenefactory.construct_network(deref(ptr.get()))
            nc.cnetworks.push_back(ptr)

        return nc

    property environments:
        def __get__(self):
            if self._environments is not None:
                return self._environments

            cdef vector[cChannelState].iterator i =  self.cfactory.environments.begin()
            envs = []
            while i != self.cfactory.environments.end():
                p = ChannelStateFrozen()
                p.init(self.cfactory_ptr, deref(i))
                envs.append(p)
                preinc(i)

            self._environments = envs
            return envs

    # Some routines just for testing, mainly for testing (NOT fast)
    def seed_random_engine(self, int s):
        self.cfactory.random_engine.seed(s)

    def get_random_double(self, double low, double high):
        return self.cfactory.get_random_double(low, high)

    def get_random_int(self, int low, int high):
        return self.cfactory.get_random_int(low, high)


cdef class Target:
    cdef:
        cTarget *ctarget
        readonly:
            Factory factory

    def __cinit__(self, Factory f, init_func):
        self.factory = f
        self.ctarget = new cTarget(f.cfactory)

        # Slow and cumbersome, but it doesn't matter
        for i, e in enumerate(f.environments):
            # TODO: Clean up the refs here
            outputs = init_func(e.as_array()[1:f.params.cue_channels+1])
            if len(outputs) != f.params.out_channels:
                raise RuntimeError

            for j, val in enumerate(outputs):
                self.ctarget.optimal_rates[i][j] = float(val)

    def __dealloc__(self):
        del self.ctarget

    def assess(self, Network net):
        assert net.factory is self.factory
        return self.ctarget.assess(deref(net.cnetwork));

    def as_array(self):
        return numpy.array(self.ctarget.optimal_rates)


cdef class Network:
    cdef:
        readonly:
            Factory factory
            bint ready
            bint dirty

        # Because we hold a reference to the shared_ptr, we know we can always
        # safely access the ACTUAL pointer. We keep the pointer around too, as
        # it makes our life easier. The cost is a tiny bit of space.
        cNetwork_ptr ptr
        cNetwork *cnetwork
        object _genes, _attractors, _rates

    # Networks take two stages. You need to create the python object and then
    # bind it to a cNetwork_ptr.
    def __cinit__(self, Factory factory):
        self.factory = factory
        self.ready = False
        self.dirty = False
        self._genes = None

    cdef bind_to(self, cNetwork_ptr &ptr):
        self.ptr = ptr
        self.cnetwork = self.ptr.get()
        self.cnetwork.pyobject = <void *>(self)
        self.ready = True

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
            return self.cnetwork.genes.size()

    property genes:
        def __get__(self):
            if self._genes is None:
                self._genes = tuple(Gene(self, i) for i in range(self.gene_count))
            return self._genes

    def cycle(self, ChannelState c):
        self.cnetwork.cycle(c.cchannel_state)

    def mutated(self, size_t nmutations):
        """Return a mutated version..."""
        cdef:
            cGeneFactory *gf = self.factory.cgenefactory
            cNetwork_ptr new_ptr

        new_ptr = gf.copy_and_mutate_network(self.ptr, nmutations)
        new_net = Network(self.factory)
        new_net.bind_to(new_ptr)
        return new_net

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
                    c.init(self.factory.cfactory_ptr, deref(cstate_iter))
                    attr.append(c)
                    preinc(cstate_iter)

                # Make everything a tuple -- you shouldn't be able to mess
                # with it!
                attrs.append(tuple(attr))
                preinc(cattr_iter)

            # Again: tuples, so you can't mess with it.
            self._attractors = tuple(attrs)
            self.dirty = False
            return self._attractors

    property rates:
        """Return a readonly numpy array of the rates"""
        def __get__(self):
            # Lazy evaluation
            if self._rates is not None:
                return self._rates

            # Construct the numpy array via python
            r = numpy.zeros((self.factory.cfactory.environments.size(),
                         self.factory.cfactory.out_channels))

            cdef:
                np.npy_double[:,:] c_r = r
                size_t i, j

            # Copy in the goods using a memory array
            for i in range(self.factory.cfactory.environments.size()):
                for j in range(self.factory.cfactory.out_channels):
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

    def testing(self):
        cdef cGene *g = self.cnetwork.genes[0];
        cdef cCisModule *m = g.modules[0];
        print 'one'
        cdef cCisModuleLogic2 *ml = dynamic_cast_cCisModuleLogic2(m)
        print ml.op


cdef class Gene:
    """A proxy to a gene.
    """
    cdef:
        # Assumption: Networks CANNOT mess with genes number once a network
        # has been established (You must copy and mutate a network).
        cGene *cgene
        object _modules

        readonly:
            # By holding this ref, we ensure the pointer is always valid (pace
            # what I said above. DON'T mess with the genes!)
            Network network
            size_t gene_number

    def __cinit__(self, Network n, size_t g):
        """This should never be called publicly"""
        self.network = n
        self.gene_number = g
        self.cgene = self.network.cnetwork.genes[g]

    property sequence:
        def __get__(self):
            return self.cgene.sequence

    property pub:
        def __get__(self):
            return self.cgene.pub

    property module_count:
        def __get__(self):
            return self.cgene.modules.size()

    property modules:
        def __get__(self):
            # Lazy construction
            if self._modules is None:
                self._modules = tuple(CisModuleLogic2(self, i) 
                                      for i in range(self.module_count))
            return self._modules

    def __repr__(self):
        p = self.network.factory.params
        return "<Gene[{}]: {}>".format(self.sequence, p.name_for_channel(self.pub))


cdef class CisModule:
    """A proxy to a CisModule.
    """
    cdef:
        cCisModule *ccismodule
    
        readonly:
            Gene gene

    def __cinit__(self, Gene g, size_t i):
        self.gene = g
        assert i < g.cgene.modules.size()
        self.ccismodule = g.cgene.modules[i]

    def site_count(self):
        return self.ccismodule.site_count()

    def get_site(self, size_t i):
        assert i < self.ccismodule.site_count()
        return self.ccismodule.get_site_channel(i)

    cdef reset_network(self):
        self.gene.network._attractors = None
        self.gene.network._rates = None
        self.gene.network.cnetwork.calc_attractors()

    # property op:
    #     def __get__(self):
    #         return self.ccismodule.op
    #
    # property sub1:
    #     def __get__(self):
    #         return self.ccismodule.sub1
    #     def __set__(self, size_t c):
    #         self.ccismodule.sub1 = c
    #         self.reset_network()
    #
    # property sub2:
    #     def __get__(self):
    #         return self.ccismodule.sub2
    #     def __set__(self, size_t c):
    #         self.ccismodule.sub2 = c
    #         self.reset_network()

    # def __cmp__(self, CisModule other):
    #     cdef int retval 
    #     retval = c_cmp(self.ccismodule.op, other.ccismodule.op)
    #     if retval != 0:
    #         return retval
    #     retval = c_cmp(self.ccismodule.sub1, other.ccismodule.sub1)
    #     if retval != 0:
    #         return retval
    #     return c_cmp(self.ccismodule.sub2, other.ccismodule.sub2)

    # def test(self, unsigned int a, unsigned int b):
    #     return self.ccismodule.test(a, b)

    def is_active(self, ChannelState p):
        return self.ccismodule.is_active(p.cchannel_state)

    
cdef class CisModuleLogic2(CisModule):
    cdef:
        cCisModuleLogic2 *logic2

    def __cinit__(self, Gene g, size_t i):
        self.logic2 = dynamic_cast_cCisModuleLogic2(self.ccismodule)

    property op:
        def __get__(self):
            return self.logic2.op

    property sub1:
        def __get__(self):
            return self.logic2.channels[0]

    property sub2:
        def __get__(self):
            return self.logic2.channels[1]

    def __str__(self):
        p = self.gene.network.factory.params
        return "{}({}, {})".format(
            Operand(self.op).name,
            p.name_for_channel(self.logic2.channels[0]),
            p.name_for_channel(self.logic2.channels[1]),
        )

    def __repr__(self):
        return "<CisModule: {}>".format(self.__str__())

cdef class NetworkAnalysis:
    cdef:
        cNetworkAnalysis *canalysis
        readonly:
            Network network

    def __cinit__(self, Network net):
        self.network = net
        self.canalysis = new cNetworkAnalysis(net.ptr)
        # self.canalysis.find_knockouts()

    def __dealloc__(self):
        del self.canalysis

    # property knockouts:
    #     def __get__(self):
    #         cdef:
    #             cSiteIndex *si
    #             vector[cSiteIndex].iterator i = self.canalysis.knockouts.begin()
    #
    #         k = []
    #         while i != self.canalysis.knockouts.end():
    #             si = &deref(i)
    #             k.append((si.gene(), si.cis(), si.site()))
    #             preinc(i)
            # return k

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

cdef class NetworkCollection:
    cdef:
        readonly: 
            Factory factory
        cNetworkVector cnetworks

    def __cinit__(self, Factory factory):
        self.factory = factory

    def add(self, Network n):
        assert n.factory is self.factory
        self.cnetworks.push_back(n.ptr)

    cdef object get_at(self, size_t i):
        if i >= self.cnetworks.size():
            raise IndexError

        cdef:
            cNetwork_ptr ptr = self.cnetworks[i]
            cNetwork *net = ptr.get()

        # Is there an existing python object?
        # Note: cython automatically increments the reference count when we do
        # this (which is what we want)
        if net.pyobject:
            return <object>(net.pyobject)

        # Nope. So we need to create a new python wrapper object
        n = Network(self.factory)
        n.bind_to(ptr)
        return n

    property size:
        def __get__(self):
            return self.cnetworks.size()

    # def get_max_fitness(self):
    #     return max_fitness(self.cnetworks);

    def __getitem__(self, size_t i):
        return self.get_at(i)

    def __iter__(self):
        for i in range(self.size):
            yield self.get_at(i)

    def __repr__(self):
        return "<NetworkCollection: {}>".format(self.size)

    def mutate(self):
        cdef cIndexes mutated
        self.factory.cgenefactory.mutate_collection(self.cnetworks, mutated)

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

        sm = new cSelectionModel(self.factory.cfactory_ptr)
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

        This is used for testing.
        """
        cdef:
            cSelectionModel *sm
            bint sel
            cIndexes indexes
            cNetworkVector new_networks

        sm = new cSelectionModel(self.factory.cfactory_ptr)

        sel = sm.select(self.cnetworks, deref(target.ctarget), 
                     self.cnetworks.size(), indexes)

        return indexes

