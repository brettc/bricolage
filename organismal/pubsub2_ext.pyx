# distutils: include_dirs = NUMPY_PATH
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

import cython
import numpy
cimport numpy as np
from cython.operator import dereference as deref, preincrement as preinc
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.utility cimport pair

from operand import Operand

cdef extern from "<random>" namespace "std":
    cdef cppclass mt19937:
        # mt19937(size_t seed)
        void seed(size_t s)
        unsigned int operator()()

    cdef cppclass uniform_int_distribution[T]:
        uniform_int_distribution(T, T)
        T operator()(mt19937)


cdef extern from "<boost/shared_ptr.hpp>" namespace "boost":
    cdef cppclass shared_ptr[T]:
        shared_ptr()
        shared_ptr(T* ptr)
        shared_ptr(shared_ptr& r)
        T* get()


cdef extern from "<boost/dynamic_bitset.hpp>" namespace "boost":
    cdef cppclass dynamic_bitset[T]:
        dynamic_bitset()
        void resize(size_t)
        size_t size()
        bint test(size_t)
        void set(size_t)
        void reset(size_t)
        void flip(size_t)

    cdef void to_string(dynamic_bitset[size_t], string s)


cdef extern from "pubsub2_c.h" namespace "pubsub2":
    ctypedef unsigned int signal_t
    ctypedef unsigned int operand_t
    ctypedef unsigned int sequence_t

    ctypedef vector[operand_t] cOperands

    cdef cppclass cFactory:
        cFactory(size_t seed)
        mt19937 random_engine
        sequence_t next_identifier
        size_t pop_count, gene_count, cis_count
        cOperands operands
        pair[size_t, size_t] sub_range
        pair[size_t, size_t] pub_range

    ctypedef shared_ptr[cFactory] cFactory_ptr

    ctypedef vector[cCisModule] cCisModules

    ctypedef dynamic_bitset[size_t] cProducts

    cdef cppclass cProductsSequence:
        ProductStates()
        void init(size_t np)
        void push_back(dynamic_bitset[size_t] p)
        size_t size() 
        size_t products_size() 
        bint get(size_t i, size_t j)
        void set(size_t i, size_t j, bint b)

    cdef cppclass cCisModule:
        bint test(unsigned int a, unsigned int b)
        bint active(dynamic_bitset[size_t] s)
        signal_t op, sub1, sub2

    ctypedef vector[cCisModule] cCisModules

    cdef cppclass cGene:
        sequence_t sequence;
        cCisModules modules;
        signal_t pub

    cdef cppclass cNetwork:
        cNetwork(cFactory_ptr)
        void *pyobject
        vector[cGene] genes
        sequence_t identifier
        size_t gene_count

    ctypedef shared_ptr[cNetwork] cNetwork_ptr
    ctypedef vector[cNetwork_ptr] cNetworkVector


cdef class Products:
    cdef cProducts cproducts;

    def __cinit__(self, size_t size):
        self.cproducts.resize(size)

    def set(self, size_t i):
        self.cproducts.set(i)

    def reset(self, size_t i):
        self.cproducts.reset(i)

    def flip(self, size_t i):
        self.cproducts.flip(i)

    def test(self, size_t i):
        return self.cproducts.test(i)

    def __str__(self):
        cdef string s
        to_string(self.cproducts, s)
        return s

    property size:
        def __get__(self):
            return self.cproducts.size()


cdef class ProductStates:
    cdef cProductsSequence cstates

    def __cinit__(self):
        pass

    cdef init_from(self, cProductsSequence *s):
        self.cstates = deref(s)


cdef class Factory:
    cdef:
        cFactory_ptr cfactory_ptr
        cFactory * cfactory
        readonly:
            object params

    def __cinit__(self, params):
        self.params = params
        self.cfactory_ptr = cFactory_ptr(new cFactory(params.seed))
        self.cfactory = self.cfactory_ptr.get()

        # Now translate everything to cpp 
        # There is automatic conversion for most of this (stl::vector etc)
        self.cfactory.gene_count = params.gene_count
        self.cfactory.cis_count = params.cis_count
        self.cfactory.operands = params.operands
        self.cfactory.sub_range = params.sub_range
        self.cfactory.pub_range = params.pub_range

    def create_network(self):
        cdef cNetwork_ptr ptr = cNetwork_ptr(new cNetwork(self.cfactory_ptr))
        n = Network(self)
        n.bind_to(ptr)
        return n

    def test_states(self):
        self.states.init(5)
        # self.states.push_back()
        # self.states.push_back()
        # self.states.push_back()
        # self.states.set(0, 3, True)
        # self.states.set(1, 1, True)
        # cdef size_t i, j
        # for i in range(self.states.size()):
        #     print '--'
        #     for j in range(self.states.products_size()):
        #         print self.states.get(i, j)


cdef class Network:
    cdef:
        readonly:
            Factory factory
            bint ready

        # Because we hold a reference to the shared_ptr, we know we can always
        # safely access the ACTUAL pointer. We keep a reference to it too, as
        # it makes our life easier. The cost is a tiny bit of space.
        cNetwork_ptr ptr
        cNetwork *cnetwork
        object _genes

    # Networks take two stages. You need to create the python object and then
    # bind it to a cNetwork_ptr.
    def __cinit__(self, Factory factory):
        self.factory = factory
        self.ready = False
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

    property gene_count:
        def __get__(self):
            return self.cnetwork.genes.size()

    property genes:
        def __get__(self):
            if self._genes is None:
                self._genes = tuple(Gene(self, i) for i in range(self.gene_count))
            return self._genes

    # def export_genes(self):
    #     output = numpy.zeros((self.factory.gene_count, 3), dtype=numpy.uint8)
    #
    #     cdef:
    #         cGene *g
    #         np.npy_ubyte[:,:] output_c = output
    #         size_t i = 0
    #         vector[cGene].iterator gene_i = self.cnetwork.genes.begin()
    #
    #     while gene_i != self.cnetwork.genes.end():
    #         g = &deref(gene_i)
    #         output_c[i, 0] = g.sub1
    #         output_c[i, 1] = g.sub2
    #         output_c[i, 2] = g.pub
    #         preinc(gene_i)
    #         i += 1
    #
    #     return output
    #
    # def import_genes(self, np.npy_ubyte[:, :] input_c):
    #
    #     cdef:
    #         vector[cGene].iterator gene_i
    #         cGene *g
    #         size_t i = 0
    #     
    #     gene_i = self.cnetwork.genes.begin()
    #
    #     while gene_i != self.cnetwork.genes.end():
    #         g = &deref(gene_i)
    #         g.sub1 = input_c[i, 0]
    #         g.sub2 = input_c[i, 1]
    #         g.pub = input_c[i, 2]
    #         preinc(gene_i)
    #         i += 1


cdef class Gene:
    """A proxy to a gene.
    """
    cdef:
        # Assumption: Networks CANNOT mess with their genes once established
        # (You must copy and mutate a network)
        cGene *cgene
        object _modules

        readonly:
            Network network
            size_t gene_number

    def __cinit__(self, Network n, size_t g):
        self.network = n
        self.gene_number = g
        self.cgene = &self.network.cnetwork.genes[g]

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
                self._modules = tuple(CisModule(self, i) for i in range(self.module_count))
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
        self.ccismodule = &g.cgene.modules[i]

    property sub1:
        def __get__(self):
            return self.ccismodule.sub1
        def __set__(self, signal_t val):
            self.ccismodule.sub1 = val

    property sub2:
        def __get__(self):
            return self.ccismodule.sub2
        def __set__(self, signal_t val):
            self.ccismodule.sub2 = val

    property op:
        def __get__(self):
            return self.ccismodule.op
        def __set__(self, operand_t val):
            self.ccismodule.op = val

    def test(self, unsigned int a, unsigned int b):
        return self.ccismodule.test(a, b)

    def active(self, Products p):
        return self.ccismodule.active(p.cproducts)

    def __repr__(self):
        p = self.gene.network.factory.params
        return "<CisModule: {}, {}, {}>".format(
            Operand(self.op).name,
            p.name_for_channel(self.sub1),
            p.name_for_channel(self.sub2),
        )
    

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

        # Ok, so we need to create a python wrapper object
        n = Network(self.factory)
        n.bind_to(ptr)
        return n

    property size:
        def __get__(self):
            return self.cnetworks.size()

    def __getitem__(self, size_t i):
        return self.get_at(i)

    def __iter__(self):
        for i in range(self.size):
            yield self.get_at(i)

    def __repr__(self):
        return "<NetworkCollection: {}>".format(self.size)

def make_population(NetworkCollection nc):
    cdef size_t i
    for i in range(nc.factory.params.population_size):
        nc.cnetworks.push_back(
            cNetwork_ptr(new cNetwork(nc.factory.cfactory_ptr)))

