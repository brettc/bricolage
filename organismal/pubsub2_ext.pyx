# distutils: include_dirs = NUMPY_PATH
# distutils: language = c++
# distutils: sources = organismal/pubsub2_c.cpp
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

import cython
import numpy
cimport numpy as np
from cython.operator import dereference as deref, preincrement as preinc
from cpython cimport Py_DECREF, Py_INCREF
from libcpp.vector cimport vector
from libcpp.string cimport string

ctypedef np.npy_byte byte_t
ctypedef np.npy_int int_t


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

    cdef cppclass cFactory:
        cFactory(size_t seed)
        mt19937 random_engine
        size_t pop_count, gene_count, cis_count

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
        bint test(byte_t a, byte_t b)
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

    cdef cppclass cPopulation:
        vector[cNetwork_ptr] networks


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
        self.cfactory_ptr = cFactory_ptr(new cFactory(1))
        self.cfactory = self.cfactory_ptr.get()

        # Now translate everything to C
        self.cfactory.gene_count = params.gene_count

    # property gene_count:
    #     def __get__(self):
    #         return self.cfactory.gene_count

    def create_network(self):
        cdef cNetwork_ptr ptr = cNetwork_ptr(new cNetwork(self.cfactory_ptr))
        n = Network(self)
        n.bind_to(ptr)
        return n

    def create_population(self):
        p = Population(self)
        return p

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
            cdef size_t i
            if self._genes is None:
                self._genes = []
                for i in range(self.cnetwork.genes.size()):
                    self._genes.append(Gene(self, i))
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
        # Assumption: Networks don't change their gene size
        cGene *cgene

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

    def __getitem__(self, size_t i):
        return CisModule(self, i)

    def __iter__(self):
        cdef size_t i
        for i in range(self.cgene.modules.size()):
            yield CisModule(self, i)


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
        def __set__(self, np.npy_ubyte val):
            self.ccismodule.sub1 = val

    property sub2:
        def __get__(self):
            return self.ccismodule.sub2
        def __set__(self, np.npy_ubyte val):
            self.ccismodule.sub2 = val

    property op:
        def __get__(self):
            return self.ccismodule.op
        def __set__(self, np.npy_ubyte val):
            self.ccismodule.op = val

    def test(self, np.npy_ubyte a, np.npy_ubyte b):
        return self.ccismodule.test(a, b)

    def active(self, Products p):
        return self.ccismodule.active(p.cproducts)


cdef class Population:
    cdef:
        readonly: 
            bint ready

        Factory factory
        cPopulation cpop

    def __cinit__(self, Factory factory):
        self.factory = factory
        self.ready = False

        cdef size_t i
        for i in range(factory.params.population_size):
            self.cpop.networks.push_back(
                cNetwork_ptr(new cNetwork(factory.cfactory_ptr)))


    def add(self, Network n):
        # Note that we add the c++ shared pointer, rather than the python
        # object. The original python Network holding this can go out of
        # scope, and that is fine.
        #
        # TODO: should check the factory is the same.
        self.cpop.networks.push_back(n.ptr)

    def get(self, size_t i):
        cdef:
            cNetwork_ptr ptr = self.cpop.networks[i]
            cNetwork *net = ptr.get()

        # Is there an existing python object?
        if net.pyobject:
            print 'already have an object, returning it'
            return <object>(net.pyobject)

        # Ok, so we need to create a python wrapper object
        print 'creating a new object'
        n = Network(self.factory)
        n.bind_to(ptr)
        return n

