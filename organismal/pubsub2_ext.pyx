# distutils: include_dirs = NUMPY_PATH
# distutils: language = c++
# distutils: sources = organismal/pubsub2_c.cpp
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

# from cpython cimport PyObject, Py_DECREF, Py_INCREF
import numpy
cimport numpy as np
from cython.operator import dereference as deref, preincrement as preinc

from libcpp.vector cimport vector
# from libc.math cimport log2, fabs
# ctypedef np.int_t int_type
ctypedef np.int_t int_type
ctypedef np.uint8_t np_uint8
# ctypedef np.float64_t float_type
# ctypedef int bitset_type
# ctypedef vector[int_type] attractor_type
# ctypedef vector[int_type].iterator attractor_iterator_type
# ctypedef vector[attractor_type] attractor_individual_type
#
# ctypedef vector[cGene] GeneVector

cdef extern from "pubsub2_c.h" namespace "pubsub2":
    cdef cppclass cGene:
        np_uint8 sub1, sub2, pub

    cdef cppclass cNetwork:
        cNetwork()
        void init(int_type ident, size_t size)
        vector[cGene] genes

    cdef cppclass cNetworkSharedPtr:
        cNetworkSharedPtr()
        cNetworkSharedPtr(cNetwork *g)
        cNetwork *get()

    cdef cppclass cPopulation:
        vector[cNetworkSharedPtr] networks


cdef class Factory:
    cdef:
        readonly:
            object params
            size_t gene_count

        int_type next_identifier

    def __cinit__(self, params):
        self.params = params

        # Now translate everything to C
        self.gene_count = params.gene_count
        self.next_identifier = 0

    def create_network(self):
        n = Network(self)
        n.create(self.next_identifier)
        self.next_identifier += 1
        return n

    def create_population(self):
        p = Population(self)
        return p


cdef class Network:
    cdef:
        readonly:
            Factory factory
            bint ready

        # Because we hold a reference to the shared_ptr, we know we can always
        # safely access the ACTUAL pointer. We keep a reference to it too, as
        # it makes our life easier. The cost is a tiny bit of space.
        cNetworkSharedPtr cnetwork_ptr
        cNetwork *cnetwork

    def __cinit__(self, Factory factory):
        self.factory = factory
        self.ready = False

    cdef create(self, int_type ident):
        self.cnetwork_ptr = cNetworkSharedPtr(new cNetwork())
        self.cnetwork = self.cnetwork_ptr.get()
        self.cnetwork.init(ident, self.factory.gene_count)
        self.ready = True

    cdef copy(self, cNetworkSharedPtr w):
        self.cnetwork_ptr = w
        self.cnetwork = self.cnetwork_ptr.get()
        self.ready = True

    def export_genes(self):
        output = numpy.zeros((self.factory.gene_count, 3), dtype=numpy.uint8)

        cdef:
            cGene *g
            np_uint8[:,:] output_c = output
            size_t i = 0
            vector[cGene].iterator gene_i = self.cnetwork.genes.begin()

        while gene_i != self.cnetwork.genes.end():
            g = &deref(gene_i)
            output_c[i, 0] = g.sub1
            output_c[i, 1] = g.sub2
            output_c[i, 2] = g.pub
            preinc(gene_i)
            i += 1

        return output

    def import_genes(self, np_uint8[:, :] input_c):

        cdef:
            vector[cGene].iterator gene_i
            cGene *g
            size_t i = 0
        
        gene_i = self.cnetwork.genes.begin()

        while gene_i != self.cnetwork.genes.end():
            g = &deref(gene_i)
            g.sub1 = input_c[i, 0]
            g.sub2 = input_c[i, 1]
            g.pub = input_c[i, 2]
            preinc(gene_i)
            i += 1

    def __getitem__(self, i):
        cdef int_type index = i
        return Gene(self, index)


cdef class Gene:
    """A proxy to a gene. Clumsy but useful.
    """
    cdef:
        readonly:
            Network network
            int_type gene_number

    def __cinit__(self, Network n, int_type g):
        self.network = n
        self.gene_number = g

    property pub:
        def __get__(self):
            return self.network.cnetwork.genes[self.gene_number].pub
        def __set__(self, np_uint8 val):
            self.network.cnetwork.genes[self.gene_number].pub = val


cdef class Population:
    cdef:
        readonly: 
            bint ready

        Factory factory
        cPopulation cpop

    def __cinit__(self, Factory factory):
        self.factory = factory
        self.ready = False

    def add(self, Network n):
        self.cpop.networks.push_back(n.cnetwork_ptr)

    def get(self, size_t i):
        cdef:
            cNetworkSharedPtr ptr = self.cpop.networks[i]
            Network n = Network(self.factory)

        n.copy(ptr)
        return n


# cdef class NetworkRepr



# cdef class Networks:
#     pass


