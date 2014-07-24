# distutils: include_dirs = NUMPY_PATH
# distutils: language = c++
# distutils: sources = organismal/pubsub2_c.cpp
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

import numpy
cimport numpy as np
from cython.operator import dereference as deref, preincrement as preinc
from libcpp.vector cimport vector

# from numpy.pxd
# ctypedef signed char      npy_byte
# ctypedef signed int       npy_int
# ctypedef signed long      npy_long
# ctypedef signed long long npy_longlong
#
# ctypedef unsigned char      npy_ubyte
# ctypedef unsigned long      npy_ulong
# ctypedef unsigned long long npy_ulonglong
#
# ctypedef float        npy_float
# ctypedef double       npy_double
# ctypedef long double  npy_longdouble
# ctypedef Py_intptr_t npy_intp
# ctypedef size_t npy_uintp


cdef extern from "pubsub2_c.h" namespace "pubsub2":
    cdef cppclass cGene:
        np.npy_ubyte sub1, sub2, pub

    cdef cppclass cNetwork:
        cNetwork()
        void init(np.npy_int, size_t size)
        vector[cGene] genes
        np.npy_int identifier

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

        np.npy_int next_identifier

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

    property identifier:
        def __get__(self):
            return self.cnetwork.identifier

    cdef create(self, np.npy_int ident):
        self.cnetwork_ptr = cNetworkSharedPtr(new cNetwork())
        self.cnetwork = self.cnetwork_ptr.get()
        self.cnetwork.init(ident, self.factory.gene_count)
        self.ready = True

    cdef copy(self, cNetworkSharedPtr w):
        self.cnetwork_ptr = w
        self.cnetwork = self.cnetwork_ptr.get()
        self.ready = True

    def export_genes(self):
        charoutput = numpy.zeros((self.factory.gene_count, 3), dtype=numpy.uint8)

        cdef:
            cGene *g
            np.npy_ubyte[:,:] output_c = output
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

    def import_genes(self, np.npy_ubyte[:, :] input_c):

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
        cdef np.npy_int index = i
        return Gene(self, index)

    def __iter__(self):
        count = self.factory.gene_count
        for i in range(count):
            yield Gene(self, i)


cdef class Gene:
    """A proxy to a gene. Clumsy but useful.
    """
    cdef:
        readonly:
            Network network
            np.npy_int gene_number

    def __cinit__(self, Network n, np.npy_int g):
        self.network = n
        self.gene_number = g

    property pub:
        def __get__(self):
            return self.network.cnetwork.genes[self.gene_number].pub
        def __set__(self, np.npy_ubyte val):
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


