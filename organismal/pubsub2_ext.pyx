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


cdef extern from "<boost/shared_ptr.hpp>" namespace "boost":
    cdef cppclass shared_ptr[T]:
        shared_ptr()
        shared_ptr(T* ptr)
        shared_ptr(shared_ptr& r)
        T* get()


cdef extern from "pubsub2_c.h" namespace "pubsub2":
    cdef cppclass cGene:
        np.npy_ubyte sub1, sub2, pub

    cdef cppclass cNetwork:
        cNetwork()
        void init(np.npy_int, size_t size)
        void test()
        vector[cGene] genes
        np.npy_int identifier
        np.npy_int gene_count
        np.npy_byte *gene_data()

    ctypedef shared_ptr[cNetwork] cNetwork_ptr

    cdef cppclass cPopulation:
        vector[cNetwork_ptr] networks


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
        n.create_new(self.next_identifier)
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
        cNetwork_ptr ptr
        cNetwork *cnetwork

    # Networks should only be created by Factory. They have a two-stage
    # creation. Depending on what we're doing, either `create_new` or
    # `create_copy` needs to be called
    def __cinit__(self, Factory factory):
        self.factory = factory
        self.ready = False

    cdef create_new(self, np.npy_int ident):
        self.ptr = cNetwork_ptr(new cNetwork())
        self.cnetwork = self.ptr.get()
        self.cnetwork.init(ident, self.factory.gene_count)
        self.ready = True

    cdef create_copy(self, cNetwork_ptr w):
        self.ptr = w
        self.cnetwork = self.ptr.get()
        self.ready = True

    # Will this crash?
    property identifier:
        def __get__(self):
            return self.cnetwork.identifier

    def gene_array(self):
        cdef:
            np.npy_byte[:,:] copied

        copied = <np.npy_byte[:self.cnetwork.gene_count, :3]>self.cnetwork.gene_data()
        return numpy.asarray(copied)

    def export_genes(self):
        output = numpy.zeros((self.factory.gene_count, 3), dtype=numpy.uint8)

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

    def test(self):
        self.cnetwork.test()

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

    property sub1:
        def __get__(self):
            return self.network.cnetwork.genes[self.gene_number].sub1
        def __set__(self, np.npy_ubyte val):
            self.network.cnetwork.genes[self.gene_number].sub1 = val



# cdef class GeneSub:
#     cdef:
#         readonly:
#             Gene gene
#             np.npy_int gene_number
#
#     def __cinit__(self, Gene g, np.npy_int n):
#         self.gene = g
#         self.sub_number = n
#
#     def __getitem__(self, np.npy_int i):
#         return self.gene.network.


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
        # Note that we add the c++ shared pointer, rather than the python
        # object. The original python Network holding this can go out of
        # scope, and that is fine.
        self.cpop.networks.push_back(n.ptr)

    def get(self, size_t i):
        cdef:
            cNetwork_ptr ptr = self.cpop.networks[i]
            Network n = Network(self.factory)

        # We need to construct a new python Network to wrap the existing
        # cNetwork. 
      
        # NOTE: We could get clever here, and find a way to see if there is
        # already an existing python object for this cNetwork. But the
        # overhead and the pointer fandangling required is not worth it.
        n.create_copy(ptr)
        return n


# cdef class NetworkRepr


# cdef class Networks:
#     pass


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
