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
# from libc.math cimport log2, fabs
# ctypedef np.int_t int_type
ctypedef np.uint8_t np_uint8
# ctypedef np.float64_t float_type
# ctypedef int bitset_type
# ctypedef vector[int_type] attractor_type
# ctypedef vector[int_type].iterator attractor_iterator_type
# ctypedef vector[attractor_type] attractor_individual_type
#
# ctypedef vector[Gene] GeneVector

cdef extern from "pubsub2_c.h" namespace "pubsub2":
    # ctypedef unsigned char uint8
    # ctypedef unsigned int uint32
    # ctypedef float float32
    # ctypedef short int16

    cdef cppclass Gene:
        np_uint8 sub1, sub2, pub

    cdef cppclass Genome:
        Genome()
        void init(size_t size)
        # void fill(int i)
        vector[Gene] genes

cdef class cParameters:
    cdef:
        readonly:
            object params
            size_t gene_count

    def __cinit__(self, params):
        self.params = params
        self.gene_count = params.gene_count

        # Now translate everything to C


cdef class NetworkFactory:
    cdef:
        readonly:
            cParameters cparams

    def __cinit__(self, params):
        self.cparams = cParameters(params)

    def get_network(self):
        n = Network(self.cparams)
        n.ready = True
        return n


cdef class Network:
    cdef:
        readonly:
            cParameters cparams
            bint ready

        Genome genome

    def __cinit__(self, cParameters cparams):
        self.cparams = cparams
        self.ready = False
        self.genome.init(self.cparams.gene_count)

    def export_genes(self):
        output = numpy.zeros((self.cparams.gene_count, 3), dtype=numpy.uint8)

        cdef:
            vector[Gene].iterator gene_i
            Gene *g
            np_uint8[:,:] output_c = output
            size_t i = 0
        
        gene_i = self.genome.genes.begin()

        while gene_i != self.genome.genes.end():
            g = &deref(gene_i)
            output_c[i, 0] = g.sub1
            output_c[i, 1] = g.sub2
            output_c[i, 2] = g.pub
            preinc(gene_i)
            i += 1

        return output

# cdef class Population:
# cdef class NetworkRepr



# cdef class Networks:
#     pass


