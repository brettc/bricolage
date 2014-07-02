# distutils: include_dirs = NUMPY_PATH
# distutils: language = c++
# distutils: sources = organismal/func.cpp
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

import numpy
cimport numpy as np
from cython.operator import dereference as deref, preincrement as preinc

# from libcpp.vector cimport vector
# from libc.math cimport log2, fabs
# ctypedef np.int_t int_type
# ctypedef np.int8_t tiny_type
# ctypedef np.float64_t float_type
# ctypedef int bitset_type
# ctypedef vector[int_type] attractor_type
# ctypedef vector[int_type].iterator attractor_iterator_type
# ctypedef vector[attractor_type] attractor_individual_type
#

# cdef extern from "func.h" namespace "func":
#     cdef cppclass tester:
#         int i, j, k
#         int bob()

cdef class NetworkFactory:
    cdef:
        object params
        int population_size

    def __cinit__(self, params):
        self.params = params
        self.population_size = params.population_size
        # Translate some values

    def get_network(self):
        n = Network(self.params)


cdef class Network:
    cdef:
        object params
        readonly bint not_ready

    def __cinit__(self, params):
        self.params = params
        self.not_ready = True


# cdef class Networks:
#     pass


