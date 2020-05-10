"""
Definitions for c++ or boost classes
"""
from libcpp.vector cimport vector
from libcpp.set cimport set as std_set
from libcpp.map cimport map as std_map
from libcpp.string cimport string
from libcpp.utility cimport pair
from libcpp.cast import dynamic_cast, static_cast


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
