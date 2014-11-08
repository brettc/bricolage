"""
Definitions for c++ or boost classes
"""
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.utility cimport pair


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

