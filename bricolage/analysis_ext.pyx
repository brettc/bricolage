# distutils: include_dirs = NUMPY_PATH
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# distutils: define_macros=NPY_NO_DEPRECATED_API

import cython
import numpy

# cimports
cimport numpy as np
import numpy
from cython.operator import dereference as deref, preincrement as preinc

cdef class NetworkAnalysis:
    def __cinit__(self, Network net):
        self.network = net
        self._this = new cNetworkAnalysis(net._shared)

    def __dealloc__(self):
        del self._this

    def get_edges(self):
        cdef:
            cEdgeList edges
        self._this.make_edges(edges)
        return edges

    def get_active_edges(self):
        cdef:
            cEdgeList edges
        self._this.make_active_edges(edges)
        return edges

cdef class InfoE:
    def __cinit__(self, World w, init_func):
        self.world = w
        a, b = self.world._this.cue_range
        cats = []
        for i, e in enumerate(self.world.environments):
            cats.append(init_func(*e.as_array()[a:b]))

        # Make sure the categories are consecutive 0, 1, 2 ...
        catset = set(cats)
        ncats = len(catset)
        if catset != set(range(ncats)):
            raise ValueError("Categories must be consecutively numbered"
                             " integers from 0 to N")
        # and there are at least 2
        if ncats < 2:
            raise ValueError("There must be at least two categories")

        self._this = new cInfoE(w._shared, ncats)
        self._this.categories = cats

    property categories:
        def __get__(self):
            return self._this.categories

    def network_probs(self, Network net):
        cdef:
            size_t b=0, c=0, d=0

        self._this.get_extents(b, c, d)
        np_arr = numpy.zeros((b, c, d))

        cdef double [:, :, :] np_arr_view = np_arr
        self._this.network_probs(&np_arr_view[0][0][0],
                                   deref(net._this))
        return np_arr

    def collection_probs(self, CollectionBase coll):
        cdef:
            size_t a, b=0, c=0, d=0

        a = coll._collection.size()
        self._this.get_extents(b, c, d)
        np_arr = numpy.zeros((a, b, c, d))

        cdef double [:, :, :, :] np_arr_view = np_arr
        self._this.collection_probs(&np_arr_view[0][0][0][0],
                                   deref(coll._collection))
        return np_arr

    def collection_info(self, CollectionBase coll):
        cdef size_t a, b=0, c=0, d=0
        a = coll._collection.size()
        self._this.get_extents(b, c, d)
        np_arr = numpy.zeros((a, b))
        cdef double [:, :] np_arr_view = np_arr
        self._this.collection_info(&np_arr_view[0][0], deref(coll._collection))
        return np_arr

    def __dealloc__(self):
        del self._this

cdef class JointProbabilities:
    def __cinit__(self, World w):
        self.world = w
        # self.size = size
        self._this = NULL

    cdef bind(self, cJointProbabilities *j):
        self._this = j

    def __dealloc__(self):
        if self._this != NULL:
            del self._this

    def __getbuffer__(self, Py_buffer *buffer, int flags):
        if self._this == NULL:
            raise BufferError

        # We're just copying stuff here. This is not expensive and avoids
        # trying to access the const * items inside multi_array. Also, the
        # strides here are in bytes, not elements sizes as they are in
        # multi-array
        cdef size_t i
        for i in range(self._this.dimensions()):
            self.shape.push_back(self._this.shape_n(i))
            self.strides.push_back(self._this.stride_n(i) * self._this.element_size())

        buffer.buf = self._this.data()
        # You must set this to the right thing. Look in python "struct"
        # modules for the flags
        buffer.format = 'd'
        buffer.internal = NULL
        buffer.itemsize = self._this.element_size()
        buffer.len = self._this.total_size()
        buffer.ndim = self._this.dimensions()
        buffer.obj = self
        buffer.readonly = 0
        buffer.shape = self.shape.data()
        buffer.strides = self.strides.data()
        buffer.suboffsets = NULL

    def __releasebuffer__(self, Py_buffer *buffer):
        pass

cdef class Information:
    def __cinit__(self, JointProbabilities joint):
        self._this = new cInformation(deref(joint._this))

    def __dealloc__(self):
        if self._this != NULL:
            del self._this

    def __getbuffer__(self, Py_buffer *buffer, int flags):
        # We're just copying stuff here. This is not expensive and avoids
        # trying to access the const * items inside multi_array. Also, the
        # strides here are in bytes, not elements sizes as they are in
        # multi-array
        cdef size_t i
        for i in range(self._this.dimensions()):
            self.shape.push_back(self._this.shape_n(i))
            self.strides.push_back(self._this.stride_n(i) * self._this.element_size())

        buffer.buf = self._this.data()
        # You must set this to the right thing. Look in python "struct"
        # modules for the flags
        buffer.format = 'd'
        buffer.internal = NULL
        buffer.itemsize = self._this.element_size()
        buffer.len = self._this.total_size()
        buffer.ndim = self._this.dimensions()
        buffer.obj = self
        buffer.readonly = 0
        buffer.shape = self.shape.data()
        buffer.strides = self.strides.data()
        buffer.suboffsets = NULL

    def __releasebuffer__(self, Py_buffer *buffer):
        pass


cdef class CausalFlowAnalyzer:
    def __cinit__(self, World w, cRates rates):
        assert len(rates) == w.out_channels
        self.world = w
        self._this = new cCausalFlowAnalyzer(w._shared, rates)

    def __dealloc__(self):
        if self._this != NULL:
            del self._this

    def analyse_network(self, Network n):
        j = JointProbabilities(self.world)
        cdef cJointProbabilities *c_joint = NULL 
        c_joint = self._this.analyse_network(deref(n._this))
        j.bind(c_joint)
        return j

    def analyse_collection(self, CollectionBase coll):
        j = JointProbabilities(self.world)
        cdef cJointProbabilities *c_joint = NULL 
        c_joint = self._this.analyse_collection(deref(coll._collection))
        j.bind(c_joint)
        return j
