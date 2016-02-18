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

    property modified: 
        def __get__(self):
            cdef cNetwork *ptr = self._this.modified.get()

            # Is there an existing python object?
            # Note: cython automatically increments the reference count when we do
            # this (which is what we want). This move just means we get object
            # identity, at least while one python reference continues to exist.
            if ptr.pyobject:
                return <object>(ptr.pyobject)

            # Nope. We need to create a new python wrapper object
            cdef Factory f = self.network.factory
            cdef Network n = f.network_class(f, f._secret_key)
            n.bind_to(self._this.modified)
            return n

    property potential_bindings: 
        def __get__(self):
            return self._this.potential_bindings

    property active_bindings: 
        def __get__(self):
            return self._this.active_bindings

    property world: 
        def __get__(self):
            return self.modified.factory.world

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
    def __cinit__(self, JointProbabilities joint=None):
        if joint is None:
            # Should be bound later!
            self._this = NULL
        else:
            self._this = new cInformation(deref(joint._this))

    cdef bind(self, cInformation *i):
        self._this = i

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
    def __cinit__(self, World w):
        self.world = w
        self._this = new cCausalFlowAnalyzer(w._shared)

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

    def numpy_info_from_collection(self, CollectionBase coll):
        i = Information(self.analyse_collection(coll))
        return numpy.asarray(i)

    def numpy_info_from_network(self, Network n):
        i = Information(self.analyse_network(n))
        return numpy.asarray(i)


cdef class AverageControlAnalyzer:
    def __cinit__(self, World w):
        self.world = w
        self._this = new cAverageControlAnalyzer(w._shared)

    def __dealloc__(self):
        if self._this != NULL:
            del self._this

    def analyse_network(self, Network n):
        info = Information()
        cdef cInformation *c_info= NULL 
        c_info = self._this.analyse_network(deref(n._this))
        info.bind(c_info)
        return info

    def analyse_collection(self, CollectionBase coll):
        info = Information()
        cdef cInformation *c_info= NULL 
        c_info = self._this.analyse_collection(deref(coll._collection))
        info.bind(c_info)
        return info

    def numpy_info_from_collection(self, CollectionBase coll):
        i = self.analyse_collection(coll)
        return numpy.asarray(i)

    def numpy_info_from_network(self, Network n):
        i = self.analyse_network(n)
        return numpy.asarray(i)


cdef class MutualInfoAnalyzer:
    def __cinit__(self, World w, cIndexes categories):
        self.world = w
        assert categories.size() == self.world._this.environments.size()

        catset = set(list(categories))
        ncats = len(catset)
        if ncats < 2:
            raise ValueError("There must be at least two categories")

        self._this = new cMutualInfoAnalyzer(w._shared, categories)

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

    def numpy_info_from_network(self, Network n):
        i = Information(self.analyse_network(n))
        return numpy.asarray(i)

    def numpy_info_from_collection(self, CollectionBase coll):
        i = Information(self.analyse_collection(coll))
        return numpy.asarray(i)

cdef class OutputControlAnalyzer:
    def __cinit__(self, World w, cRatesVector tr = []):
        self.world = w
        self._this = new cOutputControlAnalyzer(w._shared, tr)

    def __dealloc__(self):
        if self._this != NULL:
            del self._this

    def analyse_network(self, Network n):
        info = Information()
        cdef cInformation *c_info= NULL 
        c_info = self._this.analyse_network(deref(n._this))
        info.bind(c_info)
        return info

    def analyse_collection(self, CollectionBase coll):
        info = Information()
        cdef cInformation *c_info= NULL 
        c_info = self._this.analyse_collection(deref(coll._collection))
        info.bind(c_info)
        return info

    def numpy_info_from_collection(self, CollectionBase coll):
        i = self.analyse_collection(coll)
        return numpy.asarray(i)

    def numpy_info_from_network(self, Network n):
        i = self.analyse_network(n)
        return numpy.asarray(i)


cdef class RelevantControlAnalyzer:
    def __cinit__(self, World w, cRatesVector tr = []):
        self.world = w
        self._this = new cRelevantControlAnalyzer(w._shared, tr)

    def __dealloc__(self):
        if self._this != NULL:
            del self._this

    def analyse_network(self, Network n):
        info = Information()
        cdef cInformation *c_info= NULL 
        c_info = self._this.analyse_network(deref(n._this))
        info.bind(c_info)
        return info

    def analyse_collection(self, CollectionBase coll):
        info = Information()
        cdef cInformation *c_info= NULL 
        c_info = self._this.analyse_collection(deref(coll._collection))
        info.bind(c_info)
        return info

    def numpy_info_from_collection(self, CollectionBase coll):
        i = self.analyse_collection(coll)
        arr = numpy.asarray(i)
        arr.shape = arr.shape[:2]
        return arr

    def numpy_info_from_network(self, Network n):
        i = self.analyse_network(n)
        arr = numpy.asarray(i)
        arr.shape = arr.shape[1],
        return arr

def _get_max_category_size():
    return get_max_category_size()

def _set_max_category_size(size_t m):
    set_max_category_size(m)
