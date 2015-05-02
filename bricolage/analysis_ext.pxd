from utility cimport *
from core cimport *
from core_ext cimport *

cdef class NetworkAnalysis:
    cdef:
        cNetworkAnalysis *_this
        readonly:
            Network network

cdef class InfoE:
    cdef:
        cInfoE *_this
        readonly:
            World world

cdef class JointProbabilities:
    cdef:
        readonly:
            size_t size
        cJointProbabilities *_this
        vector[Py_ssize_t] shape
        vector[Py_ssize_t] strides
        World world
    cdef bind(self, cJointProbabilities *j)

cdef class Information:
    cdef:
        readonly:
            size_t size
        cInformation *_this
        vector[Py_ssize_t] shape
        vector[Py_ssize_t] strides
        World world

cdef class CausalFlowAnalyzer:
    cdef:
        cCausalFlowAnalyzer *_this;
        World world
