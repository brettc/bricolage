from utility cimport *
from core cimport *
from targets cimport *
from analysis cimport *
from core_ext cimport *


cdef class NetworkAnalysis:
    cdef:
        cNetworkAnalysis *_this
        readonly:
            Network network


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
    cdef bind(self, cInformation *i)


cdef class CausalFlowAnalyzer:
    cdef:
        cCausalFlowAnalyzer *_this;
        World world


cdef class AverageControlAnalyzer:
    cdef:
        cAverageControlAnalyzer *_this;
        World world


cdef class MutualInfoAnalyzer:
    cdef:
        cMutualInfoAnalyzer *_this;
        World world


cdef class OutputControlAnalyzer:
    cdef:
        cOutputControlAnalyzer *_this;
        World world


cdef class RelevantControlAnalyzer:
    cdef:
        cRelevantControlAnalyzer *_this;
        World world


cdef class MIAnalyzer:
    cdef:
        cMIAnalyzer *_this;
        World world

cdef class WCAnalyzer:
    cdef:
        cWCAnalyzer *_this;
        World world

cdef class FastCAnalyzer:
    cdef:
        cFastCAnalyzer *_this;
        World world
