from utility cimport *
from core cimport *
from selection cimport *
from analysis cimport *
from core_ext cimport *

cdef class BaseTarget:
    cdef:
        cBaseTarget *_base
        readonly:
            World world

cdef class DefaultTarget(BaseTarget):
    cdef:
        cDefaultTarget *_this

cdef class NoisyTarget(BaseTarget):
    cdef:
        cNoisyTarget *_this

