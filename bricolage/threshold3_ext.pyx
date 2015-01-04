# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

from cython.operator import dereference as deref, preincrement as preinc
# from grn cimport *
from threshold_cpp cimport *
cimport core_ext

from .logic_tools import boolean_func_from_coop_binding

cdef class Constructor(core_ext.Constructor):
    def __cinit__(self, core_ext.World w, params):
        self._shared = grn.cConstructor_ptr(new cConstructor(
            w._shared,
            params.cis_count,
        ))
        self._this = self._shared.get()

        # Specialise the python classes
        self.module_class = CisModule

        
cdef class CisModule(core_ext.CisModule):
    # def __cinit__(self, core_ext.Gene g, size_t i):
    #     pass

    property bindings:
        def __get__(self):
            cdef cCisModule *cm = dynamic_cast_cCisModule(self._this) 
            return tuple(cm.binding[i] for i in range(3))

    property operation:
        def __get__(self):
            w = self.gene.network.constructor.world
            return boolean_func_from_coop_binding(w, self.channels, self.bindings)

    def is_active(self, core_ext.ChannelState c):
        cdef cCisModule *cm = dynamic_cast_cCisModule(self._this) 
        return cm.is_active(c._this)

    def _evil_set_binding(self, size_t i, int w):
        assert i <= 3
        assert -3 <= w <= 3
        cdef cCisModule *cm = dynamic_cast_cCisModule(self._this) 
        cm.binding[i] = w
        self.gene.network._invalidate_cached()

    def __str__(self):
        return self.operation

    def __repr__(self):
        return "<CisT3: {}>".format(self.__str__())

