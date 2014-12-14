# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

from cython.operator import dereference as deref, preincrement as preinc
from grn cimport *
cimport core_ext

from .logic_tools import boolean_func_from_coop_binding

cdef class Factory(core_ext.Factory):
    def __cinit__(self, params):
        # Create a mutator
        self.cfactory.constructor = new cConstructorThreshold3(
            deref(self.cfactory),
            params.gene_count,
            params.cis_count,
        )
        self.cis_class = CisModule

    
cdef class CisModule(core_ext.CisModule):
    cdef:
        cCisModuleThreshold3 *threshold3

    def __cinit__(self, core_ext.Gene g, size_t i):
        self.threshold3 = dynamic_cast_cCisModuleThreshold3(self.ccismodule)

    property bindings:
        def __get__(self):
            return tuple(self.threshold3.binding[i] for i in range(3))

    property operation:
        def __get__(self):
            params = self.gene.network.factory.params
            return boolean_func_from_coop_binding(params, self.channels, self.bindings)

    def __str__(self):
        return self.operation
        # p = self.gene.network.factory.params
        # names = [p.name_for_channel(c) for c in self.channels]
        # return "+".join(["{}*{}".format(*_) for _ in zip(names, self.bindings)])

    def __repr__(self):
        return "<threshold3.CisModule: {}>".format(self.__str__())

