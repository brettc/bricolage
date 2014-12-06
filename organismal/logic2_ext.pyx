# distutils: include_dirs = NUMPY_PATH
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

from .operand import Operand
from cython.operator import dereference as deref, preincrement as preinc
from grn cimport *
cimport core_ext

cdef class Factory(core_ext.Factory):
    def __cinit__(self, params):

        # Create a mutator
        self.cgenefactory = new cGeneFactoryLogic2(self.cfactory, 
                                         params.gene_mutation_rate)
        self.cis_class = CisModule

    
cdef class CisModule(core_ext.CisModule):
    cdef:
        cCisModuleLogic2 *logic2

    def __cinit__(self, core_ext.Gene g, size_t i):
        self.logic2 = dynamic_cast_cCisModuleLogic2(self.ccismodule)

    property op:
        def __get__(self):
            return self.logic2.op

    property sub1:
        def __get__(self):
            return self.logic2.channels[0]

    property sub2:
        def __get__(self):
            return self.logic2.channels[1]

    def __str__(self):
        p = self.gene.network.factory.params
        return "{}({}, {})".format(
            Operand(self.op).name,
            p.name_for_channel(self.logic2.channels[0]),
            p.name_for_channel(self.logic2.channels[1]),
        )

    def __repr__(self):
        return "<CisModule: {}>".format(self.__str__())

