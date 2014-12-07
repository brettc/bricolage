# distutils: include_dirs = NUMPY_PATH
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

from cython.operator import dereference as deref, preincrement as preinc
from grn cimport *
cimport core_ext

cdef class Factory(core_ext.Factory):
    def __cinit__(self, params):

        # Create a mutator
        self.cgenefactory = new cGeneFactoryThreshold3(self.cfactory, 
                                         params.gene_mutation_rate)
        self.cis_class = CisModule

    
cdef class CisModule(core_ext.CisModule):
    cdef:
        cCisModuleThreshold3 *threshold3

    def __cinit__(self, core_ext.Gene g, size_t i):
        self.threshold3 = dynamic_cast_cCisModuleThreshold3(self.ccismodule)

    def __str__(self):
        p = self.gene.network.factory.params
        return ",".join([str(self.threshold3.binding[i]) for i in range(3)])

    def __repr__(self):
        return "<CisModule: {}>".format(self.__str__())

