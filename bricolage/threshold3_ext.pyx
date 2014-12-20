# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

from cython.operator import dereference as deref, preincrement as preinc
# from grn cimport *
from threshold_cpp cimport *
cimport core_ext

from .logic_tools import boolean_func_from_coop_binding

cdef class World(core_ext.World):
    def __cinit__(self, params):
        # Create a mutator
        cdef cConstructor *c = new cConstructor(
            self.cworld_ptr,
            params.gene_count,
            params.cis_count,
        )
        self.cworld.constructor = c
        # print '----'
        # for i in range(10):
        #     print c.xxx()
        # print '----'
        self.gene_class = Gene
        self.module_class = CisModule


cdef class Network(core_ext.Network):
    def __cinit__(self, World w):
        cdef grn.cNetwork_ptr ptr = grn.cNetwork_ptr(
            new cNetwork(deref(w.cworld.constructor)))
        self.bind_to(ptr)

        
cdef class Gene(core_ext.Gene):
    def __cinit__(self):
        pass

    def test(self):
        cdef cGene *g = dynamic_cast_cGene(self.cgene)
        return g.modules.size()

        
cdef class CisModule(core_ext.CisModule):
    # cdef:
    #     cCisModuleThreshold3 *threshold3

    def __cinit__(self, Gene g, size_t i):
        pass
#
#     property bindings:
#         def __get__(self):
#             return tuple(self.threshold3.binding[i] for i in range(3))
#
#     property operation:
#         def __get__(self):
#             factory = self.gene.network.factory
#             return boolean_func_from_coop_binding(factory, self.channels, self.bindings)
#
#     def __str__(self):
#         return self.operation
#         # p = self.gene.network.factory.params
#         # names = [p.name_for_channel(c) for c in self.channels]
#         # return "+".join(["{}*{}".format(*_) for _ in zip(names, self.bindings)])
#
#     def __repr__(self):
#         return "<threshold3.CisModule: {}>".format(self.__str__())
#
