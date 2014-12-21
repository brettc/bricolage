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
            params.cis_count,
        )
        self.cworld.constructor = c

        # Specialise the python classes
        # self.gene_class = Gene
        self.module_class = CisModule


cdef class Network(core_ext.Network):
    def __cinit__(self, World w): # io_init=True):
        cdef grn.cNetwork_ptr ptr = w.cworld.constructor.construct()
        self.bind_to(ptr)

    # def clone(self):
    #     cdef grn.cNetwork_ptr ptr = self.cnetwork.clone()
    #     copy = Network(self.world, do_init=False)
    #     copy.bind_to(ptr)
        
        
# No need?
# cdef class Gene(core_ext.Gene):
#     def __cinit__(self):
#         pass
#
#     def test(self):
#         cdef cGene *g = dynamic_cast_cGene(self.cgene)
#         return g.modules.size()

        
cdef class CisModule(core_ext.CisModule):
    # def __cinit__(self, core_ext.Gene g, size_t i):
    #     pass

    property bindings:
        def __get__(self):
            cdef cCisModule *cm = dynamic_cast_cCisModule(self.ccismodule) 
            return tuple(cm.binding[i] for i in range(3))

    property operation:
        def __get__(self):
            w = self.gene.network.world
            return boolean_func_from_coop_binding(w, self.channels, self.bindings)

    def is_active(self, core_ext.ChannelState c):
        cdef cCisModule *cm = dynamic_cast_cCisModule(self.ccismodule) 
        return cm.is_active(c.cchannel_state)

    def _evil_set_binding(self, size_t i, int w):
        assert i <= 3
        assert -3 <= w <= 3
        cdef cCisModule *cm = dynamic_cast_cCisModule(self.ccismodule) 
        cm.binding[i] = w
        self.gene.network._invalidate_cached()

    def __str__(self):
        return self.operation

    def __repr__(self):
        return "<CisT3: {}>".format(self.__str__())

