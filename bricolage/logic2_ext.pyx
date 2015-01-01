# distutils: include_dirs = NUMPY_PATH
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

from .operand import Operand
from logic2_cpp cimport *
cimport core_ext

cdef class Constructor(core_ext.Constructor):
    def __cinit__(self, core_ext.World w, params):
        self._shared = grn.cConstructor_ptr(new cConstructor(
            w.cworld_ptr,
            params.cis_count,
            params.operands,
        ))
        self._this = self._shared.get()

        # Specialise the python classes
        self.module_class = CisModule

    property gene_count:
        def __get__(self):
            cdef cConstructor *c = dynamic_cast_cConstructor(self._this) 
            return c.gene_count
    
    property module_count:
        def __get__(self):
            cdef cConstructor *c = dynamic_cast_cConstructor(self._this) 
            return c.module_count

    property operands:
        def __get__(self):
            cdef cConstructor *c = dynamic_cast_cConstructor(self._this) 
            return [Operand(_) for _ in c.operands]
    
cdef class CisModule(core_ext.CisModule):
    property op:
        def __get__(self):
            cdef cCisModule *cm = dynamic_cast_cCisModule(self.ccismodule) 
            return Operand(cm.op)

    property sub1:
        def __get__(self):
            return self.ccismodule.channels[0]

    property sub2:
        def __get__(self):
            return self.ccismodule.channels[1]

    def __str__(self):
        w = self.gene.network.constructor.world
        return "{}({}, {})".format(
            self.op.name,
            w.name_for_channel(self.ccismodule.channels[0]),
            w.name_for_channel(self.ccismodule.channels[1]),
        )

    def __repr__(self):
        return "<CisModule: {}>".format(self.__str__())

