# distutils: include_dirs = NUMPY_PATH
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

from .operand import Operand
from logic2_cpp cimport *
cimport core_ext
from cython.operator import dereference as deref, preincrement as preinc

cdef class Constructor(core_ext.Constructor):
    def __cinit__(self, core_ext.World w, params):
        self._shared = grn.cConstructor_ptr(new cConstructor(
            w._shared,
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

    property bindings:
        def __get__(self):
            cdef cConstructor *c = dynamic_cast_cConstructor(self._this) 
            return c.bindings
    

cdef class CisModule(core_ext.CisModule):
    property op:
        def __get__(self):
            cdef cCisModule *cm = dynamic_cast_cCisModule(self._this) 
            return Operand(cm.op)
        def __set__(self, size_t op):
            cdef cCisModule *cm = dynamic_cast_cCisModule(self._this) 
            operand = Operand(op)
            cm.op = op

    def test(self, bint a, bint b):
        cdef cCisModule *cm = dynamic_cast_cCisModule(self._this) 
        return cm.test(a, b)

    def is_active(self, core_ext.ChannelState state):
        cdef cCisModule *cm = dynamic_cast_cCisModule(self._this) 
        return cm.is_active(state._this)

    def mutate(self):
        cdef cCisModule *cm = dynamic_cast_cCisModule(self._this) 
        cm.mutate(deref(self.gene.network._this.constructor.get()))

    def __str__(self):
        w = self.gene.network.constructor.world
        return "{}({}, {})".format(
            self.op.name,
            w.name_for_channel(self._this.channels[0]),
            w.name_for_channel(self._this.channels[1]),
        )

    def __repr__(self):
        return "<CisModule: {}>".format(self.__str__())

