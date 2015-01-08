# distutils: include_dirs = NUMPY_PATH
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

cimport core_ext
from logic2 cimport *
from .operand import Operand
from cython.operator import dereference as deref

# TODO: cleanup the numpy definitions here
cimport numpy as np
import numpy
ctypedef np.int_t int_type
ctypedef np.int8_t tiny_type

cdef class Constructor(core_ext.Constructor):
    def __cinit__(self, core_ext.World w, params):
        self._shared = core.cConstructor_ptr(new cConstructor(
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

    # Only used for some mutation models
    property bindings:
        def __get__(self):
            cdef cConstructor *c = dynamic_cast_cConstructor(self._this) 
            return c.bindings

    def dtype(self):
        cdef cConstructor *c = dynamic_cast_cConstructor(self._this) 
        return numpy.dtype([
            ('id', int),
            ('parent', int),
            ('generation', int), # note this is not filled in below
            ('pub', numpy.int8, (c.gene_count)),
            ('op', numpy.int8, (c.gene_count, c.module_count)),
            ('sub', numpy.int8, (c.gene_count, c.module_count, 2)),
        ])

    def to_numpy(self, core_ext.Population p, bint mutations_only=False):
        cdef:
            size_t i, j, k, count, net_i
            cNetwork *net
            cConstructor *c = dynamic_cast_cConstructor(self._this) 

        if mutations_only:
            count = p._this.mutated.size()
        else:
            count = p._this.networks.size()

        output = numpy.zeros(count, dtype=self.dtype())

        if count == 0:
            return output

        # Get the arrays
        cdef: 
            int_type[:] n_id = output['id']
            int_type[:] n_parent = output['parent']
            int_type[:] n_generation = output['generation']
            tiny_type[:,:] n_pub = output['pub']
            tiny_type[:,:,:] n_op = output['op']
            tiny_type[:,:,:,:] n_sub = output['sub']

        for i in range(count):
            if mutations_only:
                # Either get the networks directly...
                net_i = p._this.mutated[i]
            else:
                # Or indirectly via the mutated count
                net_i = i

            net = p._this.networks[net_i].get()
            n_id[i] = net.identifier
            n_parent[i] = net.parent_identifier
            n_generation[i] = net.generation
            for j in range(c.gene_count):
                n_pub[i, j] = net.genes[j].pub
                for k in range(c.module_count):
                    n_op[i, j, k] = net.genes[j].modules[k].op
                    n_sub[i, j, k, 0] = net.genes[j].modules[k].channels[0]
                    n_sub[i, j, k, 1] = net.genes[j].modules[k].channels[1]

        return output

    def from_numpy(self, output, core_ext.Population p):
        cdef: 
            int_type[:] n_id = output['id']
            int_type[:] n_parent = output['parent']
            int_type[:] n_generation = output['generation']
            tiny_type[:,:] n_pub = output['pub']
            tiny_type[:,:,:] n_op = output['op']
            tiny_type[:,:,:,:] n_sub = output['sub']
            size_t i, j, k
            cNetwork *net
            core.cNetwork_ptr ptr
            cConstructor *c = dynamic_cast_cConstructor(self._this) 

        p._this.networks.clear()

        assert n_sub.shape[1] == c.gene_count
        assert n_sub.shape[2] == c.module_count

        for i in range(output.shape[0]):
            ptr = self._this.construct()
            p._this.networks.push_back(ptr)
            net = ptr.get()
            net.identifier = n_id[i]
            net.parent_identifier= n_parent[i]
            net.generation = n_generation[i]
            for j in range(c.gene_count):
                net.genes[j].pub = n_pub[i, j]
                for k in range(c.module_count):
                    net.genes[j].modules[k].op = n_op[i, j, k]
                    net.genes[j].modules[k].channels[0] = n_sub[i, j, k, 0] 
                    net.genes[j].modules[k].channels[1] = n_sub[i, j, k, 1] 
            # TODO: Make a "construct blank"
            net.calc_attractors()

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

