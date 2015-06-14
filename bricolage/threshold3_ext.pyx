# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

cimport core_ext
from threshold3 cimport *
from cython.operator import dereference as deref, preincrement as preinc

from .logic_tools import boolean_func_from_coop_binding

cimport numpy as np
import numpy
ctypedef np.int_t int_type
ctypedef np.int8_t tiny_type

def _construct_factory(core_ext.World w):
    return Constructor(w)

cdef class Constructor(core_ext.Constructor):
    def __cinit__(self, core_ext.World w):
        # Hack for old stored classes
        mt = w._params.__dict__.get('mutate_type', 0)

        self._shared = core.cConstructor_ptr(new cConstructor(
            w._shared,
            w._params.cis_count,
            mt,
        ))
        self._this = self._shared.get()

        # Specialise the python classes
        self.module_class = CisModule

    def __init__(self, core_ext.World w):
        addz = w._params.__dict__.get('add_zeros', 0)
        if addz > 0:
            self.draw_from_subs += [0] * addz

    property draw_from_subs:
        def __get__(self):
            cdef cConstructor *c = dynamic_cast_cConstructor(self._this) 
            return c.draw_from_subs
        def __set__(self, core.cIndexes dsubs):
            cdef cConstructor *c = dynamic_cast_cConstructor(self._this) 
            c.set_draw_from_subs(dsubs)

    def __reduce__(self):
        return _construct_factory, (self.world, )

    def dtype(self):
        cdef cConstructor *c = dynamic_cast_cConstructor(self._this) 
        return numpy.dtype([
            ('id', int),
            ('parent', int),
            ('generation', int), # note this is not filled in below
            ('sub', numpy.int8, (c.gene_count, c.module_count, 3)),
            ('binding', numpy.int8, (c.gene_count, c.module_count, 3)),
        ])

    cdef _to_numpy(self, core.cNetworkVector *networks, core.cIndexes *indexes):
        cdef:
            size_t i, j, k, l, count, net_i
            cNetwork *net
            cConstructor *c = dynamic_cast_cConstructor(self._this) 

        if indexes == NULL:
            count = networks.size()
        else:
            count = indexes.size()

        output = numpy.zeros(count, dtype=self.dtype())

        if count == 0:
            return output

        # Get the arrays
        cdef: 
            int_type[:] n_id = output['id']
            int_type[:] n_parent = output['parent']
            int_type[:] n_generation = output['generation']
            tiny_type[:,:,:,:] n_sub = output['sub']
            tiny_type[:,:,:,:] n_bnd = output['binding']

        for i in range(count):
            if indexes == NULL:
                # Either get the networks directly...
                net_i = i
            else:
                # Or indirectly via the indexes
                net_i = deref(indexes)[i]

            net = deref(networks)[net_i].get()
            n_id[i] = net.identifier
            n_parent[i] = net.parent_identifier
            n_generation[i] = net.generation
            for j in range(c.gene_count):
                for k in range(c.module_count):
                    for l in range(3):
                        n_sub[i, j, k, l] = net.genes[j].modules[k].channels[l]
                        n_bnd[i, j, k, l] = net.genes[j].modules[k].binding[l]

        return output

    def to_numpy(self, core_ext.CollectionBase c):
        return self._to_numpy(c._collection, NULL)

    def pop_to_numpy(self, core_ext.Population p, bint mutations_only=False):
        cdef:
            core.cIndexes *indexes = NULL
            core.cNetworkVector *networks = &p._this.networks
        if mutations_only:
            indexes = &p._this.mutated
        return self._to_numpy(networks, indexes)

    cdef _from_numpy(self, output, core.cNetworkVector *networks):
        cdef: 
            int_type[:] n_id = output['id']
            int_type[:] n_parent = output['parent']
            int_type[:] n_generation = output['generation']
            tiny_type[:,:,:,:] n_sub = output['sub']
            tiny_type[:,:,:,:] n_bnd = output['binding']
            size_t i, j, k, l
            cNetwork *net
            core.cNetwork_ptr ptr
            cConstructor *c = dynamic_cast_cConstructor(self._this) 

        networks.clear()

        assert n_sub.shape[1] == c.gene_count
        assert n_sub.shape[2] == c.module_count

        for i in range(output.shape[0]):
            ptr = self._this.construct(False)
            networks.push_back(ptr)
            net = ptr.get()
            net.identifier = n_id[i]
            net.parent_identifier= n_parent[i]
            net.generation = n_generation[i]
            for j in range(c.gene_count):
                for k in range(c.module_count):
                    for l in range(3):
                        net.genes[j].modules[k].channels[l] = n_sub[i, j, k, l] 
                        net.genes[j].modules[k].binding[l] = n_bnd[i, j, k, l] 
            net.calc_attractors()

    def from_numpy(self, output, core_ext.CollectionBase c):
        self._from_numpy(output, c._collection)
        
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

    def same_as(self, other):
        return \
            self.channels == other.channels and\
            self.bindings == other.bindings
