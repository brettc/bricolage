cimport grn

cdef extern from "<src/threshold3.hpp>" namespace "thresh3":
    cdef cppclass cConstructor(grn.cConstructor):
        cConstructor(grn.cWorld_ptr &w, size_t cc)

    cdef cppclass cCisModule(grn.cCisModule):
        cCisModule(cConstructor &c); 
        void mutate(cConstructor &c)
        bint is_active(grn.cChannelState &state)
        grn.int_t binding[3];
        grn.signal_t channels[3];

    cdef cppclass cGene(grn.cGene):
        cGene(grn.sequence_t sequence, grn.signal_t p)
        grn.vector[cCisModule] modules;

    cdef cppclass cNetwork(grn.cNetwork):
        cNetwork(cConstructor &c)
        grn.vector[cGene] genes

    cGene * dynamic_cast_cGene \
        "dynamic_cast<thresh3::cGene *>" (grn.cGene *) except NULL

    cCisModule * dynamic_cast_cCisModule \
        "dynamic_cast<thresh3::cCisModule *>" (grn.cCisModule *) except NULL
