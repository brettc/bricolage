cimport core

cdef extern from "<src/threshold3.hpp>" namespace "thresh3":
    cdef cppclass cConstructor(core.cConstructor):
        cConstructor(core.cWorld_ptr &w, size_t cc)
        size_t gene_count, module_count

    cdef cppclass cCisModule(core.cCisModule):
        cCisModule(cConstructor &c); 
        void mutate(cConstructor &c)
        bint is_active(core.cChannelState &state)
        core.int_t binding[3];
        core.signal_t channels[3];

    cdef cppclass cGene(core.cGene):
        cGene(core.sequence_t sequence, core.signal_t p)
        core.vector[cCisModule] modules;

    cdef cppclass cNetwork(core.cNetwork):
        cNetwork(cConstructor &c)
        core.vector[cGene] genes

    cConstructor* dynamic_cast_cConstructor \
        "dynamic_cast<thresh3::cConstructor*>" (core.cConstructor *) except NULL

    cNetwork * dynamic_cast_cNetwork \
        "dynamic_cast<thresh3::cNetwork *>" (core.cNetwork *) except NULL

    cGene * dynamic_cast_cGene \
        "dynamic_cast<thresh3::cGene *>" (core.cGene *) except NULL

    cCisModule * dynamic_cast_cCisModule \
        "dynamic_cast<thresh3::cCisModule *>" (core.cCisModule *) except NULL
