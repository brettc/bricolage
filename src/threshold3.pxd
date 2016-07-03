cimport core

cdef extern from "<src/threshold3.hpp>" namespace "thresh3":
    cdef enum MutateType:
        JUMP=0
        PROGRESSIVE=1

    cdef cppclass cFactory(core.cFactory):
        cFactory(core.cWorld_ptr &w, size_t cc, MutateType)

    cdef cppclass cCisModule(core.cCisModule):
        cCisModule(cFactory &c); 
        void mutate(cFactory &c)
        bint is_active(core.cChannels &state)
        core.int_t binding[3];
        core.signal_t channels[3];

    cdef cppclass cGene(core.cGene):
        cGene(core.sequence_t sequence, core.signal_t p)
        core.vector[cCisModule] modules;

    cdef cppclass cNetwork(core.cNetwork):
        cNetwork(cFactory &c)
        core.vector[cGene] genes

    cFactory* dynamic_cast_cFactory \
        "dynamic_cast<thresh3::cFactory*>" (core.cFactory *) except NULL

    cNetwork * dynamic_cast_cNetwork \
        "dynamic_cast<thresh3::cNetwork *>" (core.cNetwork *) except NULL

    cGene * dynamic_cast_cGene \
        "dynamic_cast<thresh3::cGene *>" (core.cGene *) except NULL

    cCisModule * dynamic_cast_cCisModule \
        "dynamic_cast<thresh3::cCisModule *>" (core.cCisModule *) except NULL
