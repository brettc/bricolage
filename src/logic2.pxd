cimport core

cdef extern from "<src/logic2.hpp>" namespace "logic2":
    ctypedef unsigned int operand_t
    ctypedef core.vector[operand_t] cOperands
    ctypedef core.pair[core.signal_t, core.signal_t] signal_pair_t
    # ctypedef core.std_map[signal_pair_t, operand_t] binding_map_t

    cdef cppclass cFactory(core.cFactory):
        cFactory(core.cWorld_ptr &w, size_t cc, cOperands &ops)
        size_t gene_count, module_count
        cOperands operands
        # binding_map_t bindings

    cdef cppclass cCisModule(core.cCisModule):
        cCisModule(cFactory &c)
        void mutate(cFactory &c)
        bint is_active(core.cChannels &state)
        bint test(unsigned int a, unsigned int b)
        operand_t op

    cdef cppclass cGene(core.cGene):
        cGene(core.sequence_t sequence, core.signal_t p)
        core.vector[cCisModule] modules;

    cdef cppclass cNetwork(core.cNetwork):
        cNetwork(core.cFactory_ptr &c)
        core.vector[cGene] genes

    cFactory* dynamic_cast_cFactory \
        "dynamic_cast<logic2::cFactory*>" (core.cFactory *) except NULL

    cNetwork * dynamic_cast_cNetwork \
        "dynamic_cast<logic2::cNetwork *>" (core.cNetwork *) except NULL

    cGene * dynamic_cast_cGene \
        "dynamic_cast<logic2::cGene *>" (core.cGene *) except NULL

    cCisModule * dynamic_cast_cCisModule \
        "dynamic_cast<logic2::cCisModule *>" (core.cCisModule *) except NULL
