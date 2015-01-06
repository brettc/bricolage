cimport grn

cdef extern from "<src/logic2.hpp>" namespace "logic2":
    ctypedef unsigned int operand_t
    ctypedef grn.vector[operand_t] cOperands
    ctypedef grn.pair[grn.signal_t, grn.signal_t] signal_pair_t
    ctypedef grn.std_map[signal_pair_t, operand_t] binding_map_t

    cdef cppclass cConstructor(grn.cConstructor):
        cConstructor(grn.cWorld_ptr &w, size_t cc, cOperands &ops)
        size_t gene_count, module_count
        cOperands operands
        binding_map_t bindings

    cdef cppclass cCisModule(grn.cCisModule):
        cCisModule(cConstructor &c)
        void mutate(cConstructor &c)
        bint is_active(grn.cChannelState &state)
        bint test(unsigned int a, unsigned int b)
        operand_t op

    cdef cppclass cGene(grn.cGene):
        cGene(grn.sequence_t sequence, grn.signal_t p)
        grn.vector[cCisModule] modules;

    cdef cppclass cNetwork(grn.cNetwork):
        cNetwork(grn.cConstructor_ptr &c)
        grn.vector[cGene] genes

    cConstructor* dynamic_cast_cConstructor \
        "dynamic_cast<logic2::cConstructor*>" (grn.cConstructor *) except NULL

    cNetwork * dynamic_cast_cNetwork \
        "dynamic_cast<logic2::cNetwork *>" (grn.cNetwork *) except NULL

    cGene * dynamic_cast_cGene \
        "dynamic_cast<logic2::cGene *>" (grn.cGene *) except NULL

    cCisModule * dynamic_cast_cCisModule \
        "dynamic_cast<logic2::cCisModule *>" (grn.cCisModule *) except NULL
