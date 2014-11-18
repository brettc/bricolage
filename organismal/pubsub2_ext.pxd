from _cpp cimport *

cdef extern from "pubsub2_c.h" namespace "pubsub2":
    ctypedef unsigned int signal_t
    ctypedef unsigned int operand_t
    ctypedef unsigned int sequence_t
    ctypedef int int_t
    ctypedef random_engine_t
    ctypedef uniform_int_distribution[size_t] randint_t

    ctypedef vector[operand_t] cOperands

    ctypedef dynamic_bitset[size_t] cChannelState
    ctypedef vector[cChannelState] cChannelStateVector
    ctypedef vector[cChannelStateVector] cAttractors
    ctypedef vector[size_t] cIndexes
    ctypedef vector[double] cRates
    ctypedef vector[cRates] cRatesVector

    cdef int bitset_cmp(cChannelState &, cChannelState &)
    cdef int c_sgn(int)
    cdef int c_cmp(int, int)

    cdef cppclass cFactory:
        cFactory(size_t seed)
        mt19937 random_engine
        sequence_t next_identifier
        size_t pop_count, gene_count, cis_count
        size_t cue_channels, reg_channels, out_channels
        size_t total_channels
        cOperands operands
        pair[size_t, size_t] sub_range
        pair[size_t, size_t] pub_range
        pair[size_t, size_t] out_range

        void init_environments()
        cChannelStateVector environments

    ctypedef shared_ptr[cFactory] cFactory_ptr

    ctypedef vector[cCisModule] cCisModules

    cdef cppclass cCisModule:
        bint test(unsigned int a, unsigned int b)
        bint is_active(dynamic_bitset[size_t] s)
        operand_t op
        signal_t sub1, sub2
        bint silenced

    ctypedef vector[cCisModule] cCisModules

    cdef cppclass cGene:
        sequence_t sequence;
        cCisModules modules;
        signal_t pub

    cdef cppclass cNetwork:
        cNetwork(cFactory_ptr)
        void *pyobject
        vector[cGene] genes
        int_t identifier, parent_identifier
        size_t gene_count
        void cycle(cChannelState c)
        cAttractors attractors
        cRatesVector rates

    ctypedef shared_ptr[cNetwork] cNetwork_ptr
    ctypedef vector[cNetwork_ptr] cNetworkVector

    cdef cppclass cMutationModel:
        cMutationModel(cFactory *f, double rate)
        cFactory *factory

        void construct_cis(cCisModule &m);
        void construct_network(cNetwork &network);

        void mutate_cis(cCisModule &m)
        void mutate_gene(cGene &g)
        void mutate_network(cNetwork_ptr &n, size_t n)

        cNetwork_ptr copy_and_mutate_network(cNetwork_ptr &n, size_t mutations)
        void mutate_collection(cNetworkVector &networks, cIndexes &mutated)
