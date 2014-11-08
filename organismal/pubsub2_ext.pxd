from _cpp cimport *

cdef extern from "pubsub2_c.h" namespace "pubsub2":
    ctypedef unsigned int signal_t
    ctypedef unsigned int operand_t
    ctypedef unsigned int sequence_t
    ctypedef random_engine_t
    ctypedef uniform_int_distribution[size_t] randint_t

    ctypedef vector[operand_t] cOperands

    ctypedef dynamic_bitset[size_t] cChannelState
    ctypedef vector[cChannelState] cChannelStateVector
    ctypedef vector[cChannelStateVector] cAttractors

    cdef int bitset_cmp(cChannelState &, cChannelState &)

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

        void construct_random(cNetwork &network)

        void init_environments()
        cChannelStateVector environments

    ctypedef shared_ptr[cFactory] cFactory_ptr

    ctypedef vector[cCisModule] cCisModules

    cdef cppclass cCisModule:
        bint test(unsigned int a, unsigned int b)
        bint active(dynamic_bitset[size_t] s)
        signal_t op, sub1, sub2

    ctypedef vector[cCisModule] cCisModules

    cdef cppclass cGene:
        sequence_t sequence;
        cCisModules modules;
        signal_t pub

    cdef cppclass cNetwork:
        cNetwork(cFactory_ptr)
        void *pyobject
        vector[cGene] genes
        sequence_t identifier
        size_t gene_count
        void cycle(cChannelState c)
        cAttractors attractors

    ctypedef shared_ptr[cNetwork] cNetwork_ptr
    ctypedef vector[cNetwork_ptr] cNetworkVector

    cdef cppclass cGeneMutator:
        cGeneMutator(cFactory *f, double rate)
        cFactory *factory;
        randint_t r_sub, r_oper, r_cis
        void mutate_gene(cGene &g)
        void mutate_network(cNetwork_ptr &n, size_t n)
        cNetwork_ptr copy_and_mutate_network(cNetwork_ptr &n, size_t mutations)
        void mutate_collection(cNetworkVector &networks)
