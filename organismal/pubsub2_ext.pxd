from _cpp cimport *

cdef extern from "<pubsub2_c.h>" namespace "pubsub2":
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
        double get_random_double(double low, double high)
        double get_random_int(int low, int high)

        sequence_t next_identifier
        size_t pop_count, gene_count, cis_count
        size_t cue_channels, reg_channels, out_channels
        size_t total_channels
        cOperands operands
        pair[size_t, size_t] sub_range
        pair[size_t, size_t] pub_range
        pair[size_t, size_t] out_range

        cChannelStateVector environments
        void init_environments()

    ctypedef shared_ptr[cFactory] cFactory_ptr

    ctypedef vector[cCisModule] cCisModules

    cdef cppclass cCisModule:
        # bint test(unsigned int a, unsigned int b)
        bint is_active(dynamic_bitset[size_t] s)
        signal_t get_site_channel(size_t index)
        signal_t set_site_channel(size_t index, signal_t channel);
        size_t site_count()
        
    ctypedef vector[cCisModule *] cCisModules

    cdef cppclass cGene:
        sequence_t sequence;
        cCisModules modules;
        signal_t pub

    cdef cppclass cNetwork
    ctypedef shared_ptr[cNetwork] cNetwork_ptr

    cdef cppclass cNetwork:
        cNetwork(cFactory_ptr)
        void *pyobject
        vector[cGene *] genes
        int_t identifier, parent_identifier
        size_t gene_count
        void cycle(cChannelState c)
        cAttractors attractors
        cRatesVector rates
        int_t target
        double fitness
        void calc_attractors()
        bint is_detached()

    # cdef cppclass cSiteIndex:
    #     int_t gene()
    #     int_t cis()
    #     int_t site()
    #
    # ctypedef vector[cSiteIndex] cSiteLocations
    ctypedef pair[char, size_t] Node_t
    ctypedef pair[Node_t, Node_t] Edge_t
    ctypedef set[Edge_t] cEdgeList

    cdef cppclass cNetworkAnalysis:
        cNetworkAnalysis(const cNetwork_ptr &n)
        void make_active_edges(cEdgeList e)
        void make_edges(cEdgeList e)

    # ctypedef shared_ptr[cNetwork] cNetwork_ptr
    # cNetwork_ptr get_detached_copy(cNetwork_ptr)
    # ctypedef shared_ptr[const cNetwork] cConstNetwork_ptr
    # ctypedef vector[cConstNetwork_ptr] cNetworkVector
    ctypedef vector[cNetwork_ptr] cNetworkVector

    cdef cppclass cTarget:
        cTarget(cFactory *factory)
        double assess(cNetwork &net)
        cFactory *factory
        int_t identifier
        cRatesVector optimal_rates

    cdef cppclass cGeneFactory:
        cGeneFactory()
        cGeneFactory(cFactory *f, double rate)
        cFactory *factory

        void construct_network(cNetwork &network);
        void mutate_network(cNetwork_ptr &n, size_t n)

        cNetwork_ptr copy_and_mutate_network(cNetwork_ptr &n, size_t mutations)
        void mutate_collection(cNetworkVector &networks, cIndexes &mutated)

    cdef cppclass cSelectionModel:
        cSelectionModel(cFactory_ptr &factory)
        cFactory_ptr factory

        bint select(
            const cNetworkVector &networks, cTarget &target, 
            size_t number, cIndexes &selected)

        void copy_using_indexes(
            const cNetworkVector &fr, cNetworkVector &to, const cIndexes &selected)

    cdef cppclass cGeneFactoryLogic2(cGeneFactory):
        cGeneFactoryLogic2(cFactory *f, double rate)

    cdef cppclass cCisModuleLogic2(cCisModule):
        operand_t op
        signal_t *channels

    cCisModuleLogic2 * dynamic_cast_cCisModuleLogic2 \
        "dynamic_cast<pubsub2::cCisModuleLogic2 *>" (cCisModule *) except NULL

