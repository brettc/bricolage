# Define everything in the external library
from utility cimport *

cdef extern from "<src/core.hpp>" namespace "pubsub2":
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

    cdef cppclass cConstructor
    cdef cppclass cNetwork
    ctypedef shared_ptr[cNetwork] cNetwork_ptr
    ctypedef vector[cNetwork_ptr] cNetworkVector
        
    cdef cppclass cFactory:
        cFactory(size_t seed, size_t cue, size_t reg, size_t out)
        mt19937 rand
        size_t cue_channels, reg_channels, out_channels
        size_t channel_count
        pair[size_t, size_t] cue_range
        pair[size_t, size_t] out_range
        pair[size_t, size_t] reg_range
        pair[size_t, size_t] sub_range
        pair[size_t, size_t] pub_range
        cChannelStateVector environments
        cConstructor *constructor
        double get_random_double(double low, double high)
        double get_random_int(int low, int high)

    ctypedef shared_ptr[cFactory] cFactory_ptr

    cdef cppclass cConstructor:
        cConstructor()
        cConstructor(cFactory &f, size_t gene_count_, size_t cis_count_)
        void construct_network(cNetwork &network);
        void mutate_network(cNetwork_ptr &n, size_t n)
        cNetwork_ptr copy_and_mutate_network(
            cNetwork_ptr &n, size_t mutations)
        void mutate_collection(
            cNetworkVector &networks, cIndexes &mutated, double site_rate)
        size_t gene_count, cis_count;
        randint_t r_gene, r_cis, r_sub;

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

    ctypedef pair[char, size_t] Node_t
    ctypedef pair[Node_t, Node_t] Edge_t
    ctypedef std_set[Edge_t] cEdgeList

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

    cdef cppclass cSelectionModel:
        cSelectionModel(cFactory_ptr &factory)
        cFactory_ptr factory

        bint select(
            const cNetworkVector &networks, cTarget &target, 
            size_t number, cIndexes &selected)

        void copy_using_indexes(
            const cNetworkVector &fr, cNetworkVector &to, const cIndexes &selected)

cdef extern from "<src/logic2.hpp>" namespace "pubsub2":
    cdef cppclass cConstructorLogic2(cConstructor):
        cConstructorLogic2(cFactory &f, size_t gc, size_t cc, cOperands &ops)

    cdef cppclass cCisModuleLogic2(cCisModule):
        operand_t op
        signal_t *channels

    cCisModuleLogic2 * dynamic_cast_cCisModuleLogic2 \
        "dynamic_cast<pubsub2::cCisModuleLogic2 *>" (cCisModule *) except NULL

cdef extern from "<src/threshold3.hpp>" namespace "pubsub2":
    cdef cppclass cConstructorThreshold3(cConstructor):
        cConstructorThreshold3(cFactory &f, size_t gc, size_t cc)

    cdef cppclass cCisModuleThreshold3(cCisModule):
        signal_t channels[3]
        int_t binding[3]

    cCisModuleThreshold3 * dynamic_cast_cCisModuleThreshold3 \
        "dynamic_cast<pubsub2::cCisModuleThreshold3*>" (cCisModule *) except NULL

