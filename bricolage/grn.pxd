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
        
    cdef cppclass cWorld:
        cWorld(size_t seed, size_t cue, size_t reg, size_t out)
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

    ctypedef shared_ptr[cWorld] cWorld_ptr

    cdef cppclass cConstructor:
        cConstructor()
        cNetwork_ptr construct()
        void mutate_collection(cNetworkVector &networks, 
                               cIndexes &mutated, double site_rate)
        cWorld world
    
    cdef cppclass cCisModule:
        signal_t get_site(size_t index)
        signal_t set_site(size_t index, signal_t channel);
        size_t site_count()
        # InterventionState intervene;

    cdef cppclass cGene:
        size_t module_count()
        const cCisModule *get_module(size_t i)
        sequence_t sequence
        signal_t pub
        # InterventionState intervene;

    cdef cppclass cNetwork:
        cNetwork(cWorld_ptr)
        void *pyobject
        int_t identifier, parent_identifier
        void cycle(cChannelState c)
        void cycle_with_intervention(cChannelState c)
        size_t gene_count()
        void mutate(size_t)
        cGene *get_gene(size_t)
        cAttractors attractors
        cRatesVector rates
        int_t target
        double fitness
        void calc_attractors()
        # bint is_detached()
        
        # TODO: Needed for cython bug, never used
        # See https://groups.google.com/forum/#!topic/cython-users/ko7X_fQ0n9Q
        cNetwork() 

    # ctypedef pair[char, size_t] Node_t
    # ctypedef pair[Node_t, Node_t] Edge_t
    # ctypedef std_set[Edge_t] cEdgeList

    # cdef cppclass cNetworkAnalysis:
    #     cNetworkAnalysis(const cNetwork_ptr &n)
    #     void make_active_edges(cEdgeList e)
    #     void make_edges(cEdgeList e)
    #
    # # ctypedef shared_ptr[cNetwork] cNetwork_ptr
    # # cNetwork_ptr get_detached_copy(cNetwork_ptr)
    # # ctypedef shared_ptr[const cNetwork] cConstNetwork_ptr
    # # ctypedef vector[cConstNetwork_ptr] cNetworkVector
    # ctypedef vector[cNetwork_ptr] cNetworkVector
    #
    cdef cppclass cTarget:
        cTarget(cWorld_ptr &w)
        double assess(cNetwork &net)
        cWorld *factory
        int_t identifier
        cRatesVector optimal_rates

    cdef cppclass cSelectionModel:
        cSelectionModel(cWorld_ptr &factory)
        cWorld_ptr factory

        bint select(
            const cNetworkVector &networks, cTarget &target, 
            size_t number, cIndexes &selected)

        void copy_using_indexes(
            const cNetworkVector &fr, cNetworkVector &to, const cIndexes &selected)

# cdef extern from "<src/logic2.hpp>" namespace "pubsub2":
#     cdef cppclass cConstructorLogic2(cConstructor):
#         cConstructorLogic2(cWorld &f, size_t gc, size_t cc, cOperands &ops)
#
#     cdef cppclass cCisModuleLogic2(cCisModule):
#         operand_t op
#         signal_t *channels
#
#     cCisModuleLogic2 * dynamic_cast_cCisModuleLogic2 \
#         "dynamic_cast<pubsub2::cCisModuleLogic2 *>" (cCisModule *) except NULL



