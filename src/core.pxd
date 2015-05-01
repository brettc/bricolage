# Define everything in the external library
from utility cimport *

cdef extern from "<src/core.hpp>" namespace "pubsub2":
    cdef enum:
        MAX_CIS_CHANNELS = 4

    ctypedef unsigned int signal_t

    ctypedef unsigned int sequence_t
    ctypedef int int_t
    ctypedef random_engine_t
    ctypedef uniform_int_distribution[size_t] randint_t

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
    ctypedef shared_ptr[cConstructor] cConstructor_ptr
    ctypedef vector[cNetwork_ptr] cNetworkVector
        
    cdef enum InterventionState:
        INTERVENE_NONE = 0
        INTERVENE_ON = 1
        INTERVENE_OFF = 2

    cdef cppclass cWorld:
        cWorld(size_t seed, size_t cue, size_t reg, size_t out)
        mt19937 rand
        sequence_t next_network_identifier, next_target_identifier
        size_t cue_channels, reg_channels, out_channels
        size_t channel_count
        pair[size_t, size_t] cue_range
        pair[size_t, size_t] out_range
        pair[size_t, size_t] reg_range
        pair[size_t, size_t] sub_range
        pair[size_t, size_t] pub_range
        cChannelStateVector environments
        double get_random_double(double low, double high)
        double get_random_int(int low, int high)
        string get_random_state()
        void set_random_state(string &s)

    ctypedef shared_ptr[cWorld] cWorld_ptr

    cdef cppclass cConstructor:
        cConstructor()
        cNetwork_ptr construct(bint fill)
        size_t site_count(cNetworkVector &networks)
        void mutate_collection(cNetworkVector &networks, 
                               cIndexes &mutated, double site_rate)
        cWorld_ptr world
    
    cdef cppclass cCisModule:
        signal_t get_site(size_t index)
        signal_t set_site(size_t index, signal_t channel);
        size_t site_count()
        InterventionState intervene
        signal_t channels[MAX_CIS_CHANNELS]

    cdef cppclass cGene:
        size_t module_count()
        const cCisModule *get_module(size_t i)
        sequence_t sequence
        signal_t pub
        InterventionState intervene

    cdef cppclass cNetwork:
        cNetwork(cConstructor_ptr &)

        void cycle(cChannelState c)
        void cycle_with_intervention(cChannelState c)
        size_t gene_count()
        void mutate(size_t)
        cGene *get_gene(size_t)
        void calc_attractors()
        void calc_attractors_with_intervention()
        
        void *pyobject
        cConstructor_ptr constructor
        cWorld_ptr world
        int_t identifier, parent_identifier, generation
        cAttractors attractors
        cRatesVector rates
        int_t target
        double fitness
        
        # TODO: Needed for cython bug, never used
        # See https://groups.google.com/forum/#!topic/cython-users/ko7X_fQ0n9Q
        cNetwork() 

    ctypedef pair[char, size_t] Node_t
    ctypedef pair[Node_t, Node_t] Edge_t
    ctypedef std_set[Edge_t] cEdgeList

    cdef cppclass cNetworkAnalysis:
        cNetworkAnalysis(const cNetwork_ptr &n)
        void make_active_edges(cEdgeList e)
        void make_edges(cEdgeList e)

    # # ctypedef shared_ptr[cNetwork] cNetwork_ptr
    # # cNetwork_ptr get_detached_copy(cNetwork_ptr)
    # # ctypedef shared_ptr[const cNetwork] cConstNetwork_ptr
    # # ctypedef vector[cConstNetwork_ptr] cNetworkVector
    # ctypedef vector[cNetwork_ptr] cNetworkVector
    #
    cdef cppclass cTarget:
        cTarget(cWorld_ptr &w, string name, int_t ident)
        double assess(cNetwork &net)
        cWorld *factory
        int_t identifier
        string name
        cRatesVector optimal_rates

    cdef cppclass cSelectionModel:
        cSelectionModel(cWorld_ptr &factory)
        cWorld_ptr factory

        bint select(
            const cNetworkVector &networks, size_t number, cIndexes &selected)

    cdef cppclass cPopulation:
        cPopulation(const cConstructor_ptr &c, size_t n)
        size_t mutate(double site_rate, int_t generation)
        void assess(const cTarget &target)
        bint select(const cSelectionModel &sm, size_t size)
        pair[double, double] worst_and_best()
        void best_indexes(cIndexes &best)
        cConstructor_ptr constructor
        cWorld_ptr world

        sequence_t next_id
        size_t generation
        size_t size

        cIndexes selected, mutated
        cNetworkVector networks

    cdef cppclass cInfoE:
        cInfoE(const cWorld_ptr &world, size_t ncategories)
        cWorld_ptr world
        size_t category_count
        cIndexes categories
        void get_extents(size_t &channels, size_t &categories, size_t &on_off)
        void network_probs(double *data, cNetwork &net)
        void collection_probs(double *data, cNetworkVector &networks)
        void collection_info(double *data, cNetworkVector &networks)

    cdef cppclass cInformation:
        cInformation(const cWorld_ptr &w, size_t networks)
        cWorld_ptr world
        size_t stride_n(size_t n)
        size_t shape_n(size_t n)
        size_t dimensions()
        size_t element_size()
        size_t total_size()
        void *data()

    cdef cppclass cJointProbabilities:
        cJointProbabilities(const cWorld_ptr &w, size_t networks)
        bint calc_information(cInformation &information)
        cWorld_ptr world
        size_t stride_n(size_t n)
        size_t shape_n(size_t n)
        size_t dimensions()
        size_t element_size()
        size_t total_size()
        void *data()
    
    ctypedef vector[double] category_test
    ctypedef vector[category_test] categorizer_list
    cdef cppclass cCausalFlowAnalyzer:
        cCausalFlowAnalyzer(const cWorld_ptr& world, cRates rates);
        cRates rates
        cRates natural_probabilities
        bint analyse_network(cNetwork &net, cJointProbabilities &joint)
        bint analyse_collection(const cNetworkVector &networks, cJointProbabilities &joint)
