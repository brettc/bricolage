# Define everything in the external library
from utility cimport *

cdef extern from "<src/core.hpp>" namespace "bricolage":
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

    cdef cppclass cFactory
    cdef cppclass cNetwork
    ctypedef shared_ptr[cNetwork] cNetwork_ptr
    ctypedef shared_ptr[cFactory] cFactory_ptr
    ctypedef vector[cNetwork_ptr] cNetworkVector
        
    cdef enum InterventionState:
        INTERVENE_NONE = 0
        INTERVENE_ON = 1
        INTERVENE_OFF = 2

    cdef enum InputType:
        INPUT_CONSTANT = 0
        INPUT_PULSE = 1

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
        InputType input_type
        double get_random_double(double low, double high)
        double get_random_int(int low, int high)
        string get_random_state()
        void set_random_state(string &s)

    ctypedef shared_ptr[cWorld] cWorld_ptr

    cdef cppclass cFactory:
        cFactory()
        cNetwork_ptr construct(bint fill)
        size_t site_count(cNetworkVector &networks)
        cNetwork_ptr clone_and_mutate_network(
            cNetwork_ptr &n, size_t mutations, int_t generation)
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
        cNetwork(cFactory_ptr &)

        void cycle(cChannelState c)
        void cycle_with_intervention(cChannelState c)
        size_t gene_count()
        void mutate(size_t)
        cNetwork_ptr clone()
        cGene *get_gene(size_t)
        void calc_attractors()
        void calc_attractors_with_intervention()
        void calc_perturbation()
        
        void *pyobject
        cFactory_ptr factory
        cWorld_ptr world
        int_t identifier, parent_identifier, generation
        cAttractors attractors, pert_attractors
        cRatesVector rates, pert_rates
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
        size_t calc_active_bindings()
        cNetwork_ptr original
        cNetwork_ptr modified
        size_t active_bindings, potential_bindings

    # # ctypedef shared_ptr[cNetwork] cNetwork_ptr
    # # cNetwork_ptr get_detached_copy(cNetwork_ptr)
    # # ctypedef shared_ptr[const cNetwork] cConstNetwork_ptr
    # # ctypedef vector[cConstNetwork_ptr] cNetworkVector
    # ctypedef vector[cNetwork_ptr] cNetworkVector
    
    cdef enum ScoringMethod:
        SCORE_LINEAR = 0
        SCORE_EXPONENTIAL = 1

    cdef cppclass cBaseTarget:
        cBaseTarget(cWorld_ptr &w, string name, int_t ident,
                       ScoringMethod meth, double strength)

        void assess_networks(cNetworkVector &networks,
                             cRates &fitnesses)
        double assess(cNetwork &net)
        void set_weighting(const cRates &w);
        cWorld *factory
        int_t identifier
        string name
        cRates weighting
        cRatesVector optimal_rates
        ScoringMethod scoring_method
        double strength

    cdef cppclass cDefaultTarget(cBaseTarget):
        cDefaultTarget(cWorld_ptr &w, string name, int_t ident,
                       ScoringMethod meth, double strength)

    cdef cppclass cNoisyTarget(cBaseTarget):
        cNoisyTarget(cWorld_ptr &w, string name, int_t ident, 
                     ScoringMethod meth, double strength, 
                     int_t perturb_count)
        size_t perturb_count

    cdef cppclass cSelectionModel:
        cSelectionModel(cWorld_ptr &factory)
        cWorld_ptr factory

        bint select(
            const cNetworkVector &networks, size_t number, 
            cIndexes &selected, const cBaseTarget &target)

    cdef cppclass cPopulation:
        cPopulation(const cFactory_ptr &c, size_t n)
        size_t mutate(double site_rate, int_t generation)
        void assess(const cBaseTarget &target)
        bint select(const cSelectionModel &sm, size_t size,
                    const cBaseTarget &target)
        pair[double, double] worst_and_best()
        void best_indexes(cIndexes &best)
        cFactory_ptr factory
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

    cdef cppclass cJointProbabilities

    cdef cppclass cInformation:
        cInformation(const cJointProbabilities &joint)
        cWorld_ptr world
        size_t stride_n(size_t n)
        size_t shape_n(size_t n)
        size_t dimensions()
        size_t element_size()
        size_t total_size()
        void *data()

    cdef cppclass cJointProbabilities:
        cJointProbabilities(const cWorld_ptr &w, size_t network_size, 
                        size_t per_network, size_t per_channel)
        bint calc_information(cInformation &information)
        cWorld_ptr world
        size_t stride_n(size_t n)
        size_t shape_n(size_t n)
        size_t dimensions()
        size_t element_size()
        size_t total_size()
        void *data()

    cdef cppclass cCausalFlowAnalyzer:
        cCausalFlowAnalyzer(const cWorld_ptr& world)
        cRates natural_probabilities
        cJointProbabilities *analyse_network(cNetwork &net) except +
        cJointProbabilities *analyse_collection(const cNetworkVector &networks) except +

    cdef cppclass cAverageControlAnalyzer:
        cAverageControlAnalyzer(const cWorld_ptr& world)
        cRates natural_probabilities
        cInformation *analyse_network(cNetwork &net) except +
        cInformation *analyse_collection(const cNetworkVector &networks) except +

    cdef cppclass cMutualInfoAnalyzer:
        cMutualInfoAnalyzer(const cWorld_ptr& world, const cIndexes categories);
        cIndexes categories
        cJointProbabilities *analyse_network(cNetwork &net) except +
        cJointProbabilities *analyse_collection(const cNetworkVector &networks) except +

    cdef cppclass cOutputControlAnalyzer:
        cOutputControlAnalyzer(const cWorld_ptr& world, cRatesVector)
        cRates natural_probabilities
        cInformation *analyse_network(cNetwork &net) except +
        cInformation *analyse_collection(const cNetworkVector &networks) except +

    cdef cppclass cRelevantControlAnalyzer:
        cRelevantControlAnalyzer(const cWorld_ptr& world, cRatesVector)
        cRates natural_probabilities
        cInformation *analyse_network(cNetwork &net) except +
        cInformation *analyse_collection(const cNetworkVector &networks) except +

cdef extern from "<src/core.hpp>" namespace "bricolage::cBaseCausalAnalyzer":
    # Hack for allowing access to static class functions
    size_t get_max_category_size()
    void set_max_category_size(size_t)

