# From core stuff
from utility cimport *

cdef extern from "<src/defines.hpp>" namespace "bricolage":
    ctypedef unsigned int signal_t

    ctypedef unsigned int sequence_t
    ctypedef int int_t
    ctypedef random_engine_t
    ctypedef uniform_int_distribution[size_t] randint_t
    ctypedef unsigned long bits_t
    ctypedef unsigned char index_t

    cdef int c_sgn(int)
    cdef int c_cmp(int, int)


cdef extern from "<src/channels.hpp>" namespace "bricolage":

    cdef cppclass cChannels:
        cChannels()
        bint test(index_t i, index_t sz) except +
        void set(index_t i, index_t sz) except +
        void clear(index_t i, index_t sz) except +
        void flip(index_t i, index_t sz) except +
        string to_string(index_t size) except +
        void unchecked_union(cChannels &other)
        void unchecked_intersection(cChannels &other)
        index_t max() 

        bits_t bits

    cdef int bitset_cmp(cChannels &, cChannels &)


cdef extern from "<src/core.hpp>" namespace "bricolage":
    cdef enum:
        MAX_CIS_CHANNELS = 4


    ctypedef vector[cChannels] cAttractor
    ctypedef vector[bits_t] cAttractorBits

    # Cast between these two
    cAttractorBits* to_cAttractorBits\
        "reinterpret_cast<bricolage::cAttractorBits *>" (cAttractor *) except NULL

    ctypedef vector[cAttractor] cAttractorSet
    ctypedef vector[cAttractorBits] cAttractorSetBits

    cAttractorSetBits* to_cAttractorSetBits\
        "reinterpret_cast<bricolage::cAttractorSetBits *>" (cAttractorSet *) except NULL


    ctypedef vector[size_t] cIndexes
    ctypedef vector[double] cRates

    ctypedef std_map[cChannels, cRates] cChannelsRatesMap
    ctypedef std_map[bits_t, cRates] cChannelsRatesMapBits

    cChannelsRatesMapBits* to_cChannelRatesMapBits\
        "reinterpret_cast<bricolage::cChannelsRatesMapBits *>" (cChannelsRatesMap *) except NULL
    ctypedef vector[cRates] cRatesVector

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
        cWorld(size_t seed, size_t cue, size_t reg, size_t out, size_t reg_genes)
        mt19937 rand
        sequence_t next_network_identifier, next_target_identifier
        size_t reg_gene_count
        size_t cue_channels, reg_channels, out_channels
        size_t channel_count
        pair[size_t, size_t] cue_range
        pair[size_t, size_t] out_range
        pair[size_t, size_t] reg_range
        pair[size_t, size_t] sub_range
        pair[size_t, size_t] pub_range
        cAttractor environments
        InputType input_type
        size_t pulse_for
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
            cNetwork_ptr &n, size_t n_sub, size_t n_pub, size_t dups, int_t generation)
        void set_draw_from_subs(cIndexes &d)
        void set_draw_from_regs(cIndexes &d)
        cIndexes draw_from_subs, draw_from_regs
        size_t gene_count, module_count
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

        void cycle(cChannels c)
        void cycle_with_intervention(cChannels c)
        size_t gene_count()
        void mutate(size_t, size_t)
        void duplicate(size_t)
        cNetwork_ptr clone()
        cGene *get_gene(size_t)
        void calc_attractors()
        void calc_attractors_with_intervention()
        void stabilise(cAttractor &, bint intervene, cAttractor &, cAttractor &,
                       cRates &)
        void get_rates(cChannels &initial, cRates &rates, bint use_cache)
        void clear_rate_cache()
        double attractor_robustness()

        void *pyobject
        cFactory_ptr factory
        cWorld_ptr world
        int_t identifier, parent_identifier, generation
        cAttractorSet attractors
        cRatesVector rates, transients
        cChannelsRatesMap cached_mappings
        int_t target
        double fitness
        
        # TODO: Needed for cython bug, never used
        # See https://groups.google.com/forum/#!topic/cython-users/ko7X_fQ0n9Q
        cNetwork() 


    # Forward declarations for selection.pxd
    cdef cppclass cBaseTarget
    cdef cppclass cSelectionModel

    cdef cppclass cPopulation:
        cPopulation(const cFactory_ptr &c, size_t n)
        size_t mutate(double cis_rate, double trans_rate, double dup_rate, 
                      int_t generation)
        void assess(const cBaseTarget &target)
        bint select(const cSelectionModel &sm, size_t size)
        pair[double, double] worst_and_best()
        void best_indexes(cIndexes &best)
        cFactory_ptr factory
        cWorld_ptr world

        sequence_t next_id
        size_t generation
        size_t size

        cIndexes selected, mutated
        cNetworkVector networks
        cRates fitnesses

