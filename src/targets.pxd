# Define everything in the external library
from utility cimport *
from core cimport *

cdef extern from "<src/targets.hpp>" namespace "bricolage":
    
    cdef enum ScoringMethod:
        SCORE_LINEAR = 0
        SCORE_EXPONENTIAL = 1
        SCORE_EXPONENTIAL_VEC = 2
        SCORE_POWER_HAMMING = 3

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
                     int_t perturb_count, double perturb_prop, bint env_only)
        size_t perturb_count
        double perturb_prop
        bint env_only

    cdef cppclass cTransNoisyTarget(cBaseTarget):
        cTransNoisyTarget(cWorld_ptr &w, string name, int_t ident, 
                     ScoringMethod meth, double strength)

    cdef cppclass cMultiTarget(cBaseTarget):
        cMultiTarget(cWorld_ptr &w, string name, int_t ident, 
                     ScoringMethod meth, double strength)
        cAttractor pulses

    cdef cppclass cSelectionModel:
        cSelectionModel(cWorld_ptr &factory)
        cWorld_ptr factory

        bint select(
            const cNetworkVector &networks, size_t number, 
            cIndexes &selected, const cBaseTarget &target)


