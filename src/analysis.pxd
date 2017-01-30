# From core stuff
from utility cimport *
from core cimport *

cdef extern from "<src/analysis.hpp>" namespace "bricolage":

    cdef enum NodeType:
        NT_GENE = 0
        NT_MODULE = 1
        NT_CHANNEL = 2

    ctypedef pair[char, size_t] Node_t
    ctypedef pair[Node_t, Node_t] Edge_t
    ctypedef std_set[Edge_t] cEdgeList

    size_t encode_module_id(size_t i, size_t j)

    cdef cppclass cNetworkAnalysis:
        cNetworkAnalysis(const cNetwork_ptr &n)
        void make_active_edges(cEdgeList e)
        void make_edges(cEdgeList e)
        size_t calc_active_bindings()
        cNetwork_ptr original
        cNetwork_ptr modified
        size_t active_bindings, potential_bindings

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
        cRates intervention_probs
        cJointProbabilities *analyse_network(cNetwork &net) except +
        cJointProbabilities *analyse_collection(const cNetworkVector &networks) except +

    cdef cppclass cAverageControlAnalyzer:
        cAverageControlAnalyzer(const cWorld_ptr& world)
        cRates intervention_probs
        cInformation *analyse_network(cNetwork &net) except +
        cInformation *analyse_collection(const cNetworkVector &networks) except +

    cdef cppclass cMutualInfoAnalyzer:
        cMutualInfoAnalyzer(const cWorld_ptr& world, const cIndexes categories);
        cIndexes categories
        cJointProbabilities *analyse_network(cNetwork &net) except +
        cJointProbabilities *analyse_collection(const cNetworkVector &networks) except +

    cdef cppclass cOutputControlAnalyzer:
        cOutputControlAnalyzer(const cWorld_ptr& world, cRatesVector)
        cRates intervention_probs
        cInformation *analyse_network(cNetwork &net) except +
        cInformation *analyse_collection(const cNetworkVector &networks) except +

    cdef cppclass cRelevantControlAnalyzer:
        cRelevantControlAnalyzer(const cWorld_ptr& world, cRatesVector, bint use_natural)
        bint use_natural
        cRates intervention_probs
        cInformation *analyse_network(cNetwork &net) except +
        cInformation *analyse_collection(const cNetworkVector &networks) except +

    cdef cppclass cMIAnalyzer:
        cMIAnalyzer(const cWorld_ptr& world, const cIndexes categories)
        cIndexes categories
        cInformation *analyse_network(cNetwork &net) except +
        cInformation *analyse_collection(const cNetworkVector &networks) except +

    cdef cppclass cWCAnalyzer:
        cWCAnalyzer(const cWorld_ptr& world, const cIndexes ind, 
                    const cRates t1, const cRates t2, double weighting)
        cJointProbabilities *get_joint(cNetwork &net) except +
        cInformation *analyse_network(cNetwork &net) except +
        cInformation *analyse_collection(const cNetworkVector &networks) except +

    cdef cppclass cFastCAnalyzer:
        cFastCAnalyzer(const cWorld_ptr& world, const cIndexes ind, 
                    const cRates t1, const cRates t2)
        cInformation *analyse_network(cNetwork &net) except +
        cInformation *analyse_collection(const cNetworkVector &networks) except +

    cdef cppclass cFastCandBAnalyzer:
        cFastCandBAnalyzer(const cWorld_ptr& world, const cIndexes ind, 
                    const cRates t1, const cRates t2)
        cInformation *analyse_network(cNetwork &net) except +
        cInformation *analyse_collection(const cNetworkVector &networks) except +


cdef extern from "<src/core.hpp>" namespace "bricolage::cBaseCausalAnalyzer":
    # Hack for allowing access to static class functions
    size_t get_max_category_size()
    void set_max_category_size(size_t)

