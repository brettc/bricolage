#pragma once
#include "core.hpp"

namespace bricolage
{


enum NodeType { NT_GENE=0, NT_MODULE, NT_CHANNEL };
typedef std::pair<NodeType, size_t> Node_t;
typedef std::pair<Node_t, Node_t> Edge_t;
typedef std::set<Edge_t> cEdgeList;

// TODO: Should really make this into a union, but not sure if it makes life
// harder...
inline size_t make_module_node_id(size_t g, size_t m)
{
    // Combine gene and module id
    return g << 8 | m;
}

struct cNetworkAnalysis
{
    cNetworkAnalysis(cNetwork_ptr &n);
    cNetwork_ptr original;
    cNetwork_ptr modified;
    size_t potential_bindings, active_bindings;

    void make_edges(cEdgeList &edges);
    void make_active_edges(cEdgeList &edges);
    size_t calc_active_bindings();
};

struct cJointProbabilities;

typedef boost::multi_array<double, 3> info_array_type;
struct cInformation
{
    cInformation(const cJointProbabilities &jp);
    cInformation(const cWorld_ptr &w, size_t network_size, size_t per_channel_size);

    cWorld_ptr world;
    info_array_type _array;

    // These provide information for constructing a buffer interface from
    // python, allowing us to operate on the array via numpy
    size_t stride_n(size_t n) { return _array.strides()[n]; }
    size_t shape_n(size_t n) { return _array.shape()[n]; }
    size_t dimensions() { return _array.num_dimensions(); }
    size_t element_size() { return sizeof(info_array_type::element); }
    size_t total_size() { return element_size() * _array.num_elements(); }
    void *data() { return _array.data(); }
};

typedef boost::multi_array<double, 5> joint_array_type;
struct cJointProbabilities
{
    cJointProbabilities(const cWorld_ptr &w, size_t network_size,
                        size_t per_network, size_t per_channel);
    cWorld_ptr world;
    joint_array_type _array;

    void calc_information(cInformation &information) const;

    // These provide information for constructing a buffer interface from
    // python, allowing us to operate on the array via numpy
    size_t stride_n(size_t n) { return _array.strides()[n]; }
    size_t shape_n(size_t n) { return _array.shape()[n]; }
    size_t dimensions() { return _array.num_dimensions(); }
    size_t element_size() { return sizeof(joint_array_type::element); }
    size_t total_size() { return element_size() * _array.num_elements(); }
    void *data() { return _array.data(); }
};

struct cBaseCausalAnalyzer
{
    cBaseCausalAnalyzer(cWorld_ptr &world);
    cWorld_ptr world;
    cRates intervention_probs;

    // We only allocate this many categories. More than this and we're screwed.
    static size_t max_category;
    static size_t get_max_category_size();
    static void set_max_category_size(size_t);

    void _calc_natural(cNetwork &net);
};

struct cRateCategorizer
{
    // We only allocate this many categories. More than this and we're screwed.
    size_t next_category;
    std::map<double, int> rate_categories;
    std::vector<double> category_probabilities;

    cRateCategorizer() : next_category(0) {}
    size_t get_category(double rate, double prob);
    void clear();
};


struct cCausalFlowAnalyzer : public cBaseCausalAnalyzer
{
    cCausalFlowAnalyzer(cWorld_ptr &world);

    // Note you need to delete the return values from these!
    cJointProbabilities *analyse_network(cNetwork &net);
    cJointProbabilities *analyse_collection(const cNetworkVector &networks);

    void _analyse(cNetwork &net, joint_array_type::reference sub);
};

struct cAverageControlAnalyzer : public cBaseCausalAnalyzer
{
    boost::multi_array<cRateCategorizer, 2> categorizers;
    cJointProbabilities joint_over_envs;

    cAverageControlAnalyzer(cWorld_ptr &world);

    // Note you need to delete the return values from these!
    cInformation *analyse_network(cNetwork &net);
    cInformation *analyse_collection(const cNetworkVector &networks);

    void _analyse(cNetwork &net, info_array_type::reference sub);
    void _clear();
};

struct cMutualInfoAnalyzer
{
    cMutualInfoAnalyzer(cWorld_ptr &world, const cIndexes &cats);
    cWorld_ptr world;
    cIndexes categories;
    size_t max_category;

    // Note you need to delete the return values from these!
    cJointProbabilities *analyse_network(cNetwork &net);
    cJointProbabilities *analyse_collection(const cNetworkVector &networks);

    void _analyse(cNetwork &net, joint_array_type::reference sub);
};

//-----------------------------------------------------------------------------
// Keep a map of the unique OUTPUTS (all rates) and assign them to persistent
// categories within a network. We'll use these categories to calculate the
// information. We also keep the total probabilities of all of these outputs,
// as this then allows us to calculate the per-regulatory ENTROPY of the
// outputs.
struct cOutputCategorizer
{
    size_t next_category;
    std::map<std::vector<double>, size_t> rate_categories;
    std::vector<double> category_probabilities;
    const cRatesVector &target_rates;
    std::vector<double> targets_hit_in_env;

    cOutputCategorizer(const cRatesVector &tr, size_t env_size);
    size_t get_category(const cRates &rates, double prob, size_t env);
    void clear();
};

struct cOutputControlAnalyzer : public cBaseCausalAnalyzer
{
    cOutputControlAnalyzer(cWorld_ptr &world, const cRatesVector &tr);
    std::vector<cOutputCategorizer> categorizers;
    cJointProbabilities joint_over_envs;
    cRatesVector target_rates;

    // Note you need to delete the return values from these!
    cInformation *analyse_network(cNetwork &net);
    cInformation *analyse_collection(const cNetworkVector &networks);

    void _clear();
    void _analyse(cNetwork &net, info_array_type::reference sub);
};

struct cRelevantControlAnalyzer : public cBaseCausalAnalyzer
{
    cRelevantControlAnalyzer(cWorld_ptr &world, const cRatesVector &tr, bool use_natural_);
    cRatesVector target_rates;
    boost::multi_array<int, 2> categories;
    boost::multi_array<double, 2> info;
    bool use_natural;

    // Note you need to delete the return values from these!
    int categorize(const cRates rates);
    cInformation *analyse_network(cNetwork &net);
    cInformation *analyse_collection(const cNetworkVector &networks);


    void _clear();
    void _analyse(cNetwork &net, info_array_type::reference sub);
};

//---------------------------------------------------------------------------
// New information measures that handle multiple categories
//
struct cMIAnalyzer
{
    cMIAnalyzer(cWorld_ptr &world, const cIndexes &cats);
    cWorld_ptr world;
    cIndexes categories;
    size_t num_categories;

    // Note you need to delete the return values from these!
    cInformation *analyse_network(cNetwork &net);
    cInformation *analyse_collection(const cNetworkVector &networks);

    void _analyse(cNetwork &net, joint_array_type::reference sub);
};


// Weighted control
struct cWCAnalyzer
{
    // Construct using two targets
    cWCAnalyzer(cWorld_ptr &world, 
                const cIndexes &indexes,
                const cRates &t1, 
                const cRates &t2,
                double weighting);

    // Keep this around for general sizes
    const cWorld_ptr world_ptr;

    // Convenience
    const cWorld &world; 
    
    // Which particular rates do we want?
    cIndexes indexes;

    // What are the values of the targets?
    const cRates target1, target2;

    // Weighting value 
    double weighting;

    // What values have we found?
    cRatesVector rates_found;

    size_t processing_index;

    // What probabilities do these rates have 
    typedef std::array<double, 2> off_on_type;
    typedef std::vector<off_on_type> off_on_vector_type;
    off_on_vector_type empty_probs;

    struct RateDetail
    {
        RateDetail(size_t enc, size_t idx, size_t sz) 
            : encountered(enc)
            , used_for(idx)
            , probs(sz, {{0.0, 0.0}})
        {}

        size_t encountered; // the order it was encountered
        size_t used_for; // indicates if it is in use for this index;
        off_on_vector_type probs; // OFF/ON per reg channel
        std::array<double, 2> similarity; // target1, target2
    };

    typedef std::map<cRates, RateDetail> rate_detail_map_type;
    rate_detail_map_type rate_detail_map;

    // Getting actual values back.
    // Note you need to delete the return values from these!
    cInformation *analyse_network(cNetwork &net);
    cInformation *analyse_collection(const cNetworkVector &networks);

    cJointProbabilities *get_joint(cNetwork &net);

    // The internal functions
    // void _calc_weightings();
    void _analyse(cNetwork &net, info_array_type::reference sub);
    inline double _similarity(const cRates &a, const cRates &b);
    void _wiggle(cNetwork &net);
    inline void _add_probability(const cRates &rates, size_t gindex, size_t is_on, double pr);
    inline RateDetail& _match(const cRates &rates);
    // void _clear();

};


} // end namespace bricolage
// vim: path=.,/usr/include/c++/4.2.1,/usr/include/c++/4.2.1/tr1,/usr/local/include fdm=syntax
