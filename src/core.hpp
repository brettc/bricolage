// vim: path=.,/usr/include/c++/4.2.1,/usr/include/c++/4.2.1/tr1,/usr/local/include fdm=syntax
#pragma once

#include <cstdint>
#include <vector>
#include <set>
#include <map>
#include <memory>

// #include <tr1/memory> Another shared_ptr?
// #include <boost/shared_ptr.hpp>
//
#include <boost/dynamic_bitset.hpp>
#include <boost/multi_array.hpp>
#include <tuple>
#include <random>

// Failed attempt to templatize. Let's deal with this later.
// template <typename T, std::size_t NumDims>
// struct multi_array_wrapper
// {
//     multi_array_wrapper::multi_array_wrapper(const extent_gen &ranges)
//     typedef boost::multi_array<T, NumDims> array_type;
//     array_type _array;
//
//     // These provide information for constructing a buffer interface from
//     // python, allowing us to operate on the array via numpy
//     size_t stride_n(size_t n) { return _array.strides()[n]; }
//     size_t shape_n(size_t n) { return _array.shape()[n]; }
//     size_t dimensions() { return _array.num_dimensions(); }
//     size_t element_size() { return sizeof(array_type::element); }
//     size_t total_size() { return element_size() * _array.num_elements(); }
//     void *data() { return _array.data(); }
// };

namespace bricolage
{

// http://stackoverflow.com/questions/1903954/
// is-there-a-standard-sign-function-signum-sgn-in-c-c
template <typename T> inline int c_sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

template <typename T> inline int c_cmp(T a, T b)
{
    return c_sgn(a - b);
}

typedef int_fast8_t signal_t;
typedef signed int int_t;
typedef unsigned int sequence_t;

// TODO: Maybe change dynamic_bitset for something a little lighter (Do we really
// need anything bigger than 64bits?)
typedef boost::dynamic_bitset<size_t> cChannelState;
typedef std::vector<cChannelState> cChannelStateVector;
typedef std::vector<cChannelStateVector> cAttractors;
typedef std::vector<size_t> cIndexes;
typedef std::vector<double> cRates;
typedef std::vector<cRates> cRatesVector;

// TODO: Do these need implementation? I forget my C++.
const size_t reserved_channels = 2;
const size_t off_channel = 0;
const size_t on_channel = 1;

// Set this globally
const size_t MAX_CHANNELS_PER_CIS = 4;

// Allows us to intervene on the state of genes and modules, forcing them on or
// off to ascertain what causal role they have.
enum InterventionState
{
    INTERVENE_NONE = 0,
    INTERVENE_ON = 1,
    INTERVENE_OFF = 2,
};

inline int bitset_cmp(cChannelState &a, cChannelState &b)
{
    if (a < b) return -1;
    if (a == b) return 0;
    return 1;
}

class cCisModule
{
public:
    cCisModule() : intervene(INTERVENE_NONE) {}
    virtual ~cCisModule() {}
    bricolage::signal_t get_site(size_t i) const { return channels[i]; }
    bricolage::signal_t set_site(size_t i, bricolage::signal_t c)
    {
        bricolage::signal_t old = channels[i];
        channels[i] = c;
        return old;
    }
    // Defines how many channels you'll actually use.
    virtual size_t site_count() const = 0;
    InterventionState intervene;
    signal_t channels[MAX_CHANNELS_PER_CIS];
};

struct cGene
{
    cGene(sequence_t sequence, signal_t p);
    virtual ~cGene() {}

    virtual size_t module_count() const=0;
    virtual cCisModule *get_module(size_t i)=0;

    sequence_t sequence;
    signal_t pub;
    InterventionState intervene;
};


class cNetwork;
typedef std::shared_ptr<cNetwork> cNetwork_ptr;
typedef std::vector<cNetwork_ptr> cNetworkVector;

typedef std::shared_ptr<const cNetwork> cConstNetwork_ptr;
// typedef std::vector<cConstNetwork_ptr> cNetworkVector;

typedef std::mt19937 random_engine_t;
typedef std::uniform_int_distribution<size_t> randint_t;
typedef std::uniform_real_distribution<> randreal_t;

class cFactory;

class cWorld //: std::enable_shared_from_this<cWorld>
{
public:
    cWorld(size_t seed, size_t cue, size_t reg, size_t out);
    ~cWorld();

    sequence_t get_next_network_ident() { return next_network_identifier++; }
    sequence_t get_next_target_ident() { return next_target_identifier++; }

    // Counters
    sequence_t next_network_identifier, next_target_identifier;

    // Channel definitions
    size_t cue_channels, reg_channels, out_channels;
    size_t channel_count;

    size_t gene_count;

    std::pair<size_t, size_t> cue_range;
    std::pair<size_t, size_t> reg_range;
    std::pair<size_t, size_t> out_range;

    std::pair<size_t, size_t> sub_range;
    std::pair<size_t, size_t> pub_range;

    // Once all the values are in place we call this
    cChannelStateVector environments;

    // Randomising stuff
    random_engine_t rand;

    std::string get_random_state();
    void set_random_state(const std::string &s);

private:
    void init_channels();
    void init_environments();
public:
    // Extra cruft down here-----
    // Mainly just for testing
    double get_random_double(double low, double high) {
        return std::uniform_real_distribution<double>(low, high)(rand); }
    double get_random_int(int low, int high) {
        return std::uniform_int_distribution<int>(low, high)(rand); }

};
typedef std::shared_ptr<cWorld> cWorld_ptr;

// Base class for overriding world operations
class cFactory : public std::enable_shared_from_this<cFactory>
{
public:
    cFactory(const cWorld_ptr &w);
    virtual ~cFactory() {};

    // This is where the factorys and mutators for networks live
    virtual cNetwork_ptr construct(bool fill)=0;
    virtual size_t site_count(cNetworkVector &networks)=0;

    cNetwork_ptr clone_and_mutate_network(
        cNetwork_ptr &n, size_t mutations, int_t generation);

    cWorld_ptr world;
};

typedef std::shared_ptr<cFactory> cFactory_ptr;

class cNetwork
{
public:
    cNetwork(const cFactory_ptr &c);
    virtual ~cNetwork() {}

    virtual cNetwork_ptr clone() const=0;
    virtual void mutate(size_t nmutations)=0;
    virtual size_t gene_count() const=0;
    virtual cGene *get_gene(size_t i)=0;

    virtual void cycle(cChannelState &c) const=0;
    virtual void cycle_with_intervention(cChannelState &c) const=0;

    void _calc_attractors(bool intervention);
    void calc_attractors() { _calc_attractors(false); }
    void calc_attractors_with_intervention() { _calc_attractors(true); }

    cFactory_ptr factory;
    cWorld_ptr world;
    int_t identifier, parent_identifier;

    // Optional -- the generation that this was created (default: 0)
    int_t generation;

    // Calculated attractor and rates
    cAttractors attractors;
    cRatesVector rates;

    // Record the fitness and the target against which it was calculated.
    // This means we don't need to recalc the fitness if the target has not
    // changed.
    mutable int_t target;
    mutable double fitness;

    // This is where we store the Python object for easy recall.
    mutable void *pyobject;

private:
    // Don't allow copying
    cNetwork(const cNetwork &);
};

struct cTarget
{
    cTarget(const cWorld_ptr &world, const std::string &name, int_t id=-1);
    cWorld_ptr world;
    std::string name;
    int_t identifier;
    cRatesVector optimal_rates;
    cRates weighting;

    // TODO: per env weighting
    // TODO: per output weighting
    double assess(const cNetwork &net) const;
    void assess_networks(cNetworkVector &networks) const;
    void set_weighting(const cRates &w);
};


struct cSelectionModel
{
    cSelectionModel(cWorld_ptr &world);
    cWorld_ptr world;

    // TODO: Make this virtual -- come up with different selection models
    bool select(const cNetworkVector &networks,
                size_t number, cIndexes &selected) const;
};

class cPopulation
{
public:
    cPopulation(const cFactory_ptr &c, size_t n);

    size_t mutate(double site_rate, int_t generation);
    void assess(const cTarget &target) const;
    bool select(const cSelectionModel &sm, size_t size);
    std::pair<double, double> worst_and_best() const;
    void best_indexes(cIndexes &best) const;

    // cConstNetwork_ptr get_network(size_t index) const;
    size_t get_generation() const { return generation; }

// protected:
    cFactory_ptr factory;
    cWorld_ptr world;

    sequence_t next_id;
    size_t generation;

    cIndexes selected;
    cIndexes mutated;

    cNetworkVector networks;
};

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

    void make_edges(cEdgeList &edges);
    void make_active_edges(cEdgeList &edges);
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
    cRates natural_probabilities;

    // We only allocate this many categories. More than this and we're screwed.
    static size_t max_category;
    static size_t get_max_category_size();
    static void set_max_category_size(size_t);

    void _calc_natural(cNetwork &net);
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
    cAverageControlAnalyzer(cWorld_ptr &world);

    // Note you need to delete the return values from these!
    cInformation *analyse_network(cNetwork &net);
    cInformation *analyse_collection(const cNetworkVector &networks);

    void _analyse(cNetwork &net, info_array_type::reference sub);
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

    cOutputCategorizer() : next_category(0) {}
    size_t get_category(const cRates &rates, double prob);
    void clear();
};

struct cOutputControlAnalyzer : public cBaseCausalAnalyzer
{
    cOutputControlAnalyzer(cWorld_ptr &world);
    std::vector<cOutputCategorizer> categorizers;
    cJointProbabilities joint_over_envs;

    // Note you need to delete the return values from these!
    cInformation *analyse_network(cNetwork &net);
    cInformation *analyse_collection(const cNetworkVector &networks);

    void _clear();
    void _analyse(cNetwork &net, info_array_type::reference sub);
};

} // end namespace bricolage
