// vim: path=.,/usr/include/c++/4.2.1,/usr/include/c++/4.2.1/tr1,/usr/local/include fdm=syntax
#pragma once

#include <cstdint>
#include <vector>
#include <set>
#include <memory>

// #include <tr1/memory> Another shared_ptr?
// #include <boost/shared_ptr.hpp>
//
#include <boost/dynamic_bitset.hpp>
#include <boost/multi_array.hpp>
#include <tuple>
#include <random>

namespace pubsub2
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
    pubsub2::signal_t get_site(size_t i) const { return channels[i]; }
    pubsub2::signal_t set_site(size_t i, pubsub2::signal_t c) 
    { 
        pubsub2::signal_t old = channels[i];
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

class cConstructor;

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
class cConstructor : public std::enable_shared_from_this<cConstructor>
{
public:
    cConstructor(const cWorld_ptr &w);
    virtual ~cConstructor() {};

    // This is where the constructors and mutators for networks live
    virtual cNetwork_ptr construct(bool fill)=0;
    virtual size_t site_count(cNetworkVector &networks)=0;

    cNetwork_ptr clone_and_mutate_network(
        cNetwork_ptr &n, size_t mutations, int_t generation);

    cWorld_ptr world;
};

typedef std::shared_ptr<cConstructor> cConstructor_ptr;

class cNetwork
{
public:
    cNetwork(const cConstructor_ptr &c);
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
    
    cConstructor_ptr constructor;
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

    // TODO: per env weighting
    // TODO: per output weighting
    double assess(const cNetwork &net) const;
};

typedef boost::multi_array<double, 3> freqs_t;
typedef boost::multi_array<double, 4> networks_probs_t;
typedef boost::multi_array_ref<double, 4> networks_probs_ref_t;
struct cEnvironmentI
{
    cEnvironmentI(const cWorld_ptr &world, size_t nc);
    ~cEnvironmentI();
    cWorld_ptr world;
    size_t category_count;
    freqs_t *frequencies;
    cIndexes categories;
    void copy_frequencies(double *view) const;
    void calculate(const cNetwork &net);

    void get_extents(size_t &channels, size_t &categories, size_t &on_off);
    void calc_collection(double *data, const cNetworkVector &networks);
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
    cPopulation(const cConstructor_ptr &c, size_t n);

    size_t mutate(double site_rate, int_t generation);
    void assess(const cTarget &target) const;
    bool select(const cSelectionModel &sm, size_t size);
    std::pair<double, double> worst_and_best() const;
    void best_indexes(cIndexes &best) const;

    // cConstNetwork_ptr get_network(size_t index) const;
    size_t get_generation() const { return generation; }

// protected:
    cConstructor_ptr constructor;
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


} // end namespace pubsub2
