// vim: path=.,/usr/include/c++/4.2.1,/usr/include/c++/4.2.1/tr1,/usr/local/include fdm=syntax
#pragma once

#include <cstdint>
#include <vector>
#include <set>

// #include <tr1/memory> Another shared_ptr?
#include <boost/shared_ptr.hpp>
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
typedef int_fast16_t operand_t;
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
    virtual signal_t get_site(size_t i) const=0;
    virtual signal_t set_site(size_t i, signal_t c)=0;
    virtual size_t site_count() const = 0;

    InterventionState intervene;
};

// typedef std::vector<cCisModule *> cCisModules;

struct cGene
{
    cGene(sequence_t sequence, signal_t p);
    virtual ~cGene();

    virtual size_t module_count() const=0;
    virtual const cCisModule *get_module(size_t i)=0;

    sequence_t sequence;
    InterventionState intervene;
    signal_t pub;
};


class cNetwork;
typedef boost::shared_ptr<cNetwork> cNetwork_ptr;
typedef std::vector<cNetwork_ptr> cNetworkVector;

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

    std::pair<size_t, size_t> cue_range;
    std::pair<size_t, size_t> reg_range;
    std::pair<size_t, size_t> out_range;

    std::pair<size_t, size_t> sub_range;
    std::pair<size_t, size_t> pub_range;

    // Once all the values are in place we call this
    cChannelStateVector environments;

    // Randomising stuff
    random_engine_t rand;

    // The specialised constructor for networks
    cConstructor *constructor;
protected:
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
typedef boost::shared_ptr<cWorld> cWorld_ptr;

// Base class for overriding world operations
class cConstructor
{
public:
    cConstructor(const cWorld_ptr &w);
    virtual ~cConstructor() {};

    // This is where the constructors and mutators for networks live
    virtual cNetwork *construct()=0;
    // virtual void mutate_network(cNetwork_ptr *n, size_t mutations)=0;
    // cNetwork_ptr copy_and_mutate_network(cNetwork_ptr *n, size_t mutations);

    // Mutates the gene in place
    // virtual void mutate_cis(cCisModule *m) = 0;
    // void mutate_gene(cGene *g);
    // void mutate_network(cNetwork_ptr &n, size_t mutations);
    // void mutate_collection(cNetworkVector &networks, cIndexes &mutated, double site_rate);

    cWorld_ptr world;
    // Basic parameters
    // size_t gene_count, cis_count;
    // randint_t r_gene, r_cis, r_sub;
    // randreal_t r_uniform_01;
};

class cNetwork
{
public:
    cNetwork(const cWorld_ptr &f, bool no_ident=false);
    virtual ~cNetwork() {}

    // size_t gene_count() const { return genes.size(); }
    virtual cNetwork *clone() const=0;
    virtual void cycle(cChannelState &c) const=0;

    void calc_attractors();
    // void clone_genes(cGeneVector &gv) const;
    // cNetwork_ptr get_detached_copy() const;
    // bool is_detached() const { return identifier < 0; }
    
    cWorld_ptr world;
    int_t identifier, parent_identifier;

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

// cNetwork_ptr get_detached_copy(cNetwork_ptr original);
// inline cNetwork_ptr get_detached_copy(cNetwork_ptr original)
// {
//     cNetwork *copy = new cNetwork(original->world, true);
//     copy->parent_identifier = original->identifier;
//     copy->genes = original->genes;
//
//     return cNetwork_ptr(copy);
// }

// enum NodeType { NT_GENE=0, NT_MODULE, NT_CHANNEL };
// typedef std::pair<NodeType, size_t> Node_t;
// typedef std::pair<Node_t, Node_t> Edge_t;
// typedef std::set<Edge_t> cEdgeList;
//
// // TODO: Should really make this into a union, but not sure if it makes life
// // harder...
// inline size_t make_module_node_id(size_t g, size_t m)
// {
//     // Combine gene and module id
//     return g << 8 | m;
// }
//
// struct cNetworkAnalysis
// {
//     cNetworkAnalysis(const cNetwork_ptr &n);
//     cNetwork_ptr original;
//     cNetwork modified;
//
//     void make_edges(cEdgeList &edges);
//     void make_active_edges(cEdgeList &edges);
// };
//
// struct cTarget
// {
//     cTarget(cWorld *world);
//     cWorld *world;
//     int_t identifier;
//     cRatesVector optimal_rates;
//
//     // TODO: per env weighting
//     // TODO: per output weighting
//     double assess(const cNetwork &net);
// };
//
// struct cSelectionModel
// {
//     cSelectionModel(cWorld_ptr &world);
//     cWorld_ptr world;
//
//     bool select(
//         const cNetworkVector &networks, cTarget &target, 
//         size_t number, cIndexes &selected);
//
//     void copy_using_indexes(
//         const cNetworkVector &from, cNetworkVector &to, const cIndexes &selected);
// };
//
} // end namespace pubsub2
