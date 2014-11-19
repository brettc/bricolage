// vim: path=.,/usr/include/c++/4.2.1,/usr/include/c++/4.2.1/tr1,/usr/local/include
#pragma once

#include <cstdint>
#include <vector>

// #include <tr1/memory> Another shared_ptr?
#include <boost/shared_ptr.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/multi_array.hpp>
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

// TODO: Maybe change dynamic_bitset for something a little light (Do we really
// need anything bigger than 64bits?)
typedef boost::dynamic_bitset<size_t> cChannelState;
typedef std::vector<cChannelState> cChannelStateVector;
typedef std::vector<cChannelStateVector> cAttractors;
typedef std::vector<size_t> cIndexes;
typedef std::vector<double> cRates;
typedef std::vector<cRates> cRatesVector;

inline int bitset_cmp(cChannelState &a, cChannelState &b)
{
    if (a < b) return -1;
    if (a == b) return 0;
    return 1;
}


struct cCisModule
{
    // Default constructor is fine
    operand_t op;
    signal_t sub1, sub2;
    bool silenced;

    // Inline this stuff. It won't change.
    inline bool test(unsigned int a, unsigned int b) const 
    { 
        // Note: C++ standard guarantees integral conversion from bool results
        // in 0 or 1.
        return op & (8 >> ((a << 1) | b)); 
    }

    inline bool is_active(cChannelState const &state) const 
    {
        if (silenced) return false;

        return test(state.test(sub1), state.test(sub2));
    }
};

typedef std::vector<cCisModule> cCisModules;

struct cGene
{
    cGene(sequence_t sequence, signal_t p);

    sequence_t sequence;
    cCisModules modules;
    signal_t pub;
};

typedef std::vector<cGene> cGeneVector;
typedef std::vector<operand_t> cOperands;

struct cNetwork;
typedef boost::shared_ptr<cNetwork> cNetwork_ptr;
typedef std::vector<cNetwork_ptr> cNetworkVector;

typedef std::mt19937 random_engine_t;
typedef std::uniform_int_distribution<size_t> randint_t;
typedef std::uniform_real_distribution<> randreal_t;

struct cFactory
{
    cFactory(size_t seed);

    sequence_t get_next_network_ident() { return next_network_identifier++; }
    sequence_t get_next_target_ident() { return next_target_identifier++; }

    random_engine_t random_engine;
    double get_random_double(double low, double high) {
        return std::uniform_real_distribution<double>(low, high)(random_engine); }
    double get_random_int(int low, int high) {
        return std::uniform_int_distribution<int>(low, high)(random_engine); }


    sequence_t next_network_identifier, next_target_identifier;
    size_t gene_count, cis_count;
    size_t reg_gene_count;
    size_t cue_channels, reg_channels, out_channels;
    std::pair<size_t, size_t> sub_range;
    std::pair<size_t, size_t> pub_range;
    std::pair<size_t, size_t> out_range;

    cOperands operands;

    // Once all the values are in place we call this
    void init_environments();
    cChannelStateVector environments;
    size_t total_channels;

};

typedef boost::shared_ptr<cFactory> cFactory_ptr;

struct cNetwork
{
    cNetwork(cFactory_ptr &f);

    void *pyobject;
    
    // Set this when the network has been finalised
    // bool _locked;

    cFactory_ptr factory;
    int_t identifier, parent_identifier;
    cGeneVector genes;

    void cycle(cChannelState &c);
    void calc_attractors();

    // Calculated attractor and rates
    cAttractors attractors;
    cRatesVector rates;

    // Record the fitness and the target against which it was calculated. 
    // This means we don't need to recalc the fitness if the target has not
    // changed.
    mutable int_t target;
    mutable double fitness;

};

struct cTarget
{
    cTarget(cFactory *factory);
    cFactory *factory;
    int_t identifier;
    cRatesVector optimal_rates;

    // TODO: per env weighting
    // TODO: per output weighting
    double assess(const cNetwork &net);

};

// TODO: Eventually this will be a base class, then have different types of
// generators.
struct cMutationModel
{
    cMutationModel(cFactory *f, double rate_per_gene_);
    cFactory *factory;
    double rate_per_gene;

    // Not sure whether this is worth it, but still...
    randint_t r_gene, r_mod;
    randint_t r_sub, r_oper, r_cis;

    // default constructed [0, 1]
    randreal_t r_uniform_01;

    // This is where the constructors and mutators for networks live
    void construct_cis(cCisModule &m);
    void construct_network(cNetwork &network);

    // Mutates the gene in place
    void mutate_cis(cCisModule &m);
    void mutate_gene(cGene &g);
    void mutate_network(cNetwork_ptr &n, size_t mutations);

    cNetwork_ptr copy_and_mutate_network(cNetwork_ptr &n, size_t mutations);
    void mutate_collection(cNetworkVector &networks, cIndexes &mutated);
};

struct cSelectionModel
{
    cSelectionModel(cFactory_ptr &factory);
    cFactory_ptr factory;

    bool select(
        const cNetworkVector &networks, cTarget &target, 
        size_t number, cIndexes &selected);

    void copy_using_indexes(
        const cNetworkVector &from, cNetworkVector &to, const cIndexes &selected);
};

// struct cSimpleMutationModel : public cMutationModel
// {
//
// };
//

} // end namespace pubsub2
