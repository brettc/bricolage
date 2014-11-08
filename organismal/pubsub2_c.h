// vim: path=.,/usr/include/c++/4.2.1,/usr/include/c++/4.2.1/tr1,/usr/local/include
#pragma once

#include <cstdint>
#include <vector>

// #include <tr1/memory> Another shared_ptr?
#include <boost/shared_ptr.hpp>
#include <boost/dynamic_bitset.hpp>
#include <random>

namespace pubsub2
{

typedef int_fast8_t signal_t;
typedef int_fast16_t operand_t;
typedef signed int int_t;
typedef unsigned int sequence_t;

typedef boost::dynamic_bitset<size_t> cChannelState;
typedef std::vector<cChannelState> cChannelStateVector;
typedef std::vector<cChannelStateVector> cAttractors;

inline int bitset_cmp(cChannelState &a, cChannelState &b)
{
    if (a < b) return -1;
    if (a == b) return 0;
    return 1;
}


struct cCisModule
{
    cCisModule(operand_t op_, signal_t sub1_, signal_t sub2_);
    operand_t op;
    signal_t sub1, sub2;

    // Inline this stuff. It won't change.
    inline bool test(unsigned int a, unsigned int b) const 
    { 
        return op & (8 >> ((a << 1) | b)); 
    }

    inline bool active(cChannelState const &state) const 
    {
        unsigned int a = state[sub1] ? 1: 0;
        unsigned int b = state[sub2] ? 1: 0;
        return test(a, b);
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

struct cFactory
{
    cFactory(size_t seed);
    sequence_t get_next_ident() { return next_identifier++; }
    void construct_random(cNetwork &network);
    void construct_population(cNetworkVector &nv);

    random_engine_t random_engine;
    sequence_t next_identifier;
    size_t pop_count, gene_count, cis_count;
    size_t reg_gene_count;
    size_t cue_channels, reg_channels, out_channels;
    std::pair<size_t, size_t> sub_range;
    std::pair<size_t, size_t> pub_range;

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

    cFactory_ptr factory;
    sequence_t identifier;
    cGeneVector genes;

    void cycle(cChannelState &c);
    void calc_attractors();

    // Calculated attractor
    cAttractors attractors;

};

// TODO: Pbly should be a base class, then have different types of generators.
struct cGeneMutator
{
    cGeneMutator(cFactory *f, double rate);
    cFactory *factory;
    double rate_per_gene;
    randint_t r_gene;
    randint_t r_sub, r_oper, r_cis;

    // Mutates the gene in place
    void mutate_gene(cGene &g);
    void mutate_network(cNetwork_ptr &n, size_t mutations);
    cNetwork_ptr copy_and_mutate_network(cNetwork_ptr &n, size_t mutations);
    void mutate_collection(cNetworkVector &networks);
};


} // end namespace pubsub2
