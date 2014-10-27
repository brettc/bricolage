// vim: path=.,/usr/include/c++/4.2.1,/usr/include/c++/4.2.1/tr1
#pragma once

#include <cstdint>
#include <vector>
// #include <tr1/memory> Another shared_ptr?
#include <boost/shared_ptr.hpp>
#include <boost/dynamic_bitset.hpp>
#include <random>

namespace pubsub2
{

// I'm simply copying these here from npy_common to avoid including them from
// npy_common
typedef int_fast8_t signal_t;
typedef int_fast16_t operand_t;
typedef signed int int_t;
typedef unsigned int sequence_t;

typedef boost::dynamic_bitset<size_t> cProducts;
typedef std::vector<cProducts> cProductsVector;

struct cProductsSequence
{
    cProductsSequence() { num_products = 0; }

    size_t num_products;
    cProductsVector products;

    void init(size_t np) { num_products = np; }
    void push_back(cProducts &p) { products.push_back(p); }

    size_t size() const { return products.size(); }
    size_t products_size() const { return num_products; }

    bool get(size_t i, size_t j) { return products[i][j]; }
    void set(size_t i, size_t j, bool b) { products[i][j] = b; }
};

typedef std::vector<std::string> cNames;

struct cCisModule
{
    cCisModule();
    operand_t op;
    signal_t sub1, sub2;

    // Inline this stuff. It won't change.
    inline bool test(operand_t a, operand_t b) const 
    { 
        return op & (1 << ((a << 1) | b)); 
    }

    inline bool active(cProducts const &products) const 
    {
        operand_t a = products[sub1] ? 1: 0;
        operand_t b = products[sub2] ? 1: 0;
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

struct cFactory
{
    cFactory(size_t seed);
    sequence_t get_next_ident() { return next_identifier++; }

    std::mt19937 random_engine;

    sequence_t next_identifier;
    size_t pop_count, gene_count, cis_count;

};

typedef boost::shared_ptr<cFactory> cFactory_ptr;

struct cNetwork
{
    cNetwork(cFactory_ptr &f);
    ~cNetwork();

    cFactory_ptr factory;
    sequence_t identifier;
    size_t gene_count;
    cGeneVector genes;
    // void *object_ptr; // Back ptr to python object

};

typedef boost::shared_ptr<cNetwork> cNetwork_ptr;
typedef std::vector<cNetwork_ptr> cNetworkVector;

struct cPopulation
{
    cNetworkVector networks;
};




} // end namespace pubsub2
