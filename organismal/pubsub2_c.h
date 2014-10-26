#pragma once

#include <vector>
// #include <tr1/memory>
#include <boost/shared_ptr.hpp>
#include <boost/dynamic_bitset.hpp>
#include <random>

// #include "Python.h"
// #include "numpy/npy_common.h"
namespace pubsub2
{

// I'm simply copying these here from npy_common to avoid including them from
// npy_common
typedef signed char byte_t;
typedef unsigned char npy_ubyte;
typedef signed int int_t;
typedef unsigned int uint_t;

typedef unsigned int npy_uint;
typedef signed long npy_long;
typedef unsigned long npy_ulong;


struct cRandom
{
    cRandom(int seed=0);
    void reseed(int seed) { engine.seed(seed); }
    double get_uniform() { return uniform(engine); }

    std::mt19937 engine;
    std::uniform_real_distribution<double> uniform;
    std::uniform_int_distribution<double> uniform;
};

typedef boost::dynamic_bitset<size_t> cProducts;
typedef std::vector<cProducts> cProductsVector;

struct cProductStates
{
    cProductStates() { num_products = 0; }

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
typedef std::vector<uint_t> cResults;

struct cBinaryOps
{
    cBinaryOps() { num_ops = 0; }

    uint_t add_op(std::string name, bool offoff, bool offon, bool onoff, bool onon);
    bool result_from_op(size_t opn, bool a, bool b);
    size_t num_ops;
    cNames names;
    cResults results;
};

struct cGene
{
    cGene();
    byte_t op, sub1, sub2, pub;

    // Inline this stuff. It won't change.
    inline bool test(byte_t a, byte_t b) const 
    { 
        return op & (1 << ((a << 1) | b)); 
    }

    inline bool active(cProducts const &products) const 
    {
        byte_t a = products[sub1] ? 1: 0;
        byte_t b = products[sub2] ? 1: 0;
        return test(a, b);
    }

};

typedef std::vector<cGene> cGeneVector;

struct cNetwork
{
    int_t identifier;
    int_t gene_count;
    cGeneVector genes;
    // void *object_ptr; // Back ptr to python object

    cNetwork(int_t gene_count);
    ~cNetwork();
    void init(npy_long ident, size_t size);
};


typedef boost::shared_ptr<cNetwork> cNetwork_ptr;
typedef std::vector<cNetwork_ptr> cNetworkVector;

struct cPopulation
{
    cNetworkVector networks;
};



} // end namespace pubsub2
