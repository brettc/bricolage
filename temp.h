#pragma once

#include <vector>
// #include <tr1/memory>
#include <boost/shared_ptr.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/multi_array.hpp>
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
typedef unsigned int npy_uint;
typedef signed long npy_long;
typedef unsigned long npy_ulong;

typedef boost::multi_array<byte_t, 2> garray_type;

typedef boost::dynamic_bitset<> cProducts;
typedef std::vector<cProducts> cProductsVector;

struct cProductStates
{
    cProductStates() { num_products = 0; }

    size_t num_products;
    cProductsVector products;

    void init(size_t np) { num_products = np; }
    void push_back() { products.push_back(cProducts(num_products)); }

    size_t size() const { return products.size(); }
    size_t products_size() const { return num_products; }

    bool get(size_t i, size_t j) { return products[i][j]; }
    void set(size_t i, size_t j, bool b) { products[i][j] = b; }
};

struct cGene
{
    cGene();
    byte_t sub1, sub2, pub;
};

typedef std::vector<cGene> cGeneVector;

struct cNetwork
{
    int_t identifier;
    int_t gene_count;
    cGeneVector genes;
    garray_type &garray;

    // void *object_ptr; // Back ptr to python object

    cNetwork(int_t gene_count);
    ~cNetwork();
    void init(npy_long ident, size_t size);
    void test();
    byte_t *gene_data(); 
};


typedef boost::shared_ptr<cNetwork> cNetwork_ptr;
typedef std::vector<cNetwork_ptr> cNetworkVector;

struct cPopulation
{
    cNetworkVector networks;
};

struct cRandom
{
    cRandom(int seed=0);
    void reseed(int seed) { engine.seed(seed); }
    double get_uniform() { return uniform(engine); }

    std::mt19937 engine;
    std::uniform_real_distribution<double> uniform;
};


} // end namespace pubsub2
