#pragma once

#include <vector>
// #include <tr1/memory>
#include <boost/shared_ptr.hpp>
#include <boost/multi_array.hpp>

// #include "Python.h"
// #include "numpy/npy_common.h"

namespace pubsub2
{

// I'm simply copying these here from npy_common to avoid including them from
// npy_common
typedef signed char npy_byte;
typedef unsigned char npy_ubyte;
typedef signed int npy_int;
typedef unsigned int npy_uint;
typedef signed long npy_long;
typedef unsigned long npy_ulong;

typedef boost::multi_array<npy_byte, 2> garray_type;

struct cGene
{
    cGene();
    npy_byte sub1, sub2, pub;
};

typedef std::vector<cGene> cGeneVector;

struct cNetwork
{
    npy_int identifier;
    npy_int gene_count;
    cGeneVector genes;
    garray_type *garray;

    // void *object_ptr; // Back ptr to python object

    cNetwork();
    ~cNetwork();
    void init(npy_long ident, size_t size);
    void test();
    npy_byte *gene_data(); 
};

typedef boost::shared_ptr<cNetwork> cNetwork_ptr;
typedef std::vector<cNetwork_ptr> cNetworkVector;

struct cPopulation
{
    cNetworkVector networks;
};


} // end namespace pubsub2
