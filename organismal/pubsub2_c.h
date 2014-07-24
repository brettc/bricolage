#ifndef FUNC_H_SK3NECRB
#define FUNC_H_SK3NECRB

#include <vector>
// #include <tr1/memory>
#include <boost/shared_ptr.hpp>

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

struct cGene
{
    cGene();
    npy_byte sub1, sub2, pub;
};

typedef std::vector<cGene> cGeneVector;

struct cNetwork
{
    npy_long identifier;
    cGeneVector genes;
    // void *object_ptr; // Back ptr to python object

    cNetwork();
    ~cNetwork();
    void init(npy_long ident, size_t size);
    // void fill(int i);
};

typedef boost::shared_ptr<cNetwork> cNetwork_ptr;

struct cNetworkSharedPtr 
{
    cNetwork_ptr ptr;
    cNetworkSharedPtr(cNetwork *n=0) : ptr(n) {}
    cNetwork *get() { return this->ptr.get(); }
};

typedef std::vector<cNetworkSharedPtr> cNetworkVector;

struct cPopulation
{
    cNetworkVector networks;
};


} 

#endif /* end of include guard: FUNC_H_SK3NECRB */

