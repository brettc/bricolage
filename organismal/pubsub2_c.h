#ifndef FUNC_H_SK3NECRB
#define FUNC_H_SK3NECRB

#include <vector>
// #include <tr1/memory>
#include <boost/shared_ptr.hpp>

// #include "pyconfig.h"
// #include "Python.h"
// #include "numpy/npy_common.h"

namespace pubsub2
{

typedef unsigned char uint8;
typedef unsigned int uint32;
typedef long int32;

struct cGene
{
    cGene();
    uint8 sub1, sub2, pub;
};

typedef std::vector<cGene> cGeneVector;

struct cNetwork
{
    int32 identifier;
    cGeneVector genes;
    // void *object_ptr; // Back ptr to python object

    cNetwork();
    ~cNetwork();
    void init(int32 ident, size_t size);
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

