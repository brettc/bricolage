#ifndef FUNC_H_SK3NECRB
#define FUNC_H_SK3NECRB

#include <vector>
// #include <tr1/memory>
// #include <boost/shared_ptr.hpp>

// #include "pyconfig.h"
// #include "Python.h"
// #include "numpy/npy_common.h"

namespace pubsub2
{

typedef unsigned char uint8;
typedef unsigned int uint32;

struct Gene
{
    Gene();
    uint8 sub1, sub2, pub;
};

typedef std::vector<Gene> GeneVector;

struct Genome
{
    GeneVector genes;

    Genome();
    void init(size_t size);
    // void fill(int i);
};


// typedef boost::shared_ptr<boop> bpp;
//     
// struct tester 
// {
//     tester();
//     int i, j, k;
//     // std::shared_ptr<boop> bp;
//     // std::tr1::shared_ptr<boop> bp;
//     // boost::shared_ptr<boop> bp;
//     bpp bp;
//     int bob();
//
// };
//
// typedef std::vector<tester> tester_vec;

} 

#endif /* end of include guard: FUNC_H_SK3NECRB */

