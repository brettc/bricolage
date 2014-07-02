#ifndef FUNC_H_SK3NECRB
#define FUNC_H_SK3NECRB

#include <vector>
// #include <tr1/memory>
#include <boost/shared_ptr.hpp>
// #include <memory>

namespace func
{
struct boop
{
    int q;
};

typedef boost::shared_ptr<boop> bpp;
    
struct tester 
{
    tester();
    int i, j, k;
    // std::shared_ptr<boop> bp;
    // std::tr1::shared_ptr<boop> bp;
    // boost::shared_ptr<boop> bp;
    bpp bp;
    int bob();

};

typedef std::vector<tester> tester_vec;

} 

#endif /* end of include guard: FUNC_H_SK3NECRB */

