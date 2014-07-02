#include "func.h"

using namespace func;

tester::tester()
    : i(10), j(3), k(33), bp(new boop())
{
    this->bp->q = 5;
}


int tester::bob()
{
    return  this->i * this->j * 100;
}
