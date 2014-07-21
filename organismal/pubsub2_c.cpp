#include "pubsub2_c.h"

using namespace pubsub2;

Gene::Gene()
    // : sub1(0), sub2(0), pub(0)
{
}

// void Genes::test()
// {
//     for (auto &i: this->genes)
//         i.pub = 1;
//
// }

Genome::Genome()
{
}

void Genome::init(size_t size)
{
    for (size_t i=0; i < size; ++i)
        this->genes.push_back(Gene());

}

// tester::tester()
//     : i(10), j(3), k(33), bp(new boop())
// {
//     this->bp->q = 5;
// }
//
//
// int tester::bob()
// {
//     return  this->i * this->j * 100;
// }
