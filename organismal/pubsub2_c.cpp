#include "pubsub2_c.h"
#include <iostream>
// #include "Python.h"

using namespace pubsub2;

cGene::cGene()
    : sub1(0), sub2(0), pub(0)
{
}

// void Genes::test()
// {
//     for (auto &i: this->genes)
//         i.pub = 1;
//
// }

cNetwork::cNetwork()
    : identifier(-1)
{
    std::cout << "CREATED!" << std::endl;
}

cNetwork::~cNetwork()
{
    // if (this->object_ptr != 0)
    //     Py_XDECREF(this->object_ptr);

    std::cout << "DYING!!" << std::endl;
}

void cNetwork::init(int32 ident, size_t size)
{
    this->identifier = ident;
    for (size_t i=0; i < size; ++i)
        this->genes.push_back(cGene());

}
