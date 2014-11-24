#include "pubsub2_c.h"
// #include <iostream>
// #include "Python.h"

using namespace pubsub2;

cGene::cGene()
    : sub1(0), sub2(0), pub(0)
{
}

void cNetwork::test()
{
    // for (auto &g: genes)
    //     ipub = 1;
}

cNetwork::cNetwork()
    : identifier(-1), garray(0), gene_count(0)
{
    // std::cout << "CREATED!" << std::endl;
}

cNetwork::~cNetwork()
{
    if (garray != 0)
        delete garray;

    // if (this->object_ptr != 0)
    //     Py_XDECREF(this->object_ptr);

    // std::cout << "DYING!!" << std::endl;
}

void cNetwork::init(npy_long ident, size_t size)
{
    gene_count = size;
    identifier = ident;
    for (size_t i=0; i < size; ++i)
        genes.push_back(cGene());

    garray = new garray_type(boost::extents[size][3]); 

    (*garray)[0][2] = 55;
    (*garray)[1][1] = 66;
}

npy_byte *cNetwork::gene_data()
{
    return garray->data();
}
