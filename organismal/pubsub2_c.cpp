#include "pubsub2_c.h"
// #include <iostream>
// #include "Python.h"

using namespace pubsub2;

cGene::cGene()
    : sub1(0), sub2(0), pub(0)
{
}

cRandom::cRandom(int seed) 
    : engine(seed), uniform(0, 1) 
{
}


void cNetwork::test()
{
    for (auto g: garray)
        g[0] += 1;
}

cNetwork::cNetwork(int_t gc_)
    : identifier(-1), gene_count(gc_), 
    garray(*(new garray_type(boost::extents[gc_][3])))
{
    // std::cout << "CREATED!" << std::endl;
}

cNetwork::~cNetwork()
{
    // if (garray != 0)
    delete &garray;

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

    // garray = new garray_type(boost::extents[size][3]); 
    //
    garray[0][2] = 55;
    garray[1][1] = 66;
}

byte_t *cNetwork::gene_data()
{
    return garray.data();
}
