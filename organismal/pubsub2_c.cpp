#include "pubsub2_c.h"
// #include <iostream>
// #include "Python.h"

using namespace pubsub2;

cRandom::cRandom(int seed) 
    : engine(seed), uniform(0, 1) 
{
}


uint_t cBinaryOps::add_op(std::string name, bool offoff, bool offon, 
                          bool onoff, bool onon)
{
    uint_t ret = names.size();
    names.push_back(name);

    // Store everything in binary
    uint_t state = 0;
    if (offoff) state &= 1;
    if (offon) state &= 1 << 1;
    if (onoff) state &= 1 << 2;
    if (onon) state &= 1 << 3;

    results.push_back(state);

    // Return the number of operation 
    return ret;
}

bool cBinaryOps::result_from_op(size_t opn, bool a, bool b)
{
    // Which operation do we want?
    uint op = results[opn];

    // Which bit do we want?
    uint shift = (a & 1) | ((b & 1) << 1);

    // Get that bit out
    uint r = op & (1 << shift);

    return r != 0;
};

cGene::cGene()
    : sub1(0), sub2(0), pub(0)
{
}

cNetwork::cNetwork(int_t gc_)
    : identifier(-1), gene_count(gc_)
{
    // std::cout << "CREATED!" << std::endl;
}

cNetwork::~cNetwork()
{
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
}
