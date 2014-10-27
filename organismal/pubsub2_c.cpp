#include "pubsub2_c.h"
// #include <iostream>
// #include "Python.h"

using namespace pubsub2;

// cRandom::cRandom(int seed) 
//     : engine(seed), uniform(0, 1) 
// {
// }


// uint_t cBinaryOps::add_op(std::string name, bool offoff, bool offon, 
//                           bool onoff, bool onon)
// {
//     uint_t ret = names.size();
//     names.push_back(name);
//
//     // Store everything in binary
//     uint_t state = 0;
//     if (offoff) state &= 1;
//     if (offon) state &= 1 << 1;
//     if (onoff) state &= 1 << 2;
//     if (onon) state &= 1 << 3;
//
//     results.push_back(state);
//
//     // Return the number of operation 
//     return ret;
// }
//
// bool cBinaryOps::result_from_op(size_t opn, bool a, bool b)
// {
//     // Which operation do we want?
//     uint op = results[opn];
//
//     // Which bit do we want?
//     uint shift = (a & 1) | ((b & 1) << 1);
//
//     // Get that bit out
//     uint r = op & (1 << shift);
//
//     return r != 0;
// };

cGene::cGene(sequence_t seq, signal_t p)
    : sequence(seq)
    , pub(p)
{
}

cNetwork::cNetwork(cFactory_ptr &f)
    : factory(f)
{
    // We initialize
    auto &rng = factory->random_engine;
    std::uniform_int_distribution<> rpub(0, 5);
    
    identifier = factory->get_next_ident();
    for (size_t i=0; i < factory->gene_count; ++i)
        genes.emplace_back(i, rpub(rng));


    // std::cout << "CREATED!" << std::endl;
}

cNetwork::~cNetwork()
{
    // if (this->object_ptr != 0)
    //     Py_XDECREF(this->object_ptr);
    // std::cout << "DYING!!" << std::endl;
}

cFactory::cFactory(size_t seed)
    : random_engine(seed), next_identifier(0), pop_count(1), gene_count(4), cis_count(3)
{
}
