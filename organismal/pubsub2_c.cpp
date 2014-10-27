// vim: path=.,/usr/include/c++/4.2.1,/usr/include/c++/4.2.1/tr1
#include "pubsub2_c.h"
// #include <iostream>
// #include "Python.h"

using namespace pubsub2;

cGene::cGene(sequence_t seq, signal_t p)
    : sequence(seq)
    , pub(p)
{
}

cNetwork::cNetwork(cFactory_ptr &f)
    : pyobject(0)
    , factory(f)
{
    // We initialize
    auto &rng = factory->random_engine;
    std::uniform_int_distribution<> rpub(0, 5);
    
    identifier = factory->get_next_ident();
    for (size_t i=0; i < factory->gene_count; ++i)
        genes.push_back(cGene(i, rpub(rng)));
        // genes.emplace_back(i, rpub(rng));
    // std::cout << "CREATED!" << std::endl;
    //
}

cNetwork::~cNetwork()
{
    // if (this->object_ptr != 0)
    //     Py_XDECREF(this->object_ptr);
    // std::cout << "DYING!!" << std::endl;
}

cFactory::cFactory(size_t seed)
    : pyobject(0)
    , random_engine(seed)
    , next_identifier(0)
    , pop_count(1)
    , gene_count(4)
    , cis_count(3)
{
}
