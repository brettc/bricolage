// vim: path=.,/usr/include/c++/4.2.1,/usr/include/c++/4.2.1/tr1
#include "pubsub2_c.h"
// #include <iostream>
// #include "Python.h"

using namespace pubsub2;

cCisModule::cCisModule(operand_t op_, signal_t sub1_, signal_t sub2_)
    : op(op_)
    , sub1(sub1_)
    , sub2(sub2_)
{
}

cGene::cGene(sequence_t seq, signal_t p)
    : sequence(seq)
    , pub(p)
{
}

cNetwork::cNetwork(cFactory_ptr &f)
    : pyobject(0)
    , factory(f)
{
    identifier = factory->get_next_ident();

    // We initialize
    auto &rng = factory->random_engine;
    auto &operands = factory->operands;
    auto &sub = factory->sub_range;
    auto pub_base = factory->pub_range.first;

    std::uniform_int_distribution<size_t> r_sub(sub.first, sub.second-1);
    std::uniform_int_distribution<size_t> r_oper(0, operands.size()-1);
    
    for (size_t i=0; i < factory->gene_count; ++i)
    {
        genes.push_back(cGene(i, pub_base + i));
        // genes.emplace_back(i, rpub(rng));
        // Get a reference to it
        cGene &g = genes.back();
        for (size_t j=0; j < factory->cis_count; ++j)
        {
            g.modules.push_back(
                cCisModule(operands[r_oper(rng)], 
                           r_sub(rng), 
                           r_sub(rng)
                        ));
        }
    }
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
