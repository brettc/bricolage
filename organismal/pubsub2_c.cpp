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

void cNetwork::cycle(cChannelState &c)
{
    cChannelState next(c.size());
    
    for (auto &g: genes)
    {
        for (auto &cis: g.modules)
        {
            if (cis.active(c))
            {
                next.set(g.pub);
                // just one is sufficient
                break;
            }
        }
    }
    // "return" this value
    c.swap(next);
}

cFactory::cFactory(size_t seed)
    : random_engine(seed)
    , next_identifier(0)
    , pop_count(1)
    , gene_count(4)
    , cis_count(3)
{
}

void cFactory::init_environments()
{
    total_channels = cue_channels + reg_channels + out_channels;
    // Environments
    size_t env_count = 1 << cue_channels;
    for (size_t i = 0; i < env_count; ++i)
    {
        cChannelState c = cChannelState(total_channels, i);
        environments.push_back(c);
    }
}

