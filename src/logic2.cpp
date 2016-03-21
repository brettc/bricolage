#include "logic2.hpp"
#include "algorithm.hpp"
#include <cmath>
#include <stdexcept>
#include <iostream>

using namespace logic2;

// --------------------------------------------------------------------
// cFactory
cFactory::cFactory(const bricolage::cWorld_ptr &w, size_t cc, 
                           const cOperands &ops)
    : bricolage::cFactory(w)
    , gene_count(w->reg_channels + w->out_channels)
    , module_count(cc)
    , operands(ops)
    , r_gene(random_int_range(0, gene_count, w))
    , r_module(random_int_range(0, module_count, w))
    , r_operand(random_int_range(0, ops.size(), w))
    , r_site(random_int_range(0, 2, w))
    , r_input(random_int_range(w->sub_range.first, w->sub_range.second, w))
{

    // Randomly allocate operands to binding pairs
    // std::pair<size_t, size_t> &r = w->sub_range;
    // for (bricolage::signal_t a = r.first; a < r.second; ++a)
    //     for (bricolage::signal_t b = r.first; b < r.second; ++b)
    //         bindings[std::make_pair(a, b)] = r_operand();
}

bricolage::cNetwork_ptr cFactory::construct(bool fill)
{
    bricolage::cFactory_ptr p = shared_from_this();
    cNetwork *net = new cNetwork(p);
    if (fill)
        net->identifier = world->get_next_network_ident();

    for (size_t pub = world->reg_range.first, gindex = 0; 
         pub < world->pub_range.second; ++pub, ++gindex)
    {
        net->genes.emplace_back(gindex, pub);
        auto &g = net->genes.back();

        for (size_t j=0; j < module_count; ++j)
        {
            g.modules.emplace_back(cCisModule(*this));
            if (fill)
            {
                auto &m = g.modules.back();
                m.op = operands[r_operand()];
                m.channels[0] = r_input();
                m.channels[1] = r_input();
            }
        }
    }

    // Calculate the attractors
    if (fill)
        net->calc_attractors();
    return bricolage::cNetwork_ptr(net);
}

size_t cFactory::site_count(bricolage::cNetworkVector &networks)
{
    return gene_count * module_count * 2 * networks.size();
}

// --------------------------------------------------------------------
// cNetwork
// TODO: Make this a template
cNetwork::cNetwork(const bricolage::cFactory_ptr &c)
    : bricolage::cNetwork(c)
{
    // nothing special here...
}

void cNetwork::mutate(size_t nmutations)
{
    auto &ctor = static_cast<const cFactory &>(*factory);

    // Select the genes that should be mutated
    while (nmutations > 0)
    {
        // Choose a gene and mutate it
        size_t i = ctor.r_gene();
        auto &g = genes[i];
        
        size_t j = ctor.r_module();
        auto &c = g.modules[j];
        c.mutate(ctor);
        --nmutations;
    }
}

bricolage::cNetwork_ptr cNetwork::clone() const
{
    // We don't use the construct here -- as we're copying
    cNetwork *copy = new cNetwork(factory);

    // This is the only extra things that needs copying.
    // Everything magically works here.
    copy->genes = genes;

    // We also don't calculate attractors or anything, as we might be doing
    // some mutating first.
    return bricolage::cNetwork_ptr(copy);
}

// This is the outer-inner loop!
// TODO: This should be part of a template
void cNetwork::cycle(bricolage::cChannels &c) const
{
    static algo::Cycle<cNetwork> cycler;
    cycler.cycle(*this, c);
}

// A slower version with the ability intervene
void cNetwork::cycle_with_intervention(bricolage::cChannels &c) const
{
    static algo::Cycle<cNetwork> cycler;
    cycler.cycle_with_intervention(*this, c);
}

// --------------------------------------------------------------------
// cGene
cGene::cGene(bricolage::sequence_t sequence, bricolage::signal_t p)
    : bricolage::cGene(sequence, p)
{
}

// --------------------------------------------------------------------
// CisModule
cCisModule::cCisModule(const cFactory &c)
{
}

// This is where the action really is.
void cCisModule::mutate(const cFactory &c)
{
    // Pick a channel
    size_t i = c.r_site();
    channels[i] = c.r_input();
    // signal_pair_t cc = std::make_pair(channels[0], channels[1]);
    // op = c.bindings.at(cc);
    op = c.operands[c.r_operand()];
}
