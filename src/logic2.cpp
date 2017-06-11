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
    : bricolage::cFactory(w, cc)
    , operands(ops)
    , r_operand(random_int_range(0, ops.size(), w))
    , r_site(random_int_range(0, 2, w))
{
}

bricolage::cNetwork_ptr cFactory::construct(bool create)
{
    bricolage::cFactory_ptr p = shared_from_this();
    cNetwork *net = new cNetwork(p);
    if (create)
        net->identifier = world->get_next_network_ident();

    size_t pub;
    for (size_t gindex = 0; gindex < gene_count; ++gindex)
    {
        // Figure out the next publish value 
        if (gindex < world->reg_gene_count)
        {
            pub = world->reg_range.first + (gindex % world->reg_channels);
        }
        else
        {
            pub = gindex - world->reg_gene_count + world->out_range.first;
        }

        net->genes.emplace_back(gindex, pub);
        auto &g = net->genes.back();

        for (size_t j=0; j < module_count; ++j)
        {
            g.modules.emplace_back(cCisModule(*this));
            if (create)
            {
                auto &m = g.modules.back();
                m.op = operands[r_operand()];
                m.channels[0] = draw_from_subs[r_sub()];
                m.channels[1] = draw_from_subs[r_sub()];
            }
        }

    }

    // Calculate the attractors
    if (create)
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

void cNetwork::mutate(size_t n_cis, size_t n_trans)
{
    static algo::Mutator<cNetwork, cFactory> mutator;
    mutator.mutate(*this, n_cis, n_trans);
}

void cNetwork::duplicate(size_t n_dups)
{
    static algo::Mutator<cNetwork, cFactory> mutator;
    mutator.duplicate(*this, n_dups);
}

bricolage::cNetwork_ptr cNetwork::clone() const
{
    // We don't use the construct here -- as we're copying
    cNetwork *copy = new cNetwork(factory);

    // This is the only extra things that needs copying.
    // Everything magically works here.
    copy->genes = genes;
    // copy->identifier = world->get_next_network_ident();
    // copy->parent_identifier = identifier;

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
    channels[i] = c.draw_from_subs[c.r_sub()];
    op = c.operands[c.r_operand()];
}
