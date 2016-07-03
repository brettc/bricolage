#include "threshold3.hpp"
#include "algorithm.hpp"

#include <iostream>

using namespace thresh3;

cFactory::cFactory(const bricolage::cWorld_ptr &w, size_t cc,
                           const MutateType mtype)
    : bricolage::cFactory(w, cc)
    , mutate_type(mtype)
    , r_site(random_int_range(0, 3, w))
    , r_binding(random_int_range(-3, 4, w))
    , r_direction(random_int_range(0, 2, w))
{
}

// Construct a brand new Network with random stuff.
bricolage::cNetwork_ptr cFactory::construct(bool fill)
{
    bricolage::cFactory_ptr p = shared_from_this();
    cNetwork *net = new cNetwork(p);
    if (fill)
        net->identifier = world->get_next_network_ident();

    for (size_t pub = world->reg_range.first, gindex = 0;
         pub < world->pub_range.second; ++pub, ++gindex)
    {
        net->genes.push_back(cGene(gindex, pub));
        auto &g = net->genes.back();

        for (size_t j=0; j < module_count; ++j)
        {
            g.modules.push_back(cCisModule(*this));
            if (fill)
            {
                auto &m = g.modules.back();
                for (size_t i = 0; i < 3; ++i)
                {
                    m.channels[i] = draw_from_subs[r_sub()];
                    m.binding[i] = r_binding();
                }
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
    // TODO: should multiple by 3!!!
    // Just keeping it this way for comparison
    return gene_count * module_count * networks.size();
}

void cNetwork::mutate(size_t nmutations)
{
    auto &fy = static_cast<const cFactory &>(*factory);

    // Select the genes that should be mutated
    while (nmutations > 0)
    {
        // Choose a gene and mutate it
        size_t i = fy.r_gene();
        auto &g = genes[i];

        size_t j = fy.r_module();
        auto &c = g.modules[j];
        c.mutate(fy, g);
        --nmutations;
    }
}

void cNetwork::duplicate(size_t nmutations)
{
    auto &ctor = static_cast<const cFactory &>(*factory);

    while (nmutations > 0)
    {
        // Choose a gene and mutate it
        size_t i, j;
        i = ctor.r_regulatory();
        do 
        {
            j = ctor.r_regulatory();
        } while (i == j);

        genes[j] = genes[i];
        --nmutations;
    }
}

cCisModule::cCisModule(const cFactory &fy)
{
}

// This is where the action really is.
void cCisModule::mutate(const cFactory &fy, const cGene &gene)
{
    size_t site = fy.r_site();
    switch (fy.mutate_type)
    {
    case JUMP:
        {
            channels[site] = fy.draw_from_subs[fy.r_sub()];
            binding[site] = fy.r_binding();
        }
        break;
    case PROGRESSIVE:
        {
            bricolage::int_t current = binding[site];

            bricolage::int_t mutate;
            if (current == 3)
                mutate = -1;
            else if (current == -3)
                mutate = 1;
            else
                // Make it either -1 or +1
                mutate = fy.r_direction() * 2 - 1;

            // If we're at zero, possibly change into a different binding
            if (current == 0)
                // Randomly draw from the provided list
                channels[site] = fy.draw_from_subs[fy.r_sub()];

            binding[site] += mutate;
        }
        break;
    }
}

// This is where the action really is.
// void cCisModule::mutate(const cFactory &fy)
// {
//     // TODO: mutation size is ....
//     // 1 + poisson something?
// }

bool cCisModule::is_active(bricolage::cChannels const &state) const
{
    // Calculate the weighted sum. Unrolled.
    bricolage::int_t sum = 0;
    if (state.unchecked_test(channels[0]))
        sum += binding[0];
    if (state.unchecked_test(channels[1]))
        sum += binding[1];
    if (state.unchecked_test(channels[2]))
        sum += binding[2];
    // Thresholded
    return sum >= 3;
}

cNetwork::cNetwork(const bricolage::cFactory_ptr &fy)
    : bricolage::cNetwork(fy)
{
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

cGene::cGene(bricolage::sequence_t sequence, bricolage::signal_t p)
    : bricolage::cGene(sequence, p)
{
}

// void cFactory::mutate_gene(cGene *g)
// {
//     if (g->modules.size() == 0)
//         throw std::runtime_error("no cis modules");
//
//     // First, decide what cis module we're using. Account for the fact that
//     // there might be only one cis module.
//     size_t cis_i = 0;
//
//     // More than? We need to pick one...
//     if (cis_count > 0)
//         cis_i = r_cis(world.rand);
//
//     // Grab this and mutate it.
//     cCisModule *m = g->modules[cis_i];
//     mutate_cis(m);
// }
//

