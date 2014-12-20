#include "threshold3.hpp"
#include "algorithm.hpp"

using namespace thresh3;

cConstructor::cConstructor(const pubsub2::cWorld_ptr &w, size_t gc, size_t cc)
    : pubsub2::cConstructor(w)
    , gene_count(gc)
    , module_count(cc)
    , r_binding(-3, 3)
    , r_direction(0, 1)
    , r_site(0, 2)
    , r_input(0, w->sub_range.second-1)
    , xxx(boost::bind(pubsub2::randint_t(-3, 3), w->rand))
{
}

pubsub2::cNetwork *cConstructor::construct()
{
    pubsub2::cNetwork *net = new cNetwork(*this);
    return net;
}

cCisModule::cCisModule(const cConstructor &c)
{
    for (size_t i = 0; i < 3; ++i)
    {
        channels[i] = c.r_input(c.world->rand);
        binding[i] = c.r_binding(c.world->rand);
    }
}

// This is where the action really is.
void cCisModule::mutate(const cConstructor &c)
{
    size_t site = c.r_site(c.world->rand);
    pubsub2::int_t current = binding[site];

    pubsub2::int_t mutate;
    if (current == 3)
        mutate = -1;
    else if (current == -3)
        mutate = 1;
    else
        mutate = c.r_direction(c.world->rand) * 2 - 1;

    // If we're at zero, possibly change into a different binding 
    if (current == 0)
        channels[site] = c.r_input(c.world->rand);

    binding[site] += mutate;
}

bool cCisModule::is_active(pubsub2::cChannelState const &state) const 
{
    // Calculate the weighted sum
    pubsub2::int_t sum = 0;
    for (size_t i = 0; i < site_count(); ++i)
    {
        if (state[channels[i]])
            sum += binding[i];
    }
    // Thresholded
    return sum >= 3;
}

cNetwork::cNetwork(const cConstructor &c)
    : pubsub2::cNetwork(c.world)
    , constructor(c)
{
    for (size_t i=0; i < c.gene_count; ++i)
    {
        genes.emplace_back(i, c.world->pub_range.first + i);
        auto &g = genes.back();

        for (size_t j=0; j < c.module_count; ++j)
            g.modules.emplace_back(c);
    }

    // Calculate the attractors
    calc_attractors();
}

pubsub2::cNetwork *cNetwork::clone() const
{
    return 0;
}

// This is the inner-inner loop!
// TODO: Maybe this could be moved down the the CIS level to prevent constant
// calling of virtual function. It would have to be Factory level call:
// cGeneFactory::cycle(Network &, ChannelState &). But this would mean building 
// the cis action into the world somehow (or static_cast-ing the CIS which
// would, I guess, be safe.
void cNetwork::cycle(pubsub2::cChannelState &c) const
{
    algo::Cycle<cNetwork> cycler;
    cycler.cycle(*this, c);
}

// A slower version with the ability intervene
void cNetwork::cycle_with_intervention(pubsub2::cChannelState &c) const
{
    algo::Cycle<cNetwork> cycler;
    cycler.cycle_with_intervention(*this, c);
}
   
cGene::cGene(pubsub2::sequence_t sequence, pubsub2::signal_t p)
    : pubsub2::cGene(sequence, p)
{
}
