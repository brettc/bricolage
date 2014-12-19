#include "threshold3.hpp"

using namespace thresh3;

cConstructor::cConstructor(const pubsub2::cWorld_ptr &w, size_t gc, size_t cc)
    : pubsub2::cConstructor(w)
    , r_binding(-3, 3)
    , r_direction(0, 1)
    , r_site(0, 2)
    , r_input(0, w->sub_range.second-1)
    , xxx(boost::bind(pubsub2::randint_t(-3, 3), w->rand))
{
}

pubsub2::cNetwork *cConstructor::construct()
{
    pubsub2::cNetwork *net = new cNetwork(world);
    // net->calc_attractors();
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

cNetwork::cNetwork(const pubsub2::cWorld_ptr &w, bool no_ident)
    : pubsub2::cNetwork(w, no_ident)
{
}

pubsub2::cNetwork *cNetwork::clone() const
{
    return 0;
}

void cNetwork::cycle(pubsub2::cChannelState &c) const
{
}
    
