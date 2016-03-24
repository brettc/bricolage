#include "core.hpp"
#include "algorithm.hpp"
#include <cmath>
#include <stdexcept>

using namespace bricolage;

cGene::cGene(sequence_t seq, signal_t p)
    : sequence(seq)
    , pub(p)
    , intervene(INTERVENE_NONE)
{
}

cNetwork_ptr cFactory::clone_and_mutate_network(
    cNetwork_ptr &n, size_t nmutations, int_t generation)
{
    cNetwork_ptr copy(n->clone());
    copy->identifier = world->get_next_network_ident();
    copy->parent_identifier = n->identifier;
    copy->generation = generation;
    copy->mutate(nmutations);
    copy->calc_attractors();
    return copy;
}

cNetwork::cNetwork(const cFactory_ptr &c)
    : factory(c)
    , world(c->world)
    , identifier(-1)
    , parent_identifier(-1)
    , generation(0)
    , target(-1)
    , pyobject(0)
{
    // identifier = world->get_next_network_ident();
}

// This is the outer-inner loop, where we find the attractors.
void cNetwork::_calc_attractors(bool intervention)
{
    transients.clear();
    attractors.clear();
    rates.clear();

    // Go through each environment.
    for (auto &env : world->environments)
    {
        // Set the state to current environment, and set it as the start of the
        // path to the attractor.
        // Add a new attractor for this environment.
        attractors.emplace_back();
        cAttractor &this_attr = attractors.back();
        rates.emplace_back();
        cRates &this_rate = rates.back();
        transients.emplace_back();
        cAttractor &this_trans = transients.back();

        stabilise(env, intervention, this_attr, this_trans, this_rate);
    }
}



// This is the core calculation for everything (at the moment)
void cNetwork::stabilise(const cChannels &initial,
                         bool intervention,
                         cAttractor &attractor_,
                         cAttractor &transient_,
                         cRates &rates_) const
{
    cAttractor path;
    cChannels current = initial;

    // Make sure the on channel is on
    current.unchecked_set(on_channel);

    path.push_back(current);

    size_t attractor_begins_at;
    bool found;

    for (;;)
    {
        // Update the current state.
        if (!intervention)
            cycle(current);
        else
            cycle_with_intervention(current);

        // Make sure the on channel is on
        current.unchecked_set(on_channel);

        // Have we already seen this?
        attractor_begins_at = 0;
        found = false;
        for (cChannels &prev : path)
        {
            if (prev == current)
            {
                found = true;
                break;
            }
            attractor_begins_at++;
        }

        // If we have seen this state, we've found the attractor.
        if (found)
            break;

        // Add the current to our attractor.
        path.push_back(current);
    }

    // Add a new attractor for this environment.
    rates_.clear();
    for (size_t i = 0; i < world->out_channels; ++i)
        rates_.push_back(0.0);

    // Copy the attractor and transient
    for (size_t copy_at = 0; copy_at < path.size(); ++copy_at)
    {
        cChannels &c = path[copy_at];
        if (copy_at >= attractor_begins_at)
        {
            attractor_.push_back(c);
            // We construct the rates at the same time
            for (size_t i = 0; i < world->out_channels; ++i)
                if (c.unchecked_test(i + world->out_range.first))
                    rates_[i] += 1.0;
        }
        else
            transient_.push_back(c);
    }

    // Now normalise the rates. We want the *average* rate across the
    // attractor; the stable state expression rate.
    double norm = 1.0 / double(attractor_.size());
    for (size_t i = 0; i < world->out_channels; ++i)
        rates_[i] *= norm;
    
    // Fill cache (if we're not interventing) -- everything in transient and
    // attractor map to rates
    if (!intervention)
    {
        for (auto &ch: path)
        {
            cached_mappings[ch] = rates_;
            // auto ret = cached_mappings.insert(std::make_pair(initial,
            // rates_));

        }
    }
}


void cNetwork::get_rates(const cChannels &initial, cRates &rates, bool use_cache) const
{
    if (use_cache)
    {
        auto found = cached_mappings.find(initial);
        if (found != cached_mappings.end())
        {
            rates = found->second;
            return;
        }
    }

    cAttractor attr, trans;
    stabilise(initial, false, attr, trans, rates);
}
