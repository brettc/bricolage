#include "core.hpp"
#include "algorithm.hpp"
#include <cmath>
#include <stdexcept>
#include <iostream>

using namespace bricolage;

cGene::cGene(sequence_t seq, signal_t p)
    : sequence(seq)
    , pub(p)
    , intervene(INTERVENE_NONE)
{
}

cNetwork_ptr cFactory::clone_and_mutate_network(
    cNetwork_ptr &n, size_t n_sub, size_t n_pub, size_t dups, int_t generation)
{
    cNetwork_ptr copy(n->clone());
    copy->identifier = world->get_next_network_ident();
    copy->parent_identifier = n->identifier;
    copy->generation = generation;
    if (n_sub || n_pub)
        copy->mutate(n_sub, n_pub);
    if (dups)
        copy->duplicate(dups);
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
    for (auto const &external : world->inputs)
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

        stabilise(external, intervention, this_attr, this_trans, this_rate);
    }
}


// This is the core calculation for everything (at the moment)
void cNetwork::stabilise(const cAttractor &external,
                         bool intervention,
                         cAttractor &attractor_,
                         cAttractor &transient_,
                         cRates &rates_) const
{
    cAttractor path;
    size_t external_index = 0;
    size_t external_size = external.size();

    // Initialise the states
    cChannels current = external[external_index];
    external_index++;

    // Make sure the on channel is on
    current.unchecked_set(on_channel);

    path.push_back(current);

    size_t attractor_begins_at;
    bool found;

    // This is a bit repetitive
    for (;;)
    {
        // Update the current state.
        if (!intervention)
            cycle(current);
        else
            cycle_with_intervention(current);

        // Make sure the on channel is always on
        current.unchecked_set(on_channel);

        // Mix in the current external states 
        if (external_index < external_size)
        {
            current.unchecked_union(external[external_index]);
            external_index++;
        }
        else
        {
            // Okay, all external changes have been seen. Now we can look for
            // the attractor. Only begin looking at the last time we injected
            // something.
            attractor_begins_at = external_size - 1;
            found = false;

            // We're looking to see if the very same set of channels has
            // occured. If so, we about to into a loop.
            while (attractor_begins_at < path.size())
            {
                cChannels &prev = path[attractor_begins_at];
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
        }

        // Add the current to our attractor.
        path.push_back(current);
    }

    // Add a new attractor for this environment.
    attractor_.clear();
    transient_.clear();
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
    
    // Fill cache (if we're not intervening) -- everything in transient and
    // attractor map to rates
    if (!intervention)
    {
        for (auto &ch: path)
            cached_mappings.emplace(ch, rates_);
    }
}


void cNetwork::get_rates(const cChannels &initial, cRates &rates, bool use_cache) const
{
    cAttractor external;
    std::fill_n(std::back_inserter(external), world->pulse_for, initial);
    
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
    stabilise(external, false, attr, trans, rates);
}

double cNetwork::attractor_robustness() const
{
    cRates new_rates;
    cAttractor new_attr, new_trans;
    double success = 0.0;
    double per_channel = 1.0 / 
        double((world->cue_channels + world->reg_channels) * 
        attractors.size());
   
    for (size_t i = 0; i < attractors.size(); ++i)
    {
        const auto &attr = attractors[i];
        const auto &orig_rates = rates[i];
        double per_test = per_channel / double(attr.size());
        for (const auto &attr_state: attr)
        {
            for (size_t j = world->cue_range.first; 
                 j < world->reg_range.second;
                 ++j)
            {
                cChannels c(attr_state);
                c.unchecked_flip(j);
                get_rates(c, new_rates, true);
                bool same = true;
                for (size_t k = 0; k < orig_rates.size(); ++k)
                {
                    if (!is_close(orig_rates[k], new_rates[k]))
                    {
                        same = false;
                        break;
                    }
                }
                if (same)
                    success += per_test;
            }
        }
    }

    return success;
}

