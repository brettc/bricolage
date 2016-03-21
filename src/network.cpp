#include "core.hpp"
#include "algorithm.hpp"
#include <cmath>
#include <stdexcept>
#include <boost/function.hpp>
#include <boost/bind.hpp>

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

// TODO: move to main file
typedef std::function<int()> random_int_t;
inline random_int_t random_int_range(int a, int b, const bricolage::cWorld_ptr &w)
{
    return std::bind(std::uniform_int_distribution<>(a, b-1), std::ref(w->rand));
}

void cNetwork::calc_perturbation(cDynamics &dynamics, bool env_only) const
{
    // We should already have initial attractors and rates.  This is NOT
    // deterministic, unlike the calculation of attractors, as we randomize
    // the inputs. It also assumes we'd doing *pulses* of inputs rather than
    // constant inputs.
    dynamics.clear();
    random_int_t r_env_state(random_int_range(0, 2, world));
    random_int_t r_reg_channel(random_int_range(world->cue_range.first,
                                                world->reg_range.second,
                                                world));

    // Go through each environment.
    for (auto &attr : attractors)
    {
        // For each environmental attractor:
        // 1. Randomly select one state of the attractor.
        // 2. Introduce a pulse of random states (maybe just in the env)
        // 3. Find new attractor
        //
        // 1. Get starting point
        size_t state_i = 0;
        if (attr.size() > 1)
        {
            random_int_t r_state(random_int_range(0, attr.size(), world));
            state_i = r_state();
        }
        cChannels start_state = attr[state_i];

        // 2. Random state
        if (env_only)
        {
            for (size_t k = world->cue_range.first; k < world->cue_range.second; ++k)
                if (r_env_state())
                    start_state.unchecked_set(k);
        } else {
            start_state.unchecked_flip(r_reg_channel());
        }

        // Add a new attractor for this environment.
        dynamics.add();
        auto &this_attractor = dynamics.attractors.back();
        auto &this_transient = dynamics.transients.back();
        auto &this_rate = dynamics.rates.back();

        stabilise(start_state, false, this_attractor, this_transient, this_rate);
    }
}

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
}

