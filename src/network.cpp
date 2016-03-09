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
    attractors.clear();
    rates.clear();

    InputType it = world->input_type;

    size_t attractor_begins_at;
    bool found;

    // Go through each environment.
    for (auto &env : world->environments)
    {
        // Set the state to current environment, and set it as the start of the
        // path to the attractor.
        cChannelStateVector path;
        cChannelState current = env;
        path.push_back(current);

        for (;;)
        {
            // Update the current state.
            if (!intervention)
                cycle(current);
            else
                cycle_with_intervention(current);

            // Put back the environment if it is constant
            if (it == INPUT_CONSTANT)
                current |= env;

            // Have we already seen this?
            attractor_begins_at = 0;
            found = false;
            for (cChannelState &prev : path)
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
        attractors.emplace_back();
        cChannelStateVector &this_attr = attractors.back();

        rates.emplace_back();
        cRates &this_rate = rates.back();
        for (size_t i=0; i < world->out_channels; ++i)
            this_rate.push_back(0.0);

        // Copy the part the path that is the attractor, ignoring the transient
        for (size_t copy_at=attractor_begins_at; copy_at < path.size(); ++copy_at)
        {
            cChannelState &c = path[copy_at];
            this_attr.push_back(c);

            // We construct the rates at the same time
            for (size_t i=0; i < world->out_channels; ++i)
                this_rate[i] += double(c[i + world->out_range.first]);

        }
        // Now normalise the rates
        for (size_t i=0; i < world->out_channels; ++i)
            this_rate[i] /= double(this_attr.size());

    }
}

typedef std::function<int()> random_int_t;
inline random_int_t random_int_range(int a, int b, const bricolage::cWorld_ptr &w)
{
    return std::bind(std::uniform_int_distribution<>(a, b-1), std::ref(w->rand));
}

void cNetwork::calc_perturbation()
{
    // We should already have initial attractors and rates.  This is NOT
    // deterministic (unlike the calculation of attractors, as we randomize
    // the inputs. It also assumes we'd doing *pulses* of inputs rather than
    // constant inputs.
    pert_attractors.clear();
    pert_rates.clear();
    random_int_t r_env_state(random_int_range(0, 2, world));

    size_t attractor_begins_at;
    bool found;

    // Go through each environment.
    for (auto &attr : attractors)
    {
        cChannelStateVector path;
        // For each environmental attractor:
        // 1. Randomly select one state of the attractor.
        // 2. Introduce a pulse of random env states.
        // 3. Iterator until we get another atttractor.
        //
        // 1.
        size_t state_i = 0;
        if (attr.size() > 1)
        {
            random_int_t r_state(random_int_range(0, attr.size(), world));
            state_i = r_state();
        }
        cChannelState current = attr[state_i];

        // 2.
        for (size_t k = world->cue_range.first; k < world->cue_range.second; ++k)
            if (r_env_state())
                current.set(k);

        path.push_back(current);

        for (;;)
        {
            cycle(current);

            // Have we already seen this?
            attractor_begins_at = 0;
            found = false;
            for (cChannelState &prev : path)
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
        pert_attractors.emplace_back();
        cChannelStateVector &this_attr = pert_attractors.back();

        pert_rates.emplace_back();
        cRates &this_rate = pert_rates.back();
        for (size_t i = 0; i < world->out_channels; ++i)
            this_rate.push_back(0.0);

        // Copy the part the path that is the attractor, ignoring the transient
        for (size_t copy_at=attractor_begins_at; copy_at < path.size(); ++copy_at)
        {
            cChannelState &c = path[copy_at];
            this_attr.push_back(c);

            // We construct the rates at the same time
            for (size_t i=0; i < world->out_channels; ++i)
                this_rate[i] += double(c[i + world->out_range.first]);

        }
        // Now normalise the rates
        for (size_t i=0; i < world->out_channels; ++i)
            this_rate[i] /= double(this_attr.size());

    }
}


