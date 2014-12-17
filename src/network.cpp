#include "core.hpp"
#include "algorithm.hpp"
#include <cmath>
#include <stdexcept>

using namespace pubsub2;

signal_t cCisModule::get_site_channel(size_t index) const
{
    // TODO: Fix runtime error?
    if (index >= site_count())
        throw std::runtime_error("illegal site index");

    return _channels[index];
}

signal_t cCisModule::set_site_channel(size_t index, signal_t c)
{
    // TODO: Fix runtime error?
    if (index >= site_count())
        throw std::runtime_error("illegal site index");

    signal_t old = _channels[index];
    _channels[index] = c;
    return old;
}

cGene::cGene(sequence_t seq, signal_t p)
    : sequence(seq)
    , intervene(INTERVENE_NONE)
    , pub(p)
{
}

cGene *cGene::clone()
{
    cGene *g = new cGene(sequence, pub);
    for (auto m : modules)
        g->modules.push_back(m->clone());
    return g;
}

cGene::~cGene()
{
    for (auto m : modules)
        delete m;
}

cNetwork::cNetwork(const cFactory_ptr &f, bool no_ident)
    : factory(f)
    , parent_identifier(-1)
    , target(-1)
    , pyobject(0)
{
    if (!no_ident)
        identifier = factory->get_next_network_ident();
    else
        // A network that is not part of a lineage, for analysis only.
        // We'll call this a "detached" network
        identifier = -1;
}

cNetwork::~cNetwork()
{
    for (auto g : genes)
        delete g;
}

void cNetwork::clone_genes(cGeneVector &gv) const
{
    for (auto g : genes)
        gv.push_back(g->clone());
}

// This is the inner-inner loop!
// TODO: Maybe this could be moved down the the CIS level to prevent constant
// calling of virtual function. It would have to be Factory level call:
// cGeneFactory::cycle(Network &, ChannelState &). But this would mean building 
// the cis action into the factory somehow (or static_cast-ing the CIS which
// would, I guess, be safe.
void cNetwork::cycle(cChannelState &c) const
{
    Cycle<cNetwork> cycler;
    // cycler.cycle(*this, c);
    cycler.cycle_with_intervention(*this, c);
}

// This is the outer-inner loop, where we find the attractors. 
// TODO: Maybe template-ize this so that it runs faster without constract calls
// to the virtual function. Like this maybe:
// template class attractor_calc<Cis_Type, Factory_Type>
// {
// }
// But it would still need a dynamic runtime selector to call the right one...
// Hmmm.
void cNetwork::calc_attractors()
{
    attractors.clear();
    rates.clear();

    size_t attractor_begins_at;
    bool found;

    // Go through each environment.
    for (auto &env : factory->environments)
    {
        // Set the state to current environment, and set it as the start of the
        // path to the attractor.
        cChannelStateVector path;
        cChannelState current = env;
        path.push_back(current);
        
        for (;;)
        {
            // Update the current state.
            cycle(current);

            // Put back the environment (as this remains constant)
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
        attractors.push_back(cChannelStateVector());
        cChannelStateVector &this_attr = attractors.back();

        rates.push_back(cRates());
        cRates &this_rate = rates.back();
        for (size_t i=0; i < factory->out_channels; ++i)
            this_rate.push_back(0.0);

        // Copy the part the path that is the attractor, ignoring the transient
        for (size_t copy_at=attractor_begins_at; copy_at < path.size(); ++copy_at)
        {
            cChannelState &c = path[copy_at];
            this_attr.push_back(c);

            // We construct the rates at the same time
            for (size_t i=0; i < factory->out_channels; ++i)
                this_rate[i] += double(c[i + factory->out_range.first]);

        }
        // Now normalise the rates
        for (size_t i=0; i < factory->out_channels; ++i)
            this_rate[i] /= double(this_attr.size());

    }
}


// cNetwork_ptr cNetwork::get_detached_copy() const
// {
//     cNetwork *copy = new cNetwork(factory, true);
//     copy->parent_identifier = identifier;
//     copy->genes = genes;
//     return cNetwork_ptr(copy);
// }
