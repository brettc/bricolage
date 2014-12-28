#include "core.hpp"
#include <cmath>
#include <stdexcept>

using namespace pubsub2;

cWorld::cWorld(size_t seed, size_t cue, size_t reg, size_t out)
    // Key stuff
    : next_network_identifier(0)
    , next_target_identifier(0)
    // Channels
    , cue_channels(cue)
    , reg_channels(reg)
    , out_channels(out)
    // Random 
    , rand(seed)
{
    init_channels();
    init_environments();
}

cWorld::~cWorld()
{
}

void cWorld::init_environments()
{
    // Number of environments is 2^cue_channels. Use some binary math...
    size_t env_count = 1 << cue_channels;

    for (size_t i = 0; i < env_count; ++i)
    {
        // Shift one, to account for channel 0
        cChannelState c = cChannelState(channel_count, i << reserved_channels);
        // Turn on bias channel
        c.set(on_channel); 
        environments.push_back(c);
    }
}

void cWorld::init_channels()
{
    // Calculate the total number of elements given the overlap
    // Example: Given cue = 2, reg = 2, out = 2
    //             0 = ALWAYS OFF
    //               1 = ALWAYS ON
    // reserved  [ 0 1 ]
    // cues          [ 2 3 ]
    // regs              [ 4 5 ] 
    // outs                  [ 6 7 8 ]
    //
    // subs          [ 2 3 4 5 ]
    // pubs              [ 4 5 6 7 8 ]
    channel_count = cue_channels + reg_channels + out_channels + reserved_channels;

    // These are python-like *ranges*, thus the interval is [first, second) or
    // first <= v < second
    cue_range.first = reserved_channels;
    cue_range.second = cue_range.first + cue_channels;

    reg_range.first = cue_range.second;
    reg_range.second = reg_range.first + reg_channels;

    out_range.first = reg_range.second;
    out_range.second = out_range.first + out_channels;

    sub_range.first = cue_range.first;
    sub_range.second = reg_range.second;

    pub_range.first = reg_range.first;
    pub_range.second = out_range.second;
}

cConstructor::cConstructor(const cWorld_ptr &w)
    : world(w)
{
}

cNetwork_ptr cConstructor::clone_and_mutate_network(cNetwork_ptr &n, size_t nmutations)
{
    cNetwork_ptr copy(n->clone());
    copy->parent_identifier = n->identifier;
    copy->mutate(nmutations);
    copy->calc_attractors();
    return copy;
}

void cConstructor::mutate_collection(
    cNetworkVector &networks, cIndexes &mutated, double site_rate)
{
    // How many mutations are we going to have? That depends on the total
    // number of sites that might mutate. Per network, this is rate_per_gene *
    // gene_count. We multiply this by the number of networks in the collection
    // to get the expected number of mutations. Then we generate a number using
    // a poisson distribution.
    double expected = site_rate * site_count(networks);
    std::poisson_distribution<> r_pop(expected);
    size_t mutations = r_pop(world->rand);

    // Clear this
    mutated.clear();

    // If we're not generating any mutations, let's just bail.
    if (mutations == 0)
        return;

    // Now we need to assign these to individual networks. Only these networks 
    // will change. The rest remain constant.
    randint_t r_network(0, networks.size()-1);
    cIndexes mutes;
    for (size_t i=0; i < mutations; ++i)
        mutes.push_back(r_network(world->rand));

    // Sort them so that any repeated networks are adjacent.
    std::sort(mutes.begin(), mutes.end());

    // We now let the networks figure out how *exactly* they will mutate.  This
    // looks longwinded, but we want to handle cases where a single network is
    // mutated more than once in one step, rather than calling multiple times
    // on the same network, so most of this is putting together multiple
    // mutations into a single call.
    auto it = mutes.begin();
    size_t network_num = *it, count = 1;
    for (;;) 
     {
        // Skip ahead to the next one
        ++it;

        // If we're at the end, or the next one is different, go ahead and
        // apply the mutations ...
        if (it == mutes.end() || *it != network_num)
        {
            networks[network_num] = clone_and_mutate_network(networks[network_num], count);

            // Add to the indexes that changed
            mutated.push_back(network_num);
            if (it == mutes.end())
                break;

            count = 1;
            network_num = *it;
        }
        // ... otherwise we have a repeat of the same network_num. Let's just
        // accumulate this into one call.
        else
        {
            count += 1;
        }
    }
}

