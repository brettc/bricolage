#include "core.hpp"
#include <cmath>
#include <stdexcept>

using namespace pubsub2;

cFactory::cFactory(size_t seed, size_t cue, size_t reg, size_t out)
    // Key stuff
    : next_network_identifier(0)
    , next_target_identifier(0)
    // Channels
    , cue_channels(cue)
    , reg_channels(reg)
    , out_channels(out)
    // Random 
    , rand(seed)
    , constructor(0)
{
    init_channels();
    init_environments();
}

cFactory::~cFactory()
{
    if (constructor)
        delete constructor;
}

void cFactory::init_environments()
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

void cFactory::init_channels()
{
    // Calculate the total number of elements given the overlap
    // Given cue = 2, reg = 2, out = 2
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

cConstructor::cConstructor(cFactory &f, size_t gene_count_, size_t cis_count_)
    : factory(f)
    , gene_count(gene_count_)
    , cis_count(cis_count_)
    , r_gene(0, gene_count-1)
    , r_cis(0, cis_count-1)
    , r_sub(f.sub_range.first, f.sub_range.second-1)
{
}

void cConstructor::construct_network(cNetwork &network)
{
    for (size_t i=0; i < gene_count; ++i)
    {
        // TODO: genes.emplace_back(...)
        cGene *g = new cGene(i, factory.pub_range.first + i);
        network.genes.push_back(g);

        for (size_t j=0; j < cis_count; ++j)
            g->modules.push_back(construct_cis());
    }

    // Calculate the attractors
    network.calc_attractors();
}

void cConstructor::mutate_network(cNetwork_ptr &n, size_t mutations)
{
    // NOTE: This done INPLACE mutation. It should never be called on a network
    // that has already had its attractors calculated! 
    
    // Select the genes that should be mutated
    while (mutations > 0)
    {
        size_t i = r_gene(factory.rand);
        mutate_gene(n->genes[i]);
        --mutations;
    }
}

cNetwork_ptr cConstructor::copy_and_mutate_network(cNetwork_ptr &n, size_t mutations)
{
    cNetwork_ptr copy(new cNetwork(n->factory));
    // The copy constructor of vector does the hard work here...
    copy->parent_identifier = n->identifier;
    n->clone_genes(copy->genes);

    // Now mutate it and calculate the attractors
    mutate_network(copy, mutations);
    copy->calc_attractors();
    return copy;
}

typedef std::poisson_distribution<> poisson_t;

void cConstructor::mutate_collection(cNetworkVector &networks, cIndexes &mutated,
                                  double site_rate)
{
    // How many mutations are we going to have? That depends on the total
    // number of sites that might mutate. Per network, this is rate_per_gene *
    // gene_count. We multiply this by the number of networks in the collection
    // to get the expected number of mutations. Then we generate a number using
    // a poisson distribution.
    // std::cout << "rate per" << rate_per_gene;
    double expected = site_rate * networks.size() * gene_count;
    poisson_t r_pop(expected);
    size_t mutations = r_pop(factory.rand);

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
        mutes.push_back(r_network(factory.rand));

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
            networks[network_num] = copy_and_mutate_network(networks[network_num], count);

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

void cConstructor::mutate_gene(cGene *g)
{
    if (g->modules.size() == 0)
        throw std::runtime_error("no cis modules");

    // First, decide what cis module we're using. Account for the fact that
    // there might be only one cis module.
    size_t cis_i = 0; 

    // More than? We need to pick one...
    if (cis_count > 0)
        cis_i = r_cis(factory.rand);

    // Grab this and mutate it.
    cCisModule *m = g->modules[cis_i];
    mutate_cis(m);
}
