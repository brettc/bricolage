#include "threshold3.hpp"
#include "algorithm.hpp"

using namespace thresh3;

cConstructor::cConstructor(const pubsub2::cWorld_ptr &w, size_t cc)
    : pubsub2::cConstructor(w)
    , gene_count(w->reg_channels + w->out_channels)
    , module_count(cc)
    , r_binding(-3, 3)
    , r_direction(0, 1)
    , r_site(0, 2)
    , r_input(0, w->sub_range.second-1)
    , r_gene(boost::bind(pubsub2::randint_t(0, gene_count-1), w->rand))
    , r_cis(boost::bind(pubsub2::randint_t(0, cc-1), w->rand))
{
}

// Construct a brand new Network with random stuff.
pubsub2::cNetwork_ptr cConstructor::construct()
{
    pubsub2::cConstructor_ptr p = shared_from_this();
    cNetwork *net = new cNetwork(p);

    for (size_t pub = world->reg_range.first, gindex = 0; 
         pub < world->pub_range.second; ++pub, ++gindex)
    {
        net->genes.push_back(cGene(gindex, pub));
        auto &g = net->genes.back();

        for (size_t j=0; j < module_count; ++j)
            g.modules.push_back(cCisModule(*this));
    }

    // Calculate the attractors
    net->calc_attractors();
    return pubsub2::cNetwork_ptr(net);
}

size_t cConstructor::site_count(pubsub2::cNetworkVector &networks)
{
    // TODO: should multiple by 3!!!
    // Just keeping it this way for comparison
    return gene_count * module_count * networks.size();
}

void cNetwork::mutate(size_t nmutations)
{
    auto &ctor = static_cast<const cConstructor &>(*constructor);

    // Select the genes that should be mutated
    while (nmutations > 0)
    {
        // Choose a gene and mutate it
        size_t i = ctor.r_gene();
        auto &g = genes[i];

        size_t j = ctor.r_cis();
        auto &c = g.modules[j];
        c.mutate(ctor);
        --nmutations;
    }
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
// void cCisModule::mutate(const cConstructor &c)
// {
//     size_t site = c.r_site(c.world->rand);
//     pubsub2::int_t current = binding[site];
//
//     pubsub2::int_t mutate;
//     if (current == 3)
//         mutate = -1;
//     else if (current == -3)
//         mutate = 1;
//     else
//         mutate = c.r_direction(c.world->rand) * 2 - 1;
//
//     // If we're at zero, possibly change into a different binding 
//     if (current == 0)
//         channels[site] = c.r_input(c.world->rand);
//
//     binding[site] += mutate;
// }

// This is where the action really is.
void cCisModule::mutate(const cConstructor &c)
{
    // TODO: mutation size is ....
    // 1 + poisson something?
    size_t i = c.r_site(c.world->rand);
    channels[i] = c.r_input(c.world->rand);
    binding[i] = c.r_binding(c.world->rand);
}

bool cCisModule::is_active(pubsub2::cChannelState const &state) const 
{
    // Calculate the weighted sum. Unrolled.
    pubsub2::int_t sum = 0;
    if (state[channels[0]])
        sum += binding[0];
    if (state[channels[1]])
        sum += binding[1];
    if (state[channels[2]])
        sum += binding[2];
    // Thresholded
    return sum >= 3;
}

cNetwork::cNetwork(const pubsub2::cConstructor_ptr &c)
    : pubsub2::cNetwork(c)
{
}

pubsub2::cNetwork_ptr cNetwork::clone() const
{
    // We don't use the construct here -- as we're copying
    cNetwork *copy = new cNetwork(constructor);

    // This is the only extra things that needs copying.
    // Everything magically works here.
    copy->genes = genes;

    // We also don't calculate attractors or anything, as we might be doing
    // some mutating first.
    return pubsub2::cNetwork_ptr(copy);
}

// This is the outer-inner loop!
void cNetwork::cycle(pubsub2::cChannelState &c) const
{
    static algo::Cycle<cNetwork> cycler;
    cycler.cycle(*this, c);
}

// A slower version with the ability intervene
void cNetwork::cycle_with_intervention(pubsub2::cChannelState &c) const
{
    static algo::Cycle<cNetwork> cycler;
    cycler.cycle_with_intervention(*this, c);
}
   
cGene::cGene(pubsub2::sequence_t sequence, pubsub2::signal_t p)
    : pubsub2::cGene(sequence, p)
{
}

// void cConstructor::mutate_gene(cGene *g)
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

