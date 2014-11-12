// vim: path=.,/usr/include/c++/4.2.1,/usr/include/c++/4.2.1/tr1 fdm=syntax
#include "pubsub2_c.h"

using namespace pubsub2;

cCisModule::cCisModule(operand_t op_, signal_t sub1_, signal_t sub2_)
    : op(op_)
    , sub1(sub1_)
    , sub2(sub2_)
{
}

cGene::cGene(sequence_t seq, signal_t p)
    : sequence(seq)
    , pub(p)
{
}

cNetwork::cNetwork(cFactory_ptr &f)
    : pyobject(0)
    , factory(f)
{
    identifier = factory->get_next_ident();
    parent_identifier = -1;
}

void cFactory::construct_random(cNetwork &network)
{
    randint_t r_sub(sub_range.first, sub_range.second-1);
    randint_t r_oper(0, operands.size()-1);
    
    for (size_t i=0; i < gene_count; ++i)
    {
        network.genes.push_back(cGene(i, pub_range.first + i));
        // genes.emplace_back(i, rpub(rng));
        // Get a reference to it
        cGene &g = network.genes.back();
        for (size_t j=0; j < cis_count; ++j)
        {
            g.modules.push_back(
                cCisModule(operands[r_oper(random_engine)], 
                           r_sub(random_engine), 
                           r_sub(random_engine)
                        ));
        }
    }

    // Calculate the attractors
    network.calc_attractors();
}

void cNetwork::cycle(cChannelState &c)
{
    cChannelState next(c.size());
    for (auto &g : genes)
        for (auto &cis : g.modules)
            if (cis.active(c))
            {
                next.set(g.pub);
                // The gene is active, no use looking further
                break;
            }

    // Update the "return" value.
    c.swap(next);
}

// This is the inner loop, where we find the attractors. We should tune this.
void cNetwork::calc_attractors()
{
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
        cChannelStateVector &attr = attractors.back();

        // Copy the part the path that is the attractor, ignoring the transient
        for (size_t copy_at=attractor_begins_at; copy_at < path.size(); ++copy_at)
        {
            // TODO: This is where we could calculate the rates, and
            // constructed them at the same time.
            attr.push_back(path[copy_at]);
        }
    }
}

cFactory::cFactory(size_t seed)
    : random_engine(seed)
    , next_identifier(0)
    , pop_count(1)
    , gene_count(4)
    , cis_count(3)
{
}

void cFactory::init_environments()
{
    total_channels = cue_channels + reg_channels + out_channels;
    // Environments
    size_t env_count = 1 << cue_channels;
    for (size_t i = 0; i < env_count; ++i)
    {
        cChannelState c = cChannelState(total_channels, i);
        environments.push_back(c);
    }
}

cGeneMutator::cGeneMutator(cFactory *f, double r)
    : factory(f)
    , rate_per_gene(r)
    , r_gene(0, f->gene_count-1)
    , r_sub(f->sub_range.first, f->sub_range.second-1)
    , r_oper(0, f->operands.size()-1)
{
}

void cGeneMutator::mutate_gene(cGene &g)
{
    // TODO: This needs fixing !!
    cCisModule &m = g.modules[0];
    m.op = factory->operands[r_oper(factory->random_engine)];
    m.sub1 = r_sub(factory->random_engine);
}

void cGeneMutator::mutate_network(cNetwork_ptr &n, size_t mutations)
{
    // NOTE: This done INPLACE mutation. It should never be called on a network
    // that has already had its attractors calculated! 
    
    // Select the genes that should be mutated
    while (mutations > 0)
    {
        size_t i = r_gene(factory->random_engine);
        mutate_gene(n->genes[i]);
        --mutations;
    }
}

cNetwork_ptr cGeneMutator::copy_and_mutate_network(cNetwork_ptr &n, size_t mutations)
{
    cNetwork_ptr copy(new cNetwork(n->factory));
    // The copy constructor of vector does the hard work here...
    copy->genes = n->genes;
    copy->parent_identifier = n->identifier;

    // Now mutate it and calculate the attractors
    mutate_network(copy, mutations);
    copy->calc_attractors();
    return copy;
}

typedef std::poisson_distribution<> poisson_t;

void cGeneMutator::mutate_collection(cNetworkVector &networks)
{
    // How many mutations are we going to have? That depends on the total
    // number of sites that might mutate. Per network, this is rate_per_gene *
    // gene_count. We multiply this by the number of networks in the collection
    // to get the expected number of mutations. Then we generate a number using
    // a poisson distribution.
    double expected = rate_per_gene * networks.size() * factory->gene_count;
    poisson_t r_pop(expected);
    size_t mutations = r_pop(factory->random_engine);
    if (mutations == 0)
        return;

    // Now we need to assign these to individual networks. Only these networks 
    // will change. The rest remain constant.
    randint_t r_network(0, networks.size()-1);
    std::vector<size_t> mutes;
    for (size_t i=0; i < mutations; ++i)
        mutes.push_back(r_network(factory->random_engine));

    // Sort them so that repeats are next to one another.
    std::sort(mutes.begin(), mutes.end());

    // We now let the networks figure out how *exactly* they will mutate.  This
    // looks longwinded, but we want to handle cases where a single network is
    // mutated more than once in one step, rather than calling multiple times
    // on the same network.
    auto it = mutes.begin();
    size_t network = *it, count = 1;
    ++it;
    for (;;) 
     {
        // If we're at the end, or the next one is different
        if (it == mutes.end() || *it != network)
        {
            networks[network] = copy_and_mutate_network(networks[network], count);
            if (it == mutes.end())
                break;
            count = 1;
            network = *it;
        }
        else
        {
            // We have another that is the same
            count += 1;
            ++it;
        }
    }

}

