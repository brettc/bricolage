// vim: path=.,/usr/include/c++/4.2.1,/usr/include/c++/4.2.1/tr1
#include "pubsub2_c.h"
// #include <iostream>
// #include "Python.h"

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
}

typedef std::uniform_int_distribution<size_t> randint;
void cFactory::construct_random(cNetwork &network)
{
    randint r_sub(sub_range.first, sub_range.second-1);
    randint r_oper(0, operands.size()-1);
    
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
            attr.push_back(path[copy_at]);
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

