#include "core.hpp"
#include <cmath>
#include <limits>
// #include <stdexcept>
// #include <iostream>

using namespace pubsub2;

cPopulation::cPopulation(const cFactory_ptr &c, size_t size)
    : factory(c)
    , world(c->world)
{
    for (size_t i = 0; i < size; ++i)
        networks.push_back(factory->construct(true));
}

std::pair<double, double> cPopulation::worst_and_best() const
{
    double best = std::numeric_limits<double>::min();
    double worst = std::numeric_limits<double>::max();
    double f;
    for (size_t i = 0; i < networks.size(); ++i)
    {
        f = networks[i]->fitness;
        if (f > best)
            best = f;
        if (f < worst)
            worst = f;
    }
    return std::make_pair(worst, best);
}

void cPopulation::best_indexes(cIndexes &best) const
{
    best.clear();
    std::pair<double, double> wb = worst_and_best();
    for (size_t i = 0; i < networks.size(); ++i)
    {
        if (networks[i]->fitness == wb.second)
            best.push_back(i);
    }
}

void cPopulation::assess(const cTarget &target) const
{
    for (size_t i = 0; i < networks.size(); ++i)
    {
        const cNetwork &net = *(networks[i]);
        target.assess(net);
    }
}

bool cPopulation::select(const cSelectionModel &sm, size_t size)
{
    bool done = sm.select(networks, size, selected);
    if (!done)
        return false;

    cNetworkVector new_networks;
    for (auto i : selected)
        new_networks.push_back(networks[i]);

    // Replace everything -- this is fast
    networks.swap(new_networks);
    return true;
}

size_t cPopulation::mutate(double site_rate, int_t generation)
{
    // How many mutations are we going to have? That depends on the total
    // number of sites that might mutate. Per network, this is rate_per_gene *
    // gene_count. We multiply this by the number of networks in the collection
    // to get the expected number of mutations. Then we generate a number using
    // a poisson distribution.
    double expected = site_rate * factory->site_count(networks);
    std::poisson_distribution<> r_pop(expected);
    size_t m_count = r_pop(world->rand);

    // Clear this
    mutated.clear();

    // If we're not generating any m_count, let's just bail.
    if (m_count == 0)
        return 0;

    // Now we need to assign these to individual networks. Only these networks 
    // will change. The rest remain constant.
    randint_t r_network(0, networks.size()-1);
    cIndexes mutes;
    for (size_t i=0; i < m_count; ++i)
        mutes.push_back(r_network(world->rand));

    // Sort them so that any repeated networks are adjacent.
    std::sort(mutes.begin(), mutes.end());

    // We now let the networks figure out how *exactly* they will mutate.  This
    // looks longwinded, but we want to handle cases where a single network is
    // mutated more than once in one step, rather than calling multiple times
    // on the same network, so most of this is putting together multiple
    // m_count into a single call.
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
            networks[network_num] = factory->clone_and_mutate_network(
                networks[network_num], count, generation);

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

    // Number of Networks mutated (not number of mutations)
    return mutated.size();
}

