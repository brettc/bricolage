#include "selection.hpp"
#include <cmath>
#include <limits>
// #include <stdexcept>
#include <iostream>

using namespace bricolage;

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

// TODO: This should construct the fitnesses
// The next step should select based on them....
// How do we do the caching...?
void cPopulation::assess(const cBaseTarget &target)
{
    // Update the fitnesses based on the networks there.
    target.assess_networks(networks, fitnesses);
}

bool cPopulation::select(const cSelectionModel &sm, size_t size)
{
    bool done = sm.select(fitnesses, size, selected);

    if (!done)
        return false;

    cNetworkVector new_networks;
    for (auto i : selected)
        new_networks.push_back(networks[i]);

    // Replace everything -- this is fast
    networks.swap(new_networks);
    return true;
}


typedef std::map<size_t, std::pair<size_t, size_t> > mutate_index;

size_t cPopulation::mutate(double site_rate, double dup_rate, int_t generation)
{
    // How many mutations are we going to have? That depends on the total
    // number of sites that might mutate. Per network, this is rate_per_gene *
    // gene_count. We multiply this by the number of networks in the collection
    // to get the expected number of mutations. Then we generate a number using
    // a poisson distribution.
    double expected = site_rate * factory->site_count(networks);
    std::poisson_distribution<> r_pop(expected);
    size_t mut_count = r_pop(world->rand);

    size_t dup_count = 0;
    if (dup_rate > 0.0)
    {
        double dup_expected = dup_rate * world->reg_channels * networks.size();
        std::poisson_distribution<> r_pop2(dup_expected);
        dup_count = r_pop2(world->rand);
        // std::cout << "m expect " << expected << " dup expect " << dup_expected << std::endl;
    }

    // Clear our record of what has been mutated
    mutated.clear();

    // If we're not generating any, let's just bail.
    if (mut_count == 0 && dup_count == 0)
        return 0;

    // std::cout << "drate " << dup_rate << " Dup count " << dup_count << std::endl;

    // First, we figure the index of the networks that are going to mutate.  It
    // is possible for the same number to come up more than once, so we need to
    // keep the *number* of mutations per network.
    randint_t r_network(0, networks.size()-1);
    mutate_index mutes;

    // First the mutations
    for (size_t i=0; i < mut_count; ++i)
    {
        auto ret = mutes.emplace(r_network(world->rand), std::make_pair(1, 0));
        if (!ret.second)
        {
            mutate_index::value_type &v = (*ret.first);
            v.second.first += 1;
        }
    }
    // Now, duplications
    for (size_t i=0; i < dup_count; ++i)
    {
        auto ret = mutes.emplace(r_network(world->rand), std::make_pair(0, 1));
        if (!ret.second)
        {
            mutate_index::value_type &v = (*ret.first);
            v.second.second += 1;
        }
    }

    // Now we can go ahead and mutate the actual networks.
    for (const auto &v : mutes)
    {
        // Replace the networks.
        networks[v.first] = factory->clone_and_mutate_network(
                networks[v.first], 
                v.second.first, 
                v.second.second,
                generation);

        // Keep a record of what we mutated.
        mutated.push_back(v.first);
    }

    // Number of Networks mutated (not number of mutations)
    return mutated.size();
}

