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

struct Mutations
{
    Mutations(size_t c, size_t t, size_t d) : cis(c), trans(t), dups(d) {}
    size_t cis, trans, dups;
};

typedef std::map<size_t, Mutations> mutate_index;

size_t cPopulation::mutate(double cis_rate,
                           double trans_rate,
                           double dup_rate, 
                           int_t generation)
{
    // How many mutations are we going to have? That depends on the total
    // number of sites that might mutate. Per network, this is rate_per_gene *
    // gene_count. We multiply this by the number of networks in the collection
    // to get the expected number of mutations. Then we generate a number using
    // a poisson distribution.
    size_t cis_count = 0;
    if (cis_rate > 0.0)
    {
        double cis_expected = cis_rate * factory->site_count(networks);
        std::poisson_distribution<> r_cis(cis_expected);
        cis_count = r_cis(world->rand);
    }

    size_t trans_count = 0;
    if (trans_rate > 0.0)
    {
        double trans_expected = trans_rate * world->reg_gene_count * networks.size();
        std::poisson_distribution<> r_trans(trans_expected);
        trans_count = r_trans(world->rand);
    }

    size_t dup_count = 0;
    if (dup_rate > 0.0)
    {
        double dup_expected = dup_rate * world->reg_gene_count * networks.size();
        std::poisson_distribution<> r_dup(dup_expected);
        dup_count = r_dup(world->rand);
    }

    // Clear our record of what has been mutated
    mutated.clear();

    // If we're not generating any, let's just bail.
    if (cis_count == 0 && trans_count == 0 && dup_count == 0)
        return 0;

    // First, we figure the index of the networks that are going to mutate.  It
    // is possible for the same number to come up more than once, and we want
    // to collate all of the mutations together.
    
    randint_t r_network(0, networks.size()-1);
    mutate_index mutes;

    // First cis mutations
    for (size_t i = 0; i < cis_count; ++i)
    {
        auto ret = mutes.emplace(r_network(world->rand), Mutations(1, 0, 0));
        if (!ret.second)
        {
            auto &v = (*ret.first);
            v.second.cis += 1;
        }
    }

    // Now, trans mutations
    for (size_t i = 0; i < trans_count; ++i)
    {
        auto ret = mutes.emplace(r_network(world->rand), Mutations(0, 1, 0));
        if (!ret.second)
        {
            auto &v = (*ret.first);
            v.second.trans += 1;
        }
    }

    // Now, duplications
    for (size_t i = 0; i < dup_count; ++i)
    {
        auto ret = mutes.emplace(r_network(world->rand), Mutations(0, 0, 1));
        if (!ret.second)
        {
            auto &v = (*ret.first);
            v.second.dups += 1;
        }
    }

    // Now we can go ahead and mutate the actual networks.
    for (const auto &v : mutes)
    {
        const size_t index = v.first;
        const Mutations &m = v.second;

        std::cout << index << ' ' << m.trans << std::endl;

        // Replace the networks with the mutated ones.
        networks[index] = factory->clone_and_mutate_network(
                networks[index], 
                m.cis,
                m.trans, 
                m.dups,
                generation);

        // Keep a record of what we mutated.
        mutated.push_back(index);
    }

    // Number of Networks mutated (not number of mutations)
    return mutated.size();
}

