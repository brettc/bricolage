#include "core.hpp"
#include <cmath>
#include <stdexcept>

using namespace pubsub2;

cTarget::cTarget(const cWorld_ptr &w)
    : world(w) 
    , identifier(w->get_next_target_ident())
{
    cRates temp;
    std::fill_n(std::back_inserter(temp), w->out_channels, 0.0);
    std::fill_n(std::back_inserter(optimal_rates), w->environments.size(), temp);
}

// TODO: per output weighting
// TODO: per env weighting
// TODO: Profiling -- make faster?
// TODO: Maybe we should be use multi_array??
double cTarget::assess(const cNetwork &net) const
{
    // We've already done it!
    if (net.target == identifier)
        return net.fitness;

    // TODO: check that factories are the same?
    size_t nsize = net.rates.size();
    size_t osize = world->out_channels;

    if (nsize != optimal_rates.size())
        throw std::runtime_error("optimal rates and networks rates differ in size");

    // We need to score each of the outputs individually
    // TODO: Maybe make this a member
    cRates scores;
    std::fill_n(std::back_inserter(scores), osize, 0.0);

    // The difference from the target reduces the score from a max of 1.0
    for (size_t i = 0; i < nsize; ++i)
        for (size_t j = 0; j < osize; ++j)
            scores[j] += ((1.0 - fabs(net.rates[i][j] - optimal_rates[i][j])) 
                          / double(nsize));

    // Now take the average
    double score = 0.0;
    for (size_t j = 0; j < osize; ++j)
    {
        scores[j] /= double(osize);
        score += scores[j];
    }

    // Must be greater 0 >= n <= 1.0
    net.fitness = score;

    // Record the target we've been used to assess. We can skip doing every
    // time then.
    net.target = identifier;
    return score;
}

cSelectionModel::cSelectionModel(cWorld_ptr &w)
    : world(w)
{
}

bool cSelectionModel::select(
    const cNetworkVector &networks, size_t number, cIndexes &selected) const
{
    selected.clear();

    std::vector<double> cum_scores;
    cIndexes indexes;
    double score, cum_score = 0.0;
    // First, score everyone
    for (size_t i = 0; i < networks.size(); ++i)
    {
        // NOTE: There must be a fitness assigned
        score = networks[i]->fitness;

        // Zero fitness has no chance
        if (score <= 0.0)
            continue;

        // TODO: Maybe add some scaling factor to the fitness here
        // ie. score * score, exp(score * y)
        cum_score += score;
        cum_scores.push_back(cum_score);
        indexes.push_back(i);
    }

    // Everyone was crap. 
    if (cum_scores.size() == 0)
        return false;

    // This is standard "roulette selection". The fitness of each individual,
    // once normalised across the population is proportional to their likelihood
    // of being selected for the next generation.
    std::uniform_real_distribution<double> wheel(0.0, cum_score);
    for (size_t i = 0; i < number; ++i)
    {
        double locator = wheel(world->rand);
        auto it = std::lower_bound(cum_scores.begin(), cum_scores.end(), locator);
        size_t found_index = it - cum_scores.begin();

        selected.push_back(indexes[found_index]);
    }

    return true;
}

cPopulation::cPopulation(const cConstructor_ptr &c, size_t size)
    : constructor(c)
    , world(c->world)
{
    for (size_t i = 0; i < size; ++i)
        networks.push_back(constructor->construct());
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
    double expected = site_rate * constructor->site_count(networks);
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
            networks[network_num] = constructor->clone_and_mutate_network(
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

