#include "core.hpp"
#include <cmath>
#include <stdexcept>
#include <iostream>

using namespace pubsub2;

cTarget::cTarget(const cWorld_ptr &w, const std::string &n, int_t id)
    : world(w) 
    , name(n)
{
    if (id == -1)
        identifier = w->get_next_target_ident();
    else
        identifier = id;

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

void cTarget::assess_networks(cNetworkVector &networks) const
{
    for (auto net: networks)
        net->fitness = assess(*net);
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

