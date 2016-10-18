#include "targets.hpp"
#include <cmath>
#include <stdexcept>
#include <iostream>

using namespace bricolage;

cBaseTarget::cBaseTarget(const cWorld_ptr &w, 
                         const std::string &n, 
                         int_t id,
                         ScoringMethod m,
                         double s)
    : world(w)
    , name(n)
    , scoring_method(m)
    , strength(s)

{
    if (id == -1)
        identifier = w->get_next_target_ident();
    else
        identifier = id;

    // Default is equal weighting
    std::fill_n(std::back_inserter(weighting), w->out_channels, 1.0 /
                double(w->out_channels));

    cRates temp;
    std::fill_n(std::back_inserter(temp), w->out_channels, 0.0);
    std::fill_n(std::back_inserter(optimal_rates), w->environments.size(), temp);
}

void cBaseTarget::set_weighting(const cRates &wghts)
{
    if (wghts.size() != weighting.size())
        throw std::runtime_error("weighting size and output size differ");

    double sum = 0.0;
    for (auto w : wghts)
        sum += w;

    // Normalize the weights
    double normfactor = 1.0 / sum;
    for (size_t i = 0; i < wghts.size(); ++i)
        weighting[i] = wghts[i] * normfactor;

}

double cBaseTarget::score_rates(const cRatesVector &rates) const
{
    size_t nsize = rates.size();
    size_t osize = world->out_channels;

    if (nsize != optimal_rates.size())
        throw std::runtime_error("optimal rates and networks rates differ in size");

    // We need to score each of the outputs individually
    double normalize = double(nsize);
    double final_score = 0.0;

    switch (scoring_method)
    {
    case SCORE_EXPONENTIAL_VEC:
        {
        for (size_t i = 0; i < nsize; ++i)
        {
            double dist = 0.0;
            for (size_t j = 0; j < osize; ++j)
            {
                double diff = rates[i][j] - optimal_rates[i][j];
                diff *= weighting[j];
                dist += diff * diff;

            }
            if (dist != 0.0)
                dist = sqrt(dist);
            final_score += exp(-fabs(dist) / strength)
                / normalize;
        }
        }
        break;

    case SCORE_EXPONENTIAL:
        {
        cRates scores;
        std::fill_n(std::back_inserter(scores), osize, 0.0);

        for (size_t i = 0; i < nsize; ++i)
            for (size_t j = 0; j < osize; ++j)
                scores[j] += exp(-fabs(rates[i][j] - optimal_rates[i][j]) / strength)
                                  / normalize;
        // Now apply the weighting
        for (size_t j = 0; j < osize; ++j)
            final_score += scores[j] * weighting[j];
        }
        break;

    default:
    case SCORE_LINEAR:
        // The difference from the target reduces the score from a max of 1.0
        {
        cRates scores;
        std::fill_n(std::back_inserter(scores), osize, 0.0);

        for (size_t i = 0; i < nsize; ++i)
            for (size_t j = 0; j < osize; ++j)
                scores[j] += (1.0 - fabs(rates[i][j] - optimal_rates[i][j]))
                            / normalize;
        // Now apply the weighting
        for (size_t j = 0; j < osize; ++j)
            final_score += scores[j] * weighting[j];
        }
    }
    return final_score;
}

void cBaseTarget::assess_networks(const cNetworkVector &networks, 
                                  std::vector<double> &scores) const
{
    scores.clear();
    for (auto net: networks)
        scores.push_back(assess(*net));
}

cDefaultTarget::cDefaultTarget(const cWorld_ptr &w, 
                               const std::string &n, int_t id,
                               ScoringMethod m, double s)
    : cBaseTarget(w, n, id, m, s)
{
}

double cDefaultTarget::assess(const cNetwork &net) const
{
    // We've already done it!
    if (net.target == identifier)
        return net.fitness;

    double score = score_rates(net.rates);

    // Must be greater 0 >= n <= 1.0
    net.fitness = score;

    // Record the target we've been used to assess. We can skip doing every
    // time then.
    net.target = identifier;
    return score;
}

cNoisyTarget::cNoisyTarget(const cWorld_ptr &w, 
                           const std::string &n, int_t id,
                           ScoringMethod m, double s,
                           size_t perturb,
                           double pert_prop,
                           bool e_only)
    : cBaseTarget(w, n, id, m, s)
    , perturb_count(perturb)
    , perturb_prop(pert_prop)
    , env_only(e_only)
{
    if (pert_prop < 0.0 or pert_prop > 1.0)
        throw std::runtime_error("perturb prop must be >= 0.0 and <= 1.0");
}

double cNoisyTarget::assess(const cNetwork &net) const
{
    double base_score = score_rates(net.rates) * (1.0 - perturb_prop);
    double score = 0.0;

    // Note: we need to reassess every time here, as there is it is not
    // deterministic
    if (perturb_count > 0)
    {
        double mult = 1.0 / double(perturb_count);
        for (size_t i = 0; i < perturb_count; ++i)
        {
            calc_perturbation(net, env_only);
            score += mult * score_rates(rates_vec);
        }
        score *= perturb_prop;
    }

    // Add proportions
    score += base_score;

    // Must be greater 0 >= n <= 1.0
    net.fitness = score;

    // Record the target we've been used to assess. 
    net.target = identifier;

    return score;
}

void cNoisyTarget::calc_perturbation(const cNetwork &net, bool env_only) const
{
    // We should already have initial attractors and rates.  This is NOT
    // deterministic, unlike the calculation of attractors, as we randomize
    // the inputs. It also assumes we'd doing *pulses* of inputs rather than
    // constant inputs.
    rates_vec.clear();
    random_int_t r_env_state(random_int_range(0, 2, world));
    random_int_t r_reg_channel(random_int_range(world->cue_range.first,
                                                world->reg_range.second,
                                                world));

    // Go through each environment.
    for (auto &attr : net.attractors)
    {
        // For each environmental attractor:
        // 1. Randomly select one state of the attractor.
        // 2. Introduce a pulse of random states (maybe just in the env)
        // 3. Find new attractor
        //
        // 1. Get starting point
        size_t state_i = 0;
        if (attr.size() > 1)
        {
            random_int_t r_state(random_int_range(0, attr.size(), world));
            state_i = r_state();
        }
        cChannels start_state = attr[state_i];

        // 2. Random state
        if (env_only)
        {
            for (size_t k = world->cue_range.first; k < world->cue_range.second; ++k)
                if (r_env_state())
                    start_state.unchecked_set(k);
        } else {
            start_state.unchecked_flip(r_reg_channel());
        }

        // Add a new attractor for this environment.
        rates_vec.emplace_back();
        auto &this_rate = rates_vec.back();

        net.get_rates(start_state, this_rate, true);
    }
}

cTransNoisyTarget::cTransNoisyTarget(const cWorld_ptr &w, 
                           const std::string &n, int_t id,
                           ScoringMethod m, double s)
    : cBaseTarget(w, n, id, m, s)
{
}

double cTransNoisyTarget::assess(const cNetwork &net) const
{
    // // We've already done it!
    if (net.target == identifier)
        return net.fitness;
    
    cRatesVector rates_vec;
    double score = 0.0;
    // Go through each environment.
    for (size_t noisy_channel = world->reg_range.first;
         noisy_channel < world->reg_range.second; ++noisy_channel)
    {
        rates_vec.clear();
        
        for (const auto &env : world->environments)
        {
            cChannels start_state = env;

            // Randomly flip a regulatory gene channel
            start_state.unchecked_set(noisy_channel);

            rates_vec.emplace_back();
            auto &this_rate = rates_vec.back();
            net.get_rates(start_state, this_rate, true);
        }

        score += score_rates(rates_vec);
    }

    score /= world->reg_channels;

    // Must be greater 0 >= n <= 1.0
    net.fitness = score;

    // Record the target we've been used to assess. We can skip doing every
    // time then.
    net.target = identifier;
    return score;
}

cMultiTarget::cMultiTarget(const cWorld_ptr &w, 
                           const std::string &n, int_t id,
                           ScoringMethod m, double s)
    : cBaseTarget(w, n, id, m, s)
{
}

void recurse_attractors(cChannels initial, cAttractor &inputs, 
                        cNetwork &net, cRates &retval)
{
    // Killed this for now
    //
    // initial.unchecked_union(inputs[0]);
    //
    // cAttractor attractor;
    // cAttractor transient;
    // cRates rates;
    //
    // if (inputs.size() == 1)
    // {
    //     // We've bottomed out. There are no more inputs to come. So this is our
    //     // final rates. We return it, and if there were multiple paths to here
    //     // then it will be averaged.
    //     net.stabilise(initial, false, attractor, transient, retval);
    //     return;
    // }
    // else
    // {
    //     net.stabilise(initial, false, attractor, transient, retval);
    // }
    //
    // // We need to send this inputs
    // cAttractor more_filter;
    // std::copy(inputs.begin() + 1, inputs.end(), std::back_inserter(more_filter));
    //
    // for (auto &state : attractor)
    // {
    //     cRates attr_rates;
    //     recurse_attractors(state, more_filter, net, attr_rates);
    //     // Average over the attractors
    // }

}

double cMultiTarget::assess(const cNetwork &net) const
{
    // Kill this for now
    double score = 0.0;
    //
    // // Only score the rates at the end
    // rates_vec.clear();
    //
    // cAttractorSet attractors;
    // cAttractorSet transients;
    // cRatesVector rates;
    //
    // // For each environment
    // for (auto &env : world->environments)
    // {
    //     // Got through and filter out only what we need
    //     cChannels env_copy(env);
    //     for (auto &pulse : pulses)
    //     {
    //         env_copy.unchecked_intersection(pulse);
    //
    //         attractors.emplace_back();
    //         cAttractor &this_attr = attractors.back();
    //         rates.emplace_back();
    //         cRates &this_rate = rates.back();
    //         transients.emplace_back();
    //         cAttractor &this_trans = transients.back();
    //
    //         net.stabilise(env_copy, false, this_attr, this_trans, this_rate);
    //
    //         // Now we need to use the attractor that comes from here. Merge it
    //         // with the pulse, and then see where we go next.
    //         env_copy = this_attr[0];
    //     }
    // }

    return score;
}


cSelectionModel::cSelectionModel(cWorld_ptr &w)
    : world(w)
{
}

bool cSelectionModel::select(const cRates &scores,
                             size_t number, 
                             cIndexes &selected) const
{
    selected.clear();

    cRates cum_scores;
    cIndexes indexes;
    double score, cum_score = 0.0;
    // First, score everyone
    for (size_t i = 0; i < scores.size(); ++i)
    {
        // NOTE: There must be a fitness assigned
        score = scores[i];

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

