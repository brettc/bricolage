#include "analysis.hpp"
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <map>
#include <exception>
#include <sstream>
#include <array>


using namespace bricolage;


double calc_entropy(double normalizer, const std::vector<double> &probs)
{
    double entropy = 0.0;
    for (auto pr : probs)
    {
        // We need to normalize this; it's not done above.
        double cond_pr = normalizer * pr;
        if (not_zeroish(cond_pr))
        {
            double result = cond_pr * -log2(cond_pr);
            entropy += result;
        }
    }
    return entropy;
}


cNetworkAnalysis::cNetworkAnalysis(cNetwork_ptr &n)
    : original(n)
    , modified(original->clone())
    , potential_bindings(0)
    , active_bindings(0)
{
}

void cNetworkAnalysis::make_edges(cEdgeList &edges)
{
    // Generate three types of edges
    // Gene -> Channel
    // Channel -> Module
    // Module -> Gene
    edges.clear();
    potential_bindings = 0;
    for (size_t i=0; i < original->gene_count(); ++i)
    {
        cGene *g = original->get_gene(i);
        Node_t gnode = std::make_pair(NT_GENE, i);
        Node_t cnode = std::make_pair(NT_CHANNEL, g->pub);
        edges.insert(Edge_t(gnode, cnode));
        for (size_t j=0; j < g->module_count(); ++j)
        {
            cCisModule *m = g->get_module(j);
            Node_t mnode = std::make_pair(NT_MODULE, make_module_node_id(i, j));
            edges.insert(Edge_t(mnode, gnode));
            for (size_t k=0; k < m->site_count(); ++k)
            {
                signal_t current = m->get_site(k);
                // 0 is NEVER active. Ever.
                if (current == 0)
                    continue;
                Node_t cnode = std::make_pair(NT_CHANNEL, current);
                edges.insert(Edge_t(cnode, mnode));
                potential_bindings++;
            }
        }
    }
}


void cNetworkAnalysis::make_active_edges(cEdgeList &edges)
{
    // Incrementally add knockouts that have no effect on the rates. Start big,
    // at the gene level and then move down.
    //
    // This isn't perfect, but it will do for now.  TODO: Fix this so it handle
    // double knockouts etc.
    edges.clear();
    active_bindings = 0;
    for (size_t i=0; i < modified->gene_count(); ++i)
    {
        cGene *g = modified->get_gene(i);

        // Can we knock this gene out?
        g->intervene = INTERVENE_OFF;
        modified->calc_attractors_with_intervention();

        // If this makes no diff, then no edges attached to this Gene make any
        // difference. Let's skip to the next one.
        if (modified->rates == original->rates)
            continue;

        // Ok, it does something ...
        g->intervene = INTERVENE_NONE;
        Node_t gnode = std::make_pair(NT_GENE, i);
        Node_t cnode = std::make_pair(NT_CHANNEL, g->pub);
        edges.insert(Edge_t(gnode, cnode));

        for (size_t j=0; j < g->module_count(); ++j)
        {
            cCisModule *m = g->get_module(j);
            m->intervene = INTERVENE_OFF;
            modified->calc_attractors_with_intervention();
            if (modified->rates == original->rates)
                continue;

            // Reset and add edge
            m->intervene = INTERVENE_NONE;
            Node_t mnode = std::make_pair(NT_MODULE, make_module_node_id(i, j));
            edges.insert(Edge_t(mnode, gnode));

            for (size_t k=0; k < m->site_count(); ++k)
            {
                // Knock this out by setting it to ZERO channel. Nothing ever
                // publishes to the ZERO channel!
                signal_t original_channel = m->set_site(k, 0);

                // It was already zero?
                if (original_channel == 0)
                    continue;

                // Otherwise, we need to test to see if it changed anything
                modified->calc_attractors_with_intervention();

                // Did it change?
                if (modified->rates != original->rates)
                {
                    // We need to reset it.
                    m->set_site(k, original_channel);

                    // We don't both including any signals that come from the
                    // permanently ON channel.
                    if (original_channel != 1)
                    {
                        // And include it in the graph.
                        Node_t cnode = std::make_pair(NT_CHANNEL, original_channel);
                        edges.insert(Edge_t(cnode, mnode));
                    }
                    active_bindings++;

                }
            }
        }
    }
}

// As above, but don't bother constructing edges
size_t cNetworkAnalysis::calc_active_bindings()
{
    active_bindings = 0;
    for (size_t i=0; i < modified->gene_count(); ++i)
    {
        cGene *g = modified->get_gene(i);

        // Can we knock this gene out?
        g->intervene = INTERVENE_OFF;
        modified->calc_attractors_with_intervention();

        // If this makes no diff, then no edges attached to this Gene make any
        // difference. Let's skip to the next one.
        if (modified->rates == original->rates)
            continue;

        // Ok, it does something ...
        g->intervene = INTERVENE_NONE;
        for (size_t j=0; j < g->module_count(); ++j)
        {
            cCisModule *m = g->get_module(j);
            m->intervene = INTERVENE_OFF;
            modified->calc_attractors_with_intervention();
            if (modified->rates == original->rates)
                continue;

            // Reset and add edge
            m->intervene = INTERVENE_NONE;
            for (size_t k=0; k < m->site_count(); ++k)
            {
                // Knock this out by setting it to ZERO channel. Nothing ever
                // publishes to the ZERO channel!
                signal_t original_channel = m->set_site(k, 0);

                // It was already zero?
                if (original_channel == 0)
                    continue;

                // Otherwise, we need to test to see if it changed anything
                modified->calc_attractors_with_intervention();

                // Did it change?
                if (modified->rates != original->rates)
                {
                    // We need to reset it.
                    m->set_site(k, original_channel);
                    active_bindings++;
                }
            }
        }
    }

    return active_bindings;
}


//-----------------------------------------------------------------------------
cInformation::cInformation(const cJointProbabilities &jp)
    : world(jp.world)
    , _array(boost::extents
             [jp._array.shape()[0]]
             [jp._array.shape()[1]]
             [jp._array.shape()[2]]
             )
{
    jp.calc_information(*this);
}

// This is used for information about output channels
cInformation::cInformation(const cWorld_ptr &w, size_t network_size,
                           size_t per_channel_size)
    : world(w)
    , _array(boost::extents
             [network_size]
             [w->reg_channels]
             [per_channel_size])
{
    // Create an empty array
    std::fill(_array.origin(),
              _array.origin() + _array.num_elements(), 0.0);
}

cJointProbabilities::cJointProbabilities(
    const cWorld_ptr &w, size_t network_size, size_t per_network,
    size_t per_channel)
    : world(w)
    , _array(boost::extents
             [network_size]
             [w->reg_channels]
             [per_network]
             [2] // This always Binary, ON/OFF
             [per_channel])
{
}

void cJointProbabilities::calc_information(cInformation &info) const
{
    auto networks_size = _array.shape()[0];
    auto channels_size = _array.shape()[1];

    // Multiple output channels, or simple one input environment
    auto per_channel_size = _array.shape()[2];

    // This is ON/OFF
    auto row_size = _array.shape()[3];

    // This is the number of categories (environments or rates)
    auto col_size = _array.shape()[4];

    typedef joint_array_type::index index;

    // Used to hold the marginals (column/row sums).
    cRates rows(row_size);
    cRates cols(col_size);

    for (index i = 0; i < networks_size; ++i)
    {
        for (index j = 0; j < channels_size; ++j)
        {
            for (index k = 0; k < per_channel_size; ++k)
            {
                // Reset marginals
                for (index ri = 0; ri < row_size; ++ri)
                    rows[ri] = 0.0;
                for (index ci = 0; ci < col_size; ++ci)
                    cols[ci] = 0.0;

                // Sum the marginals
                for (index ri = 0; ri < row_size; ++ri)
                    for (index ci = 0; ci < col_size; ++ci)
                    {
                        double val = _array[i][j][k][ri][ci];
                        rows[ri] += val;
                        cols[ci] += val;
                    }

                // Calculate the info
                double I = 0.0;
                for (index ri = 0; ri < row_size; ++ri)
                    for (index ci = 0; ci < col_size; ++ci)
                    {
                        double val = _array[i][j][k][ri][ci];
                        double denom = rows[ri] * cols[ci];
                        if (not_zeroish(val) && not_zeroish(denom))
                            I += val * log2(val / denom);
                    }
                info._array[i][j][k] = I;
            }
        }
    }
}

cBaseCausalAnalyzer::cBaseCausalAnalyzer(cWorld_ptr &w)
    : world(w)
    , intervention_probs(w->reg_channels, 0.5)
{
}

size_t cBaseCausalAnalyzer::max_category = 16;

size_t cBaseCausalAnalyzer::get_max_category_size()
{
    return max_category;
}

void cBaseCausalAnalyzer::set_max_category_size(size_t m)
{
    max_category = m;
}

// Calculate the probability of any particular signal being on when the
// attractor network is not being intervened upon
void cBaseCausalAnalyzer::_calc_natural(cNetwork &net)
{
    for (size_t i = 0; i < world->reg_channels; ++i)
        intervention_probs[i] = 0.0;

    size_t reg_base = world->reg_range.first;

    // Go through attractors
    for (auto &cv : net.attractors)
    {
        // Probability of each state depends on environments, and number states
        // in this attractor. Everything is equiprobable (for now)
        double p_state = 1.0 / net.attractors.size() / cv.size();
        // Each state in the attractor
        for (auto &cs : cv)
        {
            // See what is on or off
            for (size_t i = 0; i < world->reg_channels; ++i)
            {
                if (cs.unchecked_test(reg_base + i))
                    intervention_probs[i] += p_state;
            }
        }
    }
}

//-----------------------------------------------------------------------------
// Keep a map of the unique rates that are output and assign them to persistent
// categories within a network. We'll use these categories to calculate the
// information.
size_t cRateCategorizer::get_category(double rate, double prob)
{
    auto result = rate_categories.insert(
        std::make_pair(rate, next_category));
    if (result.second)
    {
        // Successful insert. We have a new category. Update the category.
        if (++next_category > cBaseCausalAnalyzer::max_category)
            throw std::out_of_range("Maximum categories reached!!");

        // Set the probabilities on this category
        category_probabilities.push_back(prob);
    }
    else
    {
        // Update the probabilities
        category_probabilities[(*result.first).second] += prob;
    }

    // In either case, return the value in the pair the insert points to.
    // (It will either be the new pair, or the one found).
    return (*result.first).second;
}

void cRateCategorizer::clear()
{
    next_category = 0;
    rate_categories.clear();
    category_probabilities.clear();
}

cCausalFlowAnalyzer::cCausalFlowAnalyzer(cWorld_ptr &w)
    : cBaseCausalAnalyzer(w)
{
}

cJointProbabilities *cCausalFlowAnalyzer::analyse_network(cNetwork &net)
{
    cJointProbabilities *joint =
        new cJointProbabilities(world,
                                1,
                                world->out_channels,
                                cBaseCausalAnalyzer::max_category);

    _analyse(net, joint->_array[0]);

    return joint;
}

cJointProbabilities *cCausalFlowAnalyzer::analyse_collection(
    const cNetworkVector &networks)
{
    cJointProbabilities *joint =
        new cJointProbabilities(world,
                                networks.size(),
                                world->out_channels,
                                cBaseCausalAnalyzer::max_category);

    for (size_t i = 0; i < networks.size(); ++i)
        _analyse(*networks[i], joint->_array[i]);
    return joint;
}


void cCausalFlowAnalyzer::_analyse(cNetwork &net, joint_array_type::reference sub)
{
    _calc_natural(net);

    cRateCategorizer categorizer;
    double p_env = 1.0 / net.rates.size();

    // For each channel
    for (size_t i = 0; i < world->reg_channels; ++i)
    {
        cGene *gene = net.get_gene(i);
        gene->intervene = INTERVENE_OFF;
        net.calc_attractors_with_intervention();
        double p_gene_off = p_env * (1.0 - intervention_probs[i]);

        for (size_t j = 0; j < net.rates.size(); ++j)
            for (size_t k = 0; k < world->out_channels; ++k)
            {
                size_t cat = categorizer.get_category(net.rates[j][k], p_gene_off);
                sub[i][k][0][cat] += p_gene_off;
            }

        gene->intervene = INTERVENE_ON;
        net.calc_attractors_with_intervention();
        double p_gene_on = p_env * intervention_probs[i];

        // for each environment (there are rates for each)
        for (size_t j = 0; j < net.rates.size(); ++j)
            for (size_t k = 0; k < world->out_channels; ++k)
            {
                size_t cat = categorizer.get_category(net.rates[j][k], p_gene_on);
                sub[i][k][1][cat] += p_gene_on;
            }

        // Reset
        gene->intervene = INTERVENE_NONE;
    }

    // Finally, reset everything, without using interventions.
    net.calc_attractors();
}



// --------------------------------------------------------------------------
cAverageControlAnalyzer::cAverageControlAnalyzer(cWorld_ptr &world)
    : cBaseCausalAnalyzer(world)
    , categorizers(boost::extents[world->reg_channels][world->out_channels])
    , joint_over_envs(world, world->environments.size(),
                      world->out_channels,
                      cBaseCausalAnalyzer::max_category)
{
}

// Note you need to delete the return values from these!
cInformation *cAverageControlAnalyzer::analyse_network(cNetwork &net)
{
    // Entropies and mutual info
    cInformation *info = new cInformation(world, 1, world->out_channels * 2);
    _analyse(net, info->_array[0]);
    return info;
}

cInformation *cAverageControlAnalyzer::analyse_collection(
    const cNetworkVector &networks)
{
    // Entropies and mutual info
    cInformation *info = new cInformation(world, networks.size(), world->out_channels * 2);
    for (size_t i = 0; i < networks.size(); ++i)
        _analyse(*networks[i], info->_array[i]);
    return info;
}

void cAverageControlAnalyzer::_clear()
{
    std::fill(joint_over_envs._array.origin(),
              joint_over_envs._array.origin() + joint_over_envs._array.num_elements(), 0.0);
    for (size_t i = 0; i < world->reg_channels; ++i)
        for (size_t j = 0; j < world->out_channels; ++j)
            categorizers[i][j].clear();
}

void cAverageControlAnalyzer::_analyse(
    cNetwork &net,
    info_array_type::reference sub
)
{
    // We need several probability distributions. Note that we use this
    // differently than above, as each of the major axis entries represents an
    // environment for a particular network, rather than the summed dist_n for
    // a network in population.
    auto &world = net.factory->world;

    _calc_natural(net);
    _clear();

    // Note the reversed order here (j, i, k). I've exchanged the loops because
    // it is very expensive to recalculate the attractors. So we do it as
    // little as possible.
    for (size_t j = 0; j < world->reg_channels; ++j)
    {
        cGene *gene = net.get_gene(j);
        double p_gene_on = intervention_probs[j];
        double p_gene_off = 1.0 - p_gene_on;

        // Force gene OFF
        gene->intervene = INTERVENE_OFF;
        net.calc_attractors_with_intervention();
        // for each environment and each output channel
        for (size_t i = 0; i < net.rates.size(); ++i)
            for (size_t k = 0; k < world->out_channels; ++k)
            {
                size_t cat = categorizers[j][k].get_category(net.rates[i][k], p_gene_off);
                joint_over_envs._array[i][j][k][0][cat] += p_gene_off;
            }

        // Force gene ON
        gene->intervene = INTERVENE_ON;
        net.calc_attractors_with_intervention();
        // for each environment and each output channel
        for (size_t i = 0; i < net.rates.size(); ++i)
            for (size_t k = 0; k < world->out_channels; ++k)
            {
                size_t cat = categorizers[j][k].get_category(net.rates[i][k], p_gene_on);
                joint_over_envs._array[i][j][k][1][cat] += p_gene_on;
            }

        // Reset
        gene->intervene = INTERVENE_NONE;
    }

    // Recalculate one more time without using interventions.
    net.calc_attractors();

    // Now calculate the information in each of these environments. Simply
    // using the constructor does this.
    cInformation info(joint_over_envs);

    // Now take the env weighted average of all of these, and summarize them in
    // the object that was passed.
    // Currently all environments have the same probability.
    size_t output_size = world->out_channels;
    double p_env = 1.0 / net.rates.size();
    for (size_t i = 0; i < net.rates.size(); ++i)
        for (size_t j = 0; j < world->reg_channels; ++j)
            for (size_t k = 0; k < world->out_channels; ++k)
            {
                sub[j][k] += info._array[i][j][k] * p_env;
                cRateCategorizer &categ = categorizers[j][k];
                sub[j][k + output_size] = calc_entropy(p_env, categ.category_probabilities);
            }

}


// --------------------------------------------------------------------------
// The categories define how each of the possible enviroments is mapped to a
// category.
cMutualInfoAnalyzer::cMutualInfoAnalyzer(cWorld_ptr &w, const cIndexes &cats)
    : world(w)
    , categories(cats)
    , max_category(1 + *std::max_element(categories.begin(), categories.end()))
{
}

cJointProbabilities *cMutualInfoAnalyzer::analyse_network(
    cNetwork &net)
{
    cJointProbabilities *joint =
        new cJointProbabilities(world, 1, 1, max_category);

    _analyse(net, joint->_array[0]);

    return joint;
}

cJointProbabilities *cMutualInfoAnalyzer::analyse_collection(
    const cNetworkVector &networks)
{
    cJointProbabilities *joint =
        new cJointProbabilities(world, networks.size(), 1, max_category);

    for (size_t i = 0; i < networks.size(); ++i)
        _analyse(*networks[i], joint->_array[i]);
    return joint;
}

void cMutualInfoAnalyzer::_analyse(
    cNetwork &net, joint_array_type::reference sub)
{
    size_t reg_base = world->reg_range.first;
    double normal_p_event = 1.0 / double(world->environments.size());

    // For each environment
    for (size_t ei = 0; ei < net.attractors.size(); ++ei)
    {
        auto attrs = net.attractors[ei];

        // Each environment is given a category. We then analyse the
        // information it has about that category.
        size_t cat = categories[ei];
        double p_event = normal_p_event / double(attrs.size());
        for (auto &att : attrs)
        {
            // Each attractor state...
            for (size_t ci = 0; ci < world->reg_channels; ++ci)
            {
                // Is this gene on or off?
                size_t on_off = att.unchecked_test(reg_base + ci);
                // What information does this carry about the particular
                // category assigned to this environment
                sub[ci][0][on_off][cat] += p_event;
            }
        }
    }
}

//-----------------------------------------------------------------------------
cOutputCategorizer::cOutputCategorizer(const cRatesVector &tr, size_t env_size)
    : next_category(0)
    , target_rates(tr)
{
    std::fill_n(std::back_inserter(targets_hit_in_env), env_size, 0.0);
}

size_t cOutputCategorizer::get_category(const cRates &rates, double prob, size_t env)
{
    // Record the probability that we actually hit one of the rates we wanted.
    bool match;
    for (auto &tr : target_rates)
    {
        match = true;
        for (size_t i = 0; i < tr.size(); ++i)
        {
            if (!is_close(tr[i], rates[i]))
            {
                match = false;
                break;
            }
        }
        if (match)
        {
            targets_hit_in_env[env] += prob;
            break;
        }
    }

    auto result = rate_categories.insert(
        std::make_pair(rates, next_category));
    if (result.second)
    {
        // Successful insert. We have a new category. Update the category.
        if (++next_category > cBaseCausalAnalyzer::max_category)
            throw std::out_of_range("Maximum categories reached!!");

        // Set the probabilities on this category
        category_probabilities.push_back(prob);
    }
    else
    {
        // Update the probabilities
        category_probabilities[(*result.first).second] += prob;
    }

    // In either case, return the value in the pair the insert points to.
    // (It will either be the new pair, or the one found).
    return (*result.first).second;
}

void cOutputCategorizer::clear()
{
    next_category = 0;
    rate_categories.clear();
    category_probabilities.clear();
    std::fill(targets_hit_in_env.begin(), targets_hit_in_env.end(), 0.0);
}

// --------------------------------------------------------------------------
cOutputControlAnalyzer::cOutputControlAnalyzer(cWorld_ptr &world, const cRatesVector &tr)
    : cBaseCausalAnalyzer(world)
    , joint_over_envs(world, world->environments.size(),
                      1, // only 1 per channel in this case.
                      cBaseCausalAnalyzer::max_category)
    , target_rates(tr)
{
    for (size_t i = 0; i < world->reg_channels; ++i)
        categorizers.push_back(cOutputCategorizer(target_rates, world->environments.size()));
}

// Note you need to delete the return values from these!
cInformation *cOutputControlAnalyzer::analyse_network(cNetwork &net)
{
    // One network. 2 information measures. 0 is causal power. 1 is output.
    // entropy. 2 is weighted info.
    cInformation *info = new cInformation(world, 1, 3);
    _analyse(net, info->_array[0]);
    return info;
}

cInformation *cOutputControlAnalyzer::analyse_collection(
    const cNetworkVector &networks)
{
    // Many networks. 2 information measures. 0 is causal power. 1 is output
    // entropy. 2 is weighted info.
    cInformation *info = new cInformation(world, networks.size(), 3);
    for (size_t i = 0; i < networks.size(); ++i)
        _analyse(*networks[i], info->_array[i]);
    return info;
}

void cOutputControlAnalyzer::_clear()
{
    std::fill(joint_over_envs._array.origin(),
              joint_over_envs._array.origin() + joint_over_envs._array.num_elements(), 0.0);
    for (auto &cat : categorizers)
        cat.clear();
}

void cOutputControlAnalyzer::_analyse(
    cNetwork &net,
    info_array_type::reference sub
)
{
    auto &world = net.factory->world;

    _clear();
    _calc_natural(net);

    // Note the reversed order here (j, i). I've exchanged the loops because
    // it is very expensive to recalculate the attractors. So we do it as
    // little as possible.
    for (size_t j = 0; j < world->reg_channels; ++j)
    {
        cGene *gene = net.get_gene(j);
        double p_gene_on = intervention_probs[j];
        double p_gene_off = 1.0 - p_gene_on;

        // Force gene OFF
        gene->intervene = INTERVENE_OFF;
        net.calc_attractors_with_intervention();
        // for each environment and each output channel
        if (p_gene_off != 0.0)
            for (size_t i = 0; i < net.rates.size(); ++i)
            {
                size_t cat = categorizers[j].get_category(net.rates[i], p_gene_off, i);
                joint_over_envs._array[i][j][0][0][cat] += p_gene_off;
            }

        // Force gene ON
        gene->intervene = INTERVENE_ON;
        net.calc_attractors_with_intervention();
        // for each environment and each output channel
        if (p_gene_on != 0.0)
            for (size_t i = 0; i < net.rates.size(); ++i)
            {
                size_t cat = categorizers[j].get_category(net.rates[i], p_gene_on, i);
                joint_over_envs._array[i][j][0][1][cat] += p_gene_on;
            }

        // Reset
        gene->intervene = INTERVENE_NONE;
    }

    // Recalculate one more time without using interventions. This just resets everything.
    net.calc_attractors();

    // Now calculate the information in each of these environments.
    cInformation info(world, world->environments.size(), 1);

    // Note that only the first entry will be calculated here. We'll put
    // entropy into the other.
    joint_over_envs.calc_information(info);

    // Now take the env weighted average of all of these, and summarize them in
    // the object that was passed.
    // Currently all environments have the same probability.
    double p_env = 1.0 / net.rates.size();
    for (size_t i = 0; i < net.rates.size(); ++i)
        for (size_t j = 0; j < world->reg_channels; ++j)
        {
            // Summarize info
            sub[j][0] += info._array[i][j][0] * p_env;

            // Now calc the entropy of the output from manipulating this channel
            cOutputCategorizer &categ = categorizers[j];
            sub[j][1] = calc_entropy(p_env, categ.category_probabilities);

            // Summarize info, but SCALE it, according to the targets we want
            sub[j][2] += info._array[i][j][0] * p_env * categ.targets_hit_in_env[i];
        }

}

// --------------------------------------------------------------------------
cRelevantControlAnalyzer::cRelevantControlAnalyzer(cWorld_ptr &world, const cRatesVector &tr,
                                                   bool use_natural_)
    : cBaseCausalAnalyzer(world)
    , target_rates(tr)
    , categories(boost::extents[world->environments.size()][world->reg_channels])
    , info(boost::extents[world->environments.size()][world->reg_channels])
    , use_natural(use_natural_)
{
}

cInformation *cRelevantControlAnalyzer::analyse_network(cNetwork &net)
{
    cInformation *info = new cInformation(world, 1, 1);
    _analyse(net, info->_array[0]);
    return info;
}

cInformation *cRelevantControlAnalyzer::analyse_collection(
    const cNetworkVector &networks)
{
    cInformation *info = new cInformation(world, networks.size(), 1);
    for (size_t i = 0; i < networks.size(); ++i)
        _analyse(*networks[i], info->_array[i]);
    return info;
}

void cRelevantControlAnalyzer::_clear()
{
    std::fill(categories.origin(), categories.origin() + categories.num_elements(), -1);
    std::fill(info.origin(), info.origin() + info.num_elements(), 0.0);
}

int cRelevantControlAnalyzer::categorize(const cRates rates)
{
    bool match;
    int i = 0;
    for (; i < target_rates.size(); ++i)
    {
        auto &tr = target_rates[i];
        match = true;
        for (size_t j = 0; j < tr.size(); ++j)
        {
            // Negative rates match EVERYTHING (we don't care about these outputs)
            if (tr[j] < 0.0)
                continue;

            // If any of the rates don't match, then it is NOT a match. So we
            // can quit early.
            if (!is_close(tr[j], rates[j]))
            {
                match = false;
                break;
            }
        }
        if (match)
            break;
    }
    if (match)
        return i;
    return -1;
}

void cRelevantControlAnalyzer::_analyse(cNetwork &net, info_array_type::reference sub)
{
    auto &world = net.factory->world;

    _clear();

    if (use_natural)
        _calc_natural(net);

    // Note the reversed order here (j, i). I've exchanged the loops because
    // it is very expensive to recalculate the attractors. So we do it as
    // little as possible.
    for (size_t j = 0; j < world->reg_channels; ++j)
    {
        cGene *gene = net.get_gene(j);
        double p_gene_on = intervention_probs[j];
        double p_gene_off = 1.0 - p_gene_on;

        // Force gene OFF
        gene->intervene = INTERVENE_OFF;
        net.calc_attractors_with_intervention();
        // for each environment and each output channel
        if (not_zeroish(p_gene_off))
            for (size_t i = 0; i < net.rates.size(); ++i)
            {
                int cat = categorize(net.rates[i]);
                categories[i][j] = cat;
                if (cat >= 0)
                    info[i][j] += p_gene_off * log2(1.0 / p_gene_off);
            }

        // Force gene ON
        gene->intervene = INTERVENE_ON;
        net.calc_attractors_with_intervention();
        // for each environment and each output channel
        if (not_zeroish(p_gene_on))
            for (size_t i = 0; i < net.rates.size(); ++i)
            {
                int cat = categorize(net.rates[i]);
                if (cat >= 0)
                {
                    if (cat == categories[i][j])
                        info[i][j] = 0.0;
                    else
                        info[i][j] += p_gene_on * log2(1.0 / p_gene_on);
                }
            }

        // Reset
        gene->intervene = INTERVENE_NONE;
    }

    // Recalculate one more time without using interventions. This just resets everything.
    net.calc_attractors();

    // Now take the env weighted average of all of these, and summarize them in
    // the object that was passed.
    double p_env = 1.0 / net.rates.size();
    for (size_t i = 0; i < net.rates.size(); ++i)
        for (size_t j = 0; j < world->reg_channels; ++j)
            sub[j][0] += info[i][j] * p_env;
}


//---------------------------------------------------------------------------
cMIAnalyzer::cMIAnalyzer(cWorld_ptr &w, const cIndexes &cats)
    : world(w)
    , categories(cats)
    // , p_cat({{0.0, 0.0}})
{
    // double p_cat = 1.0 / (double)categories.size();
    for (auto &c : categories)
    {
        if (c == 0)
            ;
            // p_cat[0] += p_cat;
        else if (c == 1)
            ;
            // p_cat[1] += p_cat;
        else
            throw std::out_of_range("categories must be 0 or 1");
    }
}

cInformation *cMIAnalyzer::analyse_network(
    cNetwork &net)
{
    cJointProbabilities joint(world, 1, 1, 2);
    _analyse(net, joint._array[0]);
    return new cInformation(joint);
}

cInformation *cMIAnalyzer::analyse_collection(
    const cNetworkVector &networks)
{
    cJointProbabilities joint(world, networks.size(), 1, 2);

    for (size_t i = 0; i < networks.size(); ++i)
        _analyse(*networks[i], joint._array[i]);
    return new cInformation(joint);
}

void cMIAnalyzer::_analyse(
    cNetwork &net, joint_array_type::reference sub)
{
    size_t reg_base = world->reg_range.first;
    double normal_p_event = 1.0 / double(world->environments.size());

    // For each environment
    for (size_t ei = 0; ei < net.attractors.size(); ++ei)
    {
        auto attrs = net.attractors[ei];

        // We categorise the environments into two cats, depending on
        // whether they are the focal category or not.
        size_t cat = categories[ei];
        double p_event = normal_p_event / double(attrs.size());
        for (auto &att : attrs)
        {
            // Each attractor state...
            for (size_t ci = 0; ci < world->reg_channels; ++ci)
            {
                // Is this gene on or off?
                size_t on_off = att.unchecked_test(reg_base + ci);
                // What information does this carry about the particular
                // category assigned to this environment
                sub[ci][0][on_off][cat] += p_event;
            }
        }
    }
}

cWCAnalyzer::cWCAnalyzer(cWorld_ptr &w,
                         const cIndexes &ind,
                         const cRates &t1, 
                         const cRates &t2,
                         double wght)
    : world_ptr(w)
    , world(*w)
    , indexes(ind)
    , target1(t1)
    , target2(t2)
    , weighting(wght)
    , processing_index(0)
    , empty_probs(world.reg_channels, {{0.0, 0.0}})
{
}

cInformation *cWCAnalyzer::analyse_network(
    cNetwork &net)
{
    // Entropies and mutual info
    cInformation *info = new cInformation(world_ptr, 1, 1);
    _analyse(net, info->_array[0]);
    return info;
}

cInformation *cWCAnalyzer::analyse_collection(
    const cNetworkVector &networks)
{
    cInformation *info = new cInformation(world_ptr, networks.size(), 1);
    for (size_t i = 0; i < networks.size(); ++i)
        _analyse(*networks[i], info->_array[i]);
    return info;
}

inline double cWCAnalyzer::_similarity(const cRates &a, const cRates &b)
{
    // ASSUMES: rates are same size (at least b >= a)
    double dist = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
    {
        double diff = a[i] - b[i];
        dist += diff * diff;
    }
    if (dist == 0.0)
        return 1.0;
    else
        dist = sqrt(dist);

    return exp(-fabs(dist) / weighting);
}


void cWCAnalyzer::_wiggle(cNetwork &net)
{
    processing_index++;

    double pr = 0.5 / (world.environments.size());

    // Note the reversed order here (j, i); it is very expensive to recalculate
    // the attractors. So we do it as little as possible.
    for (size_t j = 0; j < world.reg_channels; ++j)
    {
        cGene *gene = net.get_gene(j);

        // Force gene OFF
        gene->intervene = INTERVENE_OFF;
        net.calc_attractors_with_intervention();

        // For each environment and each output channel
        for (auto const &r : net.rates)
            _add_probability(r, j, 0, pr);

        // Force gene ON
        gene->intervene = INTERVENE_ON;
        net.calc_attractors_with_intervention();

        for (auto const &r : net.rates)
            _add_probability(r, j, 1, pr);

        // Reset the current gene
        gene->intervene = INTERVENE_NONE;
    }

    // Recalculate one more time without using interventions. This just resets
    // everything.
    net.calc_attractors();

}

inline void cWCAnalyzer::_add_probability(const cRates &rate, 
                                          size_t gindex, size_t is_on, double pr)
{
    RateDetail &found = _match(rate);
    found.probs[gindex][is_on] += pr;
}

// Simple linear search
inline cWCAnalyzer::RateDetail& cWCAnalyzer::_match(const cRates &rate)
{
    cRates to_match;
    for (size_t i = 0; i < indexes.size(); ++i)
        to_match.push_back(rate[indexes[i]]);

    rate_detail_map_type::iterator iter = rate_detail_map.find(to_match);
    if (iter == rate_detail_map.end())
    {
        // It is a new one. We need to initialize it and fill out the similarity metric.
        auto result = rate_detail_map.emplace(to_match, 
            RateDetail(rate_detail_map.size(), processing_index, world.reg_channels));
        iter = result.first;
        auto &new_detail = (*iter).second;
        new_detail.similarity[0] = _similarity(to_match, target1);
        new_detail.similarity[1] = _similarity(to_match, target2);
    }
    else
    {
        auto &detail = (*iter).second;
        // Reset everything to zero if we're processing something new
        if (detail.used_for != processing_index)
        {
            detail.used_for = processing_index;
            detail.probs = empty_probs;
        }
    }
    return (*iter).second;
}

void cWCAnalyzer::_analyse(cNetwork &net, info_array_type::reference sub)
{
    _wiggle(net);

    // Construct the row marginals for each regulatory channel
    off_on_vector_type row_sums(world.reg_channels, {{0.0, 0.0}});
    for (auto &rate_detail : rate_detail_map)
    {
        auto &detail = rate_detail.second;
        if (detail.used_for == processing_index)
        {
            for (size_t j = 0; j < world.reg_channels; ++j)
            {
                row_sums[j][0] += detail.probs[j][0];
                row_sums[j][1] += detail.probs[j][1];
            }
        }
    }

    // Now constructed the two weighted info sums
    off_on_vector_type info(world.reg_channels, {{0.0, 0.0}});
    for (auto &rate_detail : rate_detail_map)
    {
        auto &detail = rate_detail.second;

        // Only process currently used ones
        if (detail.used_for == processing_index)
        {
            for (size_t j = 0; j < world.reg_channels; ++j)
            {
                // calculate the column sum
                double colsum = detail.probs[j][0] + detail.probs[j][1];

                // Calculate pointwise mutual info for the two rows 
                for (size_t k = 0; k < 2; ++k)
                {
                    double denom = row_sums[j][k] * colsum;
                    double val = detail.probs[j][k];
                    if (not_zeroish(val) && not_zeroish(denom))
                    {
                        double pmi = val * log2(val / denom);

                        // We want to use the opposite weights for each index...
                        info[j][0] += pmi * detail.similarity[k % 2];
                        info[j][1] += pmi * detail.similarity[(k + 1) % 2];
                    }
                }
            }
        }
    }

    for (size_t j = 0; j < world.reg_channels; ++j)
    {
        // Summarize info into something that we can return via numpy
        sub[j][0] = std::max(info[j][0], info[j][1]);
    }
}


cJointProbabilities *cWCAnalyzer::get_joint(cNetwork &net)
{
    _wiggle(net);

    cJointProbabilities *joint =
        new cJointProbabilities(world_ptr,
                                1,
                                1,
                                rate_detail_map.size());

    for (auto &rate_detail : rate_detail_map)
    {
        auto &detail = rate_detail.second;
        if (detail.used_for == processing_index)
        {
            for (size_t j = 0; j < world.reg_channels; ++j)
            {
                joint->_array[0][j][0][0][detail.encountered] = detail.probs[j][0];
                joint->_array[0][j][0][1][detail.encountered] = detail.probs[j][1];
            }
        }
    }

    return joint;
}
