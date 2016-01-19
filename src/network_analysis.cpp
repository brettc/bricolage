#include "core.hpp"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <map>
#include <exception>
#include <sstream>


using namespace bricolage;

// Taken from the numpy docs on is_close
const double RELATIVE_TOL = 1e-05;
const double ABSOLUTE_TOL = 1e-08;
inline bool is_close(double a, double b)
{
    return fabs(a - b) <= (ABSOLUTE_TOL + RELATIVE_TOL * fabs(b));
}

inline bool not_zeroish(double a)
{
    return fabs(a) > ABSOLUTE_TOL;
}

cNetworkAnalysis::cNetworkAnalysis(cNetwork_ptr &n)
    : original(n)
    , modified(original->clone())
{
}

void cNetworkAnalysis::make_edges(cEdgeList &edges)
{
    // Generate three types of edges
    // Gene -> Channel
    // Channel -> Module
    // Module -> Gene
    edges.clear();
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

                }
            }
        }
    }
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
    , natural_probabilities(w->reg_channels, 0.0)
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
        natural_probabilities[i] = 0.0;

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
                if (cs.test(reg_base + i))
                    natural_probabilities[i] += p_state;
            }
        }
    }
}

//-----------------------------------------------------------------------------
// Keep a map of the unique rates that are output and assign them to persistent
// categories within a network. We'll use these categories to calculate the
// information.
struct RateCategorizer
{
    // We only allocate this many categories. More than this and we're screwed.
    size_t next_category;
    std::map<double, int> rate_categories;

    RateCategorizer() : next_category(2) {}
    size_t get_category(double rate)
    {
        if (rate == 0.0)
            return 0;
        if (rate == 1.0)
            return 1.0;

        auto result = rate_categories.insert(
            std::make_pair(rate, next_category));
        if (result.second)
        {
            // Successful insert. We have a new category. Update the category.
            if (++next_category > cBaseCausalAnalyzer::max_category)
                throw std::out_of_range("Maximum categories reached!!");
        }

        // In either case, return the value in the pair the insert points to.
        // (It will either be the new pair, or the one found).
        auto the_pair = *result.first;
        return the_pair.second;
    }
};




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

    RateCategorizer categorizer;
    double p_env = 1.0 / net.rates.size();

    // For each channel
    for (size_t i = 0; i < world->reg_channels; ++i)
    {
        cGene *gene = net.get_gene(i);
        gene->intervene = INTERVENE_OFF;
        net.calc_attractors_with_intervention();
        double p_gene_off = p_env * (1.0 - natural_probabilities[i]);

        for (size_t j = 0; j < net.rates.size(); ++j)
            for (size_t k = 0; k < world->out_channels; ++k)
            {
                size_t cat = categorizer.get_category(net.rates[j][k]);
                sub[i][k][0][cat] += p_gene_off;
            }

        gene->intervene = INTERVENE_ON;
        net.calc_attractors_with_intervention();
        double p_gene_on = p_env * natural_probabilities[i];

        // for each environment (there are rates for each)
        for (size_t j = 0; j < net.rates.size(); ++j)
            for (size_t k = 0; k < world->out_channels; ++k)
            {
                size_t cat = categorizer.get_category(net.rates[j][k]);
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
{
}

// Note you need to delete the return values from these!
cInformation *cAverageControlAnalyzer::analyse_network(cNetwork &net)
{
    cInformation *info = new cInformation(world, 1, world->out_channels);
    _analyse(net, info->_array[0]);
    return info;
}

cInformation *cAverageControlAnalyzer::analyse_collection(
    const cNetworkVector &networks)
{
    cInformation *info = new cInformation(world, networks.size(), world->out_channels);
    for (size_t i = 0; i < networks.size(); ++i)
        _analyse(*networks[i], info->_array[i]);
    return info;
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
    auto env_count = world->environments.size();

    // We keep a set of category assignments for the rates
    RateCategorizer categorizer;
    cJointProbabilities joint_over_envs(
        world, env_count, world->out_channels, cBaseCausalAnalyzer::max_category);

    // Effectively the same as above
    _calc_natural(net);

    // Note the reversed order here (j, i, k). I've exchanged the loops because
    // it is very expensive to recalculate the attractors. So we do it as
    // little as possible.
    for (size_t j = 0; j < world->reg_channels; ++j)
    {
        cGene *gene = net.get_gene(j);
        double p_gene_on = natural_probabilities[j];
        double p_gene_off = 1.0 - p_gene_on;

        // Force gene OFF
        gene->intervene = INTERVENE_OFF;
        net.calc_attractors_with_intervention();
        // for each environment and each output channel
        for (size_t i = 0; i < net.rates.size(); ++i)
            for (size_t k = 0; k < world->out_channels; ++k)
            {
                size_t cat = categorizer.get_category(net.rates[i][k]);
                joint_over_envs._array[i][j][k][0][cat] += p_gene_off;
            }

        // Force gene ON
        gene->intervene = INTERVENE_ON;
        net.calc_attractors_with_intervention();
        // for each environment and each output channel
        for (size_t i = 0; i < net.rates.size(); ++i)
            for (size_t k = 0; k < world->out_channels; ++k)
            {
                size_t cat = categorizer.get_category(net.rates[i][k]);
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
    double p_env = 1.0 / net.rates.size();
    for (size_t i = 0; i < net.rates.size(); ++i)
        for (size_t j = 0; j < world->reg_channels; ++j)
            for (size_t k = 0; k < world->out_channels; ++k)
                sub[j][k] += info._array[i][j][k] * p_env;

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
                size_t on_off = att.test(reg_base+ci);
                // What information does this carry about the particular
                // category assigned to this environment
                sub[ci][0][on_off][cat] += p_event;
            }
        }
    }
}

//-----------------------------------------------------------------------------
size_t cOutputCategorizer::get_category(const cRates &rates, double prob)
{
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
}

// --------------------------------------------------------------------------
cOutputControlAnalyzer::cOutputControlAnalyzer(cWorld_ptr &world)
    : cBaseCausalAnalyzer(world)
    , categorizers(world->reg_channels)
    , joint_over_envs(world, world->environments.size(),
                      1, // only 1 per channel in this case.
                      cBaseCausalAnalyzer::max_category)
{
}

// Note you need to delete the return values from these!
cInformation *cOutputControlAnalyzer::analyse_network(cNetwork &net)
{
    // One network. 2 information measures. 0 is causal power. 1 is output
    // entropy
    cInformation *info = new cInformation(world, 1, 2);
    _analyse(net, info->_array[0]);
    return info;
}

cInformation *cOutputControlAnalyzer::analyse_collection(
    const cNetworkVector &networks)
{
    // Many networks. 2 information measures. 0 is causal power. 1 is output
    // entropy
    cInformation *info = new cInformation(world, networks.size(), 2);
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
        double p_gene_on = natural_probabilities[j];
        double p_gene_off = 1.0 - p_gene_on;

        // Force gene OFF
        gene->intervene = INTERVENE_OFF;
        net.calc_attractors_with_intervention();
        // for each environment and each output channel
        if (p_gene_off != 0.0)
            for (size_t i = 0; i < net.rates.size(); ++i)
            {
                size_t cat = categorizers[j].get_category(net.rates[i], p_gene_off);
                joint_over_envs._array[i][j][0][0][cat] += p_gene_off;
            }

        // Force gene ON
        gene->intervene = INTERVENE_ON;
        net.calc_attractors_with_intervention();
        // for each environment and each output channel
        if (p_gene_on != 0.0)
            for (size_t i = 0; i < net.rates.size(); ++i)
            {
                size_t cat = categorizers[j].get_category(net.rates[i], p_gene_on);
                joint_over_envs._array[i][j][0][1][cat] += p_gene_on;
            }

        // Reset
        gene->intervene = INTERVENE_NONE;
    }

    // Recalculate one more time without using interventions. This just resets everything.
    net.calc_attractors();

    // Now calculate the information in each of these environments.
    cInformation info(world, world->environments.size(), 2);

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
            double entropy = 0.0;
            for (auto pr : categ.category_probabilities)
            {
                // We need to normalize this; it's not done above.
                double cond_pr = p_env * pr;
                if (not_zeroish(cond_pr))
                {
                    double result = cond_pr * -log2(cond_pr);
                    // if (result != result)
                    // {
                    //     std::ostringstream o;
                    //     o << "Nan:" << cond_pr;
                    //     throw std::out_of_range(o.str());
                    // }
                    entropy += result;
                }
            }
            sub[j][1] = entropy;
        }
}
