#include "core.hpp"
#include <cmath>
#include <stdexcept>
#include <algorithm>


using namespace bricolage;

// Taken from the numpy docs on isclose
const double RELATIVE_TOL = 1e-05;
const double ABSOLUTE_TOL = 1e-08;
inline bool isclose(double a, double b)
{
    return fabs(a - b) <= (ABSOLUTE_TOL + RELATIVE_TOL * fabs(b));
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
                Node_t cnode = std::make_pair(NT_CHANNEL, m->get_site(k));
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
                signal_t old = m->set_site(k, 0);
                Node_t cnode = std::make_pair(NT_CHANNEL, old);

                // It was already zero?
                if (old == 0)
                    continue;

                // Otherwise, we need to test to see if it changed anything
                modified->calc_attractors_with_intervention();

                // Did it change? 
                if (modified->rates != original->rates)
                {
                    // It did, so this channel makes a difference to the module
                    edges.insert(Edge_t(cnode, mnode));

                    // ... and we need to reset it
                    m->set_site(k, old);
                }
            }
        }
    }
}

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

cJointProbabilities::cJointProbabilities(const cWorld_ptr &w, size_t network_size,
                                         size_t per_network, size_t per_channel)
    : world(w)
    , _array(boost::extents
             [network_size]
             [w->reg_channels]
             [per_network]
             [2][per_channel])
{
}

void cJointProbabilities::calc_information(cInformation &info) const
{
    auto networks_size = _array.shape()[0];
    auto channels_size = _array.shape()[1];
    auto category_size = _array.shape()[2];
    auto row_size = _array.shape()[3];
    auto col_size = _array.shape()[4];

    typedef joint_array_type::index index;

    cRates rows(row_size);
    cRates cols(col_size);

    for (index i = 0; i < networks_size; ++i)
    {
        for (index j = 0; j < channels_size; ++j)
        {
            for (index k = 0; k < category_size; ++k)
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
                        if (val != 0.0)
                            I += val * log2(val / (rows[ri] * cols[ci]));
                    }
                info._array[i][j][k] = I;
            }
        }
    }
}

cCausalFlowAnalyzer::cCausalFlowAnalyzer(cWorld_ptr &w, const cRates &rates)
    : world(w)
    , rates(rates)
    , natural_probabilities(w->reg_channels, 0.0)
{

}

cJointProbabilities *cCausalFlowAnalyzer::analyse_network(cNetwork &net)
{
    cJointProbabilities *joint = 
        new cJointProbabilities(world, 1, world->out_channels, 2);

    _analyse(net, joint->_array[0]);

    return joint;
}

cJointProbabilities *cCausalFlowAnalyzer::analyse_collection(const cNetworkVector &networks)
{
    cJointProbabilities *joint = 
        new cJointProbabilities(world, networks.size(), world->out_channels, 2);

    for (size_t i = 0; i < networks.size(); ++i)
        _analyse(*networks[i], joint->_array[i]);
    return joint;
}

// Calculate the probability of any particular signal being on when the
// attractor network is not being intervened upon
void cCausalFlowAnalyzer::_calc_natural(cNetwork &net)
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

void cCausalFlowAnalyzer::_analyse(cNetwork &net, joint_array_type::reference sub)
{
    _calc_natural(net);

    double p_env = 1.0 / net.rates.size();
    size_t close = 0;

    // For each channel
    for (size_t i = 0; i < world->reg_channels; ++i)
    {
        cGene *gene = net.get_gene(i);
        gene->intervene = INTERVENE_OFF;
        net.calc_attractors_with_intervention();
        double p_gene_off = p_env * (1.0 - natural_probabilities[i]);

        // for each environment (there are rates for each) 
        for (size_t j = 0; j < net.rates.size(); ++j)
            // for each output channel
            for (size_t k = 0; k < world->out_channels; ++k)
            {
                if (isclose(rates[k], net.rates[j][k]))
                    close = 1;
                else
                    close = 0;
                    
                sub[i][k][0][close] += p_gene_off;
            }

        gene->intervene = INTERVENE_ON;
        net.calc_attractors_with_intervention();
        double p_gene_on = p_env * natural_probabilities[i];
        // for each environment (there are rates for each) 
        for (size_t j = 0; j < net.rates.size(); ++j)
            // for each output channel
            for (size_t k = 0; k < world->out_channels; ++k)
            {
                if (isclose(rates[k], net.rates[j][k]))
                    close = 1;
                else
                    close = 0;
                    
                sub[i][k][1][close] += p_gene_on;
            }

        // Reset
        gene->intervene = INTERVENE_NONE;
        net.calc_attractors();
    }
}


// Pass in the maximum numbers of categories
cMutualInfoAnalyzer::cMutualInfoAnalyzer(cWorld_ptr &w, const cIndexes &cats)
    : world(w)
    , categories(cats)
    , max_category(1 + *std::max_element(categories.begin(), categories.end()))
{
}

cJointProbabilities *cMutualInfoAnalyzer::analyse_network(cNetwork &net)
{
    cJointProbabilities *joint = 
        new cJointProbabilities(world, 1, 1, max_category);

    _analyse(net, joint->_array[0]);

    return joint;
}

cJointProbabilities *cMutualInfoAnalyzer::analyse_collection(const cNetworkVector &networks)
{
    cJointProbabilities *joint = 
        new cJointProbabilities(world, networks.size(), 1, max_category);

    for (size_t i = 0; i < networks.size(); ++i)
        _analyse(*networks[i], joint->_array[i]);
    return joint;
}

void cMutualInfoAnalyzer::_analyse(cNetwork &net, joint_array_type::reference sub)
{
    size_t reg_base = world->reg_range.first;
    double normal_p_event = 1.0 / double(world->environments.size());

    // For each environment
    for (size_t ei = 0; ei < net.attractors.size(); ++ei)
    {
        auto attrs = net.attractors[ei];
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



