#include "core.hpp"
#include <cmath>
#include <stdexcept>


using namespace pubsub2;

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

// TODO: Fix cut and paste, once I figure out multi_array propertly (using
// views?) 

cInfoE::cInfoE(const cWorld_ptr &w, size_t ncat)
    : world(w)
    , category_count(ncat)
{
    std::fill_n(std::back_inserter(categories), w->environments.size(), 0);
}

void cInfoE::get_extents(size_t &channels, 
                                size_t &categories,
                                size_t &on_off)
{
    channels = world->reg_channels;
    categories = category_count;
    on_off = 2;
}

typedef boost::multi_array<double, 4> prob_array_type;
// typedef boost::multi_array<double, 3> prob_array_net_type;
typedef boost::multi_array_ref<double, 4> prob_array_ref_type;
typedef boost::multi_array_ref<double, 2> info_array_ref_type;
typedef std::vector<double> sum_type;


void calc_probs(cInfoE &ei, prob_array_ref_type &arr,
                     const cNetworkVector &networks)
{
    double normal_p_event = 1.0 / double(ei.world->environments.size());
    size_t reg_base = ei.world->reg_range.first;

    for (size_t i = 0; i < networks.size(); ++i)
    {
        const cNetwork &net = *networks[i];
        // calc_probs_network(ei, net, arr[i]);
        for (size_t j = 0; j < ei.world->reg_channels; ++j)
        {
            for (size_t k = 0; k < net.attractors.size(); ++k)
            {
                double p_event = normal_p_event;
                auto attrs = net.attractors[k];
                p_event /= double(attrs.size());
                for (auto &a : attrs)
                {
                    size_t on_off = a.test(reg_base+j);
                    arr[i][j][ei.categories[k]][on_off] += p_event;
                }
            }
        }
    }
}

void cInfoE::network_probs(double *data, const cNetwork &net)
{
    typedef boost::multi_array_ref<double, 3> array_type;
    size_t chan_size, cat_size, stat_size;
    get_extents(chan_size, cat_size, stat_size);
    array_type arr(data, boost::extents[chan_size][cat_size][stat_size]);

    double normal_p_event = 1.0 / double(world->environments.size());
    size_t reg_base = world->reg_range.first;

    for (size_t j = 0; j < chan_size; ++j)
    {
        for (size_t k = 0; k < net.attractors.size(); ++k)
        {
            double p_event = normal_p_event;
            auto attrs = net.attractors[k];
            p_event /= double(attrs.size());
            for (auto &att : attrs)
            {
                size_t on_off = att.test(reg_base+j);
                arr[j][categories[k]][on_off] += p_event;
            }
        }
    }
}


void cInfoE::collection_probs(double *data,
                              const cNetworkVector &networks)
{
    size_t net_size, chan_size, cat_size, stat_size;
    net_size = networks.size();
    get_extents(chan_size, cat_size, stat_size);
    prob_array_ref_type arr(data, boost::extents[net_size][chan_size][cat_size][stat_size]);

    ::calc_probs(*this, arr, networks);
}

void calc_info(cInfoE &ei, 
               info_array_ref_type::reference i_arr,
               prob_array_ref_type::reference p_arr,
               sum_type &feat_sum, 
               sum_type &chan_sum,
               size_t chan_size, size_t cat_size, size_t stat_size
               )
{
    // iterate over channels
    for (size_t j = 0; j < chan_size; ++j)
    {
        // Reset to 0
        for (size_t k = 0; k < cat_size; ++k)
            feat_sum[k] = 0.0;
        for (size_t l = 0; l < stat_size; ++l)
            chan_sum[l] = 0.0;

        // Sum the marginals
        for (size_t k = 0; k < cat_size; ++k)
            for (size_t l = 0; l < stat_size; ++l)
            {
                double val = p_arr[j][k][l];
                feat_sum[k] += val;
                chan_sum[l] += val;
            }

        // Calculate the info
        double I = 0.0;
        for (size_t k = 0; k < cat_size; ++k)
        {
            for (size_t l = 0; l < stat_size; ++l)
            {
                double val = p_arr[j][k][l];
                if (val != 0.0)
                    I += val * log2(val / (feat_sum[k] * chan_sum[l]));
            }
        }
        i_arr[j] = I;
    }
}

void cInfoE::collection_info(double *data,
                              const cNetworkVector &networks)
{
    size_t net_size, chan_size, cat_size, stat_size;
    net_size = networks.size();
    get_extents(chan_size, cat_size, stat_size);

    // NOTE: We allocate this ourself this time, not like above
    prob_array_type arr(boost::extents[net_size][chan_size][cat_size][stat_size]);
    ::calc_probs(*this, arr, networks);

    sum_type feat_sum(cat_size, 0.0);
    sum_type chan_sum(stat_size, 0.0);

    info_array_ref_type info(data, boost::extents[net_size][chan_size]);
    for (size_t i = 0; i < networks.size(); ++i)
    {
        ::calc_info(*this, info[i], arr[i], feat_sum, chan_sum,
                    chan_size, cat_size, stat_size);
    }
}

// void cInfoE::calc_info_network(double *data, const cNetwork &network)
// void network_info(double *data, const cNetwork &network)
// {
//     size_t net_size, chan_size, cat_size, stat_size;
//     get_extents(chan_size, cat_size, stat_size);
//
//     // NOTE: We allocate this ourself this time, not like above
//     prob_array_net_type arr(boost::extents[chan_size][cat_size][stat_size]);
//     ::calc_collection(*this, arr, networks);
//
//     sum_type feat_sum(cat_size, 0.0);
//     sum_type chan_sum(stat_size, 0.0);
//
//     info_array_ref_type info(data, boost::extents[net_size][chan_size]);
//     for (size_t i = 0; i < networks.size(); ++i)
//     {
//         ::calc_info(*this, info[i], arr[i], feat_sum, chan_sum,
//                     chan_size, cat_size, stat_size);
//     }
// }

//   ACTULLY: looks possible to just use
//  prob_array_type::reference to get at the sub_array 
// void calc_network(boost::multi_array_ref<double, 4>::reference m)
// {
//     m[0][0][0] = 10.0;
// }


