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

typedef boost::multi_array_ref<double, 3> freqs_t_ref;
cEnvironmentI::cEnvironmentI(
    const cWorld_ptr &w, size_t ncat)
    : world(w)
    , category_count(ncat)
{
    frequencies = new freqs_t(boost::extents[w->reg_channels][ncat][2]);
    std::fill_n(std::back_inserter(categories), w->environments.size(), 0);
}

cEnvironmentI::~cEnvironmentI()
{
    delete frequencies;
}

void cEnvironmentI::copy_frequencies(double *view) const
{
    freqs_t_ref arr_view(view, boost::extents[world->reg_channels]
                           [category_count][2]);

    // Rely on the copy constructor
    arr_view = *frequencies;

}

void cEnvironmentI::get_extents(size_t &channels, 
                                size_t &categories,
                                size_t &on_off)
{
    channels = world->reg_channels;
    categories = category_count;
    on_off = 2;
}

void cEnvironmentI::calculate(const cNetwork &net)
{
    double normal_p_event = 1.0 / double(world->environments.size());
    size_t reg_base = world->reg_range.first;
    for (size_t i = 0; i < world->reg_channels; ++i)
    {
        for (size_t j = 0; j < net.attractors.size(); ++j)
        {
            double p_event = normal_p_event;
            auto attrs = net.attractors[j];
            p_event /= double(attrs.size());
            for (auto &a : attrs)
            {
                size_t on_off = a.test(reg_base+i);
                (*frequencies)[i][categories[j]][on_off] += p_event;
            }
        }
    }
}

void cEnvironmentI::calc_collection(double *data,
                                    const cNetworkVector &networks)
{
    typedef boost::multi_array_ref<double, 4> array_type;
    size_t a, b, c, d;
    a = networks.size();
    get_extents(b, c, d);
    array_type arr(data, boost::extents[a][b][c][d]);
    double normal_p_event = 1.0 / double(world->environments.size());
    size_t reg_base = world->reg_range.first;

    for (size_t i = 0; i < a; ++i)
    {
        const cNetwork &net = *networks[i];
        for (size_t j = 0; j < b; ++j)
        {
            for (size_t k = 0; k < c; ++k)
            {
                double p_event = normal_p_event;
                auto attrs = net.attractors[k];
                p_event /= double(attrs.size());
                for (auto &a : attrs)
                {
                    size_t on_off = a.test(reg_base+j);
                    arr[i][j][categories[k]][on_off] += p_event;
                }
            }
        }
    }
}

