#include "core.hpp"
#include <cmath>
#include <stdexcept>

using namespace pubsub2;

cNetworkAnalysis::cNetworkAnalysis(const cNetwork_ptr &n)
    : original(n)
    , modified(n->factory, true)
{
    // Just copy the genes for now. Don't bother doing anything else.
    original->clone_genes(modified.genes);
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
        const cGene *g = original->genes[i];
        Node_t gnode = std::make_pair(NT_GENE, i);
        Node_t cnode = std::make_pair(NT_CHANNEL, g->pub);
        edges.emplace(gnode, cnode);
        for (size_t j=0; j < g->module_count(); ++j)
        {
            const cCisModule *m = g->modules[j];
            Node_t mnode = std::make_pair(NT_MODULE, make_module_node_id(i, j));
            edges.emplace(mnode, gnode);
            for (size_t k=0; k < m->site_count(); ++k)
            {
                Node_t cnode = std::make_pair(NT_CHANNEL, m->get_site_channel(k));
                edges.emplace(cnode, mnode);
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
    for (size_t i=0; i < modified.gene_count(); ++i)
    {
        cGene *g = modified.genes[i];

        // Can we knock this gene out?
        g->intervene = INTERVENE_OFF;
        modified.calc_attractors();

        // If this makes no diff, then no edges attached to this Gene make any
        // difference. Let's skip to the next one.
        if (modified.rates == original->rates)
            continue;

        // Ok, it does something ...
        g->intervene = INTERVENE_NONE;
        Node_t gnode = std::make_pair(NT_GENE, i);
        Node_t cnode = std::make_pair(NT_CHANNEL, g->pub);
        edges.emplace(gnode, cnode);
        
        for (size_t j=0; j < g->module_count(); ++j)
        {
            cCisModule *m = g->modules[j];
            m->intervene = INTERVENE_OFF;
            modified.calc_attractors();
            if (modified.rates == original->rates)
                continue;

            // Reset and add edge
            m->intervene = INTERVENE_NONE;
            Node_t mnode = std::make_pair(NT_MODULE, make_module_node_id(i, j));
            edges.emplace(mnode, gnode);

            for (size_t k=0; k < m->site_count(); ++k)
            {
                // Knock this out by setting it to ZERO channel. Nothing ever
                // publishes to the ZERO channel!
                signal_t old = m->set_site_channel(k, 0);
                Node_t cnode = std::make_pair(NT_CHANNEL, old);

                // It was already zero?
                if (old == 0)
                    continue;

                // Otherwise, we need to test to see if it changed anything
                modified.calc_attractors();

                // Did it change? 
                if (modified.rates != original->rates)
                {
                    // It did, so this channel makes a difference to the module
                    edges.emplace(cnode, mnode); 

                    // ... and we need to reset it
                    m->set_site_channel(k, old);
                }
            }
        }
    }
}

