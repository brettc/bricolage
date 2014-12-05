#include "scheme_cooperative.h"

using namespace pubsub2;

cSchemeCoop3::cSchemeCoop3(cFactory *f, double rate_per_gene_)
    : cGeneFactory(f, rate_per_gene_)
    , r_binding(-3, 3)
{
}

cCisModule *cSchemeCoop3::construct_cis()
{
    cCisModuleCoop3 *cis = new cCisModuleCoop3();
    random_engine_t &re = factory->random_engine;

    for (size_t i = 0; i < cis->site_count(); ++i)
    {
        cis->channels[i] = r_sub(re);
        cis->binding[i] = r_binding(re);
    }
    return cis;
}

// This is where the action really is.
void cSchemeCoop3::mutate_cis(cCisModule *m)
{
    cCisModuleCoop3 *cis = static_cast<cCisModuleCoop3 *>(m);
}

cCisModule *cCisModuleCoop3::clone() const
{
    cCisModuleCoop3 *cis = new cCisModuleCoop3();
    for (size_t i = 0; i < cis->site_count(); ++i)
    {
        cis->channels[i] = channels[i];
        cis->binding[i] = binding[i];
    }
    return cis;
}

bool cCisModuleCoop3::is_active(cChannelState const &state) const 
{
    // Calculate the weighted sum
    int_t sum = 0;
    for (size_t i = 0; i < site_count(); ++i)
    {
        if (state[channels[i]])
            sum += binding[i];
    }
    // Thresholded
    return sum >= 3;
}
