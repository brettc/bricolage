#include "threshold3.hpp"

using namespace pubsub2;

cGeneFactoryThreshold3::cGeneFactoryThreshold3(cFactory *f, double rate_per_gene_)
    : cGeneFactory(f, rate_per_gene_)
    , r_binding(-3, 3)
    , r_direction(0, 1)
    , r_site(0, 2)
{
}

cCisModule *cGeneFactoryThreshold3::construct_cis()
{
    cCisModuleThreshold3 *cis = new cCisModuleThreshold3();
    random_engine_t &re = factory->random_engine;

    for (size_t i = 0; i < cis->site_count(); ++i)
    {
        cis->channels[i] = r_sub(re);
        cis->binding[i] = r_binding(re);
    }
    return cis;
}

// This is where the action really is.
void cGeneFactoryThreshold3::mutate_cis(cCisModule *m)
{
    random_engine_t &re = factory->random_engine;
    cCisModuleThreshold3 *cis = static_cast<cCisModuleThreshold3 *>(m);
    size_t site = r_site(re);
    int_t current = cis->binding[site];

    int_t mutate;
    if (current == 3)
        mutate = -1;
    else if (current == -3)
        mutate = 1;
    else
        mutate = r_direction(re) * 2 - 1;

    // If we're at zero, possibly change into a different binding 
    if (current == 0)
        cis->channels[site] = r_sub(re);

    cis->binding[site] += mutate;
}

cCisModule *cCisModuleThreshold3::clone() const
{
    cCisModuleThreshold3 *cis = new cCisModuleThreshold3();
    for (size_t i = 0; i < cis->site_count(); ++i)
    {
        cis->channels[i] = channels[i];
        cis->binding[i] = binding[i];
    }
    return cis;
}

bool cCisModuleThreshold3::is_active(cChannelState const &state) const 
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
