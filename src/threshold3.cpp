#include "threshold3.hpp"

using namespace pubsub2;

cConstructorThreshold3::cConstructorThreshold3(
    cFactory &f, size_t gc, size_t cc)
    : cConstructor(f, gc, cc)
    , r_binding(-3, 3)
    , r_direction(0, 1)
    , r_site(0, 2)
    , r_input(0, f.sub_range.second-1)
{
}

cCisModule *cConstructorThreshold3::construct_cis()
{
    cCisModuleThreshold3 *cis = new cCisModuleThreshold3();

    for (size_t i = 0; i < cis->site_count(); ++i)
    {
        cis->channels[i] = r_input(factory.rand);
        cis->binding[i] = r_binding(factory.rand);
    }
    return cis;
}

// This is where the action really is.
void cConstructorThreshold3::mutate_cis(cCisModule *m)
{
    cCisModuleThreshold3 *cis = static_cast<cCisModuleThreshold3 *>(m);
    size_t site = r_site(factory.rand);
    int_t current = cis->binding[site];

    int_t mutate;
    if (current == 3)
        mutate = -1;
    else if (current == -3)
        mutate = 1;
    else
        mutate = r_direction(factory.rand) * 2 - 1;

    // If we're at zero, possibly change into a different binding 
    if (current == 0)
        cis->channels[site] = r_input(factory.rand);

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
