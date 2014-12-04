#include "scheme_cooperative.h"

using namespace pubsub2;

cSchemeCoop3::cSchemeCoop3(cFactory *f, double rate_per_gene_)
    : cGeneFactory(f, rate_per_gene_)
{
}

cCisModule *cSchemeCoop3::construct_cis()
{
    cCisModuleCoop3 *cis = new cCisModuleCoop3();
    // random_engine_t &re = factory->random_engine;
    // cis->op = factory->operands[r_oper(re)];
    // cis->channels[0] = r_sub(re);
    // cis->channels[1] = r_sub(re);
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
    // cis->op = op;
    // cis->channels[0] = channels[0];
    // cis->channels[1] = channels[1];
    return cis;
}

bool cCisModuleCoop3::is_active(cChannelState const &state) const 
{
    return true;
}
