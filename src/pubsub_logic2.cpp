#include "pubsub_logic2.hpp"
#include <cmath>
#include <stdexcept>
// #include <iostream>

using namespace pubsub2;

cGeneFactoryLogic2::cGeneFactoryLogic2(cFactory *f, double rate_per_gene_)
    : cGeneFactory(f, rate_per_gene_)
    , r_oper(0, f->operands.size()-1)
{
}

cCisModule *cGeneFactoryLogic2::construct_cis()
{
    cCisModuleLogic2 *cis = new cCisModuleLogic2();
    random_engine_t &re = factory->random_engine;
    cis->op = factory->operands[r_oper(re)];
    cis->channels[0] = r_sub(re);
    cis->channels[1] = r_sub(re);
    return cis;
}

// This is where the action really is.
void cGeneFactoryLogic2::mutate_cis(cCisModule *m)
{
    cCisModuleLogic2 *cis = static_cast<cCisModuleLogic2 *>(m);
    // TODO: This is shite; just a start.
    double p = r_uniform_01(factory->random_engine);
    if (p < .5)
        cis->op = factory->operands[r_oper(factory->random_engine)];
    else if (p < .75)
        cis->channels[0] = r_sub(factory->random_engine);
    else
        cis->channels[1] = r_sub(factory->random_engine);
}

cCisModule *cCisModuleLogic2::clone() const
{
    cCisModuleLogic2 *cis = new cCisModuleLogic2();
    cis->op = op;
    cis->channels[0] = channels[0];
    cis->channels[1] = channels[1];
    return cis;
}
