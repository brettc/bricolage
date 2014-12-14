#include "logic2.hpp"
#include <cmath>
#include <stdexcept>
// #include <iostream>

using namespace pubsub2;

cConstructorLogic2::cConstructorLogic2(
    cFactory &f, size_t gc, size_t cc, const cOperands &ops
    )
    : cConstructor(f, gc, cc)
    , operands(ops)
    , r_oper(0, ops.size()-1)
{
}

cCisModule *cConstructorLogic2::construct_cis()
{
    cCisModuleLogic2 *cis = new cCisModuleLogic2();
    cis->op = operands[r_oper(factory.rand)];
    cis->channels[0] = r_sub(factory.rand);
    cis->channels[1] = r_sub(factory.rand);
    return cis;
}

// This is where the action really is.
void cConstructorLogic2::mutate_cis(cCisModule *m)
{
    cCisModuleLogic2 *cis = static_cast<cCisModuleLogic2 *>(m);
    // TODO: This is shite; just a start.
    double p = r_uniform_01(factory.rand);
    if (p < .5)
        cis->op = operands[r_oper(factory.rand)];
    else if (p < .75)
        cis->channels[0] = r_sub(factory.rand);
    else
        cis->channels[1] = r_sub(factory.rand);
}

cCisModule *cCisModuleLogic2::clone() const
{
    cCisModuleLogic2 *cis = new cCisModuleLogic2();
    cis->op = op;
    cis->channels[0] = channels[0];
    cis->channels[1] = channels[1];
    return cis;
}
