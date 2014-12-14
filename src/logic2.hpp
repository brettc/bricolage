#pragma once

#include "core.hpp"

namespace pubsub2 {

struct cConstructorLogic2 : public cConstructor
{
    cConstructorLogic2(cFactory &f, size_t gene_count_, size_t cis_count,
                       const cOperands &ops);
    cOperands operands;
    randint_t r_oper;

    // Overrides
    cCisModule *construct_cis();
    void mutate_cis(cCisModule *m);
};

class cCisModuleLogic2 : public cCisModule
{
public:
    cCisModuleLogic2() { _channels = channels; }

    size_t site_count() const { return 2; }
    virtual cCisModule* clone() const;
    void mutate();

    // Inline this stuff. It won't change.
    inline bool test(unsigned int a, unsigned int b) const 
    { 
        // Note: C++ standard guarantees integral conversion from bool results
        // in 0 or 1.
        return op & (8 >> ((a << 1) | b)); 
    }

    inline bool is_active(cChannelState const &state) const 
    {
        return test(state.test(channels[0]), state.test(channels[1]));
    }
// protected:
    // Default constructor is fine
    operand_t op;
    signal_t channels[2];
};

}
