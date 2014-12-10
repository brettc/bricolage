#include "core.hpp"

namespace pubsub2 {

struct cGeneFactoryThreshold3 : public cGeneFactory
{
    cGeneFactoryThreshold3(cFactory *f, double rate_per_gene_);

    // Randomize binding
    randint_t r_binding, r_direction, r_site, r_input;

    // virtual overrides
    cCisModule *construct_cis();
    void mutate_cis(cCisModule *m);
};

class cCisModuleThreshold3 : public cCisModule
{
public:
    cCisModuleThreshold3() { _channels = channels; }

    // Overrides
    size_t site_count() const { return 3; }
    virtual cCisModule* clone() const;
    void mutate();
    bool is_active(cChannelState const &state) const;

    int_t binding[3];
    signal_t channels[3];
};

}
