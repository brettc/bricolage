#include "pubsub2_c.h"

namespace pubsub2 {

struct cSchemeCoop3 : public cGeneFactory
{
    cSchemeCoop3(cFactory *f, double rate_per_gene_);

    // Randomize binding
    randint_t r_binding;

    // virtual overrides
    cCisModule *construct_cis();
    void mutate_cis(cCisModule *m);
};

class cCisModuleCoop3 : public cCisModule
{
public:
    cCisModuleCoop3() { _channels = channels; }

    // Overrides
    size_t site_count() const { return 3; }
    virtual cCisModule* clone() const;
    void mutate();
    bool is_active(cChannelState const &state) const;

    int_t binding[3];
    signal_t channels[3];
};

}
