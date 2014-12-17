#include "core.hpp"

namespace pubsub2 {

// This is the inner-inner loop!
// TODO: Maybe this could be moved down the the CIS level to prevent constant
// calling of virtual function. It would have to be Factory level call:
// cGeneFactory::cycle(Network &, ChannelState &). But this would mean building 
// the cis action into the factory somehow (or static_cast-ing the CIS which
// would, I guess, be safe.
template<typename Network> 
struct Cycle
{
    void cycle_with_intervention(const Network &network, 
                                 cChannelState &c) const
    {
        cChannelState next(c.size());
        for (auto g : network.genes)
            switch (g->intervene)
            {
            case INTERVENE_ON:
                // If it is forced ON, don't both checking anything
                next.set(g->pub);
                break;
            case INTERVENE_NONE:
                for (auto m : g->modules)
                    // NOTE: is_active is what does all the work here. It is
                    // virtual, and we call it a lot. Hence the ideas above about a
                    // better solution.
                    if (m->intervene == INTERVENE_ON || 
                        (m->intervene == INTERVENE_NONE && m->is_active(c)))
                    {
                        next.set(g->pub);
                        // The gene is active, no use looking further
                        break;
                    }
            case INTERVENE_OFF:
                // Do nothing 
                ;
            }
        // Update the "return" value.
        c.swap(next);
    }

    void cycle(const Network &network, 
               cChannelState &c) const
    {
        cChannelState next(c.size());
        for (auto g : network.genes)
            for (auto m : g->modules)
                if (m->is_active(c))
                    {
                        next.set(g->pub);
                        // The gene is active, no use looking further
                        break;
                    }

        // Update the "return" value.
        c.swap(next);
    }
};

}
