#include "core.hpp"

namespace algo {

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
                                 pubsub2::cChannelState &c) const
    {
        pubsub2::cChannelState next(c.size());
        for (auto &g : network.genes)
            switch (g.intervene)
            {
            case pubsub2::INTERVENE_ON:
                // If it is forced ON, don't both checking anything
                next.set(g.pub);
                break;
            case pubsub2::INTERVENE_NONE:
                for (auto &m : g.modules)
                    // NOTE: is_active is what does all the work here. It is
                    // virtual, and we call it a lot. Hence the ideas above about a
                    // better solution.
                    if (m.intervene == pubsub2::INTERVENE_ON || 
                        (m.intervene == pubsub2::INTERVENE_NONE && m.is_active(c)))
                    {
                        next.set(g.pub);
                        // The gene is active, no use looking further
                        break;
                    }
            case pubsub2::INTERVENE_OFF:
                // Do nothing 
                ;
            }
        // Update the "return" value.
        c.swap(next);
    }

    void cycle(const Network &network, 
               pubsub2::cChannelState &c) const
    {
        pubsub2::cChannelState next(c.size());
        for (auto &g : network.genes)
            for (auto &m : g.modules)
                if (m.is_active(c))
                    {
                        next.set(g.pub);
                        // The gene is active, no use looking further
                        break;
                    }

        // Update the "return" value.
        c.swap(next);
    }
};

}
