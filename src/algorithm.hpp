#include "core.hpp"

namespace algo {

// Templatize this so we can re-use it for different networks.
//
// Requirements:
// * Network must have .genes
// * Module must have .modules
// * Module must have .is_active(cChannelState &)
template<typename Network> 
struct Cycle
{
    void cycle_with_intervention(const Network &network, 
                                 bricolage::cChannelState &c) const
    {
        bricolage::cChannelState next(c.size());
        for (auto &g : network.genes)
            switch (g.intervene)
            {
            case bricolage::INTERVENE_ON:
                // If it is forced ON, don't both checking anything
                next.set(g.pub);
                break;
            case bricolage::INTERVENE_NONE:
                for (auto &m : g.modules)
                    // NOTE: is_active is what does all the work here. It is
                    // virtual, and we call it a lot. Hence the ideas above about a
                    // better solution.
                    if (m.intervene == bricolage::INTERVENE_ON || 
                        (m.intervene == bricolage::INTERVENE_NONE && m.is_active(c)))
                    {
                        next.set(g.pub);
                        // The gene is active, no use looking further
                        break;
                    }
            case bricolage::INTERVENE_OFF:
                // Do nothing 
                ;
            }
        // Update the "return" value.
        c.swap(next);
    }

    void cycle(const Network &network, 
               bricolage::cChannelState &c) const
    {
        bricolage::cChannelState next(c.size());
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
