#include "core.hpp"

namespace algo {

// Templatize this so we can re-use it for different networks.
//
// Requirements:
// * Network must have .genes
// * Module must have .modules
// * Module must have .is_active(cChannels &)
template<typename Network> 
struct Cycle
{
    void cycle_with_intervention(const Network &network, 
                                 bricolage::cChannels &c) const
    {
        bricolage::cChannels next;
        for (auto &g : network.genes)
            switch (g.intervene)
            {
            case bricolage::INTERVENE_ON:
                // If it is forced ON, don't both checking anything
                next.unchecked_set(g.pub);
                break;
            case bricolage::INTERVENE_NONE:
                for (auto &m : g.modules)
                    // NOTE: is_active is what does all the work here. It is
                    // virtual, and we call it a lot. Hence the ideas above about a
                    // better solution.
                    if (m.intervene == bricolage::INTERVENE_ON || 
                        (m.intervene == bricolage::INTERVENE_NONE && m.is_active(c)))
                    {
                        next.unchecked_set(g.pub);
                        // The gene is active, no use looking further
                        break;
                    }
            case bricolage::INTERVENE_OFF:
                // Do nothing 
                ;
            }
        // Update the "return" value.
        c.bits = next.bits;
    }

    void cycle(const Network &network, 
               bricolage::cChannels &c) const
    {
        bricolage::cChannels next;
        for (auto &g : network.genes)
            for (auto &m : g.modules)
                if (m.is_active(c))
                {
                    next.unchecked_set(g.pub);
                    // The gene is active, no use looking further
                    break;
                }

        // Update the "return" value.
        c.bits = next.bits;
    }
};

template<typename Network, typename Factory>
struct Mutator
{
    void mutate(Network &net, size_t n_cis, size_t n_trans)
    {
        auto &fty = static_cast<const Factory &>(*net.factory);

        // Select the genes that should be mutated
        while (n_cis > 0)
        {
            // Choose a gene and mutate it
            size_t i = fty.r_gene();
            auto &g = net.genes[i];
            
            // Choose a module
            size_t j = fty.r_module();
            auto &c = g.modules[j];

            // Let the module do it's own mutation
            c.mutate(fty);
            --n_cis;
        }

        while (n_trans > 0)
        {
            // Choose a gene and mutate it
            size_t i = fty.r_reg_gene();
            auto &g = net.genes[i];
            g.pub = fty.draw_from_regs[fty.r_reg()];
            
            --n_trans;
        }
    }

    void duplicate(Network &net, size_t n_dups)
    {
        auto &fty = static_cast<const Factory &>(*net.factory);

        while (n_dups > 0)
        {
            // Choose a gene and mutate it
            size_t i, j;
            i = fty.r_reg_gene();
            do 
            {
                j = fty.r_reg_gene();
            } while (i == j);

            net.genes[j] = net.genes[i];
            --n_dups;
        }
    }

};

// template<typename Network, typename Factory, typename Gene, typename CisModule>
// struct Constructor
// {
//     Network *construct(Factory &f, bool create)
//     {
//         auto p = f.shared_from_this(); 
//         cNetwork *net = new cNetwork(p);
//         if (create)
//             net->identifier = world->get_next_network_ident();
//
//         size_t gindex = 0;
//         size_t pub = world->reg_range.first;
//         // Create the regulatory genes
//         for (; gindex < f.world->reg_gene_count; ++gindex)
//         {
//             net->genes.emplace_back(gindex, pub);
//             auto &g = net->genes.back();
//
//             ++pub;
//
//             // Cycle repeatedly through these
//             if (pub == f.world->pub_range.second)
//                 pub = world->reg_range.first;
//         }
//         // Create the structural genes
//         for (size_t st_gene = 0; st_gene < f.world->out_channels; ++st_gene, ++gindex)
//         {
//
//         for (size_t pub = world->reg_range.first, gindex = 0; 
//             pub < world->pub_range.second; ++pub, ++gindex)
//         {
//             net->genes.emplace_back(gindex, pub);
//             auto &g = net->genes.back();
//
//             for (size_t j=0; j < module_count; ++j)
//             {
//                 g.modules.emplace_back(cCisModule(*this));
//                 if (create)
//                 {
//                     auto &m = g.modules.back();
//                     m.op = operands[r_operand()];
//                     m.channels[0] = draw_from_subs[r_sub()];
//                     m.channels[1] = draw_from_subs[r_sub()];
//                 }
//             }
//         }
//
//         // Calculate the attractors
//         if (create)
//             net->calc_attractors();
//     }
//     return net;
//
// };


} // end namespace
