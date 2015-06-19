#pragma once

#include <functional>
#include <map>
#include "core.hpp"

namespace logic2 {

typedef int_fast16_t operand_t;
typedef std::vector<operand_t> cOperands;
typedef std::pair<bricolage::signal_t, bricolage::signal_t> signal_pair_t;
// typedef std::map<signal_pair_t, operand_t> binding_map_t;

typedef std::function<int()> random_int_t;
inline random_int_t random_int_range(int a, int b, const bricolage::cWorld_ptr &w)
{ 
    return std::bind(std::uniform_int_distribution<>(a, b-1), std::ref(w->rand));
}

struct cFactory : public bricolage::cFactory
{
    cFactory(const bricolage::cWorld_ptr &w, size_t module_count, 
                 const cOperands &ops);
    size_t gene_count, module_count;
    cOperands operands;
    random_int_t r_gene, r_module, r_operand, r_site, r_input;
    // binding_map_t bindings;

    // Overrides
    bricolage::cNetwork_ptr construct(bool fill);
    size_t site_count(bricolage::cNetworkVector &networks);

};

class cCisModule : public bricolage::cCisModule
{
public:
    cCisModule(const cFactory &c);

    // Overrides
    size_t site_count() const { return 2; }
    void mutate(const cFactory &c);

    // Inline this stuff. It won't change.
    inline bool test(unsigned int a, unsigned int b) const 
    { 
        // Note: C++ standard guarantees that integral conversion from bool
        // will result in 0 or 1.
        return op & (8 >> ((a << 1) | b)); 
    }

    inline bool is_active(bricolage::cChannelState const &state) const 
    {
        return test(state.test(channels[0]), state.test(channels[1]));
    }

    operand_t op;
};

class cGene : public bricolage::cGene
{
public:
    cGene(bricolage::sequence_t sequence, bricolage::signal_t p);
    std::vector<cCisModule> modules;

    // Overrides
    size_t module_count() const { return modules.size(); }
    bricolage::cCisModule *get_module(size_t i) { return &modules[i]; }
};


class cNetwork : public bricolage::cNetwork
{
public:
    cNetwork(const bricolage::cFactory_ptr &c);
    std::vector<cGene> genes;

    // Overrides
    virtual size_t gene_count() const { return genes.size(); }
    virtual cGene *get_gene(size_t i) { return &genes[i]; }
    virtual void mutate(size_t nmutations);
    bricolage::cNetwork_ptr clone() const;
    void cycle(bricolage::cChannelState &c) const;
    void cycle_with_intervention(bricolage::cChannelState &c) const;
};


}
