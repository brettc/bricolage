#include "core.hpp"
#include <boost/function.hpp>
#include <boost/bind.hpp>

namespace thresh3 {

enum MutateType { JUMP=0, PROGRESSIVE=1, JUMP_LAYERED=2 };

typedef std::function<int()> random_int_t;
inline random_int_t random_int_range(int a, int b, const bricolage::cWorld_ptr &w)
{
    return std::bind(std::uniform_int_distribution<>(a, b-1), std::ref(w->rand));
}

class cFactory : public bricolage::cFactory
{
public:
    cFactory(const bricolage::cWorld_ptr &w, size_t module_count,
                 const MutateType mtype);
    size_t gene_count, module_count;
    MutateType mutate_type;
    bricolage::cIndexes draw_from_subs; // A vector of sub values that is used to draw from
    bricolage::cIndexes draw_from_regs; // A vector of only regulatory (no environmental inputs)
    random_int_t r_gene, r_module, r_site;
    random_int_t r_binding, r_direction;
    random_int_t r_sub, r_reg;

    void set_draw_from_subs(const bricolage::cIndexes &d);

    // Overrides
    bricolage::cNetwork_ptr construct(bool fill);
    size_t site_count(bricolage::cNetworkVector &networks);
};

class cGene;
class cCisModule : public bricolage::cCisModule
{
public:
    cCisModule(const cFactory &c);

    // Overrides
    size_t site_count() const { return 3; }
    void mutate(const cFactory &fc, const cGene &gene);
    bool is_active(bricolage::cChannelState const &state) const;

    bricolage::int_t binding[3];
};


class cGene : public bricolage::cGene
{
public:
    cGene(bricolage::sequence_t sequence, bricolage::signal_t p);
    std::vector<cCisModule> modules;

    // TODO: overrides
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
