#include "core.hpp"
#include <boost/function.hpp>
#include <boost/bind.hpp>

namespace thresh3 {

enum MutateType { JUMP=0, PROGRESSIVE=1 };

class cFactory : public bricolage::cFactory
{
public:
    cFactory(const bricolage::cWorld_ptr &w, size_t module_count,
                 const MutateType mtype);
    MutateType mutate_type;

    bricolage::random_int_t 
        r_site, r_binding, r_direction;

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
    void mutate(const cFactory &fc);
    bool is_active(bricolage::cChannels const &state) const;

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
    virtual void mutate(size_t n_sub, size_t n_pub);
    virtual void duplicate(size_t ndups);
    bricolage::cNetwork_ptr clone() const;
    void cycle(bricolage::cChannels &c) const;
    void cycle_with_intervention(bricolage::cChannels &c) const;
};

}
