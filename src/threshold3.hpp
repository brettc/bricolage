#include "core.hpp"
#include <boost/function.hpp>
#include <boost/bind.hpp> 

namespace thresh3 {

enum MutateType { JUMP=0, PROGRESSIVE=1 };

typedef std::function<int()> random_int_t;
inline random_int_t random_int_range(int a, int b, const pubsub2::cWorld_ptr &w)
{ 
    return std::bind(std::uniform_int_distribution<>(a, b-1), std::ref(w->rand));
}

class cFactory : public pubsub2::cFactory
{
public:
    cFactory(const pubsub2::cWorld_ptr &w, size_t module_count, 
                 const MutateType mtype);
    size_t gene_count, module_count;
    MutateType mutate_type;
    pubsub2::cIndexes draw_from_subs; // A vector of sub values that is used to draw from
    random_int_t r_gene, r_module, r_site;
    random_int_t r_binding, r_direction;
    random_int_t r_sub;

    void set_draw_from_subs(const pubsub2::cIndexes &d);

    // Overrides
    pubsub2::cNetwork_ptr construct(bool fill);
    size_t site_count(pubsub2::cNetworkVector &networks);
};

class cCisModule : public pubsub2::cCisModule
{
public:
    cCisModule(const cFactory &c); 

    // Overrides
    size_t site_count() const { return 3; }
    void mutate(const cFactory &c);
    bool is_active(pubsub2::cChannelState const &state) const;

    pubsub2::int_t binding[3];
};


class cGene : public pubsub2::cGene
{
public:
    cGene(pubsub2::sequence_t sequence, pubsub2::signal_t p);
    std::vector<cCisModule> modules;

    // TODO: overrides
    size_t module_count() const { return modules.size(); }
    pubsub2::cCisModule *get_module(size_t i) { return &modules[i]; }
};


class cNetwork : public pubsub2::cNetwork
{
public:
    cNetwork(const pubsub2::cFactory_ptr &c);
    std::vector<cGene> genes;

    // Overrides
    virtual size_t gene_count() const { return genes.size(); }
    virtual cGene *get_gene(size_t i) { return &genes[i]; }
    virtual void mutate(size_t nmutations);
    pubsub2::cNetwork_ptr clone() const;
    void cycle(pubsub2::cChannelState &c) const;
    void cycle_with_intervention(pubsub2::cChannelState &c) const;
};

}
