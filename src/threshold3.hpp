#include "core.hpp"
#include <boost/function.hpp>
#include <boost/bind.hpp> 

namespace thresh3 {

typedef boost::function<int()> random_int_t;


class cConstructor : public pubsub2::cConstructor
{
public:
    cConstructor(const pubsub2::cWorld_ptr &w, size_t cis_count);

    pubsub2::cNetwork_ptr construct();
    size_t site_count(pubsub2::cNetworkVector &networks);

    // Randomize binding
    size_t gene_count, module_count;
    mutable pubsub2::randint_t r_binding, r_direction, r_site, r_input;
    mutable random_int_t r_gene, r_cis;
};

class cCisModule : public pubsub2::cCisModule
{
public:
    cCisModule(const cConstructor &c); 

    // Overrides
    size_t site_count() const { return 3; }
    pubsub2::signal_t get_site(size_t i) const { return channels[i]; }
    pubsub2::signal_t set_site(size_t i, pubsub2::signal_t c) 
    { 
        pubsub2::signal_t old = channels[i];
        channels[i] = c; 
        return old;
    }

    void mutate(const cConstructor &c);
    bool is_active(pubsub2::cChannelState const &state) const;

    pubsub2::int_t binding[3];
    pubsub2::signal_t channels[3];
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
    cNetwork(const cConstructor &c);
    const cConstructor &constructor;
    std::vector<cGene> genes;

    // Overrides
    virtual size_t gene_count() { return genes.size(); }
    virtual cGene *get_gene(size_t i) { return &genes[i]; }
    virtual void mutate(size_t nmutations);
    pubsub2::cNetwork_ptr clone() const;
    void cycle(pubsub2::cChannelState &c) const;
    void cycle_with_intervention(pubsub2::cChannelState &c) const;
};

}
