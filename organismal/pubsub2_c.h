// vim: path=.,/usr/include/c++/4.2.1,/usr/include/c++/4.2.1/tr1,/usr/local/include fdm=syntax
#pragma once

#include <cstdint>
#include <vector>

// #include <tr1/memory> Another shared_ptr?
#include <boost/shared_ptr.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/multi_array.hpp>
#include <random>

namespace pubsub2
{

// http://stackoverflow.com/questions/1903954/
// is-there-a-standard-sign-function-signum-sgn-in-c-c
template <typename T> inline int c_sgn(T val) 
{
    return (T(0) < val) - (val < T(0));
}

template <typename T> inline int c_cmp(T a, T b) 
{
    return c_sgn(a - b);
}

typedef int_fast8_t signal_t;
typedef int_fast16_t operand_t;
typedef signed int int_t;
typedef unsigned int sequence_t;

// TODO: Maybe change dynamic_bitset for something a little lighter (Do we really
// need anything bigger than 64bits?)
typedef boost::dynamic_bitset<size_t> cChannelState;
typedef std::vector<cChannelState> cChannelStateVector;
typedef std::vector<cChannelStateVector> cAttractors;
typedef std::vector<size_t> cIndexes;
typedef std::vector<double> cRates;
typedef std::vector<cRates> cRatesVector;

inline int bitset_cmp(cChannelState &a, cChannelState &b)
{
    if (a < b) return -1;
    if (a == b) return 0;
    return 1;
}

class cCisModule
{
public:
    signal_t get_site_channel(size_t index) const;
    signal_t set_site_channel(size_t index, signal_t channel);

    virtual size_t site_count() const = 0;
    virtual bool is_active(cChannelState const &state) const = 0;
    virtual cCisModule* clone() const = 0;
    virtual ~cCisModule() {}
protected:
    signal_t *_channels;
};

typedef std::vector<cCisModule *> cCisModules;


struct cGene
{
    cGene(sequence_t sequence, signal_t p);
    ~cGene();
    cGene *clone();

    size_t module_count() { return modules.size(); }

    sequence_t sequence;
    cCisModules modules;
    signal_t pub;
};

typedef std::vector<cGene *> cGeneVector;
typedef std::vector<operand_t> cOperands;

class cNetwork;
typedef boost::shared_ptr<cNetwork> cNetwork_ptr;
typedef std::vector<cNetwork_ptr> cNetworkVector;

typedef std::mt19937 random_engine_t;
typedef std::uniform_int_distribution<size_t> randint_t;
typedef std::uniform_real_distribution<> randreal_t;

struct cFactory
{
    cFactory(size_t seed);

    random_engine_t random_engine;
    sequence_t next_network_identifier, next_target_identifier;

    sequence_t get_next_network_ident() { return next_network_identifier++; }
    sequence_t get_next_target_ident() { return next_target_identifier++; }

    double get_random_double(double low, double high) {
        return std::uniform_real_distribution<double>(low, high)(random_engine); }
    double get_random_int(int low, int high) {
        return std::uniform_int_distribution<int>(low, high)(random_engine); }

    size_t gene_count, cis_count;
    size_t reg_gene_count;
    size_t cue_channels, reg_channels, out_channels;
    std::pair<size_t, size_t> sub_range;
    std::pair<size_t, size_t> pub_range;
    std::pair<size_t, size_t> out_range;

    cOperands operands;

    // Once all the values are in place we call this
    void init_environments();
    cChannelStateVector environments;
    size_t total_channels;
};

typedef boost::shared_ptr<cFactory> cFactory_ptr;

struct cSiteIndex 
{
    cSiteIndex(int_t g=0, int_t c=0, int_t s=0)
        : gene(g), cismod(c), site(s) {}
    int_t gene, cismod, site;
};
typedef std::vector<cSiteIndex> cSiteLocations;

// struct cSiteLocator
// {
//     cSiteLocator(const cFactory_ptr &f);
//     bool next_site(Network_ptr, cSiteIndex &x, bool init=false);
// };

class cNetwork
{
public:
    cNetwork(const cFactory_ptr &f, bool no_ident=false);
    ~cNetwork();

    size_t gene_count() { return genes.size(); }
    void cycle(cChannelState &c) const;
    void calc_attractors();
    void clone_genes(cGeneVector &gv) const;
    bool is_detached() const { return identifier < 0; }
    cNetwork_ptr get_detached_copy() const;

    mutable void *pyobject;
    
    cFactory_ptr factory;
    int_t identifier, parent_identifier;
    cGeneVector genes;

    // Calculated attractor and rates
    cAttractors attractors;
    cRatesVector rates;

    // Record the fitness and the target against which it was calculated. 
    // This means we don't need to recalc the fitness if the target has not
    // changed.
    mutable int_t target;
    mutable double fitness;

private:
    // Don't allow copying
    cNetwork(const cNetwork &);

};

// cNetwork_ptr get_detached_copy(cNetwork_ptr original);
// inline cNetwork_ptr get_detached_copy(cNetwork_ptr original)
// {
//     cNetwork *copy = new cNetwork(original->factory, true);
//     copy->parent_identifier = original->identifier;
//     copy->genes = original->genes;
//
//     return cNetwork_ptr(copy);
// }

struct cNetworkAnalysis
{
    cNetworkAnalysis(const cNetwork_ptr &n);
    cNetwork_ptr original;
    cNetwork modified;
    cSiteLocations knockouts;

    // Go through and change everything.
    void find_knockouts();
};

struct cTarget
{
    cTarget(cFactory *factory);
    cFactory *factory;
    int_t identifier;
    cRatesVector optimal_rates;

    // TODO: per env weighting
    // TODO: per output weighting
    double assess(const cNetwork &net);

};

// TODO: Eventually this will be a base class, then have different types of
// generators.
struct cGeneFactory
{
    cGeneFactory(cFactory *f, double rate_per_gene_);
    virtual ~cGeneFactory() {};
    cFactory *factory;
    double rate_per_gene;

    // Not sure whether this is worth it, but still...
    randint_t r_gene, r_mod, r_sub;

    // default constructed [0, 1]
    randreal_t r_uniform_01;

    // This is where the constructors and mutators for networks live
    virtual cCisModule *construct_cis() = 0;
    void construct_network(cNetwork &network);

    // Mutates the gene in place
    virtual void mutate_cis(cCisModule *m) = 0;
    void mutate_gene(cGene *g);
    void mutate_network(cNetwork_ptr &n, size_t mutations);

    cNetwork_ptr copy_and_mutate_network(cNetwork_ptr &n, size_t mutations);
    void mutate_collection(cNetworkVector &networks, cIndexes &mutated);
};

struct cSelectionModel
{
    cSelectionModel(cFactory_ptr &factory);
    cFactory_ptr factory;

    bool select(
        const cNetworkVector &networks, cTarget &target, 
        size_t number, cIndexes &selected);

    void copy_using_indexes(
        const cNetworkVector &from, cNetworkVector &to, const cIndexes &selected);
};


struct cGeneFactoryLogic2 : public cGeneFactory
{
    cGeneFactoryLogic2(cFactory *f, double rate_per_gene_);
    randint_t r_oper;

    // Overrides
    cCisModule *construct_cis();
    void mutate_cis(cCisModule *m);
};

class cCisModuleLogic2 : public cCisModule
{
public:
    cCisModuleLogic2() { _channels = channels; }

    size_t site_count() const { return 2; }
    virtual cCisModule* clone() const;
    void mutate();

    // Inline this stuff. It won't change.
    inline bool test(unsigned int a, unsigned int b) const 
    { 
        // Note: C++ standard guarantees integral conversion from bool results
        // in 0 or 1.
        return op & (8 >> ((a << 1) | b)); 
    }

    inline bool is_active(cChannelState const &state) const 
    {
        return test(state.test(channels[0]), state.test(channels[1]));
    }
// protected:
    // Default constructor is fine
    operand_t op;
    signal_t channels[2];

    friend cGeneFactoryLogic2;
};

// struct cSimpleMutationModel : public cMutationModel
// {
//
// };
//

} // end namespace pubsub2
