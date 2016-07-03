#pragma once

#include <vector>
#include <set>
#include <map>
#include <memory>

#include <boost/multi_array.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <tuple>
#include <random>

#include "defines.hpp"
#include "channels.hpp"

namespace bricolage
{

typedef std::vector<cChannels> cAttractor;
typedef std::vector<bits_t> cAttractorBits;

typedef std::vector<cAttractor> cAttractorSet;
typedef std::vector<cAttractorBits> cAttractorSetBits;

typedef std::vector<size_t> cIndexes;
typedef std::vector<double> cRates;
typedef std::vector<cRates> cRatesVector;

typedef std::map<cChannels, cRates> cChannelsRatesMap;
typedef std::map<bits_t, cRates> cChannelsRatesMapBits;

// TODO: Do these need implementation? I forget my C++.
const size_t reserved_channels = 2;
const size_t off_channel = 0;
const size_t on_channel = 1;

// Set this globally
const size_t MAX_CHANNELS_PER_CRM = 4;

// Allows us to intervene on the state of genes and modules, forcing them on or
// off to ascertain what causal role they have.
enum InterventionState
{
    INTERVENE_NONE = 0,
    INTERVENE_ON = 1,
    INTERVENE_OFF = 2,
};

enum InputType
{
    INPUT_CONSTANT = 0,
    INPUT_PULSE = 1,
};


class cCisModule
{
public:
    cCisModule() : intervene(INTERVENE_NONE) {}
    virtual ~cCisModule() {}
    bricolage::signal_t get_site(size_t i) const { return channels[i]; }
    bricolage::signal_t set_site(size_t i, bricolage::signal_t c)
    {
        bricolage::signal_t old = channels[i];
        channels[i] = c;
        return old;
    }
    // Defines how many channels you'll actually use.
    virtual size_t site_count() const = 0;
    InterventionState intervene;
    signal_t channels[MAX_CHANNELS_PER_CRM];
};

struct cGene
{
    cGene(sequence_t sequence, signal_t p);
    virtual ~cGene() {}

    virtual size_t module_count() const=0;
    virtual cCisModule *get_module(size_t i)=0;

    sequence_t sequence;
    signal_t pub;
    InterventionState intervene;
};


class cNetwork;
typedef std::shared_ptr<cNetwork> cNetwork_ptr;
typedef std::vector<cNetwork_ptr> cNetworkVector;

// typedef std::shared_ptr<const cNetwork> cConstNetwork_ptr;
// typedef std::vector<cConstNetwork_ptr> cNetworkVector;


class cFactory;

class cWorld //: std::enable_shared_from_this<cWorld>
{
public:
    cWorld(size_t seed, size_t cue, size_t reg, size_t out);
    ~cWorld();

    sequence_t get_next_network_ident() { return next_network_identifier++; }
    sequence_t get_next_target_ident() { return next_target_identifier++; }

    // Counters
    sequence_t next_network_identifier, next_target_identifier;

    // Channel definitions
    size_t cue_channels, reg_channels, out_channels;
    size_t channel_count;

    size_t gene_count;

    std::pair<size_t, size_t> cue_range;
    std::pair<size_t, size_t> reg_range;
    std::pair<size_t, size_t> out_range;

    std::pair<size_t, size_t> sub_range;
    std::pair<size_t, size_t> pub_range;

    // Once all the values are in place we call this
    cAttractor environments;

    // Randomising stuff
    random_engine_t rand;

    // How the input is handled
    InputType input_type;

    std::string get_random_state();
    void set_random_state(const std::string &s);

private:
    void init_channels();
    void init_environments();
public:
    // Extra cruft down here-----
    // Mainly just for testing
    double get_random_double(double low, double high) {
        return std::uniform_real_distribution<double>(low, high)(rand); }
    double get_random_int(int low, int high) {
        return std::uniform_int_distribution<int>(low, high)(rand); }

};
typedef std::shared_ptr<cWorld> cWorld_ptr;


// Helper function for generating random ranges
typedef std::function<int()> random_int_t;
inline random_int_t random_int_range(int a, int b, const bricolage::cWorld_ptr &w)
{
    return std::bind(std::uniform_int_distribution<>(a, b-1), std::ref(w->rand));
}

// Base class for overriding world operations
class cFactory : public std::enable_shared_from_this<cFactory>
{
public:
    cFactory(const cWorld_ptr &w, size_t cc);
    virtual ~cFactory() {};

    cWorld_ptr world;
    size_t gene_count, module_count;
    bricolage::random_int_t r_gene, r_regulatory, r_module;

    // This is where the factorys and mutators for networks live
    virtual cNetwork_ptr construct(bool fill)=0;
    virtual size_t site_count(cNetworkVector &networks)=0;

    cNetwork_ptr clone_and_mutate_network(
        cNetwork_ptr &n, size_t mutations, size_t dups, int_t generation);

};

typedef std::shared_ptr<cFactory> cFactory_ptr;

class cNetwork
{
public:
    cNetwork(const cFactory_ptr &c);
    virtual ~cNetwork() {}

    virtual cNetwork_ptr clone() const=0;
    virtual void mutate(size_t nmutations)=0;
    virtual void duplicate(size_t ndups)=0;
    virtual size_t gene_count() const=0;
    virtual cGene *get_gene(size_t i)=0;

    virtual void cycle(cChannels &c) const=0;
    virtual void cycle_with_intervention(cChannels &c) const=0;

    void _calc_attractors(bool intervention);
    void calc_attractors() { _calc_attractors(false); }
    void calc_attractors_with_intervention() { _calc_attractors(true); }

    void get_rates(const cChannels &initial, cRates &rates, bool use_cache) const;
    void clear_rate_cache() const { cached_mappings.clear(); }
    void stabilise(const cChannels &initial,
                   bool intervention,
                   cAttractor &attractor_,
                   cAttractor &transient_,
                   cRates &rates_) const;

    double attractor_robustness() const;

    cFactory_ptr factory;
    cWorld_ptr world;
    int_t identifier, parent_identifier;

    // Optional -- the generation that this was created (default: 0)
    int_t generation;

    // Calculated attractor and rates
    cAttractorSet attractors;
    cAttractorSet transients;
    cRatesVector rates;

    // This is a cached map of all computed channels and the resulting rates
    // they lead to (anything in a transient or attractor -> rates).
    mutable cChannelsRatesMap cached_mappings;

    // Record the fitness and the target against which it was calculated.
    // This means we don't need to recalc the fitness if the target has not
    // changed (ideally). 
    mutable int_t target;
    mutable double fitness;

    // This is where we store the Python object for easy recall. This ensures
    // some degree of python object consistence.
    mutable void *pyobject;

private:
    // Don't allow copying
    cNetwork(const cNetwork &);
};


struct cBaseTarget;
struct cSelectionModel;

class cPopulation
{
public:
    cPopulation(const cFactory_ptr &c, size_t n);

    size_t mutate(double site_rate, double dup_rate, int_t generation);
    void assess(const cBaseTarget &target);
    bool select(const cSelectionModel &sm, size_t size);

    std::pair<double, double> worst_and_best() const;
    void best_indexes(cIndexes &best) const;
    size_t get_generation() const { return generation; }

    cFactory_ptr factory;
    cWorld_ptr world;

    sequence_t next_id;
    size_t generation;

    cIndexes selected;
    cIndexes mutated;

    cNetworkVector networks;
    cRates fitnesses;
};



} // end namespace bricolage
// vim: path=.,/usr/include/c++/4.2.1,/usr/include/c++/4.2.1/tr1,/usr/local/include fdm=syntax
