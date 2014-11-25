// vim: path=.,/usr/include/c++/4.2.1,/usr/include/c++/4.2.1/tr1 fdm=syntax
#include "pubsub2_c.h"
#include <cmath>
#include <stdexcept>
// #include <iostream>

using namespace pubsub2;

signal_t cCisModule::get_site_channel(size_t index) const
{
    // TODO: Fix runtime error?
    if (index >= site_count())
        throw std::runtime_error("illegal site index");

    return _channels[index];
}

signal_t cCisModule::set_site_channel(size_t index, signal_t c)
{
    // TODO: Fix runtime error?
    if (index >= site_count())
        throw std::runtime_error("illegal site index");

    signal_t old = _channels[index];
    _channels[index] = c;
    return old;
}

cGene::cGene(sequence_t seq, signal_t p)
    : sequence(seq)
    , pub(p)
{
}

cGene *cGene::clone()
{
    cGene *g = new cGene(sequence, pub);
    for (auto m : modules)
        g->modules.push_back(m->clone());
    return g;
}

cGene::~cGene()
{
    for (auto m : modules)
        delete m;
}

cNetwork::cNetwork(const cFactory_ptr &f, bool no_ident)
    : pyobject(0)
    , factory(f)
    , parent_identifier(-1)
    , target(-1)
{
    if (!no_ident)
        identifier = factory->get_next_network_ident();
    else
        // A network that is not part of a lineage, for analysis only.
        // We'll call this a "detached" network
        identifier = -1;
}

cNetwork::~cNetwork()
{
    for (auto g : genes)
        delete g;
}

void cNetwork::clone_genes(cGeneVector &gv) const
{
    for (auto g : genes)
        gv.push_back(g->clone());
}

// This is the inner-inner loop!
// TODO: Maybe this could be moved down the the CIS level to prevent constant
// calling of virtual function. It would have to be Factory level call:
// cGeneFactory::cycle(Network &, ChannelState &). But this would mean building 
// the cis action into the factory somehow (or static_cast-ing the CIS which
// would, I guess, be safe.
void cNetwork::cycle(cChannelState &c) const
{
    cChannelState next(c.size());
    for (auto g : genes)
        for (auto cis : g->modules)
            if (cis->is_active(c))
            {
                next.set(g->pub);
                // The gene is active, no use looking further
                break;
            }
    // Update the "return" value.
    c.swap(next);
}

// This is the outer-inner loop, where we find the attractors. 
// TODO: Maybe template-ize this so that it runs faster without constract calls
// to the virtual function. Like this maybe:
// template class attractor_calc<Cis_Type, Factory_Type>
// {
// }
// But it would still need a dynamic runtime selector to call the right one...
// Hmmm.
void cNetwork::calc_attractors()
{
    attractors.clear();
    rates.clear();

    size_t attractor_begins_at;
    bool found;

    // Go through each environment.
    for (auto &env : factory->environments)
    {
        // Set the state to current environment, and set it as the start of the
        // path to the attractor.
        cChannelStateVector path;
        cChannelState current = env;
        path.push_back(current);
        
        for (;;)
        {
            // Update the current state.
            cycle(current);

            // Put back the environment (as this remains constant)
            current |= env;

            // Have we already seen this?
            attractor_begins_at = 0;
            found = false;
            for (cChannelState &prev : path)
            {
                if (prev == current)
                {
                    found = true;
                    break;
                }
                attractor_begins_at++;
            }

            // If we have seen this state, we've found the attractor.
            if (found)
                break;

            // Add the current to our attractor.
            path.push_back(current);
        }
        // Add a new attractor for this environment.
        attractors.push_back(cChannelStateVector());
        cChannelStateVector &this_attr = attractors.back();

        rates.push_back(cRates());
        cRates &this_rate = rates.back();
        for (size_t i=0; i < factory->out_channels; ++i)
            this_rate.push_back(0.0);

        // Copy the part the path that is the attractor, ignoring the transient
        for (size_t copy_at=attractor_begins_at; copy_at < path.size(); ++copy_at)
        {
            cChannelState &c = path[copy_at];
            this_attr.push_back(c);

            // We construct the rates at the same time
            for (size_t i=0; i < factory->out_channels; ++i)
                this_rate[i] += double(c[i + factory->out_range.first]);

        }
        // Now normalise the rates
        for (size_t i=0; i < factory->out_channels; ++i)
            this_rate[i] /= double(this_attr.size());

    }
}

// void cNetwork::make_graph_edges(cEdgeList &edges)
// {
//     // We're building a bipartite graph 
//     // Gene -> Channel
//     // Channel -> Gene
//
//     for g in genes:
//         add g -> pub
//
//         for m in mod:
//             for c in mod:
//                 add c -> g (Label with index = CISMOD/SITE)
//
//     // Now if required FILTER By find_knockouts
// }
//
// void cNetwork::find_knockouts()
// {
//     // return positions GENE -> CISMOD -> SITE
// }

// cNetwork_ptr cNetwork::get_detached_copy() const
// {
//     cNetwork *copy = new cNetwork(factory, true);
//     copy->parent_identifier = identifier;
//     copy->genes = genes;
//     return cNetwork_ptr(copy);
// }

cNetworkAnalysis::cNetworkAnalysis(const cNetwork_ptr &n)
    : original(n)
    , modified(n->factory, true)
{
    // Just copy the genes for now. Don't bother doing anything else.
    original->clone_genes(modified.genes);
}

void cNetworkAnalysis::find_knockouts()
{
    // Incrementally add site knockouts that have no effect. This isn't
    // perfect, but it will do for now.
    for (size_t i=0; i < modified.gene_count(); ++i)
    {
        cGene *g = modified.genes[i];
        for (size_t j=0; j < g->module_count(); ++j)
        {
            cCisModule *m = g->modules[j];
            for (size_t k=0; k < m->site_count(); ++k)
            {
                // Knock this out by setting it to ZERO channel
                signal_t old = m->set_site_channel(k, 0);
                modified.calc_attractors();

                // Did it change? 
                if (modified.attractors == original->attractors)
                    // Record it
                    knockouts.push_back(cSiteIndex(i, j, k));
                else
                    // Reset it
                    m->set_site_channel(k, old);
            }
        }
    }
}

cFactory::cFactory(size_t seed)
    : random_engine(seed)
    , next_network_identifier(0)
    , next_target_identifier(0)
    , gene_count(4)
    , cis_count(3)
{
}

void cFactory::init_environments()
{   
    // Remember: Channel 0 is reserved, so we need to add 1.
    total_channels = cue_channels + reg_channels + out_channels + 1;

    // Number of environments is 2^cue_channels. Use some binary math...
    size_t env_count = 1 << cue_channels;
    for (size_t i = 0; i < env_count; ++i)
    {
        // Shift one, to account for channel 0
        cChannelState c = cChannelState(total_channels, i << 1);
        environments.push_back(c);
    }
}

cTarget::cTarget(cFactory *f)
    : factory(f) 
    , identifier(f->get_next_target_ident())
{
    cRates temp;
    std::fill_n(std::back_inserter(temp), f->out_channels, 0.0);
    std::fill_n(std::back_inserter(optimal_rates), f->environments.size(), temp);
}

// TODO: per env weighting
// TODO: per output weighting
// TODO: Profiling -- make faster?
// TODO: Maybe we should be use multi_array??
double cTarget::assess(const cNetwork &net)
{
    // We've already done it!
    if (net.target == identifier)
        return net.fitness;

    // TODO: check that factories are the same?
    size_t nsize = net.rates.size();
    size_t osize = factory->out_channels;

    if (nsize != optimal_rates.size())
        throw std::runtime_error("optimal rates and networks rates differ in size");

    // We need to score each of the outputs individually
    cRates scores;
    std::fill_n(std::back_inserter(scores), osize, 0.0);

    // The difference from the target reduces the score from a max of 1.0
    for (size_t i = 0; i < nsize; ++i)
        for (size_t j = 0; j < osize; ++j)
            scores[j] += ((1.0 - fabs(net.rates[i][j] - optimal_rates[i][j])) 
                          / double(nsize));

    // Now take the average
    double score = 0.0;
    for (size_t j = 0; j < osize; ++j)
    {
        scores[j] /= double(osize);
        score += scores[j];
    }

    // Must be greater 0 >= n <= 1.0
    net.fitness = score;
    net.target = identifier;
    return score;
}

cGeneFactory::cGeneFactory(cFactory *f, double rate_per_gene_)
    : factory(f)
    , rate_per_gene(rate_per_gene_)
    // Construct a bunch of useful random number generators
    , r_gene(0, f->gene_count-1)
    , r_mod(0, f->cis_count-1)
    , r_sub(f->sub_range.first, f->sub_range.second-1)
{
}

void cGeneFactory::construct_network(cNetwork &network)
{
    for (size_t i=0; i < factory->gene_count; ++i)
    {
        // TODO: genes.emplace_back(...)
        cGene *g = new cGene(i, factory->pub_range.first + i);
        network.genes.push_back(g);

        for (size_t j=0; j < factory->cis_count; ++j)
        {
            g->modules.push_back(construct_cis());
        }
    }

    // Calculate the attractors
    network.calc_attractors();
}

void cGeneFactory::mutate_network(cNetwork_ptr &n, size_t mutations)
{
    // NOTE: This done INPLACE mutation. It should never be called on a network
    // that has already had its attractors calculated! 
    
    // Select the genes that should be mutated
    while (mutations > 0)
    {
        size_t i = r_gene(factory->random_engine);
        mutate_gene(n->genes[i]);
        --mutations;
    }
}

cNetwork_ptr cGeneFactory::copy_and_mutate_network(cNetwork_ptr &n, size_t mutations)
{
    cNetwork_ptr copy(new cNetwork(n->factory));
    // The copy constructor of vector does the hard work here...
    n->clone_genes(copy->genes);
    copy->parent_identifier = n->identifier;

    // Now mutate it and calculate the attractors
    mutate_network(copy, mutations);
    copy->calc_attractors();
    return copy;
}

typedef std::poisson_distribution<> poisson_t;

void cGeneFactory::mutate_collection(cNetworkVector &networks,
                                     cIndexes &mutated)
{
    // How many mutations are we going to have? That depends on the total
    // number of sites that might mutate. Per network, this is rate_per_gene *
    // gene_count. We multiply this by the number of networks in the collection
    // to get the expected number of mutations. Then we generate a number using
    // a poisson distribution.
    // std::cout << "rate per" << rate_per_gene;
    double expected = rate_per_gene * networks.size() * factory->gene_count;
    poisson_t r_pop(expected);
    size_t mutations = r_pop(factory->random_engine);

    // Clear this
    mutated.clear();
    
    // If we're not generating any mutations, let's just bail.
    if (mutations == 0)
        return;


    // Now we need to assign these to individual networks. Only these networks 
    // will change. The rest remain constant.
    randint_t r_network(0, networks.size()-1);
    cIndexes mutes;
    for (size_t i=0; i < mutations; ++i)
        mutes.push_back(r_network(factory->random_engine));

    // Sort them so that any repeated networks are adjacent.
    std::sort(mutes.begin(), mutes.end());

    // We now let the networks figure out how *exactly* they will mutate.  This
    // looks longwinded, but we want to handle cases where a single network is
    // mutated more than once in one step, rather than calling multiple times
    // on the same network, so most of this is putting together multiple
    // mutations into a single call.
    auto it = mutes.begin();
    size_t network_num = *it, count = 1;
    for (;;) 
     {
        // Skip ahead to the next one
        ++it;

        // If we're at the end, or the next one is different, go ahead and
        // apply the mutations ...
        if (it == mutes.end() || *it != network_num)
        {
            networks[network_num] = copy_and_mutate_network(networks[network_num], count);

            // Add to the indexes that changed
            mutated.push_back(network_num);
            if (it == mutes.end())
                break;

            count = 1;
            network_num = *it;
        }
        // ... otherwise we have a repeat of the same network_num. Let's just
        // accumulate this into one call.
        else
        {
            count += 1;
        }
    }
}

void cGeneFactory::mutate_gene(cGene *g)
{
    if (g->modules.size() == 0)
        throw std::runtime_error("no cis modules");

    // First, decide what cis module we're using.
    size_t cis_i = 0; 

    // Account for the fact that there might be only one cis module.
    if (factory->cis_count > 0)
        cis_i = r_mod(factory->random_engine);

    // Grab this and mutate it.
    cCisModule *m = g->modules[cis_i];
    mutate_cis(m);
}


cSelectionModel::cSelectionModel(cFactory_ptr &f)
    : factory(f)
{
}

bool cSelectionModel::select(
    const cNetworkVector &networks, cTarget &target,
    size_t number, cIndexes &selected)
{
    selected.clear();

    std::vector<double> cum_scores;
    cIndexes indexes;
    double score, cum_score = 0.0;
    // First, score everyone
    for (size_t i = 0; i < networks.size(); ++i)
    {
        const cNetwork &net = *(networks[i]);
        score = target.assess(net);

        // Zero fitness has no chance
        if (score <= 0.0)
            continue;

        // TODO: Maybe add some scaling factor to the fitness here
        // ie. score * score, exp(score * y)
        cum_score += score;
        cum_scores.push_back(cum_score);
        indexes.push_back(i);
    }

    // Everyone was crap. 
    if (cum_scores.size() == 0)
        return false;

    // Let's make things a little tidier
    random_engine_t &engine = factory->random_engine;

    // This is standard "roulette selection". The fitness of each individual,
    // once normalised across the population is proportional to their likelihood
    // of being selected for the next generation.
    std::uniform_real_distribution<double> wheel(0.0, cum_score);
    for (size_t i = 0; i < number; ++i)
    {
        double locator = wheel(engine);
        auto it = std::lower_bound(cum_scores.begin(), cum_scores.end(), locator);
        size_t found_index = it - cum_scores.begin();

        selected.push_back(indexes[found_index]);
    }

    return true;
}

void cSelectionModel::copy_using_indexes(
        const cNetworkVector &from, cNetworkVector &to, const cIndexes &selected)
{
    for (auto i : selected)
        // We could check ...
        to.push_back(from[i]);
}

cGeneFactoryLogic2::cGeneFactoryLogic2(cFactory *f, double rate_per_gene_)
    : cGeneFactory(f, rate_per_gene_)
    , r_oper(0, f->operands.size()-1)
{
}

cCisModule *cGeneFactoryLogic2::construct_cis()
{
    cCisModuleLogic2 *cis = new cCisModuleLogic2();
    random_engine_t &re = factory->random_engine;
    cis->op = factory->operands[r_oper(re)];
    cis->channels[0] = r_sub(re);
    cis->channels[1] = r_sub(re);
    return cis;
}

// This is where the action really is.
void cGeneFactoryLogic2::mutate_cis(cCisModule *m)
{
    cCisModuleLogic2 *cis = static_cast<cCisModuleLogic2 *>(m);
    // TODO: This is shite; just a start.
    double p = r_uniform_01(factory->random_engine);
    if (p < .5)
        cis->op = factory->operands[r_oper(factory->random_engine)];
    else if (p < .75)
        cis->channels[0] = r_sub(factory->random_engine);
    else
        cis->channels[1] = r_sub(factory->random_engine);
}

cCisModule *cCisModuleLogic2::clone() const
{
    cCisModuleLogic2 *cis = new cCisModuleLogic2();
    cis->op = op;
    cis->channels[0] = channels[0];
    cis->channels[1] = channels[1];
    return cis;
}
