// vim: path=.,/usr/include/c++/4.2.1,/usr/include/c++/4.2.1/tr1 fdm=syntax
#include "pubsub2_c.h"
#include <cmath>
#include <stdexcept>

using namespace pubsub2;

cGene::cGene(sequence_t seq, signal_t p)
    : sequence(seq)
    , pub(p)
{
}

cNetwork::cNetwork(cFactory_ptr &f)
    : pyobject(0)
    , factory(f)
    , parent_identifier(-1)
    , target(-1)
{
    identifier = factory->get_next_network_ident();
}

void cNetwork::cycle(cChannelState &c)
{
    cChannelState next(c.size());
    for (auto &g : genes)
        for (auto &cis : g.modules)
            if (cis.is_active(c))
            {
                next.set(g.pub);
                // The gene is active, no use looking further
                break;
            }

    // Update the "return" value.
    c.swap(next);
}

// This is the inner loop, where we find the attractors. We should tune this.
void cNetwork::calc_attractors()
{
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
            // TODO: This is where we could calculate the rates, and
            // constructed them at the same time.
            cChannelState &c = path[copy_at];
            this_attr.push_back(c);

            for (size_t i=0; i < factory->out_channels; ++i)
                this_rate[i] += double(c[i + factory->out_range.first]);

        }
        for (size_t i=0; i < factory->out_channels; ++i)
            this_rate[i] /= double(this_attr.size());

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
    total_channels = cue_channels + reg_channels + out_channels;
    // Numbor of environments is 2^cue_channels. Use some binary math...
    size_t env_count = 1 << cue_channels;
    for (size_t i = 0; i < env_count; ++i)
    {
        cChannelState c = cChannelState(total_channels, i);
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

cMutationModel::cMutationModel(cFactory *f, double rate_per_gene_)
    : factory(f)
    , rate_per_gene(rate_per_gene_)
    // Construct a bunch of useful random number generators
    , r_gene(0, f->gene_count-1)
    , r_mod(0, f->cis_count-1)
    , r_sub(f->sub_range.first, f->sub_range.second-1)
    , r_oper(0, f->operands.size()-1)
{
}

void cMutationModel::construct_network(cNetwork &network)
{
    for (size_t i=0; i < factory->gene_count; ++i)
    {
        network.genes.push_back(cGene(i, factory->pub_range.first + i));
        // genes.emplace_back(i, rpub(rng));
        // Get a reference to it
        cGene &g = network.genes.back();
        for (size_t j=0; j < factory->cis_count; ++j)
        {
            g.modules.push_back(cCisModule());
            construct_cis(g.modules.back());
        }
    }

    // Calculate the attractors
    network.calc_attractors();
}

void cMutationModel::construct_cis(cCisModule &cis)
{
    random_engine_t &re = factory->random_engine;
    cis.op = factory->operands[r_oper(re)];
    cis.sub1 = r_sub(re);
    cis.sub2 = r_sub(re);
    cis.silenced = false;
    // cis.silenced = bool(r_uniform_01(re) > p_silenced);
}

void cMutationModel::mutate_network(cNetwork_ptr &n, size_t mutations)
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

cNetwork_ptr cMutationModel::copy_and_mutate_network(cNetwork_ptr &n, size_t mutations)
{
    cNetwork_ptr copy(new cNetwork(n->factory));
    // The copy constructor of vector does the hard work here...
    copy->genes = n->genes;
    copy->parent_identifier = n->identifier;

    // Now mutate it and calculate the attractors
    mutate_network(copy, mutations);
    copy->calc_attractors();
    return copy;
}

typedef std::poisson_distribution<> poisson_t;

void cMutationModel::mutate_collection(cNetworkVector &networks,
                                     cIndexes &mutated)
{
    // How many mutations are we going to have? That depends on the total
    // number of sites that might mutate. Per network, this is rate_per_gene *
    // gene_count. We multiply this by the number of networks in the collection
    // to get the expected number of mutations. Then we generate a number using
    // a poisson distribution.
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

    // Sort them so that repeats are next to one another.
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

void cMutationModel::mutate_gene(cGene &g)
{
    // First, decide what cis module we're using.
    size_t cis_i = 0; 

    // Account for the fact that there might be only one cis module.
    if (factory->cis_count > 0)
        cis_i = r_mod(factory->random_engine);

    // Grab this and mutate it.
    cCisModule &m = g.modules[cis_i];
    mutate_cis(m);
}

// This is where the action really is.
void cMutationModel::mutate_cis(cCisModule &m)
{
    // otherwise, let's change something
    double p = r_uniform_01(factory->random_engine);
    if (p < .5)
        m.op = factory->operands[r_oper(factory->random_engine)];
    else if (p < .75)
        m.sub1 = r_sub(factory->random_engine);
    else
        m.sub2 = r_sub(factory->random_engine);
}

cSelectionModel::cSelectionModel(cFactory_ptr &f):
    factory(f)
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
