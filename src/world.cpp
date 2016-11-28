#include "core.hpp"
#include <cmath>
#include <stdexcept>
#include <sstream>

using namespace bricolage;

std::string cChannels::to_string(index_t sz)
{
    // Note: This is the OPPOSITE of what might be expected, as the least
    // significant bits are first. But this is the way that it makes sense to
    // write our channels.
    _check_size(sz);
    std::string rep(sz, '0');
    for (size_t i = 0; i < sz; ++i)
        if (test(i, sz))
            rep[i] = '1';

    return rep;
}

cWorld::cWorld(size_t seed, size_t cue, size_t reg, size_t out, size_t reg_gene_count)
    // Key stuff
    : next_network_identifier(0)
    , next_target_identifier(0)
    // Genes
    , reg_gene_count(reg_gene_count)
    // Channels
    , cue_channels(cue)
    , reg_channels(reg)
    , out_channels(out)
    // Random
    , rand(seed)
    , input_type(INPUT_PULSE)
    , pulse_for(1)
{
    if (reg_gene_count > 0 && reg == 0)
        throw std::runtime_error("no regulatory channels!");

    init_channels();
    init_environments();
}

cWorld::~cWorld()
{
}


void cWorld::init_channels()
{
    // Calculate the total number of elements given the overlap
    // Example: Given cue = 2, reg = 2, out = 2
    //             0 = ALWAYS OFF
    //               1 = ALWAYS ON
    // reserved  [ 0 1 ]
    // cues          [ 2 3 ]
    // regs              [ 4 5 ]
    // outs                  [ 6 7 8 ]
    //
    // subs      [ 0 1 2 3 4 5 ]
    // pubs              [ 4 5 6 7 8 ]
    channel_count = cue_channels + reg_channels + out_channels + reserved_channels;

    // These are python-like *ranges*, thus the interval is [first, second) or
    // first <= v < second
    cue_range.first = reserved_channels;
    cue_range.second = cue_range.first + cue_channels;

    reg_range.first = cue_range.second;
    reg_range.second = reg_range.first + reg_channels;

    out_range.first = reg_range.second;
    out_range.second = out_range.first + out_channels;

    sub_range.first = 0;
    sub_range.second = reg_range.second;

    pub_range.first = reg_range.first;
    pub_range.second = out_range.second;
}

void cWorld::init_environments()
{
    // Number of environments is 2^cue_channels. Use some binary math...
    size_t env_count = 1 << cue_channels;

    for (size_t i = 0; i < env_count; ++i)
    {
        cChannels c;
        // Shift to account for reserved channels
        c.bits = i << reserved_channels;
        // Turn on bias channel
        c.unchecked_set(on_channel);
        environments.push_back(c);
    }
}

std::string cWorld::get_random_state()
{
    std::stringstream ostr;
    ostr << rand;
    return ostr.str();
}

void cWorld::set_random_state(const std::string &s)
{
    std::stringstream istr(s);
    istr >> rand;
}

cFactory::cFactory(const cWorld_ptr &w, size_t cc)
    : world(w)
    , gene_count(w->reg_gene_count + w->out_channels)
    , module_count(cc)
    , r_gene(random_int_range(0, gene_count, w))
    , r_reg_gene(random_int_range(0, w->reg_gene_count, w))
    , r_module(random_int_range(0, module_count, w))
{
    cIndexes subs, regs;
    for (size_t i = w->sub_range.first; i < w->sub_range.second; ++i)
        subs.push_back(i);
    set_draw_from_subs(subs);

    for (size_t i = w->reg_range.first; i < w->reg_range.second; ++i)
        regs.push_back(i);
    set_draw_from_regs(regs);
}

void cFactory::set_draw_from_subs(const cIndexes &dsubs)
{
    if (dsubs.size() == 0)
        throw std::runtime_error("VARIABLE mutation requires subs");
    for (auto s: dsubs)
        if (s < world->sub_range.first || s >= world->sub_range.second)
            throw std::runtime_error("VARIABLE mutation has invalid subs");

    draw_from_subs = dsubs;
    r_sub = random_int_range(0, draw_from_subs.size(), world);
}

void cFactory::set_draw_from_regs(const cIndexes &dregs)
{
    // It could be the case that we have NO regulatory genes
    // if (dregs.size() == 0)
    //     throw std::runtime_error("VARIABLE mutation requires regs");
    for (auto s: dregs)
        if (s < world->reg_range.first || s >= world->reg_range.second)
            throw std::runtime_error("VARIABLE mutation has invalid regs");

    draw_from_regs = dregs;
    r_reg = random_int_range(0, draw_from_regs.size(), world);
}
