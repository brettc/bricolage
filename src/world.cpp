#include "core.hpp"
#include <cmath>
#include <stdexcept>
#include <sstream>

using namespace bricolage;

cWorld::cWorld(size_t seed, size_t cue, size_t reg, size_t out)
    // Key stuff
    : next_network_identifier(0)
    , next_target_identifier(0)
    // Channels
    , cue_channels(cue)
    , reg_channels(reg)
    , out_channels(out)
    // Random 
    , rand(seed)
{
    init_channels();
    init_environments();
}

cWorld::~cWorld()
{
}

void cWorld::init_environments()
{
    // Number of environments is 2^cue_channels. Use some binary math...
    size_t env_count = 1 << cue_channels;

    for (size_t i = 0; i < env_count; ++i)
    {
        // Shift one, to account for channel 0
        cChannelState c = cChannelState(channel_count, i << reserved_channels);
        // Turn on bias channel
        c.set(on_channel); 
        environments.push_back(c);
    }
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

cFactory::cFactory(const cWorld_ptr &w)
    : world(w)
{
}

