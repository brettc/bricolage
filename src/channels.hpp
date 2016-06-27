#pragma once

#include "defines.hpp"

namespace bricolage
{

struct cChannels
{
    bits_t bits;

    cChannels()
        : bits(0)
    {
    }

    cChannels(cChannels const &other)
        : bits(other.bits)
    {
    }

    void operator=(const cChannels &other)
    {
        bits = other.bits;
    }

    bool operator<(const cChannels &other) const
    {
        return bits < other.bits;
    }

    bool operator==(const cChannels &other) const
    {
        return bits == other.bits;
    }

    index_t max() 
    {
        return MAX_CHANNELS;
    }

    std::string to_string(index_t with_size);

    void _check_size(index_t sz) const
    {
        if (sz >= MAX_CHANNELS)
            throw std::runtime_error("channel size is too big");
    }

    void _check_index(index_t i, index_t size) const
    {
        _check_size(size);
        if (i >= size)
            throw std::runtime_error("channel index is too big");
    }

    void set(index_t i, index_t size)
    {
        _check_index(i, size);
        unchecked_set(i);
    }

    void clear(index_t i, index_t size)
    {
        _check_index(i, size);
        unchecked_clear(i);
    }

    void flip(index_t i, index_t size)
    {
        _check_index(i, size);
        unchecked_flip(i);
    }

    bool test(index_t i, index_t size) const
    {
        _check_index(i, size);
        return unchecked_test(i);
    }

    // The ones to use if you know what you're doing
    void unchecked_set(index_t i)
    {
        bits |= 1 << i;
    }

    void unchecked_clear(index_t i)
    {
        bits &= ~(1 << i);
    }

    void unchecked_flip(index_t i)
    {
        bits ^= (1 << i);
    }

    bool unchecked_test(index_t i) const
    {
        return bits & (1 << i);
    }

    void unchecked_union(const cChannels &other)
    {
        bits |= other.bits;
    }

    void unchecked_intersection(const cChannels &other)
    {
        bits &= other.bits;
    }
};

inline int bitset_cmp(cChannels &a, cChannels &b)
{
    if (a < b) return -1;
    if (a == b) return 0;
    return 1;
}


} // end namespace bricolage
