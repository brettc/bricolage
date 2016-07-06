#pragma once

#include <cstdint>
#include <climits> // For CHAR_BIT

namespace bricolage
{

typedef int_fast8_t signal_t;
typedef int_fast8_t index_t;
typedef signed int int_t;
typedef unsigned int sequence_t;
typedef unsigned long bits_t;

#define MAX_CHANNELS (sizeof(bits_t) * CHAR_BIT)

typedef std::mt19937 random_engine_t;
typedef std::uniform_int_distribution<size_t> randint_t;
typedef std::uniform_real_distribution<> randreal_t;

// Taken from the numpy docs on is_close
const double RELATIVE_TOL = 1e-05;
const double ABSOLUTE_TOL = 1e-08;
inline bool is_close(double a, double b)
{
    return fabs(a - b) <= (ABSOLUTE_TOL + RELATIVE_TOL * fabs(b));
}

inline bool not_zeroish(double a)
{
    return fabs(a) > ABSOLUTE_TOL;
}

// http://stackoverflow.com/questions/1903954/
// is-there-a-standard-sign-function-signum-sgn-in-c-c
template <typename T> inline int c_sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

// python-like comparison function which returns -1, 0, 1
template <typename T> inline int c_cmp(T a, T b)
{
    return c_sgn(a - b);
}

} // end namespace bricolage

