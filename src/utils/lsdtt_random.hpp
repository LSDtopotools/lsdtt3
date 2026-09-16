// lsdtt_random.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Fully specified random numbers: splitmix64.
//
// Why this exists
// ---------------
// The C++ standard fixes the output of engines such as std::mt19937, but NOT of
// std::uniform_int_distribution, std::uniform_real_distribution or
// std::shuffle: libstdc++, libc++ and MSVC produce different values from the
// same seed. Anything that must give identical results on Linux, macOS and
// native Windows (for example the QGIS backend) must not use them. Everything
// here is written out in full, so a seed produces the same sequence on every
// compiler and standard library.
//
// splitmix64 is Steele, Lea & Flood (2014), "Fast splittable pseudorandom number
// generators", with the constants of Vigna's public-domain reference
// implementation (https://prng.di.unimi.it/splitmix64.c).
//
// Relationship to older copies
// ----------------------------
// splitmix64() below is the same function as the static copies in
// lsdtt_hillslope_cellular_cpu.cpp and lsdtt_raster_synthetic_surfaces_cpu.cpp,
// and splitmix64_stream::next() gives the same sequence as the inline generator
// in lsdtt_reference_field_ladder_cpu.cpp. Those copies are left in place;
// migrating them is a separate change that must be shown bit-identical.
//
// Header-only; no dependencies beyond the standard library.
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace lsdtt {
namespace random {

inline constexpr std::uint64_t k_splitmix64_gamma = 0x9E3779B97F4A7C15ULL;

// -----------------------------------------------------------------------------
// splitmix64
//
// Stateless mixing function: adds the golden-ratio increment to x and applies
// the splitmix64 finaliser. Useful as a hash of (seed, index) when values must
// not depend on evaluation order.
// -----------------------------------------------------------------------------
constexpr std::uint64_t splitmix64(std::uint64_t x) noexcept
{
    x += k_splitmix64_gamma;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

// -----------------------------------------------------------------------------
// splitmix64_stream
//
// Sequential generator. The k-th output (k = 1, 2, ...) is
// splitmix64(seed + (k-1) * gamma), which is the sequence of the reference
// implementation started from state = seed.
//
//   next()        - next 64-bit value
//   bounded(n)    - uniform integer in [0, n); n must be > 0. Exactly uniform:
//                   draws below 2^64 mod n are rejected, then the remainder is
//                   taken. Uses only 64-bit arithmetic (no __int128), so it
//                   compiles identically on MSVC.
//   uniform01()   - double in [0, 1) from the top 53 bits
// -----------------------------------------------------------------------------
class splitmix64_stream
{
public:
    explicit constexpr splitmix64_stream(std::uint64_t seed) noexcept
        : state_(seed)
    {
    }

    constexpr std::uint64_t next() noexcept
    {
        const std::uint64_t out = splitmix64(state_);
        state_ += k_splitmix64_gamma;
        return out;
    }

    constexpr std::uint64_t bounded(std::uint64_t n) noexcept
    {
        // 2^64 mod n, computed without 128-bit arithmetic. Draws in
        // [0, threshold) are rejected so every residue is equally likely.
        const std::uint64_t threshold = (0ULL - n) % n;
        for (;;)
        {
            const std::uint64_t x = next();
            if (x >= threshold)
            {
                return x % n;
            }
        }
    }

    double uniform01() noexcept
    {
        return static_cast<double>(next() >> 11) * (1.0 / 9007199254740992.0);
    }

    constexpr std::uint64_t state() const noexcept { return state_; }

private:
    std::uint64_t state_;
};

// -----------------------------------------------------------------------------
// shuffle
//
// Fisher-Yates shuffle driven by splitmix64_stream::bounded. Deterministic for a
// given stream state on every platform (unlike std::shuffle).
// -----------------------------------------------------------------------------
template <typename T>
void shuffle(std::vector<T>& values, splitmix64_stream& rng)
{
    const std::size_t n = values.size();
    if (n < 2)
    {
        return;
    }
    for (std::size_t i = n - 1; i > 0; --i)
    {
        const std::size_t j = static_cast<std::size_t>(
            rng.bounded(static_cast<std::uint64_t>(i) + 1));
        using std::swap;
        swap(values[i], values[j]);
    }
}

} // namespace random
} // namespace lsdtt
