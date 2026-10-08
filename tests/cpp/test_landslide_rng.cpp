// test_landslide_rng.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Characterisation test for random::xorshift64star_uniform01
// (src/utils/lsdtt_random.hpp; consolidation L3, Simon's D6 approval).
//
// lsdtt_landslide_event_cpu and lsdtt_landslide_placement_cpu each had an
// identical xorshift64* next_uniform(state). Both are kept below verbatim and
// compared with the shared function: the returned doubles (by bit pattern) and
// the state left behind, for the first 10000 draws from several seeds,
// including the seeds the landslide tests use.
//
// Usage: lsdtt-test-landslide-rng
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>

#include "../../src/utils/lsdtt_random.hpp"

using namespace lsdtt;

namespace
{

// src/cpu/lsdtt_landslide_event_cpu.cpp:26-33 (as of 2026-10-08)
inline double event_next_uniform(std::uint64_t& state)
{
    state ^= state >> 12;
    state ^= state << 25;
    state ^= state >> 27;
    const std::uint64_t v = state * 2685821657736338717ULL;
    return static_cast<double>(v >> 11) * (1.0 / 9007199254740992.0);
}

// src/cpu/lsdtt_landslide_placement_cpu.cpp:28-36 (as of 2026-10-08)
inline double placement_next_uniform(std::uint64_t& state)
{
    state ^= state >> 12;
    state ^= state << 25;
    state ^= state >> 27;
    const std::uint64_t v = state * 2685821657736338717ULL;
    // 53 significant bits, in [0,1)
    return static_cast<double>(v >> 11) * (1.0 / 9007199254740992.0);
}

bool same_bits(double a, double b)
{
    std::uint64_t x, y;
    std::memcpy(&x, &a, sizeof x);
    std::memcpy(&y, &b, sizeof y);
    return x == y;
}

int g_failures = 0;

void check(const std::string& what, bool ok)
{
    std::cout << (ok ? "PASS " : "FAIL ") << what << "\n";
    if (!ok) ++g_failures;
}

} // namespace

int main()
{
    std::cout << "== L3 landslide xorshift64* copies vs lsdtt_random ==\n";
    bool ok_event = true, ok_place = true;
    for (std::uint64_t seed : {1ULL, 7ULL, 42ULL, 12345ULL, 20261008ULL, 0x9E3779B97F4A7C15ULL, ~0ULL})
    {
        std::uint64_t se = seed, sp = seed, sn1 = seed, sn2 = seed;
        for (int i = 0; i < 10000; ++i)
        {
            const double e = event_next_uniform(se);
            const double p = placement_next_uniform(sp);
            const double n1 = random::xorshift64star_uniform01(sn1);
            const double n2 = random::xorshift64star_uniform01(sn2);
            ok_event = ok_event && same_bits(e, n1) && se == sn1;
            ok_place = ok_place && same_bits(p, n2) && sp == sn2;
        }
    }
    check("landslide_event next_uniform == xorshift64star_uniform01 (values and state, 7 seeds x 10000)", ok_event);
    check("landslide_placement next_uniform == xorshift64star_uniform01 (values and state, 7 seeds x 10000)", ok_place);

    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
