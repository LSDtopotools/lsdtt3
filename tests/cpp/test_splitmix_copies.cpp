// test_splitmix_copies.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Characterisation test for the splitmix64 copies replaced by
// src/utils/lsdtt_random.hpp (consolidation C15).
//
//   lsdtt_raster_synthetic_surfaces_cpu  private splitmix64(x) -> random::splitmix64
//   lsdtt_hillslope_cellular_cpu         private splitmix64(x) -> random::splitmix64
//                                        (stage 10, C15-H)
//   lsdtt_reference_field_ladder_cpu     inline "state += gamma; mix(state)"
//                                        lambda -> random::splitmix64_stream
//
// lsdtt_random.hpp's own header says migrating a copy "must be shown
// bit-identical"; this is that test. Both removed copies are kept below
// verbatim. The landslide xorshift64* pair (stage 10, L3) has its own test.
//
// Usage: lsdtt-test-splitmix-copies
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "../../src/utils/lsdtt_random.hpp"

using namespace lsdtt;

namespace
{

// src/cpu/lsdtt_raster_synthetic_surfaces_cpu.cpp:115-122 (as of 2026-10-08)
static inline std::uint64_t old_surfaces_splitmix64(std::uint64_t x)
{
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    x =  x ^ (x >> 31);
    return x;
}

// src/cpu/lsdtt_hillslope_cellular_cpu.cpp:55-61 (as of 2026-10-08), stage 10
// (C15-H, Simon's D6 approval); feeds particle_uniform, which stays local.
std::uint64_t old_cellular_splitmix64(std::uint64_t x)
{
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

// src/cpu/lsdtt_reference_field_ladder_cpu.cpp:627-636 (as of 2026-10-08)
std::vector<std::uint64_t> old_ladder_sequence(std::uint64_t seed, int n)
{
    std::uint64_t state = seed ? seed : 1u;
    auto next = [&state]() -> std::uint64_t
    {
        state += 0x9E3779B97F4A7C15ull;
        std::uint64_t z = state;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        return z ^ (z >> 31);
    };
    std::vector<std::uint64_t> out;
    for (int i = 0; i < n; ++i) out.push_back(next());
    return out;
}

// As the ladder now does it.
std::vector<std::uint64_t> new_ladder_sequence(std::uint64_t seed, int n)
{
    random::splitmix64_stream rng(seed ? seed : 1u);
    std::vector<std::uint64_t> out;
    for (int i = 0; i < n; ++i) out.push_back(rng.next());
    return out;
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
    std::cout << "== C15 splitmix64 copies vs lsdtt_random ==\n";

    // Stateless mixer: chained inputs, edge values, and the hash-key pattern
    // synthetic surfaces feeds it.
    bool ok = true;
    std::uint64_t x = 0;
    for (int i = 0; i < 100000; ++i)
    {
        ok = ok && old_surfaces_splitmix64(x) == random::splitmix64(x);
        x = old_surfaces_splitmix64(x) ^ static_cast<std::uint64_t>(i);
    }
    for (std::uint64_t e : {0ULL, 1ULL, ~0ULL, 0x8000000000000000ULL, 0x9E3779B97F4A7C15ULL})
        ok = ok && old_surfaces_splitmix64(e) == random::splitmix64(e);
    check("synthetic_surfaces splitmix64 == random::splitmix64 (100005 inputs)", ok);

    bool ok_c = true;
    std::uint64_t y = 12345;
    for (int i = 0; i < 100000; ++i)
    {
        ok_c = ok_c && old_cellular_splitmix64(y) == random::splitmix64(y);
        y = old_cellular_splitmix64(y) + static_cast<std::uint64_t>(i);
    }
    for (std::uint64_t e : {0ULL, 1ULL, ~0ULL, 0x632be59bd9b4e019ULL, 0x8cb92baa3f3d8dd7ULL, 0xdb4f0b9175ae2165ULL})
        ok_c = ok_c && old_cellular_splitmix64(e) == random::splitmix64(e);
    check("hillslope_cellular splitmix64 == random::splitmix64 (100006 inputs)", ok_c);

    bool ok2 = true;
    for (std::uint64_t seed : {0ULL, 1ULL, 2ULL, 33ULL, 20261008ULL, ~0ULL})
        ok2 = ok2 && old_ladder_sequence(seed, 1000) == new_ladder_sequence(seed, 1000);
    check("ladder bootstrap sequence == splitmix64_stream (first 1000, 6 seeds incl. 0 -> 1)", ok2);

    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
