// test_flat_resolution.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Unit tests for Barnes et al. (2014b) flat resolution
// (src/cpu/lsdtt_raster_flat_resolution_cpu.hpp).
//
// Hermetic: all inputs are built in-process, nothing is read from disk.
//
// Checks:
//   F1  step_up_float32(z, k) is BITWISE equal to k looped
//       std::nextafter(z, +inf) calls, including the sign of zero, the
//       denormals, the crossing from negative to positive values, and the
//       documented overflow / NaN behaviour.
//
// Usage: lsdtt-test-flat-resolution
// Prints one PASS/FAIL line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "lsdtt_raster_flat_resolution_cpu.hpp"

namespace
{

// =============================================================================
// Tiny check harness (same shape as tests/cpp/test_flow_properties.cpp)
// =============================================================================
int g_n_checks = 0;
int g_n_failed = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!ok && !detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << "\n";
}

std::uint32_t bits_of(float x)
{
    std::uint32_t b;
    std::memcpy(&b, &x, sizeof b);
    return b;
}

float float_of(std::uint32_t b)
{
    float x;
    std::memcpy(&x, &b, sizeof x);
    return x;
}

std::string hex_bits(float x)
{
    std::ostringstream os;
    os << "0x" << std::hex << std::setw(8) << std::setfill('0') << bits_of(x);
    return os.str();
}

// The reference: k successive nextafter(z, +inf) calls.
float looped_nextafter(float z, std::int32_t k)
{
    const float up = std::numeric_limits<float>::infinity();
    for (std::int32_t i = 0; i < k; ++i) z = std::nextafter(z, up);
    return z;
}

// Compare step_up_float32 with the loop for every z in zs and every k in ks.
// Reports one line per k, listing the first mismatch if there is one.
void check_against_loop(const std::string& label,
                        const std::vector<float>& zs,
                        const std::vector<std::int32_t>& ks)
{
    for (std::int32_t k : ks) {
        int n_bad = 0;
        std::string first;
        for (float z : zs) {
            const float got = lsdtt::step_up_float32(z, k);
            const float ref = looped_nextafter(z, k);
            if (bits_of(got) != bits_of(ref)) {
                if (n_bad == 0) {
                    std::ostringstream os;
                    os << "z=" << hex_bits(z) << " got " << hex_bits(got)
                       << " want " << hex_bits(ref);
                    first = os.str();
                }
                ++n_bad;
            }
        }
        std::ostringstream name;
        name << "F1 " << label << " k=" << k << " (" << zs.size()
             << " values) bitwise == looped nextafter";
        check(name.str(), n_bad == 0,
              std::to_string(n_bad) + " mismatches; first: " + first);
    }
}

}  // namespace

int main()
{
    std::cout << "=== F1: step_up_float32 vs looped std::nextafter ===\n";

    const float inf        = std::numeric_limits<float>::infinity();
    const float denorm_min = std::numeric_limits<float>::denorm_min();
    auto ulp_below = [](float x) { return std::nextafter(x, 0.0f); };

    // Named values: magnitudes around powers of two (where the ULP changes),
    // zeros of both signs, denormals, and the normal/denormal boundary.
    const std::vector<float> named = {
        1.0f, ulp_below(2.0f), ulp_below(4096.0f), 1024.0f,
        0.0f, -0.0f, denorm_min, -denorm_min, -3.0f * denorm_min,
        -1.0f, 8848.0f, FLT_MIN, -FLT_MIN, ulp_below(FLT_MIN),
        ulp_below(-1024.0f),            // -1024 + ulp (magnitude just below 1024)
        ulp_below(-1.0f), ulp_below(-2.0f), ulp_below(-4096.0f),
        -1024.0f, -2.0f, 100.0f, -50.0f, 4000.0f};

    // Fixed-seed random values in [-9000, 9000].
    std::mt19937 rng(20141009u);
    std::uniform_real_distribution<float> dist(-9000.0f, 9000.0f);
    std::vector<float> randoms(200);
    for (float& v : randoms) v = dist(rng);

    const std::vector<std::int32_t> small_ks = {0, 1, 2, 3, 1000};
    check_against_loop("named", named, small_ks);
    check_against_loop("random", randoms, small_ks);

    // k = 1<<20 costs ~1M nextafter calls per value, ~220M in all here.
    // MEASURED 2026-10-09: the whole driver ran in 0.76 s (8-core laptop),
    // so every value gets the large step.
    const std::vector<std::int32_t> big_k = {1 << 20};
    check_against_loop("named", named, big_k);
    check_against_loop("random", randoms, big_k);

    // The far end of the range: FLT_MAX, the infinities, and -FLT_MAX. The loop
    // saturates at +inf, as step_up_float32 is documented to.
    const std::vector<float> extremes = {FLT_MAX, ulp_below(FLT_MAX), -FLT_MAX,
                                         inf, -inf};
    check_against_loop("extremes", extremes, {0, 1, 2, 3, 1000});

    std::cout << "\n=== F1: named zero / denormal paths ===\n";

    check("F1 -denorm_min k=1 is -0.0 (bits 0x80000000)",
          bits_of(lsdtt::step_up_float32(-denorm_min, 1)) == 0x80000000u,
          hex_bits(lsdtt::step_up_float32(-denorm_min, 1)));
    check("F1 -denorm_min k=2 is +denorm_min",
          bits_of(lsdtt::step_up_float32(-denorm_min, 2)) == bits_of(denorm_min),
          hex_bits(lsdtt::step_up_float32(-denorm_min, 2)));
    check("F1 -0.0 k=1 is +denorm_min",
          bits_of(lsdtt::step_up_float32(-0.0f, 1)) == bits_of(denorm_min),
          hex_bits(lsdtt::step_up_float32(-0.0f, 1)));
    check("F1 +0.0 k=1 is +denorm_min",
          bits_of(lsdtt::step_up_float32(0.0f, 1)) == bits_of(denorm_min),
          hex_bits(lsdtt::step_up_float32(0.0f, 1)));
    check("F1 -0.0 k=0 keeps bits 0x80000000",
          bits_of(lsdtt::step_up_float32(-0.0f, 0)) == 0x80000000u,
          hex_bits(lsdtt::step_up_float32(-0.0f, 0)));
    {
        // Crossing zero: -3*denorm_min -> -2 -> -1 -> -0.0 -> +1 -> +2 denorm.
        const float z   = -3.0f * denorm_min;
        const float got = lsdtt::step_up_float32(z, 5);
        const float ref = looped_nextafter(z, 5);
        check("F1 -3*denorm_min k=5 crosses zero, matches loop (+2*denorm_min)",
              bits_of(got) == bits_of(ref) && bits_of(got) == 2u,
              "got " + hex_bits(got) + " want " + hex_bits(ref));
    }
    {
        // A wider crossing, from -1 to past +1: k = 2 * bits(1.0f) + 1.
        const std::int32_t k = 2 * static_cast<std::int32_t>(bits_of(1.0f)) + 1;
        const float got = lsdtt::step_up_float32(-1.0f, k);
        check("F1 -1.0 stepped by 2*bits(1.0)+1 is nextafter(1.0, +inf)",
              bits_of(got) == bits_of(std::nextafter(1.0f, inf)),
              "got " + hex_bits(got));
    }

    std::cout << "\n=== F1: documented overflow / NaN behaviour ===\n";

    check("F1 FLT_MAX k=1 is +inf",
          bits_of(lsdtt::step_up_float32(FLT_MAX, 1)) == bits_of(inf));
    check("F1 FLT_MAX k=1<<20 is +inf (saturates, no wrap)",
          bits_of(lsdtt::step_up_float32(FLT_MAX, 1 << 20)) == bits_of(inf));
    check("F1 +inf k=INT32_MAX is +inf",
          bits_of(lsdtt::step_up_float32(inf, std::numeric_limits<std::int32_t>::max()))
              == bits_of(inf));
    check("F1 -inf k=1 is -FLT_MAX",
          bits_of(lsdtt::step_up_float32(-inf, 1)) == bits_of(-FLT_MAX));
    {
        // -FLT_MAX key is -0x7F7FFFFF; + INT32_MAX (0x7FFFFFFF) = 0x00800000,
        // i.e. FLT_MIN.
        const float got = lsdtt::step_up_float32(
            -FLT_MAX, std::numeric_limits<std::int32_t>::max());
        check("F1 -FLT_MAX k=INT32_MAX is FLT_MIN (no int32 overflow)",
              bits_of(got) == bits_of(FLT_MIN), "got " + hex_bits(got));
    }
    {
        const float qnan = std::numeric_limits<float>::quiet_NaN();
        const float snan_payload = float_of(0x7FC01234u);
        const float neg_nan      = float_of(0xFFC00001u);
        check("F1 quiet NaN k=1 unchanged bits",
              bits_of(lsdtt::step_up_float32(qnan, 1)) == bits_of(qnan));
        check("F1 NaN with payload k=1000 unchanged bits",
              bits_of(lsdtt::step_up_float32(snan_payload, 1000))
                  == bits_of(snan_payload));
        check("F1 negative NaN k=1 unchanged bits",
              bits_of(lsdtt::step_up_float32(neg_nan, 1)) == bits_of(neg_nan));
    }

    std::cout << "\n" << g_n_checks << " checks, " << g_n_failed << " failed\n";
    std::cout << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";
    return g_n_failed == 0 ? 0 : 1;
}
