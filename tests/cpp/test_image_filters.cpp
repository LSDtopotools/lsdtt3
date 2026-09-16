// test_image_filters.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Unit tests for src/utils/lsdtt_random.hpp and src/utils/lsdtt_image_filters.*
// (OpenCV removal plan, Phase 1). No OpenCV needed.
//
// Every expected value comes from something other than the code under test:
//   - published splitmix64 outputs, and an independent Python implementation
//     (for the shuffle and bounded() sequences);
//   - hand-written tables (border mapping, Sobel kernels);
//   - brute-force reference implementations written directly from the
//     definitions (box mean, separable filter, Sobel, Otsu);
//   - exhaustive property checks (the fixed-point box-mean division equals
//     true rounding for every possible sum).
//
//   A  random numbers          E  separable float filter
//   B  border mapping          F  Otsu threshold
//   C  Gaussian kernels        G  Sobel
//   D  box mean                H  Bresenham lines
// -----------------------------------------------------------------------------

#include "../../src/utils/lsdtt_image_filters.hpp"
#include "../../src/utils/lsdtt_random.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using lsdtt::image::border_mode;
namespace img = lsdtt::image;
namespace rnd = lsdtt::random;

namespace {

int g_n_checks = 0;
int g_n_failed = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << "\n";
}

std::string hex(std::uint64_t v)
{
    std::ostringstream s;
    s << "0x" << std::hex << v;
    return s.str();
}

// Deterministic test images from splitmix64 (not from the code under test's
// stream class, so a broken stream cannot hide a broken filter).
std::vector<std::uint8_t> random_u8(int rows, int cols, std::uint64_t seed)
{
    std::vector<std::uint8_t> v(static_cast<std::size_t>(rows) * cols);
    for (std::size_t i = 0; i < v.size(); ++i)
        v[i] = static_cast<std::uint8_t>(rnd::splitmix64(seed * 1000003ULL + i) & 0xFF);
    return v;
}

const border_mode k_borders[] = {border_mode::replicate, border_mode::reflect_101};
const char* border_name(border_mode b) { return b == border_mode::replicate ? "replicate" : "reflect_101"; }

// Border mapping written independently of border_index, by explicit
// reflection steps rather than a period formula.
int ref_border(int i, int n, border_mode b)
{
    if (b == border_mode::replicate) return std::max(0, std::min(n - 1, i));
    if (n == 1) return 0;
    while (i < 0 || i >= n)
    {
        if (i < 0)  i = -i;
        if (i >= n) i = 2 * (n - 1) - i;
    }
    return i;
}


// =============================================================================
// A - random numbers
// =============================================================================
void test_random()
{
    std::cout << "--- A: splitmix64 ---\n";

    // Published reference outputs (Vigna's splitmix64.c, state = seed).
    {
        rnd::splitmix64_stream s0(0);
        const std::uint64_t a = s0.next(), b = s0.next(), c = s0.next();
        check("A1 seed 0 first three outputs match the reference implementation",
              a == 0xe220a8397b1dcdafULL && b == 0x6e789e6aa1b965f4ULL &&
              c == 0x06c45d188009454fULL,
              hex(a) + " " + hex(b) + " " + hex(c));
        rnd::splitmix64_stream s1(1);
        const std::uint64_t d = s1.next(), e = s1.next();
        check("A2 seed 1 first two outputs match the reference implementation",
              d == 0x910a2dec89025cc1ULL && e == 0xbeeb8da1658eec67ULL,
              hex(d) + " " + hex(e));
    }

    // Stateless function equals the stream's first draw, and the stream equals
    // the inline generator in lsdtt_reference_field_ladder_cpu.cpp.
    {
        bool same = true;
        std::uint64_t state = 123456789;
        rnd::splitmix64_stream s(123456789);
        for (int i = 0; i < 1000; ++i)
        {
            state += 0x9E3779B97F4A7C15ull;
            std::uint64_t z = state;
            z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
            z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
            z = z ^ (z >> 31);
            if (s.next() != z) { same = false; break; }
        }
        check("A3 stream sequence equals the reference_field_ladder generator (1000 draws)", same);
        check("A4 splitmix64(seed) equals the first stream output",
              rnd::splitmix64(987) == rnd::splitmix64_stream(987).next());
    }

    // bounded(): values from the independent Python implementation.
    {
        rnd::splitmix64_stream s(7);
        std::vector<std::uint64_t> got;
        for (int i = 0; i < 5; ++i) got.push_back(s.bounded(1000));
        const std::vector<std::uint64_t> want = {487, 804, 346, 203, 674};
        check("A5 bounded(1000) sequence for seed 7 matches independent implementation", got == want);
    }

    // bounded(): range, and every residue reachable, for awkward n.
    {
        bool in_range = true, all_hit = true;
        for (std::uint64_t n : {1ULL, 2ULL, 3ULL, 7ULL, 10ULL, 255ULL, 1000ULL})
        {
            rnd::splitmix64_stream s(n * 31 + 5);
            std::vector<int> hits(static_cast<std::size_t>(n), 0);
            for (std::uint64_t i = 0; i < n * 200; ++i)
            {
                const std::uint64_t v = s.bounded(n);
                if (v >= n) { in_range = false; break; }
                ++hits[static_cast<std::size_t>(v)];
            }
            for (int h : hits) if (h == 0) all_hit = false;
        }
        check("A6 bounded(n) never returns >= n", in_range);
        check("A7 bounded(n) reaches every residue (n in {1,2,3,7,10,255,1000})", all_hit);

        // Rejection threshold: with n = 2^63 + 1 the threshold is 2^63 - 1, so
        // roughly half of all raw draws must be rejected.
        const std::uint64_t n = (1ULL << 63) + 1;
        rnd::splitmix64_stream s(99);
        for (int i = 0; i < 2000; ++i) (void)s.bounded(n);
        // Count raw draws by stepping a second stream until its state matches
        // (state differences wrap modulo 2^64, so they cannot be divided).
        rnd::splitmix64_stream probe(99);
        double draws = 0.0;
        while (probe.state() != s.state() && draws < 1e6) { (void)probe.next(); draws += 1.0; }
        check("A8 bounded(2^63+1) rejects about half the raw draws",
              draws > 3600.0 && draws < 4400.0, "raw draws for 2000 values: " + std::to_string(draws));
    }

    // shuffle: independent Python implementation, and permutation property.
    {
        rnd::splitmix64_stream s(42);
        std::vector<int> v = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
        rnd::shuffle(v, s);
        const std::vector<int> want = {0, 9, 5, 8, 6, 4, 7, 2, 1, 3};
        check("A9 shuffle of 0..9 with seed 42 matches independent implementation", v == want);

        std::vector<int> big(5000);
        for (int i = 0; i < 5000; ++i) big[static_cast<std::size_t>(i)] = i;
        rnd::splitmix64_stream s2(3);
        rnd::shuffle(big, s2);
        std::vector<int> sorted = big;
        std::sort(sorted.begin(), sorted.end());
        bool perm = true;
        for (int i = 0; i < 5000; ++i) if (sorted[static_cast<std::size_t>(i)] != i) perm = false;
        int fixed_points = 0;
        for (int i = 0; i < 5000; ++i) if (big[static_cast<std::size_t>(i)] == i) ++fixed_points;
        check("A10 shuffle is a permutation", perm);
        check("A11 shuffle actually moves elements (fixed points of 5000 < 20)",
              fixed_points < 20, std::to_string(fixed_points) + " fixed points");
    }
}


// =============================================================================
// B - border mapping
// =============================================================================
void test_border()
{
    std::cout << "\n--- B: border_index ---\n";

    // Hand-written tables for n = 5, i = -7..11.
    const std::vector<int> rep  = {0,0,0,0,0,0,0, 0,1,2,3,4, 4,4,4,4,4,4,4};
    const std::vector<int> r101 = {1,2,3,4,3,2,1, 0,1,2,3,4, 3,2,1,0,1,2,3};
    bool ok_rep = true, ok_101 = true;
    for (int i = -7; i <= 11; ++i)
    {
        if (img::border_index(i, 5, border_mode::replicate)   != rep [static_cast<std::size_t>(i + 7)]) ok_rep = false;
        if (img::border_index(i, 5, border_mode::reflect_101) != r101[static_cast<std::size_t>(i + 7)]) ok_101 = false;
    }
    check("B1 replicate matches hand table (n=5, i=-7..11)", ok_rep);
    check("B2 reflect_101 matches hand table (n=5, i=-7..11)", ok_101);

    bool ok_small = true;
    for (int n : {1, 2, 3})
        for (int i = -60; i <= 60; ++i)
            for (border_mode b : k_borders)
                if (img::border_index(i, n, b) != ref_border(i, n, b)) ok_small = false;
    check("B3 matches step-by-step reflection for n in {1,2,3}, |i| <= 60 (kernel wider than image)", ok_small);
}


// =============================================================================
// C - Gaussian kernels
// =============================================================================
void test_gaussian()
{
    std::cout << "\n--- C: Gaussian kernels ---\n";

    auto sum = [](const std::vector<double>& k) { double s = 0; for (double v : k) s += v; return s; };
    auto symmetric = [](const std::vector<double>& k)
    {
        for (std::size_t i = 0; i < k.size() / 2; ++i) if (k[i] != k[k.size() - 1 - i]) return false;
        return true;
    };

    // Fixed tables for sigma <= 0 (OpenCV getGaussianKernelBitExact).
    check("C1 ksize 3 sigma 0 is {0.25, 0.5, 0.25}",
          img::gaussian_kernel_opencv(3, 0.0) == std::vector<double>{0.25, 0.5, 0.25});
    check("C2 ksize 7 sigma 0 is the OpenCV table",
          img::gaussian_kernel_opencv(7, 0.0) ==
          std::vector<double>{0.03125, 0.109375, 0.21875, 0.28125, 0.21875, 0.109375, 0.03125});
    check("C3 ksize 9 sigma -1 uses the table (4,13,30,51,60,...)/256",
          img::gaussian_kernel_opencv(9, -1.0)[4] == 60.0 / 256 &&
          img::gaussian_kernel_opencv(9, -1.0)[0] == 4.0 / 256);

    // Computed kernels: independent formula exp(-d^2 / 2 sigma^2), normalised.
    double worst = 0.0;
    bool sums_ok = true, sym_ok = true;
    for (int ksize : {11, 15, 51})
    {
        for (double sigma : {0.0, 0.8, 2.5})
        {
            const std::vector<double> k = img::gaussian_kernel_opencv(ksize, sigma);
            const double s_x = sigma > 0 ? sigma : 0.3 * ((ksize - 1) * 0.5 - 1) + 0.8;
            std::vector<double> ref(static_cast<std::size_t>(ksize));
            double rs = 0;
            for (int i = 0; i < ksize; ++i)
            {
                const double d = i - (ksize - 1) / 2.0;
                ref[static_cast<std::size_t>(i)] = std::exp(-d * d / (2 * s_x * s_x));
                rs += ref[static_cast<std::size_t>(i)];
            }
            for (int i = 0; i < ksize; ++i)
                worst = std::max(worst, std::fabs(k[static_cast<std::size_t>(i)] - ref[static_cast<std::size_t>(i)] / rs));
            if (std::fabs(sum(k) - 1.0) > 1e-12) sums_ok = false;
            if (!symmetric(k)) sym_ok = false;
        }
    }
    check("C4 computed kernels match exp(-d^2/2s^2)/sum, with s = 0.3((k-1)/2-1)+0.8 when sigma<=0",
          worst < 1e-15, "max abs diff " + std::to_string(worst));
    check("C5 computed kernels sum to 1 within 1e-12", sums_ok);
    check("C6 computed kernels are exactly symmetric", sym_ok);

    // Truncated kernel (moved from the Frangi filter): radius ceil(3 sigma).
    const std::vector<double> t = img::gaussian_kernel_truncated(1.2);
    check("C7 truncated kernel radius is ceil(3 sigma) (sigma 1.2 -> 9 taps)", t.size() == 9);
    check("C8 truncated kernel sums to 1", std::fabs(sum(t) - 1.0) < 1e-12);
    check("C9 truncated kernel has radius >= 1 for tiny sigma", img::gaussian_kernel_truncated(0.01).size() == 3);
    bool threw = false;
    try { (void)img::gaussian_kernel_truncated(0.0); } catch (const std::runtime_error&) { threw = true; }
    check("C10 truncated kernel throws std::runtime_error for sigma 0 (callers catch runtime_error)", threw);
}


// =============================================================================
// D - box mean
// =============================================================================

// Brute force from the definition: neighbourhood sum via border mapping, then
// the rounding rule stated independently as true rounding (odd areas have no
// ties). D2/D3 prove that equals OpenCV's fixed-point formula.
std::vector<std::uint8_t> ref_box_mean(const std::vector<std::uint8_t>& src, int rows, int cols,
                                       int k, border_mode b)
{
    const int r = k / 2;
    std::vector<std::uint8_t> out(src.size());
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < cols; ++x)
        {
            long s = 0;
            for (int dy = -r; dy <= r; ++dy)
                for (int dx = -r; dx <= r; ++dx)
                    s += src[static_cast<std::size_t>(ref_border(y + dy, rows, b)) * cols + ref_border(x + dx, cols, b)];
            const long area = static_cast<long>(k) * k;
            out[static_cast<std::size_t>(y) * cols + x] = static_cast<std::uint8_t>((2 * s + area) / (2 * area));
        }
    return out;
}

void test_box_mean()
{
    std::cout << "\n--- D: box_mean_u8 ---\n";

    // D1: brute force on a non-square image smaller than the largest block.
    {
        const int rows = 23, cols = 37;
        const auto src = random_u8(rows, cols, 1);
        int n_bad_cases = 0;
        std::string first_bad;
        for (border_mode b : k_borders)
            for (int k = 3; k <= 51; k += 2)
            {
                std::vector<std::uint8_t> got;
                img::box_mean_u8(src, rows, cols, k, b, got);
                if (got != ref_box_mean(src, rows, cols, k, b))
                {
                    if (n_bad_cases++ == 0) first_bad = std::string(border_name(b)) + " k=" + std::to_string(k);
                }
            }
        check("D1 equals brute force for every odd block 3..51, both borders, 23x37 random image",
              n_bad_cases == 0, n_bad_cases ? "first mismatch " + first_bad : "");
    }

    // D2/D3: the OpenCV rounding rules equal true rounding for EVERY possible
    // sum, derived here independently of the implementation's constants.
    {
        bool fixed_ok = true;
        std::string where;
        for (int k = 3; k * k <= 256; k += 2)
        {
            const int d = k * k;
            double scalef = static_cast<double>(1 << 23) / d;
            int ds = static_cast<int>(std::floor(scalef));
            scalef -= ds;
            int dd = d / 2;
            if (scalef < 0.5) ++dd; else ++ds;
            for (long s = 0; s <= 255L * d; ++s)
            {
                const long fx = ((s + dd) * static_cast<long>(ds)) >> 23;
                const long tr = (2 * s + d) / (2 * d);
                if (fx != tr) { fixed_ok = false; where = "k=" + std::to_string(k) + " s=" + std::to_string(s); break; }
            }
            if (!fixed_ok) break;
        }
        check("D2 OpenCV fixed-point division ((s+dd)*ds)>>23 == round(s/area) for all s, blocks 3..15",
              fixed_ok, where);

        bool float_ok = true;
        for (int k = 17; k <= 51 && float_ok; k += 2)
        {
            const long area = static_cast<long>(k) * k;
            const double inv = 1.0 / static_cast<double>(area);
            for (long s = 0; s <= 255L * area; ++s)
            {
                const double v = static_cast<double>(s) * inv;
                const long nearest = static_cast<long>(std::nearbyint(v));
                if (nearest != (2 * s + area) / (2 * area)) { float_ok = false; where = "k=" + std::to_string(k); break; }
            }
        }
        check("D3 OpenCV float scale round(s*(1/area)) == round(s/area) for all s, blocks 17..51",
              float_ok, where);
    }

    // D4: constant image stays constant; ksize 1 is identity.
    {
        const int rows = 10, cols = 9;
        std::vector<std::uint8_t> c(static_cast<std::size_t>(rows) * cols, 173), got;
        img::box_mean_u8(c, rows, cols, 51, border_mode::reflect_101, got);
        check("D4 constant image is unchanged by a 51x51 mean", got == c);
        const auto src = random_u8(rows, cols, 9);
        img::box_mean_u8(src, rows, cols, 1, border_mode::replicate, got);
        check("D5 ksize 1 returns the input", got == src);
    }

    // D6: invalid parameters.
    {
        std::vector<std::uint8_t> v(12, 0), out;
        bool even = false, size = false;
        try { img::box_mean_u8(v, 3, 4, 4, border_mode::replicate, out); } catch (const std::invalid_argument&) { even = true; }
        try { img::box_mean_u8(v, 3, 5, 3, border_mode::replicate, out); } catch (const std::invalid_argument&) { size = true; }
        check("D6 even ksize and wrong image size throw std::invalid_argument", even && size);
    }
}


// =============================================================================
// E - separable float filter
// =============================================================================
void test_separable()
{
    std::cout << "\n--- E: separable_filter_f32 ---\n";

    // E1: impulse response shows correlation orientation (no flip).
    {
        const int rows = 7, cols = 7;
        std::vector<float> src(49, 0.0f), dst;
        src[3 * 7 + 3] = 1.0f;
        const std::vector<double> kx = {1.0, 2.0, 3.0}, ky = {10.0, 20.0, 30.0};
        img::separable_filter_f32(src, rows, cols, kx, ky, border_mode::reflect_101, dst);
        // Correlation: out(y,x) = sum ky[j] kx[i] src(y+j-1, x+i-1); the impulse at
        // (3,3) is seen from (y,x) = (3-(j-1), 3-(i-1)), so out(2,2) = ky[2]*kx[2].
        check("E1 correlation orientation: out(2,2) = kx[2]*ky[2] = 90, out(4,4) = kx[0]*ky[0] = 10",
              dst[2 * 7 + 2] == 90.0f && dst[4 * 7 + 4] == 10.0f && dst[3 * 7 + 3] == 40.0f,
              "out(2,2)=" + std::to_string(dst[16]) + " out(4,4)=" + std::to_string(dst[32]));
    }

    // E2: brute-force double convolution with borders.
    {
        const int rows = 19, cols = 26;
        const auto u = random_u8(rows, cols, 5);
        std::vector<float> src(u.begin(), u.end()), dst;
        double worst = 0.0;
        for (border_mode b : k_borders)
            for (int k : {3, 7, 21, 51})
            {
                const std::vector<double> kx = img::gaussian_kernel_opencv(k, 0.0);
                const std::vector<double> ky = img::gaussian_kernel_opencv(k, 1.7);
                img::separable_filter_f32(src, rows, cols, kx, ky, b, dst);
                const int r = k / 2;
                for (int y = 0; y < rows; ++y)
                    for (int x = 0; x < cols; ++x)
                    {
                        double acc = 0.0;
                        for (int j = -r; j <= r; ++j)
                            for (int i = -r; i <= r; ++i)
                                acc += ky[static_cast<std::size_t>(j + r)] * kx[static_cast<std::size_t>(i + r)] *
                                       src[static_cast<std::size_t>(ref_border(y + j, rows, b)) * cols + ref_border(x + i, cols, b)];
                        worst = std::max(worst, std::fabs(acc - dst[static_cast<std::size_t>(y) * cols + x]));
                    }
            }
        check("E2 matches brute-force 2-D convolution within float precision (max |diff| < 1e-3 on 0..255 data)",
              worst < 1e-3, "max |diff| " + std::to_string(worst));
    }
}


// =============================================================================
// E3-E6 - fixed-point 8-bit Gaussian blur
// =============================================================================
void test_fixed_blur()
{
    std::cout << "\n--- E: gaussian_blur_u8 (fixed point) ---\n";

    bool sums = true, symm = true;
    for (int k = 3; k <= 21; k += 2)
        for (double sigma : {0.5, 1.0, 2.3, 5.0})
        {
            const auto kr = img::gaussian_kernel_fixed_8_8(k, sigma);
            long t = 0;
            for (auto v : kr) t += v;
            if (t != 256) sums = false;
            for (std::size_t i = 0; i < kr.size() / 2; ++i) if (kr[i] != kr[kr.size() - 1 - i]) symm = false;
        }
    check("E3 fixed-point kernels sum to exactly 256", sums);
    check("E4 fixed-point kernels are symmetric", symm);

    // Error diffusion, from the definition: for sigma 2.3, ksize 15 the
    // running error must stay within +-0.5 and the centre absorbs the rest.
    {
        const auto kd = img::gaussian_kernel_opencv(15, 2.3);
        const auto kr = img::gaussian_kernel_fixed_8_8(15, 2.3);
        double err = 0.0, worst = 0.0;
        for (int i = 0; i < 7; ++i)
        {
            err += kd[static_cast<std::size_t>(i)] * 256.0 - kr[static_cast<std::size_t>(i)];
            worst = std::max(worst, std::fabs(err));
        }
        check("E5 error-diffused kernel: cumulative rounding error stays within 0.5", worst <= 0.5,
              "max |cumulative error| " + std::to_string(worst));
    }

    // Brute force: exact integer convolution then (v + 2^15) >> 16.
    {
        const int rows = 17, cols = 31;
        const auto src = random_u8(rows, cols, 44);
        int bad = 0;
        for (border_mode b : k_borders)
            for (double sigma : {0.7, 1.5, 3.1})
            {
                const int k = 2 * static_cast<int>(std::ceil(3 * sigma)) + 1;
                const auto kr = img::gaussian_kernel_fixed_8_8(k, sigma);
                std::vector<std::uint8_t> got;
                img::gaussian_blur_u8(src, rows, cols, k, sigma, b, got);
                const int r = k / 2;
                for (int y = 0; y < rows; ++y)
                    for (int x = 0; x < cols; ++x)
                    {
                        std::uint64_t acc = 0;
                        for (int j = -r; j <= r; ++j)
                            for (int i = -r; i <= r; ++i)
                                acc += static_cast<std::uint64_t>(kr[static_cast<std::size_t>(j + r)]) * kr[static_cast<std::size_t>(i + r)] *
                                       src[static_cast<std::size_t>(ref_border(y + j, rows, b)) * cols + ref_border(x + i, cols, b)];
                        if (got[static_cast<std::size_t>(y) * cols + x] != static_cast<std::uint8_t>((acc + 32768) >> 16)) ++bad;
                    }
            }
        check("E6 fixed-point blur equals brute-force 2-D integer convolution (kernel wider than image, both borders)",
              bad == 0, std::to_string(bad) + " pixels differ");
        std::vector<std::uint8_t> c(static_cast<std::size_t>(rows) * cols, 91), out;
        img::gaussian_blur_u8(c, rows, cols, 19, 3.1, border_mode::reflect_101, out);
        check("E7 constant image is unchanged by the fixed-point blur", out == c);
    }
}


// =============================================================================
// F - Otsu
// =============================================================================

// Brute force from the definition: between-class variance for "foreground =
// value > t", evaluated for every t, first strict maximum. Returns -1 if no t
// separates two non-empty classes.
int ref_otsu(const std::array<std::uint64_t, 256>& h, double* best_margin)
{
    double total = 0, sum_all = 0;
    for (int i = 0; i < 256; ++i) { total += static_cast<double>(h[static_cast<std::size_t>(i)]); sum_all += i * static_cast<double>(h[static_cast<std::size_t>(i)]); }
    double best = -1, second = -1;
    int best_t = -1;
    double w0 = 0, s0 = 0;
    for (int t = 0; t < 256; ++t)
    {
        w0 += static_cast<double>(h[static_cast<std::size_t>(t)]);
        s0 += t * static_cast<double>(h[static_cast<std::size_t>(t)]);
        const double w1 = total - w0;
        if (w0 == 0 || w1 == 0) continue;
        const double m0 = s0 / w0, m1 = (sum_all - s0) / w1;
        const double var = (w0 / total) * (w1 / total) * (m0 - m1) * (m0 - m1);
        if (var > best) { second = best; best = var; best_t = t; }
        else if (var > second) second = var;
    }
    if (best_margin) *best_margin = best > 0 ? (best - second) / best : 0;
    return best_t;
}

void test_otsu()
{
    std::cout << "\n--- F: otsu_threshold_u8 ---\n";

    // F1: two separated spikes: every t in [50, 199] gives the same variance, so
    // OpenCV's first-maximum rule returns 50.
    {
        std::array<std::uint64_t, 256> h{};
        h[50] = 1000; h[200] = 3000;
        check("F1 two spikes at 50 and 200: first maximum is 50", img::otsu_threshold_u8(h) == 50,
              "got " + std::to_string(img::otsu_threshold_u8(h)));
    }
    // F2: a single occupied bin cannot be split.
    {
        std::array<std::uint64_t, 256> h{};
        h[128] = 777;
        check("F2 single-valued histogram returns 0", img::otsu_threshold_u8(h) == 0);
        std::array<std::uint64_t, 256> e{};
        check("F3 empty histogram returns 0", img::otsu_threshold_u8(e) == 0);
    }
    // F4: two overlapping Gaussians: threshold between the modes.
    {
        std::array<std::uint64_t, 256> h{};
        for (int i = 0; i < 256; ++i)
            h[static_cast<std::size_t>(i)] = static_cast<std::uint64_t>(
                5000 * std::exp(-std::pow(i - 70, 2) / 400.0) + 3000 * std::exp(-std::pow(i - 180, 2) / 900.0));
        const int t = img::otsu_threshold_u8(h);
        check("F4 bimodal (modes 70 and 180): threshold lies between the modes", t > 70 && t < 180,
              "got " + std::to_string(t));
    }
    // F5: brute force on random histograms where the maximum is unambiguous.
    {
        int compared = 0, mismatches = 0;
        std::string first;
        for (std::uint64_t seed = 0; seed < 400; ++seed)
        {
            std::array<std::uint64_t, 256> h{};
            for (int i = 0; i < 256; ++i)
                h[static_cast<std::size_t>(i)] = rnd::splitmix64(seed * 7919 + static_cast<std::uint64_t>(i)) % 1000;
            double margin = 0;
            const int want = ref_otsu(h, &margin);
            if (margin < 1e-9) continue;
            ++compared;
            const int got = img::otsu_threshold_u8(h);
            if (got != want && mismatches++ == 0)
                first = "seed " + std::to_string(seed) + " got " + std::to_string(got) + " want " + std::to_string(want);
        }
        check("F5 matches brute-force between-class variance on random histograms",
              mismatches == 0 && compared > 300,
              std::to_string(compared) + " compared, " + std::to_string(mismatches) + " mismatches " + first);
    }
    // F7: a class holding less than FLT_EPSILON of the pixels is never split
    // off. h[0] = h[10] = 1 and h[200] = 1e8 - 2: every candidate threshold
    // 0..199 has q1 <= 2e-8 < FLT_EPSILON, so no bin qualifies and OpenCV's
    // rule returns 0. Without the skip, t = 10 would win (larger variance than
    // t = 0). Bin 0 is occupied on purpose: an empty leading bin makes q1 = 0,
    // and the resulting 0/0 would hide a missing skip.
    {
        std::array<std::uint64_t, 256> h{};
        h[0] = 1; h[10] = 1; h[200] = 100000000ULL - 2;
        check("F7 classes smaller than FLT_EPSILON of the pixels are not split off (returns 0)",
              img::otsu_threshold_u8(h) == 0, "got " + std::to_string(img::otsu_threshold_u8(h)));
    }
    // F6: histogram_u8 counts.
    {
        const auto src = random_u8(31, 17, 12);
        const auto h = img::histogram_u8(src);
        std::uint64_t total = 0;
        for (auto v : h) total += v;
        std::uint64_t n200 = 0;
        for (auto v : src) if (v == 200) ++n200;
        check("F6 histogram_u8 totals the pixel count and counts a value correctly",
              total == src.size() && h[200] == n200);
    }
}


// =============================================================================
// G - Sobel
// =============================================================================
void test_sobel()
{
    std::cout << "\n--- G: Sobel ---\n";

    check("G1 ksize 3 kernels {1,2,1} and {-1,0,1}",
          img::sobel_kernel(3, 0) == std::vector<int>{1, 2, 1} &&
          img::sobel_kernel(3, 1) == std::vector<int>{-1, 0, 1});
    check("G2 ksize 5 kernels {1,4,6,4,1} and {-1,-2,0,2,1}",
          img::sobel_kernel(5, 0) == std::vector<int>{1, 4, 6, 4, 1} &&
          img::sobel_kernel(5, 1) == std::vector<int>{-1, -2, 0, 2, 1});
    check("G3 ksize 7 kernels {1,6,15,20,15,6,1} and {-1,-4,-5,0,5,4,1}",
          img::sobel_kernel(7, 0) == std::vector<int>{1, 6, 15, 20, 15, 6, 1} &&
          img::sobel_kernel(7, 1) == std::vector<int>{-1, -4, -5, 0, 5, 4, 1});
    check("G4 ksize 5 second derivative {1,0,-2,0,1}",
          img::sobel_kernel(5, 2) == std::vector<int>{1, 0, -2, 0, 1});

    // G5: brute force with borders, every aperture, both derivatives, both divisors.
    {
        const int rows = 21, cols = 29;
        const auto src = random_u8(rows, cols, 33);
        int bad = 0;
        std::string first;
        for (border_mode b : k_borders)
            for (int k : {3, 5, 7})
                for (int dir = 0; dir < 2; ++dir)
                    for (int divisor : {1, 16})
                    {
                        const int dx = dir == 0 ? 1 : 0, dy = 1 - dx;
                        std::vector<std::int16_t> got;
                        img::sobel_s16(src, rows, cols, dx, dy, k, divisor, b, got);
                        const auto kx = img::sobel_kernel(k, dx), ky = img::sobel_kernel(k, dy);
                        const int r = k / 2;
                        for (int y = 0; y < rows; ++y)
                            for (int x = 0; x < cols; ++x)
                            {
                                long acc = 0;
                                for (int j = -r; j <= r; ++j)
                                    for (int i = -r; i <= r; ++i)
                                        acc += ky[static_cast<std::size_t>(j + r)] * kx[static_cast<std::size_t>(i + r)] *
                                               static_cast<long>(src[static_cast<std::size_t>(ref_border(y + j, rows, b)) * cols + ref_border(x + i, cols, b)]);
                                // Ties to even, written with floating point here.
                                double v = std::nearbyint(static_cast<double>(acc) / divisor);
                                v = std::max(-32768.0, std::min(32767.0, v));
                                if (got[static_cast<std::size_t>(y) * cols + x] != static_cast<std::int16_t>(v) && bad++ == 0)
                                    first = std::string(border_name(b)) + " k=" + std::to_string(k) +
                                            " dx=" + std::to_string(dx) + " div=" + std::to_string(divisor);
                            }
                    }
        check("G5 equals brute force (apertures 3/5/7, dx and dy, divisor 1 and 16, both borders)",
              bad == 0, bad ? std::to_string(bad) + " pixels, first " + first : "");
    }

    // G6: ramp of slope a along x: dx = a * sum(i*kx[i]) * sum(ky) = 8a for ksize 3.
    {
        const int rows = 9, cols = 12;
        std::vector<std::uint8_t> ramp(static_cast<std::size_t>(rows) * cols);
        for (int y = 0; y < rows; ++y) for (int x = 0; x < cols; ++x) ramp[static_cast<std::size_t>(y) * cols + x] = static_cast<std::uint8_t>(3 * x);
        std::vector<std::int16_t> gx, gy;
        img::sobel_s16(ramp, rows, cols, 1, 0, 3, 1, border_mode::replicate, gx);
        img::sobel_s16(ramp, rows, cols, 0, 1, 3, 1, border_mode::replicate, gy);
        bool ok = gx[4 * cols + 5] == 24;
        for (auto v : gy) if (v != 0) ok = false;
        check("G6 ramp slope 3: interior dx = 24, dy = 0 everywhere", ok,
              "dx(4,5)=" + std::to_string(gx[4 * cols + 5]));
    }

    // G7: saturation and tie rounding.
    {
        const int rows = 9, cols = 9;
        std::vector<std::uint8_t> step(81, 0);
        for (int y = 0; y < rows; ++y) for (int x = 5; x < cols; ++x) step[static_cast<std::size_t>(y) * cols + x] = 255;
        std::vector<std::int16_t> g;
        img::sobel_s16(step, rows, cols, 1, 0, 7, 1, border_mode::replicate, g);
        check("G7 ksize 7 on a 0/255 step saturates to 32767", g[4 * cols + 4] == 32767,
              "got " + std::to_string(g[4 * cols + 4]));

        // A single 1 at (4,4), dx, ksize 3: raw values are ky*kx in {+-1, +-2}.
        // With divisor 4: 2/4 = 0.5 -> 0 (even), -2/4 -> 0; 1/4 -> 0.
        // Use divisor 2 via a scaled pixel: value 3 -> raw 6 at the centre row
        // (2*3=6 -> 3), raw 3 at off rows (3/2 = 1.5 -> 2, -1.5 -> -2).
        std::vector<std::uint8_t> dot(81, 0);
        dot[4 * cols + 4] = 3;
        img::sobel_s16(dot, rows, cols, 1, 0, 3, 2, border_mode::replicate, g);
        // out(3,3) = ky[2]*kx[2]*3 = 1*1*3 = 3 -> 3/2 = 1.5 -> 2; out(4,3) = 2*1*3 = 6 -> 3;
        // out(3,5) = 1*(-1)*3 = -3 -> -1.5 -> -2.
        check("G8 divisor rounds ties to even (1.5 -> 2, -1.5 -> -2, 3 -> 3)",
              g[3 * cols + 3] == 2 && g[4 * cols + 3] == 3 && g[3 * cols + 5] == -2,
              std::to_string(g[3 * cols + 3]) + " " + std::to_string(g[4 * cols + 3]) + " " + std::to_string(g[3 * cols + 5]));
        img::sobel_s16(dot, rows, cols, 1, 0, 3, 4, border_mode::replicate, g);
        // out(3,3) raw 3 -> 0.75 -> 1; out(4,3) raw 6 -> 1.5 -> 2; dot of 2 would tie at 0.5.
        std::vector<std::uint8_t> dot2(81, 0);
        dot2[4 * cols + 4] = 1;
        std::vector<std::int16_t> g2;
        img::sobel_s16(dot2, rows, cols, 1, 0, 3, 4, border_mode::replicate, g2);
        // out(4,3) raw 2 -> 0.5 -> 0 (even); out(4,5) raw -2 -> -0.5 -> 0.
        check("G9 divisor 4: 0.5 -> 0 and -0.5 -> 0 (ties to even)",
              g2[4 * cols + 3] == 0 && g2[4 * cols + 5] == 0 && g[4 * cols + 3] == 2);
    }
}


// =============================================================================
// H - Bresenham
// =============================================================================
// Verbatim copy of rasterise_line_bresenham_ as it stood in
// lsdtt_raster_skeleton_graph_cpu.cpp before it moved to utils (git HEAD at
// 2026-09-14). The moved routine must reproduce it exactly, tie-breaking
// included, so skeleton bridging output does not change.
std::vector<std::pair<int,int>> pre_move_bresenham(int r0, int c0, int r1, int c1)
{
    std::vector<std::pair<int,int>> pixels;
    int x0 = c0, y0 = r0;
    const int x1 = c1, y1 = r1;
    const int dx = std::abs(x1 - x0);
    const int sx = (x0 < x1) ? 1 : -1;
    const int dy = -std::abs(y1 - y0);
    const int sy = (y0 < y1) ? 1 : -1;
    int err = dx + dy;
    while (true)
    {
        pixels.emplace_back(y0, x0);
        if (x0 == x1 && y0 == y1) break;
        const int e2 = 2 * err;
        if (e2 >= dy) { err += dy; x0 += sx; }
        if (e2 <= dx) { err += dx; y0 += sy; }
    }
    return pixels;
}

void test_bresenham()
{
    std::cout << "\n--- H: bresenham_line ---\n";
    bool endpoints = true, connected = true, length = true, unique = true;
    const int ends[][4] = {
        {0, 0, 0, 10}, {0, 0, 10, 0}, {0, 0, 7, 7}, {5, 5, 0, 0},
        {2, 3, 9, 20}, {9, 20, 2, 3}, {0, 10, 6, 0}, {6, 0, 0, 10},
        {3, 3, 3, 3}, {10, 0, 0, 3}, {-4, -2, 5, 8},
    };
    for (const auto& e : ends)
    {
        const auto p = img::bresenham_line(e[0], e[1], e[2], e[3]);
        if (p.front() != std::make_pair(e[0], e[1]) || p.back() != std::make_pair(e[2], e[3])) endpoints = false;
        for (std::size_t i = 1; i < p.size(); ++i)
            if (std::abs(p[i].first - p[i - 1].first) > 1 || std::abs(p[i].second - p[i - 1].second) > 1)
                connected = false;
        const std::size_t want = static_cast<std::size_t>(std::max(std::abs(e[2] - e[0]), std::abs(e[3] - e[1]))) + 1;
        if (p.size() != want) length = false;
        std::map<std::pair<int, int>, int> seen;
        for (const auto& q : p) if (++seen[q] > 1) unique = false;
    }
    check("H1 both endpoints included, in order", endpoints);
    check("H2 consecutive pixels are 8-connected", connected);
    check("H3 pixel count is max(|dr|,|dc|) + 1 (no redundant steps)", length);
    check("H4 no pixel repeated", unique);
    check("H5 zero-length line is one pixel", img::bresenham_line(4, 4, 4, 4).size() == 1);

    bool same = true;
    for (int r1 = -6; r1 <= 6 && same; ++r1)
        for (int c1 = -6; c1 <= 6 && same; ++c1)
            if (img::bresenham_line(0, 0, r1, c1) != pre_move_bresenham(0, 0, r1, c1)) same = false;
    check("H6 identical to the pre-move skeleton-graph routine for every line from (0,0) to [-6,6]^2", same);
}

} // namespace

int main()
{
    test_random();
    test_border();
    test_gaussian();
    test_box_mean();
    test_separable();
    test_fixed_blur();
    test_otsu();
    test_sobel();
    test_bresenham();

    std::cout << "\n" << (g_n_checks - g_n_failed) << "/" << g_n_checks << " checks passed\n";
    return g_n_failed == 0 ? 0 : 1;
}
