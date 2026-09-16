// image_ops_oracle_inputs.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Deterministic synthetic rasters for the image-ops oracle.
//
// These generated the inputs from which OpenCV 4.13's outputs were captured
// (2026-09-14, before OpenCV was removed; the capture tool was deleted with
// it), and the native image-ops tests regenerate them to compare against those
// fixtures, so both sides use identical inputs. The capture manifest also
// records a hash of every input; a test that finds a different hash knows the
// generators here have drifted and the fixtures are stale.
//
// Every value is computed from closed-form expressions and a fully specified
// splitmix64 hash, never from std::*_distribution, so the inputs are identical
// on every compiler and standard library.
//
// Header-only; no OpenCV dependency.
// -----------------------------------------------------------------------------

#pragma once

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "../../src/cpu/lsdtt_raster_cpu.hpp"

namespace lsdtt_oracle {

struct oracle_input
{
    std::string            name;
    lsdtt::lsdtt_raster_cpu raster;
};

// A planted straight segment in the lines input, in pixel coordinates.
struct planted_segment
{
    int c0, r0, c1, r1;
};

inline std::uint64_t splitmix64(std::uint64_t x)
{
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

// Uniform [0, 1) from (seed, row, col), independent of evaluation order.
inline double hash_uniform(std::uint64_t seed, int r, int c)
{
    std::uint64_t h = splitmix64(seed ^ splitmix64(static_cast<std::uint64_t>(r) * 0x100000001B3ULL
                                                   + static_cast<std::uint64_t>(c)));
    return static_cast<double>(h >> 11) * (1.0 / 9007199254740992.0);
}

inline lsdtt::lsdtt_raster_meta oracle_meta(int n, int n_cols = 0)
{
    lsdtt::lsdtt_raster_meta m;
    m.n_rows        = n;
    m.n_cols        = n_cols > 0 ? n_cols : n;
    m.cell_size     = 1.0;
    m.x_min         = 0.0;
    m.y_min         = 0.0;
    m.no_data_value = -9999.0;
    return m;
}

// Nodata pattern: a 6-pixel border plus an interior disc of radius n/12
// centred at (row n/3, col 5n/8).
inline bool oracle_is_nodata(int n, int r, int c, int n_cols = 0)
{
    const int border = 6;
    const int nc = n_cols > 0 ? n_cols : n;
    if (r < border || c < border || r >= n - border || c >= nc - border) return true;
    const double rr = r - n / 3.0;
    const double cc = c - 5.0 * n / 8.0;
    const double rad = n / 12.0;
    return rr * rr + cc * cc <= rad * rad;
}

inline lsdtt::lsdtt_raster_cpu make_oracle_raster(
    int n, bool with_nodata, double (*z)(int, int, int), int n_cols = 0)
{
    const int nc = n_cols > 0 ? n_cols : n;
    const lsdtt::lsdtt_raster_meta meta = oracle_meta(n, nc);
    std::vector<float>   data(static_cast<std::size_t>(n) * nc);
    std::vector<uint8_t> mask(static_cast<std::size_t>(n) * nc, 1);
    for (int r = 0; r < n; ++r)
    {
        for (int c = 0; c < nc; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(r) * nc + c;
            if (with_nodata && oracle_is_nodata(n, r, c, nc))
            {
                data[i] = static_cast<float>(meta.no_data_value);
                mask[i] = 0;
            }
            else
            {
                data[i] = static_cast<float>(z(n, r, c));
            }
        }
    }
    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

// --- field definitions --------------------------------------------------------
// Values are deliberately positive and far from zero for the DEM-like fields,
// so the legacy normalisation defect (nodata zeros included in min/max) is
// exercised by the with-nodata variants.

inline double z_ramp(int, int r, int c)  { return 1000.0 + 5.0 * c + 2.0 * r; }

inline double z_peaks(int n, int r, int c)
{
    auto g = [&](double r0, double c0, double s, double a)
    {
        const double dr = r - r0 * n, dc = c - c0 * n;
        return a * std::exp(-(dr * dr + dc * dc) / (2.0 * s * s * n * n));
    };
    return 500.0 + g(0.3, 0.3, 0.12, 300.0) + g(0.7, 0.6, 0.08, 220.0)
                 - g(0.5, 0.8, 0.10, 150.0);
}

inline double z_noise(int, int r, int c) { return 100.0 + 50.0 * hash_uniform(7, r, c); }

inline double z_steps(int, int r, int c)
{
    return 200.0 + 10.0 * std::floor((c + r) / 16.0);
}

// Signed curvature-like field spanning zero (VEC-like).
inline double z_signed(int n, int r, int c)
{
    const double k = 2.0 * 3.14159265358979323846 / n;
    return 0.02 * std::sin(3.0 * k * c) * std::cos(2.0 * k * r)
         + 0.005 * (hash_uniform(11, r, c) - 0.5);
}

inline double z_flat(int, int, int) { return 5.0; }

// High-entropy field with fractional values: every 8-bit rounding decision is
// exercised, so float accumulation order, float-vs-double scaling and the
// fixed-point kernel construction all change the output if done differently
// from OpenCV (measured 2026-09-14 against OpenCV 4.13 while it was still
// installed). Non-square (97 x 113) to catch row/column mix-ups.
inline double z_entropy(int, int r, int c)
{
    return 1000.0 + 0.9 * c + 0.4 * r + 60.0 * hash_uniform(101, r, c);
}

inline const std::vector<planted_segment>& lines_planted_segments()
{
    static const std::vector<planted_segment> segs = {
        { 20,  30, 230,  30},   // horizontal
        { 40,  20,  40, 220},   // vertical
        { 60,  60, 220, 220},   // 45 degrees
        { 30, 240, 200, 150},   // shallow, falling
        {120,  40, 150, 200},   // steep
        {180,  20, 240,  45},   // short-ish (length ~65)
    };
    return segs;
}

// Low-amplitude noise with planted bright 1-px segments and one bright arc.
inline double z_lines(int n, int r, int c)
{
    double v = 10.0 * hash_uniform(23, r, c);
    for (const auto& s : lines_planted_segments())
    {
        // Distance from pixel centre to the segment.
        const double vx = s.c1 - s.c0, vy = s.r1 - s.r0;
        const double wx = c - s.c0,   wy = r - s.r0;
        const double len2 = vx * vx + vy * vy;
        double t = len2 > 0.0 ? (wx * vx + wy * vy) / len2 : 0.0;
        t = t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t);
        const double dx = wx - t * vx, dy = wy - t * vy;
        if (dx * dx + dy * dy <= 0.25) v = 100.0;
    }
    const double dr = r - 0.62 * n, dc = c - 0.30 * n;
    const double rho = std::sqrt(dr * dr + dc * dc);
    if (std::fabs(rho - 0.18 * n) <= 0.5 && dr > 0.0) v = 100.0;
    return v;
}

// --- the input set -----------------------------------------------------------

inline std::vector<oracle_input> make_oracle_inputs()
{
    struct field { const char* name; int n; double (*z)(int, int, int); int n_cols; };
    const field fields[] = {
        {"ramp",    128, z_ramp,    0},
        {"peaks",   128, z_peaks,   0},
        {"noise",   128, z_noise,   0},
        {"steps",   128, z_steps,   0},
        {"signed",  128, z_signed,  0},
        {"flat",    128, z_flat,    0},
        {"lines",   256, z_lines,   0},
        {"entropy",  97, z_entropy, 113},
    };
    std::vector<oracle_input> out;
    for (const auto& f : fields)
    {
        out.push_back({std::string(f.name) + "_full",   make_oracle_raster(f.n, false, f.z, f.n_cols)});
        out.push_back({std::string(f.name) + "_nodata", make_oracle_raster(f.n, true,  f.z, f.n_cols)});
    }
    return out;
}

// --- primitive-level oracle input ---------------------------------------------
// A large high-entropy image for checking the filtering primitives' own 8-bit
// output (not a thresholded result) against OpenCV. At this size a single
// different rounding decision anywhere changes the hash. Float variant for the
// normalisation check.
inline constexpr int k_primitive_rows = 512;
inline constexpr int k_primitive_cols = 487;

inline std::vector<std::uint8_t> make_primitive_u8()
{
    std::vector<std::uint8_t> v(static_cast<std::size_t>(k_primitive_rows) * k_primitive_cols);
    for (int r = 0; r < k_primitive_rows; ++r)
        for (int c = 0; c < k_primitive_cols; ++c)
            v[static_cast<std::size_t>(r) * k_primitive_cols + c] =
                static_cast<std::uint8_t>(splitmix64(0xA5A5ULL + static_cast<std::uint64_t>(r) * 100003ULL + c) & 0xFF);
    return v;
}

inline std::vector<float> make_primitive_f32()
{
    std::vector<float> v(static_cast<std::size_t>(k_primitive_rows) * k_primitive_cols);
    for (int r = 0; r < k_primitive_rows; ++r)
        for (int c = 0; c < k_primitive_cols; ++c)
            v[static_cast<std::size_t>(r) * k_primitive_cols + c] =
                static_cast<float>(1000.0 + 0.37 * c - 0.21 * r + 61.3 * hash_uniform(202, r, c));
    return v;
}

// FNV-1a 64-bit over raw bytes: a tripwire hash for inputs and records.
inline std::uint64_t fnv1a64(const unsigned char* p, std::size_t n)
{
    std::uint64_t h = 0xcbf29ce484222325ULL;
    for (std::size_t i = 0; i < n; ++i)
    {
        h ^= p[i];
        h *= 0x100000001b3ULL;
    }
    return h;
}

inline std::uint64_t hash_raster(const lsdtt::lsdtt_raster_cpu& r)
{
    const auto& d = r.get_data();
    const auto& m = r.get_mask();
    std::uint64_t h = fnv1a64(reinterpret_cast<const unsigned char*>(d.data()),
                              d.size() * sizeof(float));
    return h ^ splitmix64(fnv1a64(m.data(), m.size()));
}

} // namespace lsdtt_oracle
