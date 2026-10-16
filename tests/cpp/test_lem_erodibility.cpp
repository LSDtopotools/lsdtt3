// test_lem_erodibility.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Spatially varying erodibility (K) in the LEM.
//
// Modes:
//   lsdtt-test-lem-erodibility [--suite fast|analytic [work_dir]]
//       Run the test suite: all of it by default; "fast" skips the analytic
//       steady states (unit_lem_erodibility), "analytic" runs only them
//       (rigorous_lem_erodibility_analytic). Prints one line per check and a
//       final "RESULT: PASS|FAIL".
//
//   Kernels and units:
//     K1  D8 n = 1, 2: scalar K and a uniform field of the same K, bitwise.
//     K2  D-inf n = 1, 2, 1 and 4 threads: scalar vs uniform field, bitwise.
//     K2b D-inf, non-uniform field: 1, 2, 4 and 0 threads identical.
//     K3  D8 n = 1, random per-pixel K: the closed form with each pixel's own K,
//         recomputed here, agrees to 1 ulp.
//     K4  A field of size N-1 throws, naming both sizes.
//     K5  D8 n = 2, block K: the Newton residual with the pixel's own K is
//         within tolerance; with the receiver's K it is not.
//     K6  D-inf two-receiver pixels (n = 1, 2) and single-receiver Newton
//         pixels (n = 2): as K5. Two-receiver pixels that end below their
//         cardinal receiver are excluded and counted: the kernel's Newton is
//         ill-conditioned there whatever K is.
//     (K5/K6 discrimination uses crossing pixels incising >= 0.01 m, 100x the
//     tolerance; smaller incisions cannot separate the two K values.)
//     U1  describe_erodibility_units strings.
//     U2  erodibility_area_conversion_factor identities.
//     U3  parse_erodibility_area_convention.
//   Calibration:
//     V11 m and n: equal, unstated (warning), different (throws).
//     V12 area convention and calibration cell size (D3b).
//   Steady state (K-weighted chi):
//     S1  Uniform field == opts.K: z bitwise, K_used and scale (D8, D-inf,
//         K-mode and relief mode).
//     S2  200-node chain on the tree core, K2 = K1/4 upstream of node 100:
//         z against a double reference, and U = K_i a^m S^n on every link.
//     S3  z does not depend on the reference K.
//     S4  Relief mode: relief on target, scale = K_used/K, and K-mode with the
//         scaled field reproduces z.
//     S5  D-inf, two-region field: the steady law on single-receiver cells.
//     S5b D-inf two-receiver cells: full-gradient identity, or cardinal fallback.
//     S6  K <= 0 or NaN on a valid pixel throws naming it; on nodata it does not.
//   Values, strict alignment, scaling:
//     V1  Same grid, all valid: values copied, summary exact.
//     V2-V5 Rows, origin, cell size or CRS differ: refused (label, both grids,
//         the placement hint). V6 no CRS: accepted with one warning.
//     V7  V2-V6 via check_lem_field_grid_alignment("lithology code raster").
//     V8  K nodata on a DEM-valid pixel: refused with (r, c) and the count.
//     V9  K nodata on a DEM-nodata pixel: accepted, stored as 0.
//     V10 K <= 0, NaN or inf on a valid pixel: refused.
//     V13 Converted values equal float(double(K) * factor); factor 1 is bitwise.
//     V14 The validator returns warnings and prints nothing.
//     F1  scale_lem_erodibility_field identities on 1e5 random floats.
//     F2  General formula, refusals, overflow and underflow.
//   run_lem_step:
//     L1  5 steps, uniform field == scalar (D8, D-inf, n = 1, 2), bitwise.
//     L2  Rewriting the field between steps equals fresh calls; clearing it
//         gives the scalar step.
//     L3  Wrong size throws naming run_lem_step; K = 0 on a valid pixel throws,
//         on a nodata pixel it does not.
//     L4  K rewritten from z every step for 50 steps (K1 above z_c, K1/4 below):
//         each step equals a fresh kernel call with that step's field, and (D8)
//         the closed form with that step's K.
//     F3  Forcing x field: uniform base == scalar run; two-region field scaled
//         per F2; a K step at t = 500 changes the applied field exactly then.
//   Model grids and input-raster preparation (fixtures written to work_dir):
//     G1  Mode 2 grid from a raster; geographic, US-foot and no-CRS refused.
//     G2  Mode 1 (a): source dims on the synthetic georeference; values bitwise.
//     G3  Centre clip: every window id exact (even and odd sizes; grid lines).
//     G4  Lower-left clip: corner on a grid corner, and inside a pixel.
//     G5  A window off the source on any side, or a point outside, is refused.
//     G6  4 m source on a 10 m model: near, ids of the containing pixels;
//         categorical always near.
//     G7  warp_to_model_grid onto a 0.3-cell shift; no overlap / no CRS refused.
//     G8  as_georeferenced on the shifted grid refused.  G9 categorical bilinear
//         warp refused.  G10 relabel_* on a DEM grid refused.
//     G11 Every placement's output passes the strict check.
//     G12 The in-memory warp leaves no /vsimem/ file.
//     G13 Public clip-window functions (A9), including a 3-band cube grid.
//     G14 A rotated file is refused through the LEM path; empty grids refused.
//     G15 Distortion rule: UTM and LAEA pass; Web Mercator and UTM far
//         outside its zone are refused, in mode 2 and in clips.
//     G16 Mode 2 footprint: source nodata becomes DEM nodata; boundary flats.
//   Analytic steady states (--suite analytic):
//     A1  80x3 strip, m = 0: S_i = (U/K_i)^(1/n) either side of the contact.
//     A1b 60x9, K1 west / K2 east, m = 0.5: K_i a^m S^n = U on trunk pixels.
//     A2  30x30 red noise, two-region field, D8 and D-inf, n = 1 and 2.
//     A3  Snap with a field: the first step after it moves the surface as much as
//         after a scalar snap (p999 within 1%); 200-step values recorded.
//   (A1-A2 run to max |dz| < 1e-3 U dt on every pixel and check the steady law
//   on the incision the step solved, so |q - 1| <= 1e-3/n plus rounding.)
//
//   lsdtt-test-lem-erodibility --write-fixtures <dir>
//       Write the lsdtt-lem CLI test fixtures (DEM, K rasters, placement
//       sources, expected placement results, channel masks, forcing CSVs).
//   lsdtt-test-lem-erodibility --compare-rasters <a> <b>
//       Pass if dims, cell size, origin, CRS, mask and every valid data value
//       (bitwise) agree. --compare-raster-values: dims, mask and data only.
//
//   lsdtt-test-lem-erodibility --print-scalar-hashes
//       Print FNV-1a hashes of the elevation field produced by scalar-K runs:
//       five run_lem_step steps (D8 / D-inf, n = 1 / 2) and the analytical
//       steady state (D8 / D-inf, K-mode / relief mode). These are same-machine
//       before/after references for changes that must leave scalar-K runs
//       bitwise unchanged. They are never compared against stored values.
//   lsdtt-test-lem-erodibility --print-raster-data-hash <raster>...
//       FNV-1a hash of each raster's pixel data and mask only. The georeference
//       is excluded, so a deliberate CRS change does not change the hash.
//   lsdtt-test-lem-erodibility --expect-raster-epsg <raster> <epsg>
//                              [--expect-centre-near <tol_m>]
//       Pass if the raster's CRS identifies as EPSG:<epsg>. With
//       --expect-centre-near, also require the grid centre to be within tol_m
//       of where set_default_synthetic_georeference() would centre a grid of
//       the same dimensions (the Watopia anchor), computed with the same
//       transform rather than hard-coded.
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <random>
#include <stdexcept>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <cpl_string.h>
#include <cpl_vsi.h>
#include <gdal_priv.h>
#include <ogr_spatialref.h>

#include "../../src/core/lsdtt_data_loader.hpp"
#include "../../src/core/lsdtt_raster_meta.hpp"
#include "../../src/core/lsdtt_raster_warp.hpp"
#include "../../src/cpu/lsdtt_flowgraph_cpu.hpp"
#include "../../src/cpu/lsdtt_flowinfo_cpu.hpp"
#include "../../src/cpu/lsdtt_lem_cpu.hpp"
#include "../../src/cpu/lsdtt_lem_input_grid_cpu.hpp"
#include "../../src/cpu/lsdtt_lem_parameter_rasters_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_synthetic_surfaces_cpu.hpp"
#include "../../src/cpu/lsdtt_steady_state_cpu.hpp"
#include "../../src/cpu/lsdtt_stream_power_cpu.hpp"
#include "../../src/utils/lsdtt_grid_neighbours.hpp"

using namespace lsdtt;

namespace
{

std::uint64_t fnv1a(const std::vector<float>& v)
{
    std::uint64_t h = 1469598103934665603ULL;
    const auto* p = reinterpret_cast<const unsigned char*>(v.data());
    for (std::size_t i = 0; i < v.size() * sizeof(float); ++i)
    {
        h ^= p[i];
        h *= 1099511628211ULL;
    }
    return h;
}

std::string hex(std::uint64_t h)
{
    std::ostringstream s;
    s << std::hex << std::setw(16) << std::setfill('0') << h;
    return s.str();
}

// 40 rows x 30 cols of red noise on a northward tilt, pit-filled so that the
// steady-state flow objects (which do not allow pits) can route it.
lsdtt_raster_cpu make_surface(double dx = 30.0)
{
    const int nr = 40, nc = 30;
    lsdtt_raster_cpu z = generate_red_noise_surface(nr, nc, dx, 7, 0.0, 2.0, 6, 1);
    std::vector<float> d = z.get_data();
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
            d[static_cast<std::size_t>(r) * nc + c] += static_cast<float>(0.01 * dx * r);
    std::vector<std::uint8_t> msk = z.get_mask();
    lsdtt_raster_cpu tilted(std::move(d), std::move(msk), z.get_meta());
    return fill_pits(tilted, 1e-4f, 1);
}

lsdtt_lem_options base_options()
{
    lsdtt_lem_options o;
    o.timestep = 200.0f;
    o.erodibility = 2e-5f;
    o.m_exponent = 0.45f;
    o.incision_threshold_area = 0.0f;
    o.enable_hillslope_transport = true;
    o.hillslope_transport_model = "cidre";
    o.hillslope_diffusivity = 5e-3f;
    o.uplift_rate = 1e-4f;
    o.boundary_conditions = "oooo";
    o.n_threads = 1;
    return o;
}

std::vector<float> run_steps(const lsdtt_raster_cpu& dem0, const lsdtt_lem_options& o, int n)
{
    const std::vector<std::int64_t> flats = build_boundary_flats(dem0, o.boundary_conditions);
    lsdtt_raster_cpu dem = dem0;
    for (int i = 0; i < n; ++i) dem = run_lem_step(dem, o, flats).first;
    return dem.get_data();
}

int print_scalar_hashes()
{
    const lsdtt_raster_cpu dem = make_surface();

    struct step_case { const char* name; bool dinf; float n; };
    const step_case cases[] = {
        {"lem_step d8   n=1", false, 1.0f},
        {"lem_step d8   n=2", false, 2.0f},
        {"lem_step dinf n=1", true, 1.0f},
        {"lem_step dinf n=2", true, 2.0f},
    };
    for (const auto& sc : cases)
    {
        lsdtt_lem_options o = base_options();
        o.enable_d8_fluvial_erosion = !sc.dinf;
        o.enable_dinf_fluvial_erosion = sc.dinf;
        o.n_exponent = sc.n;
        std::cout << sc.name << " " << hex(fnv1a(run_steps(dem, o, 5))) << "\n";
    }

    lsdtt_flowinfo_options fo;
    fo.boundary_conditions = "oooo";
    fo.n_threads = 1;
    const lsdtt_flowinfo d8(dem.get_data_array(), dem.get_mask_array(), dem.get_meta(), fo);
    const lsdtt_flowgraph_cpu dinf(dem.get_data_array(), dem.get_mask_array(), dem.get_meta());
    for (int relief = 0; relief < 2; ++relief)
    {
        lsdtt_fluvial_steady_state_options so;
        so.K = 2e-5f;
        so.U = 1e-4f;
        so.relief_mode = relief != 0;
        so.target_relief = 400.0f;
        so.n_threads = 1;
        const auto a = compute_d8_fluvial_steady_state(d8, so);
        const auto b = compute_dinf_fluvial_steady_state(dinf, so);
        const char* mode = relief ? "relief" : "K     ";
        std::cout << "steady d8   " << mode << " " << hex(fnv1a(a.elevation.get_data()))
                  << " K_used " << std::setprecision(9) << a.K_used << "\n";
        std::cout << "steady dinf " << mode << " " << hex(fnv1a(b.elevation.get_data()))
                  << " K_used " << std::setprecision(9) << b.K_used << "\n";
    }
    return 0;
}

// ---------------------------------------------------------------------------
// Test harness
// ---------------------------------------------------------------------------
int g_failures = 0;

void check(const std::string& what, bool ok, const std::string& detail = "")
{
    std::cout << (ok ? "PASS " : "FAIL ") << what;
    if (!detail.empty()) std::cout << " -- " << detail;
    std::cout << "\n";
    if (!ok) ++g_failures;
}

bool same_bits(const std::vector<float>& a, const std::vector<float>& b)
{
    return a.size() == b.size() &&
           std::memcmp(a.data(), b.data(), a.size() * sizeof(float)) == 0;
}

bool throws_containing(const std::function<void()>& f, const std::vector<std::string>& needles,
                       std::string* what_out = nullptr)
{
    try { f(); }
    catch (const std::exception& e)
    {
        if (what_out) *what_out = e.what();
        const std::string w = e.what();
        for (const auto& n : needles)
            if (w.find(n) == std::string::npos) return false;
        return true;
    }
    return false;
}

float ulp(float x)
{
    const float a = std::fabs(x);
    return std::nextafter(a, std::numeric_limits<float>::infinity()) - a;
}

std::string fmt(double v)
{
    std::ostringstream s;
    s << std::setprecision(6) << v;
    return s.str();
}

std::size_t n_cells(const lsdtt_raster_cpu& r)
{
    return static_cast<std::size_t>(r.get_meta().n_rows) * r.get_meta().n_cols;
}

// Checkerboard of 4x4-pixel blocks alternating Ka and Kb, so that many pixels
// drain into a block with a different K.
std::vector<float> block_field(const lsdtt_raster_cpu& dem, float Ka, float Kb)
{
    const int nr = dem.get_meta().n_rows, nc = dem.get_meta().n_cols;
    std::vector<float> k(static_cast<std::size_t>(nr) * nc);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
            k[static_cast<std::size_t>(r) * nc + c] = ((r / 4 + c / 4) % 2) ? Ka : Kb;
    return k;
}

lsdtt_flowinfo d8_flow(const lsdtt_raster_cpu& dem)
{
    lsdtt_flowinfo_options fo;
    fo.boundary_conditions = "oooo";
    fo.n_threads = 1;
    return lsdtt_flowinfo(dem.get_data_array(), dem.get_mask_array(), dem.get_meta(), fo);
}

// ---------------------------------------------------------------------------
// K1-K4
// ---------------------------------------------------------------------------
void test_kernel_bitwise()
{
    const lsdtt_raster_cpu dem = make_surface();
    const float K = 2e-5f;
    const std::vector<float> uniform(n_cells(dem), K);
    const lsdtt_flowinfo d8 = d8_flow(dem);
    const lsdtt_flowgraph_cpu dinf(dem.get_data_array(), dem.get_mask_array(), dem.get_meta());

    for (float n : {1.0f, 2.0f})
    {
        const auto a = compute_d8_stream_power_incision(dem, d8, K, 0.45f, n, 500.0f, 0.0f,
                                                         true, 1e-5f).get_data();
        const auto b = compute_d8_stream_power_incision(dem, d8, K, 0.45f, n, 500.0f, 0.0f,
                                                         true, 1e-5f, 1, "", uniform).get_data();
        check("K1 D8 n=" + fmt(n) + ": uniform field == scalar K, bitwise", same_bits(a, b));
    }
    for (float n : {1.0f, 2.0f})
        for (int t : {1, 4})
        {
            const auto a = compute_dinf_stream_power_incision(dem, dinf, K, 0.45f, n, 500.0f,
                                                               0.0f, true, 1e-5f, t).get_data();
            const auto b = compute_dinf_stream_power_incision(dem, dinf, K, 0.45f, n, 500.0f,
                                                               0.0f, true, 1e-5f, t, 1, "",
                                                               uniform).get_data();
            check("K2 D-inf n=" + fmt(n) + " threads=" + std::to_string(t) +
                      ": uniform field == scalar K, bitwise", same_bits(a, b));
        }

    {
        const std::vector<float> field = block_field(dem, 8e-5f, 2e-5f);
        std::vector<std::vector<float>> res;
        for (int t : {1, 2, 4, 0})
            res.push_back(compute_dinf_stream_power_incision(dem, dinf, K, 0.45f, 2.0f, 500.0f,
                                                              0.0f, true, 1e-5f, t, 1, "",
                                                              field).get_data());
        bool ok = true;
        for (const auto& r : res) ok = ok && same_bits(r, res.front());
        const bool moved = !same_bits(res.front(), dem.get_data());
        check("K2b D-inf non-uniform field: threads 1, 2, 4, 0 identical (and the field erodes)",
              ok && moved);
    }

    // K3: D8 n = 1 closed form with each pixel's own K.
    {
        std::mt19937 rng(11);
        std::uniform_real_distribution<float> U(1e-5f, 1e-3f);
        std::vector<float> field(n_cells(dem));
        for (float& k : field) k = U(rng);
        const float m = 0.45f, dt = 200.0f;
        const auto out = compute_d8_stream_power_incision(dem, d8, K, m, 1.0f, dt, 0.0f, true,
                                                           1e-5f, 1, "", field).get_data();
        const auto& z_old = dem.get_data();
        const std::vector<double> area = d8.accumulate();
        const float dx = static_cast<float>(dem.get_meta().cell_size);
        int n_checked = 0, n_bad = 0;
        for (std::int64_t i = 0; i < static_cast<std::int64_t>(field.size()); ++i)
        {
            if (!d8.is_valid_cell(i) || d8.is_outlet_unchecked(i)) continue;
            const std::int64_t r = d8.receiver_of_unchecked(i);
            const float dx_eff = (d8.receiver_dir_of_unchecked(i) & 1u) ? dx * grid::k_sqrt2_f : dx;
            const float zr = out[static_cast<std::size_t>(r)];
            const float zo = z_old[static_cast<std::size_t>(i)];
            if (zr >= zo) continue;
            const double kam = static_cast<double>(field[static_cast<std::size_t>(i)]) *
                               std::pow(area[static_cast<std::size_t>(i)] / dx,
                                        static_cast<double>(m));
            const double spf = kam * (dt / dx_eff);
            const float expect = static_cast<float>((zo + zr * spf) / (1.0 + spf));
            const float got = out[static_cast<std::size_t>(i)];
            ++n_checked;
            if (std::fabs(got - expect) > ulp(expect)) ++n_bad;
        }
        check("K3 D8 n=1 random K: closed form with each pixel's own K, within 1 ulp",
              n_checked > 500 && n_bad == 0,
              std::to_string(n_checked) + " pixels checked, " + std::to_string(n_bad) + " off");
    }

    {
        const std::vector<float> short_field(n_cells(dem) - 1, K);
        const std::string N = std::to_string(n_cells(dem)), N1 = std::to_string(n_cells(dem) - 1);
        check("K4 D8: field of size N-1 throws naming both sizes",
              throws_containing([&] { (void)compute_d8_stream_power_incision(
                                    dem, d8, K, 0.45f, 1.0f, 500.0f, 0.0f, true, 1e-5f, 1, "",
                                    short_field); }, {N, N1}));
        check("K4 D-inf: field of size N-1 throws naming both sizes",
              throws_containing([&] { (void)compute_dinf_stream_power_incision(
                                    dem, dinf, K, 0.45f, 1.0f, 500.0f, 0.0f, true, 1e-5f, 1, 1,
                                    "", short_field); }, {N, N1}));
    }
}

// ---------------------------------------------------------------------------
// K5: D8 n = 2 Newton residual uses the pixel's own K
// ---------------------------------------------------------------------------
double median_incision(const std::vector<float>& a, const std::vector<float>& b)
{
    std::vector<double> d;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (a[i] - b[i] > 0.0f) d.push_back(a[i] - b[i]);
    if (d.empty()) return 0.0;
    std::nth_element(d.begin(), d.begin() + d.size() / 2, d.end());
    return d[d.size() / 2];
}

void test_k5_d8_newton()
{
    const lsdtt_raster_cpu dem = make_surface(10.0);
    const float m = 0.45f, n = 2.0f, slope_clamp = 1e-5f;
    const float dx = 10.0f;
    const double tol = static_cast<double>(slope_clamp) * dx;   // 1e-4 m
    const float Kb = 2e-4f, Ka = 4.0f * Kb;
    const std::vector<float> field = block_field(dem, Ka, Kb);
    const lsdtt_flowinfo d8 = d8_flow(dem);

    float dt = 100.0f;
    std::vector<float> out;
    for (int i = 0; i < 20; ++i, dt *= 2.0f)
    {
        out = compute_d8_stream_power_incision(dem, d8, Kb, m, n, dt, 0.0f, true, slope_clamp,
                                               1, "", field).get_data();
        if (median_incision(dem.get_data(), out) >= 0.05) break;
    }
    const auto& z_old = dem.get_data();
    const std::vector<double> area = d8.accumulate();

    auto residual = [&](std::int64_t i, float Kx, float z, float zr, float dx_eff) {
        const double kam = static_cast<double>(Kx) *
                           std::pow(area[static_cast<std::size_t>(i)] / dx, static_cast<double>(m));
        const double s = (z - zr) / dx_eff;
        const double f = z - z_old[static_cast<std::size_t>(i)] + kam * dt * std::pow(s, n);
        const double fp = 1.0 + kam * dt * (n / dx_eff) * std::pow(s, n - 1.0);
        return std::fabs(f / fp);
    };

    int n_checked = 0, n_bad = 0, n_cross = 0, n_disc = 0;
    double worst = 0.0;
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(field.size()); ++i)
    {
        if (!d8.is_valid_cell(i) || d8.is_outlet_unchecked(i)) continue;
        const std::int64_t r = d8.receiver_of_unchecked(i);
        const float dx_eff = (d8.receiver_dir_of_unchecked(i) & 1u) ? dx * grid::k_sqrt2_f : dx;
        const float z = out[static_cast<std::size_t>(i)];
        const float zr = out[static_cast<std::size_t>(r)];
        const float zo = z_old[static_cast<std::size_t>(i)];
        if (zr >= zo || z == zo) continue;                                   // untouched / reverted
        if (z == static_cast<float>(zr + slope_clamp * dx_eff) || z <= zr) continue;   // clamped
        const double res = residual(i, field[static_cast<std::size_t>(i)], z, zr, dx_eff);
        const double lim = tol + 2.0 * ulp(z);
        worst = std::max(worst, res);
        ++n_checked;
        if (res > lim) ++n_bad;
        if (field[static_cast<std::size_t>(i)] != field[static_cast<std::size_t>(r)] &&
            zo - z >= 0.01f)
        {
            ++n_cross;
            if (residual(i, field[static_cast<std::size_t>(r)], z, zr, dx_eff) > 10.0 * tol)
                ++n_disc;
        }
    }
    check("K5 D8 n=2: Newton residual with the pixel's own K within tolerance",
          n_checked > 100 && n_bad == 0,
          std::to_string(n_checked) + " pixels, worst |f/f'| " + fmt(worst) + " m, tolerance " +
              fmt(tol) + " m + 2 ulp(z), dt " + fmt(dt) + " yr");
    check("K5 D8 n=2: the receiver's K gives a residual > 10x tolerance on >= 90% of crossings "
          "incising >= 0.01 m",
          n_cross >= 50 && n_disc >= 0.9 * n_cross,
          std::to_string(n_disc) + " of " + std::to_string(n_cross) + " crossing pixels");
}

// ---------------------------------------------------------------------------
// K6: D-inf two-receiver (n = 1, 2) and single-receiver Newton (n = 2)
// ---------------------------------------------------------------------------
void test_k6_dinf_newton()
{
    const lsdtt_raster_cpu dem = make_surface(10.0);
    const float m = 0.45f, slope_clamp = 1e-5f, dx = 10.0f;
    const double tol = static_cast<double>(slope_clamp) * dx;
    const float Kb = 2e-4f, Ka = 4.0f * Kb;
    const std::vector<float> field = block_field(dem, Ka, Kb);
    const lsdtt_flowgraph_cpu fg(dem.get_data_array(), dem.get_mask_array(), dem.get_meta());
    const std::vector<double> area = fg.accumulate();
    const auto& z_old = dem.get_data();

    for (float n : {1.0f, 2.0f})
    {
        float dt = 100.0f;
        std::vector<float> out;
        for (int i = 0; i < 20; ++i, dt *= 2.0f)
        {
            out = compute_dinf_stream_power_incision(dem, fg, Kb, m, n, dt, 0.0f, true,
                                                     slope_clamp, 1, 1, "", field).get_data();
            if (median_incision(z_old, out) >= 0.05) break;
        }
        int two_checked = 0, two_bad = 0, two_clamped = 0, two_below_cardinal = 0;
        int one_checked = 0, one_bad = 0;
        int n_cross = 0, n_disc = 0;
        double worst = 0.0;
        for (std::int64_t i = 0; i < static_cast<std::int64_t>(field.size()); ++i)
        {
            if (fg.is_outlet_flat(i) || !fg.is_valid_flat(i)) continue;
            const std::int64_t nr = fg.n_receivers(i);
            const float w0 = nr > 0 ? fg.receiver_weight_at(i, 0) : 0.0f;
            const float w1 = nr > 1 ? fg.receiver_weight_at(i, 1) : 0.0f;
            const bool only_cardinal = (nr == 1) || (w0 >= 1.0f - 1e-6f);
            const bool only_diagonal = (nr > 1) && (w1 >= 1.0f - 1e-6f);
            const bool two = (nr > 1) && !only_cardinal && !only_diagonal;
            const float z = out[static_cast<std::size_t>(i)];
            const float zo = z_old[static_cast<std::size_t>(i)];
            const double a = area[static_cast<std::size_t>(i)] / dx;
            const float Ki = field[static_cast<std::size_t>(i)];
            if (z == zo) continue;
            if (two)
            {
                const std::int64_t r1 = fg.receiver_index_at(i, 0), r2 = fg.receiver_index_at(i, 1);
                const float ze1 = out[static_cast<std::size_t>(r1)];
                const float ze2 = out[static_cast<std::size_t>(r2)];
                if (z == std::min(ze1, ze2) + slope_clamp * dx || z <= std::min(ze1, ze2))
                {
                    ++two_clamped;
                    continue;
                }
                // Ended below the cardinal receiver: f' changes sign there, so
                // the kernel's Newton is ill-conditioned and stops at its
                // iteration cap with or without a K field. Not a K property.
                if (z < ze1)
                {
                    ++two_below_cardinal;
                    continue;
                }
                auto res2 = [&](float Kx) {
                    const double kam = static_cast<double>(Kx) * std::pow(a, static_cast<double>(m));
                    const double sbt = std::pow((ze1 - z) / dx, 2.0) + std::pow((ze2 - ze1) / dx, 2.0);
                    const double f = z - zo + kam * std::pow(sbt, 0.5 * n) * dt;
                    const double fp = 1.0 - kam * n * (ze1 - z) / (dx * dx) *
                                                std::pow(sbt, 0.5 * n - 1.0) * dt;
                    return std::fabs(f / fp);
                };
                const double res = res2(Ki);
                worst = std::max(worst, res);
                ++two_checked;
                if (res > tol + 2.0 * ulp(z)) ++two_bad;
                const float Kr = field[static_cast<std::size_t>(r1)];
                if (Kr != Ki && zo - z >= 0.01f)
                {
                    ++n_cross;
                    if (res2(Kr) > 10.0 * tol) ++n_disc;
                }
            }
            else if (n != 1.0f)
            {
                const std::int64_t r = only_diagonal ? fg.receiver_index_at(i, 1)
                                                     : fg.receiver_index_at(i, 0);
                const float dx_eff = only_diagonal ? dx * grid::k_sqrt2_f : dx;
                const float zr = out[static_cast<std::size_t>(r)];
                if (zr >= zo || z <= zr || z == static_cast<float>(zr + slope_clamp * dx_eff))
                    continue;
                const double kam = static_cast<double>(Ki) * std::pow(a, static_cast<double>(m));
                const double s = (z - zr) / dx_eff;
                const double f = z - zo + kam * dt * std::pow(s, n);
                const double fp = 1.0 + kam * dt * (n / dx_eff) * std::pow(s, n - 1.0);
                ++one_checked;
                if (std::fabs(f / fp) > tol + 2.0 * ulp(z)) ++one_bad;
            }
        }
        const std::string tag = "K6 D-inf n=" + fmt(n);
        check(tag + ": two-receiver Newton residual with the pixel's own K within tolerance",
              two_checked >= 20 && two_bad == 0,
              std::to_string(two_checked) + " pixels (excluded: " + std::to_string(two_clamped) +
                  " clamped, " + std::to_string(two_below_cardinal) +
                  " below the cardinal receiver), worst |f/f'| " + fmt(worst) + " m, dt " + fmt(dt));
        check(tag + ": the receiver's K gives a residual > 10x tolerance on >= 90% of crossings "
                    "incising >= 0.01 m",
              n_cross >= 10 && n_disc >= 0.9 * n_cross,
              std::to_string(n_disc) + " of " + std::to_string(n_cross) + " crossing pixels");
        if (n != 1.0f)
            check(tag + ": single-receiver Newton residual within tolerance",
                  one_checked >= 50 && one_bad == 0,
                  std::to_string(one_checked) + " pixels, " + std::to_string(one_bad) + " off");
    }
}

// ---------------------------------------------------------------------------
// U1-U3, V11-V12
// ---------------------------------------------------------------------------
void test_units_and_calibration()
{
    check("U1 specific-area units string",
          describe_erodibility_units(0.45f, true) ==
              "m^(0.55) yr^-1 (specific area a = A/dx; K is tied to m and to the cell size)",
          describe_erodibility_units(0.45f, true));
    check("U1 map-area units string",
          describe_erodibility_units(0.45f, false) == "m^(0.1) yr^-1 (map area A; K is tied to m)",
          describe_erodibility_units(0.45f, false));

    using C = erodibility_area_convention;
    const double m = static_cast<double>(0.45f);
    check("U2 map -> map is exactly 1",
          erodibility_area_conversion_factor(0.45f, C::map_area, -1.0, false, 10.0) == 1.0);
    check("U2 specific(10) -> specific(10) is exactly 1",
          erodibility_area_conversion_factor(0.45f, C::specific_area, 10.0, true, 10.0) == 1.0);
    check("U2 map -> specific(10) == pow(10, m)",
          erodibility_area_conversion_factor(0.45f, C::map_area, -1.0, true, 10.0) ==
              std::pow(10.0, m));
    check("U2 specific(5) -> specific(10) == pow(2, m)",
          erodibility_area_conversion_factor(0.45f, C::specific_area, 5.0, true, 10.0) ==
              std::pow(2.0, m));
    {
        const double there = erodibility_area_conversion_factor(0.45f, C::map_area, -1.0, true, 7.3);
        const double back = erodibility_area_conversion_factor(0.45f, C::specific_area, 7.3, false, 7.3);
        const double rt = there * back;
        check("U2 map -> specific -> map round trip within 2 ulp(double)",
              std::fabs(rt - 1.0) <= 2.0 * std::numeric_limits<double>::epsilon(), fmt(rt - 1.0));
    }
    check("U2 a needed cell size <= 0 throws",
          throws_containing([] { (void)erodibility_area_conversion_factor(
                                0.45f, C::map_area, -1.0, true, 0.0); }, {"cell size"}) &&
          throws_containing([] { (void)erodibility_area_conversion_factor(
                                0.45f, C::specific_area, -1.0, true, 10.0); }, {"cell size"}));

    check("U3 parse accepts map_area and specific_area",
          parse_erodibility_area_convention("map_area") == C::map_area &&
              parse_erodibility_area_convention("specific_area") == C::specific_area);
    check("U3 parse of 'map' throws naming both accepted values",
          throws_containing([] { (void)parse_erodibility_area_convention("map"); },
                            {"map_area", "specific_area"}));

    lsdtt_erodibility_run_convention run;
    run.m_exponent = 0.45f;
    run.n_exponent = 1.0f;
    run.divide_area_by_dx = true;
    run.cell_size = 10.0;
    const std::string label = "erodibility raster test.tif";

    auto has = [](const std::vector<std::string>& v, const std::string& needle) {
        for (const auto& s : v)
            if (s.find(needle) != std::string::npos) return true;
        return false;
    };

    {
        lsdtt_erodibility_calibration st;
        st.m_exponent = 0.45f;
        st.n_exponent = 1.0f;
        st.area_convention = C::specific_area;
        st.cell_size = 10.0;
        const auto r = check_lem_erodibility_calibration(st, run, label);
        check("V11 m and n equal, specific(10) on 10 m: factor 1, no notes or warnings",
              r.conversion_factor == 1.0 && r.notes.empty() && r.warnings.empty());

        lsdtt_erodibility_calibration un = st;
        un.m_exponent = -1.0f;
        un.n_exponent = -1.0f;
        const auto r2 = check_lem_erodibility_calibration(un, run, label);
        check("V11 m and n unstated: one warning each, naming the key",
              r2.warnings.size() == 2 && has(r2.warnings, "m_exponent") &&
                  has(r2.warnings, "n_exponent"));

        lsdtt_erodibility_calibration bm = st;
        bm.m_exponent = 0.5f;
        check("V11 different m throws containing the source label",
              throws_containing([&] { (void)check_lem_erodibility_calibration(bm, run, label); },
                                {label, "m_exponent"}));
        lsdtt_erodibility_calibration bn = st;
        bn.n_exponent = 2.0f;
        check("V11 different n throws containing the source label",
              throws_containing([&] { (void)check_lem_erodibility_calibration(bn, run, label); },
                                {label, "n_exponent"}));
    }
    {
        lsdtt_erodibility_calibration st;
        st.m_exponent = 0.45f;
        st.n_exponent = 1.0f;
        st.area_convention = C::map_area;
        const auto r = check_lem_erodibility_calibration(st, run, label);
        check("V12 map stated on a 10 m specific run: factor pow(10, m), note with 2.81838",
              r.conversion_factor == std::pow(10.0, m) && has(r.notes, "2.81838"),
              r.notes.empty() ? "" : r.notes.front());

        lsdtt_erodibility_calibration un = st;
        un.area_convention.reset();
        const auto r2 = check_lem_erodibility_calibration(un, run, label);
        check("V12 convention unstated: warning with dx^m = 2.81838, factor 1",
              r2.conversion_factor == 1.0 && has(r2.warnings, "dx^m = 2.81838"));

        lsdtt_erodibility_calibration s5 = st;
        s5.area_convention = C::specific_area;
        s5.cell_size = 5.0;
        const auto r3 = check_lem_erodibility_calibration(s5, run, label);
        check("V12 specific(5) on 10 m: factor pow(2, m), note 'calibration cell size 5 to 10'",
              r3.conversion_factor == std::pow(2.0, m) &&
                  has(r3.notes, "calibration cell size 5 to 10"),
              r3.notes.empty() ? "" : r3.notes.front());

        lsdtt_erodibility_calibration snc = st;
        snc.area_convention = C::specific_area;
        snc.cell_size = -1.0;
        const auto r4 = check_lem_erodibility_calibration(snc, run, label);
        check("V12 specific with no cell size on a specific run: warning, factor 1",
              r4.conversion_factor == 1.0 && has(r4.warnings, "calibration_cell_size"));

        lsdtt_erodibility_run_convention map_run = run;
        map_run.divide_area_by_dx = false;
        check("V12 specific with no cell size on a map-area run throws",
              throws_containing([&] { (void)check_lem_erodibility_calibration(snc, map_run, label); },
                                {label, "calibration cell size"}));
    }
}

// ---------------------------------------------------------------------------
// S1-S6: steady state with a K field (K-weighted chi)
// ---------------------------------------------------------------------------
float max_valid(const std::vector<float>& z, const std::vector<std::uint8_t>& mask)
{
    float m = -std::numeric_limits<float>::infinity();
    for (std::size_t i = 0; i < z.size(); ++i)
        if (mask[i]) m = std::max(m, z[i]);
    return m;
}

// Two regions: K1 north of the middle row, K2 south of it.
std::vector<float> two_region_field(const lsdtt_raster_cpu& dem, float K1, float K2)
{
    const int nr = dem.get_meta().n_rows, nc = dem.get_meta().n_cols;
    std::vector<float> k(static_cast<std::size_t>(nr) * nc);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
            k[static_cast<std::size_t>(r) * nc + c] = (r < nr / 2) ? K1 : K2;
    return k;
}

void test_steady_state()
{
    const lsdtt_raster_cpu dem = make_surface();
    const lsdtt_flowinfo d8 = d8_flow(dem);
    const lsdtt_flowgraph_cpu dinf(dem.get_data_array(), dem.get_mask_array(), dem.get_meta());
    const float K = 2e-5f, U = 1e-4f;

    // S1: uniform field equal to opts.K.
    for (int relief = 0; relief < 2; ++relief)
    {
        lsdtt_fluvial_steady_state_options o;
        o.K = K;
        o.U = U;
        o.relief_mode = relief != 0;
        o.target_relief = 400.0f;
        o.n_threads = 1;
        lsdtt_fluvial_steady_state_options of = o;
        of.erodibility_field.assign(n_cells(dem), K);
        const std::string mode = relief ? "relief" : "K-mode";
        {
            const auto a = compute_d8_fluvial_steady_state(d8, o);
            const auto b = compute_d8_fluvial_steady_state(d8, of);
            const float want_scale = relief ? b.K_used / o.K : 1.0f;
            check("S1 D8 " + mode + ": uniform field == scalar K (z bitwise, K_used, scale)",
                  same_bits(a.elevation.get_data(), b.elevation.get_data()) &&
                      a.K_used == b.K_used && b.erodibility_field_scale == want_scale &&
                      a.erodibility_field_scale == 1.0f);
        }
        {
            const auto a = compute_dinf_fluvial_steady_state(dinf, o);
            const auto b = compute_dinf_fluvial_steady_state(dinf, of);
            const float want_scale = relief ? b.K_used / o.K : 1.0f;
            check("S1 D-inf " + mode + ": uniform field == scalar K (z bitwise, K_used, scale)",
                  same_bits(a.elevation.get_data(), b.elevation.get_data()) &&
                      a.K_used == b.K_used && b.erodibility_field_scale == want_scale);
        }
    }

    // S2: a 200-node chain on the tree core, K2 = K1/4 upstream of node 100.
    for (float n : {1.0f, 2.0f})
    {
        const int N = 200;
        const double dx = 10.0;
        const float m = 0.5f, K1 = 1e-4f, K2 = K1 / 4.0f;
        std::vector<std::int64_t> receiver(N), topo(N);
        std::vector<double> area(N);
        std::vector<float> dx_eff(N, static_cast<float>(dx)), field(N);
        std::vector<std::uint8_t> valid(N, 1);
        for (int i = 0; i < N; ++i)
        {
            receiver[i] = (i == 0) ? 0 : i - 1;
            area[i] = (N - i) * dx * dx;
            field[i] = (i >= 100) ? K2 : K1;
            topo[i] = N - 1 - i;   // upstream first
        }
        lsdtt_fluvial_steady_state_options o;
        o.U = K1;
        o.K = K1;
        o.m = m;
        o.n = n;
        o.erodibility_field = field;
        const auto r = compute_fluvial_steady_state_on_tree(receiver, area, dx_eff, valid, topo,
                                                            dx, -9999.0, o);
        // Reference in double.
        std::vector<double> zref(N, 0.0);
        const double f = std::pow(static_cast<double>(o.U) / o.K, 1.0 / n);
        double chi = 0.0;
        for (int i = 1; i < N; ++i)
        {
            const double a = area[i] / dx;
            const double w = std::pow(static_cast<double>(o.K) / field[i], 1.0 / n);
            chi += dx * w * std::pow(a, -static_cast<double>(m) / n);
            zref[i] = f * chi;
        }
        const float zmax = *std::max_element(r.elevation.begin(), r.elevation.end());
        double worst = 0.0;
        for (int i = 0; i < N; ++i) worst = std::max(worst, std::fabs(r.elevation[i] - zref[i]));
        const double ztol = 200.0 * ulp(zmax);
        check("S2 chain n=" + fmt(n) + ": z within 200 ulp(z_max) of the double reference",
              worst <= ztol, "worst " + fmt(worst) + " m, tolerance " + fmt(ztol) + " m");

        double dz_min = std::numeric_limits<double>::infinity();
        for (int i = 1; i < N; ++i)
            dz_min = std::min(dz_min, static_cast<double>(r.elevation[i] - r.elevation[i - 1]));
        const double qtol = 4.0 * ulp(zmax) / dz_min + 1e-6;
        double qworst = 0.0;
        for (int i = 1; i < N; ++i)
        {
            const double S = (r.elevation[i] - r.elevation[i - 1]) / dx;
            const double q = std::pow(static_cast<double>(field[i]), 1.0 / n) * S *
                             std::pow(area[i] / dx, static_cast<double>(m) / n) /
                             std::pow(static_cast<double>(o.U), 1.0 / n);
            qworst = std::max(qworst, std::fabs(q - 1.0));
        }
        check("S2 chain n=" + fmt(n) + ": U = K_i a^m S^n on every link (q = 1)",
              qworst <= qtol && qtol < 0.01,
              "worst |q-1| " + fmt(qworst) + ", tolerance " + fmt(qtol));
    }

    const std::vector<float> two = two_region_field(dem, 8e-5f, 2e-5f);

    // S3: the reference K does not change z.
    {
        lsdtt_fluvial_steady_state_options a, b;
        a.U = b.U = U;
        a.K = 2e-5f;
        b.K = 7e-5f;
        a.erodibility_field = b.erodibility_field = two;
        const auto ra = compute_d8_fluvial_steady_state(d8, a);
        const auto rb = compute_d8_fluvial_steady_state(d8, b);
        const auto& za = ra.elevation.get_data();
        const auto& zb = rb.elevation.get_data();
        const float zmax = max_valid(za, ra.elevation.get_mask());
        double worst = 0.0;
        for (std::size_t i = 0; i < za.size(); ++i)
            if (ra.elevation.get_mask()[i]) worst = std::max(worst, double(std::fabs(za[i] - zb[i])));
        check("S3 D8: z independent of the reference K (within 200 ulp(z_max))",
              worst <= 200.0 * ulp(zmax), "worst " + fmt(worst) + " m");
    }

    // S4: relief mode with a field.
    for (int dinf_path = 0; dinf_path < 2; ++dinf_path)
    {
        lsdtt_fluvial_steady_state_options o;
        o.U = U;
        o.K = 2e-5f;
        o.relief_mode = true;
        o.target_relief = 400.0f;
        o.erodibility_field = two;
        o.n_threads = 1;
        const auto r = dinf_path ? compute_dinf_fluvial_steady_state(dinf, o)
                                 : compute_d8_fluvial_steady_state(d8, o);
        const std::string tag = dinf_path ? "S4 D-inf" : "S4 D8";
        check(tag + ": relief within 4 ulp of the target",
              std::fabs(r.computed_relief - o.target_relief) <= 4.0f * ulp(o.target_relief),
              fmt(r.computed_relief));
        check(tag + ": erodibility_field_scale == K_used / opts.K",
              r.erodibility_field_scale == r.K_used / o.K, fmt(r.erodibility_field_scale));
        lsdtt_fluvial_steady_state_options k = o;
        k.relief_mode = false;
        for (float& v : k.erodibility_field) v *= r.erodibility_field_scale;
        const auto rk = dinf_path ? compute_dinf_fluvial_steady_state(dinf, k)
                                  : compute_d8_fluvial_steady_state(d8, k);
        const auto& z1 = r.elevation.get_data();
        const auto& z2 = rk.elevation.get_data();
        const float zmax = max_valid(z1, r.elevation.get_mask());
        double worst = 0.0;
        for (std::size_t i = 0; i < z1.size(); ++i)
            if (r.elevation.get_mask()[i]) worst = std::max(worst, double(std::fabs(z1[i] - z2[i])));
        check(tag + ": K-mode with scale * field reproduces the relief-mode z",
              worst <= 200.0 * ulp(zmax), "worst " + fmt(worst) + " m");
    }

    // S5 / S5b: D-inf, two-region field: the steady law holds cell by cell.
    {
        lsdtt_fluvial_steady_state_options o;
        o.U = U;
        o.K = 2e-5f;
        o.erodibility_field = two;
        o.n_threads = 1;
        const auto r = compute_dinf_fluvial_steady_state(dinf, o);
        const auto& z = r.elevation.get_data();
        const std::vector<double> area = dinf.accumulate();
        const float dx = static_cast<float>(dem.get_meta().cell_size);
        const float zmax = max_valid(z, r.elevation.get_mask());
        const double rel_m = static_cast<double>(o.m) / o.n;
        int one = 0, one_bad = 0, two_full = 0, two_fallback = 0, two_bad = 0;
        double worst1 = 0.0, worst2 = 0.0;
        for (std::int64_t i = 0; i < static_cast<std::int64_t>(z.size()); ++i)
        {
            if (dinf.is_outlet_flat(i) || !dinf.is_valid_flat(i)) continue;
            const std::int64_t nr = dinf.n_receivers(i);
            const float w0 = nr > 0 ? dinf.receiver_weight_at(i, 0) : 0.0f;
            const float w1 = nr > 1 ? dinf.receiver_weight_at(i, 1) : 0.0f;
            const bool only_diag = nr > 1 && w1 >= 1.0f - 1e-6f;
            const bool two_r = nr > 1 && !(w0 >= 1.0f - 1e-6f) && !only_diag;
            const double a = area[static_cast<std::size_t>(i)] / dx;
            const double kfac = std::pow(two[static_cast<std::size_t>(i)] / o.U, 1.0 / o.n);
            if (two_r)
            {
                const float ze1 = z[static_cast<std::size_t>(dinf.receiver_index_at(i, 0))];
                const float ze2 = z[static_cast<std::size_t>(dinf.receiver_index_at(i, 1))];
                const double zi = z[static_cast<std::size_t>(i)];
                const double dz = std::min(zi - ze1, std::fabs(double(ze1) - ze2) + (zi - ze1));
                const double tol = 4.0 * ulp(zmax) / std::max(dz, 1e-12) + 1e-6;
                const double g = std::hypot((zi - ze1) / dx, (double(ze1) - ze2) / dx);
                const double q_full = g * std::pow(a, rel_m) * kfac;
                const double q_card = (zi - ze1) / dx * std::pow(a, rel_m) * kfac;
                if (std::fabs(q_full - 1.0) <= tol) { ++two_full; worst2 = std::max(worst2, std::fabs(q_full - 1.0)); }
                else if (std::fabs(q_card - 1.0) <= tol) ++two_fallback;
                else ++two_bad;
            }
            else
            {
                const int k0 = only_diag ? 1 : 0;
                const std::int64_t rr = dinf.receiver_index_at(i, k0);
                const int dr = std::abs(dinf.row_from_flat(rr) - dinf.row_from_flat(i));
                const int dc = std::abs(dinf.col_from_flat(rr) - dinf.col_from_flat(i));
                const double dxe = (dr == 1 && dc == 1) ? dx * grid::k_sqrt2_f : dx;
                const double dz = double(z[static_cast<std::size_t>(i)]) - z[static_cast<std::size_t>(rr)];
                const double tol = 4.0 * ulp(zmax) / dz + 1e-6;
                const double q = dz / dxe * std::pow(a, rel_m) * kfac;
                worst1 = std::max(worst1, std::fabs(q - 1.0));
                ++one;
                if (std::fabs(q - 1.0) > tol) ++one_bad;
            }
        }
        check("S5 D-inf two-region field: single-receiver cells satisfy U = K_i a^m S^n",
              one >= 50 && one_bad == 0,
              std::to_string(one) + " cells, worst |q-1| " + fmt(worst1));
        check("S5b D-inf two-receiver cells: |grad z| identity (>= 100 cells) or cardinal fallback",
              two_full >= 100 && two_bad == 0,
              std::to_string(two_full) + " full-gradient, " + std::to_string(two_fallback) +
                  " fallback, " + std::to_string(two_bad) + " neither; worst |q-1| " + fmt(worst2));
    }

    // S6: bad K on a valid pixel throws, on a nodata pixel does not.
    {
        lsdtt_raster_cpu holed = dem;
        const int nc = dem.get_meta().n_cols;
        const std::size_t hole = static_cast<std::size_t>(20) * nc + 15;
        holed.get_mask_array().data()[hole] = 0;
        holed.get_data_array().data()[hole] = static_cast<float>(dem.get_meta().no_data_value);
        lsdtt_flowinfo_options fo;
        fo.boundary_conditions = "oooo";
        const lsdtt_flowinfo hd8(holed.get_data_array(), holed.get_mask_array(), holed.get_meta(), fo);
        const lsdtt_flowgraph_cpu hdinf(holed.get_data_array(), holed.get_mask_array(), holed.get_meta());
        const std::size_t good = static_cast<std::size_t>(10) * nc + 7;
        for (float bad : {0.0f, -1e-5f, std::numeric_limits<float>::quiet_NaN()})
        {
            lsdtt_fluvial_steady_state_options o;
            o.U = U;
            o.K = K;
            o.erodibility_field.assign(n_cells(dem), K);
            o.erodibility_field[good] = bad;
            const std::string v = fmt(bad);
            check("S6 K = " + v + " on a valid pixel: D8 and D-inf throw 'K must be > 0' at (10, 7)",
                  throws_containing([&] { (void)compute_d8_fluvial_steady_state(hd8, o); },
                                    {"K must be > 0", "(10, 7)"}) &&
                      throws_containing([&] { (void)compute_dinf_fluvial_steady_state(hdinf, o); },
                                        {"K must be > 0", "(10, 7)"}));
            o.erodibility_field[good] = K;
            o.erodibility_field[hole] = bad;
            bool ok = true;
            try
            {
                (void)compute_d8_fluvial_steady_state(hd8, o);
                (void)compute_dinf_fluvial_steady_state(hdinf, o);
            }
            catch (const std::exception&) { ok = false; }
            check("S6 K = " + v + " on a nodata pixel: accepted", ok);
        }
    }
}

// ---------------------------------------------------------------------------
// V1-V10, V13, V14, F1, F2: values, strict alignment, scaling
// ---------------------------------------------------------------------------
lsdtt_raster_meta utm_meta(int nr, int nc, double dx, int epsg)
{
    lsdtt_raster_meta m;
    m.n_rows = nr;
    m.n_cols = nc;
    m.cell_size = dx;
    m.x_min = 500000.0;
    m.y_min = 8700000.0;
    m.no_data_value = -9999.0;
    if (epsg > 0) m.set_crs_from_epsg(epsg);
    return m;
}

lsdtt_raster_cpu raster_on(const lsdtt_raster_meta& m, float value)
{
    const std::size_t n = static_cast<std::size_t>(m.n_rows) * m.n_cols;
    return lsdtt_raster_cpu(std::vector<float>(n, value), std::vector<std::uint8_t>(n, 1), m);
}

void test_values_and_scaling()
{
    const lsdtt_raster_meta gm = utm_meta(6, 5, 10.0, 32758);
    const lsdtt_raster_cpu dem = raster_on(gm, 100.0f);
    lsdtt_lem_erodibility_field_options opts;
    opts.reference_erodibility = 1e-4f;
    opts.stated.m_exponent = 0.45f;
    opts.stated.n_exponent = 1.0f;
    opts.stated.area_convention = erodibility_area_convention::specific_area;
    opts.stated.cell_size = 10.0;
    opts.source_label = "erodibility";

    // V1
    {
        lsdtt_raster_cpu k = raster_on(gm, 1e-4f);
        k.get_data_array().data()[7] = 3e-4f;
        const auto f = validate_lem_erodibility_field(k, dem, opts);
        check("V1 same grid, all valid: values copied bitwise, summary exact",
              same_bits(f.values, k.get_data()) && f.n_valid == 30 && f.min_value == 1e-4f &&
                  f.max_value == 3e-4f && f.conversion_factor == 1.0 && f.warnings.empty());
    }

    // V2-V7: each mismatch throws through the strict check.
    struct mismatch { std::string what; lsdtt_raster_meta meta; std::vector<std::string> needles; };
    lsdtt_raster_meta rows = gm; rows.n_rows = 7;
    lsdtt_raster_meta origin = gm; origin.x_min += 5.0;
    lsdtt_raster_meta cell = gm; cell.cell_size = 12.0;
    const lsdtt_raster_meta crs = utm_meta(6, 5, 10.0, 32601);
    const std::vector<mismatch> cases = {
        {"V2 rows differ", rows, {"raster grid", "7x5", "6x5", "input_raster_placement"}},
        {"V3 origin differs by half a cell", origin, {"raster grid", "6x5", "input_raster_placement"}},
        {"V4 cell size differs", cell, {"raster grid", "cell 12", "input_raster_placement"}},
        {"V5 CRS differs (UTM 58S vs 1N)", crs, {"raster grid", "CRS", "input_raster_placement"}},
    };
    for (const auto& c : cases)
    {
        const lsdtt_raster_cpu k = raster_on(c.meta, 1e-4f);
        std::string what;
        const bool ok = throws_containing([&] { (void)validate_lem_erodibility_field(k, dem, opts); },
                                          c.needles, &what);
        check(c.what + ": refused with label, both grids and the placement hint", ok, what);
        std::vector<std::string> lith = c.needles;
        lith.push_back("lithology code raster");
        check("V7 " + c.what + " via check_lem_field_grid_alignment('lithology code raster')",
              throws_containing([&] { (void)check_lem_field_grid_alignment(
                                    c.meta, gm, "lithology code raster"); }, lith));
    }
    {
        lsdtt_raster_meta nocrs = gm;
        nocrs.crs_wkt.clear();
        const lsdtt_raster_cpu k = raster_on(nocrs, 1e-4f);
        const auto f = validate_lem_erodibility_field(k, dem, opts);
        int n_nocrs = 0;
        for (const auto& w : f.warnings)
            if (w.find("no CRS") != std::string::npos) ++n_nocrs;
        check("V6 no CRS, same geometry: accepted with one 'no CRS' warning", n_nocrs == 1);
        const auto w = check_lem_field_grid_alignment(nocrs, gm, "lithology code raster");
        check("V7 no CRS via check_lem_field_grid_alignment: one warning", w.size() == 1);
    }

    // V8 / V9: K nodata.
    {
        lsdtt_raster_cpu k = raster_on(gm, 1e-4f);
        k.get_mask_array().data()[13] = 0;   // (2, 3)
        k.get_mask_array().data()[21] = 0;
        std::string what;
        check("V8 K nodata on DEM-valid pixels: throws with (r, c) and the count",
              throws_containing([&] { (void)validate_lem_erodibility_field(k, dem, opts); },
                                {"nodata on 2 DEM-valid pixels", "(2, 3)"}, &what), what);
        lsdtt_raster_cpu holed = dem;
        holed.get_mask_array().data()[13] = 0;
        holed.get_mask_array().data()[21] = 0;
        const auto f = validate_lem_erodibility_field(k, holed, opts);
        check("V9 K nodata on DEM-nodata pixels: accepted, stored as 0",
              f.values[13] == 0.0f && f.values[21] == 0.0f && f.n_valid == 28);
    }

    // V10: K <= 0, NaN, inf.
    for (float bad : {-1e-4f, 0.0f, std::numeric_limits<float>::quiet_NaN(),
                      std::numeric_limits<float>::infinity()})
    {
        lsdtt_raster_cpu k = raster_on(gm, 1e-4f);
        k.get_data_array().data()[13] = bad;
        check("V10 K = " + fmt(bad) + " on a valid pixel: throws 'K must be > 0' at (2, 3)",
              throws_containing([&] { (void)validate_lem_erodibility_field(k, dem, opts); },
                                {"K must be > 0", "(2, 3)"}));
    }

    // V13: conversion of the values.
    {
        std::mt19937 rng(5);
        std::uniform_real_distribution<float> U(1e-6f, 1e-3f);
        std::vector<float> v(30);
        for (float& x : v) x = U(rng);
        lsdtt_lem_erodibility_field_options mo = opts;
        mo.stated.area_convention = erodibility_area_convention::map_area;
        const auto f = validate_lem_erodibility_values(v, {}, dem, mo);
        bool ok = f.conversion_factor == std::pow(10.0, static_cast<double>(0.45f));
        for (std::size_t i = 0; i < v.size(); ++i)
            ok = ok && f.values[i] == static_cast<float>(static_cast<double>(v[i]) * f.conversion_factor);
        check("V13 map-area values: each equals float(double(K) * factor)", ok);
        const auto g = validate_lem_erodibility_values(v, {}, dem, opts);
        check("V13 factor 1: values bitwise equal to the input", same_bits(g.values, v));
    }

    // V14: returns its warnings, prints nothing.
    {
        std::ostringstream cap_out, cap_err;
        std::streambuf* o = std::cout.rdbuf(cap_out.rdbuf());
        std::streambuf* e = std::cerr.rdbuf(cap_err.rdbuf());
        lsdtt_lem_erodibility_field_options uo = opts;
        uo.stated = lsdtt_erodibility_calibration{};   // everything unstated -> warnings
        const auto f = validate_lem_erodibility_values(std::vector<float>(30, 1.0f), {}, dem, uo);
        std::cout.rdbuf(o);
        std::cerr.rdbuf(e);
        check("V14 validate_lem_erodibility_values returns warnings and prints nothing",
              f.warnings.size() >= 3 && cap_out.str().empty() && cap_err.str().empty(),
              std::to_string(f.warnings.size()) + " warnings");
    }

    // F1: identities over random floats.
    {
        std::mt19937 rng(9);
        std::uniform_real_distribution<float> U(-30.0f, 0.0f);
        std::vector<float> base(100000);
        for (float& b : base) b = std::pow(10.0f, U(rng) / 3.0f);
        const float ref = 3.7e-5f;
        std::vector<float> out;
        scale_lem_erodibility_field(base, ref, ref, out);
        const bool id1 = same_bits(out, base);
        std::vector<float> refs(base.size(), ref);
        scale_lem_erodibility_field(refs, ref, 9.1e-4f, out);
        bool id2 = true;
        for (float x : out) id2 = id2 && x == 9.1e-4f;
        check("F1 applied == reference: out == base bitwise (1e5 random floats)", id1);
        check("F1 base == reference: out == applied bitwise", id2);
    }

    // F2: the general formula, and the refusals.
    {
        const std::vector<float> base = {1e-4f, 0.0f, 2.5e-5f};
        std::vector<float> out;
        scale_lem_erodibility_field(base, 2e-5f, 7e-5f, out);
        bool ok = out[1] == 0.0f;
        for (std::size_t i : {std::size_t(0), std::size_t(2)})
            ok = ok && out[i] == static_cast<float>((static_cast<double>(base[i]) * 7e-5f) / 2e-5f);
        check("F2 general formula exact; nodata (0) stays 0", ok);
        check("F2 reference <= 0 throws",
              throws_containing([&] { scale_lem_erodibility_field(base, 0.0f, 1e-5f, out); },
                                {"reference K must be > 0"}));
        check("F2 applied <= 0 throws",
              throws_containing([&] { scale_lem_erodibility_field(base, 1e-5f, -1e-5f, out); },
                                {"applied K must be > 0"}));
        bool overflow = false, underflow = false;
        try { scale_lem_erodibility_field({1e-4f, 3e38f}, 1.0f, 10.0f, out); }
        catch (const std::range_error& e) { overflow = std::string(e.what()).find("pixel 1") != std::string::npos; }
        try { scale_lem_erodibility_field({1e-38f}, 1.0f, 1e-10f, out); }
        catch (const std::range_error& e) { underflow = std::string(e.what()).find("pixel 0") != std::string::npos; }
        check("F2 3e38 x 10 throws range_error naming the pixel", overflow);
        check("F2 1e-38 x 1e-10 (underflow to 0) throws range_error naming the pixel", underflow);
    }

    // The per-step check used by run_lem_step.
    {
        std::vector<float> f(30, 1e-4f);
        f[13] = 0.0f;
        check("check_lem_erodibility_field: K = 0 on a valid pixel throws at (2, 3)",
              throws_containing([&] { check_lem_erodibility_field(f, dem); },
                                {"K must be > 0", "(2, 3)"}));
        check("check_lem_erodibility_field: wrong size throws naming both sizes",
              throws_containing([&] { check_lem_erodibility_field(std::vector<float>(29, 1e-4f), dem); },
                                {"29", "30"}));
    }
}

// ---------------------------------------------------------------------------
// P4: run_lem_step with a field (L1-L4), analytic tests (A1, A1b, A2, A3), F3
// ---------------------------------------------------------------------------
lsdtt_lem_options fluvial_only(bool dinf, float n, float K, float U, float dt, float m = 0.45f)
{
    lsdtt_lem_options o;
    o.timestep = dt;
    o.erodibility = K;
    o.m_exponent = m;
    o.n_exponent = n;
    o.incision_threshold_area = 0.0f;
    o.enable_hillslope_transport = false;
    o.enable_d8_fluvial_erosion = !dinf;
    o.enable_dinf_fluvial_erosion = dinf;
    o.uplift_rate = U;
    o.enable_uplift = U > 0.0f;
    o.boundary_conditions = "oooo";
    o.n_threads = 1;
    return o;
}

// The D8 flow object run_lem_step builds for its fluvial step.
lsdtt_flowinfo lem_d8_flow(const lsdtt_raster_cpu& dem, const std::string& bc)
{
    lsdtt_flowinfo_options fo;
    fo.boundary_conditions = bc;
    fo.interior_nodata_is_outlet = true;
    fo.allow_pits_as_outlets = true;
    return lsdtt_flowinfo(dem.get_data_array(), dem.get_mask_array(), dem.get_meta(), fo);
}

lsdtt_flowgraph_cpu lem_dinf_flow(const lsdtt_raster_cpu& dem, const std::string& bc)
{
    lsdtt_flowgraph_options fg;
    fg.boundary_conditions = bc;
    fg.interior_nodata_is_outlet = true;
    fg.method = lsdtt_flowgraph_method::d_infinity;
    fg.n_threads = 1;
    return lsdtt_flowgraph_cpu(dem.get_data_array(), dem.get_mask_array(), dem.get_meta(), fg);
}

void test_run_lem_step()
{
    const lsdtt_raster_cpu dem = make_surface();
    const float K = 2e-5f;

    // L1
    for (bool dinf : {false, true})
        for (float n : {1.0f, 2.0f})
        {
            lsdtt_lem_options a = base_options();
            a.enable_d8_fluvial_erosion = !dinf;
            a.enable_dinf_fluvial_erosion = dinf;
            a.n_exponent = n;
            a.erodibility = K;
            lsdtt_lem_options b = a;
            b.erodibility_field.assign(n_cells(dem), K);
            check(std::string("L1 ") + (dinf ? "D-inf" : "D8") + " n=" + fmt(n) +
                      ": 5 steps, uniform field == scalar, bitwise",
                  same_bits(run_steps(dem, a, 5), run_steps(dem, b, 5)));
        }

    // L2
    {
        const std::vector<std::int64_t> flats = build_boundary_flats(dem, "oooo");
        lsdtt_lem_options o = base_options();
        lsdtt_raster_cpu z1 = dem, z2 = dem;
        bool same = true;
        for (int step = 0; step < 4; ++step)
        {
            const std::vector<float> f = block_field(dem, K * (1.0f + step), K / (1.0f + step));
            o.erodibility_field = f;                       // rewritten in place
            lsdtt_lem_options fresh = base_options();      // a new options object
            fresh.erodibility_field = f;
            z1 = run_lem_step(z1, o, flats).first;
            z2 = run_lem_step(z2, fresh, flats).first;
            same = same && same_bits(z1.get_data(), z2.get_data());
        }
        check("L2 per-step rewrite of the field equals fresh calls, bitwise", same);
        o.erodibility_field.clear();
        o.erodibility = K;
        lsdtt_lem_options scalar = base_options();
        scalar.erodibility = K;
        check("L2 clearing the field gives the scalar step",
              same_bits(run_lem_step(z1, o, flats).first.get_data(),
                        run_lem_step(z1, scalar, flats).first.get_data()));
    }

    // L3
    {
        const std::vector<std::int64_t> flats = build_boundary_flats(dem, "oooo");
        lsdtt_lem_options o = base_options();
        o.erodibility_field.assign(n_cells(dem) - 1, K);
        check("L3 wrong size throws naming run_lem_step and both sizes",
              throws_containing([&] { (void)run_lem_step(dem, o, flats); },
                                {"run_lem_step", std::to_string(n_cells(dem)),
                                 std::to_string(n_cells(dem) - 1)}));
        o.erodibility_field.assign(n_cells(dem), K);
        o.erodibility_field[35] = 0.0f;
        check("L3 K = 0 on a valid pixel throws 'K must be > 0'",
              throws_containing([&] { (void)run_lem_step(dem, o, flats); },
                                {"run_lem_step", "K must be > 0", "(1, 5)"}));
        lsdtt_raster_cpu holed = dem;
        holed.get_mask_array().data()[35] = 0;
        bool ok = true;
        try { (void)run_lem_step(holed, o, build_boundary_flats(holed, "oooo")); }
        catch (const std::exception&) { ok = false; }
        check("L3 K = 0 on a nodata pixel does not throw", ok);
    }

    // L4: elevation-dependent rewrite every step (the lithology mechanism).
    for (bool dinf : {false, true})
    {
        lsdtt_lem_options o = fluvial_only(dinf, 1.0f, 1e-4f, 0.0f, 200.0f);
        const std::vector<std::int64_t> flats = build_boundary_flats(dem, "oooo");
        const float K1 = 1e-4f, K2 = K1 / 4.0f;
        lsdtt_raster_cpu z = dem;
        const float zc = 0.5f * (*std::max_element(dem.get_data().begin(), dem.get_data().end()));
        bool bitwise_ok = true, closed_ok = true;
        int n_checked = 0, n_k2 = 0, n_k1 = 0;
        for (int step = 0; step < 50; ++step)
        {
            std::vector<float> f(n_cells(dem));
            for (std::size_t i = 0; i < f.size(); ++i) f[i] = z.get_data()[i] > zc ? K1 : K2;
            o.erodibility_field = f;
            const lsdtt_raster_cpu next = run_lem_step(z, o, flats).first;
            // Fresh recomputation with the K written for this step.
            lsdtt_raster_cpu fresh = dinf
                ? compute_dinf_stream_power_incision(z, lem_dinf_flow(z, "oooo"), K1, o.m_exponent,
                                                     1.0f, o.timestep, 0.0f, true, o.slope_clamp,
                                                     1, 1, "oooo", f)
                : compute_d8_stream_power_incision(z, lem_d8_flow(z, "oooo"), K1, o.m_exponent,
                                                   1.0f, o.timestep, 0.0f, true, o.slope_clamp,
                                                   1, "oooo", f);
            std::vector<float> fz = fresh.get_data();
            for (std::size_t k = 0; k < flats.size(); ++k)   // boundary restore
                fz[static_cast<std::size_t>(flats[k])] = z.get_data()[static_cast<std::size_t>(flats[k])];
            bitwise_ok = bitwise_ok && same_bits(next.get_data(), fz);
            if (!dinf)
            {
                const lsdtt_flowinfo fl = lem_d8_flow(z, "oooo");
                const std::vector<double> area = fl.accumulate();
                const float dx = static_cast<float>(dem.get_meta().cell_size);
                std::vector<std::uint8_t> fixed(f.size(), 0);
                for (std::int64_t b : flats) fixed[static_cast<std::size_t>(b)] = 1;
                for (std::int64_t i = 0; i < static_cast<std::int64_t>(f.size()); ++i)
                {
                    if (!fl.is_valid_cell(i) || fl.is_outlet_unchecked(i)) continue;
                    if (fixed[static_cast<std::size_t>(i)]) continue;   // restored after the step
                    const std::int64_t r = fl.receiver_of_unchecked(i);
                    const float zo = z.get_data()[static_cast<std::size_t>(i)];
                    // The receiver as the kernel saw it: boundary pixels can erode
                    // in the sweep and are only restored after it.
                    const float zr = fresh.get_data()[static_cast<std::size_t>(r)];
                    if (zr >= zo) continue;
                    const float dxe = (fl.receiver_dir_of_unchecked(i) & 1u) ? dx * grid::k_sqrt2_f : dx;
                    const double kam = static_cast<double>(f[static_cast<std::size_t>(i)]) *
                                       std::pow(area[static_cast<std::size_t>(i)] / dx,
                                                static_cast<double>(o.m_exponent));
                    const double spf = kam * (o.timestep / dxe);
                    const float expect = static_cast<float>((zo + zr * spf) / (1.0 + spf));
                    if (std::fabs(next.get_data()[static_cast<std::size_t>(i)] - expect) > ulp(expect))
                        closed_ok = false;
                    ++n_checked;
                    (f[static_cast<std::size_t>(i)] == K2 ? n_k2 : n_k1)++;
                }
            }
            z = next;
        }
        const std::string tag = dinf ? "L4 D-inf" : "L4 D8";
        check(tag + ": 50 steps with K rewritten from z each step; each step equals a fresh "
                    "kernel call with that step's field", bitwise_ok);
        if (!dinf)
            check(tag + ": per-pixel lowering matches the closed form with that step's K (1 ulp)",
                  closed_ok && n_k1 > 1000 && n_k2 > 1000,
                  std::to_string(n_checked) + " pixel-steps (" + std::to_string(n_k1) + " K1, " +
                      std::to_string(n_k2) + " K2)");
    }
}

// Run fluvial-only steps until max |dz| < 1e-3 U dt on EVERY pixel. Returns
// the step count, or -1 if the cap is hit (a failure). The maximum, not p999,
// because the identities below are checked on every pixel: at the stop each
// pixel's erosion rate is within 1e-3 of U, so the identity
// q = (K_i a^m S^n / U)^(1/n) holds to 1e-3/n plus float rounding.
int run_to_steady(lsdtt_raster_cpu& z, const lsdtt_lem_options& o, int cap)
{
    const std::vector<std::int64_t> flats = build_boundary_flats(z, o.boundary_conditions);
    const float stop = 1e-3f * o.uplift_rate * o.timestep;
    for (int i = 0; i < cap; ++i)
    {
        auto r = run_lem_step(z, o, flats);
        z = std::move(r.first);
        if (r.second.max_local_elevation_change < stop)
        {
            // Precondition: the stop is meaningful only above float noise.
            const float zmax = *std::max_element(z.get_data().begin(), z.get_data().end());
            if (!(stop >= 16.0f * ulp(zmax)))
                throw std::logic_error("run_to_steady: 1e-3 U dt is below 16 ulp(z_max)");
            return i + 1;
        }
    }
    return -1;
}

lsdtt_raster_cpu strip(int nr, int nc, double dx, float slope)
{
    lsdtt_raster_meta m;
    m.n_rows = nr;
    m.n_cols = nc;
    m.cell_size = dx;
    m.no_data_value = -9999.0;
    std::vector<float> z(static_cast<std::size_t>(nr) * nc);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
            z[static_cast<std::size_t>(r) * nc + c] =
                static_cast<float>(slope * dx * (nr - 1 - r) + 0.01 * std::abs(c - nc / 2));
    return lsdtt_raster_cpu(std::move(z), std::vector<std::uint8_t>(z.size(), 1), m);
}

// The incision the implicit solver produces from surface z in one step, before
// uplift and before the boundary restore (the kernel's own output, with the
// flow object run_lem_step builds). At steady state z - zp = U dt on every
// pixel that is not a fixed boundary pixel, and the solver's equation then
// gives K_i a^m S'^n = U with S' measured on zp. Measuring S on zp rather than
// on z matters next to the boundary: boundary pixels that are not flow outlets
// erode within the sweep and are restored only afterwards.
std::vector<double> post_incision_surface(const lsdtt_raster_cpu& z, const lsdtt_lem_options& o,
                                          const std::vector<float>& field)
{
    const lsdtt_raster_cpu zp = o.enable_dinf_fluvial_erosion
        ? compute_dinf_stream_power_incision(z, lem_dinf_flow(z, o.boundary_conditions),
                                             o.erodibility, o.m_exponent, o.n_exponent,
                                             o.timestep, 0.0f, true, o.slope_clamp, 1, 1,
                                             o.boundary_conditions, field)
        : compute_d8_stream_power_incision(z, lem_d8_flow(z, o.boundary_conditions),
                                           o.erodibility, o.m_exponent, o.n_exponent,
                                           o.timestep, 0.0f, true, o.slope_clamp, 1,
                                           o.boundary_conditions, field);
    return std::vector<double>(zp.get_data().begin(), zp.get_data().end());
}

void test_analytic()
{
    // A1: 80x3 strip, m = 0, contact half way: S_i = (U/K_i)^(1/n).
    for (float n : {1.0f, 2.0f})
    {
        const int nr = 80, nc = 3;
        const double dx = 10.0;
        const float U = 1e-3f, K1 = 1e-2f, K2 = K1 / 4.0f, dt = 1000.0f;
        lsdtt_raster_cpu z = strip(nr, nc, dx, 0.02f);
        lsdtt_lem_options o = fluvial_only(false, n, K1, U, dt, 0.0f);
        o.boundary_conditions = "nnon";
        std::vector<float> f(n_cells(z));
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
                f[static_cast<std::size_t>(r) * nc + c] = (r < nr / 2) ? K1 : K2;   // K1 upstream
        o.erodibility_field = f;
        const int steps = run_to_steady(z, o, 20000);
        const auto& d = z.get_data();
        auto slope_at = [&](int r) {   // middle column, to the pixel downstream (south)
            return (d[static_cast<std::size_t>(r) * nc + 1] - d[static_cast<std::size_t>(r + 1) * nc + 1]) / dx;
        };
        const int contact = nr / 2;
        const double S1 = slope_at(contact - 3), S2 = slope_at(contact + 2);
        const double S1e = std::pow(U / K1, 1.0 / n), S2e = std::pow(U / K2, 1.0 / n);
        const float zmax = *std::max_element(d.begin(), d.end());
        const double tol = 4.0 * ulp(zmax) / (S1 * dx) + 1e-6 + 1e-3 * U * dt / (n * S1 * dx);
        check("A1 strip n=" + fmt(n) + ": reached steady state", steps > 0,
              std::to_string(steps) + " steps");
        check("A1 strip n=" + fmt(n) + ": S_i = (U/K_i)^(1/n) either side of the contact and "
              "S1/S2 = (K2/K1)^(1/n)",
              std::fabs(S1 / S1e - 1.0) <= tol && std::fabs(S2 / S2e - 1.0) <= tol &&
                  std::fabs((S1 / S2) / std::pow(K2 / K1, 1.0 / n) - 1.0) <= 2.0 * tol && tol < 0.01,
              "S1 " + fmt(S1) + " (" + fmt(S1e) + "), S2 " + fmt(S2) + " (" + fmt(S2e) +
                  "), tolerance " + fmt(tol));
    }

    // Per-pixel steady identity on a D8 surface: K_i a^m S^n = U.
    auto d8_identity = [&](const lsdtt_raster_cpu& z, const lsdtt_lem_options& o,
                           const std::vector<float>& f, double min_area_px, int& n_out,
                           double& worst) {
        const lsdtt_flowinfo fl = lem_d8_flow(z, o.boundary_conditions);
        const std::vector<double> area = fl.accumulate();
        const double dx = z.get_meta().cell_size;
        const std::vector<double> d = post_incision_surface(z, o, f);
        std::vector<std::uint8_t> fixed(d.size(), 0);
        for (std::int64_t b : build_boundary_flats(z, o.boundary_conditions))
            fixed[static_cast<std::size_t>(b)] = 1;
        const float zmax = static_cast<float>(*std::max_element(d.begin(), d.end()));
        bool ok = true;
        n_out = 0;
        worst = 0.0;
        for (std::int64_t i = 0; i < static_cast<std::int64_t>(f.size()); ++i)
        {
            if (!fl.is_valid_cell(i) || fl.is_outlet_unchecked(i)) continue;
            if (fixed[static_cast<std::size_t>(i)]) continue;
            if (area[static_cast<std::size_t>(i)] < min_area_px * dx * dx) continue;
            const std::int64_t r = fl.receiver_of_unchecked(i);
            const double dxe = (fl.receiver_dir_of_unchecked(i) & 1u) ? dx * std::sqrt(2.0) : dx;
            const double dz = d[static_cast<std::size_t>(i)] - d[static_cast<std::size_t>(r)];
            const double q = std::pow(f[static_cast<std::size_t>(i)] / o.uplift_rate, 1.0 / o.n_exponent) *
                             (dz / dxe) *
                             std::pow(area[static_cast<std::size_t>(i)] / dx,
                                      static_cast<double>(o.m_exponent) / o.n_exponent);
            const double tol = 4.0 * ulp(zmax) / dz + 1e-6 + 1e-3 / o.n_exponent;
            worst = std::max(worst, std::fabs(q - 1.0));
            if (std::fabs(q - 1.0) > tol || tol >= 0.01) ok = false;
            ++n_out;
        }
        return ok;
    };

    // A1b: 60x9 with K1 west and K2 east, m = 0.5, nnon.
    {
        const int nr = 60, nc = 9;
        lsdtt_raster_cpu z = strip(nr, nc, 10.0, 0.02f);
        lsdtt_lem_options o = fluvial_only(false, 1.0f, 1e-3f, 1e-3f, 500.0f, 0.5f);
        o.boundary_conditions = "nnon";
        std::vector<float> f(n_cells(z));
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
                f[static_cast<std::size_t>(r) * nc + c] = (c < nc / 2) ? 1e-3f : 2.5e-4f;
        o.erodibility_field = f;
        const int steps = run_to_steady(z, o, 40000);
        int n_px = 0;
        double worst = 0.0;
        const bool ok = d8_identity(z, o, f, 5.0, n_px, worst);
        check("A1b 60x9 K1 west / K2 east: steady, and K_i a^m S^n = U on trunk pixels",
              steps > 0 && ok && n_px >= 20,
              std::to_string(steps) + " steps, " + std::to_string(n_px) + " pixels, worst |q-1| " +
                  fmt(worst));
    }

    // A2: 30x30 red noise, two-region field, D8 and D-inf.
    for (bool dinf : {false, true})
        for (float n : {1.0f, 2.0f})
        {
            lsdtt_raster_cpu z0 = generate_red_noise_surface(30, 30, 10.0, 3, 0.0, 1.0, 6, 1);
            z0 = fill_pits(z0, 1e-4f, 1);
            lsdtt_lem_options o = fluvial_only(dinf, n, 1e-3f, 1e-3f, 500.0f, 0.5f);
            const std::vector<float> f = two_region_field(z0, 1e-3f, 2.5e-4f);
            o.erodibility_field = f;
            lsdtt_raster_cpu z = z0;
            const int steps = run_to_steady(z, o, 40000);
            const std::string tag = std::string("A2 ") + (dinf ? "D-inf" : "D8") + " n=" + fmt(n);
            if (!dinf)
            {
                int n_px = 0;
                double worst = 0.0;
                const bool ok = d8_identity(z, o, f, 1.0, n_px, worst);
                check(tag + ": steady, and K_i a^m S^n = U on every pixel",
                      steps > 0 && ok && n_px > 500,
                      std::to_string(steps) + " steps, " + std::to_string(n_px) +
                          " pixels, worst |q-1| " + fmt(worst));
                continue;
            }
            const lsdtt_flowgraph_cpu fg = lem_dinf_flow(z, "oooo");
            const std::vector<double> area = fg.accumulate();
            const std::vector<double> d = post_incision_surface(z, o, f);
            std::vector<std::uint8_t> fixed(d.size(), 0);
            for (std::int64_t b : build_boundary_flats(z, "oooo")) fixed[static_cast<std::size_t>(b)] = 1;
            const double dx = 10.0;
            const float zmax = static_cast<float>(*std::max_element(d.begin(), d.end()));
            int one = 0, two = 0, bad = 0, clamped = 0;
            double worst = 0.0;
            for (std::int64_t i = 0; i < static_cast<std::int64_t>(f.size()); ++i)
            {
                if (fg.is_outlet_flat(i) || !fg.is_valid_flat(i)) continue;
                if (fixed[static_cast<std::size_t>(i)]) continue;
                const std::int64_t nrcv = fg.n_receivers(i);
                const float w0 = nrcv > 0 ? fg.receiver_weight_at(i, 0) : 0.0f;
                const float w1 = nrcv > 1 ? fg.receiver_weight_at(i, 1) : 0.0f;
                const bool only_diag = nrcv > 1 && w1 >= 1.0f - 1e-6f;
                const bool tw = nrcv > 1 && !(w0 >= 1.0f - 1e-6f) && !only_diag;
                const double a = area[static_cast<std::size_t>(i)] / dx;
                const double zi = d[static_cast<std::size_t>(i)];
                double lhs, dz;
                if (tw)
                {
                    const double ze1 = d[static_cast<std::size_t>(fg.receiver_index_at(i, 0))];
                    const double ze2 = d[static_cast<std::size_t>(fg.receiver_index_at(i, 1))];
                    if (zi <= std::min(ze1, ze2) + 2.0 * o.slope_clamp * dx) { ++clamped; continue; }
                    const double sbt = std::pow((ze1 - zi) / dx, 2.0) + std::pow((ze2 - ze1) / dx, 2.0);
                    lhs = f[static_cast<std::size_t>(i)] * std::pow(a, 0.5) * std::pow(sbt, 0.5 * n);
                    dz = zi - std::min(ze1, ze2);
                    ++two;
                }
                else
                {
                    const int k0 = only_diag ? 1 : 0;
                    const std::int64_t rr = fg.receiver_index_at(i, k0);
                    const double dxe = only_diag ? dx * std::sqrt(2.0) : dx;
                    dz = zi - d[static_cast<std::size_t>(rr)];
                    if (dz <= 0.0) { ++clamped; continue; }
                    lhs = f[static_cast<std::size_t>(i)] * std::pow(a, 0.5) * std::pow(dz / dxe, n);
                    ++one;
                }
                const double q = std::pow(lhs / o.uplift_rate, 1.0 / n);
                const double tol = 4.0 * ulp(zmax) / dz + 1e-6 + 1e-3 / n;
                worst = std::max(worst, std::fabs(q - 1.0));
                if (std::fabs(q - 1.0) > tol) ++bad;
            }
            check(tag + ": steady; single-receiver S^n and two-receiver sbt^(n/2) identities hold",
                  steps > 0 && bad == 0 && one + two > 500,
                  std::to_string(steps) + " steps, " + std::to_string(one) + " single, " +
                      std::to_string(two) + " two-receiver, " + std::to_string(clamped) +
                      " clamped (excluded), worst |q-1| " + fmt(worst));
        }

    // A3: snap with the field, then time-step; compare with the scalar run.
    //
    // The snap integrates chi on the routing of the surface it is given, and the
    // snapped surface routes differently (about 330 of 840 D8 receivers change
    // here, with scalar K too), so it is not a fixed point of the transient
    // model, and on this surface both runs then settle into a limit cycle with
    // max |dz| of a few cm per step that never meets 1e-3 U dt. So the check is
    // that the field snap is as close to steady as the scalar snap: their first
    // step after snapping moves the surface by the same p999 within 1%. The
    // 200-step values are printed as a record, not asserted (the K/4 region
    // settles four times more slowly).
    {
        lsdtt_raster_cpu z0 = generate_red_noise_surface(30, 30, 10.0, 3, 0.0, 1.0, 6, 1);
        z0 = fill_pits(z0, 1e-4f, 1);
        const float U = 1e-3f, Kr = 1e-3f;
        auto snap_and_run = [&](const std::vector<float>& field, float& p999_first) {
            lsdtt_fluvial_steady_state_options so;
            so.U = U;
            so.K = Kr;
            so.m = 0.5f;
            so.n = 1.0f;
            so.erodibility_field = field;
            const lsdtt_flowinfo fl = d8_flow(z0);
            lsdtt_raster_cpu z = compute_d8_fluvial_steady_state(fl, so).elevation;
            lsdtt_lem_options o = fluvial_only(false, 1.0f, Kr, U, 500.0f, 0.5f);
            o.erodibility_field = field;
            const std::vector<std::int64_t> flats = build_boundary_flats(z, "oooo");
            float p999 = 0.0f;
            for (int i = 0; i < 200; ++i)
            {
                auto r = run_lem_step(z, o, flats);
                z = std::move(r.first);
                p999 = r.second.p999_local_elevation_change;
                if (i == 0) p999_first = p999;
            }
            return p999;
        };
        float f_first = 0.0f, s_first = 0.0f;
        const float p_field = snap_and_run(two_region_field(z0, 1e-3f, 2.5e-4f), f_first);
        const float p_scalar = snap_and_run(std::vector<float>(n_cells(z0), Kr), s_first);
        check("A3 first step after a snap: p999 with the field within 1% of the scalar snap's",
              std::fabs(f_first - s_first) <= 0.01f * s_first,
              "first-step p999 " + fmt(f_first) + " m vs " + fmt(s_first) +
                  " m; after 200 steps (recorded, not asserted) " + fmt(p_field) + " m vs " +
                  fmt(p_scalar) + " m; 1e-3 U dt = " + fmt(1e-3f * U * 500.0f) + " m");
    }
}

void test_forcing_scaling()
{
    const lsdtt_raster_cpu dem = make_surface();
    const float Kp = 4e-5f;
    // F3 (i): a uniform base at K_param/2 with K_hist == K_param equals the scalar K_param/2 run.
    {
        std::vector<float> base(n_cells(dem), Kp / 2.0f), applied;
        scale_lem_erodibility_field(base, Kp, Kp, applied);
        lsdtt_lem_options a = base_options();
        a.erodibility = Kp / 2.0f;
        lsdtt_lem_options b = base_options();
        b.erodibility = Kp;
        b.erodibility_field = applied;
        check("F3 (i) uniform base K_param/2, forcing K = K_param: bitwise equal to the scalar run",
              same_bits(run_steps(dem, a, 3), run_steps(dem, b, 3)));
    }
    // F3 (ii) and (iii): a two-region base under a K step at t = 500.
    {
        const lsdtt_lem_forcing_history h = [] {
            lsdtt_lem_forcing_history x;
            lsdtt_lem_forcing_point p0, p1;
            p0.time_yr = 0.0f;   p0.erodibility = 4e-5f;
            p1.time_yr = 500.0f; p1.erodibility = 2e-5f;
            x.points = {p0, p1};
            return x;
        }();
        const std::vector<float> base = two_region_field(dem, 8e-5f, 2e-5f);
        std::vector<float> applied, previous;
        int changed_at = -1;
        bool formula = true;
        for (int step = 0; step < 10; ++step)
        {
            const float t = 100.0f * step;
            const float Kh = h.at(t).erodibility;
            scale_lem_erodibility_field(base, Kp, Kh, applied);
            for (std::size_t i = 0; i < base.size(); ++i)
                formula = formula && applied[i] == static_cast<float>(
                                         (static_cast<double>(base[i]) * Kh) / static_cast<double>(Kp));
            if (step > 0 && !same_bits(applied, previous) && changed_at < 0) changed_at = step;
            previous = applied;
        }
        check("F3 (ii) two-region field scaled per F2 at every step", formula);
        check("F3 (iii) a K step at t = 500 changes the applied field exactly then (step 5)",
              changed_at == 5, "changed at step " + std::to_string(changed_at));
    }
}

// ---------------------------------------------------------------------------
// G1-G16: model grids and input-raster preparation
// ---------------------------------------------------------------------------
std::string g_work_dir = ".";

lsdtt_raster_meta ids_meta(double dx, int epsg = 32632)
{
    lsdtt_raster_meta m;
    m.n_rows = static_cast<int>(std::lround(1000.0 / dx));
    m.n_cols = static_cast<int>(std::lround(800.0 / dx));
    m.cell_size = dx;
    m.x_min = 500000.0;
    m.y_min = 5000000.0;
    m.no_data_value = -9999.0;
    if (epsg > 0) m.set_crs_from_epsg(epsg);
    return m;
}

float id_of(int r, int c) { return static_cast<float>(1 + 1000 * r + c); }

// Writes a source raster of ids (1 + 1000 r + c, top-row r) and returns its path.
std::string write_ids(const std::string& name, const lsdtt_raster_meta& m)
{
    std::vector<float> v(static_cast<std::size_t>(m.n_rows) * m.n_cols);
    for (int r = 0; r < m.n_rows; ++r)
        for (int c = 0; c < m.n_cols; ++c) v[static_cast<std::size_t>(r) * m.n_cols + c] = id_of(r, c);
    const std::string path = g_work_dir + "/" + name;
    lsdtt_data_loader w;
    w.write_full(v, m, path, raster_output_format::gtiff);
    return path;
}

lsdtt_raster_meta synthetic_model(int nr, int nc, double dx)
{
    return make_lem_synthetic_grid(nr, nc, dx, "");
}

lsdtt_lem_input_raster_placement_options placement(lem_input_raster_placement p, double E = NAN,
                                                   double N = NAN)
{
    lsdtt_lem_input_raster_placement_options o;
    o.placement = p;
    o.clip_easting = E;
    o.clip_northing = N;
    return o;
}

// Every window pixel (r, c) holds id(r0 + r, c0 + c).
bool window_ids(const lsdtt_raster_cpu& r, int r0, int c0, std::string* first_bad = nullptr)
{
    const int nr = r.get_meta().n_rows, nc = r.get_meta().n_cols;
    for (int i = 0; i < nr; ++i)
        for (int j = 0; j < nc; ++j)
        {
            const float got = r.get_data()[static_cast<std::size_t>(i) * nc + j];
            if (got != id_of(r0 + i, c0 + j))
            {
                if (first_bad)
                    *first_bad = "(" + std::to_string(i) + ", " + std::to_string(j) + ") = " +
                                 fmt(got) + ", expected " + fmt(id_of(r0 + i, c0 + j));
                return false;
            }
        }
    return true;
}

void test_grid_modes()
{
    std::filesystem::create_directories(g_work_dir);
    const lsdtt_raster_meta S = ids_meta(10.0);
    const std::string src = write_ids("lemk_src_ids.tif", S);
    const std::string src4 = write_ids("lemk_src_ids_4m.tif", ids_meta(4.0));
    const std::string src_geo = write_ids("lemk_src_geo.tif", [] {
        lsdtt_raster_meta m;
        m.n_rows = 10; m.n_cols = 10; m.cell_size = 0.001; m.x_min = 9.0; m.y_min = 45.0;
        m.no_data_value = -9999.0; m.set_crs_from_epsg(4326); return m; }());
    const std::string src_ft = write_ids("lemk_src_ft.tif", [] {
        lsdtt_raster_meta m;
        m.n_rows = 10; m.n_cols = 10; m.cell_size = 30.0; m.x_min = 6000000.0; m.y_min = 2000000.0;
        m.no_data_value = -9999.0; m.set_crs_from_epsg(2227); return m; }());
    const std::string src_nocrs = write_ids("lemk_src_nocrs.tif", ids_meta(10.0, 0));
    const int R = S.n_rows;
    auto source_of = [](const std::string& path) {
        lsdtt_data_loader l;
        l.open_raster(path);
        return describe_lem_source_grid(l, "test");
    };

    // G1
    {
        const lsdtt_raster_meta g = make_lem_model_grid_from_raster(source_of(src), "erodibility");
        check("G1 mode 2: the model grid is the source's grid (same_grid + CRS)",
              g.n_rows == S.n_rows && g.n_cols == S.n_cols && g.cell_size == S.cell_size &&
                  g.x_min == S.x_min && g.y_min == S.y_min && crs_equal_semantically(g.crs_wkt, S.crs_wkt));
        check("G1 geographic source throws 'geographic'",
              throws_containing([&] { (void)make_lem_model_grid_from_raster(source_of(src_geo), "erodibility"); }, {"geographic"}));
        check("G1 US-foot source throws 'metre'",
              throws_containing([&] { (void)make_lem_model_grid_from_raster(source_of(src_ft), "erodibility"); }, {"metre"}));
        check("G1 source with no CRS throws 'CRS'",
              throws_containing([&] { (void)make_lem_model_grid_from_raster(source_of(src_nocrs), "erodibility"); }, {"CRS"}));
    }

    // G2
    {
        const lsdtt_raster_meta g = make_lem_model_grid_from_raster_dimensions(source_of(src), "");
        const lsdtt_raster_meta w = synthetic_model(S.n_rows, S.n_cols, S.cell_size);
        check("G2 mode 1 (a): rows, cols, cell size from the source; CRS and origin synthetic",
              g.n_rows == S.n_rows && g.n_cols == S.n_cols && g.cell_size == S.cell_size &&
                  g.x_min == w.x_min && g.y_min == w.y_min && crs_equal_semantically(g.crs_wkt, w.crs_wkt) &&
                  !crs_equal_semantically(g.crs_wkt, S.crs_wkt));
        const auto p = prepare_lem_input_raster(src, g, placement(lem_input_raster_placement::relabel_whole_raster), "erodibility");
        lsdtt_data_loader l;
        l.open_raster(src);
        check("G2 relabel_whole_raster: values bitwise equal to the source, note on the georeference",
              same_bits(p.raster.get_data(), l.read_full().data) && p.relabelled && !p.notes.empty());
    }

    const double dx = 10.0;
    auto px_centre_x = [&](int c) { return S.x_min + (c + 0.5) * dx; };
    auto px_centre_y = [&](int r) { return S.y_min + (R - r - 0.5) * dx; };

    // G3 centre clip
    for (int odd = 0; odd < 2; ++odd)
    {
        const int nr = 20 + odd, nc = 10 + odd;
        const lsdtt_raster_meta model = synthetic_model(nr, nc, dx);
        std::string bad;
        auto p = prepare_lem_input_raster(src, model, placement(lem_input_raster_placement::relabel_centre_clip,
                                                                px_centre_x(30), px_centre_y(40)), "erodibility");
        check("G3 centre clip " + std::to_string(nr) + "x" + std::to_string(nc) +
                  ": window (r, c) = id(40 - " + std::to_string(nr / 2) + " + r, 30 - " + std::to_string(nc / 2) + " + c)",
              window_ids(p.raster, 40 - nr / 2, 30 - nc / 2, &bad), bad);
    }
    {
        const lsdtt_raster_meta model = synthetic_model(20, 10, dx);
        std::string bad;
        auto p = prepare_lem_input_raster(src, model, placement(lem_input_raster_placement::relabel_centre_clip,
                                                                S.x_min + 30 * dx, px_centre_y(40)), "erodibility");
        check("G3 centre on a vertical grid line selects the pixel to the east",
              window_ids(p.raster, 30, 25, &bad), bad);
        auto q = prepare_lem_input_raster(src, model, placement(lem_input_raster_placement::relabel_centre_clip,
                                                                px_centre_x(30), S.y_min + (R - 40) * dx), "erodibility");
        check("G3 centre on a horizontal grid line selects the pixel to the north",
              window_ids(q.raster, 29, 25, &bad), bad);
    }

    // G4 lower-left clip
    {
        const int nr = 20, nc = 10, k = 50;   // corner at the line y_min + k dx
        const lsdtt_raster_meta model = synthetic_model(nr, nc, dx);
        const int top = R - 1 - k - (nr - 1);
        std::string bad;
        auto a = prepare_lem_input_raster(src, model, placement(lem_input_raster_placement::relabel_lower_left_clip,
                                                                S.x_min + 30 * dx, S.y_min + k * dx), "erodibility");
        check("G4 lower-left corner on a grid corner: bottom-left pixel is the one to its north-east",
              window_ids(a.raster, top, 30, &bad), bad);
        auto b = prepare_lem_input_raster(src, model, placement(lem_input_raster_placement::relabel_lower_left_clip,
                                                                S.x_min + 30.4 * dx, S.y_min + (k + 0.6) * dx), "erodibility");
        check("G4 lower-left corner inside a pixel: that pixel", window_ids(b.raster, top, 30, &bad), bad);
    }

    // G5 off the source
    {
        const lsdtt_raster_meta model = synthetic_model(20, 10, dx);
        struct off { const char* side; double E, N; };
        const off cases[] = {
            {"west", S.x_min + 2 * dx, px_centre_y(50)}, {"east", S.x_min + 78 * dx, px_centre_y(50)},
            {"north", px_centre_x(40), px_centre_y(3)}, {"south", px_centre_x(40), px_centre_y(96)},
            {"outside", S.x_min - 500.0, S.y_min - 500.0}};
        for (const auto& o : cases)
            check(std::string("G5 centre clip off the ") + o.side + ": throws 'runs off the source'",
                  throws_containing([&] { (void)prepare_lem_input_raster(src, model,
                                        placement(lem_input_raster_placement::relabel_centre_clip, o.E, o.N), "erodibility"); },
                                    {"runs off the source"}));
    }

    // G6 different cell size, resampled with near
    for (int cat = 0; cat < 2; ++cat)
    {
        const lsdtt_raster_meta S4 = ids_meta(4.0);
        const int nr = 20, nc = 10;
        const lsdtt_raster_meta model = synthetic_model(nr, nc, dx);
        lsdtt_lem_input_raster_placement_options o =
            placement(lem_input_raster_placement::relabel_lower_left_clip, S.x_min + 102.0, S.y_min + 202.0);
        o.categorical = cat != 0;
        if (cat) o.resampling_method = "bilinear";
        const auto p = prepare_lem_input_raster(src4, model, o, "erodibility");
        bool ok = p.resampled;
        for (int r = 0; r < nr && ok; ++r)
            for (int c = 0; c < nc && ok; ++c)
            {
                const double xc = S.x_min + 102.0 + (c + 0.5) * dx;
                const double yc = S.y_min + 202.0 + (nr - r - 0.5) * dx;
                const int sc = static_cast<int>(std::floor((xc - S4.x_min) / 4.0));
                const int sr = S4.n_rows - 1 - static_cast<int>(std::floor((yc - S4.y_min) / 4.0));
                ok = p.raster.get_data()[static_cast<std::size_t>(r) * nc + c] == id_of(sr, sc);
            }
        bool note = false;
        for (const auto& n : p.notes) note = note || (n.find("cell size 4 to 10") != std::string::npos);
        check(std::string("G6 4 m source on a 10 m model, ") + (cat ? "categorical with bilinear asked (near used)" : "near") +
                  ": each value is the source pixel containing the target centre; note gives both cell sizes",
              ok && note);
    }

    // G7 warp_to_model_grid, G8 as_georeferenced refusal
    {
        lsdtt_raster_meta shifted = S;
        shifted.x_min += 3.0;
        shifted.y_min += 3.0;
        std::string bad;
        const auto p = prepare_lem_input_raster(src, shifted, placement(lem_input_raster_placement::warp_to_model_grid), "erodibility");
        check("G7 warp_to_model_grid onto a 0.3-cell shift (near): ids of the containing pixels",
              p.warped && window_ids(p.raster, 0, 0, &bad), bad);
        lsdtt_raster_meta far = S;
        far.x_min += 1e6;
        check("G7 warp with no overlap throws",
              throws_containing([&] { (void)prepare_lem_input_raster(src, far, placement(lem_input_raster_placement::warp_to_model_grid), "erodibility"); }, {"overlap"}));
        check("G7 warp of a source with no CRS throws",
              throws_containing([&] { (void)prepare_lem_input_raster(src_nocrs, shifted, placement(lem_input_raster_placement::warp_to_model_grid), "erodibility"); }, {"CRS"}));
        check("G8 as_georeferenced onto the 0.3-cell shift throws naming input_raster_placement",
              throws_containing([&] { (void)prepare_lem_input_raster(src, shifted, placement(lem_input_raster_placement::as_georeferenced), "erodibility"); },
                                {"input_raster_placement"}));
        lsdtt_lem_input_raster_placement_options cat = placement(lem_input_raster_placement::warp_to_model_grid);
        cat.categorical = true;
        cat.resampling_method = "bilinear";
        check("G9 categorical with bilinear under warp_to_model_grid throws naming near and mode",
              throws_containing([&] { (void)prepare_lem_input_raster(src, shifted, cat, "lithology code raster"); }, {"near", "mode"}));
    }

    // G10
    for (auto pl : {lem_input_raster_placement::relabel_whole_raster, lem_input_raster_placement::relabel_centre_clip,
                    lem_input_raster_placement::relabel_lower_left_clip})
    {
        lsdtt_lem_input_raster_placement_options o = placement(pl, px_centre_x(40), px_centre_y(50));
        o.model_grid_is_synthetic = false;
        check("G10 " + lem_input_raster_placement_name(pl) + " on a DEM grid throws",
              throws_containing([&] { (void)prepare_lem_input_raster(src, S, o, "erodibility"); }, {"synthetic model grid"}));
    }

    // G11: every placement yields a raster that passes the strict check.
    {
        const lsdtt_raster_meta clipm = synthetic_model(20, 10, dx);
        lsdtt_raster_meta shifted = S;
        shifted.x_min += 3.0;
        const lsdtt_raster_meta whole = make_lem_model_grid_from_raster_dimensions(source_of(src), "");
        struct c { lem_input_raster_placement p; lsdtt_raster_meta m; std::string path; };
        const std::vector<c> cases = {
            {lem_input_raster_placement::as_georeferenced, S, src},
            {lem_input_raster_placement::warp_to_model_grid, shifted, src},
            {lem_input_raster_placement::relabel_whole_raster, whole, src},
            {lem_input_raster_placement::relabel_centre_clip, clipm, src},
            {lem_input_raster_placement::relabel_lower_left_clip, clipm, src4}};
        bool ok = true;
        for (const auto& x : cases)
        {
            const auto p = prepare_lem_input_raster(x.path, x.m, placement(x.p, px_centre_x(40), px_centre_y(50)), "erodibility");
            try { (void)check_lem_field_grid_alignment(p.raster.get_meta(), x.m, "erodibility"); }
            catch (const std::exception&) { ok = false; }
        }
        check("G11 every placement's prepared raster passes check_lem_field_grid_alignment", ok);
    }

    // G12
    {
        auto vsimem_count = [] {
            char** list = VSIReadDir("/vsimem/");
            const int n = CSLCount(list);
            CSLDestroy(list);
            return n;
        };
        const int before = vsimem_count();
        lsdtt_raster_meta shifted = S;
        shifted.x_min += 3.0;
        (void)warp_raster_to_grid_in_memory(src, shifted, lsdtt_warp_options{});
        check("G12 warp_raster_to_grid_in_memory leaves no /vsimem/ file behind", vsimem_count() == before);
    }

    // G13 (A9) public window functions
    {
        const lsdtt_lem_source_grid g = source_of(src);
        const auto w = compute_lem_clip_window(g, lem_input_raster_placement::relabel_centre_clip,
                                               px_centre_x(30), px_centre_y(40), 20, 10, dx, "test");
        check("G13 compute_lem_clip_window (G3 inputs): pixel-aligned at row 30, col 25",
              w.pixel_aligned && w.row0 == 30 && w.col0 == 25);
        const auto wl = compute_lem_clip_window(g, lem_input_raster_placement::relabel_lower_left_clip,
                                                S.x_min + 30 * dx, S.y_min + 50 * dx, 20, 10, dx, "test");
        check("G13 compute_lem_clip_window (G4 inputs): row R-1-50-19, col 30",
              wl.pixel_aligned && wl.row0 == R - 1 - 50 - 19 && wl.col0 == 30);
        lsdtt_lem_source_grid cube = g;
        cube.n_bands = 3;
        const auto wc = compute_lem_clip_window(cube, lem_input_raster_placement::relabel_centre_clip,
                                                px_centre_x(30), px_centre_y(40), 20, 10, dx, "test");
        check("G13 a 3-band cube description gives the same window",
              wc.row0 == w.row0 && wc.col0 == w.col0 && wc.pixel_aligned);
        const auto w6 = compute_lem_clip_window(source_of(src4), lem_input_raster_placement::relabel_lower_left_clip,
                                                S.x_min + 102.0, S.y_min + 202.0, 20, 10, dx, "test");
        check("G13 (G6 inputs): not pixel-aligned, extent from compute_lem_clip_extent",
              !w6.pixel_aligned && w6.extent.x_min == S.x_min + 102.0 && w6.extent.x_max == S.x_min + 102.0 + 100.0 &&
                  w6.extent.y_min == S.y_min + 202.0 && w6.extent.y_max == S.y_min + 202.0 + 200.0);
        const auto ec = compute_lem_clip_extent(lem_input_raster_placement::relabel_centre_clip, 1000.0, 2000.0, 20, 10, 10.0);
        const auto el = compute_lem_clip_extent(lem_input_raster_placement::relabel_lower_left_clip, 1000.0, 2000.0, 20, 10, 10.0);
        check("G13 compute_lem_clip_extent (b) and (c) equal the hand-computed extents",
              ec.x_min == 950.0 && ec.x_max == 1050.0 && ec.y_min == 1900.0 && ec.y_max == 2100.0 &&
                  el.x_min == 1000.0 && el.x_max == 1100.0 && el.y_min == 2000.0 && el.y_max == 2200.0);
        const lsdtt_lem_clip_extent se{0.0, 0.0, 100.0, 100.0};
        bool sides = true;
        for (int i = 0; i < 4; ++i)
        {
            lsdtt_lem_clip_extent e = se;
            if (i == 0) e.x_min -= 10.0;
            if (i == 1) e.x_max += 10.0;
            if (i == 2) e.y_min -= 10.0;
            if (i == 3) e.y_max += 10.0;
            sides = sides && throws_containing([&] { check_lem_clip_extent_within(e, se, "test"); }, {"runs off the source"});
        }
        bool equal_ok = true;
        try { check_lem_clip_extent_within(se, se, "test"); } catch (const std::exception&) { equal_ok = false; }
        check("G13 check_lem_clip_extent_within: one cell off on each side throws, an equal extent passes",
              sides && equal_ok);
    }

    // G14
    {
        const std::string rot = g_work_dir + "/lemk_rotated.tif";
        {
            GDALDriver* drv = GetGDALDriverManager()->GetDriverByName("GTiff");
            GDALDataset* ds = drv->Create(rot.c_str(), 5, 5, 1, GDT_Float32, nullptr);
            double gt[6] = {500000.0, 10.0, 0.5, 5000050.0, 0.0, -10.0};
            ds->SetGeoTransform(gt);
            GDALClose(ds);
        }
        check("G14 a rotated GeoTIFF given to prepare_lem_input_raster throws 'rotated'",
              throws_containing([&] { (void)prepare_lem_input_raster(rot, S, placement(lem_input_raster_placement::as_georeferenced), "erodibility"); }, {"rotated"}));
        lsdtt_raster_meta empty = S;
        empty.n_rows = 0;
        lsdtt_raster_meta zero = S;
        zero.cell_size = 0.0;
        check("G14 make_lem_source_grid with 0 rows or cell size 0 throws naming the label",
              throws_containing([&] { (void)make_lem_source_grid(empty, 1, "erodibility"); }, {"erodibility", "empty"}) &&
                  throws_containing([&] { (void)make_lem_source_grid(zero, 1, "erodibility"); }, {"erodibility", "cell size"}));
    }

    // G15 distortion
    {
        auto grid = [](int epsg, double x, double y, double dxg) {
            lsdtt_raster_meta m;
            m.n_rows = 100; m.n_cols = 100; m.cell_size = dxg; m.x_min = x; m.y_min = y;
            m.set_crs_from_epsg(epsg);
            return m;
        };
        // UTM 32N central meridian 9E; Web Mercator y for 48N and 2N; LAEA (EPSG:3035) near its centre.
        const lsdtt_raster_meta utm_ok = grid(32632, 495000.0, 4980000.0, 100.0);
        const lsdtt_raster_meta wm48 = grid(3857, 1000000.0, 6106854.0, 30.0);
        const lsdtt_raster_meta wm2 = grid(3857, 1000000.0, 222684.0, 30.0);
        const lsdtt_raster_meta utm_far = grid(32632, 1450000.0, 5000000.0, 30.0);
        const lsdtt_raster_meta laea = grid(3035, 4321000.0, 3210000.0, 100.0);
        bool ok = true;
        try { check_lem_projection_distortion(utm_ok, "erodibility"); } catch (const std::exception&) { ok = false; }
        check("G15 UTM 32N at its central meridian passes", ok, "scale error " + fmt(utm_ok.max_linear_scale_error()));
        std::string what;
        const bool wm_ok = throws_containing([&] { check_lem_projection_distortion(wm48, "erodibility"); },
                                             {"Web Mercator", "scale factor 1.4"}, &what);
        check("G15 EPSG:3857 at 48 N throws 'Web Mercator' with scale factor 1.4x", wm_ok, what);
        check("G15 EPSG:3857 at 2 N throws by name",
              throws_containing([&] { check_lem_projection_distortion(wm2, "erodibility"); }, {"Web Mercator"}));
        const bool far_ok = throws_containing([&] { check_lem_projection_distortion(utm_far, "erodibility"); },
                                              {"distorts distances", "%"}, &what);
        check("G15 UTM 32N ~12 degrees east of its central meridian throws 'distorts distances' with a percentage",
              far_ok, what + " (scale error " + fmt(utm_far.max_linear_scale_error()) + ")");
        ok = true;
        try { check_lem_projection_distortion(laea, "erodibility"); } catch (const std::exception&) { ok = false; }
        check("G15 an equal-area LAEA 10 km grid near its centre passes", ok, "scale error " + fmt(laea.max_linear_scale_error()));
        const std::string wm_path = write_ids("lemk_src_webmerc.tif", wm48);
        check("G15 mode 2 and compute_lem_clip_window both raise it",
              throws_containing([&] { (void)make_lem_model_grid_from_raster(source_of(wm_path), "erodibility"); }, {"Web Mercator"}) &&
                  throws_containing([&] { (void)compute_lem_clip_window(source_of(wm_path), lem_input_raster_placement::relabel_centre_clip,
                                                                         1001500.0, 6108354.0, 10, 10, 30.0, "erodibility"); }, {"Web Mercator"}));
    }

    // G16 footprint
    {
        const int nr = 12, nc = 10;
        lsdtt_raster_meta m;
        m.n_rows = nr; m.n_cols = nc; m.cell_size = 10.0; m.no_data_value = -9999.0;
        lsdtt_raster_cpu dem(std::vector<float>(nr * nc, 5.0f), std::vector<std::uint8_t>(nr * nc, 1), m);
        std::vector<std::uint8_t> src_valid(nr * nc, 1);
        int west = 0;
        for (int r = 0; r < nr && west < 20; ++r)          // an irregular staircase on the west edge
            for (int c = 0; c <= (r % 3) && west < 20; ++c) { src_valid[r * nc + c] = 0; ++west; }
        for (int r = 5; r < 7; ++r)
            for (int c = 5; c < 7; ++c) src_valid[r * nc + c] = 0;   // 2x2 interior hole
        const std::size_t removed = apply_lem_source_footprint(dem, src_valid, "erodibility");
        bool same = true, nd = true;
        for (int i = 0; i < nr * nc; ++i)
        {
            same = same && dem.get_mask()[i] == src_valid[i];
            if (!src_valid[i]) nd = nd && dem.get_data()[i] == -9999.0f;
        }
        check("G16 apply_lem_source_footprint returns 24; DEM mask equals the source mask; removed pixels hold nodata",
              removed == 24 && same && nd, std::to_string(removed) + " removed");
        lsdtt_raster_cpu dem2(std::vector<float>(nr * nc, 5.0f), std::vector<std::uint8_t>(nr * nc, 1), m);
        check("G16 an all-nodata source throws 'no valid pixels'",
              throws_containing([&] { (void)apply_lem_source_footprint(dem2, std::vector<std::uint8_t>(nr * nc, 0), "erodibility"); }, {"no valid pixels"}));
        check("G16 a size mismatch throws",
              throws_containing([&] { (void)apply_lem_source_footprint(dem2, std::vector<std::uint8_t>(nr * nc - 1, 1), "erodibility"); }, {"pixels"}));
        const std::vector<std::int64_t> flats = build_boundary_flats(dem, "oooo");
        std::vector<std::uint8_t> is_flat(nr * nc, 0);
        for (auto f : flats) is_flat[static_cast<std::size_t>(f)] = 1;
        bool all_adjacent = true;
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
            {
                if (!src_valid[r * nc + c]) continue;
                bool touches = false;
                for (int dr = -1; dr <= 1; ++dr)
                    for (int dc = -1; dc <= 1; ++dc)
                    {
                        const int rr = r + dr, cc = c + dc;
                        if (rr >= 0 && rr < nr && cc >= 0 && cc < nc && !src_valid[rr * nc + cc]) touches = true;
                    }
                if (touches && !is_flat[r * nc + c]) all_adjacent = false;
            }
        check("G16 with oooo every valid pixel 8-adjacent to the removed region is a boundary pixel", all_adjacent);
        // With "nooo", a pixel whose only removed neighbours are to its north is not fixed.
        lsdtt_raster_cpu dem3(std::vector<float>(nr * nc, 5.0f), std::vector<std::uint8_t>(nr * nc, 1), m);
        std::vector<std::uint8_t> hole(nr * nc, 1);
        hole[5 * nc + 5] = 0;
        (void)apply_lem_source_footprint(dem3, hole, "erodibility");
        const auto f3 = build_boundary_flats(dem3, "nooo");
        const bool south_of_hole_fixed = std::find(f3.begin(), f3.end(), 6 * nc + 5) != f3.end();
        check("G16 with nooo a pixel whose only removed neighbour is to its north is not fixed",
              !south_of_hole_fixed);
    }
}

// ---------------------------------------------------------------------------
// --write-fixtures: the CLI test fixtures (generated, never committed)
// ---------------------------------------------------------------------------
void write_raster(const std::string& path, std::vector<float> v, const lsdtt_raster_meta& m,
                  const std::vector<std::uint8_t>& valid = {})
{
    if (!valid.empty())
        for (std::size_t i = 0; i < v.size(); ++i)
            if (!valid[i]) v[i] = static_cast<float>(m.no_data_value);
    lsdtt_data_loader w;
    w.write_full(v, m, path, raster_output_format::gtiff);
}

void write_text(const std::string& path, const std::string& body) { std::ofstream(path) << body; }

int write_fixtures(const std::string& dir)
{
    std::filesystem::create_directories(dir);
    auto P = [&](const std::string& f) { return dir + "/" + f; };

    // The DEM grid: 40 rows x 20 cols at 10 m, UTM 32N.
    lsdtt_raster_meta dm;
    dm.n_rows = 40; dm.n_cols = 20; dm.cell_size = 10.0;
    dm.x_min = 500000.0; dm.y_min = 5000000.0; dm.no_data_value = -9999.0;
    dm.set_crs_from_epsg(32632);
    const int nr = dm.n_rows, nc = dm.n_cols;
    const std::size_t n = static_cast<std::size_t>(nr) * nc;
    std::vector<float> z(n);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
            z[static_cast<std::size_t>(r) * nc + c] = static_cast<float>(
                2.0 * (nr - 1 - r) + 3.0 * std::fabs(std::sin(c * 0.37)) +
                0.3 * std::sin(r * 1.1 + c * 0.7));
    write_raster(P("lemk_dem.tif"), z, dm);

    const float Kref = 1e-4f;
    write_raster(P("lemk_k_uniform.tif"), std::vector<float>(n, Kref), dm);
    std::vector<float> two(n), half(n);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
            two[static_cast<std::size_t>(r) * nc + c] = (r < nr / 2) ? 1e-3f : 2.5e-4f;
    for (std::size_t i = 0; i < n; ++i)
        half[i] = static_cast<float>((static_cast<double>(two[i]) * 5e-5f) / static_cast<double>(Kref));
    write_raster(P("lemk_k_two_region.tif"), two, dm);
    write_raster(P("lemk_k_two_region_half.tif"), half, dm);

    lsdtt_raster_meta shifted = dm;
    shifted.x_min += 5.0;
    write_raster(P("lemk_k_shifted.tif"), std::vector<float>(n, Kref), shifted);
    lsdtt_raster_meta offset = dm;
    offset.x_min += 3.0;
    offset.y_min += 3.0;
    std::vector<float> off_ids(n);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c) off_ids[static_cast<std::size_t>(r) * nc + c] = id_of(r, c) * 1e-7f;
    write_raster(P("lemk_k_offset.tif"), off_ids, offset);
    // Near warp of the 0.3-cell offset onto the DEM grid: every DEM centre lies
    // in the offset pixel with the same (r, c).
    write_raster(P("lemk_expect_warped.tif"), off_ids, dm);
    lsdtt_raster_meta wrong = dm;
    wrong.n_rows = 41;
    write_raster(P("lemk_k_wrong_dims.tif"), std::vector<float>(41 * 20, Kref), wrong);

    auto one_bad = [&](const std::string& name, float v, bool nodata) {
        std::vector<float> k(n, Kref);
        std::vector<std::uint8_t> valid;
        k[static_cast<std::size_t>(10) * nc + 7] = v;
        if (nodata) { valid.assign(n, 1); valid[static_cast<std::size_t>(10) * nc + 7] = 0; }
        write_raster(P(name), k, dm, valid);
    };
    one_bad("lemk_k_nodata.tif", Kref, true);
    one_bad("lemk_k_negative.tif", -1e-4f, false);
    one_bad("lemk_k_zero.tif", 0.0f, false);
    one_bad("lemk_k_huge.tif", 3e38f, false);

    // Placement sources: ids x 1e-6, 100 x 80 at 10 m (and 4 m), UTM 32N.
    auto ids_k = [](const lsdtt_raster_meta& m) {
        std::vector<float> v(static_cast<std::size_t>(m.n_rows) * m.n_cols);
        for (int r = 0; r < m.n_rows; ++r)
            for (int c = 0; c < m.n_cols; ++c)
                v[static_cast<std::size_t>(r) * m.n_cols + c] = id_of(r, c) * 1e-6f;
        return v;
    };
    const lsdtt_raster_meta S = ids_meta(10.0), S4 = ids_meta(4.0);
    const std::vector<float> sv = ids_k(S), sv4 = ids_k(S4);
    write_raster(P("lemk_k_src_ids.tif"), sv, S);
    write_raster(P("lemk_k_src_ids_4m.tif"), sv4, S4);
    {
        lsdtt_raster_meta g;
        g.n_rows = 10; g.n_cols = 10; g.cell_size = 0.001; g.x_min = 9.0; g.y_min = 45.0;
        g.no_data_value = -9999.0;
        g.set_crs_from_epsg(4326);
        write_raster(P("lemk_k_geographic.tif"), std::vector<float>(100, Kref), g);
    }
    {
        // Irregular footprint: a staircase along the west and south edges plus
        // an interior lake, 2000 of the 8000 pixels exactly (75 % valid).
        std::vector<std::uint8_t> valid(sv.size(), 1);
        int removed = 0;
        for (int r = 0; r < S.n_rows; ++r)
            for (int c = 0; c < 3 + (r % 7); ++c) { valid[static_cast<std::size_t>(r) * S.n_cols + c] = 0; ++removed; }
        for (int r = 90; r < S.n_rows; ++r)
            for (int c = 0; c < 40 + 5 * (r % 3); ++c)
                if (valid[static_cast<std::size_t>(r) * S.n_cols + c]) { valid[static_cast<std::size_t>(r) * S.n_cols + c] = 0; ++removed; }
        for (int r = 30; r < 80 && removed < 2000; ++r)
            for (int c = 30; c < 70 && removed < 2000; ++c)
                if (valid[static_cast<std::size_t>(r) * S.n_cols + c]) { valid[static_cast<std::size_t>(r) * S.n_cols + c] = 0; ++removed; }
        write_raster(P("lemk_k_src_footprint.tif"), sv, S, valid);
    }

    // Expected clips onto a synthetic 20 x 30 grid at 10 m.
    {
        const int mr = 20, mc = 30;
        std::vector<float> centre(mr * mc), ll(mr * mc), res(mr * mc);
        for (int r = 0; r < mr; ++r)
            for (int c = 0; c < mc; ++c)
            {
                centre[r * mc + c] = id_of(50 - mr / 2 + r, 40 - mc / 2 + c) * 1e-6f;   // centre of pixel (50, 40)
                ll[r * mc + c] = id_of(S.n_rows - 1 - 40 - (mr - 1) + r, 30 + c) * 1e-6f; // corner (30, 40) lines
                const double xc = S.x_min + 102.0 + (c + 0.5) * 10.0;
                const double yc = S.y_min + 202.0 + (mr - r - 0.5) * 10.0;
                const int sc = static_cast<int>(std::floor((xc - S4.x_min) / 4.0));
                const int sr = S4.n_rows - 1 - static_cast<int>(std::floor((yc - S4.y_min) / 4.0));
                res[r * mc + c] = id_of(sr, sc) * 1e-6f;
            }
        lsdtt_raster_meta em;
        em.n_rows = mr; em.n_cols = mc; em.cell_size = 10.0; em.no_data_value = -9999.0;
        write_raster(P("lemk_expect_centre_clip.tif"), centre, em);
        write_raster(P("lemk_expect_lower_left_clip.tif"), ll, em);
        write_raster(P("lemk_expect_resampled_clip.tif"), res, em);
    }

    // Channel masks: the middle column on the DEM grid, and the same shifted.
    std::vector<float> mask(n, 0.0f);
    for (int r = 0; r < nr; ++r) mask[static_cast<std::size_t>(r) * nc + nc / 2] = 1.0f;
    write_raster(P("lemk_mask_ondem.tif"), mask, dm);
    write_raster(P("lemk_mask_shifted.tif"), mask, shifted);

    // Forcing histories (uplift_rate written explicitly: an absent column is 0).
    write_text(P("lemk_forcing_K.csv"), "time_yr,uplift_rate,erodibility\n0,1e-4,1e-4\n1000,1e-4,5e-5\n");
    write_text(P("lemk_forcing_K_x10.csv"), "time_yr,uplift_rate,erodibility\n0,1e-4,1e-3\n");
    write_text(P("lemk_forcing_K_t0.csv"), "time_yr,uplift_rate,erodibility\n0,1e-4,5e-4\n");
    write_text(P("lemk_forcing_K_zero.csv"), "time_yr,uplift_rate,erodibility\n0,1e-4,1e-4\n500,1e-4,0\n");
    std::cout << "fixtures written to " << dir << "\n";
    return 0;
}

// ---------------------------------------------------------------------------
// --compare-rasters / --compare-raster-values
// ---------------------------------------------------------------------------
int compare_rasters(const std::string& a, const std::string& b, bool full)
{
    lsdtt_data_loader la, lb;
    try
    {
        la.open_raster(a);
        lb.open_raster(b);
    }
    catch (const std::exception& e)
    {
        std::cout << "FAIL " << a << " vs " << b << " -- " << e.what() << "\n";
        return 1;
    }
    const lsdtt_raster_meta ma = la.get_meta(), mb = lb.get_meta();
    std::vector<std::string> diffs;
    if (ma.n_rows != mb.n_rows || ma.n_cols != mb.n_cols)
        diffs.push_back("dims " + std::to_string(ma.n_rows) + "x" + std::to_string(ma.n_cols) +
                        " vs " + std::to_string(mb.n_rows) + "x" + std::to_string(mb.n_cols));
    if (full)
    {
        if (ma.cell_size != mb.cell_size) diffs.push_back("cell size");
        if (ma.x_min != mb.x_min || ma.y_min != mb.y_min) diffs.push_back("origin");
        if (!crs_equal_semantically(ma.crs_wkt, mb.crs_wkt)) diffs.push_back("CRS");
    }
    if (diffs.empty())
    {
        const lsdtt_raster_tile ta = la.read_full(), tb = lb.read_full();
        if (ta.mask != tb.mask) diffs.push_back("mask");
        std::size_t n_diff = 0, first = 0;
        for (std::size_t i = 0; i < ta.data.size(); ++i)
            if (ta.mask[i] && std::memcmp(&ta.data[i], &tb.data[i], sizeof(float)) != 0)
            {
                if (n_diff == 0) first = i;
                ++n_diff;
            }
        if (n_diff > 0)
            diffs.push_back(std::to_string(n_diff) + " data values differ, first at flat " +
                            std::to_string(first) + ": " + fmt(ta.data[first]) + " vs " +
                            fmt(tb.data[first]));
    }
    std::cout << (diffs.empty() ? "PASS " : "FAIL ") << a << " vs " << b
              << (full ? " (full)" : " (values)");
    for (const auto& d : diffs) std::cout << " -- " << d;
    std::cout << "\n";
    return diffs.empty() ? 0 : 1;
}

// fast: everything except the analytic steady states; analytic: those only.
int run_suite(bool fast, bool analytic)
{
    std::cout << "== lem_erodibility ==\n";
    if (!fast)
    {
        test_analytic();
        std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
        return g_failures == 0 ? 0 : 1;
    }
    test_kernel_bitwise();
    test_k5_d8_newton();
    test_k6_dinf_newton();
    test_units_and_calibration();
    test_steady_state();
    test_values_and_scaling();
    test_run_lem_step();
    test_forcing_scaling();
    test_grid_modes();
    if (analytic) test_analytic();
    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}

int print_raster_data_hashes(int argc, char* argv[])
{
    for (int i = 2; i < argc; ++i)
    {
        lsdtt_data_loader loader;
        loader.open_raster(argv[i]);
        const lsdtt_raster_tile t = loader.read_full();
        std::vector<float> both = t.data;
        for (std::uint8_t m : t.mask) both.push_back(static_cast<float>(m));
        std::cout << hex(fnv1a(both)) << " " << argv[i] << "\n";
    }
    return 0;
}

int expect_raster_epsg(int argc, char* argv[])
{
    if (argc < 4)
    {
        std::cerr << "usage: --expect-raster-epsg <raster> <epsg> [--expect-centre-near <tol_m>]\n";
        return 2;
    }
    const std::string path = argv[2];
    const int want = std::stoi(argv[3]);

    lsdtt_data_loader loader;
    loader.open_raster(path);
    const lsdtt_raster_meta meta = loader.get_meta();

    int got = 0;
    OGRSpatialReference srs;
    if (!meta.crs_wkt.empty() && srs.importFromWkt(meta.crs_wkt.c_str()) == OGRERR_NONE)
    {
        srs.AutoIdentifyEPSG();
        const char* code = srs.GetAuthorityCode(nullptr);
        if (code != nullptr) got = std::atoi(code);
    }
    bool ok = got == want;
    std::cout << (ok ? "PASS " : "FAIL ") << path << ": EPSG " << got << " (expected " << want
              << ")\n";

    if (argc >= 6 && std::string(argv[4]) == "--expect-centre-near")
    {
        const double tol = std::stod(argv[5]);
        lsdtt_raster_meta anchor;
        anchor.n_rows = meta.n_rows;
        anchor.n_cols = meta.n_cols;
        anchor.cell_size = meta.cell_size;
        anchor.set_default_synthetic_georeference();
        const double ex = anchor.x_min + 0.5 * anchor.n_cols * anchor.cell_size;
        const double ey = anchor.y_min + 0.5 * anchor.n_rows * anchor.cell_size;
        const double gx = meta.x_min + 0.5 * meta.n_cols * meta.cell_size;
        const double gy = meta.y_min + 0.5 * meta.n_rows * meta.cell_size;
        const double d = std::hypot(gx - ex, gy - ey);
        const bool c_ok = !anchor.crs_wkt.empty() && d <= tol;
        std::cout << (c_ok ? "PASS " : "FAIL ") << path << ": grid centre (" << std::fixed
                  << std::setprecision(3) << gx << ", " << gy << ") is " << d
                  << " m from the synthetic anchor (" << ex << ", " << ey << "), tolerance "
                  << tol << " m\n";
        ok = ok && c_ok;
    }
    return ok ? 0 : 1;
}

} // namespace

int main(int argc, char* argv[])
{
    const std::string mode = (argc > 1) ? argv[1] : "";
    if (mode == "--suite" && argc > 3) g_work_dir = argv[3];
    if (mode.empty()) return run_suite(true, true);
    if (mode == "--suite" && argc > 2 && std::string(argv[2]) == "fast") return run_suite(true, false);
    if (mode == "--suite" && argc > 2 && std::string(argv[2]) == "analytic")
        return run_suite(false, true);
    if (mode == "--print-scalar-hashes") return print_scalar_hashes();
    if (mode == "--print-raster-data-hash") return print_raster_data_hashes(argc, argv);
    if (mode == "--expect-raster-epsg") return expect_raster_epsg(argc, argv);
    if (mode == "--write-fixtures" && argc > 2) return write_fixtures(argv[2]);
    if (mode == "--compare-rasters" && argc > 3) return compare_rasters(argv[2], argv[3], true);
    if (mode == "--compare-raster-values" && argc > 3)
        return compare_rasters(argv[2], argv[3], false);
    std::cerr << "usage: lsdtt-test-lem-erodibility --print-scalar-hashes\n"
                 "       lsdtt-test-lem-erodibility --print-raster-data-hash <raster>...\n"
                 "       lsdtt-test-lem-erodibility --expect-raster-epsg <raster> <epsg>"
                 " [--expect-centre-near <tol_m>]\n"
                 "       lsdtt-test-lem-erodibility --write-fixtures <dir>\n"
                 "       lsdtt-test-lem-erodibility --compare-rasters|--compare-raster-values <a> <b>\n";
    return 2;
}
