// test_hydro_roughness.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Unit tests for src/cpu/lsdtt_hydro_roughness_cpu -- Phase 3 of the
// distributed-hydrology work.
//
// What is under test
// ------------------
//   the seven channel-n modes against their closed forms, including the
//   published Strickler and Jarrett coefficients and Ferguson's VPE asymptotics
//
//   the roughness PARTITION: that channel n comes from the bed and hillslope
//   and floodplain n come from land cover, and -- the assertion that matters
//   for scenario work -- that changing land cover leaves channel n untouched
//
//   downstream behaviour: n falls downstream where the mechanism is downstream
//   fining, and tracks steepness where the mechanism is ksn
//
// A NOTE ON WHAT "n DECLINES DOWNSTREAM" ACTUALLY CLAIMS
// ------------------------------------------------------
// The gate is `rigorous_channel_n_declines_downstream`. That is true of
// `drainage_area` and `grain_size`, where the mechanism is Sternberg fining, but
// it is NOT a prediction of `ksn` mode. In a graded profile ksn is roughly
// constant, so ksn-mode n is roughly constant too: that mode says n tracks
// STEEPNESS, not position. Testing it as a downstream decline would be testing a
// claim the physics does not make, so it is tested as "n rises with ksn"
// instead, which is what the mode is for.
//
// Usage: lsdtt-test-hydro-roughness
// Prints one PASS/FAIL line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "lsdtt_hydro_roughness_cpu.hpp"
#include "lsdtt_hydro_velocity_field_cpu.hpp"

namespace
{

int g_n_checks = 0;
int g_n_failed = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << std::endl;
}

bool close_rel(double a, double b, double tol)
{
    const double d = std::max(std::fabs(a), std::fabs(b));
    if (d < 1e-300) return std::fabs(a - b) < tol;
    return std::fabs(a - b) / d < tol;
}

// Rasters store float, so a value read back out of one and compared against a
// double literal can only agree to about 1e-7 relative. Two checks here were
// first written at 1e-9 and failed on exactly that -- 0.2f is 0.20000000298 --
// which is the same float/double trap the Phase 0 TWI work hit. Anything that
// has been through a raster uses this.
constexpr double k_float_tol = 1e-6;

constexpr double k_nodata = -9999.0;

lsdtt::lsdtt_raster_meta make_meta(int n_rows, int n_cols)
{
    lsdtt::lsdtt_raster_meta m;
    m.n_rows = n_rows; m.n_cols = n_cols; m.cell_size = 10.0;
    m.no_data_value = k_nodata;
    m.set_default_synthetic_georeference();
    return m;
}

lsdtt::lsdtt_raster_cpu make_raster(const lsdtt::lsdtt_raster_meta& meta,
                                    const std::vector<float>& v)
{
    lsdtt::lsdtt_array2d_f  d(meta.n_rows, meta.n_cols, static_cast<float>(k_nodata));
    lsdtt::lsdtt_array2d_u8 m(meta.n_rows, meta.n_cols, std::uint8_t(0));
    for (int r = 0; r < meta.n_rows; ++r)
        for (int c = 0; c < meta.n_cols; ++c)
        {
            const float x = v[static_cast<std::size_t>(r) * meta.n_cols + c];
            d(r, c) = x;
            m(r, c) = (x == static_cast<float>(k_nodata)) ? std::uint8_t(0)
                                                          : std::uint8_t(1);
        }
    return lsdtt::lsdtt_raster_cpu(std::move(d), std::move(m), meta);
}

lsdtt::lsdtt_raster_cpu uniform(const lsdtt::lsdtt_raster_meta& meta, float v)
{
    return make_raster(meta, std::vector<float>(
        static_cast<std::size_t>(meta.n_rows) * meta.n_cols, v));
}

// =============================================================================
// 1. The closed forms
// =============================================================================
void test_closed_forms()
{
    // Strickler: n = 0.047 D50^(1/6). A 50 mm gravel bed gives ~0.0285, which
    // is the right order for a gravel-bed river.
    const double n_gravel = lsdtt::strickler_roughness(0.05, 0.047);
    check("strickler: n = a D50^(1/6)",
          close_rel(n_gravel, 0.047 * std::pow(0.05, 1.0 / 6.0), 1e-12),
          "n(D50=50mm) = " + std::to_string(n_gravel));
    check("strickler: a 64-fold coarsening doubles n",
          close_rel(lsdtt::strickler_roughness(0.05 * 64.0, 0.047),
                    2.0 * n_gravel, 1e-12));
    check("strickler: non-positive D50 returns zero",
          lsdtt::strickler_roughness(0.0, 0.047) == 0.0);

    // Jarrett (1984): n = 0.39 S^0.38 R^(-0.16).
    const double n_j = lsdtt::jarrett_roughness(0.01, 0.5);
    check("jarrett: n = 0.39 S^0.38 R^(-0.16)",
          close_rel(n_j, 0.39 * std::pow(0.01, 0.38) * std::pow(0.5, -0.16), 1e-12),
          "n(S=0.01,R=0.5) = " + std::to_string(n_j));
    check("jarrett: n rises with slope",
          lsdtt::jarrett_roughness(0.03, 0.5) > n_j);
    check("jarrett: n falls as the channel deepens",
          lsdtt::jarrett_roughness(0.01, 2.0) < n_j);
}

// =============================================================================
// 2. Ferguson's variable power equation
// =============================================================================
void test_ferguson()
{
    const double a1 = 6.5, a2 = 2.5, d84 = 0.1, S = 0.01;

    // Asymptote 1: at large relative submergence the VPE collapses to
    // Manning-Strickler, (8/f)^(1/2) -> a1 (R/D84)^(1/6).
    {
        const double R  = d84 * 1.0e6;
        const double got = lsdtt::ferguson_sqrt_8_over_f(R, d84, a1, a2);
        const double want = a1 * std::pow(R / d84, 1.0 / 6.0);
        check("ferguson: large R/D84 tends to Manning-Strickler",
              close_rel(got, want, 1e-4),
              "got " + std::to_string(got) + " want " + std::to_string(want));
    }

    // Asymptote 2: at small relative submergence it tends to a2 (R/D84), the
    // roughness-layer regime where a constant Manning n is the wrong shape.
    {
        const double R = d84 * 1.0e-6;
        const double got = lsdtt::ferguson_sqrt_8_over_f(R, d84, a1, a2);
        check("ferguson: small R/D84 tends to the roughness-layer limit",
              close_rel(got, a2 * (R / d84), 1e-4),
              "got " + std::to_string(got));
    }

    // The equivalent n must reproduce the VPE velocity exactly through Manning.
    {
        const double R = 0.3;
        const double v_vpe = lsdtt::ferguson_vpe_velocity(R, d84, S, a1, a2);
        const double n_eq  = lsdtt::ferguson_equivalent_manning_n(R, d84, S, a1, a2);
        const double v_man = (1.0 / n_eq) * std::pow(R, 2.0 / 3.0) * std::sqrt(S);
        check("ferguson: the equivalent n reproduces the VPE velocity exactly",
              close_rel(v_man, v_vpe, 1e-12),
              "manning " + std::to_string(v_man) + " vpe " + std::to_string(v_vpe));
    }

    // And the equivalent n is NOT a constant: it rises steeply as the channel
    // shallows, which is the whole reason headwaters are slow.
    {
        const double n_deep    = lsdtt::ferguson_equivalent_manning_n(2.0, d84, S, a1, a2);
        const double n_shallow = lsdtt::ferguson_equivalent_manning_n(0.15, d84, S, a1, a2);
        check("ferguson: equivalent n rises sharply at low submergence",
              n_shallow > 1.5 * n_deep,
              "n(R=2m) = " + std::to_string(n_deep) +
              ", n(R=0.15m) = " + std::to_string(n_shallow));
    }
}

// =============================================================================
// 3. Each mode over a grid, against its closed form
// =============================================================================
void test_modes_on_grid()
{
    const int n = 5;
    const auto meta = make_meta(1, n);
    const auto area  = make_raster(meta, {1e4f, 1e5f, 1e6f, 1e7f, 1e8f});
    const auto slope = make_raster(meta, {0.10f, 0.05f, 0.02f, 0.01f, 0.005f});
    const auto d50   = make_raster(meta, {0.20f, 0.10f, 0.05f, 0.02f, 0.01f});
    const auto d84   = make_raster(meta, {0.40f, 0.20f, 0.10f, 0.04f, 0.02f});
    const auto rad   = make_raster(meta, {0.10f, 0.30f, 0.60f, 1.00f, 1.50f});

    lsdtt::roughness_field_inputs in;
    in.drainage_area    = &area;
    in.slope            = &slope;
    in.d50              = &d50;
    in.d84              = &d84;
    in.hydraulic_radius = &rad;

    lsdtt::channel_roughness_options o;
    o.channel_roughness = 0.035;
    o.n_min = 1e-4; o.n_max = 10.0;

    // constant
    o.mode = lsdtt::channel_roughness_mode::constant;
    {
        const auto r = lsdtt::compute_channel_roughness(in, meta, o);
        bool all = true;
        for (int i = 0; i < n; ++i)
            if (!close_rel(r.get_data_array()(0, i), 0.035, 1e-6)) all = false;
        check("grid: constant mode gives n_ref everywhere", all);
    }

    // drainage_area
    o.mode = lsdtt::channel_roughness_mode::drainage_area;
    o.area_exponent = 0.08; o.reference_area_m2 = 1e6;
    {
        const auto r = lsdtt::compute_channel_roughness(in, meta, o);
        bool ok = true;
        for (int i = 0; i < n; ++i)
        {
            const double a = area.get_data_array()(0, i);
            const double want = 0.035 * std::pow(a / 1e6, -0.08);
            if (!close_rel(r.get_data_array()(0, i), want, 1e-5)) ok = false;
        }
        check("grid: drainage_area matches n_ref (A/A_ref)^(-e_A)", ok);
        check("grid: at A = A_ref it returns n_ref exactly",
              close_rel(r.get_data_array()(0, 2), 0.035, 1e-6));
    }

    // constant is the e_A = 0 limit of drainage_area
    o.area_exponent = 0.0;
    {
        const auto r = lsdtt::compute_channel_roughness(in, meta, o);
        bool all = true;
        for (int i = 0; i < n; ++i)
            if (!close_rel(r.get_data_array()(0, i), 0.035, k_float_tol)) all = false;
        check("grid: constant is the zero-exponent limit of drainage_area", all);
    }
    o.area_exponent = 0.08;

    // grain_size
    o.mode = lsdtt::channel_roughness_mode::grain_size;
    {
        const auto r = lsdtt::compute_channel_roughness(in, meta, o);
        bool ok = true;
        for (int i = 0; i < n; ++i)
        {
            const double want = lsdtt::strickler_roughness(
                d50.get_data_array()(0, i), o.strickler_coefficient);
            if (!close_rel(r.get_data_array()(0, i), want, 1e-6)) ok = false;
        }
        check("grid: grain_size matches Strickler on the D50 raster", ok);
    }

    // jarrett
    o.mode = lsdtt::channel_roughness_mode::jarrett;
    {
        const auto r = lsdtt::compute_channel_roughness(in, meta, o);
        bool ok = true;
        for (int i = 0; i < n; ++i)
        {
            const double want = lsdtt::jarrett_roughness(
                std::max<double>(slope.get_data_array()(0, i), o.min_slope),
                rad.get_data_array()(0, i));
            if (!close_rel(r.get_data_array()(0, i), want, 1e-6)) ok = false;
        }
        check("grid: jarrett matches its closed form", ok);
    }

    // ferguson_vpe
    o.mode = lsdtt::channel_roughness_mode::ferguson_vpe;
    {
        const auto r = lsdtt::compute_channel_roughness(in, meta, o);
        bool ok = true;
        for (int i = 0; i < n; ++i)
        {
            const double want = lsdtt::ferguson_equivalent_manning_n(
                rad.get_data_array()(0, i), d84.get_data_array()(0, i),
                std::max<double>(slope.get_data_array()(0, i), o.min_slope),
                o.ferguson_a1, o.ferguson_a2);
            if (!close_rel(r.get_data_array()(0, i), want, 1e-6)) ok = false;
        }
        check("grid: ferguson_vpe matches its equivalent-n form", ok);
    }

    // slope_area: the median reference centres the field on n_ref
    o.mode = lsdtt::channel_roughness_mode::slope_area;
    o.reference_ksn = 0.0;
    {
        const auto r = lsdtt::compute_channel_roughness(in, meta, o);
        std::vector<double> v;
        for (int i = 0; i < n; ++i) v.push_back(r.get_data_array()(0, i));
        std::sort(v.begin(), v.end());
        check("grid: slope_area centres on n_ref via the median ksn",
              close_rel(v[v.size() / 2], 0.035, 1e-6),
              "median n = " + std::to_string(v[v.size() / 2]));
    }

    // clamping
    o.mode = lsdtt::channel_roughness_mode::constant;
    o.channel_roughness = 5.0; o.n_max = 0.2;
    {
        const auto r = lsdtt::compute_channel_roughness(in, meta, o);
        check("grid: n is clamped to n_max",
              close_rel(r.get_data_array()(0, 0), 0.2, k_float_tol));
    }
}

// =============================================================================
// 4. Downstream behaviour -- the Phase 3 gate
// =============================================================================
void test_downstream_behaviour()
{
    const int n = 8;
    const auto meta = make_meta(1, n);

    std::vector<float> a(n), d50(n), ksn(n);
    for (int i = 0; i < n; ++i)
    {
        const double A = 1e3 * std::pow(10.0, 4.0 * i / (n - 1));   // 1e3 -> 1e7
        a[static_cast<std::size_t>(i)] = static_cast<float>(A);
        // Sternberg fining: D50 ~ A^-0.4
        d50[static_cast<std::size_t>(i)] =
            static_cast<float>(0.2 * std::pow(A / 1e3, -0.4));
        // A graded profile: ksn roughly constant, with a steep reach in the middle
        ksn[static_cast<std::size_t>(i)] = (i == 3) ? 40.0f : 20.0f;
    }
    const auto area   = make_raster(meta, a);
    const auto d50_r  = make_raster(meta, d50);
    const auto ksn_r  = make_raster(meta, ksn);

    lsdtt::roughness_field_inputs in;
    in.drainage_area = &area;
    in.d50 = &d50_r;
    in.ksn = &ksn_r;

    lsdtt::channel_roughness_options o;
    o.channel_roughness = 0.035;
    o.n_min = 1e-4; o.n_max = 10.0;

    // drainage_area: strictly falling downstream
    o.mode = lsdtt::channel_roughness_mode::drainage_area;
    {
        const auto r = lsdtt::compute_channel_roughness(in, meta, o);
        bool falling = true;
        for (int i = 1; i < n; ++i)
            if (!(r.get_data_array()(0, i) < r.get_data_array()(0, i - 1)))
                falling = false;
        check("downstream: drainage_area n falls monotonically downstream",
              falling,
              "head " + std::to_string(r.get_data_array()(0, 0)) +
              " mouth " + std::to_string(r.get_data_array()(0, n - 1)));
    }

    // grain_size with Sternberg fining: also falling, and by the predicted
    // amount -- n ~ D50^(1/6) and D50 ~ A^-0.4 gives n ~ A^-0.0667.
    o.mode = lsdtt::channel_roughness_mode::grain_size;
    {
        const auto r = lsdtt::compute_channel_roughness(in, meta, o);
        bool falling = true;
        for (int i = 1; i < n; ++i)
            if (!(r.get_data_array()(0, i) < r.get_data_array()(0, i - 1)))
                falling = false;
        check("downstream: grain_size n falls with Sternberg fining", falling);

        const double ratio = r.get_data_array()(0, n - 1) / r.get_data_array()(0, 0);
        const double predicted = std::pow(1e4, -0.4 / 6.0);   // (A ratio)^(-0.4/6)
        check("downstream: and by the predicted D50^(1/6) amount",
              close_rel(ratio, predicted, 1e-5),
              "measured " + std::to_string(ratio) +
              " predicted " + std::to_string(predicted));
    }

    // ksn: n tracks STEEPNESS, not position. The steep reach must be the
    // roughest cell, and the graded reaches must all share one value.
    o.mode = lsdtt::channel_roughness_mode::ksn;
    o.reference_ksn = 20.0;
    {
        const auto r = lsdtt::compute_channel_roughness(in, meta, o);
        double steep = r.get_data_array()(0, 3);
        bool graded_flat = true;
        for (int i = 0; i < n; ++i)
        {
            if (i == 3) continue;
            if (!close_rel(r.get_data_array()(0, i), 0.035, 1e-6)) graded_flat = false;
            if (!(steep > r.get_data_array()(0, i))) graded_flat = false;
        }
        check("downstream: ksn n tracks steepness, not position",
              graded_flat,
              "steep reach n = " + std::to_string(steep) + " vs graded 0.035");
        check("downstream: the steep reach matches (ksn/ksn_ref)^e_k",
              close_rel(steep, 0.035 * std::pow(2.0, o.ksn_exponent), 1e-6));
    }
}

// =============================================================================
// 5. The partition -- the assertion that matters for scenarios
// =============================================================================
void test_partition()
{
    const int n = 6;
    const auto meta = make_meta(1, n);
    //                                  hill  hill  chan  chan  flood flood
    const auto chan  = make_raster(meta, {0.f,  0.f,  1.f,  1.f,  0.f,  0.f});
    const auto surf  = make_raster(meta, {0.f,  0.f,  0.f,  0.f,  1.f,  1.f});
    const auto area  = uniform(meta, 1.0e6f);
    const auto lc_a  = uniform(meta, 0.060f);      // grazing
    const auto lc_b  = uniform(meta, 0.300f);      // woodland

    lsdtt::roughness_partition_options o;
    o.channel.mode = lsdtt::channel_roughness_mode::constant;
    o.channel.channel_roughness = 0.035;
    o.channel.n_min = 1e-4; o.channel.n_max = 10.0;
    o.floodplain_class = 1.0;

    lsdtt::roughness_field_inputs in;
    in.drainage_area = &area;
    in.channel_extent  = &chan;
    in.surface_class = &surf;

    in.landcover_roughness = &lc_a;
    const auto ra = lsdtt::compute_roughness_field(in, meta, o);

    in.landcover_roughness = &lc_b;
    const auto rb = lsdtt::compute_roughness_field(in, meta, o);

    check("partition: hillslope and floodplain take land-cover n",
          close_rel(ra.roughness.get_data_array()(0, 0), 0.060, 1e-6) &&
          close_rel(ra.roughness.get_data_array()(0, 4), 0.060, 1e-6));
    check("partition: channel takes the channel law, not land cover",
          close_rel(ra.roughness.get_data_array()(0, 2), 0.035, 1e-6));

    // THE assertion for scenario work: planting a field changes the bank, not
    // the bed.
    bool channel_unchanged = true;
    bool non_channel_moved = true;
    for (int i = 0; i < n; ++i)
    {
        const double x = ra.roughness.get_data_array()(0, i);
        const double y = rb.roughness.get_data_array()(0, i);
        const bool is_chan = (chan.get_data_array()(0, i) != 0.0f);
        if (is_chan  && !close_rel(x, y, k_float_tol)) channel_unchanged = false;
        if (!is_chan && close_rel(x, y, k_float_tol))  non_channel_moved = false;
    }
    check("partition: a land-cover change leaves channel n UNCHANGED",
          channel_unchanged);
    check("partition: and does change hillslope and floodplain n",
          non_channel_moved);

    check("partition: the regime counts are right",
          ra.stats.n_channel == 2 && ra.stats.n_floodplain == 2 &&
          ra.stats.n_hillslope == 2,
          "channel " + std::to_string(ra.stats.n_channel) +
          " floodplain " + std::to_string(ra.stats.n_floodplain) +
          " hillslope " + std::to_string(ra.stats.n_hillslope));

    // Jarrett's calibration envelope is counted, not silently accepted.
    {
        const auto slope = make_raster(meta, {0.01f, 0.01f, 0.5f, 0.01f, 0.01f, 0.01f});
        const auto rad   = uniform(meta, 0.5f);
        lsdtt::roughness_partition_options oj = o;
        oj.channel.mode = lsdtt::channel_roughness_mode::jarrett;
        lsdtt::roughness_field_inputs ij = in;
        ij.slope = &slope; ij.hydraulic_radius = &rad;
        const auto rj = lsdtt::compute_roughness_field(ij, meta, oj);
        check("partition: cells outside Jarrett's fitted range are counted",
              rj.stats.n_outside_jarrett_range == 1,
              std::to_string(rj.stats.n_outside_jarrett_range) + " of 2 channel cells");
    }
}

// =============================================================================
// 5b. End to end: the roughness field drives the velocity field
//
// This is the Phase 3 payoff and the reason the two modules share
// hydraulic_radius_from_area(). It also settles a design question that was left
// open: the velocity module does NOT need a separate `roughness_raster` mode.
// Composing compute_roughness_field() into compute_velocity_field()'s existing
// roughness input already does everything that mode was for, so it is not
// added -- a redundant enum value is a thing to keep correct forever.
// =============================================================================
void test_roughness_drives_velocity()
{
    const int n = 4;
    const auto meta = make_meta(1, n);
    //                                   hill  chan  chan  hill
    const auto chan  = make_raster(meta, {0.f,  1.f,  1.f,  0.f});
    const auto area  = uniform(meta, 1.0e6f);
    const auto slope = uniform(meta, 0.02f);
    const auto graze = uniform(meta, 0.060f);
    const auto wood  = uniform(meta, 0.300f);

    lsdtt::roughness_partition_options po;
    po.channel.mode = lsdtt::channel_roughness_mode::constant;
    po.channel.channel_roughness = 0.035;
    po.channel.n_min = 1e-4; po.channel.n_max = 10.0;

    lsdtt::roughness_field_inputs ri;
    ri.drainage_area = &area;
    ri.channel_extent  = &chan;

    lsdtt::hydro_velocity_options vo;
    vo.hydraulic_radius_coeff = 0.5;
    vo.v_min_metres_per_second = 1e-6;
    vo.v_max_metres_per_second = 1e6;

    auto velocity_with = [&](const lsdtt::lsdtt_raster_cpu& landcover) {
        ri.landcover_roughness = &landcover;
        const auto rf = lsdtt::compute_roughness_field(ri, meta, po);
        lsdtt::hydro_velocity_inputs vi;
        vi.drainage_area = &area;
        vi.slope         = &slope;
        vi.roughness     = &rf.roughness;
        return lsdtt::compute_velocity_field(vi, meta, vo);
    };

    const auto v_graze = velocity_with(graze);
    const auto v_wood  = velocity_with(wood);

    check("end to end: channel cells are faster than hillslope cells",
          v_graze.get_data_array()(0, 1) > v_graze.get_data_array()(0, 0),
          "hillslope " + std::to_string(v_graze.get_data_array()(0, 0)) +
          " channel " + std::to_string(v_graze.get_data_array()(0, 1)));

    // Manning is linear in 1/n, so the velocity ratio must be exactly the
    // inverse roughness ratio at equal A and S.
    check("end to end: the velocity ratio is exactly the inverse n ratio",
          close_rel(v_graze.get_data_array()(0, 1) / v_graze.get_data_array()(0, 0),
                    0.060 / 0.035, k_float_tol),
          "measured " + std::to_string(v_graze.get_data_array()(0, 1) /
                                       v_graze.get_data_array()(0, 0)));

    // The scenario result: planting woodland slows the hillslope and leaves the
    // channel alone. This is the whole mechanism by which an NFM scenario acts.
    check("end to end: a woodland scenario slows hillslope flow",
          v_wood.get_data_array()(0, 0) < v_graze.get_data_array()(0, 0),
          "grazing " + std::to_string(v_graze.get_data_array()(0, 0)) +
          " woodland " + std::to_string(v_wood.get_data_array()(0, 0)));
    check("end to end: and leaves channel velocity untouched",
          close_rel(v_wood.get_data_array()(0, 1),
                    v_graze.get_data_array()(0, 1), k_float_tol));
}

// =============================================================================
// 6. Validation
// =============================================================================
void test_validation()
{
    const auto meta = make_meta(2, 2);
    const auto area = uniform(meta, 1.0e6f);

    auto throws = [](auto&& fn) {
        try { fn(); } catch (const std::invalid_argument&) { return true; }
        catch (...) { return false; }
        return false;
    };

    lsdtt::roughness_field_inputs in;
    lsdtt::channel_roughness_options o;
    o.mode = lsdtt::channel_roughness_mode::grain_size;
    check("validation: grain_size without a D50 raster throws",
          throws([&]{ (void)lsdtt::compute_channel_roughness(in, meta, o); }));

    o.mode = lsdtt::channel_roughness_mode::ferguson_vpe;
    in.slope = &area; in.hydraulic_radius = &area;
    check("validation: ferguson_vpe without a D84 raster throws",
          throws([&]{ (void)lsdtt::compute_channel_roughness(in, meta, o); }));

    o.mode = lsdtt::channel_roughness_mode::jarrett;
    in.hydraulic_radius = nullptr;
    check("validation: jarrett without a hydraulic radius throws",
          throws([&]{ (void)lsdtt::compute_channel_roughness(in, meta, o); }));

    o.mode = lsdtt::channel_roughness_mode::drainage_area;
    in.drainage_area = &area;
    o.n_min = 5.0; o.n_max = 1.0;
    check("validation: n_min above n_max throws",
          throws([&]{ (void)lsdtt::compute_channel_roughness(in, meta, o); }));

    o.n_min = 0.005; o.n_max = 1.0;
    const auto wrong = uniform(make_meta(3, 3), 1.0f);
    in.drainage_area = &wrong;
    check("validation: a raster on the wrong grid throws",
          throws([&]{ (void)lsdtt::compute_channel_roughness(in, meta, o); }));
}

} // namespace

int main()
{
    std::cout << "lsdtt-test-hydro-roughness\n"
              << "--------------------------\n";
    try
    {
        test_closed_forms();
        test_ferguson();
        test_modes_on_grid();
        test_downstream_behaviour();
        test_partition();
        test_roughness_drives_velocity();
        test_validation();
    }
    catch (const std::exception& e)
    {
        std::cout << "FAIL  unexpected exception: " << e.what() << std::endl;
        ++g_n_failed; ++g_n_checks;
    }

    std::cout << "--------------------------\n"
              << g_n_checks - g_n_failed << " / " << g_n_checks << " checks passed\n"
              << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << std::endl;
    return g_n_failed == 0 ? 0 : 1;
}
