// test_channel_extent.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Unit tests for src/cpu/lsdtt_channel_extent_cpu -- Phase 6.
//
// What is under test
// ------------------
//   all five ladder rungs against their definitions
//   the emergent channel area, and that it agrees with where the velocity
//   closure's hydraulic radius actually crosses the sheet-flow floor
//   vector burn-and-snap, including rejection by area
//
//   THE PHASE 6 GATE: that the mask is a REFINEMENT, not a prerequisite --
//   supplying one, or changing which one, must not move the velocity field's
//   geometry at all. Only roughness may respond.
//
// WHY THE GATE IS PHRASED THAT WAY
// --------------------------------
// The gate name is `rigorous_channel_transition_continuity`, and Phase 1
// already measures the continuity itself (largest step-to-step velocity change
// 0.92% across six decades of area). What Phase 6 adds is the claim that
// introducing a MASK cannot break that continuity, because the mask is wired to
// roughness and not to geometry. That is the thing this phase could plausibly
// get wrong, so that is what is tested here.
//
// Usage: lsdtt-test-channel-mask
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "lsdtt_channel_extent_cpu.hpp"
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

constexpr double k_nodata = -9999.0;
constexpr double k_float_tol = 1e-6;   // rasters store float; see Phase 3

lsdtt::lsdtt_raster_meta make_meta(int nr, int nc, double cell = 10.0)
{
    lsdtt::lsdtt_raster_meta m;
    m.n_rows = nr; m.n_cols = nc; m.cell_size = cell;
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
// 1. The emergent channel head
// =============================================================================
void test_emergent()
{
    lsdtt::channel_extent_options o;
    o.hydraulic_radius_coeff      = 0.5;
    o.hydraulic_radius_exponent   = 0.4;
    o.reference_area_m2           = 1.0e6;
    o.sheet_flow_min_depth_metres = 0.01;

    const double A_head = lsdtt::emergent_channel_area_m2(o);
    const double want = 1.0e6 * std::pow(0.01 / 0.5, 1.0 / 0.4);
    check("emergent: A_head = A_ref (h_sheet/a_R)^(1/c_R)",
          close_rel(A_head, want, 1e-12),
          "A_head = " + std::to_string(A_head) + " m2");

    // The definition that matters: at exactly A_head the hydraulic radius the
    // VELOCITY closure uses equals the sheet-flow floor. If these two ever
    // disagree, the mask describes a different channel head from the one the
    // model actually has.
    lsdtt::hydro_velocity_options vo;
    vo.hydraulic_radius_coeff      = o.hydraulic_radius_coeff;
    vo.hydraulic_radius_exponent   = o.hydraulic_radius_exponent;
    vo.reference_area_m2           = o.reference_area_m2;
    vo.sheet_flow_min_depth_metres = o.sheet_flow_min_depth_metres;

    const double R_at   = lsdtt::hydraulic_radius_from_area(A_head, vo);
    const double R_below = lsdtt::hydraulic_radius_from_area(A_head * 0.5, vo);
    const double R_above = lsdtt::hydraulic_radius_from_area(A_head * 2.0, vo);

    check("emergent: at A_head the velocity closure's R equals h_sheet",
          close_rel(R_at, o.sheet_flow_min_depth_metres, 1e-9),
          "R = " + std::to_string(R_at));
    check("emergent: below it the floor binds",
          close_rel(R_below, o.sheet_flow_min_depth_metres, 1e-12));
    // Above the head the geometry binds, and by the EXACT amount the exponent
    // dictates: R(A) = h_sheet (A/A_head)^c_R. Asserting a round-number ratio
    // instead is how this check first failed -- 2^0.4 is 1.32, not the 1.5 I had
    // guessed at.
    check("emergent: above it R follows h_sheet (A/A_head)^c_R exactly",
          close_rel(R_above,
                    o.sheet_flow_min_depth_metres *
                        std::pow(2.0, o.hydraulic_radius_exponent), 1e-9),
          "R = " + std::to_string(R_above) + ", want " +
          std::to_string(o.sheet_flow_min_depth_metres *
                         std::pow(2.0, o.hydraulic_radius_exponent)));

    // Raising the sheet-flow depth pushes the channel head downstream.
    auto o2 = o; o2.sheet_flow_min_depth_metres = 0.05;
    check("emergent: a deeper sheet-flow floor moves the head downstream",
          lsdtt::emergent_channel_area_m2(o2) > A_head);

    bool threw = false;
    try { auto b = o; b.hydraulic_radius_coeff = 0.0;
          (void)lsdtt::emergent_channel_area_m2(b); }
    catch (const std::invalid_argument&) { threw = true; }
    check("emergent: a non-positive coefficient throws", threw);
}

// =============================================================================
// 2. The ladder
// =============================================================================
void test_ladder()
{
    const int n = 7;
    const auto meta = make_meta(1, n);
    //                                   1e2  1e3  1e4  1e5  1e6  1e7  1e8
    const auto area = make_raster(meta, {1e2f, 1e3f, 1e4f, 1e5f, 1e6f, 1e7f, 1e8f});
    const auto slope = make_raster(meta, {0.5f, 0.3f, 0.2f, 0.1f, 0.05f, 0.02f, 0.01f});

    lsdtt::channel_extent_options o;

    // none
    o.mode = lsdtt::channel_extent_mode::none;
    o.hydraulic_radius_coeff = 0.5; o.sheet_flow_min_depth_metres = 0.01;
    {
        const auto r = lsdtt::compute_channel_extent(area, nullptr, nullptr, meta, o);
        const double head = lsdtt::emergent_channel_area_m2(o);
        bool right = true;
        for (int i = 0; i < n; ++i)
        {
            const bool want = area.get_data_array()(0, i) >= head;
            if ((r.mask.get_data_array()(0, i) != 0.0f) != want) right = false;
        }
        check("ladder: `none` marks exactly the cells above the emergent head",
              right, "head at " + std::to_string(head) + " m2, " +
              std::to_string(r.n_channel) + " of " + std::to_string(r.n_valid));
    }

    // area_threshold
    o.mode = lsdtt::channel_extent_mode::area_threshold;
    o.channel_area_threshold_m2 = 1.0e5;
    {
        const auto r = lsdtt::compute_channel_extent(area, nullptr, nullptr, meta, o);
        check("ladder: `area_threshold` is a clean cut at the threshold",
              r.n_channel == 4 &&
              r.mask.get_data_array()(0, 2) == 0.0f &&
              r.mask.get_data_array()(0, 3) == 1.0f,
              std::to_string(r.n_channel) + " channel cells");
        check("ladder: and the channel fraction is reported",
              close_rel(r.channel_fraction(), 4.0 / 7.0, 1e-12));
    }

    // area_slope_threshold: steeper cells become channel at SMALLER area, which
    // is the whole point of Montgomery & Dietrich over a plain area cut.
    o.mode = lsdtt::channel_extent_mode::area_slope_threshold;
    o.area_slope_theta = 2.0;
    o.area_slope_threshold = 1.0e2;
    {
        const auto r = lsdtt::compute_channel_extent(area, &slope, nullptr, meta, o);
        bool right = true;
        for (int i = 0; i < n; ++i)
        {
            const double A = area.get_data_array()(0, i);
            const double S = slope.get_data_array()(0, i);
            const bool want = A * S * S >= 1.0e2;
            if ((r.mask.get_data_array()(0, i) != 0.0f) != want) right = false;
        }
        check("ladder: `area_slope_threshold` matches A S^theta >= threshold",
              right, std::to_string(r.n_channel) + " channel cells");

        // A steep small-area cell must qualify where a gentle one of the same
        // area does not.
        const auto a2 = make_raster(make_meta(1, 2), {1e3f, 1e3f});
        const auto s2 = make_raster(make_meta(1, 2), {0.5f, 0.01f});
        const auto r2 = lsdtt::compute_channel_extent(a2, &s2, nullptr,
                                                    make_meta(1, 2), o);
        check("ladder: at equal area, the steeper cell channelises and the gentle one does not",
              r2.mask.get_data_array()(0, 0) == 1.0f &&
              r2.mask.get_data_array()(0, 1) == 0.0f);
    }

    // extraction: any non-zero value in the supplied raster counts
    o.mode = lsdtt::channel_extent_mode::extraction;
    {
        const auto sup = make_raster(meta, {0.f, 1.f, 0.f, 3.f, 0.f, 1.f, 0.f});
        const auto r = lsdtt::compute_channel_extent(area, nullptr, &sup, meta, o);
        check("ladder: `extraction` takes any non-zero value as channel",
              r.n_channel == 3 &&
              r.mask.get_data_array()(0, 3) == 1.0f &&
              r.mask.get_data_array()(0, 4) == 0.0f,
              std::to_string(r.n_channel) + " channel cells");
    }
}

// =============================================================================
// 3. THE GATE: the mask is a refinement, not a prerequisite
// =============================================================================
void test_mask_is_a_refinement()
{
    const int n = 40;
    const auto meta = make_meta(1, n);
    std::vector<float> a(n);
    for (int i = 0; i < n; ++i)
        a[static_cast<std::size_t>(i)] =
            static_cast<float>(std::pow(10.0, 2.0 + 6.0 * i / (n - 1)));
    const auto area  = make_raster(meta, a);
    const auto slope = uniform(meta, 0.02f);

    lsdtt::hydro_velocity_options vo;
    vo.hydraulic_radius_coeff      = 0.5;
    vo.sheet_flow_min_depth_metres = 0.01;
    vo.v_min_metres_per_second     = 1e-9;
    vo.v_max_metres_per_second     = 1e9;
    vo.roughness_default           = 0.05;

    lsdtt::hydro_velocity_inputs vi;
    vi.drainage_area = &area;
    vi.slope         = &slope;

    // No mask at all: the model must run and produce a continuous field.
    const auto v_nomask = lsdtt::compute_velocity_field(vi, meta, vo);

    double worst_jump = 0.0;
    for (int i = 1; i < n; ++i)
    {
        const double p = v_nomask.get_data_array()(0, i - 1);
        const double q = v_nomask.get_data_array()(0, i);
        worst_jump = std::max(worst_jump, std::fabs(q - p) / std::max(p, q));
    }
    check("gate: with NO mask the model runs and the velocity field is continuous",
          worst_jump < 0.5,
          "largest step-to-step change " + std::to_string(100.0 * worst_jump) + "%");

    // Now supply three different masks. The GEOMETRY must not move at all: the
    // mask is wired to roughness, never to the hydraulic radius.
    lsdtt::channel_extent_options mo;
    mo.hydraulic_radius_coeff      = vo.hydraulic_radius_coeff;
    mo.sheet_flow_min_depth_metres = vo.sheet_flow_min_depth_metres;

    auto velocity_with_mask = [&](lsdtt::channel_extent_mode mode, double thresh) {
        mo.mode = mode;
        mo.channel_area_threshold_m2 = thresh;
        const auto cm = lsdtt::compute_channel_extent(area, &slope, nullptr, meta, mo);

        lsdtt::roughness_partition_options po;
        po.channel.mode = lsdtt::channel_roughness_mode::constant;
        po.channel.channel_roughness = 0.05;     // SAME as the hillslope value,
        po.channel.n_min = 1e-4; po.channel.n_max = 10.0;
        po.landcover_roughness_default = 0.05;   // so only geometry can differ

        lsdtt::roughness_field_inputs ri;
        ri.drainage_area = &area;
        ri.channel_extent  = &cm.mask;
        const auto rf = lsdtt::compute_roughness_field(ri, meta, po);

        lsdtt::hydro_velocity_inputs vi2 = vi;
        vi2.roughness = &rf.roughness;
        return lsdtt::compute_velocity_field(vi2, meta, vo);
    };

    const auto v_a = velocity_with_mask(lsdtt::channel_extent_mode::none, 0.0);
    const auto v_b = velocity_with_mask(lsdtt::channel_extent_mode::area_threshold, 1e4);
    const auto v_c = velocity_with_mask(lsdtt::channel_extent_mode::area_threshold, 1e6);

    bool identical = true;
    double worst = 0.0;
    for (int i = 0; i < n; ++i)
    {
        const double x = v_nomask.get_data_array()(0, i);
        for (const auto* v : {&v_a, &v_b, &v_c})
        {
            const double y = v->get_data_array()(0, i);
            if (!close_rel(x, y, k_float_tol)) identical = false;
            worst = std::max(worst, std::fabs(x - y) / std::max(x, y));
        }
    }
    check("gate: which mask is used does NOT change the velocity geometry",
          identical,
          "largest difference across three masks and no mask: " +
          std::to_string(100.0 * worst) + "%");

    // And the mask DOES do its job when roughness differs between regimes, or
    // the previous check would be vacuous -- it would pass for a mask that was
    // simply ignored everywhere.
    {
        mo.mode = lsdtt::channel_extent_mode::area_threshold;
        mo.channel_area_threshold_m2 = 1e5;
        const auto cm = lsdtt::compute_channel_extent(area, &slope, nullptr, meta, mo);

        lsdtt::roughness_partition_options po;
        po.channel.mode = lsdtt::channel_roughness_mode::constant;
        po.channel.channel_roughness = 0.035;    // now DIFFERENT
        po.channel.n_min = 1e-4; po.channel.n_max = 10.0;
        po.landcover_roughness_default = 0.30;

        lsdtt::roughness_field_inputs ri;
        ri.drainage_area = &area;
        ri.channel_extent  = &cm.mask;
        const auto rf = lsdtt::compute_roughness_field(ri, meta, po);

        lsdtt::hydro_velocity_inputs vi2 = vi;
        vi2.roughness = &rf.roughness;
        const auto v = lsdtt::compute_velocity_field(vi2, meta, vo);

        bool differs = false;
        for (int i = 0; i < n; ++i)
            if (!close_rel(v.get_data_array()(0, i),
                           v_nomask.get_data_array()(0, i), 1e-3)) differs = true;
        check("gate: but the mask IS load-bearing when roughness differs",
              differs, "otherwise the check above would be vacuous");
    }
}

// =============================================================================
// 4. Vector burn and snap
// =============================================================================
void test_vector_mode(const std::string& out_dir)
{
    const int n = 9;
    const auto meta = make_meta(3, n);

    // Drainage area: a valley along the middle row.
    std::vector<float> a(static_cast<std::size_t>(3) * n, 1.0e3f);
    for (int c = 0; c < n; ++c)
        a[static_cast<std::size_t>(1) * n + c] = static_cast<float>(1.0e6 * (c + 1));
    const auto area = make_raster(meta, a);

    // A GeoJSON line along the TOP row, one cell off the modelled flow path --
    // which is the realistic case: the survey and the DEM disagree by a pixel.
    const lsdtt::lsdtt_raster_meta& m = area.get_meta();
    const double x0 = m.x_min + 0.5 * m.cell_size;
    const double y_top = m.y_max() - 0.5 * m.cell_size;
    const std::string path = out_dir + "/channel_line.geojson";
    {
        std::ofstream f(path);
        // The raster carries a projected CRS (set_default_synthetic_georeference
        // uses EPSG:32758, WGS 84 / UTM zone 58S). GeoJSON defaults to WGS84
        // lat/lon, so without an explicit CRS member GDAL reads these UTM
        // metres as degrees, PROJ rejects them ("Invalid latitude"), and the
        // burn silently produces nothing.
        f << "{\"type\":\"FeatureCollection\","
          << "\"crs\":{\"type\":\"name\",\"properties\":"
          << "{\"name\":\"urn:ogc:def:crs:EPSG::32758\"}},"
          << "\"features\":[{\"type\":\"Feature\","
          << "\"properties\":{},\"geometry\":{\"type\":\"LineString\","
          << "\"coordinates\":[";
        for (int c = 0; c < n; ++c)
        {
            if (c) f << ",";
            f << "[" << (x0 + c * m.cell_size) << "," << y_top << "]";
        }
        f << "]}}]}";
    }

    lsdtt::channel_extent_options o;
    o.mode = lsdtt::channel_extent_mode::vector;
    o.channel_vector_fname = path;

    // No snapping: the mask sits on the top row, where the line was drawn.
    o.snap_search_radius_pixels = 0;
    {
        const auto r = lsdtt::compute_channel_extent(area, nullptr, nullptr, meta, o);
        check("vector: the line burns in where it was drawn",
              r.n_burned > 0 && r.n_snapped == 0,
              std::to_string(r.n_burned) + " burned, " +
              std::to_string(r.n_snapped) + " snapped");
        int on_top = 0, on_mid = 0;
        for (int c = 0; c < n; ++c)
        {
            if (r.mask.get_data_array()(0, c) != 0.0f) ++on_top;
            if (r.mask.get_data_array()(1, c) != 0.0f) ++on_mid;
        }
        check("vector: without snapping it stays off the flow path",
              on_top > 0 && on_mid == 0,
              std::to_string(on_top) + " on the drawn row, " +
              std::to_string(on_mid) + " on the valley axis");
    }

    // With snapping: it must move onto the valley, where the water is.
    o.snap_search_radius_pixels = 2;
    {
        const auto r = lsdtt::compute_channel_extent(area, nullptr, nullptr, meta, o);
        int on_mid = 0;
        for (int c = 0; c < n; ++c)
            if (r.mask.get_data_array()(1, c) != 0.0f) ++on_mid;
        check("vector: snapping moves the burn onto the valley axis",
              on_mid > 0 && r.n_snapped > 0,
              std::to_string(on_mid) + " cells on the axis, " +
              std::to_string(r.n_snapped) + " snapped");
    }

    // The area gate drops a line that cannot find a real channel nearby.
    {
        auto o2 = o;
        o2.snap_min_area_m2 = 1.0e12;    // nothing can satisfy this
        const auto r = lsdtt::compute_channel_extent(area, nullptr, nullptr, meta, o2);
        check("vector: snap_min_area_m2 rejects burns with no channel nearby",
              r.n_channel == 0 && r.n_rejected_by_area == r.n_burned,
              std::to_string(r.n_rejected_by_area) + " of " +
              std::to_string(r.n_burned) + " rejected");
    }

    // A missing file must fail loudly, not silently produce an empty mask.
    {
        auto o3 = o;
        o3.channel_vector_fname = out_dir + "/does_not_exist.geojson";
        bool threw = false;
        try { (void)lsdtt::compute_channel_extent(area, nullptr, nullptr, meta, o3); }
        catch (const std::exception&) { threw = true; }
        check("vector: a missing vector file throws rather than silently emptying",
              threw);
    }
}

// =============================================================================
// 5. Validation
// =============================================================================
void test_validation()
{
    const auto meta = make_meta(3, 3);
    const auto area = uniform(meta, 1.0e6f);

    auto throws = [](auto&& fn) {
        try { fn(); } catch (const std::invalid_argument&) { return true; }
        catch (...) { return false; }
        return false;
    };

    lsdtt::channel_extent_options o;
    o.mode = lsdtt::channel_extent_mode::area_slope_threshold;
    check("validation: area_slope_threshold without slope throws",
          throws([&]{ (void)lsdtt::compute_channel_extent(area, nullptr, nullptr, meta, o); }));

    o.mode = lsdtt::channel_extent_mode::extraction;
    check("validation: extraction without a supplied mask throws",
          throws([&]{ (void)lsdtt::compute_channel_extent(area, nullptr, nullptr, meta, o); }));

    o.mode = lsdtt::channel_extent_mode::vector;
    check("validation: vector mode without a filename throws",
          throws([&]{ (void)lsdtt::compute_channel_extent(area, nullptr, nullptr, meta, o); }));

    o.mode = lsdtt::channel_extent_mode::area_threshold;
    o.channel_area_threshold_m2 = 0.0;
    check("validation: a non-positive area threshold throws",
          throws([&]{ (void)lsdtt::compute_channel_extent(area, nullptr, nullptr, meta, o); }));

    o.channel_area_threshold_m2 = 1e5;
    const auto wrong = uniform(make_meta(4, 4), 1.0f);
    check("validation: a raster on the wrong grid throws",
          throws([&]{ (void)lsdtt::compute_channel_extent(wrong, nullptr, nullptr, meta, o); }));
}

} // namespace

int main(int argc, char** argv)
{
    const std::string out_dir = (argc > 1) ? argv[1] : ".";

    std::cout << "lsdtt-test-channel-mask\n"
              << "-----------------------\n";
    try
    {
        test_emergent();
        test_ladder();
        test_mask_is_a_refinement();
        test_vector_mode(out_dir);
        test_validation();
    }
    catch (const std::exception& e)
    {
        std::cout << "FAIL  unexpected exception: " << e.what() << std::endl;
        ++g_n_failed; ++g_n_checks;
    }

    std::cout << "-----------------------\n"
              << g_n_checks - g_n_failed << " / " << g_n_checks << " checks passed\n"
              << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << std::endl;
    return g_n_failed == 0 ? 0 : 1;
}
