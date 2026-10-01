// test_subbasin_partition.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Tests for lsdtt_subbasin_partition_cpu (S-HYPE port, WP6 gate).
//
// Every DEM is built in memory and is small enough that its D8 routing, its
// contributing areas and its flow distances can be written out by hand. The
// expected values below are those hand calculations; the routing each one rests
// on is asserted first (the "audit" checks), so a change in D8 tie-breaking
// shows up as a routing failure rather than as a mysterious wrong number.
//
// Cell size is 10 m throughout. A cardinal step is 10 m, a diagonal step
// 10*sqrt(2) m.
//
// TOLERANCES
//   * Quantities accumulated in double from exactly representable inputs
//     (pixel counts, areas, means of integer elevations, trunk lengths summed
//     step by step) are checked to 1e-9.
//   * The gradient raster is stored as float, so per-pixel gradients are
//     checked to 1e-6 (float storage). On the ramp used here every gradient
//     happens to round to exactly 0.5f.
//   * loc_rivlen is built from the float flow-distance raster, so it carries
//     float rounding of FD (FD < 200 m here, float spacing < 2e-5 m); it is
//     checked to 1e-4 m.
//
// FIXTURES
//   A  planar ramp           gradient and slope_mean_percent, elev mean/std
//   B  V-shaped valley pair  two disjoint sub-basins, unassigned count,
//                            duplicate outlets, nodata in the DEM, classes
//   C  nested tree           user outlets and Strahler-order outlets, chains,
//                            rivlen/loc_rivlen with and without upstream inflow
//   D  diagonal valley       rivlen/loc_rivlen with diagonal steps, and how the
//                            channel threshold moves the trunk head
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_flow_base.hpp"
#include "lsdtt_junction_network_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_surface_metrics_cpu.hpp"
#include "lsdtt_subbasin_partition_cpu.hpp"

namespace
{

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

void check_near(const std::string& name, double got, double want, double tol)
{
    const bool ok = std::isfinite(got) && std::abs(got - want) <= tol;
    check(name, ok,
          "got " + std::to_string(got) + ", want " + std::to_string(want) +
          ", tol " + std::to_string(tol));
}

template <typename T>
std::string vec_str(const std::vector<T>& v)
{
    std::string s = "{";
    for (std::size_t i = 0; i < v.size(); ++i)
    {
        if (i) s += ",";
        s += std::to_string(v[i]);
    }
    return s + "}";
}

constexpr double k_dx    = 10.0;
constexpr double k_sqrt2 = 1.41421356237309504880;
constexpr double k_tight = 1e-9;
constexpr double k_fd    = 1e-4;

lsdtt::lsdtt_raster_cpu make_raster(
    int n_rows, int n_cols, const std::function<double(int, int)>& z)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows        = n_rows;
    meta.n_cols        = n_cols;
    meta.cell_size     = k_dx;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();

    const std::size_t n = static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_cols);
    std::vector<float>        data(n, 0.0f);
    std::vector<std::uint8_t> mask(n, 1u);
    for (int r = 0; r < n_rows; ++r)
        for (int c = 0; c < n_cols; ++c)
            data[static_cast<std::size_t>(r) * n_cols + c] = static_cast<float>(z(r, c));

    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

std::unique_ptr<lsdtt::lsdtt_flow_base> make_flow(const lsdtt::lsdtt_raster_cpu& dem)
{
    lsdtt::lsdtt_flowinfo_options opts;
    opts.boundary_conditions     = "oooo";
    opts.build_s_vector          = true;
    opts.build_topological_order = true;
    return lsdtt::make_flow_base(
        dem.get_data_array(), dem.get_mask_array(), dem.get_meta(), opts);
}

// Receiver of (r, c) is (rr, cc).
void check_receiver(const lsdtt::lsdtt_flow_base& flow, const std::string& tag,
                    int r, int c, int rr, int cc)
{
    const std::int64_t got  = flow.receiver_of(flow.flat_index(r, c));
    const std::int64_t want = flow.flat_index(rr, cc);
    check(tag + " receiver of (" + std::to_string(r) + "," + std::to_string(c) +
              ") is (" + std::to_string(rr) + "," + std::to_string(cc) + ")",
          got == want, "got flat " + std::to_string(got));
}

std::int64_t sum_of(const std::vector<std::int64_t>& v)
{
    std::int64_t s = 0;
    for (const std::int64_t x : v) s += x;
    return s;
}

// =============================================================================
// A. Planar ramp, 6 x 6:  z = 1000 - 3 r - 4 c   (metres)
//
// Per 10 m cell the surface drops 3 m southward and 4 m eastward, so the
// gradient is (0.3, 0.4) and its magnitude exactly 0.5 m/m: a 3-4-5 triangle.
// The centred 3x3 stencil is exact on a plane, so every interior pixel
// (rows 1..4, cols 1..4: 16 pixels) has gradient 0.5; edge pixels are nodata.
//   slope_mean_percent = 100 * 0.5 = 50
//
// Routing: interior pixels go SE (7 m over 14.14 m = 0.495 beats E 0.4 and
// S 0.3); bottom row goes E, right column goes S; (5,5) is the only base-level
// outlet, so one sub-basin holds all 36 pixels.
//
// Elevation: z = 1000 - (3 r + 4 c), r and c independent and uniform on 0..5.
//   mean = 1000 - 7 * 2.5 = 982.5
//   population var = (9 + 16) * var(0..5) = 25 * (36 - 1)/12 = 25 * 35/12
//   std = 5 * sqrt(35/12) = 8.539125638...
// The sample (n-1) std would be larger by sqrt(36/35); the check pins the
// population definition.
// =============================================================================
void test_ramp()
{
    std::cout << "\n--- A: planar ramp ---\n";
    const int nr = 6, nc = 6;
    const lsdtt::lsdtt_raster_cpu dem = make_raster(nr, nc, [](int r, int c) {
        return 1000.0 - 3.0 * r - 4.0 * c;
    });
    const auto flow = make_flow(dem);

    check_receiver(*flow, "A audit", 2, 2, 3, 3);
    check_receiver(*flow, "A audit", 5, 0, 5, 1);
    check_receiver(*flow, "A audit", 0, 5, 1, 5);

    const lsdtt::lsdtt_raster_cpu grad = lsdtt::compute_gradient_magnitude_3x3_square(dem);

    std::int64_t n_valid_grad = 0;
    bool all_half = true;
    double sum_grad = 0.0;
    for (std::size_t i = 0; i < grad.get_data().size(); ++i)
    {
        if (grad.get_mask()[i] == 0) continue;
        ++n_valid_grad;
        sum_grad += static_cast<double>(grad.get_data()[i]);
        if (std::abs(static_cast<double>(grad.get_data()[i]) - 0.5) > 1e-6) all_half = false;
    }
    check("A interior gradient pixels = 16", n_valid_grad == 16,
          "got " + std::to_string(n_valid_grad));
    check("A every interior gradient = 0.5 m/m (float storage, tol 1e-6)", all_half);

    const auto part = lsdtt::partition_incremental_subbasins(
        *flow, {flow->flat_index(5, 5)});
    check("A one sub-basin of 36 pixels, none unassigned",
          part.n_subbasins() == 1 && part.n_pixels[0] == 36 && part.n_unassigned == 0);

    const auto geo = lsdtt::compute_subbasin_geodata(part, dem, grad);
    check("A gradient pixels used = 16", geo.n_valid_gradient[0] == 16);

    const double mean_grad = sum_grad / static_cast<double>(n_valid_grad);
    check_near("A mean gradient = 0.5 (analytic)", mean_grad, 0.5, 1e-6);
    check_near("A slope_mean_percent = 50 (analytic)", geo.slope_mean_percent[0], 50.0, 1e-4);
    check_near("A slope_mean_percent = 100 x mean gradient of the raster",
               geo.slope_mean_percent[0], 100.0 * mean_grad, k_tight);
    check_near("A area = 36 * 100 m2", geo.area_m2[0], 3600.0, k_tight);
    check_near("A elev_mean = 982.5", geo.elev_mean_m[0], 982.5, k_tight);
    check_near("A elev_std (population) = 5 sqrt(35/12)",
               geo.elev_std_m[0], 5.0 * std::sqrt(35.0 / 12.0), 1e-9);
}

// =============================================================================
// B. V-shaped valley pair, 6 rows x 6 cols:
//      z = 10 * d(c) + (5 - r),   d = {1, 0, 1, 1, 0, 1} for c = 0..5
//
// Two valleys, cols 1 and 4, both running south (1 m per cell). Hillslope
// pixels step sideways into their valley: the cardinal drop is 10 m over 10 m
// (slope 1.0), the diagonal into the valley 11 m over 14.14 m (0.78). Cols 2
// and 3 are at equal height, so nothing crosses between the valleys. The
// valley pixels at the foot, (5,1) and (5,4), are the only base-level outlets.
//
//   west sub-basin = cols 0..2 = 18 pixels = 1800 m2
//   east sub-basin = cols 3..5 = 18 pixels = 1800 m2
//
// West elevations: col 0 and col 2 are 10 + (5 - r) = 15,14,13,12,11,10;
// col 1 is 5,4,3,2,1,0.
//   sum = 75 + 15 + 75 = 165, mean = 165/18 = 55/6
//   sum of squares = 955 + 55 + 955 = 1965
//   population var = 1965/18 - (55/6)^2 = 905/36,  std = sqrt(905)/6
// The east side is the mirror image and has the same values.
// =============================================================================
lsdtt::lsdtt_raster_cpu make_v_pair()
{
    return make_raster(6, 6, [](int r, int c) {
        static const int d[6] = {1, 0, 1, 1, 0, 1};
        return 10.0 * d[c] + (5 - r);
    });
}

void test_v_pair()
{
    std::cout << "\n--- B: V-shaped valley pair ---\n";
    const int nr = 6, nc = 6;
    lsdtt::lsdtt_raster_cpu dem = make_v_pair();
    const auto flow = make_flow(dem);

    check_receiver(*flow, "B audit", 2, 0, 2, 1);
    check_receiver(*flow, "B audit", 2, 2, 2, 1);
    check_receiver(*flow, "B audit", 2, 3, 2, 4);
    check_receiver(*flow, "B audit", 2, 1, 3, 1);
    check("B audit (5,1) and (5,4) are base-level outlets",
          flow->is_outlet(flow->flat_index(5, 1)) && flow->is_outlet(flow->flat_index(5, 4)));

    const std::int64_t o_w = flow->flat_index(5, 1);
    const std::int64_t o_e = flow->flat_index(5, 4);

    // ---- two disjoint sub-basins ----
    {
        const auto part = lsdtt::partition_incremental_subbasins(*flow, {o_w, o_e});

        bool ids_ok = true;
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
                if (part.subbasin_id[static_cast<std::size_t>(r * nc + c)] != (c <= 2 ? 0 : 1))
                    ids_ok = false;
        check("B every pixel assigned once, to the valley it drains to", ids_ok);
        check("B pixel counts 18 + 18, none unassigned",
              part.n_pixels == std::vector<std::int64_t>{18, 18} && part.n_unassigned == 0,
              vec_str(part.n_pixels));
        check("B both downstream links are -1",
              part.downstream == std::vector<std::int32_t>{-1, -1}, vec_str(part.downstream));

        const auto grad = lsdtt::compute_gradient_magnitude_3x3_square(dem);
        const auto geo  = lsdtt::compute_subbasin_geodata(part, dem, grad);
        check_near("B areas sum to the modelled area 36 * 100 m2",
                   geo.area_m2[0] + geo.area_m2[1], 3600.0, k_tight);
        check_near("B west area 1800 m2", geo.area_m2[0], 1800.0, k_tight);
        check_near("B west elev_mean = 55/6", geo.elev_mean_m[0], 55.0 / 6.0, k_tight);
        check_near("B west elev_std = sqrt(905)/6", geo.elev_std_m[0], std::sqrt(905.0) / 6.0, 1e-9);
        check_near("B east elev_mean = 55/6", geo.elev_mean_m[1], 55.0 / 6.0, k_tight);

        const auto id_raster = lsdtt::get_subbasin_id_raster(part);
        check("B id raster carries ids where assigned",
              id_raster.mask[0] == 1u && id_raster.data[0] == 0.0f &&
              id_raster.data[static_cast<std::size_t>(nc - 1)] == 1.0f);
    }

    // ---- an outlet part-way down the west valley ----
    // Outlet (3,1): the west sub-basin is rows 0..3 of cols 0..2 = 12 pixels.
    // Rows 4..5 of cols 0..2 (6 pixels) reach base level at (5,1) without
    // passing an outlet. (3,1)'s receiver (4,1) is unassigned, so maindown -1.
    {
        const auto part = lsdtt::partition_incremental_subbasins(
            *flow, {flow->flat_index(3, 1), o_e});
        check("B mid-valley outlet: counts 12 + 18",
              part.n_pixels == std::vector<std::int64_t>{12, 18}, vec_str(part.n_pixels));
        check("B mid-valley outlet: 6 unassigned pixels", part.n_unassigned == 6,
              "got " + std::to_string(part.n_unassigned));
        check("B mid-valley outlet: assigned + unassigned = valid",
              sum_of(part.n_pixels) + part.n_unassigned == part.n_valid_pixels);
        check("B mid-valley outlet: maindown -1 (receiver unassigned)",
              part.downstream[0] == -1);
        check("B unassigned pixel (5,0) has id -1",
              part.subbasin_id[static_cast<std::size_t>(5 * nc + 0)] == -1);
    }

    // ---- duplicate and invalid outlets ----
    {
        const auto part = lsdtt::partition_incremental_subbasins(*flow, {o_w, o_e, o_w});
        check("B duplicate outlet dropped (first kept), counted",
              part.n_subbasins() == 2 && part.n_duplicate_outlets == 1 &&
              part.outlet_flat == std::vector<std::int64_t>{o_w, o_e});

        bool threw = false;
        try { (void)lsdtt::partition_incremental_subbasins(*flow, {o_w, 36}); }
        catch (const std::invalid_argument&) { threw = true; }
        check("B out-of-range outlet throws invalid_argument", threw);
    }

    // ---- nodata in the DEM is excluded from elevation statistics ----
    // Mask (0,0) (z = 15) in the DEM handed to the statistics only. Area is
    // unchanged (the partition still holds 18 pixels); the mean uses 17:
    //   (165 - 15) / 17 = 150/17
    {
        const auto part = lsdtt::partition_incremental_subbasins(*flow, {o_w, o_e});
        lsdtt::lsdtt_raster_cpu holed = make_v_pair();
        holed.get_mask_array().get_vector()[0] = 0u;
        const auto grad = lsdtt::compute_gradient_magnitude_3x3_square(holed);
        const auto geo  = lsdtt::compute_subbasin_geodata(part, holed, grad);
        check("B nodata pixel excluded: 17 elevations used", geo.n_valid_elev[0] == 17);
        check_near("B nodata pixel excluded: mean = 150/17", geo.elev_mean_m[0], 150.0 / 17.0, k_tight);
        check_near("B nodata pixel does not reduce area", geo.area_m2[0], 1800.0, k_tight);
    }

    // ---- class fractions and elevation offsets ----
    // West (cols 0..2): class 1 = col 1 (the valley floor, z 5..0, mean 2.5);
    //   class 2 = cols 0,2 rows 0..2 (z 15,14,13, mean 14);
    //   class 3 = cols 0,2 rows 3..5 (z 12,11,10, mean 11).
    //   fractions 6/18 each; elev_mean 55/6
    //   dh1 = 2.5 - 55/6 = -20/3,  dh2 = 14 - 55/6 = 29/6,  dh3 = 11 - 55/6 = 11/6
    //   (6 dh1 + 6 dh2 + 6 dh3) / 18 = 0
    // East (cols 3..5): class 5 at (0,3) only (z 15); (5,5) has NO class
    //   (nodata, z 10); every other pixel class 1. Valid class pixels 17:
    //   fraction5 = 1/17, fraction1 = 16/17.
    //   elev_mean still over all 18 DEM pixels = 55/6.
    //   class 1 mean = (165 - 15 - 10)/16 = 8.75,  dh1 = 8.75 - 55/6 = -5/12
    //   dh5 = 15 - 55/6 = 35/6
    {
        const auto part = lsdtt::partition_incremental_subbasins(*flow, {o_w, o_e});
        lsdtt::lsdtt_raster_cpu cls = make_raster(nr, nc, [](int r, int c) {
            if (c <= 2) return (c == 1) ? 1.0 : (r <= 2 ? 2.0 : 3.0);
            return (r == 0 && c == 3) ? 5.0 : 1.0;
        });
        cls.get_mask_array().get_vector()[static_cast<std::size_t>(5 * nc + 5)] = 0u;

        const auto cf = lsdtt::compute_subbasin_class_fractions(part, cls, dem);

        check("B classes present = {1,2,3,5}",
              cf.class_codes == std::vector<std::int32_t>{1, 2, 3, 5}, vec_str(cf.class_codes));
        check("B offsets = {0,3,5}", cf.offset == std::vector<std::int64_t>{0, 3, 5},
              vec_str(cf.offset));

        if (cf.entries.size() == 5)
        {
            const auto& w1 = cf.entries[0];
            const auto& w2 = cf.entries[1];
            const auto& w3 = cf.entries[2];
            const auto& e1 = cf.entries[3];
            const auto& e5 = cf.entries[4];
            check("B west entries are classes 1,2,3 of sub-basin 0",
                  w1.class_code == 1 && w2.class_code == 2 && w3.class_code == 3 &&
                  w1.subbasin == 0 && w3.subbasin == 0);
            check_near("B west fraction 1 = 1/3", w1.fraction, 1.0 / 3.0, k_tight);
            check_near("B west fraction 2 = 1/3", w2.fraction, 1.0 / 3.0, k_tight);
            check_near("B west fraction 3 = 1/3", w3.fraction, 1.0 / 3.0, k_tight);
            check_near("B west dh1 = -20/3", w1.dh_m, -20.0 / 3.0, k_tight);
            check_near("B west dh2 = 29/6", w2.dh_m, 29.0 / 6.0, k_tight);
            check_near("B west dh3 = 11/6", w3.dh_m, 11.0 / 6.0, k_tight);
            check_near("B west sum fraction*dh = 0",
                       w1.fraction * w1.dh_m + w2.fraction * w2.dh_m + w3.fraction * w3.dh_m,
                       0.0, k_tight);

            check("B east entries are classes 1,5 of sub-basin 1",
                  e1.class_code == 1 && e5.class_code == 5 && e1.subbasin == 1);
            check("B east class-valid pixels = 17", cf.n_valid_class[1] == 17);
            check_near("B east fraction 1 = 16/17", e1.fraction, 16.0 / 17.0, k_tight);
            check_near("B east fraction 5 = 1/17", e5.fraction, 1.0 / 17.0, k_tight);
            check_near("B east fractions sum to 1", e1.fraction + e5.fraction, 1.0, k_tight);
            check_near("B east dh1 = -5/12", e1.dh_m, -5.0 / 12.0, k_tight);
            check_near("B east dh5 = 35/6", e5.dh_m, 35.0 / 6.0, k_tight);
        }
        else
        {
            check("B five (sub-basin, class) entries", false,
                  "got " + std::to_string(cf.entries.size()));
        }
    }
}

// =============================================================================
// C. Nested tree, 8 rows x 10 cols:
//      rows 1..7:  z = 10 |c - 5| + (7 - r)
//      row 0:      z = z(row 1) + 30
//
// The trunk runs south down col 5 from (1,5) to (7,5), the only base-level
// outlet. Every pixel in rows 1..7 steps sideways (cardinal, slope 1.0) to
// col 5. Row 0 is a 30 m step above row 1 and drains straight south (slope 3
// beats any diagonal, at most 2.83). Row 1 therefore carries two tributaries:
//
//   West, cols 0..4 of rows 0..1 (10 pixels), flowing east into (1,5).
//     contributing pixels along row 1 (include self): 2, 4, 6, 8, 10
//   East, cols 6..9 of rows 0..1 (8 pixels), flowing west into (1,5).
//     contributing pixels along row 1, from col 9: 2, 4, 6, 8
//
// Side pixels in rows 2..7 have at most 5 contributing pixels, so with a
// channel threshold of 6 the channel network is:
//   west tributary (1,2)-(1,4), source (1,2), order 1
//   east tributary (1,7)-(1,6), source (1,7), order 1
//   trunk (1,5)-(7,5), confluence at (1,5), order 2
// Junctions (row-major ids): (1,2), (1,5), (1,7), and (7,5) because base-level
// channel pixels are junctions.
//
// Junction basin outlets (pixel just above the next junction downstream):
//   (1,2) -> (1,4);  (1,7) -> (1,6);  (1,5) -> (6,5);  (7,5) -> (7,5) itself.
//
// USER OUTLETS (1,4), (1,6), (7,5):
//   sub-basin 0 (west)  10 px, maindown 2
//   sub-basin 1 (east)   8 px, maindown 2
//   sub-basin 2 (trunk) 80 - 18 = 62 px, maindown -1  (incremental area only)
//
// STRAHLER OUTLETS:
//   order 1 exactly: (1,4), (1,6)       -> 10 + 8 assigned, 62 unassigned
//   order 2 exactly: (6,5), (7,5)       -> (7,5) gets row 7 only (10 px, every
//                                          row-7 pixel steps sideways into it);
//                                          (6,5) gets the other 70; chain 70->10
//   orders 1..2:     all four           -> 10, 8, 52, 10; chain
//                                          west/east -> (6,5) -> (7,5) -> -1
// =============================================================================
lsdtt::lsdtt_raster_cpu make_tree()
{
    return make_raster(8, 10, [](int r, int c) {
        const double z1 = 10.0 * std::abs(c - 5) + (7 - std::max(r, 1));
        return (r == 0) ? z1 + 30.0 : z1;
    });
}

void test_nested_tree()
{
    std::cout << "\n--- C: nested three-outlet tree ---\n";
    const int nc = 10;
    const lsdtt::lsdtt_raster_cpu dem = make_tree();
    const auto flow = make_flow(dem);
    const lsdtt::lsdtt_flow_base& fb = *flow;

    auto F = [&](int r, int c) { return fb.flat_index(r, c); };

    check_receiver(fb, "C audit", 0, 0, 1, 0);
    check_receiver(fb, "C audit", 0, 5, 1, 5);
    check_receiver(fb, "C audit", 1, 4, 1, 5);
    check_receiver(fb, "C audit", 1, 6, 1, 5);
    check_receiver(fb, "C audit", 1, 5, 2, 5);
    check_receiver(fb, "C audit", 4, 0, 4, 1);
    check_receiver(fb, "C audit", 4, 9, 4, 8);
    check_receiver(fb, "C audit", 6, 4, 6, 5);
    check_receiver(fb, "C audit", 7, 4, 7, 5);
    check("C audit contributing pixels (1,4)=10, (1,6)=8, (7,5)=80",
          fb.s_vector_size(F(1, 4)) == 10 && fb.s_vector_size(F(1, 6)) == 8 &&
          fb.s_vector_size(F(7, 5)) == 80);

    constexpr int k_threshold = 6;
    lsdtt::lsdtt_junction_network_options jopts;
    lsdtt::lsdtt_junction_network_cpu jn;
    jn.create(fb, fb.get_sources_index_threshold(k_threshold), jopts);
    check("C audit junctions = 4", jn.n_junctions() == 4,
          "got " + std::to_string(jn.n_junctions()));

    // ---------------- user outlets ----------------
    const std::vector<std::int64_t> user_outlets = lsdtt::subbasin_outlets_from_user_flats(
        fb, jn, {F(1, 4), F(1, 6), F(7, 5)});
    check("C user outlets pass through the junction basin manager unchanged",
          user_outlets == std::vector<std::int64_t>{F(1, 4), F(1, 6), F(7, 5)},
          vec_str(user_outlets));

    const auto part = lsdtt::partition_incremental_subbasins(fb, user_outlets);
    check("C user: counts {10, 8, 62}",
          part.n_pixels == std::vector<std::int64_t>{10, 8, 62}, vec_str(part.n_pixels));
    check("C user: tributaries link to the trunk sub-basin, trunk to -1",
          part.downstream == std::vector<std::int32_t>{2, 2, -1}, vec_str(part.downstream));
    check("C user: none unassigned, counts sum to the 80 valid pixels",
          part.n_unassigned == 0 && sum_of(part.n_pixels) == 80);
    check("C user: upstream-to-downstream order {0, 1, 2}",
          part.upstream_to_downstream == std::vector<std::int32_t>{0, 1, 2},
          vec_str(part.upstream_to_downstream));

    bool ids_ok = true;
    for (int r = 0; r < 8; ++r)
        for (int c = 0; c < nc; ++c)
        {
            const int want = (r <= 1 && c <= 4) ? 0 : (r <= 1 && c >= 6) ? 1 : 2;
            if (part.subbasin_id[static_cast<std::size_t>(r * nc + c)] != want) ids_ok = false;
        }
    check("C user: every pixel has the hand-assigned sub-basin", ids_ok);

    {
        const auto grad = lsdtt::compute_gradient_magnitude_3x3_square(dem);
        const auto geo  = lsdtt::compute_subbasin_geodata(part, dem, grad);
        check_near("C user: areas sum to 80 * 100 m2",
                   geo.area_m2[0] + geo.area_m2[1] + geo.area_m2[2], 8000.0, k_tight);
        check_near("C user: trunk area is the incremental 6200 m2", geo.area_m2[2], 6200.0, k_tight);
    }

    // ---------------- Strahler-order outlets ----------------
    auto sorted = [](std::vector<std::int64_t> v) { std::sort(v.begin(), v.end()); return v; };

    {
        const auto o1 = lsdtt::subbasin_outlets_by_stream_order(fb, jn, 1, 1);
        check("C order 1: outlets {(1,4), (1,6)}",
              sorted(o1) == std::vector<std::int64_t>{F(1, 4), F(1, 6)}, vec_str(o1));
        const auto p1 = lsdtt::partition_incremental_subbasins(fb, sorted(o1));
        check("C order 1: counts {10, 8}, 62 unassigned",
              p1.n_pixels == std::vector<std::int64_t>{10, 8} && p1.n_unassigned == 62,
              vec_str(p1.n_pixels) + " unassigned " + std::to_string(p1.n_unassigned));
        check("C order 1: maindown -1 (receiver (1,5) unassigned)",
              p1.downstream == std::vector<std::int32_t>{-1, -1}, vec_str(p1.downstream));
    }
    {
        const auto o2 = lsdtt::subbasin_outlets_by_stream_order(fb, jn, 2, 2);
        check("C order 2: outlets {(6,5), (7,5)}",
              sorted(o2) == std::vector<std::int64_t>{F(6, 5), F(7, 5)}, vec_str(o2));
        const auto p2 = lsdtt::partition_incremental_subbasins(fb, sorted(o2));
        check("C order 2: counts {70, 10}, chain 0 -> 1 -> -1",
              p2.n_pixels == std::vector<std::int64_t>{70, 10} &&
              p2.downstream == std::vector<std::int32_t>{1, -1} && p2.n_unassigned == 0,
              vec_str(p2.n_pixels) + " " + vec_str(p2.downstream));
    }
    {
        const auto o12 = lsdtt::subbasin_outlets_by_stream_order(fb, jn, 1, 2);
        const std::vector<std::int64_t> want = {F(1, 4), F(1, 6), F(6, 5), F(7, 5)};
        check("C orders 1..2: all four outlets", sorted(o12) == want, vec_str(o12));
        const auto p = lsdtt::partition_incremental_subbasins(fb, sorted(o12));
        check("C orders 1..2: counts {10, 8, 52, 10}",
              p.n_pixels == std::vector<std::int64_t>{10, 8, 52, 10}, vec_str(p.n_pixels));
        check("C orders 1..2: chain {2, 2, 3, -1}",
              p.downstream == std::vector<std::int32_t>{2, 2, 3, -1}, vec_str(p.downstream));
        check("C orders 1..2: upstream-to-downstream {0, 1, 2, 3}",
              p.upstream_to_downstream == std::vector<std::int32_t>{0, 1, 2, 3},
              vec_str(p.upstream_to_downstream));
    }

    // ---------------- river lengths (user outlets, threshold 6) ----------------
    //
    // West (headwater). Walk up from (1,4): in-sub-basin donors (1,3) [8] and
    // (0,4) [1] -> (1,3); then (1,2) [6] -> (1,2); then (1,1) [4] and (0,2) [1]
    // are below 6 -> stop. p_in = (1,2). rivlen = 2 cardinal steps = 20 m.
    // Distance to the trunk (1,2)-(1,4):
    //   row 1: (1,0) 20, (1,1) 10, trunk 0,0,0                       = 30
    //   row 0: (0,0) 10+20, (0,1) 10+10, (0,2),(0,3),(0,4) 10 each    = 80
    //   loc_rivlen = 110 / 10 = 11 m
    //
    // East (headwater). From (1,6): (1,7) [6] -> (1,7); then (1,8) [4] stops.
    // p_in = (1,7), rivlen = 10 m. Distances:
    //   row 1: 0, 0, (1,8) 10, (1,9) 20                              = 30
    //   row 0: (0,6) 10, (0,7) 10, (0,8) 20, (0,9) 30                 = 70
    //   loc_rivlen = 100 / 8 = 12.5 m
    //
    // Trunk (one upstream inflow chosen: west, 10 px, beats east, 8 px).
    // p_in = receiver of (1,4) = (1,5). Trunk (1,5)..(7,5): 6 cardinal steps,
    // rivlen = 60 m. Distances over its 62 pixels:
    //   (0,5) 10; trunk pixels 0;
    //   rows 2..7, west side 50+40+30+20+10 = 150, east side 10+20+30+40 = 100,
    //   6 rows * 250 = 1500
    //   loc_rivlen = 1510 / 62 = 24.354838...
    //
    // The in-sub-basin walk agrees here: from (7,5) the largest in-sub-basin
    // donor is always the next trunk pixel up, and at (1,5) the only
    // in-sub-basin donor is (0,5) with 1 contributing pixel < 6, so it stops
    // exactly at the pixel that receives the west outlet.
    {
        const auto rl = lsdtt::compute_subbasin_river_lengths(part, fb, k_threshold);
        const auto fd = fb.get_distance_from_outlet_raster();
        auto FD = [&](std::int64_t f) { return static_cast<double>(fd.data[static_cast<std::size_t>(f)]); };

        check("C west trunk head is (1,2)", rl.trunk_top_flat[0] == F(1, 2));
        check("C east trunk head is (1,7)", rl.trunk_top_flat[1] == F(1, 7));
        check("C trunk sub-basin head is (1,5), the receiver of the west outlet",
              rl.trunk_top_flat[2] == F(1, 5) && fb.receiver_of(F(1, 4)) == F(1, 5));
        check("C trunk pixel counts {3, 2, 7}",
              rl.n_trunk_pixels == std::vector<std::int64_t>{3, 2, 7}, vec_str(rl.n_trunk_pixels));

        check_near("C west rivlen = 20 m", rl.rivlen_m[0], 20.0, k_tight);
        check_near("C east rivlen = 10 m", rl.rivlen_m[1], 10.0, k_tight);
        check_near("C trunk rivlen = 60 m", rl.rivlen_m[2], 60.0, k_tight);
        check_near("C trunk rivlen = FD(p_in) - FD(outlet)",
                   rl.rivlen_m[2], FD(F(1, 5)) - FD(F(7, 5)), k_fd);
        check_near("C west rivlen = FD(p_in) - FD(outlet)",
                   rl.rivlen_m[0], FD(F(1, 2)) - FD(F(1, 4)), k_fd);

        check_near("C west loc_rivlen = 11 m", rl.loc_rivlen_m[0], 11.0, k_fd);
        check_near("C east loc_rivlen = 12.5 m", rl.loc_rivlen_m[1], 12.5, k_fd);
        check_near("C trunk loc_rivlen = 1510/62 m", rl.loc_rivlen_m[2], 1510.0 / 62.0, k_fd);

        // The in-sub-basin walk stops at (1,5) too: no in-sub-basin donor of
        // (1,5) reaches the threshold.
        bool walk_stops = true;
        for (int dr = -1; dr <= 1; ++dr)
            for (int dc = -1; dc <= 1; ++dc)
            {
                if (dr == 0 && dc == 0) continue;
                const std::int64_t nb = F(1 + dr, 5 + dc);
                if (fb.receiver_of(nb) != F(1, 5)) continue;
                if (part.subbasin_id[static_cast<std::size_t>(nb)] != 2) continue;
                if (fb.s_vector_size(nb) >= static_cast<std::uint64_t>(k_threshold))
                    walk_stops = false;
            }
        check("C in-sub-basin largest-area walk also ends at (1,5)", walk_stops);
    }

    // Threshold moves the headwater trunk head. With threshold 3, (1,1) [4]
    // qualifies and (1,0) [2] does not: p_in = (1,1), rivlen = 30 m.
    {
        const auto rl = lsdtt::compute_subbasin_river_lengths(part, fb, 3);
        check("C threshold 3: west trunk head (1,1)", rl.trunk_top_flat[0] == F(1, 1));
        check_near("C threshold 3: west rivlen = 30 m", rl.rivlen_m[0], 30.0, k_tight);
        check("C threshold 3: trunk sub-basin head unchanged (inflow-defined)",
              rl.trunk_top_flat[2] == F(1, 5));
    }
}

// =============================================================================
// D. Diagonal valley, 5 x 5:  z = 10 |r - c| - (r + c)
//
// The valley floor is the diagonal r = c, descending SE to (4,4), the only
// base-level outlet. With d = |r - c|:
//   d = 0  on the floor: SE to the next floor pixel (diagonal step).
//   d = 1  cardinal step onto the floor (drop 11 m over 10 m = 1.1, against
//          the diagonal along-slope option, 0.14).
//   d >= 2 diagonal step two places nearer the floor (drop 20 m over
//          14.14 m = 1.41, beats the cardinal 1.1).
//
// Upper triangle (c > r), and its mirror below:
//   (0,1)->(1,1)  (1,2)->(2,2)  (2,3)->(3,3)  (3,4)->(4,4)      cardinal
//   (0,2)->(1,1)  (1,3)->(2,2)  (2,4)->(3,3)                    diagonal
//   (0,3)->(1,2)->(2,2)   (1,4)->(2,3)->(3,3)                   diag + card
//   (0,4)->(1,3)->(2,2)                                         2 diagonals
// Contributing pixels on the floor: (0,0) 1, (1,1) 6, (2,2) 15, (3,3) 22,
// (4,4) 25.
//
// Threshold 6: walk (4,4) -> (3,3) [22] -> (2,2) [15] -> (1,1) [6]; (0,0) [1]
// stops it. Trunk (1,1)..(4,4), three diagonal steps: rivlen = 30 sqrt 2.
// Distance to the trunk, upper triangle:
//   d=1: 4 x 10                     = 40
//   d=2: 3 x 10 sqrt2               = 30 sqrt2
//   d=3: 2 x (10 sqrt2 + 10)        = 20 + 20 sqrt2
//   d=4: 20 sqrt2                   = 20 sqrt2
//   total 60 + 70 sqrt2; the lower triangle is the same; (0,0) adds 10 sqrt2.
//   loc_rivlen = (120 + 150 sqrt2) / 25 = 4.8 + 6 sqrt2
//
// Threshold 7: (1,1) [6] no longer qualifies, p_in = (2,2), rivlen = 20 sqrt2.
// The six pixels that ended at (1,1) -- (1,1), (0,0), (0,1), (1,0), (0,2),
// (2,0) -- each travel one more diagonal step:
//   loc_rivlen = (120 + 210 sqrt2) / 25 = 4.8 + 8.4 sqrt2
// =============================================================================
void test_diagonal_valley()
{
    std::cout << "\n--- D: diagonal valley ---\n";
    const lsdtt::lsdtt_raster_cpu dem = make_raster(5, 5, [](int r, int c) {
        return 10.0 * std::abs(r - c) - (r + c);
    });
    const auto flow = make_flow(dem);
    const lsdtt::lsdtt_flow_base& fb = *flow;
    auto F = [&](int r, int c) { return fb.flat_index(r, c); };

    check_receiver(fb, "D audit", 0, 0, 1, 1);
    check_receiver(fb, "D audit", 0, 1, 1, 1);
    check_receiver(fb, "D audit", 0, 2, 1, 1);
    check_receiver(fb, "D audit", 0, 3, 1, 2);
    check_receiver(fb, "D audit", 0, 4, 1, 3);
    check_receiver(fb, "D audit", 4, 0, 3, 1);
    check_receiver(fb, "D audit", 3, 4, 4, 4);
    check("D audit contributing pixels on the floor 1, 6, 15, 22, 25",
          fb.s_vector_size(F(0, 0)) == 1 && fb.s_vector_size(F(1, 1)) == 6 &&
          fb.s_vector_size(F(2, 2)) == 15 && fb.s_vector_size(F(3, 3)) == 22 &&
          fb.s_vector_size(F(4, 4)) == 25);

    const auto part = lsdtt::partition_incremental_subbasins(fb, {F(4, 4)});
    check("D one sub-basin of 25 pixels", part.n_pixels == std::vector<std::int64_t>{25});

    {
        const auto rl = lsdtt::compute_subbasin_river_lengths(part, fb, 6);
        check("D threshold 6: trunk head (1,1), 4 trunk pixels",
              rl.trunk_top_flat[0] == F(1, 1) && rl.n_trunk_pixels[0] == 4);
        check_near("D threshold 6: rivlen = 30 sqrt2", rl.rivlen_m[0], 30.0 * k_sqrt2, k_tight);
        check_near("D threshold 6: loc_rivlen = 4.8 + 6 sqrt2",
                   rl.loc_rivlen_m[0], 4.8 + 6.0 * k_sqrt2, k_fd);
    }
    {
        const auto rl = lsdtt::compute_subbasin_river_lengths(part, fb, 7);
        check("D threshold 7: trunk head (2,2), 3 trunk pixels",
              rl.trunk_top_flat[0] == F(2, 2) && rl.n_trunk_pixels[0] == 3);
        check_near("D threshold 7: rivlen = 20 sqrt2", rl.rivlen_m[0], 20.0 * k_sqrt2, k_tight);
        check_near("D threshold 7: loc_rivlen = 4.8 + 8.4 sqrt2",
                   rl.loc_rivlen_m[0], 4.8 + 8.4 * k_sqrt2, k_fd);
    }
}

} // namespace

int main()
{
    try
    {
        test_ramp();
        test_v_pair();
        test_nested_tree();
        test_diagonal_valley();
    }
    catch (const std::exception& e)
    {
        std::cerr << "test_subbasin_partition: unexpected exception: " << e.what() << "\n";
        return 2;
    }

    std::cout << "\n" << (g_n_checks - g_n_failed) << " / " << g_n_checks << " checks passed\n";
    return (g_n_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
