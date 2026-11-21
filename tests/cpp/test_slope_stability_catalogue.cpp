// test_slope_stability_catalogue.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Unit test for the landslide catalogues (src/cpu/lsdtt_slope_stability_catalogue_cpu).
//
// Hand-built rupture results pin each grouping rule exactly:
//   * the 1-D rule: a single flow line to D is dropped at 2 and kept at 0; two
//     converging paths are kept; a STABLE neighbour still counts; D on a grid
//     corner and edge is handled;
//   * "add outlets": a daylight point joins its own landslide when its own F
//     <= 1, and no longer counts towards its own downslope daylight point;
//   * amalgamation joins diagonally touching landslides; n_pieces counts the
//     8-connected pieces of one daylight point's cells;
//   * per-landslide statistics on a 3 x 7 rectangle against closed forms
//     (area, depth, volume, the regionprops axes 4 sqrt(n^2/12) b, perimeter,
//     centroid, the eq.-7 landslide safety factor).
// Then on a rough synthetic DEM through compute_rupture_planes, properties every
// catalogue must satisfy: every landslide cell has F <= 1 and either drains
// to its daylight point or is it; every kept daylight point has at least the
// required labelled neighbours (recomputed here); areas and cell counts agree
// between the two catalogues; volume = area x mean depth.
//
// Usage: lsdtt-test-slope-stability-catalogue
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_flow_base.hpp"
#include "../../src/core/lsdtt_raster_meta.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_slope_stability_catalogue_cpu.hpp"
#include "../../src/cpu/lsdtt_slope_stability_rupture_cpu.hpp"
#include "../../src/utils/lsdtt_array2d.hpp"

using namespace lsdtt;

namespace
{

int g_failures = 0;
int g_checks = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_checks;
    if (!ok) ++g_failures;
    std::cout << (ok ? "  PASS  " : "  FAIL  ") << name;
    if (!detail.empty()) std::cout << "  (" << detail << ")";
    std::cout << "\n";
}

std::string num(double v)
{
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.8g", v);
    return buf;
}

lsdtt_raster_meta make_meta(int n_rows, int n_cols, double cs)
{
    lsdtt_raster_meta meta;
    meta.n_rows = n_rows;
    meta.n_cols = n_cols;
    meta.cell_size = cs;
    meta.x_min = 500000.0;
    meta.y_min = 4000000.0;
    meta.no_data_value = -9999.0;
    return meta;
}

// A hand-built rupture result: set() makes a cell a source.
struct fake
{
    lsdtt_raster_meta meta;
    std::size_t n;
    std::vector<float> F, H, S;
    std::vector<std::uint8_t> m;
    std::vector<std::int64_t> day;
    fake(int nr, int nc, double cs) : meta(make_meta(nr, nc, cs)), n(static_cast<std::size_t>(nr) * nc),
        F(n, -9999.0f), H(n, -9999.0f), S(n, -9999.0f), m(n, 0), day(n, -1) {}
    std::size_t at(int r, int c) const { return static_cast<std::size_t>(r) * meta.n_cols + c; }
    void set(int r, int c, int dr, int dc, float f, float h = 1.0f, float s = 0.5f)
    {
        const std::size_t i = at(r, c);
        F[i] = f; H[i] = h; S[i] = s; m[i] = 1;
        day[i] = static_cast<std::int64_t>(at(dr, dc));
    }
    rupture_rasters build() const
    {
        rupture_rasters r;
        auto mk = [&](const std::vector<float>& v) { return lsdtt_raster_cpu(std::vector<float>(v), std::vector<std::uint8_t>(m), meta); };
        r.safety_factor = mk(F);
        r.depth = mk(H);
        r.rupture_slope = mk(S);
        r.daylight = day;
        return r;
    }
};

const std::map<std::int32_t, rupture_parameters> k_params{{0, rupture_parameters{}}};

landslide_catalogues run(const fake& f, int k)
{
    return build_landslide_catalogues(f.build(), f.meta, {}, k_params, k);
}

} // namespace

int main()
{
    std::cout << "lsdtt-test-slope-stability-catalogue\n";

    // ---------------------------------------------------------- 1-D rule
    {
        fake f(10, 10, 10.0);
        f.set(2, 5, 5, 5, 0.8f);
        f.set(3, 5, 5, 5, 0.8f);
        f.set(4, 5, 5, 5, 0.8f);        // only (4,5) is a neighbour of D
        const landslide_catalogues k2 = run(f, 2), k0 = run(f, 0);
        check("1-D rule: a single flow line to D is dropped at 2 and kept at 0",
              k2.by_daylight.empty() && k2.n_daylight_points_removed == 1
              && k0.by_daylight.size() == 1 && k0.by_daylight[0].n_cells == 3);
    }
    {
        fake f(10, 10, 10.0);
        f.set(4, 4, 5, 5, 0.8f);
        f.set(4, 6, 5, 5, 0.8f);
        const landslide_catalogues k2 = run(f, 2);
        check("1-D rule: two converging paths are kept",
              k2.by_daylight.size() == 1 && k2.by_daylight[0].n_cells == 2 && k2.n_daylight_points_removed == 0);
    }
    {
        fake f(10, 10, 10.0);
        f.set(4, 4, 5, 5, 0.9f);
        f.set(4, 6, 5, 5, 1.5f);         // stable, but still counts
        const landslide_catalogues k2 = run(f, 2);
        check("1-D rule: a stable neighbour still counts, and is not in the landslide",
              k2.by_daylight.size() == 1 && k2.by_daylight[0].n_cells == 1);
    }
    {
        fake f(6, 6, 10.0);
        f.set(0, 1, 0, 0, 0.5f);
        f.set(1, 0, 0, 0, 0.5f);         // D on the corner
        f.set(1, 2, 0, 3, 0.5f);
        f.set(1, 4, 0, 3, 0.5f);         // D on the north edge
        const landslide_catalogues k2 = run(f, 2);
        check("1-D rule: daylight points on a grid corner and edge are handled", k2.by_daylight.size() == 2);
    }

    // ------------------------------------------------------- add outlets
    {
        fake f(10, 10, 10.0);
        f.set(4, 4, 5, 5, 0.8f);
        f.set(4, 5, 5, 5, 0.8f);
        f.set(4, 6, 5, 5, 0.8f);
        f.set(5, 5, 6, 5, 0.7f);         // D is itself a source, draining to E
        f.set(6, 4, 6, 5, 0.7f);         // E's only other donor
        const landslide_catalogues k2 = run(f, 2);
        // D keeps (3 neighbours) and includes itself (its own F 0.7); E has only
        // (6,4) labelled E, because (5,5) was relabelled to itself.
        check("add outlets: D joins its own landslide; it no longer counts towards E, so E is dropped",
              k2.by_daylight.size() == 1 && k2.by_daylight[0].n_cells == 4
              && k2.by_daylight[0].daylight_flat == static_cast<std::int64_t>(f.at(5, 5))
              && k2.n_daylight_points_removed == 1);
    }

    // -------------------------------------------- amalgamation and pieces
    {
        fake f(12, 12, 10.0);
        f.set(2, 2, 3, 3, 0.6f);
        f.set(2, 3, 3, 3, 0.6f);         // landslide A -> (3,3)
        f.set(3, 4, 4, 5, 0.6f);         // touches A diagonally
        f.set(3, 5, 4, 5, 0.6f);         // landslide B -> (4,5)
        f.set(8, 1, 9, 4, 0.6f);         // far away, same D as the next two
        f.set(8, 3, 9, 4, 0.6f);
        f.set(8, 4, 9, 4, 0.6f);
        f.set(8, 8, 9, 4, 0.6f);
        const landslide_catalogues k2 = run(f, 2);
        int pieces = 0;
        for (const auto& r : k2.by_daylight)
            if (r.daylight_flat == static_cast<std::int64_t>(f.at(9, 4))) pieces = r.n_pieces;
        int merged = 0;
        for (const auto& r : k2.amalgamated) merged = std::max(merged, r.n_daylight_points);
        check("amalgamation: diagonally touching landslides merge (8-connected)",
              k2.by_daylight.size() == 3 && merged == 2);
        check("pieces: one daylight point's cells in three 8-connected pieces", pieces == 3);
        check("labels: by-daylight ids increase with the daylight point's flat index",
              k2.by_daylight[0].daylight_flat < k2.by_daylight[1].daylight_flat
              && k2.by_daylight[1].daylight_flat < k2.by_daylight[2].daylight_flat);
    }

    // ------------------------------------------------ statistics on a rectangle
    {
        const double cs = 10.0;
        fake f(10, 10, cs);
        const float h = 2.0f, s = 0.4f, fs = 0.5f;
        for (int r = 2; r <= 4; ++r)
            for (int c = 1; c <= 7; ++c) f.set(r, c, 5, 4, fs, h, s);
        const landslide_catalogues k2 = run(f, 2);
        const landslide_record& L = k2.by_daylight.at(0);
        const rupture_parameters p;
        const double t = std::atan(static_cast<double>(s));
        const double rho_g = p.density_kg_m3 * p.gravity_m_s2;
        const double fl = (p.cohesion_pa + p.friction_coefficient * rho_g * h * std::cos(t) * std::cos(t))
                        / (rho_g * h * std::cos(t) * std::sin(t));
        const double length = 4.0 * std::sqrt(49.0 / 12.0) * cs, width = 4.0 * std::sqrt(9.0 / 12.0) * cs;
        const lsdtt_raster_meta& M = f.meta;
        const double cx = M.x_min + (4.0 + 0.5) * cs, cy = M.y_min + (10.0 - 3.0 - 0.5) * cs;
        const bool ok = L.n_cells == 21 && std::fabs(L.area_m2 - 2100.0) < 1e-9
                     && std::fabs(L.mean_depth_m - 2.0) < 1e-12 && std::fabs(L.volume_m3 - 4200.0) < 1e-9
                     && std::fabs(L.length_m - length) < 1e-9 && std::fabs(L.width_m - width) < 1e-9
                     && std::fabs(L.perimeter_m - 200.0) < 1e-9
                     && std::fabs(L.centroid_x - cx) < 1e-9 && std::fabs(L.centroid_y - cy) < 1e-9
                     && std::fabs(L.landslide_safety_factor - fl) < 1e-9 * fl
                     && std::fabs(L.mean_safety_factor - 0.5) < 1e-12;
        check("statistics: a 3 x 7 rectangle (area, depth, volume, axes, perimeter, centroid, F_L)", ok,
              "length " + num(L.length_m) + " vs " + num(length) + ", F_L " + num(L.landslide_safety_factor) + " vs " + num(fl));
        check("statistics: the amalgamated record of a single landslide equals its by-daylight record",
              k2.amalgamated.size() == 1 && k2.amalgamated[0].n_cells == 21
              && std::fabs(k2.amalgamated[0].length_m - L.length_m) < 1e-12 && k2.amalgamated[0].n_daylight_points == 1);
    }

    // ----------------------------------------------- through the real pipeline
    {
        const int nr = 60, nc = 60;
        const double cs = 10.0;
        lsdtt_raster_meta meta = make_meta(nr, nc, cs);
        lsdtt_array2d_f d(nr, nc, 0.0f);
        lsdtt_array2d_u8 m(nr, nc, 1);
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
            {
                const double x = c * cs, y = (nr - r) * cs;
                d(r, c) = static_cast<float>(0.9 * y + 25.0 * std::sin(x / 60.0) * std::cos(y / 45.0) + 8.0 * std::sin(x / 17.0 + y / 23.0));
            }
        lsdtt_flowinfo_options o;
        o.n_threads = 1;
        const std::unique_ptr<lsdtt_flow_base> flow = make_flow_base(d, m, meta, o);
        const lsdtt_raster_cpu dem(lsdtt_array2d_f(d), lsdtt_array2d_u8(m), meta);
        const lsdtt_raster_result np = flow->get_n_contributing_pixels_raster();
        std::vector<std::uint8_t> channel(np.data.size(), 0);
        for (std::size_t i = 0; i < channel.size(); ++i) channel[i] = (np.mask[i] && np.data[i] >= 60.0f) ? 1 : 0;
        const std::map<std::int32_t, rupture_parameters> params{{0, {500.0, 0.6, 2600.0, 9.81}}};
        const rupture_rasters rup = compute_rupture_planes(*flow, dem, channel, {}, params);
        const landslide_catalogues cat = build_landslide_catalogues(rup, meta, {}, params, 2);

        // Recompute the labels after "add outlets" independently.
        const std::size_t n = channel.size();
        std::vector<std::int64_t> lab(rup.daylight);
        for (std::size_t i = 0; i < n; ++i) if (rup.daylight[i] >= 0) lab[static_cast<std::size_t>(rup.daylight[i])] = rup.daylight[i];
        bool cells_ok = true, kept_ok = true;
        std::int64_t cells = 0;
        for (std::size_t i = 0; i < n; ++i)
        {
            const std::int32_t id = cat.label[i];
            if (!id) continue;
            ++cells;
            const landslide_record& r = cat.by_daylight.at(static_cast<std::size_t>(id - 1));
            const double f = rup.safety_factor.get_data()[i];
            cells_ok &= f <= 1.0 && lab[i] == r.daylight_flat;
        }
        for (const landslide_record& r : cat.by_daylight)
        {
            const int dr = static_cast<int>(r.daylight_flat / nc), dc = static_cast<int>(r.daylight_flat % nc);
            int count = 0;
            for (int a = -1; a <= 1; ++a)
                for (int b2 = -1; b2 <= 1; ++b2)
                {
                    if (!a && !b2) continue;
                    const int rr = dr + a, cc = dc + b2;
                    if (rr < 0 || cc < 0 || rr >= nr || cc >= nc) continue;
                    if (lab[static_cast<std::size_t>(rr) * nc + cc] == r.daylight_flat) ++count;
                }
            kept_ok &= count >= 2;
        }
        double area_d = 0.0, area_a = 0.0;
        bool vol_ok = true;
        for (const auto& r : cat.by_daylight) { area_d += r.area_m2; vol_ok &= std::fabs(r.volume_m3 - r.area_m2 * r.mean_depth_m) < 1e-9 * r.volume_m3 + 1e-12; }
        for (const auto& r : cat.amalgamated) area_a += r.area_m2;
        check("pipeline: every landslide cell has F <= 1 and belongs to its daylight point",
              cells_ok && cells == cat.n_landslide_cells && !cat.by_daylight.empty(),
              std::to_string(cat.by_daylight.size()) + " landslides, " + std::to_string(cells) + " cells");
        check("pipeline: every kept daylight point has at least 2 labelled neighbours", kept_ok);
        // For one cell, eq. 7 is the cell's own eq. 2: F_L = F (float inputs).
        double worst_single = 0.0;
        int n_single = 0;
        const landslide_catalogues cat0 = build_landslide_catalogues(rup, meta, {}, params, 0);
        for (const auto& r : cat0.by_daylight)
            if (r.n_cells == 1) { ++n_single; worst_single = std::max(worst_single, std::fabs(r.landslide_safety_factor / r.mean_safety_factor - 1.0)); }
        check("pipeline: a one-cell landslide's F_L (eq. 7) equals its cell's F (eq. 2)",
              n_single > 0 && worst_single < 1e-5, std::to_string(n_single) + " one-cell landslides, max relative " + num(worst_single));
        check("pipeline: both catalogues cover the same area; volume = area x mean depth",
              std::fabs(area_d - area_a) < 1e-6 && std::fabs(area_d - cells * cs * cs) < 1e-6 && vol_ok
              && cat.amalgamated.size() <= cat.by_daylight.size());
    }

    check("arguments: min_daylight_neighbours outside 0..8 throws", [] {
        fake f(4, 4, 10.0);
        try { (void)build_landslide_catalogues(f.build(), f.meta, {}, k_params, 9); } catch (const std::invalid_argument&) { return true; }
        return false; }());

    std::cout << g_checks << " checks, " << g_failures << " failed\n";
    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
