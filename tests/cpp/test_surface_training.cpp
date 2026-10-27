// test_surface_training.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Training polygons (lsdtt_surface_training_cpu), fluvial-surface methods plan
// 7i. Hand-checkable populations: the class parser, best separation on
// disjoint and overlapping samples, a window fit with a mislabelled terrace
// number, a sloping (flagged) polygon, a high terrace kept by the relative
// limit and two polygons on one level, the too-little-training fallback, the
// reach windows, and the patch level clustering test
// (lsdtt_surface_level_clustering_cpu).
// -----------------------------------------------------------------------------

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "lsdtt_surface_level_clustering_cpu.hpp"
#include "lsdtt_surface_training_cpu.hpp"

namespace
{
int g_failures = 0;
void check(bool ok, const std::string& what)
{
    std::cout << (ok ? "PASS  " : "FAIL  ") << what << "\n";
    if (!ok) ++g_failures;
}
bool near(double a, double b) { return std::fabs(a - b) < 1e-9; }

lsdtt::lsdtt_raster_cpu row_raster(std::vector<float> v, double cell = 1.0)
{
    lsdtt::lsdtt_raster_meta m;
    m.n_rows = 1;
    m.n_cols = static_cast<int>(v.size());
    m.cell_size = cell;
    m.no_data_value = -9999.0;
    std::vector<std::uint8_t> mask(v.size(), 1);
    return lsdtt::lsdtt_raster_cpu(std::move(v), std::move(mask), m);
}
} // namespace

int main()
{
    // ---- class labels ---------------------------------------------------------
    {
        int n = -1;
        check(lsdtt::parse_training_class("fp", n) == 0 && n == 0, "'fp' is floodplain");
        check(lsdtt::parse_training_class(" Floodplain ", n) == 0, "'Floodplain' is floodplain");
        check(lsdtt::parse_training_class("t1", n) == 1 && n == 1, "'t1' is terrace 1");
        check(lsdtt::parse_training_class("Terrace 2", n) == 1 && n == 2, "'Terrace 2' is terrace 2");
        check(lsdtt::parse_training_class("terrace", n) == 1 && n == 0, "'terrace' is terrace");
        check(lsdtt::parse_training_class("terrace0", n) == 1 && n == 0,
              "'terrace0' is terrace (numbers are ignored)");
        check(lsdtt::parse_training_class("T12", n) == 1 && n == 12, "'T12' is terrace 12");
        check(lsdtt::parse_training_class("fan", n) < 0 &&
              lsdtt::parse_training_class("tx", n) < 0 &&
              lsdtt::parse_training_class("t", n) < 0, "'fan', 'tx', 't' are refused");
    }

    // ---- best separation ----------------------------------------------------
    {
        const lsdtt::training_split s = lsdtt::best_separation({1, 2, 3}, {5, 6});
        check(near(s.threshold, 4.0) && near(s.youden_j, 1.0) &&
              s.lower_below == 3 && s.upper_above == 2,
              "disjoint: midpoint of the empty interval [3, 5), J 1");
    }
    {
        // J: t=3 -> 3/4 + 4/4 - 1 = 0.75, the maximum, constant on [3, 4).
        const lsdtt::training_split s =
            lsdtt::best_separation({1, 2, 3, 6}, {4, 5, 7, 8});
        check(near(s.threshold, 3.5) && near(s.youden_j, 0.75) &&
              s.lower_below == 3 && s.upper_above == 4,
              "overlapping: threshold 3.5, J 0.75");
    }

    // ---- a window fit ---------------------------------------------------------
    // 21 cells per polygon. Linear quantiles over 21 values: p25 = the 6th
    // value, p75 = the 16th, so each core (p25..p75) is 11 cells.
    //   0 fp   1.0 .. 2.0       core 1.25 .. 1.75
    //   1 t1   4.5 .. 5.5       core 4.75 .. 5.25
    //   2 t1   9.5 .. 10.5      (labelled t1 but a level higher: numbers ignored)
    //   3 t2   3.0 .. 20.0      IQR 8.5 > 3 x gap AND > 0.25 x median 11.5: flagged
    //   4 t3   4.8 .. 5.8       median 5.3, within the gap of polygon 1: same level
    //   5 t4   100 .. 110       IQR 5 > 3 x gap but < 0.25 x 105: one level
    lsdtt::training_set set;
    std::vector<float> field;
    auto add = [&](const std::string& label, double lo, double hi)
    {
        lsdtt::training_polygon p;
        p.label = label;
        p.cls = lsdtt::parse_training_class(label, p.number);
        const std::int32_t id = static_cast<std::int32_t>(set.polygons.size());
        set.polygons.push_back(p);
        for (int k = 0; k <= 20; ++k)
        {
            field.push_back(static_cast<float>(lo + (hi - lo) * k / 20.0));
            set.polygon_of_pixel.push_back(id);
        }
    };
    add("fp", 1.0, 2.0);
    add("t1", 4.5, 5.5);
    add("t1", 9.5, 10.5);
    add("t2", 3.0, 20.0);
    add("t3", 4.8, 5.8);
    add("t4", 100.0, 110.0);
    const lsdtt::lsdtt_raster_cpu fr = row_raster(field);
    const double gap = 1.0;
    lsdtt::compute_training_polygon_stats(set, fr, 3.0 * gap, 0.25);
    check(set.polygons[3].flagged && !set.polygons[0].flagged &&
              !set.polygons[2].flagged && !set.polygons[5].flagged,
          "the sloping polygon is flagged; the high terrace is kept by the "
          "relative limit");
    check(std::fabs(set.polygons[1].median - 5.0) < 1e-6 &&
              std::fabs(set.polygons[1].p25 - 4.75) < 1e-6 &&
              std::fabs(set.polygons[1].p75 - 5.25) < 1e-6,
          "polygon median, p25 and p75");
    {
        lsdtt::training_set abs_only = set;
        lsdtt::compute_training_polygon_stats(abs_only, fr, 3.0 * gap, 0.0);
        check(abs_only.polygons[5].flagged,
              "with the relative limit off the high terrace is flagged");
    }

    std::vector<std::pair<std::int32_t, double>> cells;
    for (std::size_t i = 0; i < field.size(); ++i)
        cells.emplace_back(set.polygon_of_pixel[i], static_cast<double>(field[i]));
    const float f19 = static_cast<float>(1.75), f46 = static_cast<float>(4.75);
    const float f57 = static_cast<float>(4.8 + 1.0 * 15 / 20.0);
    const float f96 = static_cast<float>(9.75);
    const float f1025 = static_cast<float>(10.25), f1025b = static_cast<float>(102.5);
    {
        const lsdtt::training_window_fit fit =
            lsdtt::fit_training_window(cells, set, gap, 10);
        check(fit.has_bound &&
                  std::fabs(fit.bound.threshold - 0.5 * (f19 + f46)) < 1e-6 &&
                  near(fit.bound.youden_j, 1.0),
              "bound: between the floodplain core and every terrace core");
        check(fit.levels.size() == 4 && fit.levels[1].size() == 2 &&
                  fit.levels[2].size() == 1 && fit.levels[2][0] == 2 &&
                  fit.levels[3][0] == 5,
              "levels: floodplain, {1, 4} on one level, polygon 2, polygon 5");
        check(fit.breaks.size() == 3 &&
                  std::fabs(fit.breaks[0] - 0.5 * (f19 + f46)) < 1e-6 &&
                  std::fabs(fit.breaks[1] - 0.5 * (f57 + f96)) < 1e-6 &&
                  std::fabs(fit.breaks[2] - 0.5 * (f1025 + f1025b)) < 1e-4,
              "breaks: floodplain / 1, 1 / 2, 2 / 3 at the core midpoints");
    }
    {
        const lsdtt::training_window_fit fit =
            lsdtt::fit_training_window(cells, set, gap, 10);
        check(fit.has_floodplain_level && fit.polygon_core_cells.size() == 6 &&
                  fit.polygon_core_cells[0].second == 11,
              "every polygon present is counted with its 11 core cells");
        check(fit.terrace_below_floodplain.empty(),
              "no terrace below the floodplain here");
    }
    {
        // A "t9" drawn on floodplain ground: 0.5 .. 1.5 m, median 1.0 m,
        // below the floodplain polygon's median 1.5 m.
        lsdtt::training_set s2 = set;
        lsdtt::training_polygon low;
        low.label = "t9";
        low.cls = lsdtt::parse_training_class(low.label, low.number);
        low.n_with_field = 21;
        low.median = 1.0; low.p25 = 0.75; low.p75 = 1.25;
        s2.polygons.push_back(low);
        std::vector<std::pair<std::int32_t, double>> c2 = cells;
        for (int k = 0; k <= 20; ++k) c2.emplace_back(6, 0.5 + k / 20.0);
        const lsdtt::training_window_fit fit =
            lsdtt::fit_training_window(c2, s2, gap, 10);
        check(fit.terrace_below_floodplain.size() == 1 &&
                  fit.terrace_below_floodplain[0] == 6 &&
                  fit.note.find("terrace polygon 6") != std::string::npos,
              "a terrace polygon below the floodplain polygon is reported");
    }
    {
        const lsdtt::training_window_fit fit =
            lsdtt::fit_training_window(cells, set, gap, 1000);
        check(!fit.has_bound && fit.breaks.empty() && !fit.note.empty(),
              "too little training: no bound, no breaks, a note");
    }
    {
        // Floodplain cells only.
        std::vector<std::pair<std::int32_t, double>> fp_only(cells.begin(),
                                                             cells.begin() + 21);
        const lsdtt::training_window_fit fit =
            lsdtt::fit_training_window(fp_only, set, gap, 10);
        check(!fit.has_bound && fit.breaks.empty() && fit.levels.size() == 1,
              "floodplain alone: one level, no bound");
    }

    // ---- reach windows --------------------------------------------------------
    {
        std::vector<float> d(30);
        for (int i = 0; i < 30; ++i) d[static_cast<std::size_t>(i)] = static_cast<float>(i);
        const lsdtt::lsdtt_raster_cpu coord = row_raster(d);
        std::vector<std::int32_t> src(30, 0);
        src[29] = -1;
        const std::vector<std::uint8_t> defined(30, 1);
        // 10 m windows stepped at 5 m over 0..28: 28 / 5 -> 6 windows.
        const lsdtt::reach_window_set rw =
            lsdtt::assign_reach_windows(src, coord, defined, 10.0);
        check(rw.windows.size() == 6, "6 windows over 28 m at a 5 m step");
        check(rw.window_of_pixel[0] == 0 && rw.window_of_pixel[12] == 2 &&
                  rw.window_of_pixel[28] == 5 && rw.window_of_pixel[29] == -1,
              "pixels fall in (d - dmin) / step; no source -> -1");
        check(rw.windows[0].neighbours.size() == 2 &&
                  rw.windows[3].neighbours.size() == 3 &&
                  rw.windows[3].neighbours[0] == 2,
              "pooled neighbours k-1..k+1, clipped at the ends");
        check(near(rw.windows[1].centre_m, 7.5), "window centre dmin + (k + 0.5) step");
    }

    // ---- patch level clustering ----------------------------------------------
    {
        // Few patches have little power: eight values on three levels over
        // 19.1 m leave an empty fraction of 0.963, and random heights 0.90.
        const std::vector<double> lv = {1.0, 1.2, 1.1, 10.0, 10.3, 20.0, 20.1, 19.9};
        const lsdtt::level_clustering_result r = lsdtt::test_level_clustering(lv, 1.5);
        check(r.clusters.size() == 3 && near(r.clusters[0].centre, 1.1) &&
                  near(r.clusters[1].centre, 10.15) && r.clusters[2].n == 3,
              "clusters split at gaps >= the required gap; centres are medians");
        check(near(r.empty_fraction, (8.8 + 9.6) / 19.1),
              "empty fraction = large gaps / range");
        check(r.p_value > 0.05, "eight patches: not distinguishable from random");
        // Thirty patches, ten on each of three levels: p at its floor.
        std::vector<double> many;
        for (const double b : {1.0, 10.0, 20.0})
            for (int i = 0; i < 10; ++i) many.push_back(b + 0.03 * i);
        const lsdtt::level_clustering_result rm =
            lsdtt::test_level_clustering(many, 1.5);
        check(rm.clusters.size() == 3 && rm.p_value <= 0.002 &&
                  rm.null_mean_empty_fraction < 0.5,
              "thirty patches on three levels cluster (p <= 0.002)");
        check(lsdtt::test_level_clustering(many, 1.5).p_value == rm.p_value &&
                  lsdtt::test_level_clustering(lv, 1.5).p_value == r.p_value,
              "same seed, same p");
    }
    {
        std::vector<double> even;
        for (int i = 0; i < 20; ++i) even.push_back(i);
        const lsdtt::level_clustering_result r = lsdtt::test_level_clustering(even, 1.5);
        check(r.clusters.size() == 1 && r.empty_fraction == 0.0 && r.p_value == 1.0,
              "evenly spread heights: one cluster, empty fraction 0, p 1");
    }
    {
        bool threw = false;
        try { (void)lsdtt::test_level_clustering({1.0, 2.0, 3.0}, 0.0); }
        catch (const std::invalid_argument&) { threw = true; }
        check(threw, "a required gap of 0 is refused");
    }

    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
