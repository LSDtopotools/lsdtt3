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
// number, a sloping (flagged) polygon and two polygons on one level, the
// too-little-training fallback, and the reach windows.
// -----------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

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
    // 21 cells per polygon. Linear quantiles over 21 values: p10 = the 3rd
    // value, p90 = the 19th, so each core is 17 cells.
    //   0 fp   1.0 .. 2.0       core 1.1 .. 1.9
    //   1 t1   4.5 .. 5.5       core 4.6 .. 5.4
    //   2 t1   9.5 .. 10.5      (labelled t1 but a level higher: numbers ignored)
    //   3 t2   3.0 .. 20.0      spread 13.6 > 3 x gap: flagged
    //   4 t3   4.8 .. 5.8       median 5.3, within the gap of polygon 1: same level
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
    const lsdtt::lsdtt_raster_cpu fr = row_raster(field);
    const double gap = 1.0;
    lsdtt::compute_training_polygon_stats(set, fr, 3.0 * gap);
    check(set.polygons[3].flagged && !set.polygons[0].flagged &&
              !set.polygons[2].flagged,
          "the sloping polygon is flagged, the others are not");
    check(std::fabs(set.polygons[1].median - 5.0) < 1e-6 &&
              std::fabs(set.polygons[1].p10 - 4.6) < 1e-6,
          "polygon median and p10");

    std::vector<std::pair<std::int32_t, double>> cells;
    for (std::size_t i = 0; i < field.size(); ++i)
        cells.emplace_back(set.polygon_of_pixel[i], static_cast<double>(field[i]));
    const float f19 = static_cast<float>(1.9), f46 = static_cast<float>(4.6);
    const float f57 = static_cast<float>(4.8 + 1.0 * 18 / 20.0);
    const float f96 = static_cast<float>(9.6);
    {
        const lsdtt::training_window_fit fit =
            lsdtt::fit_training_window(cells, set, gap, 10);
        check(fit.has_bound &&
                  std::fabs(fit.bound.threshold - 0.5 * (f19 + f46)) < 1e-6 &&
                  near(fit.bound.youden_j, 1.0),
              "bound: between the floodplain core and every terrace core");
        check(fit.levels.size() == 3 && fit.levels[1].size() == 2 &&
                  fit.levels[2].size() == 1 && fit.levels[2][0] == 2,
              "levels: floodplain, {1, 4} on one level, then polygon 2");
        check(fit.breaks.size() == 2 &&
                  std::fabs(fit.breaks[0] - 0.5 * (f19 + f46)) < 1e-6 &&
                  std::fabs(fit.breaks[1] - 0.5 * (f57 + f96)) < 1e-6,
              "breaks: floodplain / level 1, then level 1 / level 2");
    }
    {
        const lsdtt::training_window_fit fit =
            lsdtt::fit_training_window(cells, set, gap, 100);
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

    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
