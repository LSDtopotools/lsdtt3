// test_synthetic_networks.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>

#include "lsdtt_synthetic_network_cpu.hpp"

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

} // namespace

int main(int argc, char** argv)
{
    const std::string out_dir = (argc > 1) ? argv[1] : ".";
    std::cout << "synthetic-networks scaffold output dir: " << out_dir << "\n";

    lsdtt::lsdtt_ocn_options opts;
    opts.n_rows = 8;
    opts.n_cols = 10;
    opts.cell_size = 50.0;
    opts.gamma = 0.5;
    opts.n_iter = 200;
    opts.temperature_scale = 0.0;
    opts.seed = 42;
    opts.outlet_mode = lsdtt::ocn_outlet_mode::single_south;

    lsdtt::lsdtt_synthetic_network_cpu engine(opts);

    const double e0 = engine.initial_energy();
    check("initial energy finite", std::isfinite(e0) && e0 > 0.0,
          "initial energy must be > 0");

    engine.run();

    const double e1 = engine.final_energy();
    check("final energy finite", std::isfinite(e1) && e1 > 0.0,
          "final energy must be > 0");
    check("energy does not increase", e1 <= e0 + 1e-12 * std::abs(e0),
          "anneal should not increase energy in this deterministic setup");

    const auto dem = engine.get_elevation_raster();
    const auto area = engine.get_drainage_area_cells_raster();
    const std::size_t expected =
        static_cast<std::size_t>(opts.n_rows) * static_cast<std::size_t>(opts.n_cols);
    const int n_cols = opts.n_cols;

    check("dem size matches grid", dem.data.size() == expected && dem.mask.size() == expected,
          "DEM vectors should be n_rows*n_cols");
    check("area size matches grid", area.data.size() == expected && area.mask.size() == expected,
          "area vectors should be n_rows*n_cols");

    const auto& areas = engine.areas();
    const auto& rec = engine.receivers();

    bool sea_receiver_minus_one = true;
    bool outlet_maps_to_self = true;
    bool non_outlet_land_not_to_sea = true;
    bool d8_receiver_ok = true;
    bool receiver_walk_terminates = true;
    bool area_sane = true;

    std::size_t n_land = 0;
    double outlet_area_sum = 0.0;

    for (std::size_t i = 0; i < expected; ++i)
    {
        const auto flat = static_cast<std::int64_t>(i);
        if (!engine.is_land(flat))
        {
            if (rec[i] != -1) sea_receiver_minus_one = false;
            continue;
        }

        ++n_land;

        const double a = areas[i];
        const double nearest_int = std::nearbyint(a);
        if (!std::isfinite(a) || a < 1.0 || std::abs(a - nearest_int) > 1e-9)
            area_sane = false;

        if (engine.is_outlet(flat))
        {
            if (rec[i] != flat) outlet_maps_to_self = false;
            outlet_area_sum += a;
            continue;
        }

        const std::int64_t r = rec[i];
        if (r < 0 || r >= static_cast<std::int64_t>(expected) || !engine.is_land(r))
        {
            non_outlet_land_not_to_sea = false;
            d8_receiver_ok = false;
            continue;
        }

        const int row = static_cast<int>(flat) / n_cols;
        const int col = static_cast<int>(flat) % n_cols;
        const int r_row = static_cast<int>(r) / n_cols;
        const int r_col = static_cast<int>(r) % n_cols;
        const int dr = std::abs(r_row - row);
        const int dc = std::abs(r_col - col);
        if (dr > 1 || dc > 1 || (dr == 0 && dc == 0)) d8_receiver_ok = false;
    }

    for (std::size_t i = 0; i < expected; ++i)
    {
        const auto start = static_cast<std::int64_t>(i);
        if (!engine.is_land(start)) continue;

        std::int64_t cur = start;
        std::size_t  steps = 0;
        while (true)
        {
            const std::int64_t r = rec[static_cast<std::size_t>(cur)];
            if (r < 0 || r >= static_cast<std::int64_t>(expected) || !engine.is_land(r))
            {
                receiver_walk_terminates = false;
                break;
            }
            if (r == cur)
            {
                if (!engine.is_outlet(cur)) receiver_walk_terminates = false;
                break;
            }
            ++steps;
            if (steps > expected)
            {
                receiver_walk_terminates = false;
                break;
            }
            cur = r;
        }
        if (!receiver_walk_terminates) break;
    }

    const double area_sum_tolerance = 1e-9 * std::max<std::size_t>(n_land, 1);
    const bool outlet_area_partitions_land =
        std::abs(outlet_area_sum - static_cast<double>(n_land)) <= area_sum_tolerance;

    check("sea receiver is -1", sea_receiver_minus_one,
          "all sea pixels should have receiver -1");
    check("outlet maps to self", outlet_maps_to_self,
          "all outlets should route to themselves");
    check("non-outlet land does not drain to sea", non_outlet_land_not_to_sea,
          "non-outlet land receivers must be land cells");
    check("D8 receiver neighbor constraint", d8_receiver_ok,
          "non-outlet land receivers must be one-cell D8 neighbors");
    check("receiver forest walk terminates", receiver_walk_terminates,
          "land receiver walks should terminate at an outlet self-loop");
    check("land area basic sanity", area_sane,
          "land drainage area should be finite, >= 1, and near-integer");
    check("outlet area sums to land count", outlet_area_partitions_land,
          "sum(area[outlets]) should equal total number of land cells");

    std::cout << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL")
              << "  (" << (g_n_checks - g_n_failed) << "/" << g_n_checks
              << " checks passed)\n";

    return (g_n_failed == 0) ? 0 : 1;
}
