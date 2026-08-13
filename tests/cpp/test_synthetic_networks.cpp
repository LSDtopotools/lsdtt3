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
#include <sstream>
#include <string>
#include <vector>

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

void check_case(const std::string& case_name,
                const std::string& name,
                bool               ok,
                const std::string& detail = "")
{
    check(case_name + " :: " + name, ok, detail);
}

struct synthetic_case
{
    std::string              name;
    lsdtt::lsdtt_ocn_options opts;
};

std::vector<std::uint8_t> build_island_mask(int n_rows, int n_cols)
{
    std::vector<std::uint8_t> mask(static_cast<std::size_t>(n_rows) *
                                       static_cast<std::size_t>(n_cols),
                                   static_cast<std::uint8_t>(0));
    for (int r = 1; r < n_rows - 1; ++r)
    {
        for (int c = 1; c < n_cols - 1; ++c)
        {
            mask[static_cast<std::size_t>(r) * static_cast<std::size_t>(n_cols) +
                 static_cast<std::size_t>(c)] = static_cast<std::uint8_t>(1);
        }
    }

    // Carve small coastal embayments for a simple island-like shape.
    mask[1 * n_cols + 1]                 = 0;
    mask[1 * n_cols + (n_cols - 2)]      = 0;
    mask[(n_rows - 2) * n_cols + 2]      = 0;
    mask[(n_rows - 2) * n_cols + n_cols - 3] = 0;
    return mask;
}

void run_case(const synthetic_case& sc)
{
    std::cout << "running case: " << sc.name << "\n";

    lsdtt::lsdtt_synthetic_network_cpu engine(sc.opts);

    const double e0 = engine.initial_energy();
    check_case(sc.name,
               "initial energy finite",
               std::isfinite(e0) && e0 > 0.0,
               "initial energy must be > 0");

    engine.run();

    const double e1 = engine.final_energy();
    check_case(sc.name,
               "final energy finite",
               std::isfinite(e1) && e1 > 0.0,
               "final energy must be > 0");
    if (sc.opts.temperature_scale == 0.0 && sc.opts.initial_temperature <= 0.0)
    {
        check_case(sc.name,
                   "energy does not increase",
                   e1 <= e0 + 1e-12 * std::abs(e0),
                   "anneal should not increase energy in this deterministic setup");
    }

    const auto dem = engine.get_elevation_raster();
    const auto area = engine.get_drainage_area_cells_raster();
    const std::size_t expected = static_cast<std::size_t>(sc.opts.n_rows) *
                                 static_cast<std::size_t>(sc.opts.n_cols);
    const int n_cols = sc.opts.n_cols;

    check_case(sc.name,
               "dem size matches grid",
               dem.data.size() == expected && dem.mask.size() == expected,
               "DEM vectors should be n_rows*n_cols");
    check_case(sc.name,
               "area size matches grid",
               area.data.size() == expected && area.mask.size() == expected,
               "area vectors should be n_rows*n_cols");

    const auto& areas = engine.areas();
    const auto& rec = engine.receivers();

    bool sea_receiver_minus_one = true;
    bool outlet_maps_to_self = true;
    bool non_outlet_land_not_to_sea = true;
    bool d8_receiver_ok = true;
    bool receiver_walk_terminates = true;
    bool area_sane = true;
    bool downstream_dem_strictly_decreases = true;
    bool outlet_dem_at_baselevel = true;

    std::size_t n_land = 0;
    double outlet_area_sum = 0.0;
    const double baselevel = static_cast<double>(sc.opts.baselevel_elevation);
    const double baselevel_tolerance =
        1e-9 * std::max(1.0, std::abs(baselevel));

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

            if (dem.mask[i] != 0)
            {
                const double z_outlet = static_cast<double>(dem.data[i]);
                if (!std::isfinite(z_outlet) ||
                    std::abs(z_outlet - baselevel) > baselevel_tolerance)
                {
                    outlet_dem_at_baselevel = false;
                }
            }
            continue;
        }

        const std::int64_t r = rec[i];
        if (r < 0 || r >= static_cast<std::int64_t>(expected) || !engine.is_land(r))
        {
            non_outlet_land_not_to_sea = false;
            d8_receiver_ok = false;
            downstream_dem_strictly_decreases = false;
            continue;
        }

        if (dem.mask[i] != 0 && dem.mask[static_cast<std::size_t>(r)] != 0)
        {
            const double z_i = static_cast<double>(dem.data[i]);
            const double z_r = static_cast<double>(dem.data[static_cast<std::size_t>(r)]);
            if (!std::isfinite(z_i) || !std::isfinite(z_r) || !(z_i > z_r))
                downstream_dem_strictly_decreases = false;
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

    check_case(sc.name,
               "sea receiver is -1",
               sea_receiver_minus_one,
               "all sea pixels should have receiver -1");
    check_case(sc.name,
               "outlet maps to self",
               outlet_maps_to_self,
               "all outlets should route to themselves");
    check_case(sc.name,
               "non-outlet land does not drain to sea",
               non_outlet_land_not_to_sea,
               "non-outlet land receivers must be land cells");
    check_case(sc.name,
               "D8 receiver neighbor constraint",
               d8_receiver_ok,
               "non-outlet land receivers must be one-cell D8 neighbors");
    check_case(sc.name,
               "receiver forest walk terminates",
               receiver_walk_terminates,
               "land receiver walks should terminate at an outlet self-loop");
    check_case(sc.name,
               "land area basic sanity",
               area_sane,
               "land drainage area should be finite, >= 1, and near-integer");
    check_case(sc.name,
               "outlet area sums to land count",
               outlet_area_partitions_land,
               "sum(area[outlets]) should equal total number of land cells");
    check_case(sc.name,
               "DEM decreases strictly downstream",
               downstream_dem_strictly_decreases,
               "for non-outlet land pixels, DEM must be > DEM[receiver]");
    check_case(sc.name,
               "outlet DEM equals baselevel",
               outlet_dem_at_baselevel,
               "outlet DEM should equal baselevel_elevation");
}

void run_determinism_checks(const std::string& case_name,
                            lsdtt::lsdtt_ocn_options opts)
{
    std::cout << "running case: " << case_name << "\n";

    // Same-seed repeatability check.
    lsdtt::lsdtt_synthetic_network_cpu run1(opts);
    lsdtt::lsdtt_synthetic_network_cpu run2(opts);
    run1.run();
    run2.run();

    const auto& rec1 = run1.receivers();
    const auto& rec2 = run2.receivers();
    check_case(case_name,
               "same-seed receivers exactly equal",
               rec1 == rec2,
               "receiver arrays should match exactly for identical options and seed");

    const auto& a1 = run1.areas();
    const auto& a2 = run2.areas();
    bool areas_same = (a1.size() == a2.size());
    constexpr double area_tol = 1e-12;
    if (areas_same)
    {
        for (std::size_t i = 0; i < a1.size(); ++i)
        {
            const double lhs = static_cast<double>(a1[i]);
            const double rhs = static_cast<double>(a2[i]);
            const double tol = area_tol * std::max(1.0, std::max(std::abs(lhs), std::abs(rhs)));
            if (!std::isfinite(lhs) || !std::isfinite(rhs) || std::abs(lhs - rhs) > tol)
            {
                areas_same = false;
                break;
            }
        }
    }
    check_case(case_name,
               "same-seed areas equal (tiny tolerance)",
               areas_same,
               "drainage area arrays should be identical up to tiny floating-point tolerance");

    const auto dem1 = run1.get_elevation_raster();
    const auto dem2 = run2.get_elevation_raster();
    bool dem_same_valid = dem1.mask.size() == dem2.mask.size() &&
                          dem1.data.size() == dem2.data.size() &&
                          dem1.mask == dem2.mask;
    if (dem_same_valid)
    {
        for (std::size_t i = 0; i < dem1.data.size(); ++i)
        {
            if (dem1.mask[i] == 0) continue;
            if (dem1.data[i] != dem2.data[i])
            {
                dem_same_valid = false;
                break;
            }
        }
    }
    check_case(case_name,
               "same-seed DEM data equal on valid cells",
               dem_same_valid,
               "DEM masks should match and valid-cell DEM values should be exactly equal");

    // Different-seed diversity check.
    opts.seed = 1;
    lsdtt::lsdtt_synthetic_network_cpu seed1(opts);
    seed1.run();

    opts.seed = 2;
    lsdtt::lsdtt_synthetic_network_cpu seed2(opts);
    seed2.run();

    const auto& r1 = seed1.receivers();
    const auto& r2 = seed2.receivers();
    const std::size_t n = std::min(r1.size(), r2.size());

    std::size_t land_pixels = 0;
    std::size_t different_receivers = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
        const auto idx = static_cast<std::int64_t>(i);
        if (!seed1.is_land(idx)) continue;
        ++land_pixels;
        if (r1[i] != r2[i]) ++different_receivers;
    }

    const double diff_fraction =
        (land_pixels > 0)
            ? static_cast<double>(different_receivers) / static_cast<double>(land_pixels)
            : 0.0;
    const bool diff_seed_receivers_sufficiently_different = diff_fraction >= 0.01;

    std::ostringstream diff_detail;
    diff_detail << "different receivers=" << different_receivers
                << "/" << land_pixels
                << " (fraction=" << diff_fraction << ")";

    check_case(case_name,
               "different seeds change >=1% of land receivers",
               diff_seed_receivers_sufficiently_different,
               diff_detail.str());
}

} // namespace

int main(int argc, char** argv)
{
    const std::string out_dir = (argc > 1) ? argv[1] : ".";
    std::cout << "synthetic-networks scaffold output dir: " << out_dir << "\n";

    lsdtt::lsdtt_ocn_options base_opts;
    base_opts.n_rows = 8;
    base_opts.n_cols = 10;
    base_opts.cell_size = 50.0;
    base_opts.n_iter = 200;
    base_opts.seed = 42;
    base_opts.outlet_mode = lsdtt::ocn_outlet_mode::single_south;

    std::vector<synthetic_case> cases;

    synthetic_case case_a;
    case_a.name = "A rectangle 8x10 gamma=0.5 Tscale=0";
    case_a.opts = base_opts;
    case_a.opts.gamma = 0.5;
    case_a.opts.temperature_scale = 0.0;
    cases.push_back(case_a);

    synthetic_case case_b;
    case_b.name = "B rectangle 8x10 gamma=0.75 Tscale=1";
    case_b.opts = base_opts;
    case_b.opts.gamma = 0.75;
    case_b.opts.temperature_scale = 1.0;
    cases.push_back(case_b);

    synthetic_case case_c;
    case_c.name = "C masked island-like gamma=0.5 Tscale=0";
    case_c.opts = base_opts;
    case_c.opts.gamma = 0.5;
    case_c.opts.temperature_scale = 0.0;
    case_c.opts.domain_mask = build_island_mask(base_opts.n_rows, base_opts.n_cols);
    cases.push_back(case_c);

    for (const auto& sc : cases)
    {
        run_case(sc);
    }

    lsdtt::lsdtt_ocn_options deterministic_opts = base_opts;
    deterministic_opts.gamma = 0.5;
    deterministic_opts.temperature_scale = 0.0;
    run_determinism_checks("D determinism rectangle 8x10 gamma=0.5 Tscale=0",
                           deterministic_opts);

    std::cout << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL")
              << "  (" << (g_n_checks - g_n_failed) << "/" << g_n_checks
              << " checks passed)\n";

    return (g_n_failed == 0) ? 0 : 1;
}
