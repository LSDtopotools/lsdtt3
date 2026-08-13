// test_synthetic_networks.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)

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

    check("dem size matches grid", dem.data.size() == expected && dem.mask.size() == expected,
          "DEM vectors should be n_rows*n_cols");
    check("area size matches grid", area.data.size() == expected && area.mask.size() == expected,
          "area vectors should be n_rows*n_cols");

    bool all_finite_land_energy = true;
    const auto& areas = engine.areas();
    const auto& rec = engine.receivers();
    for (std::size_t i = 0; i < expected; ++i)
    {
        const auto flat = static_cast<std::int64_t>(i);
        if (!engine.is_land(flat)) continue;

        if (!std::isfinite(areas[i]) || areas[i] < 1.0)
        {
            all_finite_land_energy = false;
            break;
        }

        if (!engine.is_outlet(flat))
        {
            const std::int64_t r = rec[i];
            if (r < 0 || r >= static_cast<std::int64_t>(expected) || !engine.is_land(r))
            {
                all_finite_land_energy = false;
                break;
            }
        }
    }
    check("land invariants smoke", all_finite_land_energy,
          "land cells need finite area and valid land receivers");

    std::cout << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL")
              << "  (" << (g_n_checks - g_n_failed) << "/" << g_n_checks
              << " checks passed)\n";

    return (g_n_failed == 0) ? 0 : 1;
}
