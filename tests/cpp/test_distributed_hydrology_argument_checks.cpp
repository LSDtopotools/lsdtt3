// test_distributed_hydrology_argument_checks.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Characterisation test for src/cpu/lsdtt_distributed_hydrology_argument_checks_cpu.hpp.
//
// require_positive_finite replaced eight character-for-character copies of
// require_positive, and require_same_grid_dims three copies of
// require_same_grid, in the distributed-hydrology kernels. No other test pins
// the text of these messages, so this one does: for every input the shared
// function must throw exactly when the removed copy threw, with the same
// exception type and the same what() text.
//
// One copy of each is kept below as the reference; the other copies were the
// same text (only brace placement differed in lsdtt_hydro_velocity_field_cpu),
// which was checked by eye when they were removed and is recorded beside each
// reference.
//
// Usage: lsdtt-test-distributed-hydrology-argument-checks
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../../src/core/lsdtt_raster_meta.hpp"
#include "../../src/cpu/lsdtt_distributed_hydrology_argument_checks_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"

using namespace lsdtt;

namespace
{

// ---- removed copies (references) --------------------------------------------

// src/cpu/lsdtt_channel_extent_cpu.cpp:28-34 (as of 2026-10-07). Identical text
// in lsdtt_canopy_interception_cpu.cpp:23, lsdtt_hydro_velocity_field_cpu.cpp:27,
// lsdtt_hydro_routing_cpu.cpp:25, lsdtt_infiltration_excess_cpu.cpp:23,
// lsdtt_hydro_roughness_cpu.cpp:28, lsdtt_topmodel_dynamic_cpu.cpp:25,
// lsdtt_topmodel_cpu.cpp:27.
void old_require_positive(double v, const char* name, const char* who)
{
    if (!(std::isfinite(v) && v > 0.0))
        throw std::invalid_argument(
            std::string(who) + ": " + name + " must be finite and > 0 (got " +
            std::to_string(v) + ")");
}

// src/cpu/lsdtt_channel_extent_cpu.cpp:36-45 (as of 2026-10-07). Identical text
// in lsdtt_hydro_roughness_cpu.cpp:36 and lsdtt_hydro_velocity_field_cpu.cpp:37.
void old_require_same_grid(const lsdtt_raster_cpu& r, const lsdtt_raster_meta& meta,
                           const char* name, const char* who)
{
    const lsdtt_raster_meta& m = r.get_meta();
    if (m.n_rows != meta.n_rows || m.n_cols != meta.n_cols)
        throw std::invalid_argument(
            std::string(who) + ": " + name + " is " + std::to_string(m.n_rows) +
            "x" + std::to_string(m.n_cols) + " but the grid is " +
            std::to_string(meta.n_rows) + "x" + std::to_string(meta.n_cols));
}

// ---- harness -----------------------------------------------------------------

int g_failures = 0;

void check(const std::string& what, bool ok, const std::string& detail = "")
{
    std::cout << (ok ? "PASS " : "FAIL ") << what;
    if (!ok && !detail.empty()) std::cout << " -- " << detail;
    std::cout << "\n";
    if (!ok) ++g_failures;
}

// Runs f and returns {threw invalid_argument, what()}. Any other exception
// type is reported as a distinct outcome so a type change is caught.
template <class F>
std::pair<std::string, std::string> outcome(F f)
{
    try { f(); return {"none", ""}; }
    catch (const std::invalid_argument& e) { return {"invalid_argument", e.what()}; }
    catch (const std::exception& e) { return {"other", e.what()}; }
}

lsdtt_raster_cpu make_raster(int n_rows, int n_cols)
{
    lsdtt_raster_meta meta;
    meta.n_rows = n_rows;
    meta.n_cols = n_cols;
    const std::size_t n = static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_cols);
    return lsdtt_raster_cpu(std::vector<float>(n, 1.0f), std::vector<uint8_t>(n, 1), meta);
}

} // namespace

int main()
{
    std::cout << "== hydrology argument checks vs the removed private copies ==\n";

    const double inf = std::numeric_limits<double>::infinity();
    const std::vector<double> values = {
        1.0, 0.5, 1e300, std::numeric_limits<double>::denorm_min(),
        std::numeric_limits<double>::min(), 0.0, -0.0, -1.0, -1e-300,
        inf, -inf, std::numeric_limits<double>::quiet_NaN(),
        123456.789, 1e-7,
    };
    const char* names[] = {"manning_n", "", "a name with spaces"};
    const char* whos[]  = {"hydro_roughness", "topmodel_dynamic", ""};

    int n_pos = 0, n_pos_bad = 0;
    std::string first_bad;
    for (double v : values)
        for (const char* name : names)
            for (const char* who : whos)
            {
                ++n_pos;
                const auto a = outcome([&] { require_positive_finite(v, name, who); });
                const auto b = outcome([&] { old_require_positive(v, name, who); });
                if (a != b)
                {
                    ++n_pos_bad;
                    if (first_bad.empty())
                        first_bad = "v=" + std::to_string(v) + " new{" + a.first + "," +
                                    a.second + "} old{" + b.first + "," + b.second + "}";
                }
            }
    check("require_positive_finite == require_positive on " + std::to_string(n_pos) +
              " cases (type and message)", n_pos_bad == 0, first_bad);

    // The expected text, spelled out once, so a change to both copies at once
    // is still caught.
    check("require_positive_finite message text",
          outcome([] { require_positive_finite(-2.0, "k", "who"); }).second ==
              "who: k must be finite and > 0 (got -2.000000)");

    const std::pair<int, int> dims[] = {{4, 5}, {5, 4}, {4, 4}, {1, 1}, {0, 0}, {100, 3}};
    int n_grid = 0, n_grid_bad = 0;
    first_bad.clear();
    for (auto [rr, rc] : dims)
        for (auto [gr, gc] : dims)
        {
            ++n_grid;
            const lsdtt_raster_cpu r = make_raster(rr, rc);
            lsdtt_raster_meta grid;
            grid.n_rows = gr;
            grid.n_cols = gc;
            const auto a = outcome([&] { require_same_grid_dims(r, grid, "slope", "hydro_routing"); });
            const auto b = outcome([&] { old_require_same_grid(r, grid, "slope", "hydro_routing"); });
            if (a != b)
            {
                ++n_grid_bad;
                if (first_bad.empty())
                    first_bad = "new{" + a.first + "," + a.second + "} old{" +
                                b.first + "," + b.second + "}";
            }
        }
    check("require_same_grid_dims == require_same_grid on " + std::to_string(n_grid) +
              " dimension pairs (type and message)", n_grid_bad == 0, first_bad);

    {
        const lsdtt_raster_cpu r = make_raster(4, 5);
        lsdtt_raster_meta grid;
        grid.n_rows = 5;
        grid.n_cols = 4;
        check("require_same_grid_dims message text",
              outcome([&] { require_same_grid_dims(r, grid, "slope", "who"); }).second ==
                  "who: slope is 4x5 but the grid is 5x4");
    }

    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
