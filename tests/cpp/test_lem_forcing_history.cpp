// test_lem_forcing_history.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Forcing-history CSV: the optional process-parameter columns.
//
// erodibility and hillslope_diffusivity let a run walk the channelization index
// through time. Their trap is the default: uplift_rate and baselevel_fall_rate
// read an absent value as 0, but a 0 diffusivity or erodibility would switch a
// process off. So an absent column or empty cell must mean "keep the parameter
// file's value", which the loader stores as a negative number.
//
//   F1  Given values are read, per row, in any column order.
//   F2  An empty cell is "not given" (negative), not 0.
//   F3  A file without the columns gives "not given" everywhere, and the
//       original uplift / base-level columns are unaffected.
//   F4  at() is piecewise constant, so a stepped D ladder holds each value for
//       its whole interval.
//
// Usage: lsdtt-test-lem-forcing-history <scratch_dir>
// -----------------------------------------------------------------------------

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "../../src/cpu/lsdtt_lem_cpu.hpp"

namespace
{
int g_checks = 0, g_failures = 0;

void check(bool ok, const std::string& what)
{
    ++g_checks;
    if (!ok) { ++g_failures; std::cout << "  FAIL: " << what << "\n"; }
    else     { std::cout << "  ok  : " << what << "\n"; }
}

bool near(float a, float b) { return std::fabs(a - b) <= 1e-6f * std::max(1.0f, std::fabs(b)); }

std::string write_file(const std::string& dir, const std::string& name,
                       const std::string& body)
{
    const std::string path = dir + "/" + name;
    std::ofstream(path) << body;
    return path;
}
} // namespace

int main(int argc, char* argv[])
{
    using namespace lsdtt;
    const std::string dir = (argc > 1) ? argv[1] : ".";
    std::cout << "\n=== LEM forcing history: erodibility / hillslope_diffusivity ===\n";

    std::cout << "\n-- F1/F2: values, column order, empty cells --\n";
    {
        const auto h = load_lem_forcing_history(write_file(dir, "forcing_kd.csv",
            "hillslope_diffusivity,time_yr,erodibility,uplift_rate\n"
            "0.02,0,6.75e-5,0.0002\n"
            ",250000,1e-4,0.0002\n"
            "0.01,500000,,0.0002\n"));
        check(h.points.size() == 3, "three rows read");
        check(near(h.points[0].hillslope_diffusivity, 0.02f) &&
              near(h.points[0].erodibility, 6.75e-5f),
              "F1: D and K read from a reordered header");
        check(h.points[1].hillslope_diffusivity < 0.0f &&
              near(h.points[1].erodibility, 1e-4f),
              "F2: empty D cell is 'not given', K still read");
        check(near(h.points[2].hillslope_diffusivity, 0.01f) &&
              h.points[2].erodibility < 0.0f,
              "F2: empty K cell is 'not given', D still read");
    }

    std::cout << "\n-- F3: legacy file without the columns --\n";
    {
        const auto h = load_lem_forcing_history(write_file(dir, "forcing_legacy.csv",
            "time_yr,uplift_rate,baselevel_fall_rate\n0,0.0001,0.0\n1000,0.0001,0.002\n"));
        check(h.points.size() == 2 &&
              h.points[0].erodibility < 0.0f && h.points[0].hillslope_diffusivity < 0.0f &&
              h.points[1].erodibility < 0.0f && h.points[1].hillslope_diffusivity < 0.0f,
              "F3: absent columns are 'not given' in every row");
        check(near(h.points[1].uplift_rate, 1e-4f) &&
              near(h.points[1].baselevel_fall_rate, 0.002f),
              "F3: uplift and base-level columns unaffected");
    }

    std::cout << "\n-- F4: piecewise-constant lookup --\n";
    {
        const auto h = load_lem_forcing_history(write_file(dir, "forcing_ladder.csv",
            "time_yr,hillslope_diffusivity\n0,0.02\n250000,0.015\n500000,0.01\n"));
        check(near(h.at(0.0f).hillslope_diffusivity, 0.02f) &&
              near(h.at(249999.0f).hillslope_diffusivity, 0.02f),
              "F4: first value held to just before the step");
        check(near(h.at(250000.0f).hillslope_diffusivity, 0.015f) &&
              near(h.at(499000.0f).hillslope_diffusivity, 0.015f),
              "F4: second value from the step onwards");
        check(near(h.at(3.0e6f).hillslope_diffusivity, 0.01f),
              "F4: last value held past the final row");
    }

    std::cout << "\n" << (g_checks - g_failures) << " / " << g_checks
              << " checks passed\n";
    return g_failures == 0 ? 0 : 1;
}
