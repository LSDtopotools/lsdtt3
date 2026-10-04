// test_distributed_hydrology_forcing.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Unit tests for src/core/lsdtt_distributed_hydrology_forcing (S-HYPE column).
//
// What is under test
// ------------------
//   load_distributed_hydrology_depth_series
//     The loader lsdtt-distributed-hydrology used to carry as a static
//     function. Metres and millimetre columns, the mm * 0.001 conversion
//     (bit-exact), the metres column winning when both are present, and the
//     rejection of negative and NaN depths with the CLI's original wording.
//
//   load_distributed_hydrology_temperature_series
//     Temperature: negative values accepted unchanged, non-finite values
//     rejected, missing column and header-only file rejected.
//
//   the HYPE oracle forcing
//     tests/expected/hype/snow/lsdtt_input/forcing.csv (columns date,
//     precipitation_mm, temperature_c) loads precipitation through the depth
//     loader and temperature through the signed loader, and every value
//     matches an independent parse of the same file. The oracle is located
//     from __FILE__ and is only read, never written.
//
// Every other input CSV is written into <output_dir> at run time; no data
// files are committed for this test.
//
// Usage: lsdtt-test-distributed-hydrology-forcing <output_dir>
// Prints one PASS/FAIL line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "lsdtt_distributed_hydrology_forcing.hpp"

namespace
{

int g_n_checks = 0;
int g_n_failed = 0;

std::string g_out_dir;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << std::endl;
}

// Write `content` to <output_dir>/<name> and return the full path.
std::string write_csv(const std::string& name, const std::string& content)
{
    const std::string path = g_out_dir + "/" + name;
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write test input " + path);
    out << content;
    return path;
}

// Run `fn`, expecting a std::runtime_error whose message contains every one
// of `needles`. Records one check.
template <typename Fn>
void expect_throw(const std::string& name, Fn fn,
                  const std::vector<std::string>& needles)
{
    try
    {
        fn();
        check(name, false, "no exception thrown");
    }
    catch (const std::runtime_error& e)
    {
        const std::string msg = e.what();
        for (const std::string& n : needles)
        {
            if (msg.find(n) == std::string::npos)
            {
                check(name, false, "message lacks '" + n + "': " + msg);
                return;
            }
        }
        check(name, true);
    }
}

bool same_vector(const std::vector<double>& a, const std::vector<double>& b)
{
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (a[i] != b[i]) return false;
    return true;
}

// -----------------------------------------------------------------------------
// Depth series
// -----------------------------------------------------------------------------
void test_depth_series()
{
    // Metres column, read as-is.
    {
        const std::string p = write_csv("depth_metres.csv",
            "timestep,rainfall_metres\n1,0.0\n2,0.0125\n3,0.002\n");
        const std::vector<double> s =
            lsdtt::load_distributed_hydrology_depth_series(p, "rainfall_metres", "rainfall_mm");
        check("depth: metres column read unscaled",
              same_vector(s, {0.0, 0.0125, 0.002}));
    }

    // Millimetre column, multiplied by 0.001 -- pinned bit-for-bit.
    {
        const std::string p = write_csv("depth_mm.csv",
            "rainfall_mm\n0\n2.5\n12.75\n0.1\n");
        const std::vector<double> s =
            lsdtt::load_distributed_hydrology_depth_series(p, "rainfall_metres", "rainfall_mm");
        check("depth: mm column scaled by exactly 0.001",
              same_vector(s, {0.0 * 0.001, 2.5 * 0.001, 12.75 * 0.001, 0.1 * 0.001}));
    }

    // Both columns present: the metres column wins and the mm column is
    // ignored entirely -- even a value in it that would otherwise be rejected.
    {
        const std::string p = write_csv("depth_both.csv",
            "rainfall_mm,rainfall_metres\n-5,0.001\n7,0.003\n");
        const std::vector<double> s =
            lsdtt::load_distributed_hydrology_depth_series(p, "rainfall_metres", "rainfall_mm");
        check("depth: metres column preferred when both present",
              same_vector(s, {0.001, 0.003}));
    }

    // Negative depth.
    {
        const std::string p = write_csv("depth_negative.csv",
            "rainfall_metres\n0.001\n-0.002\n");
        expect_throw("depth: negative value rejected",
            [&] { lsdtt::load_distributed_hydrology_depth_series(p, "rainfall_metres", "rainfall_mm"); },
            {"row 2", "'rainfall_metres' is negative or non-finite",
             "A depth per timestep cannot be either."});
    }

    // NaN depth. std::stod parses "nan", so this reaches the sign policy.
    {
        const std::string p = write_csv("depth_nan.csv",
            "rainfall_mm\n1\nnan\n");
        expect_throw("depth: NaN value rejected",
            [&] { lsdtt::load_distributed_hydrology_depth_series(p, "rainfall_metres", "rainfall_mm"); },
            {"row 2", "'rainfall_mm' is negative or non-finite",
             "A depth per timestep cannot be either."});
    }

    // Infinite depth.
    {
        const std::string p = write_csv("depth_inf.csv",
            "rainfall_metres\ninf\n");
        expect_throw("depth: infinite value rejected",
            [&] { lsdtt::load_distributed_hydrology_depth_series(p, "rainfall_metres", "rainfall_mm"); },
            {"row 1", "negative or non-finite"});
    }

    // Neither column.
    {
        const std::string p = write_csv("depth_missing.csv",
            "timestep,rain\n1,0.1\n");
        expect_throw("depth: missing column rejected",
            [&] { lsdtt::load_distributed_hydrology_depth_series(p, "rainfall_metres", "rainfall_mm"); },
            {"has no column 'rainfall_metres' and no column 'rainfall_mm'",
             "One of the two is required"});
    }

    // Header only.
    {
        const std::string p = write_csv("depth_header_only.csv",
            "rainfall_metres\n");
        expect_throw("depth: header-only file rejected",
            [&] { lsdtt::load_distributed_hydrology_depth_series(p, "rainfall_metres", "rainfall_mm"); },
            {"has a header but no data rows"});
    }
}

// -----------------------------------------------------------------------------
// Signed series
// -----------------------------------------------------------------------------
void test_signed_series()
{
    // Negative temperatures are accepted and not converted.
    {
        const std::string p = write_csv("temp_signed.csv",
            "date,temperature_c\n2001-01-01,-9\n2001-01-02,0\n2001-01-03,3.25\n"
            "2001-01-04,-0.5\n");
        const std::vector<double> s = lsdtt::load_distributed_hydrology_temperature_series(p, "temperature_c");
        check("signed: negative values accepted, unconverted",
              same_vector(s, {-9.0, 0.0, 3.25, -0.5}));
    }

    // A column with a mm-like name is still read unconverted: the signed
    // loader has no unit conversion.
    {
        const std::string p = write_csv("temp_no_scale.csv",
            "temperature_c,temperature_mm\n12,5\n");
        const std::vector<double> s = lsdtt::load_distributed_hydrology_temperature_series(p, "temperature_c");
        check("signed: named column only, no scaling", same_vector(s, {12.0}));
    }

    // Non-finite.
    {
        const std::string p = write_csv("temp_nan.csv",
            "temperature_c\n-1\nnan\n");
        expect_throw("signed: NaN rejected",
            [&] { lsdtt::load_distributed_hydrology_temperature_series(p, "temperature_c"); },
            {"row 2", "'temperature_c' is non-finite",
             "A value per timestep must be finite."});
    }
    {
        const std::string p = write_csv("temp_inf.csv",
            "temperature_c\n-inf\n");
        expect_throw("signed: -inf rejected",
            [&] { lsdtt::load_distributed_hydrology_temperature_series(p, "temperature_c"); },
            {"row 1", "'temperature_c' is non-finite"});
    }

    // Missing column.
    {
        const std::string p = write_csv("temp_missing.csv",
            "date,temperature_k\n2001-01-01,270\n");
        expect_throw("signed: missing column rejected",
            [&] { lsdtt::load_distributed_hydrology_temperature_series(p, "temperature_c"); },
            {"has no column 'temperature_c'", "It is required"});
    }

    // Header only.
    {
        const std::string p = write_csv("temp_header_only.csv",
            "temperature_c\n");
        expect_throw("signed: header-only file rejected",
            [&] { lsdtt::load_distributed_hydrology_temperature_series(p, "temperature_c"); },
            {"has a header but no data rows"});
    }
}

// -----------------------------------------------------------------------------
// The HYPE oracle forcing file
// -----------------------------------------------------------------------------
void test_oracle_forcing()
{
    const std::filesystem::path here = std::filesystem::path(__FILE__).parent_path();
    const std::string path =
        (here / ".." / "expected" / "hype" / "snow" / "lsdtt_input" / "forcing.csv")
            .lexically_normal().string();

    if (!std::filesystem::exists(path))
    {
        check("oracle: forcing.csv found", false, path);
        return;
    }

    // Independent parse: header date,precipitation_mm,temperature_c.
    std::vector<double> want_precip_m;
    std::vector<double> want_temp;
    bool any_negative_temp = false;
    {
        std::ifstream in(path);
        std::string line;
        std::getline(in, line);
        if (!line.empty() && line.back() == '\r') line.pop_back();
        check("oracle: header is date,precipitation_mm,temperature_c",
              line == "date,precipitation_mm,temperature_c", line);
        while (std::getline(in, line))
        {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty()) continue;
            std::stringstream ss(line);
            std::string date, p, t;
            std::getline(ss, date, ',');
            std::getline(ss, p, ',');
            std::getline(ss, t, ',');
            want_precip_m.push_back(std::stod(p) * 0.001);
            want_temp.push_back(std::stod(t));
            if (want_temp.back() < 0.0) any_negative_temp = true;
        }
    }

    const std::vector<double> precip =
        lsdtt::load_distributed_hydrology_depth_series(path, "precipitation_m", "precipitation_mm");
    const std::vector<double> temp =
        lsdtt::load_distributed_hydrology_temperature_series(path, "temperature_c");

    check("oracle: precipitation_mm loads via depth loader as metres",
          !precip.empty() && same_vector(precip, want_precip_m),
          std::to_string(precip.size()) + " rows");
    check("oracle: temperature_c loads via signed loader",
          !temp.empty() && same_vector(temp, want_temp),
          std::to_string(temp.size()) + " rows");
    check("oracle: file exercises negative temperatures", any_negative_temp);
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "Usage: lsdtt-test-distributed-hydrology-forcing <output_dir>\n";
        return 2;
    }
    g_out_dir = argv[1];

    try
    {
        std::filesystem::create_directories(g_out_dir);
        test_depth_series();
        test_signed_series();
        test_oracle_forcing();
    }
    catch (const std::exception& e)
    {
        std::cout << "FAIL  unexpected exception: " << e.what() << std::endl;
        ++g_n_failed; ++g_n_checks;
    }

    std::cout << "--------------------------\n"
              << g_n_checks - g_n_failed << " / " << g_n_checks << " checks passed\n"
              << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << std::endl;
    return g_n_failed == 0 ? 0 : 1;
}
