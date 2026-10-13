// test_csv_table_field_parsers.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// The strict field converters in lsdtt_csv_table (csv_parse_double_field,
// csv_parse_int_field, csv_parse_uint_field), moved from
// lsdtt_manning_hydraulics_cpu, and the two hand-written readers that now use
// them:
//
//   - lsdtt_crn_column::read_erosion_time_series (time,rate,removal), and
//   - load_hillslope_particles_csv (the particle restart file).
//
// For each reader: a well-formed file reads the same values as before (the old
// std::stod / std::stoull results), and a malformed one is reported with the
// file, the line and the offending text instead of a bare "stod" or a silent
// partial parse ("1.5e-4x" used to be read as 1.5e-4).
//
// Usage: lsdtt-test-csv-table-field-parsers <writable directory>
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_csv_table.hpp"
#include "../../src/cpu/lsdtt_crn_column.hpp"
#include "../../src/cpu/lsdtt_hillslope_particle_io.hpp"

using namespace lsdtt;

namespace
{

int g_failures = 0;

void check(const std::string& what, bool ok, const std::string& detail = "")
{
    std::cout << (ok ? "PASS " : "FAIL ") << what;
    if (!ok && !detail.empty()) std::cout << " -- " << detail;
    std::cout << "\n";
    if (!ok) ++g_failures;
}

// Runs f, which must throw std::runtime_error whose message contains every
// string in `needles`.
void check_throws(const std::string& what, const std::function<void()>& f,
                  const std::vector<std::string>& needles)
{
    try
    {
        f();
        check(what, false, "no exception");
    }
    catch (const std::runtime_error& e)
    {
        const std::string msg = e.what();
        bool ok = true;
        for (const std::string& n : needles) ok = ok && msg.find(n) != std::string::npos;
        check(what, ok, "message: " + msg);
    }
    catch (const std::exception& e)
    {
        check(what, false, std::string("wrong exception type: ") + e.what());
    }
}

void write_file(const std::string& path, const std::string& text)
{
    std::ofstream out(path, std::ios::binary);
    out << text;
}

const std::string k_particle_header =
    "id,cell,weight,source_id,age_yr,distance_cells,status,mixing_depth,seed,substep_index";

} // namespace

int main(int argc, char** argv)
{
    const std::filesystem::path dir = argc > 1 ? argv[1] : ".";
    std::filesystem::create_directories(dir);
    std::cout << "== strict CSV field parsers and the readers using them ==\n";

    // ---- the shared converters ----
    check("csv_parse_double_field reads '1e-4' as std::stod does",
          csv_parse_double_field("1e-4", "f.csv", 7) == std::stod("1e-4"));
    check("csv_parse_int_field reads '-12'", csv_parse_int_field("-12", "f.csv", 7) == -12);
    check("csv_parse_uint_field reads the largest uint64",
          csv_parse_uint_field("18446744073709551615", "f.csv", 7) ==
              std::numeric_limits<unsigned long long>::max());
    try
    {
        (void)csv_parse_double_field("1.5abc", "f.csv", 7);
        check("csv_parse_double_field rejects trailing characters", false);
    }
    catch (const std::runtime_error& e)
    {
        check("csv_parse_double_field message is the Manning one",
              std::string(e.what()) == "f.csv: line 7: cannot parse '1.5abc' as a number.", e.what());
    }
    check_throws("csv_parse_double_field rejects an empty field",
                 [] { (void)csv_parse_double_field("", "f.csv", 3); }, {"line 3", "''"});
    check_throws("csv_parse_int_field rejects '12.5'",
                 [] { (void)csv_parse_int_field("12.5", "f.csv", 3); }, {"'12.5'", "as an integer."});
    check_throws("csv_parse_uint_field rejects '-3' (stoull would wrap it)",
                 [] { (void)csv_parse_uint_field("-3", "f.csv", 3); },
                 {"'-3'", "as a non-negative integer."});

    // ---- lsdtt_crn_column::read_erosion_time_series ----
    {
        // CRLF line endings and spaces around fields, as a hand-edited file has.
        const std::string path = (dir / "crn_erosion_valid.csv").string();
        write_file(path, "time,rate,removal\r\n1000, 1e-4 ,0\r\n-9999,5e-5 , 0\r\n");
        lsdtt_crn_column col;
        col.read_erosion_time_series(path);
        const std::vector<double> t = {std::stod("1000"), std::stod("-9999")};
        const std::vector<double> r = {std::stod(" 1e-4 "), std::stod("5e-5 ")};
        const std::vector<double> m = {std::stod("0"), std::stod(" 0")};
        check("crn_column: a well-formed file (CRLF, padded fields) reads the std::stod values",
              col.time_vec() == t && col.rate_vec() == r && col.removal_vec() == m);
    }
    const auto crn_case = [&](const std::string& name, const std::string& bad_row,
                              const std::vector<std::string>& needles) {
        const std::string path = (dir / ("crn_erosion_" + name + ".csv")).string();
        write_file(path, "time,rate,removal\n1000,1e-4,0\n" + bad_row + "\n-9999,5e-5,0\n");
        std::vector<std::string> n = needles;
        n.push_back(path);
        check_throws("crn_column: " + name + " on line 3 is reported",
                     [&] { lsdtt_crn_column c; c.read_erosion_time_series(path); }, n);
    };
    crn_case("trailing junk", "2000,1e-4x,0", {"line 3", "'1e-4x'"});
    crn_case("quoted field", "\"2000\",1e-4,0", {"line 3", "'\"2000\"'"});
    crn_case("empty field", "2000,,0", {"line 3", "''"});
    crn_case("two columns", "2000,1e-4", {"line 3", "three comma-separated columns"});

    // ---- load_hillslope_particles_csv ----
    {
        const std::string path = (dir / "particles_valid.csv").string();
        write_file(path, k_particle_header + "\n4,12,1.25,3,10.5,2,0,0.75,123456789,42\n");
        const lsdtt_hillslope_particle_tracker p = load_hillslope_particles_csv(path);
        check("particle restart: a well-formed file reads the std::sto* values",
              p.id == std::vector<std::uint64_t>{std::stoull("4")} &&
                  p.cell == std::vector<std::int64_t>{std::stoll("12")} &&
                  p.weight == std::vector<double>{std::stod("1.25")} &&
                  p.source_id == std::vector<std::uint32_t>{3} &&
                  p.age == std::vector<double>{std::stod("10.5")} &&
                  p.distance_cells == std::vector<double>{std::stod("2")} &&
                  p.status == std::vector<std::uint8_t>{0} && p.mixing_depth == std::stod("0.75") &&
                  p.seed == std::stoull("123456789") && p.substep_index == 42);
    }
    {
        const std::string path = (dir / "particles_crlf.csv").string();
        write_file(path, k_particle_header + "\r\n4,12,1.25,3,10.5,2,0,0.75,123456789,42\r\n");
        check_throws("particle restart: a CRLF file is reported as CRLF",
                     [&] { (void)load_hillslope_particles_csv(path); }, {path, "CRLF"});
    }
    {
        const std::string path = (dir / "particles_bad_field.csv").string();
        write_file(path, k_particle_header + "\n4,12,1.25,3,10.5,2,0,0.75,123456789,42\n"
                                             "5,12x,1.25,3,10.5,2,0,0.75,123456789,42\n");
        check_throws("particle restart: a bad field on line 3 is reported",
                     [&] { (void)load_hillslope_particles_csv(path); }, {path, "line 3", "'12x'"});
    }
    {
        const std::string path = (dir / "particles_short_row.csv").string();
        write_file(path, k_particle_header + "\n4,12,1.25,3,10.5,2,0,0.75,123456789\n");
        check_throws("particle restart: a short row is reported with its line and field count",
                     [&] { (void)load_hillslope_particles_csv(path); },
                     {path, "line 2 has 9 fields, expected 10"});
    }

    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
