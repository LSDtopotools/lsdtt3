// test_output_file.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Characterisation test for open_output_file (src/utils/lsdtt_output_file.hpp,
// consolidation C14).
//
// 22 writers had the same three lines -- open an ofstream, test it, throw
// std::runtime_error(prefix + path) -- and now call open_output_file(path,
// prefix). The removed pattern is kept below and compared with the helper:
//
//   - an unwritable path: same exception type and exact what(), for each of
//     the 22 prefixes and for a path built as an expression (dir + "/name");
//   - a writable path: the returned stream writes, and opening truncates.
//
// Usage: lsdtt-test-output-file [scratch_dir]
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/utils/lsdtt_output_file.hpp"

using namespace lsdtt;

namespace
{

// The removed pattern, e.g. src/cpu/lsdtt_chi_network_analysis_cpu.cpp:1286-1289
// (as of 2026-10-08).
void old_open(const std::string& path, const std::string& prefix)
{
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error(prefix + path);
}

int g_failures = 0;

void check(const std::string& what, bool ok, const std::string& detail = "")
{
    std::cout << (ok ? "PASS " : "FAIL ") << what;
    if (!ok && !detail.empty()) std::cout << " -- " << detail;
    std::cout << "\n";
    if (!ok) ++g_failures;
}

template <class F>
std::string outcome(F f)
{
    try { f(); return "none"; }
    catch (const std::runtime_error& e) { return std::string("runtime_error|") + e.what(); }
    catch (const std::exception& e) { return std::string("exception|") + e.what(); }
}

} // namespace

int main(int argc, char** argv)
{
    const std::string dir = (argc > 1) ? argv[1] : "output_file_test";
    std::filesystem::create_directories(dir);
    std::cout << "== C14 open_output_file vs the removed open/test/throw lines ==\n";

    const std::vector<std::string> prefixes = {
        "lsdtt_data_loader: could not open CSV for writing: ",
        "lsdtt_param_parser::replace_and_print_parameter_file: cannot open ",
        "write_burned_stream_profiles_csv: cannot open output file: ",
        "write_junction_repair_log: cannot open ",
        "write_chi_disorder_subset_samples_csv: cannot open ",
        "write_chi_disorder_jackknife_csv: cannot open ",
        "write_chi_disorder_bootstrap_csv: cannot open ",
        "write_chi_disorder_sweep_csv: cannot open ",
        "write_chi_ksn_inversion_trials_csv: cannot open ",
        "write_chi_ksn_inversion_best_summary_csv: cannot open ",
        "write_chi_network_fit_nodes_csv: cannot open ",
        "write_chi_network_fit_segments_csv: cannot open ",
        "write_chi_ksn_basin_summary_csv: cannot open ",
        "Cannot open CRONUS v3 input file: ",
        "cannot write particle CSV: ",
        "lsdtt_hilltop_flow_routing_cpu: cannot open ",
        "cannot write ",
        "write_manning_widths_csv: could not open ",
        "write_channel_profile_csv: could not open ",
        "write_sti_lookup_csv: could not open CSV for writing: ",
        "write_macrostrat_lookup_csv: could not open CSV for writing: ",
    };
    const std::string missing_dir = dir + "/no_such_directory";
    int bad = 0;
    for (const std::string& p : prefixes)
    {
        const std::string path = missing_dir + "/out.csv";
        if (outcome([&] { old_open(path, p); }) != outcome([&] { (void)open_output_file(path, p); }))
            ++bad;
    }
    // A path built as an expression, as lsdtt_hype_raster_setup_cpu does.
    const std::string d = missing_dir;
    if (outcome([&] { old_open(d + "/GeoData.txt", "cannot write "); }) !=
        outcome([&] { (void)open_output_file(d + "/GeoData.txt", "cannot write "); }))
        ++bad;
    check("unwritable path: same exception type and message for all " +
              std::to_string(prefixes.size() + 1) + " prefixes",
          bad == 0);
    check("message text spelled out",
          outcome([&] { (void)open_output_file(missing_dir + "/x.csv", "who: cannot open "); }) ==
              "runtime_error|who: cannot open " + missing_dir + "/x.csv");

    const std::string ok_path = dir + "/written.csv";
    {
        std::ofstream pre(ok_path);
        pre << "old content that must be truncated\n";
    }
    {
        std::ofstream out = open_output_file(ok_path, "unused: ");
        out << "a,b\n1,2\n";
    }
    std::ifstream in(ok_path);
    std::stringstream ss;
    ss << in.rdbuf();
    check("writable path: stream writes and the file is truncated", ss.str() == "a,b\n1,2\n");

    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
