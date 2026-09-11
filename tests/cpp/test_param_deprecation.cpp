// test_param_deprecation.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// unit_param_parser_deprecation
//
// Phase 3b of docs/plan_docs/plan_fluvial_surfaces_cli_split.md.
//
// WHAT IS BEING PINNED, AND WHY IT NEEDED BUILDING. The typed parse functions
// iterate over a CLI's DEFAULTS MAPS and pull matching keys out of the file, so
// a key present in the file but absent from every defaults map is never looked
// at by anything. Before this facility, renaming or retiring a parameter left
// an old parameter file silently a no-op: `find_floodplain: true` would do
// nothing at all, with no error and no warning, and the user would have every
// reason to believe the analysis ran.
//
// The three arms correspond to the three outcomes that must be distinguishable:
//
//   removed    the parameter is gone and ignoring it would mislead. The exact
//              migration message is printed and the run STOPS.
//   alias      the parameter exists under a new name. The value is carried
//              across, a warning names the replacement, and the run COMPLETES
//              with the value applied - not merely accepted and dropped.
//   unknown    the key is in neither list and no defaults map declares it. It
//              is REPORTED rather than silently swallowed, but it is only a
//              warning: a false positive must never stop a run.
//
// A fourth arm checks the negative case - a clean parameter file must produce
// none of the three - because three arms that only ever assert "something was
// printed" would all pass against a parser that printed warnings constantly.
//
// The removed arm exits the process, so every arm runs in a CHILD process: the
// driver re-executes itself with --child and inspects the exit code and the
// captured output.
// -----------------------------------------------------------------------------

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_param_parser.hpp"

namespace
{

int g_n_checks = 0;
int g_n_failed = 0;

void check(bool condition, const std::string& what)
{
    ++g_n_checks;
    if (condition)
    {
        std::cout << "PASS  " << what << "\n";
    }
    else
    {
        ++g_n_failed;
        std::cout << "FAIL  " << what << "\n";
    }
}

// The exact message a removed key must produce. Kept as a constant so the test
// and the child agree on it character for character - the point of the removed
// group is that the user is told precisely where the parameter went.
const std::string k_removed_message =
    "find_floodplain has moved to lsdtt-fluvial-surfaces";

// -----------------------------------------------------------------------------
// child: build a parser with a known key policy and report what it parsed
// -----------------------------------------------------------------------------
int run_child(const std::string& path, const std::string& fname)
{
    lsdtt::lsdtt_param_parser pp(path, fname);

    pp.set_removed_keys({
        {"find_floodplain", k_removed_message},
        {"apply_plane_merges",
         "apply_plane_merges never had any effect and has been removed"} });

    pp.set_deprecated_aliases({
        {"terrace_max_distance_to_channel_metres",
         "valley_corridor_half_width_metres"},
        {"band_reach_window_metres", "ladder_reach_window_metres"},
        {"flow_dem_is_filled", "raster_is_filled"},
        {"fill_dem_for_flow_routing", "fill_before_routing"},
        {"source_pixels_threshold", "threshold_contributing_pixels"},
        {"channel_threshold_contributing_pixels",
         "threshold_contributing_pixels"},
        {"d8_channel_threshold_pixels", "threshold_contributing_pixels"},
        {"flow_accumulation_threshold_pixels",
         "threshold_contributing_pixels"} });

    std::map<std::string, float>       ff;
    std::map<std::string, int>         ii;
    std::map<std::string, bool>        bb;
    std::map<std::string, std::string> ss;

    ff["valley_corridor_half_width_metres"] = -1.0f;
    ff["ladder_reach_window_metres"]        = -1.0f;
    ii["ladder_level_count"]                = 0;
    ii["threshold_contributing_pixels"]     = 1000;
    bb["write_band_index_raster"]           = false;
    bb["raster_is_filled"]                  = false;
    bb["fill_before_routing"]               = true;
    ss["reference_field_algorithm_choice"]  = "hand";

    pp.parse_all_parameters(ff, ii, bb, ss);

    const auto f = pp.get_float_parameters();
    const auto i = pp.get_int_parameters();
    const auto b = pp.get_bool_parameters();

    // Printed so the parent can check the value actually ARRIVED under the new
    // name. A warning that names the replacement but drops the value would be
    // worse than useless: the run would proceed on a default the user never
    // chose while being told their setting had been carried across.
    std::cout << "CHILD valley_corridor_half_width_metres="
              << f.at("valley_corridor_half_width_metres") << "\n";
    std::cout << "CHILD ladder_reach_window_metres="
              << f.at("ladder_reach_window_metres") << "\n";
    std::cout << "CHILD ladder_level_count="
              << i.at("ladder_level_count") << "\n";
    std::cout << "CHILD threshold_contributing_pixels="
              << i.at("threshold_contributing_pixels") << "\n";
    std::cout << "CHILD raster_is_filled="
              << b.at("raster_is_filled") << "\n";
    std::cout << "CHILD fill_before_routing="
              << b.at("fill_before_routing") << "\n";
    std::cout << "CHILD_COMPLETED\n";
    return 0;
}

// -----------------------------------------------------------------------------
// parent helpers
// -----------------------------------------------------------------------------
struct child_run
{
    int         exit_code = 0;
    std::string output;      // stdout and stderr, merged
};

child_run run_case(const std::string& self,
                   const std::filesystem::path& dir,
                   const std::string& fname,
                   const std::string& contents)
{
    {
        std::ofstream ofs(dir / fname);
        ofs << contents;
    }
    const std::filesystem::path log = dir / (fname + ".log");

    // Quote the paths: a build directory can perfectly well contain a space.
    std::ostringstream cmd;
    cmd << "\"" << self << "\" --child \"" << dir.string()
        << "/\" \"" << fname << "\" > \"" << log.string() << "\" 2>&1";

    child_run out;
    const int rc = std::system(cmd.str().c_str());
#ifdef _WIN32
    out.exit_code = rc;
#else
    // Normalise: a child that exits 1 and a child killed by a signal are both
    // "did not complete", but they are not the same thing and the test should
    // be able to tell them apart if it ever matters.
    out.exit_code = (rc == -1) ? -1
                  : (WIFEXITED(rc) ? WEXITSTATUS(rc) : 128 + WTERMSIG(rc));
#endif

    std::ifstream ifs(log);
    std::ostringstream buf;
    buf << ifs.rdbuf();
    out.output = buf.str();
    return out;
}

bool contains(const std::string& haystack, const std::string& needle)
{
    return haystack.find(needle) != std::string::npos;
}

// The base parameter file every arm starts from: valid, minimal, no policy keys.
const std::string k_clean_params =
    "read_fname: nonexistent.tif\n"
    "write_fname: dep\n"
    "ladder_level_count: 4\n"
    "reference_field_algorithm_choice: hand\n";

int run_parent(const std::string& self, const std::filesystem::path& dir)
{
    std::cout << "lsdtt-test-param-deprecation\n"
              << "Phase 3b: deprecated, removed and unknown parameter keys\n\n";

    std::filesystem::create_directories(dir);

    // ---- arm 0: the negative control ---------------------------------------
    // Without this the other three arms could all pass against a parser that
    // warned about everything, always.
    {
        std::cout << "-- a clean parameter file --\n";
        const child_run r = run_case(self, dir, "clean.param", k_clean_params);
        check(r.exit_code == 0, "a clean file parses and the run completes");
        check(contains(r.output, "CHILD_COMPLETED"),
              "  the child reaches the end");
        check(!contains(r.output, "WARNING"),
              "  and produces NO warning of any kind");
        check(contains(r.output, "CHILD ladder_level_count=4"),
              "  the ordinary value is parsed");
    }

    // ---- arm 1: a removed key is fatal, with its exact message --------------
    {
        std::cout << "\n-- a removed key --\n";
        const child_run r = run_case(self, dir, "removed.param",
            k_clean_params + "find_floodplain: true\n");
        check(r.exit_code != 0,
              "a removed key stops the run with a non-zero exit");
        check(r.exit_code == 1,
              "  and exits 1 rather than crashing");
        check(contains(r.output, k_removed_message),
              "  the EXACT migration message is printed");
        check(!contains(r.output, "CHILD_COMPLETED"),
              "  and the run does not continue past it");
    }

    // Two removed keys should both be listed: a user migrating an old file
    // should get the whole list in one run, not one key per attempt.
    {
        std::cout << "\n-- two removed keys --\n";
        const child_run r = run_case(self, dir, "removed2.param",
            k_clean_params + "find_floodplain: true\n"
                             "apply_plane_merges: true\n");
        check(r.exit_code == 1, "two removed keys also stop the run");
        check(contains(r.output, k_removed_message) &&
              contains(r.output, "apply_plane_merges never had any effect"),
              "  BOTH messages are printed, not just the first");
    }

    // ---- arm 2: an alias warns, carries its value, and completes ------------
    {
        std::cout << "\n-- a deprecated alias --\n";
        const child_run r = run_case(self, dir, "alias.param",
            k_clean_params + "terrace_max_distance_to_channel_metres: 1234\n");
        check(r.exit_code == 0, "a deprecated key does NOT stop the run");
        check(contains(r.output, "CHILD_COMPLETED"),
              "  the run completes");
        check(contains(r.output, "terrace_max_distance_to_channel_metres") &&
              contains(r.output, "valley_corridor_half_width_metres") &&
              contains(r.output, "WARNING"),
              "  a warning names both the old key and its replacement");
        check(contains(r.output,
                       "CHILD valley_corridor_half_width_metres=1234"),
              "  and the VALUE arrives under the new name");
        check(!contains(r.output, "are not used by this tool"),
              "  the consumed alias is not then also reported as unknown");
    }

    // If both names are given the new one must win, and the collision must be
    // reported: silently preferring either would leave the user unsure which
    // number the run actually used.
    {
        std::cout << "\n-- both the old and the new name --\n";
        const child_run r = run_case(self, dir, "both.param",
            k_clean_params + "terrace_max_distance_to_channel_metres: 1234\n"
                             "valley_corridor_half_width_metres: 4321\n");
        check(r.exit_code == 0, "supplying both names does not stop the run");
        check(contains(r.output,
                       "CHILD valley_corridor_half_width_metres=4321"),
              "  the NEW name wins");
        check(contains(r.output, "was also supplied"),
              "  and the collision is reported rather than resolved silently");
    }

    // ---- arm 2b: every Phase 4 alias has the same pinned contract -----------
    const std::vector<std::string> threshold_aliases = {
        "source_pixels_threshold",
        "channel_threshold_contributing_pixels",
        "d8_channel_threshold_pixels",
        "flow_accumulation_threshold_pixels"};
    for (const std::string& old_key : threshold_aliases)
    {
        std::cout << "\n-- Phase 4 threshold alias: " << old_key << " --\n";
        const child_run r = run_case(
            self, dir, old_key + ".param",
            k_clean_params + old_key + ": 321\n");
        check(r.exit_code == 0 && contains(r.output, "CHILD_COMPLETED"),
              old_key + " completes through the alias");
        check(contains(r.output, old_key) &&
                  contains(r.output, "threshold_contributing_pixels") &&
                  contains(r.output, "WARNING"),
              "  warning names the old and canonical threshold keys");
        check(contains(r.output, "CHILD threshold_contributing_pixels=321"),
              "  old threshold value arrives under the canonical key");
        check(!contains(r.output, "are not used by this tool"),
              "  consumed threshold alias is not reported as unknown");
    }

    {
        std::cout << "\n-- Phase 4 mandatory-fill alias --\n";
        const child_run r = run_case(
            self, dir, "flow_dem_is_filled.param",
            k_clean_params + "flow_dem_is_filled: true\n");
        check(r.exit_code == 0 && contains(r.output, "CHILD_COMPLETED"),
              "flow_dem_is_filled completes through the alias");
        check(contains(r.output, "flow_dem_is_filled") &&
                  contains(r.output, "raster_is_filled") &&
                  contains(r.output, "WARNING"),
              "  mandatory-fill warning names both keys");
        check(contains(r.output, "CHILD raster_is_filled=1"),
              "  mandatory-fill value arrives under the canonical key");
    }

    {
        std::cout << "\n-- Phase 4 optional-fill alias --\n";
        const child_run r = run_case(
            self, dir, "fill_dem_for_flow_routing.param",
            k_clean_params + "fill_dem_for_flow_routing: false\n");
        check(r.exit_code == 0 && contains(r.output, "CHILD_COMPLETED"),
              "fill_dem_for_flow_routing completes through the alias");
        check(contains(r.output, "fill_dem_for_flow_routing") &&
                  contains(r.output, "fill_before_routing") &&
                  contains(r.output, "WARNING"),
              "  optional-fill warning names both keys");
        check(contains(r.output, "CHILD fill_before_routing=0"),
              "  optional-fill value arrives under the canonical key");
    }

    {
        std::cout << "\n-- Phase 4 old and canonical threshold --\n";
        const child_run r = run_case(
            self, dir, "phase4_both.param",
            k_clean_params +
                "source_pixels_threshold: 321\n"
                "threshold_contributing_pixels: 654\n");
        check(r.exit_code == 0,
              "old and canonical threshold keys parse together");
        check(contains(r.output, "CHILD threshold_contributing_pixels=654"),
              "  canonical threshold wins over its alias");
        check(contains(r.output, "was also supplied"),
              "  Phase 4 threshold collision is reported");
    }

    // ---- arm 3: an unknown key is reported, not swallowed -------------------
    {
        std::cout << "\n-- an unknown key --\n";
        const child_run r = run_case(self, dir, "unknown.param",
            k_clean_params + "ladder_levl_count: 4\n");
        check(r.exit_code == 0, "an unknown key does NOT stop the run");
        check(contains(r.output, "CHILD_COMPLETED"),
              "  the run completes");
        check(contains(r.output, "ladder_levl_count"),
              "  the unknown key is named in a warning");
        check(contains(r.output, "are not used by this tool"),
              "  and is reported as unused rather than silently dropped");
    }

    // Reserved file-I/O keys are consumed by the parser itself and appear in no
    // defaults map. Reporting them as unknown would make the warning useless.
    {
        std::cout << "\n-- reserved file I/O keys --\n";
        const child_run r = run_case(self, dir, "reserved.param",
            "read_path: .\n"
            "write_path: .\n"
            "read_fname: nonexistent.tif\n"
            "write_fname: dep\n"
            "boundary_conditions: oooo\n"
            "ladder_level_count: 2\n");
        check(r.exit_code == 0, "a file of reserved keys parses");
        check(!contains(r.output, "are not used by this tool"),
              "  and none of them is reported as unknown");
    }

    std::cout << "\n" << (g_n_checks - g_n_failed) << " / " << g_n_checks
              << " checks passed\n"
              << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";
    return g_n_failed == 0 ? 0 : 1;
}

}  // namespace

int main(int argc, char* argv[])
{
    try
    {
        if (argc == 4 && std::string(argv[1]) == "--child")
        {
            return run_child(argv[2], argv[3]);
        }

        const std::filesystem::path dir =
            (argc >= 2) ? std::filesystem::path(argv[1])
                        : std::filesystem::temp_directory_path() /
                          "lsdtt_param_deprecation";
        return run_parent(argv[0], dir);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
}
