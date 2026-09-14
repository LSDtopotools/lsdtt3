// test_param_deprecation.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
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
#include <iomanip>
#include <iostream>
#include <limits>
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
         "threshold_contributing_pixels"},
        {"flow_big_channel_threshold_pixels",
         "big_channel_threshold_contributing_pixels"},
        {"incision_threshold_area", "incision_threshold_specific_area"},
        {"channel_mask_threshold_area", "channel_mask_threshold_specific_area"},
        {"min_contributing_pixels", "basin_minimum_contributing_pixels"} });

    pp.set_deprecated_numeric_aliases({
        {"smoothing_window_diameter_metres",
         {"smoothing_window_radius_metres", 0.5, 0.0, true}},
        {"precision_diameter_metres",
         {"precision_radius_metres", 0.5, 0.0, true}} });
    pp.set_deprecated_compatibility_aliases({
        {"old_points_format",
         {"legacy_points_format", "write_vector_format"}},
        {"old_lines_format",
         {"legacy_lines_format", "write_vector_format"}} });

    std::map<std::string, float>       ff;
    std::map<std::string, int>         ii;
    std::map<std::string, bool>        bb;
    std::map<std::string, std::string> ss;
    std::map<std::string, double>      dd;

    ff["valley_corridor_half_width_metres"] = -1.0f;
    ff["ladder_reach_window_metres"]        = -1.0f;
    ff["smoothing_window_radius_metres"]    = -1.0f;
    ff["incision_threshold_specific_area"] = 0.0f;
    ff["channel_mask_threshold_specific_area"] = 0.0f;
    ii["ladder_level_count"]                = 0;
    ii["threshold_contributing_pixels"]     = 1000;
    ii["big_channel_threshold_contributing_pixels"] = 5000;
    ii["basin_minimum_contributing_pixels"] = 1000;
    bb["write_band_index_raster"]           = false;
    bb["raster_is_filled"]                  = false;
    bb["fill_before_routing"]               = true;
    ss["reference_field_algorithm_choice"]  = "hand";
    ss["write_vector_format"]               = "flatgeobuf";
    ss["legacy_points_format"]              = "same_as_write_vector_format";
    ss["legacy_lines_format"]               = "same_as_write_vector_format";
    dd["precision_radius_metres"]            = -1.0;

    pp.parse_all_parameters(ff, ii, bb, ss, dd);

    lsdtt::lsdtt_help_registry empty_help;
    pp.print_help(empty_help, pp.get_write_fname() + "_parser_io",
                  "test", "test");
    pp.print_parameters("dep_parser_params.param");

    const auto f = pp.get_float_parameters();
    const auto i = pp.get_int_parameters();
    const auto b = pp.get_bool_parameters();
    const auto d = pp.get_double_parameters();
    const auto s = pp.get_string_parameters();

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
    std::cout << "CHILD big_channel_threshold_contributing_pixels="
              << i.at("big_channel_threshold_contributing_pixels") << "\n";
    std::cout << "CHILD basin_minimum_contributing_pixels="
              << i.at("basin_minimum_contributing_pixels") << "\n";
    std::cout << "CHILD incision_threshold_specific_area="
              << f.at("incision_threshold_specific_area") << "\n";
    std::cout << "CHILD channel_mask_threshold_specific_area="
              << f.at("channel_mask_threshold_specific_area") << "\n";
    std::cout << "CHILD raster_is_filled="
              << b.at("raster_is_filled") << "\n";
    std::cout << "CHILD fill_before_routing="
              << b.at("fill_before_routing") << "\n";
    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10)
              << "CHILD smoothing_window_radius_metres="
              << f.at("smoothing_window_radius_metres") << "\n"
              << "CHILD precision_radius_metres="
              << d.at("precision_radius_metres") << "\n";
    std::cout << "CHILD write_vector_format=" << s.at("write_vector_format")
              << "\nCHILD legacy_points_format="
              << s.at("legacy_points_format")
              << "\nCHILD legacy_lines_format="
              << s.at("legacy_lines_format") << "\n";
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

    // ---- arm 2c: Phase 5 plain threshold aliases ---------------------------
    const std::vector<std::pair<std::string, std::string>> phase5_aliases = {
        {"flow_big_channel_threshold_pixels",
         "big_channel_threshold_contributing_pixels"},
        {"incision_threshold_area", "incision_threshold_specific_area"},
        {"channel_mask_threshold_area", "channel_mask_threshold_specific_area"},
        {"min_contributing_pixels", "basin_minimum_contributing_pixels"}};
    for (const auto& [old_key, new_key] : phase5_aliases)
    {
        const child_run r = run_case(
            self, dir, "phase5_" + old_key + ".param",
            k_clean_params + old_key + ": 321\n");
        check(r.exit_code == 0 && contains(r.output, "CHILD " + new_key + "=321"),
              old_key + " carries its unchanged value to " + new_key);
        check(contains(r.output, old_key) && contains(r.output, new_key) &&
                  contains(r.output, "WARNING"),
              "  Phase 5 plain-alias warning names both keys");
    }

    // ---- arm 2d: Phase 5 transformed numeric aliases -----------------------
    {
        std::cout << "\n-- transformed diameter-to-radius alias --\n";
        const child_run r = run_case(
            self, dir, "numeric_alias.param",
            k_clean_params + "smoothing_window_diameter_metres: 20.25\n");
        check(r.exit_code == 0 && contains(r.output, "CHILD_COMPLETED"),
              "a transformed numeric alias completes");
        check(contains(r.output,
                       "CHILD smoothing_window_radius_metres=10.125"),
              "  divide-by-two value arrives under the radius key");
        check(contains(r.output, "smoothing_window_diameter_metres") &&
                  contains(r.output, "smoothing_window_radius_metres") &&
                  contains(r.output,
                           "new = old * 0.5 + 0; zero remains zero") &&
                  contains(r.output, "WARNING"),
              "  warning names both keys and the exact transformation");
        check(!contains(r.output, "are not used by this tool"),
              "  transformed old key is consumed rather than unknown");
    }

    {
        std::cout << "\n-- transformed zero sentinel --\n";
        const child_run r = run_case(
            self, dir, "numeric_zero.param",
            k_clean_params + "smoothing_window_diameter_metres: 0\n");
        check(r.exit_code == 0 &&
                  contains(r.output,
                           "CHILD smoothing_window_radius_metres=0"),
              "zero sentinel remains zero");
    }

    {
        std::cout << "\n-- transformed alias and canonical key --\n";
        const child_run r = run_case(
            self, dir, "numeric_both.param",
            k_clean_params +
                "smoothing_window_diameter_metres: 20\n"
                "smoothing_window_radius_metres: 7.25\n");
        check(r.exit_code == 0 &&
                  contains(r.output,
                           "CHILD smoothing_window_radius_metres=7.25"),
              "canonical numeric key wins over transformed alias");
        check(contains(r.output, "was also supplied") &&
                  contains(r.output, "was not applied"),
              "  canonical precedence and skipped transformation are reported");
    }

    for (const std::string bad_value : {"not-a-number", "inf", "nan"})
    {
        std::cout << "\n-- invalid transformed numeric value: "
                  << bad_value << " --\n";
        const child_run r = run_case(
            self, dir, "numeric_bad_" + bad_value + ".param",
            k_clean_params + "smoothing_window_diameter_metres: " +
                bad_value + "\n");
        check(r.exit_code == 1 && !contains(r.output, "CHILD_COMPLETED"),
              bad_value + " fails before analysis");
        check(contains(r.output,
                       "FATAL: deprecated numeric parameter "
                       "'smoothing_window_diameter_metres' must be a finite "
                       "number for transformation to "
                       "'smoothing_window_radius_metres'") &&
                  contains(r.output, "got '" + bad_value + "'."),
              "  failure clearly names the key, requirement, and bad value");
    }

    {
        std::cout << "\n-- injective compatibility aliases --\n";
        const child_run r = run_case(
            self, dir, "compat_formats.param",
            k_clean_params +
                "write_vector_format: gpkg\n"
                "old_points_format: csv\n"
                "old_lines_format: flatgeobuf\n");
        check(r.exit_code == 0 &&
                  contains(r.output, "CHILD write_vector_format=gpkg") &&
                  contains(r.output, "CHILD legacy_points_format=csv") &&
                  contains(r.output, "CHILD legacy_lines_format=flatgeobuf"),
              "simultaneous old selector values remain distinct from the base");
        check(contains(r.output, "use 'write_vector_format'") &&
                   contains(r.output, "preserved for this output in this run"),
               "  compatibility warning names the public migration and preservation");

        std::ifstream params(dir / "dep_parser_params.param");
        std::ostringstream params_text;
        params_text << params.rdbuf();
        std::ifstream csv(dir / "dep_parser_io.csv");
        std::ostringstream csv_text;
        csv_text << csv.rdbuf();
        std::ifstream html(dir / "dep_parser_io.html");
        std::ostringstream html_text;
        html_text << html.rdbuf();
        check(!contains(params_text.str(), "legacy_points_format") &&
                  !contains(params_text.str(), "legacy_lines_format"),
              "  compatibility read values are absent from print_parameters");
        check(!contains(csv_text.str(), "legacy_points_format") &&
                  !contains(csv_text.str(), "legacy_lines_format") &&
                  !contains(html_text.str(), "legacy_points_format") &&
                  !contains(html_text.str(), "legacy_lines_format"),
              "  compatibility read values are absent from CSV and HTML help");
    }

    {
        std::cout << "\n-- compatibility storage defaults stay internal --\n";
        const child_run r = run_case(
            self, dir, "compat_default.param",
            k_clean_params + "old_points_format: csv\n");
        check(r.exit_code == 0 &&
                  contains(r.output, "CHILD legacy_points_format=csv") &&
                  contains(r.output,
                           "CHILD legacy_lines_format=same_as_write_vector_format"),
              "compatibility typed maps retain read and default storage values");

        std::ifstream params(dir / "dep_parser_params.param");
        std::ostringstream params_text;
        params_text << params.rdbuf();
        std::ifstream csv(dir / "dep_parser_io.csv");
        std::ostringstream csv_text;
        csv_text << csv.rdbuf();
        std::ifstream html(dir / "dep_parser_io.html");
        std::ostringstream html_text;
        html_text << html.rdbuf();
        check(!contains(params_text.str(), "legacy_points_format") &&
                  !contains(params_text.str(), "legacy_lines_format"),
              "  compatibility read/default storage is absent from print_parameters");
        check(!contains(csv_text.str(), "legacy_points_format") &&
                  !contains(csv_text.str(), "legacy_lines_format") &&
                  !contains(html_text.str(), "legacy_points_format") &&
                  !contains(html_text.str(), "legacy_lines_format"),
              "  compatibility read/default storage is absent from CSV and HTML help");
    }

    {
        std::cout << "\n-- transformed and ordinary aliases coexist --\n";
        const child_run r = run_case(
            self, dir, "alias_coexist.param",
            k_clean_params +
                "smoothing_window_diameter_metres: 8\n"
                "terrace_max_distance_to_channel_metres: 1234\n"
                "precision_diameter_metres: 0.24691357802469134\n");
        check(r.exit_code == 0 &&
                  contains(r.output,
                           "CHILD smoothing_window_radius_metres=4") &&
                  contains(r.output,
                           "CHILD valley_corridor_half_width_metres=1234"),
              "transformed and ordinary alias maps both apply");
        check(contains(r.output,
                       "CHILD precision_radius_metres=0.12345678901234566"),
              "  transformed text preserves double parsing precision");
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

        std::ifstream csv(dir / "dep_parser_io.csv");
        std::ostringstream csv_text;
        csv_text << csv.rdbuf();
        const std::string resolved_dir = dir.string() + "/";
        check(contains(csv_text.str(), "\"read_path\",\"string\",\"" +
                                      resolved_dir + "\"") &&
                  contains(csv_text.str(),
                           "\"read_fname\",\"string\",\"nonexistent.tif\"") &&
                  contains(csv_text.str(), "\"write_path\",\"string\",\"" +
                                      resolved_dir + "\"") &&
                  contains(csv_text.str(), "\"write_fname\",\"string\",\"dep\""),
              "  per-run help includes resolved parser-owned I/O without CLI declarations");
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
