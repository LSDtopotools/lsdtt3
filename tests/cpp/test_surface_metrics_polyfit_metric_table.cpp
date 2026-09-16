// test_surface_metrics_polyfit_metric_table.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Unit tests for the canonical polyfit metric name/suffix table in
// lsdtt_raster_surface_metrics_cpu.hpp.
//
// The table used to be duplicated across three CLIs and one of those copies
// carried only 21 of the 27 enum members, so a metric name accepted by two
// CLIs was rejected by the third. The completeness assertions below are the
// guard against that recurring: if a new member is added to polyfit_metric
// without a table row, this test fails.
// -----------------------------------------------------------------------------

#include <cstdio>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "../../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"

namespace {

int g_failures = 0;

void check(bool condition, const std::string& what)
{
    if (!condition)
    {
        ++g_failures;
        std::printf("FAIL: %s\n", what.c_str());
    }
}

void check_string(const std::string& got,
                  const std::string& expected,
                  const std::string& what)
{
    if (got != expected)
    {
        ++g_failures;
        std::printf("FAIL: %s -- expected '%s', got '%s'\n",
                    what.c_str(), expected.c_str(), got.c_str());
    }
}

// Every member of polyfit_metric, written out explicitly. This list is
// independent of the table so that it can be used to test the table.
const std::vector<lsdtt::polyfit_metric>& every_metric()
{
    static const std::vector<lsdtt::polyfit_metric> metrics =
    {
        lsdtt::polyfit_metric::fitted_surface,
        lsdtt::polyfit_metric::dzdx,
        lsdtt::polyfit_metric::dzdy,
        lsdtt::polyfit_metric::d2zdx2,
        lsdtt::polyfit_metric::d2zdxdy,
        lsdtt::polyfit_metric::d2zdy2,
        lsdtt::polyfit_metric::slope,
        lsdtt::polyfit_metric::slope_degrees,
        lsdtt::polyfit_metric::aspect_radians,
        lsdtt::polyfit_metric::aspect_degrees,
        lsdtt::polyfit_metric::laplacian,
        lsdtt::polyfit_metric::hessian_determinant,
        lsdtt::polyfit_metric::mean_curvature,
        lsdtt::polyfit_metric::gaussian_curvature,
        lsdtt::polyfit_metric::maximum_curvature,
        lsdtt::polyfit_metric::minimum_curvature,
        lsdtt::polyfit_metric::curvedness,
        lsdtt::polyfit_metric::shape_index,
        lsdtt::polyfit_metric::profile_curvature,
        lsdtt::polyfit_metric::plan_curvature,
        lsdtt::polyfit_metric::tangential_curvature,
        lsdtt::polyfit_metric::difference_curvature,
        lsdtt::polyfit_metric::horizontal_excess_curvature,
        lsdtt::polyfit_metric::vertical_excess_curvature,
        lsdtt::polyfit_metric::accumulation_curvature,
        lsdtt::polyfit_metric::ring_curvature,
        lsdtt::polyfit_metric::rotor
    };

    return metrics;
}

// -----------------------------------------------------------------------------
// test_table_completeness
//
// The table must carry exactly one row per enum member.
// -----------------------------------------------------------------------------
void test_table_completeness()
{
    const auto& table = lsdtt::all_polyfit_metrics();

    check(every_metric().size() == lsdtt::polyfit_metric_count,
          "enum roster size matches polyfit_metric_count (27)");

    check(table.size() == lsdtt::polyfit_metric_count,
          "table has exactly polyfit_metric_count (27) entries; got " +
          std::to_string(table.size()));

    for (lsdtt::polyfit_metric metric : every_metric())
    {
        int count = 0;
        for (const auto& entry : table)
        {
            if (entry.metric == metric)
            {
                ++count;
            }
        }

        check(count == 1,
              std::string("metric '") + lsdtt::polyfit_metric_name(metric) +
              "' appears exactly once in the table; found " +
              std::to_string(count));
    }

    // Every entry must carry a non-empty name and a non-empty suffix.
    std::set<std::string> names;
    for (const auto& entry : table)
    {
        check(entry.name != nullptr && entry.name[0] != '\0',
              "table entry has a non-empty name");
        check(entry.suffix != nullptr && entry.suffix[0] != '\0',
              std::string("table entry '") + entry.name +
              "' has a non-empty suffix");
        check(names.insert(std::string(entry.name)).second,
              std::string("table name '") + entry.name + "' is unique");
    }
}

// -----------------------------------------------------------------------------
// test_round_trip
// -----------------------------------------------------------------------------
void test_round_trip()
{
    for (lsdtt::polyfit_metric metric : every_metric())
    {
        const char* name = lsdtt::polyfit_metric_name(metric);

        check(std::string(name) != "unknown",
              "every enum value has a name (one resolved to 'unknown')");

        bool ok = false;
        try
        {
            ok = (lsdtt::polyfit_metric_from_string(name) == metric);
        }
        catch (const std::exception& e)
        {
            std::printf("FAIL: round trip threw for '%s': %s\n", name, e.what());
            ++g_failures;
            continue;
        }

        check(ok, std::string("round trip of '") + name + "'");
    }
}

// -----------------------------------------------------------------------------
// test_metrics_missing_from_basin_copy
//
// The six names that the lsdtt-basin-surface-metrics copy of the table did not
// carry. They must resolve here.
// -----------------------------------------------------------------------------
void test_metrics_missing_from_basin_copy()
{
    const std::vector<std::pair<std::string, lsdtt::polyfit_metric>> missing =
    {
        {"difference_curvature",        lsdtt::polyfit_metric::difference_curvature},
        {"horizontal_excess_curvature", lsdtt::polyfit_metric::horizontal_excess_curvature},
        {"vertical_excess_curvature",   lsdtt::polyfit_metric::vertical_excess_curvature},
        {"accumulation_curvature",      lsdtt::polyfit_metric::accumulation_curvature},
        {"ring_curvature",              lsdtt::polyfit_metric::ring_curvature},
        {"rotor",                       lsdtt::polyfit_metric::rotor}
    };

    for (const auto& item : missing)
    {
        try
        {
            check(lsdtt::polyfit_metric_from_string(item.first) == item.second,
                  "previously-missing metric '" + item.first + "' resolves");
        }
        catch (const std::exception& e)
        {
            std::printf("FAIL: previously-missing metric '%s' threw: %s\n",
                        item.first.c_str(), e.what());
            ++g_failures;
        }
    }
}

// -----------------------------------------------------------------------------
// test_unknown_name_throws
// -----------------------------------------------------------------------------
void test_unknown_name_throws()
{
    bool threw = false;
    std::string message;

    try
    {
        (void)lsdtt::polyfit_metric_from_string("nonsense");
    }
    catch (const std::runtime_error& e)
    {
        threw   = true;
        message = e.what();
    }

    check(threw, "polyfit_metric_from_string(\"nonsense\") throws");
    check_string(message, "Unknown polyfit metric: nonsense",
                 "unknown-name message is byte-identical to the CLI message");
}

// -----------------------------------------------------------------------------
// test_output_suffix
// -----------------------------------------------------------------------------
void test_output_suffix()
{
    check_string(lsdtt::polyfit_metric_output_suffix(lsdtt::polyfit_metric::slope),
                 "_wavelet_polyfit_slope.tif",
                 "default args, slope");

    check_string(lsdtt::polyfit_metric_output_suffix(
                     lsdtt::polyfit_metric::fitted_surface),
                 "_wavelet_polyfit_fitted_surface.tif",
                 "default args, fitted_surface");

    check_string(lsdtt::polyfit_metric_output_suffix(
                     lsdtt::polyfit_metric::rotor),
                 "_wavelet_polyfit_rotor.tif",
                 "default args, rotor");

    check_string(lsdtt::polyfit_metric_output_suffix(
                     lsdtt::polyfit_metric::slope, "direct", ".tif"),
                 "_direct_polyfit_slope.tif",
                 "(direct, .tif), slope");

    check_string(lsdtt::polyfit_metric_output_suffix(
                     lsdtt::polyfit_metric::plan_curvature, "direct", ".tif"),
                 "_direct_polyfit_plan_curvature.tif",
                 "(direct, .tif), plan_curvature");

    check_string(lsdtt::polyfit_metric_output_suffix(
                     lsdtt::polyfit_metric::slope, "wavelet", ".png"),
                 "_wavelet_polyfit_slope.png",
                 "(wavelet, .png), slope");

    check_string(lsdtt::polyfit_metric_output_suffix(
                     lsdtt::polyfit_metric::ring_curvature, "wavelet", ".png"),
                 "_wavelet_polyfit_ring_curvature.png",
                 "(wavelet, .png), ring_curvature");

    check_string(lsdtt::polyfit_metric_output_suffix(
                     lsdtt::polyfit_metric::slope, "direct", ".png"),
                 "_direct_polyfit_slope.png",
                 "(direct, .png), slope");

    // Every metric must produce a suffix built from its own name.
    for (const auto& entry : lsdtt::all_polyfit_metrics())
    {
        const std::string expected =
            std::string("_wavelet_polyfit_") + entry.name + ".tif";
        check_string(lsdtt::polyfit_metric_output_suffix(entry.metric),
                     expected,
                     std::string("suffix pattern for '") + entry.name + "'");
    }
}

// -----------------------------------------------------------------------------
// test_parse_list
// -----------------------------------------------------------------------------
void test_parse_list()
{
    check(lsdtt::parse_polyfit_metric_list("").empty(),
          "empty input gives an empty list");

    check(lsdtt::parse_polyfit_metric_list("   ").empty(),
          "whitespace-only input gives an empty list");

    const std::vector<lsdtt::polyfit_metric> all =
        lsdtt::parse_polyfit_metric_list("all");

    check(all.size() == lsdtt::polyfit_metric_count,
          "\"all\" expands to 27 metrics; got " + std::to_string(all.size()));

    if (all.size() == lsdtt::all_polyfit_metrics().size())
    {
        for (std::size_t i = 0; i < all.size(); ++i)
        {
            check(all[i] == lsdtt::all_polyfit_metrics()[i].metric,
                  "\"all\" preserves table order at index " +
                  std::to_string(i));
        }
    }

    const std::vector<lsdtt::polyfit_metric> dup =
        lsdtt::parse_polyfit_metric_list("slope,slope");
    check(dup.size() == 1 && dup[0] == lsdtt::polyfit_metric::slope,
          "duplicates are removed");

    const std::vector<lsdtt::polyfit_metric> two =
        lsdtt::parse_polyfit_metric_list("slope,rotor");
    check(two.size() == 2 &&
          two[0] == lsdtt::polyfit_metric::slope &&
          two[1] == lsdtt::polyfit_metric::rotor,
          "first-seen order is preserved");

    const std::vector<lsdtt::polyfit_metric> spaced =
        lsdtt::parse_polyfit_metric_list(" slope , , rotor ,");
    check(spaced.size() == 2 &&
          spaced[0] == lsdtt::polyfit_metric::slope &&
          spaced[1] == lsdtt::polyfit_metric::rotor,
          "whitespace trimmed and empty tokens dropped");

    const std::vector<lsdtt::polyfit_metric> all_then_slope =
        lsdtt::parse_polyfit_metric_list("all,slope");
    check(all_then_slope.size() == lsdtt::polyfit_metric_count,
          "\"all,slope\" does not double-count slope");

    bool threw = false;
    std::string message;
    try
    {
        (void)lsdtt::parse_polyfit_metric_list("slope,nonsense");
    }
    catch (const std::runtime_error& e)
    {
        threw   = true;
        message = e.what();
    }

    check(threw, "an unknown token throws");
    check_string(message, "Unknown polyfit metric: nonsense",
                 "unknown-token message is byte-identical to the CLI message");
}

} // namespace

int main()
{
    test_table_completeness();
    test_round_trip();
    test_metrics_missing_from_basin_copy();
    test_unknown_name_throws();
    test_output_suffix();
    test_parse_list();

    if (g_failures != 0)
    {
        std::printf("test_surface_metrics_polyfit_metric_table: %d FAILURE(S)\n",
                    g_failures);
        return 1;
    }

    std::printf("test_surface_metrics_polyfit_metric_table: all checks passed\n");
    return 0;
}
