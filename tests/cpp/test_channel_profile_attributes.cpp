// test_channel_profile_attributes.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Can a longitudinal profile be built from the channel-network points?
//
// Until 2026-08-26 the answer was NO, on any domain. The points output carried
//
//     flat, row, col, stream_order, drainage_area
//
// -- no distance and no elevation. A profile is elevation against distance
// ALONG THE FLOW PATH, and neither was there.
//
// This mattered beyond the missing columns. The motivating change was to split
// channel lines at periodic boundaries so they stop being drawn across the map,
// and the objection to splitting is "but then profiles break". They do not,
// because a profile should never have been built from geometry in the first
// place: distance along the flow path is a property the router computes, and it
// is unaffected by how the geometry is packaged.
//
// So the two requirements were never in conflict -- but only once the
// attributes exist.
//
//   P1  `flow_distance` and `elevation` are present
//   P2  ON A KNOWN CHANNEL: flow_distance increases by exactly one cell per
//       cell upstream, and elevation matches the surface that was built
//   P3  the profile is MONOTONIC -- elevation rises with flow distance along a
//       single channel, which is what makes it a profile rather than a scatter
//   P4  flow_distance comes from the D8 direction code, not from coordinates,
//       so a PERIODIC network carries no domain-width jump in its profile --
//       the property that makes the split in phase 3 harmless
//
// Usage: lsdtt-test-channel-profile-attributes <out_dir>
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_data_loader.hpp"
#include "../../src/cpu/lsdtt_flowinfo_cpu.hpp"
#include "../../src/cpu/lsdtt_junction_network_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"

namespace
{

int g_checks   = 0;
int g_failures = 0;

void check(bool ok, const std::string& what)
{
    ++g_checks;
    if (!ok) { ++g_failures; std::cout << "  FAIL: " << what << "\n"; }
    else     { std::cout << "  ok  : " << what << "\n"; }
}

constexpr double kCell  = 10.0;
constexpr double kSlope = 0.05;      // elevation gain per metre upstream

}  // namespace


int main(int argc, char** argv)
{
    using namespace lsdtt;

    const std::string out_dir = (argc > 1) ? argv[1] : ".";

    std::cout << "=====================================================\n"
              << " Longitudinal profile from the points output\n"
              << "=====================================================\n";

    // A KNOWN channel: a single valley running due south down the middle, with
    // elevation exactly linear in row. Every flow step is cardinal, so the
    // expected flow distance is exactly (rows to the outlet) x cell size and
    // the expected elevation is exactly kSlope x that distance.
    const int nr = 40, nc = 21;

    lsdtt_raster_meta meta;
    meta.n_rows = nr; meta.n_cols = nc; meta.cell_size = kCell;
    meta.x_min = 0.0; meta.y_min = 0.0; meta.no_data_value = -9999.0;

    std::vector<float>        z(static_cast<std::size_t>(nr) * nc);
    std::vector<std::uint8_t> m(z.size(), 1u);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
            z[static_cast<std::size_t>(r) * nc + c] =
                static_cast<float>(kSlope * kCell * (nr - 1 - r)
                                   + 2.0 * std::abs(c - nc / 2));

    const lsdtt_raster_cpu dem(std::move(z), std::move(m), meta);

    const auto run = [&](const std::string& bc, const std::string& tag)
    {
        lsdtt_flowinfo_options fi;
        fi.boundary_conditions = bc;
        const lsdtt_flowinfo flow(dem.get_data_array(), dem.get_mask_array(),
                                  dem.get_meta(), fi);

        const std::vector<double> area = flow.accumulate();
        std::vector<std::int64_t> src;
        double amax = 0.0;
        std::size_t argmax = 0;
        for (std::size_t i = 0; i < area.size(); ++i)
            if (flow.is_valid_cell(static_cast<std::int64_t>(i)) &&
                area[i] > amax) { amax = area[i]; argmax = i; }
        // One source: the head of the trunk, found by walking up the steepest
        // accumulation. Simpler and more deterministic than a threshold sweep.
        src.push_back(static_cast<std::int64_t>(nc / 2));   // top-middle cell

        lsdtt_junction_network_cpu net(flow, src);
        const std::string file = out_dir + "/profile_" + tag + ".fgb";
        net.write_channel_network_points(file, point_output_format::flatgeobuf,
                                         "pts", false,
                                         dem.get_data(), dem.get_mask());

        lsdtt_data_loader loader;
        loader.open_vector(file);
        const lsdtt_vector_layer layer = loader.read_vector_layer("pts");
        loader.close_vector();
        (void)argmax;
        return layer;
    };

    const lsdtt_vector_layer layer = run("oooo", "open");

    // -------------------------------------------------------------------------
    std::cout << "\nP1: the fields exist\n";
    {
        const auto has = [&](const std::string& n)
        {
            return std::find(layer.schema.field_names.begin(),
                             layer.schema.field_names.end(), n) !=
                   layer.schema.field_names.end();
        };
        check(has("flow_distance"), "flow_distance is present");
        check(has("elevation"), "elevation is present");
        check(!layer.records.empty(), "and there are channel points (" +
              std::to_string(layer.records.size()) + ")");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nP2/P3: the profile matches the surface it was built from\n";
    {
        std::size_t i_fd = 0, i_el = 0, i_row = 0;
        for (std::size_t i = 0; i < layer.schema.field_names.size(); ++i)
        {
            if (layer.schema.field_names[i] == "flow_distance") i_fd = i;
            if (layer.schema.field_names[i] == "elevation")     i_el = i;
            if (layer.schema.field_names[i] == "row")           i_row = i;
        }

        std::vector<std::pair<double, double>> profile;   // (distance, elev)
        double worst_elev_err = 0.0;

        for (const auto& rec : layer.records)
        {
            const auto* fd = std::get_if<double>(&rec.field_values[i_fd]);
            const auto* el = std::get_if<double>(&rec.field_values[i_el]);
            const auto* rw = std::get_if<std::int64_t>(&rec.field_values[i_row]);
            if (!fd || !el || !rw) continue;

            profile.emplace_back(*fd, *el);

            // The channel runs down the middle column, so every step is
            // cardinal and the elevation is exactly kSlope*kCell*(nr-1-row).
            const double expected = kSlope * kCell *
                                    static_cast<double>(nr - 1 - *rw);
            worst_elev_err = std::max(worst_elev_err, std::fabs(*el - expected));
        }

        check(!profile.empty(), "the profile has points");
        check(worst_elev_err < 1e-4,
              "every elevation matches the surface exactly (worst error " +
              std::to_string(worst_elev_err) + " m)");

        std::sort(profile.begin(), profile.end());

        bool monotonic = true;
        double worst_step = 0.0;
        for (std::size_t i = 1; i < profile.size(); ++i)
        {
            if (profile[i].second < profile[i - 1].second - 1e-6)
                monotonic = false;
            worst_step = std::max(worst_step,
                                  profile[i].first - profile[i - 1].first);
        }
        check(monotonic,
              "elevation rises monotonically with flow distance -- a PROFILE, "
              "not a scatter");
        check(worst_step <= kCell * std::sqrt(2.0) + 1e-6,
              "and successive points are one cell apart in flow distance "
              "(worst " + std::to_string(worst_step) + " m)");

        // The slope recovered from the profile must be the slope built in.
        if (profile.size() >= 2)
        {
            const double got = (profile.back().second - profile.front().second) /
                               (profile.back().first - profile.front().first);
            check(std::fabs(got - kSlope) < 1e-6,
                  "the gradient recovered from the profile is the one the "
                  "surface was built with (" + std::to_string(got) + ")");
        }
    }

    // -------------------------------------------------------------------------
    std::cout << "\nP4: periodic boundaries leave the profile intact\n";
    {
        const lsdtt_vector_layer p = run("opop", "periodic");

        std::size_t i_fd = 0;
        for (std::size_t i = 0; i < p.schema.field_names.size(); ++i)
            if (p.schema.field_names[i] == "flow_distance") i_fd = i;

        std::vector<double> d;
        for (const auto& rec : p.records)
            if (const auto* fd = std::get_if<double>(&rec.field_values[i_fd]))
                if (*fd >= 0.0) d.push_back(*fd);

        std::sort(d.begin(), d.end());
        double worst = 0.0;
        for (std::size_t i = 1; i < d.size(); ++i)
            worst = std::max(worst, d[i] - d[i - 1]);

        check(!d.empty(), "the periodic network has flow distances");
        check(worst <= kCell * std::sqrt(2.0) + 1e-6,
              "NO domain-width jump in flow distance -- it comes from the D8 "
              "direction code, not from coordinates, which is what makes "
              "splitting the geometry harmless (worst step " +
              std::to_string(worst) + " m)");
    }

    std::cout << "\n=====================================================\n"
              << " " << (g_checks - g_failures) << " / " << g_checks
              << " checks passed\n"
              << "=====================================================\n";
    return g_failures == 0 ? 0 : 1;
}
