// test_slope_stability_hazard_area.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Coseismic hazard area, SHALRUN-EQ (src/cpu/lsdtt_slope_stability_hazard_area_cpu):
// Milledge et al. (2019) Fig. 2, exact profiles, the limit of ordinary MFD
// accumulation, conservation, monotonicity in the thresholds, the Ah = 0 set
// against a brute-force search, nodata, underflow and refusals.
//
// Fig. 2e shows 0 at two cells (row 4, columns 1 and 3) that receive debris along
// a path with slopes of 13 and 11 degrees, above its 10 degree stopping angle.
// Eqs. 6-9 give 30 there, so these two cells are corrected here.
//
// Usage:
//   lsdtt-test-slope-stability-hazard-area
//   lsdtt-test-slope-stability-hazard-area --write-cli-fixture <dir>
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_data_loader.hpp"
#include "../../src/cpu/lsdtt_flowgraph_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../../src/cpu/lsdtt_slope_stability_hazard_area_cpu.hpp"

using namespace lsdtt;

namespace
{

int g_failures = 0;
int g_checks = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_checks;
    if (!ok) ++g_failures;
    std::cout << (ok ? "  PASS  " : "  FAIL  ") << name;
    if (!detail.empty()) std::cout << "  (" << detail << ")";
    std::cout << "\n";
}

lsdtt_raster_meta grid(int nr, int nc, double cell = 10.0, double nodata = -9999.0)
{
    lsdtt_raster_meta m;
    m.n_rows = nr;
    m.n_cols = nc;
    m.x_min = 500000.0;
    m.y_min = 4000000.0;
    m.cell_size = cell;
    m.no_data_value = nodata;
    m.set_crs_from_epsg(32633);
    return m;
}

lsdtt_raster_cpu raster(std::vector<float> v, const lsdtt_raster_meta& m)
{
    std::vector<std::uint8_t> mask(v.size(), 1u);
    for (std::size_t i = 0; i < v.size(); ++i)
        if (v[i] == static_cast<float>(m.no_data_value)) mask[i] = 0u;
    return lsdtt_raster_cpu(std::move(v), std::move(mask), m);
}

// Elevation depends on the column; gradient[c] is the fall from c to c + 1.
std::vector<float> profile(int nr, const std::vector<double>& gradient, double cell)
{
    const int nc = static_cast<int>(gradient.size()) + 1;
    std::vector<double> z(nc, 0.0);
    for (int c = nc - 2; c >= 0; --c) z[c] = z[c + 1] + gradient[c] * cell;
    std::vector<float> v(static_cast<std::size_t>(nr) * nc);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c) v[static_cast<std::size_t>(r) * nc + c] = static_cast<float>(z[c]);
    return v;
}

std::vector<double> fixture_gradient()
{
    std::vector<double> g;
    g.insert(g.end(), 20, 1.0);   // sources
    g.insert(g.end(), 5, 0.5);    // transit
    g.insert(g.end(), 14, 0.05);  // stopping
    return g;
}

std::vector<float> rough_surface(int nr, int nc, std::uint64_t seed, const lsdtt_raster_meta& m)
{
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> noise(0.0, 4.0);
    std::vector<float> v(static_cast<std::size_t>(nr) * nc);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
            v[static_cast<std::size_t>(r) * nc + c] = static_cast<float>(
                1000.0 - 3.0 * c + 45.0 * std::sin(r * 0.35) * std::cos(c * 0.22) + noise(rng));
    return fill_pits(raster(std::move(v), m), 0.0001f, "oooo", 1).get_data();
}

void test_published_figure()
{
    std::cout << "Milledge et al. (2019) Fig. 2\n";
    const std::vector<double> slope{19, 43, 45, 43, 19, 35, 49, 40, 49, 35, 47, 35, 9, 35, 47,
                                    40, 13, 9,  13, 40, 27, 11, 9,  11, 27};
    const std::vector<double> expected{0,  30, 30,  30, 0,  0,  30, 120, 30, 0,  30, 0, 180,
                                       0,  30, 30,  30, 0,  30, 30, 30,  30, 0,  30, 30};
    // Receivers as drawn in the figure (the last row leaves the panel).
    const std::vector<int> to{6,  7,  7,  7,  8,  11, 12, 12, 12, 13, 16, 16, 17,
                              18, 18, 20, 21, 22, 23, 24, -1, -1, -1, -1, -1};
    std::vector<std::int64_t> order(25), offsets(26, 0), receivers;
    std::vector<float> weights;
    std::vector<std::uint8_t> passes(25);
    std::vector<double> seed(25);
    for (int i = 0; i < 25; ++i)
    {
        order[i] = i;
        if (to[i] >= 0)
        {
            receivers.push_back(to[i]);
            weights.push_back(1.0f);
        }
        offsets[i + 1] = static_cast<std::int64_t>(receivers.size());
        passes[i] = slope[i] >= 10.0 ? 1 : 0;
        seed[i] = slope[i] >= 40.0 ? 30.0 : 0.0;
    }
    const std::vector<double> a = route_hazard_area(order, offsets, receivers, weights, passes, seed);
    int matches = 0;
    for (int i = 0; i < 25; ++i) matches += a[i] == expected[i] ? 1 : 0;
    check("Eqs. 6-9 on the figure's graph", matches == 25, std::to_string(matches) + "/25");
    check("the first gentle cell is hit", a[12] == 180.0);
    check("and passes nothing on", a[17] == 0.0);
}

void test_exact_profile()
{
    std::cout << "Exact profile (periodic north-south)\n";
    const lsdtt_raster_meta m = grid(20, 40);
    coseismic_hazard_area_options o;
    o.boundary_conditions = "popo";
    const auto r = compute_coseismic_hazard_area(raster(profile(20, fixture_gradient(), 10.0), m), o);
    const std::vector<float>& a = r.hazard_area.get_data();
    double worst = 0.0;
    for (int row = 0; row < 20; ++row)
        for (int c = 0; c < 40; ++c)
        {
            const double want = c < 20 ? (c + 1) * 10.0 : (c <= 25 ? 200.0 : 0.0);
            worst = std::max(worst, std::abs(a[static_cast<std::size_t>(row) * 40 + c] - want));
        }
    check("Ah by column", worst < 1e-3, "worst error " + std::to_string(worst));
    check("sources", r.n_source == 400, std::to_string(r.n_source));
    check("source area", r.source_area_m2 == 40000.0);
    check("zero-hazard cells", r.n_zero_hazard == 14 * 20, std::to_string(r.n_zero_hazard));
    check("steepest slope", std::abs(r.steepest_slope.get_data()[5] - 1.0f) < 1e-5);
    check("zone: source", r.zone.get_data()[5] == k_hazard_zone_source);
    check("zone: transit", r.zone.get_data()[22] == k_hazard_zone_transit);
    check("zone: stopping", r.zone.get_data()[30] == k_hazard_zone_stopping);
}

void test_diagonal_slope()
{
    std::cout << "Gradient uses the diagonal length\n";
    const lsdtt_raster_meta m = grid(3, 3);
    std::vector<float> z(9, 100.0f);
    z[0] = 114.2f;  // its lowest neighbour is the diagonal centre cell
    z[1] = z[3] = 110.0f;
    z[2] = z[5] = z[6] = z[7] = z[8] = 101.0f;
    z[4] = 90.0f;
    coseismic_hazard_area_options o;
    const auto r = compute_coseismic_hazard_area(raster(z, m), o);
    const double want = (114.2 - 90.0) / (10.0 * std::sqrt(2.0));
    check("steepest gradient", std::abs(r.steepest_slope.get_data()[0] - want) < 1e-4);
}

void test_limits_and_properties()
{
    std::cout << "Limits and properties on a rough surface\n";
    const lsdtt_raster_meta m = grid(70, 90);
    const lsdtt_raster_cpu surface = raster(rough_surface(70, 90, 11, m), m);

    // Every cell a source, nothing stops: ordinary MFD accumulation / cell width.
    coseismic_hazard_area_options all;
    all.initiation_slope = 2e-6;
    all.stopping_slope = 1e-6;
    const auto ra = compute_coseismic_hazard_area(surface, all);
    lsdtt_flowgraph_options fo;
    fo.method = lsdtt_flowgraph_method::mfd;
    const lsdtt_flowgraph_cpu graph(surface.get_data_array(), surface.get_mask_array(), m, fo);
    const std::vector<double> area = graph.accumulate();
    double worst = 0.0;
    for (std::size_t i = 0; i < area.size(); ++i)
    {
        if (ra.steepest_slope.get_data()[i] < 1e-6f) continue;
        worst = std::max(worst, std::abs(ra.hazard_area.get_data()[i] - area[i] / 10.0) / (area[i] / 10.0));
    }
    check("thresholds off = MFD area / cell width", worst < 1e-5, "worst relative " + std::to_string(worst));

    // Conservation: all source area ends on a cell that passes nothing on.
    coseismic_hazard_area_options o;
    o.initiation_slope = 0.9;
    o.stopping_slope = 0.2;
    const auto r = compute_coseismic_hazard_area(surface, o);
    const std::vector<std::int64_t>& off = graph.receiver_offsets();
    double ended = 0.0;
    for (std::size_t i = 0; i < area.size(); ++i)
    {
        const bool terminal = r.steepest_slope.get_data()[i] < 0.2f || off[i + 1] == off[i];
        if (terminal) ended += r.hazard_area.get_data()[i];
    }
    const double started = static_cast<double>(r.n_source) * 10.0;
    check("source area is conserved", r.n_source > 0 && std::abs(ended - started) / started < 1e-4,
          std::to_string(ended) + " of " + std::to_string(started));

    // Monotonic in both thresholds.
    coseismic_hazard_area_options higher_i = o, higher_s = o;
    higher_i.initiation_slope = 1.1;
    higher_s.stopping_slope = 0.4;
    const auto ri = compute_coseismic_hazard_area(surface, higher_i);
    const auto rs = compute_coseismic_hazard_area(surface, higher_s);
    bool mono_i = true, mono_s = true;
    for (std::size_t i = 0; i < area.size(); ++i)
    {
        const float a = r.hazard_area.get_data()[i];
        mono_i = mono_i && ri.hazard_area.get_data()[i] <= a * 1.000001f;
        mono_s = mono_s && rs.hazard_area.get_data()[i] <= a * 1.000001f;
    }
    check("raising the initiation slope never raises Ah", mono_i);
    check("raising the stopping slope never raises Ah", mono_s);

    // The Ah = 0 set against a brute-force search from every source.
    const std::vector<std::int64_t>& rec = graph.receiver_indices();
    const std::vector<float>& w = graph.receiver_weights();
    std::vector<std::uint8_t> brute(area.size(), 0), seen(area.size(), 0);
    std::vector<std::int64_t> stack;
    for (std::size_t s = 0; s < area.size(); ++s)
    {
        if (r.steepest_slope.get_data()[s] < 0.9f) continue;
        std::fill(seen.begin(), seen.end(), 0);
        stack.assign(1, static_cast<std::int64_t>(s));
        seen[s] = 1;
        while (!stack.empty())
        {
            const std::int64_t i = stack.back();
            stack.pop_back();
            brute[i] = 1;
            if (r.steepest_slope.get_data()[i] < 0.2f) continue;
            for (std::int64_t k = off[i]; k < off[i + 1]; ++k)
                if (w[k] > 0.0f && !seen[rec[k]])
                {
                    seen[rec[k]] = 1;
                    stack.push_back(rec[k]);
                }
        }
    }
    std::int64_t differ = 0;
    for (std::size_t i = 0; i < area.size(); ++i)
        differ += (r.hazard_area.get_data()[i] > 0.0f) != (brute[i] != 0) ? 1 : 0;
    check("Ah > 0 exactly where a source can reach", differ == 0, std::to_string(differ) + " differ");
}

void test_edges_nodata_and_underflow()
{
    std::cout << "Edges, nodata and underflow\n";
    {
        const lsdtt_raster_meta m = grid(9, 9);
        std::vector<float> z = profile(9, std::vector<double>(8, 1.0), 10.0);
        z[40] = -9999.0f;
        const auto r = compute_coseismic_hazard_area(raster(z, m));
        check("nodata stays nodata", r.hazard_area.get_mask()[40] == 0);
        check("nodata value of the output", r.hazard_area.get_data()[40] == -9999.0f);
        check("edge cell is flagged", r.edge_influenced.get_data()[0] == 1.0f);
        check("flag follows debris", r.n_edge_influenced > 0 && r.n_edge_influenced < r.n_valid + 1);
    }
    {
        // A DEM whose nodata value is 0 must still keep Ah = 0 as a value.
        const lsdtt_raster_meta m = grid(6, 20, 10.0, 0.0);
        std::vector<double> g(19, 0.05);
        const std::vector<float> z = profile(6, g, 10.0);
        std::vector<float> lifted(z);
        for (float& x : lifted) x += 100.0f;
        const auto r = compute_coseismic_hazard_area(raster(lifted, m));
        check("Ah = 0 is not nodata", r.hazard_area.get_meta().no_data_value == -9999.0 &&
                                          r.hazard_area.get_mask()[10] == 1 && r.hazard_area.get_data()[10] == 0.0f);
    }
    {
        // A chain whose product underflows Float64 is still "reached".
        const int n = 40;
        std::vector<std::int64_t> order(n), offsets(n + 1), receivers;
        std::vector<float> weights;
        for (int i = 0; i < n; ++i)
        {
            order[i] = i;
            if (i + 1 < n)
            {
                receivers.push_back(i + 1);
                weights.push_back(1e-30f);
            }
            offsets[i + 1] = static_cast<std::int64_t>(receivers.size());
        }
        std::vector<double> seed(n, 0.0);
        seed[0] = 30.0;
        std::vector<std::uint8_t> reached;
        const auto a = route_hazard_area(order, offsets, receivers, weights, std::vector<std::uint8_t>(n, 1), seed, &reached);
        check("underflowing chain has zero value", a[n - 1] == 0.0);
        check("but is reached", reached[n - 1] == 1);
    }
}

void test_refusals()
{
    std::cout << "Refusals\n";
    const auto refused = [](double i, double s)
    {
        coseismic_hazard_area_options o;
        o.initiation_slope = i;
        o.stopping_slope = s;
        try
        {
            validate_coseismic_hazard_area_options(o);
        }
        catch (const std::invalid_argument&)
        {
            return true;
        }
        return false;
    };
    check("stopping = 0", refused(0.8, 0.0));
    check("stopping < 0", refused(0.8, -0.1));
    check("initiation = stopping", refused(0.5, 0.5));
    check("initiation < stopping", refused(0.2, 0.5));
    check("NaN", refused(std::numeric_limits<double>::quiet_NaN(), 0.1));
    check("default accepted", !refused(0.839, 0.176));
}

void write_cli_fixture(const std::string& dir)
{
    std::filesystem::create_directories(dir);
    const lsdtt_raster_meta m = grid(20, 40);
    lsdtt_data_loader loader;
    loader.write_full(profile(20, fixture_gradient(), 10.0), m, dir + "/hazard_dem.tif");
}

}  // namespace

int main(int argc, char** argv)
{
    try
    {
        if (argc == 3 && std::string(argv[1]) == "--write-cli-fixture")
        {
            write_cli_fixture(argv[2]);
            return 0;
        }
        test_published_figure();
        test_exact_profile();
        test_diagonal_slope();
        test_limits_and_properties();
        test_edges_nodata_and_underflow();
        test_refusals();
    }
    catch (const std::exception& e)
    {
        std::cout << "  FAIL  unexpected exception: " << e.what() << "\n";
        ++g_failures;
    }
    std::cout << g_checks << " checks, " << g_failures << " failed\n";
    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
