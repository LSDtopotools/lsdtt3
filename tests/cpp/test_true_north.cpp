// test_true_north.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Grid convergence and the conversion of grid bearings to true north
// (lsdtt_true_north.hpp), the core of bearings_relative_to_true_north.
//
//   C*  convergence against pyproj 3.7 Proj.get_factors().meridian_convergence
//       (lsdtt3's c is minus pyproj's: pyproj gives grid north from true
//       north): UTM 33N at 45 N on the central meridian and both zone edges,
//       UTM 33S at 30 S, EPSG:3413 at 80 N, an oblique frame (alpha 30 at
//       38.4 N) at its centre and 50 km away; a geographic CRS is 0.
//   W*  wrapping of bearings (degrees, radians) and axial orientations;
//       non-finite values pass through.
//   P*  a plane dipping to a known true azimuth: polyfit aspect is the grid
//       bearing; converted, it is the true azimuth, in UTM and in frames at
//       alpha = -60, 0, 30, 89.
//   R*  raster conversion: lattice interpolation against an exact per-cell
//       conversion, nodata kept, axial kind.
//   M*  start-up message and converter: silent below 1 degree with the option
//       off, worded for a frame; no CRS refused.
//   N*  channel-network records: a link converted at its representative
//       point (NaN without one), directions left as grid components; a
//       junction basin's outlet -> centre-of-mass bearing at the midpoint;
//       fraction_on_d8_multiples still taken from the grid bearings.
//
// Usage: lsdtt-test-true-north
// -----------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <cpl_conv.h>
#include <ogr_spatialref.h>

#include "lsdtt_channel_bearing_statistics_cpu.hpp"
#include "lsdtt_channel_network_analysis_cpu.hpp"
#include "lsdtt_gdal_init.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_surface_metrics_cpu.hpp"
#include "lsdtt_true_north.hpp"

namespace {

int checks = 0;
int failures = 0;

void check(const std::string& name, bool ok)
{
    ++checks;
    if (!ok) ++failures;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name << '\n';
}

std::string export_wkt(OGRSpatialReference& srs)
{
    char* raw = nullptr;
    srs.exportToWkt(&raw);
    std::string wkt(raw ? raw : "");
    CPLFree(raw);
    return wkt;
}

std::string wkt_of_epsg(int code)
{
    OGRSpatialReference srs;
    srs.importFromEPSG(code);
    return export_wkt(srs);
}

std::string frame_wkt(double alpha, double lat0 = 38.4, double lonc = -110.0)
{
    OGRSpatialReference srs;
    const std::string p = "+proj=omerc +lat_0=" + std::to_string(lat0) + " +lonc="
        + std::to_string(lonc) + " +alpha=" + std::to_string(alpha)
        + " +gamma=90 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs";
    srs.importFromProj4(p.c_str());
    return export_wkt(srs);
}

// Smallest difference between two angles on a circle of the given period.
double angle_diff(double a, double b, double period = 360.0)
{
    double d = std::fmod(std::fabs(a - b), period);
    return std::min(d, period - d);
}

constexpr double deg = 3.14159265358979323846 / 180.0;

void check_convergence()
{
    using lsdtt::lsdtt_grid_convergence;
    struct ref { const char* name; int epsg; double x, y, pyproj_mc; };
    const ref refs[] = {
        {"C1 UTM 33N 45N 18E (east zone edge)", 32633, 736446.026, 4987329.505, 2.122300},
        {"C2 UTM 33N 45N 12E (west zone edge)", 32633, 263553.974, 4987329.505, -2.122300},
        {"C3 UTM 33N 45N central meridian",     32633, 500000.000, 4982950.400, 0.0},
        {"C4 UTM 33S 30S 18E",                  32733, 789409.653, 6677424.096, -1.501044},
        {"C5 EPSG:3413 80N 0E",                 3413,  767861.606, -767861.606, 45.0},
    };
    for (const ref& r : refs)
    {
        lsdtt_grid_convergence c(wkt_of_epsg(r.epsg));
        const double got = c.degrees_at(r.x, r.y);
        std::cout << "      " << r.name << ": c = " << got << " deg (pyproj "
                  << r.pyproj_mc << ")\n";
        check(r.name, std::fabs(got + r.pyproj_mc) < 1e-3);
    }
    lsdtt_grid_convergence f(frame_wkt(30.0));
    const double c0 = f.degrees_at(0.0, 0.0);
    const double c1 = f.degrees_at(41107.361, -26560.690);
    std::cout << "      frame alpha 30: c(centre) = " << c0 << ", c(-109.5, 38.6) = " << c1 << '\n';
    check("C6 frame alpha 30: true north 60 deg clockwise of grid up at the centre",
          std::fabs(c0 - 60.0) < 1e-3);
    check("C7 frame alpha 30 at (-109.5, 38.6): pyproj -59.687971",
          std::fabs(c1 - 59.687971) < 1e-3);
    check("C8 along the feature (grid 90) is true 30 at the centre",
          std::fabs(f.to_true_north(90.0, lsdtt::lsdtt_angle_kind::bearing_degrees, 0, 0) - 30.0)
              < 1e-3);
    lsdtt_grid_convergence g(wkt_of_epsg(4326));
    check("C9 geographic CRS: convergence 0", g.is_geographic() && g.degrees_at(-110, 38) == 0.0);
}

void check_wrapping()
{
    using lsdtt::grid_angle_to_true_north;
    using K = lsdtt::lsdtt_angle_kind;
    check("W1 bearing 0.5, c 1 -> 359.5",
          std::fabs(grid_angle_to_true_north(0.5, K::bearing_degrees, 1.0) - 359.5) < 1e-12);
    check("W2 bearing 350, c -20 -> 10",
          std::fabs(grid_angle_to_true_north(350.0, K::bearing_degrees, -20.0) - 10.0) < 1e-12);
    check("W3 axial 10, c 60 -> 130",
          std::fabs(grid_angle_to_true_north(10.0, K::axial_degrees, 60.0) - 130.0) < 1e-12);
    check("W4 axial 170, c -30 -> 20",
          std::fabs(grid_angle_to_true_north(170.0, K::axial_degrees, -30.0) - 20.0) < 1e-12);
    const double r = grid_angle_to_true_north(0.1, K::bearing_radians, 60.0);
    check("W5 radians 0.1, c 60 deg -> 0.1 - pi/3 + 2 pi",
          std::fabs(r - (0.1 - 60.0 * deg + 2.0 * 3.14159265358979323846)) < 1e-12);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    check("W6 NaN passes through", std::isnan(grid_angle_to_true_north(nan, K::bearing_degrees, 5.0)));
    const double tiny = grid_angle_to_true_north(-1e-17, K::bearing_degrees, 0.0);
    check("W7 a value just below 0 wraps into [0, 360)", tiny >= 0.0 && tiny < 360.0);
}

// A plane dipping to grid bearing g at slope s, on a 41 x 41 grid of 10 m
// cells centred on (cx, cy), in the CRS wkt.
lsdtt::lsdtt_raster_cpu plane(const std::string& wkt, double cx, double cy, double g)
{
    lsdtt::lsdtt_raster_meta m;
    m.n_rows = 41;
    m.n_cols = 41;
    m.cell_size = 10.0;
    m.x_min = cx - 205.0;
    m.y_min = cy - 205.0;
    m.crs_wkt = wkt;
    std::vector<float> z(41 * 41);
    for (int r = 0; r < 41; ++r)
        for (int c = 0; c < 41; ++c)
        {
            const double dx = (c - 20) * 10.0;
            const double dy = (20 - r) * 10.0;
            z[r * 41 + c] = static_cast<float>(
                1000.0 - 0.2 * (dx * std::sin(g * deg) + dy * std::cos(g * deg)));
        }
    return lsdtt::lsdtt_raster_cpu(std::move(z), std::vector<std::uint8_t>(41 * 41, 1), m);
}

void check_planes()
{
    struct site { const char* name; std::string wkt; double cx, cy; };
    std::vector<site> sites = {
        {"UTM 33N east zone edge", wkt_of_epsg(32633), 736446.026, 4987329.505},
        {"frame alpha -60", frame_wkt(-60.0), 0.0, 0.0},
        {"frame alpha 0", frame_wkt(0.0), 0.0, 0.0},
        {"frame alpha 30", frame_wkt(30.0), 0.0, 0.0},
        {"frame alpha 89", frame_wkt(89.0), 0.0, 0.0},
    };
    const double true_azimuth = 200.0;
    for (const site& s : sites)
    {
        lsdtt::lsdtt_grid_convergence conv(s.wkt);
        const double c = conv.degrees_at(s.cx, s.cy);
        const double g = std::fmod(true_azimuth + c + 360.0, 360.0);
        lsdtt::lsdtt_raster_cpu dem = plane(s.wkt, s.cx, s.cy, g);
        const lsdtt::poly_fit_coeffs fit = lsdtt::direct_polyfit(dem, 50.0, 1);
        lsdtt::polyfit_metric_request req;
        req.metrics = {lsdtt::polyfit_metric::aspect_degrees};
        req.n_threads = 1;
        lsdtt::polyfit_metric_result res = lsdtt::compute_polyfit_metrics(fit, req);
        lsdtt::lsdtt_raster_cpu aspect = res.at(lsdtt::polyfit_metric::aspect_degrees);
        const std::size_t centre = 20 * 41 + 20;
        const double grid_aspect = aspect.get_data()[centre];
        lsdtt::convert_angle_raster_to_true_north(aspect.get_data(), aspect.get_meta(),
                                                  lsdtt::lsdtt_angle_kind::bearing_degrees, conv);
        const double true_aspect = aspect.get_data()[centre];
        std::cout << "      " << s.name << ": c = " << c << ", grid aspect " << grid_aspect
                  << ", true aspect " << true_aspect << '\n';
        check(std::string("P grid aspect is the grid bearing, ") + s.name,
              angle_diff(grid_aspect, g) < 1e-3);
        check(std::string("P true aspect is 200, ") + s.name,
              angle_diff(true_aspect, true_azimuth) < 1e-3);
    }
}

void check_raster()
{
    // A 30 m frame raster 1500 x 400 cells (45 x 12 km) around the centre of
    // an alpha 30 frame, filled with grid bearing 90, one nodata cell.
    lsdtt::lsdtt_raster_meta m;
    m.n_rows = 400;
    m.n_cols = 1500;
    m.cell_size = 30.0;
    m.x_min = -22500.0;
    m.y_min = -6000.0;
    m.crs_wkt = frame_wkt(30.0);
    m.no_data_value = -9999.0;
    std::vector<float> v(m.n_pixels(), 90.0f);
    v[123 * 1500 + 456] = -9999.0f;
    lsdtt::lsdtt_grid_convergence conv(m.crs_wkt);
    std::vector<float> fast = v;
    lsdtt::convert_angle_raster_to_true_north(fast, m, lsdtt::lsdtt_angle_kind::bearing_degrees,
                                              conv);
    std::vector<float> exact = v;
    lsdtt::convert_angle_raster_to_true_north(exact, m, lsdtt::lsdtt_angle_kind::bearing_degrees,
                                              conv, 1);
    double max_err = 0.0;
    double min_v = 360.0;
    double max_v = 0.0;
    for (std::size_t i = 0; i < v.size(); ++i)
    {
        if (i == 123u * 1500u + 456u) continue;
        max_err = std::max(max_err, angle_diff(fast[i], exact[i]));
        min_v = std::min(min_v, static_cast<double>(exact[i]));
        max_v = std::max(max_v, static_cast<double>(exact[i]));
    }
    std::cout << "      45 x 12 km frame: true bearing of grid 90 spans " << min_v << " to "
              << max_v << "; lattice (32 cells) vs per cell max |diff| " << max_err << " deg\n";
    check("R1 lattice interpolation within 1e-4 deg of per-cell conversion", max_err < 1e-4);
    check("R2 nodata stays nodata", fast[123 * 1500 + 456] == -9999.0f);
    check("R3 convergence varies across the frame (span > 0.1 deg)", max_v - min_v > 0.1);
    std::vector<float> ax(m.n_pixels(), 170.0f);
    lsdtt::convert_angle_raster_to_true_north(ax, m, lsdtt::lsdtt_angle_kind::axial_degrees, conv);
    bool in_range = true;
    for (float a : ax) in_range = in_range && a >= 0.0f && a < 180.0f;
    check("R4 axial kind stays in [0, 180)", in_range);
    check("R5 axial 170 at the centre is 110",
          std::fabs(ax[200 * 1500 + 750] - 110.0f) < 0.01f);
}

void check_messages()
{
    lsdtt::lsdtt_raster_meta m;
    m.n_rows = 100;
    m.n_cols = 100;
    m.cell_size = 30.0;
    m.x_min = 500000.0 - 1500.0;
    m.y_min = 4982950.4 - 1500.0;
    m.crs_wkt = wkt_of_epsg(32633);
    check("M1 UTM central meridian, option off: no message",
          lsdtt::grid_north_message(m, false).empty());
    const std::string on = lsdtt::grid_north_message(m, true);
    std::cout << "      " << on << '\n';
    check("M2 option on: always reported", on.find("converted to true north") != std::string::npos);
    m.x_min = -1500.0;
    m.y_min = -1500.0;
    m.crs_wkt = frame_wkt(30.0);
    const std::string f = lsdtt::grid_north_message(m, false);
    std::cout << "      " << f << '\n';
    check("M3 frame alpha 30: 60.0 deg anticlockwise",
          f.find("Grid north is 60.0 deg anticlockwise of true north") == 0);
    check("M4 converter is null with the option off", !lsdtt::true_north_converter(m, false));
    m.crs_wkt.clear();
    bool refused = false;
    try { lsdtt::true_north_converter(m, true); }
    catch (const std::invalid_argument&) { refused = true; }
    check("M5 option on, no CRS: refused", refused);
    check("M6 no CRS: no message", lsdtt::grid_north_message(m, true).empty());
}

void check_channel_records()
{
    lsdtt::lsdtt_grid_convergence conv(frame_wkt(30.0));
    const double nan = std::numeric_limits<double>::quiet_NaN();

    std::vector<lsdtt::lsdtt_channel_link_bearing> links(2);
    links[0].representative_x = 0.0;
    links[0].representative_y = 0.0;
    links[0].direction_x = 1.0;
    links[0].direction_y = 0.0;
    links[0].bearing_deg = 90.0;
    links[0].bearing_rad = 90.0 * deg;
    links[1].bearing_deg = 45.0;          // no representative point
    links[1].bearing_rad = 45.0 * deg;
    const auto out = lsdtt::channel_link_bearings_to_true_north(links, conv);
    check("N1 link grid 90 at the frame centre is true 30",
          std::fabs(out[0].bearing_deg - 30.0) < 1e-3
              && std::fabs(out[0].bearing_rad - 30.0 * deg) < 1e-5);
    check("N2 direction_x / direction_y stay grid components",
          out[0].direction_x == 1.0 && out[0].direction_y == 0.0);
    check("N3 a link without a representative point gets NaN",
          std::isnan(out[1].bearing_deg) && std::isnan(out[1].bearing_rad));
    check("N4 the input records are not changed", links[0].bearing_deg == 90.0);

    std::vector<lsdtt::lsdtt_junction_com_bearing> com(2);
    com[0].outlet_x = -500.0; com[0].outlet_y = 0.0;
    com[0].com_x = 500.0;     com[0].com_y = 0.0;
    com[0].bearing_deg = 90.0;
    com[1].bearing_deg = nan;
    const auto cout_ = lsdtt::junction_com_bearings_to_true_north(com, conv);
    check("N5 COM bearing converted at the outlet / COM midpoint (the centre)",
          std::fabs(cout_[0].bearing_deg - 30.0) < 1e-3);
    check("N6 an undefined COM bearing stays NaN", std::isnan(cout_[1].bearing_deg));

    // Four links all on D8 grid directions; true bearings are 60 deg off, so
    // a fraction taken from them would be 0.
    const std::vector<double> grid = {0.0, 45.0, 90.0, 135.0};
    std::vector<double> truth;
    for (double g : grid) truth.push_back(std::fmod(g - 60.0 + 360.0, 360.0) + 7.0);
    const std::vector<double> w(4, 1.0);
    const std::vector<std::int64_t> blocks = {0, 1, 2, 3};
    lsdtt::channel_record_categories groups;
    groups.key = "all";
    groups.category = {0, 0, 0, 0};
    groups.names = {{0, "all"}};
    lsdtt::channel_bearing_statistics_options opt;
    opt.do_axial = false;
    const auto with_grid = lsdtt::summarise_channel_link_bearings_by_category(
        truth, w, blocks, groups, nullptr, opt, &grid);
    const auto without = lsdtt::summarise_channel_link_bearings_by_category(
        truth, w, blocks, groups, nullptr, opt);
    check("N7 fraction_on_d8_multiples from the grid bearings is 1",
          !with_grid.empty() && with_grid[0].summary.fraction_on_d8_multiples == 1.0);
    check("N8 the same summary without them measures the bearings given (0)",
          !without.empty() && without[0].summary.fraction_on_d8_multiples == 0.0);
    check("N9 the circular mean is the same either way",
          std::fabs(with_grid[0].summary.circular_mean_deg_count
                    - without[0].summary.circular_mean_deg_count) < 1e-12);
}

} // namespace

int main()
{
    lsdtt::gdal_init();
    check_convergence();
    check_wrapping();
    check_planes();
    check_raster();
    check_messages();
    check_channel_records();
    std::cout << (checks - failures) << "/" << checks << " checks passed\n";
    return failures == 0 ? 0 : 1;
}
