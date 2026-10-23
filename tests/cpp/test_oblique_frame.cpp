// test_oblique_frame.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Oblique analysis frames (src/core/lsdtt_oblique_frame). Geodesic positions
// are made with Vincenty's formulae on WGS84, written here independently of
// the module (which goes through PROJ's azimuthal equidistant projection).
//
//   F1  fold_axis_azimuth into (-90, 90].
//   F2  frame from two points, for centres in both hemispheres up to 78 N and
//       azimuths in all quadrants including within 0.1 deg of +-90: centre,
//       folded azimuth, length, snapped extent; the end points map to
//       (+-L/2, 0) and a point d m along the azimuth to (+d, 0), within
//       0.01 % of the distance; swapped end points, and azimuths 180 deg
//       apart, give the same CRS.
//   F3  grid convergence at the centre is 90 - azimuth (60 deg for 30 deg,
//       as measured with pyproj in the spike).
//   F4  datum follows the base CRS: NAD83 frames stay on GRS80/NAD83;
//       base_geographic_wkt_of(UTM 12N NAD83) is NAD83.
//   F5  line: a curved line, and the same line split into shuffled and
//       reversed parts, give the frame of its end points; with several parts
//       the furthest-apart part ends are used.
//   F6  polygon: a rotated rectangle returns its own azimuth, length and
//       width, end buffer on all sides; an irregular polygon matches a
//       brute-force sweep of enclosing rectangles; a square is warned.
//   F7  refusals: width <= 0, end points under 3 cells apart, collinear
//       polygon, cell size <= 0.
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "lsdtt_oblique_frame.hpp"
#include "lsdtt_raster_warp.hpp"

#include <cpl_conv.h>
#include <ogr_spatialref.h>

namespace {

int checks = 0;
int failures = 0;

void check(const std::string& name, bool ok)
{
    ++checks;
    if (!ok) ++failures;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name << '\n';
}

constexpr double pi = 3.14159265358979323846;
constexpr double deg = pi / 180.0;
constexpr double wgs_a = 6378137.0;
constexpr double wgs_f = 1.0 / 298.257223563;

using lsdtt::lsdtt_lonlat;

// Vincenty direct on WGS84: the point s metres from p along azimuth az (deg).
lsdtt_lonlat vincenty_direct(const lsdtt_lonlat& p, double az, double s)
{
    const double a = wgs_a, f = wgs_f, b = a * (1 - f);
    const double alpha1 = az * deg;
    const double sa1 = std::sin(alpha1), ca1 = std::cos(alpha1);
    const double tanU1 = (1 - f) * std::tan(p.latitude * deg);
    const double cosU1 = 1 / std::sqrt(1 + tanU1 * tanU1), sinU1 = tanU1 * cosU1;
    const double sigma1 = std::atan2(tanU1, ca1);
    const double sinAlpha = cosU1 * sa1;
    const double cosSqAlpha = 1 - sinAlpha * sinAlpha;
    const double uSq = cosSqAlpha * (a * a - b * b) / (b * b);
    const double A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)));
    const double B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)));
    double sigma = s / (b * A), sigmaP = 0, cos2SigmaM = 0, sinSigma = 0, cosSigma = 0;
    for (int i = 0; i < 200; ++i)
    {
        cos2SigmaM = std::cos(2 * sigma1 + sigma);
        sinSigma = std::sin(sigma);
        cosSigma = std::cos(sigma);
        const double dSigma = B * sinSigma * (cos2SigmaM + B / 4 * (cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM) -
            B / 6 * cos2SigmaM * (-3 + 4 * sinSigma * sinSigma) * (-3 + 4 * cos2SigmaM * cos2SigmaM)));
        sigmaP = sigma;
        sigma = s / (b * A) + dSigma;
        if (std::abs(sigma - sigmaP) < 1e-14) break;
    }
    cos2SigmaM = std::cos(2 * sigma1 + sigma);
    sinSigma = std::sin(sigma);
    cosSigma = std::cos(sigma);
    const double tmp = sinU1 * sinSigma - cosU1 * cosSigma * ca1;
    const double lat2 = std::atan2(sinU1 * cosSigma + cosU1 * sinSigma * ca1,
                                   (1 - f) * std::sqrt(sinAlpha * sinAlpha + tmp * tmp));
    const double lambda = std::atan2(sinSigma * sa1, cosU1 * cosSigma - sinU1 * sinSigma * ca1);
    const double C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha));
    const double L = lambda - (1 - C) * f * sinAlpha *
        (sigma + C * sinSigma * (cos2SigmaM + C * cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM)));
    return {p.longitude + L / deg, lat2 / deg};
}

// Vincenty inverse on WGS84: distance (m) and initial azimuth (deg).
void vincenty_inverse(const lsdtt_lonlat& p, const lsdtt_lonlat& q, double& s, double& az)
{
    const double a = wgs_a, f = wgs_f, b = a * (1 - f);
    const double L = (q.longitude - p.longitude) * deg;
    const double U1 = std::atan((1 - f) * std::tan(p.latitude * deg));
    const double U2 = std::atan((1 - f) * std::tan(q.latitude * deg));
    const double sinU1 = std::sin(U1), cosU1 = std::cos(U1), sinU2 = std::sin(U2), cosU2 = std::cos(U2);
    double lambda = L, lambdaP = 0, sinSigma = 0, cosSigma = 0, sigma = 0, cosSqAlpha = 0, cos2SigmaM = 0;
    double sinLambda = 0, cosLambda = 0;
    for (int i = 0; i < 200; ++i)
    {
        sinLambda = std::sin(lambda);
        cosLambda = std::cos(lambda);
        sinSigma = std::sqrt((cosU2 * sinLambda) * (cosU2 * sinLambda) +
                             (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda) * (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda));
        if (sinSigma == 0) { s = 0; az = 0; return; }
        cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda;
        sigma = std::atan2(sinSigma, cosSigma);
        const double sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;
        cosSqAlpha = 1 - sinAlpha * sinAlpha;
        cos2SigmaM = cosSqAlpha != 0 ? cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha : 0;
        const double C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha));
        lambdaP = lambda;
        lambda = L + (1 - C) * f * sinAlpha *
            (sigma + C * sinSigma * (cos2SigmaM + C * cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM)));
        if (std::abs(lambda - lambdaP) < 1e-14) break;
    }
    const double uSq = cosSqAlpha * (a * a - b * b) / (b * b);
    const double A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)));
    const double B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)));
    const double dSigma = B * sinSigma * (cos2SigmaM + B / 4 * (cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM) -
        B / 6 * cos2SigmaM * (-3 + 4 * sinSigma * sinSigma) * (-3 + 4 * cos2SigmaM * cos2SigmaM)));
    s = b * A * (sigma - dSigma);
    az = std::atan2(cosU2 * sinLambda, cosU1 * sinU2 - sinU1 * cosU2 * cosLambda) / deg;
}

double distance(const lsdtt_lonlat& p, const lsdtt_lonlat& q)
{
    double s = 0, az = 0;
    vincenty_inverse(p, q, s, az);
    return s;
}

// Smallest difference between two axis azimuths (mod 180), in degrees.
double axis_difference(double a, double b)
{
    double d = std::fmod(std::abs(a - b), 180.0);
    return std::min(d, 180.0 - d);
}

template <typename F>
bool throws_invalid(F&& f)
{
    try { f(); }
    catch (const std::invalid_argument&) { return true; }
    catch (...) { return false; }
    return false;
}

std::string epsg_wkt(int code)
{
    OGRSpatialReference srs;
    srs.importFromEPSG(code);
    char* raw = nullptr;
    srs.exportToWkt(&raw);
    std::string wkt(raw);
    CPLFree(raw);
    return wkt;
}

double snapped(double half, double cell) { return std::ceil(half / cell - 1e-9) * cell; }

void test_fold()
{
    using lsdtt::fold_axis_azimuth;
    check("F1 fold 135 -> -45", std::abs(fold_axis_azimuth(135.0) + 45.0) < 1e-12);
    check("F1 fold -45 stays", std::abs(fold_axis_azimuth(-45.0) + 45.0) < 1e-12);
    check("F1 fold 90 stays", fold_axis_azimuth(90.0) == 90.0);
    check("F1 fold -90 -> 90", fold_axis_azimuth(-90.0) == 90.0);
    check("F1 fold 270 -> 90", fold_axis_azimuth(270.0) == 90.0);
    check("F1 fold 200 -> 20", std::abs(fold_axis_azimuth(200.0) - 20.0) < 1e-12);
    check("F1 fold -100 -> 80", std::abs(fold_axis_azimuth(-100.0) - 80.0) < 1e-12);
    // PROJ omerc misplaces the centre by half a meridian within ~1e-6 deg of
    // alpha = -90 (PROJ 8.2.1), so that edge is taken as +90.
    check("F1 fold -89.99999999 -> 90", fold_axis_azimuth(-89.99999999) == 90.0);
    check("F1 fold 270.00000001 -> 90", fold_axis_azimuth(270.00000001) == 90.0);
    check("F1 fold -89.9998 stays", std::abs(fold_axis_azimuth(-89.9998) + 89.9998) < 1e-12);
}

void test_points()
{
    const std::vector<lsdtt_lonlat> centres = {
        {-110.03, 38.415}, {18.42, -33.92}, {25.0, 65.0}, {166.97, -11.63}, {-45.0, 78.0}};
    const std::vector<double> azimuths = {0.0, 10.0, 30.0, 45.0, 89.9, 90.0, -89.95, 135.0, -45.0};
    const double half = 10000.0;
    const double width = 2000.0;
    const double cell = 30.0;
    int bad_centre = 0, bad_az = 0, bad_len = 0, bad_ends = 0, bad_probe = 0, bad_extent = 0,
        bad_swap = 0, bad_build = 0, bad_conv = 0;
    double worst_rel = 0.0;
    for (const auto& c : centres)
        for (const double az : azimuths)
        {
            const lsdtt_lonlat e1 = vincenty_direct(c, az, half);
            const lsdtt_lonlat e2 = vincenty_direct(c, az + 180.0, half);
            try
            {
                const auto f = lsdtt::frame_from_points(e1, e2, width, -1.0, cell);
                const double folded = lsdtt::fold_axis_azimuth(az);
                if (distance(f.centre, c) > 1e-3) ++bad_centre;
                if (axis_difference(f.azimuth_degrees, folded) > 1e-5 ||
                    f.azimuth_degrees <= -90.0 || f.azimuth_degrees > 90.0) ++bad_az;
                if (std::abs(f.axis_length_metres - 2 * half) > 1e-6 * 2 * half) ++bad_len;
                if (f.x_max != snapped(half + 10 * cell, cell) || f.x_min != -f.x_max ||
                    f.y_max != snapped(0.5 * width, cell) || f.y_min != -f.y_max ||
                    f.n_cols() != static_cast<int>(std::lround(2 * f.x_max / cell)) ||
                    f.n_rows() != static_cast<int>(std::lround(2 * f.y_max / cell)))
                    ++bad_extent;

                lsdtt::lsdtt_oblique_frame_projector p(f);
                // e1 lies along az from the centre: +x when the frame azimuth
                // points the same way (at +-90 either direction is valid).
                const double sign = std::cos((az - f.azimuth_degrees) * deg) > 0.0 ? 1.0 : -1.0;
                double x = 0, y = 0;
                p.to_frame(e1.longitude, e1.latitude, x, y);
                double rel = std::max(std::abs(x - sign * half), std::abs(y)) / half;
                p.to_frame(e2.longitude, e2.latitude, x, y);
                rel = std::max(rel, std::max(std::abs(x + sign * half), std::abs(y)) / half);
                if (rel > 1e-4)
                {
                    ++bad_ends;
                    std::cout << "      ends off at (" << c.longitude << ", " << c.latitude
                              << ") az " << az << ": relative error " << rel << '\n';
                }
                worst_rel = std::max(worst_rel, rel);
                const double d = 0.5 * half;
                const lsdtt_lonlat probe = vincenty_direct(c, f.azimuth_degrees, d);
                p.to_frame(probe.longitude, probe.latitude, x, y);
                rel = std::max(std::abs(x - d), std::abs(y)) / d;
                if (rel > 1e-4) ++bad_probe;
                worst_rel = std::max(worst_rel, rel);
                const double conv = p.grid_convergence_degrees(0.0, 0.0);
                if (std::abs(conv - (90.0 - f.azimuth_degrees)) > 0.01)
                {
                    ++bad_conv;
                    std::cout << "      convergence at (" << c.longitude << ", " << c.latitude
                              << ") az " << az << ": " << conv << '\n';
                }

                const auto g = lsdtt::frame_from_points(e2, e1, width, -1.0, cell);
                if (g.to_wkt() != f.to_wkt() || g.x_max != f.x_max || g.y_max != f.y_max)
                    ++bad_swap;
            }
            catch (const std::exception& e)
            {
                ++bad_build;
                std::cout << "      build failed at (" << c.longitude << ", " << c.latitude << ") az " << az
                          << ": " << e.what() << '\n';
            }
        }
    check("F2 every frame builds (incl. azimuth 0, 90, -89.95, 78 N)", bad_build == 0);
    check("F2 centre is the geodesic midpoint (< 1 mm)", bad_centre == 0);
    check("F2 azimuth folded into (-90, 90] and correct (1e-5 deg)", bad_az == 0);
    check("F2 axis length is the geodesic length (1e-6)", bad_len == 0);
    check("F2 extent symmetric, end buffer 10 cells, snapped outward", bad_extent == 0);
    check("F2 end points map to (+-L/2, 0) within 0.01 %", bad_ends == 0);
    check("F2 a point d m along the azimuth maps to (+d, 0) within 0.01 %", bad_probe == 0);
    check("F2 swapped end points give the same frame", bad_swap == 0);
    check("F3 convergence at the centre is 90 - azimuth (0.01 deg)", bad_conv == 0);
    std::cout << "      worst relative position error " << worst_rel << '\n';

    // 135 and -45 describe the same axis.
    const lsdtt_lonlat c = centres[0];
    const auto f135 = lsdtt::frame_from_points(vincenty_direct(c, 135.0, half),
                                               vincenty_direct(c, -45.0, half), width, 0.0, cell);
    const auto fm45 = lsdtt::frame_from_points(vincenty_direct(c, -45.0, half),
                                               vincenty_direct(c, 135.0, half), width, 0.0, cell);
    check("F2 azimuths 135 and -45 give the same frame", f135.to_wkt() == fm45.to_wkt());
    check("F2 end buffer 0 keeps x_max at L/2 snapped", f135.x_max == snapped(half, cell));

    // The spike's frame: alpha = 30 at the Green River, true north 60 deg
    // clockwise of grid up.
    const auto spike = lsdtt::frame_from_points(vincenty_direct(c, 30.0, half),
                                                vincenty_direct(c, 210.0, half), width, -1.0, cell);
    lsdtt::lsdtt_oblique_frame_projector sp(spike);
    check("F3 spike frame: convergence 60.00 deg at the centre",
          std::abs(sp.grid_convergence_degrees(0.0, 0.0) - 60.0) < 0.005);
    const std::string proj = spike.to_proj();
    check("F2 PROJ string is omerc with gamma 90 on WGS84",
          proj.find("+proj=omerc") != std::string::npos &&
          proj.find("+gamma=90") != std::string::npos &&
          proj.find("+datum=WGS84") != std::string::npos);
    std::cout << "      " << proj << '\n';
}

void test_datum()
{
    const std::string nad83 = epsg_wkt(4269);
    const lsdtt_lonlat c{-110.03, 38.415};
    const auto f = lsdtt::frame_from_points(vincenty_direct(c, 30.0, 5000.0),
                                            vincenty_direct(c, 210.0, 5000.0), 1000.0, -1.0, 30.0,
                                            nad83);
    const std::string proj = f.to_proj();
    check("F4 NAD83 frame keeps NAD83",
          proj.find("+datum=NAD83") != std::string::npos ||
          proj.find("+ellps=GRS80") != std::string::npos);
    check("F4 NAD83 base CRS recorded", lsdtt::crs_equal_semantically(f.base_geographic_wkt, nad83));
    check("F4 base_geographic_wkt_of(UTM 12N NAD83) is NAD83",
          lsdtt::crs_equal_semantically(lsdtt::base_geographic_wkt_of(epsg_wkt(26912)), nad83));
    check("F4 base_geographic_wkt_of(\"\") is WGS84",
          lsdtt::crs_equal_semantically(lsdtt::base_geographic_wkt_of(""), epsg_wkt(4326)));
    std::cout << "      " << proj << '\n';
}

void test_lines()
{
    const lsdtt_lonlat c{-110.03, 38.415};
    const double az = 30.0;
    const lsdtt_lonlat a = vincenty_direct(c, az + 180.0, 8000.0);
    const lsdtt_lonlat b = vincenty_direct(c, az, 8000.0);
    const auto reference = lsdtt::frame_from_points(a, b, 1500.0, -1.0, 30.0);

    // A curved line from a to b: an arc bulging up to 2 km to one side.
    std::vector<lsdtt_lonlat> curve;
    for (int i = 0; i <= 40; ++i)
    {
        const double t = i / 40.0;
        const lsdtt_lonlat on_axis = vincenty_direct(a, az, 16000.0 * t);
        curve.push_back(vincenty_direct(on_axis, az + 90.0, 2000.0 * std::sin(pi * t)));
    }
    curve.front() = a;
    curve.back() = b;
    const auto curved = lsdtt::frame_from_line_parts({curve}, 1500.0, -1.0, 30.0);
    check("F5 curved line uses only its end points",
          curved.to_wkt() == reference.to_wkt() && curved.x_max == reference.x_max &&
          curved.definition == lsdtt::oblique_frame_definition::line);

    // The same line split into parts, shuffled and some reversed.
    std::vector<std::vector<lsdtt_lonlat>> parts = {
        {curve.begin(), curve.begin() + 14},
        {curve.begin() + 14, curve.begin() + 28},
        {curve.begin() + 28, curve.end()}};
    std::reverse(parts[1].begin(), parts[1].end());
    std::swap(parts[0], parts[2]);
    std::reverse(parts[0].begin(), parts[0].end());
    const auto split = lsdtt::frame_from_line_parts(parts, 1500.0, -1.0, 30.0);
    check("F5 shuffled, reversed parts give the same frame", split.to_wkt() == reference.to_wkt());

    // A short separate part further out than b extends the axis to its far end.
    const lsdtt_lonlat far_start = vincenty_direct(c, az, 9000.0);
    const lsdtt_lonlat far_end = vincenty_direct(c, az, 11000.0);
    const auto extended = lsdtt::frame_from_line_parts({curve, {far_start, far_end}}, 1500.0, -1.0, 30.0);
    const auto expected = lsdtt::frame_from_points(a, far_end, 1500.0, -1.0, 30.0);
    check("F5 multipart uses the furthest-apart part ends",
          std::abs(extended.axis_length_metres - 19000.0) < 0.01 &&
          extended.to_wkt() == expected.to_wkt());
}

std::vector<lsdtt_lonlat> rectangle(const lsdtt_lonlat& c, double az, double half_l, double half_w)
{
    std::vector<lsdtt_lonlat> v;
    for (const double sl : {1.0, -1.0})
        for (const double sw : {1.0, -1.0})
        {
            const lsdtt_lonlat along = vincenty_direct(c, sl > 0 ? az : az + 180.0, half_l);
            double s = 0, az_here = 0;
            vincenty_inverse(along, c, s, az_here);
            // Perpendicular at the end of the axis (azimuth back to the centre + 90).
            v.push_back(vincenty_direct(along, az_here + (sw * sl > 0 ? 90.0 : -90.0), half_w));
        }
    std::swap(v[2], v[3]);   // walk the corners in ring order
    return v;
}

void test_polygons()
{
    const double cell = 30.0;
    const lsdtt_lonlat c{-110.03, 38.415};
    int bad = 0;
    for (const double az : {20.0, 75.0, -40.0, 120.0})
    {
        const auto f = lsdtt::frame_from_polygon_vertices(rectangle(c, az, 6000.0, 1500.0), 60.0, cell);
        if (axis_difference(f.azimuth_degrees, az) > 0.02 ||
            std::abs(f.axis_length_metres - 12000.0) > 12.0 ||
            std::abs(f.width_metres - 3000.0) > 3.0 || !f.buffer_all_sides ||
            f.definition != lsdtt::oblique_frame_definition::polygon ||
            f.y_max != snapped(0.5 * f.width_metres + 60.0, cell) ||
            f.x_max != snapped(0.5 * f.axis_length_metres + 60.0, cell) || !f.warnings.empty())
        {
            ++bad;
            std::cout << "      az " << az << " gave " << f.azimuth_degrees << ", L "
                      << f.axis_length_metres << ", W " << f.width_metres << '\n';
        }
    }
    check("F6 rotated rectangles return their azimuth, length and width", bad == 0);

    // Irregular polygon against a brute-force sweep in the same local plane
    // (distances and azimuths from the vertex mean, by Vincenty).
    std::mt19937 rng(20261023);
    std::uniform_real_distribution<double> u_az(0.0, 360.0), u_r(0.0, 1.0);
    std::vector<lsdtt_lonlat> v;
    for (int i = 0; i < 40; ++i)
    {
        const double a = u_az(rng);
        // Elongated cloud: radius larger along 55 deg.
        const double r = 1000.0 + 4000.0 * u_r(rng) * std::abs(std::cos((a - 55.0) * deg));
        v.push_back(vincenty_direct(c, a, r));
    }
    lsdtt_lonlat mean;
    for (const auto& p : v) { mean.longitude += p.longitude / v.size(); mean.latitude += p.latitude / v.size(); }
    std::vector<double> px, py;
    for (const auto& p : v)
    {
        double s = 0, a = 0;
        vincenty_inverse(mean, p, s, a);
        px.push_back(s * std::sin(a * deg));
        py.push_back(s * std::cos(a * deg));
    }
    double best_area = 1e300, best_long_az = 0, best_long = 0, best_short = 0;
    for (int k = 0; k < 18000; ++k)
    {
        const double t = k * 0.01 * deg;   // direction azimuth of the first side
        const double ux = std::sin(t), uy = std::cos(t);
        double lo_u = 1e300, hi_u = -1e300, lo_n = 1e300, hi_n = -1e300;
        for (std::size_t i = 0; i < px.size(); ++i)
        {
            const double a = px[i] * ux + py[i] * uy;
            const double b = -px[i] * uy + py[i] * ux;
            lo_u = std::min(lo_u, a); hi_u = std::max(hi_u, a);
            lo_n = std::min(lo_n, b); hi_n = std::max(hi_n, b);
        }
        const double area = (hi_u - lo_u) * (hi_n - lo_n);
        if (area < best_area)
        {
            best_area = area;
            const bool along = hi_u - lo_u >= hi_n - lo_n;
            best_long_az = along ? k * 0.01 : k * 0.01 + 90.0;
            best_long = std::max(hi_u - lo_u, hi_n - lo_n);
            best_short = std::min(hi_u - lo_u, hi_n - lo_n);
        }
    }
    const auto f = lsdtt::frame_from_polygon_vertices(v, 0.0, cell);
    const double frame_area = f.axis_length_metres * f.width_metres;
    std::cout << "      sweep: az " << lsdtt::fold_axis_azimuth(best_long_az) << ", " << best_long
              << " x " << best_short << "; frame: az " << f.azimuth_degrees << ", "
              << f.axis_length_metres << " x " << f.width_metres << '\n';
    check("F6 irregular polygon: axis matches the brute-force sweep (0.1 deg)",
          axis_difference(f.azimuth_degrees, best_long_az) < 0.1);
    check("F6 irregular polygon: area within 0.1 % of the sweep minimum",
          std::abs(frame_area - best_area) < 1e-3 * best_area);

    const auto square = lsdtt::frame_from_polygon_vertices(rectangle(c, 10.0, 2000.0, 1950.0), 0.0, cell);
    check("F6 near-square polygon is warned as ill-defined", square.warnings.size() == 1);
}

void test_refusals()
{
    const lsdtt_lonlat c{-110.03, 38.415};
    const lsdtt_lonlat a = vincenty_direct(c, 30.0, 5000.0);
    const lsdtt_lonlat b = vincenty_direct(c, 210.0, 5000.0);
    check("F7 width 0 refused", throws_invalid([&] { lsdtt::frame_from_points(a, b, 0.0, -1.0, 30.0); }));
    check("F7 width < 0 refused", throws_invalid([&] { lsdtt::frame_from_points(a, b, -5.0, -1.0, 30.0); }));
    check("F7 cell size 0 refused", throws_invalid([&] { lsdtt::frame_from_points(a, b, 100.0, -1.0, 0.0); }));
    check("F7 end points 2 cells apart refused", throws_invalid([&] {
        lsdtt::frame_from_points(c, vincenty_direct(c, 30.0, 60.0), 100.0, -1.0, 30.0);
    }));
    check("F7 single-vertex line refused",
          throws_invalid([&] { lsdtt::frame_from_line_parts({{a}}, 100.0, -1.0, 30.0); }));
    check("F7 collinear polygon refused", throws_invalid([&] {
        lsdtt::frame_from_polygon_vertices({a, c, b}, 0.0, 30.0);
    }));
    check("F7 two-vertex polygon refused",
          throws_invalid([&] { lsdtt::frame_from_polygon_vertices({a, b}, 0.0, 30.0); }));
}

} // namespace

int main()
{
    test_fold();
    test_points();
    test_datum();
    test_lines();
    test_polygons();
    test_refusals();
    std::cout << "RESULT: " << (failures == 0 ? "PASS" : "FAIL") << " (" << checks - failures
              << "/" << checks << ")\n";
    return failures == 0 ? 0 : 1;
}
