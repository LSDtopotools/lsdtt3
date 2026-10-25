// lsdtt_true_north.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// See lsdtt_true_north.hpp.

#include "lsdtt_true_north.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <stdexcept>

#include <ogr_spatialref.h>

#include "lsdtt_gdal_init.hpp"

namespace lsdtt {

namespace {

constexpr double pi = 3.14159265358979323846;
constexpr double deg = pi / 180.0;

// Convergence below this (degrees) is not reported at start-up.
constexpr double grid_north_message_threshold_degrees = 1.0;

struct ct_deleter
{
    void operator()(OGRCoordinateTransformation* ct) const
    {
        OGRCoordinateTransformation::DestroyCT(ct);
    }
};
using ct_ptr = std::unique_ptr<OGRCoordinateTransformation, ct_deleter>;

void transform(OGRCoordinateTransformation& ct, double& x, double& y)
{
    if (!ct.Transform(1, &x, &y))
        throw std::runtime_error("lsdtt_grid_convergence: coordinate transform failed");
}

double wrap(double v, double period)
{
    double w = std::fmod(v, period);
    if (w < 0.0) w += period;
    // fmod of a value just below 0 can round up to the period itself.
    if (w >= period) w -= period;
    return w;
}

} // namespace

struct lsdtt_grid_convergence::impl
{
    bool geographic = false;
    ct_ptr to_geog;
    ct_ptr to_grid;
};

lsdtt_grid_convergence::lsdtt_grid_convergence(const std::string& crs_wkt)
    : impl_(std::make_unique<impl>())
{
    if (crs_wkt.empty())
        throw std::invalid_argument(
            "lsdtt_grid_convergence: the raster has no CRS, so grid north cannot be "
            "related to true north");
    gdal_init();
    OGRSpatialReference grid;
    if (grid.importFromWkt(crs_wkt.c_str()) != OGRERR_NONE)
        throw std::runtime_error("lsdtt_grid_convergence: could not read the raster CRS");
    grid.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    if (grid.IsGeographic())
    {
        impl_->geographic = true;
        return;
    }
    OGRSpatialReference geog;
    if (geog.CopyGeogCSFrom(&grid) != OGRERR_NONE)
        throw std::runtime_error("lsdtt_grid_convergence: the raster CRS has no geographic CRS");
    geog.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    impl_->to_geog.reset(OGRCreateCoordinateTransformation(&grid, &geog));
    impl_->to_grid.reset(OGRCreateCoordinateTransformation(&geog, &grid));
    if (!impl_->to_geog || !impl_->to_grid)
        throw std::runtime_error("lsdtt_grid_convergence: could not create a CRS transform");
}

lsdtt_grid_convergence::~lsdtt_grid_convergence() = default;

bool lsdtt_grid_convergence::is_geographic() const
{
    return impl_->geographic;
}

double lsdtt_grid_convergence::degrees_at(double x, double y) const
{
    if (impl_->geographic) return 0.0;
    double lon = x;
    double lat = y;
    transform(*impl_->to_geog, lon, lat);
    // Probe a short step north (south within a degree of the north pole; the
    // sign of the step sets the direction measured).
    const double step = lat > 89.0 ? -1e-4 : 1e-4;
    double x2 = lon;
    double y2 = lat + step;
    transform(*impl_->to_grid, x2, y2);
    double dx = x2 - x;
    double dy = y2 - y;
    if (step < 0.0) { dx = -dx; dy = -dy; }
    return std::atan2(dx, dy) / deg;
}

double lsdtt_grid_convergence::to_true_north(
    double angle, lsdtt_angle_kind kind, double x, double y) const
{
    if (!std::isfinite(angle)) return angle;
    return grid_angle_to_true_north(angle, kind, degrees_at(x, y));
}

double grid_angle_to_true_north(double angle, lsdtt_angle_kind kind,
                                double convergence_degrees)
{
    if (!std::isfinite(angle)) return angle;
    switch (kind)
    {
    case lsdtt_angle_kind::bearing_degrees:
        return wrap(angle - convergence_degrees, 360.0);
    case lsdtt_angle_kind::bearing_radians:
        return wrap(angle - convergence_degrees * deg, 2.0 * pi);
    case lsdtt_angle_kind::axial_degrees:
        return wrap(angle - convergence_degrees, 180.0);
    }
    return angle;
}

void convert_angle_raster_to_true_north(
    std::vector<float>&           data,
    const lsdtt_raster_meta&      meta,
    lsdtt_angle_kind              kind,
    const lsdtt_grid_convergence& convergence,
    int                           lattice_step)
{
    const int nr = meta.n_rows;
    const int nc = meta.n_cols;
    if (!meta.is_valid() || data.size() != static_cast<std::size_t>(meta.n_pixels()))
        throw std::invalid_argument(
            "convert_angle_raster_to_true_north: data do not match the raster metadata");
    if (lattice_step < 1)
        throw std::invalid_argument("convert_angle_raster_to_true_north: lattice_step must be >= 1");

    // Lattice rows and columns: every lattice_step cells, plus the last one.
    auto lattice = [lattice_step](int n) {
        std::vector<int> v;
        for (int i = 0; i < n; i += lattice_step) v.push_back(i);
        if (v.back() != n - 1) v.push_back(n - 1);
        return v;
    };
    const std::vector<int> lr = lattice(nr);
    const std::vector<int> lc = lattice(nc);
    const double cs = meta.cell_size;
    const double ph = meta.pixel_height();

    std::vector<double> c(lr.size() * lc.size());
    for (std::size_t i = 0; i < lr.size(); ++i)
        for (std::size_t j = 0; j < lc.size(); ++j)
        {
            const double x = meta.x_min + (lc[j] + 0.5) * cs;
            const double y = meta.y_min + (nr - lr[i] - 0.5) * ph;
            c[i * lc.size() + j] = convergence.degrees_at(x, y);
        }

    const int nlc = static_cast<int>(lc.size());
    for (int r = 0; r < nr; ++r)
    {
        const std::size_t i0 = static_cast<std::size_t>(
            std::upper_bound(lr.begin(), lr.end(), r) - lr.begin()) - 1;
        const std::size_t i1 = std::min(i0 + 1, lr.size() - 1);
        const double fr = i1 == i0 ? 0.0
            : static_cast<double>(r - lr[i0]) / static_cast<double>(lr[i1] - lr[i0]);
        std::size_t j0 = 0;
        for (int col = 0; col < nc; ++col)
        {
            while (j0 + 1 < lc.size() && lc[j0 + 1] <= col) ++j0;
            float& v = data[static_cast<std::size_t>(r) * nc + col];
            if (meta.is_nodata(v) || !std::isfinite(v)) continue;
            const std::size_t j1 = std::min<std::size_t>(j0 + 1, lc.size() - 1);
            const double fc = j1 == j0 ? 0.0
                : static_cast<double>(col - lc[j0]) / static_cast<double>(lc[j1] - lc[j0]);
            const double c00 = c[i0 * nlc + j0];
            const double c01 = c[i0 * nlc + j1];
            const double c10 = c[i1 * nlc + j0];
            const double c11 = c[i1 * nlc + j1];
            const double cv = (1.0 - fr) * ((1.0 - fc) * c00 + fc * c01)
                            + fr * ((1.0 - fc) * c10 + fc * c11);
            v = static_cast<float>(grid_angle_to_true_north(v, kind, cv));
        }
    }
}

std::string grid_north_message(const lsdtt_raster_meta& meta,
                               bool bearings_relative_to_true_north)
{
    if (meta.crs_wkt.empty() || !meta.is_valid()) return "";
    double c = 0.0;
    try
    {
        lsdtt_grid_convergence conv(meta.crs_wkt);
        if (conv.is_geographic()) return "";
        c = conv.degrees_at(meta.x_min + 0.5 * (meta.x_max() - meta.x_min),
                            meta.y_min + 0.5 * (meta.y_max() - meta.y_min));
    }
    catch (const std::exception&)
    {
        return "";
    }
    if (!bearings_relative_to_true_north
        && std::fabs(c) <= grid_north_message_threshold_degrees)
        return "";
    std::ostringstream os;
    os << std::fixed << std::setprecision(1)
       << "Grid north is " << std::fabs(c) << " deg "
       << (c >= 0.0 ? "anticlockwise" : "clockwise")
       << " of true north at the grid centre; ";
    if (bearings_relative_to_true_north)
        os << "bearings are converted to true north at each output's location "
              "(hillshade illumination stays relative to grid north).";
    else
        os << "bearings are relative to grid north "
              "(bearings_relative_to_true_north is false).";
    return os.str();
}

std::unique_ptr<lsdtt_grid_convergence> true_north_converter(
    const lsdtt_raster_meta& meta,
    bool bearings_relative_to_true_north)
{
    if (!bearings_relative_to_true_north) return nullptr;
    return std::make_unique<lsdtt_grid_convergence>(meta.crs_wkt);
}

std::unique_ptr<lsdtt_grid_convergence> start_true_north_bearings(
    const lsdtt_raster_meta& meta,
    bool bearings_relative_to_true_north,
    std::ostream& out)
{
    const std::string msg = grid_north_message(meta, bearings_relative_to_true_north);
    if (!msg.empty()) out << msg << '\n';
    return true_north_converter(meta, bearings_relative_to_true_north);
}

} // namespace lsdtt
