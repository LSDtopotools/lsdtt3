// lsdtt_vertical_datum.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_vertical_datum.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>

#include <ogr_spatialref.h>

#include "../utils/lsdtt_raster_geometry.hpp"
#include "lsdtt_gdal_init.hpp"
#include "lsdtt_proj_grid_check.hpp"

#include <cstdlib>
#include <filesystem>
#include <sstream>

namespace lsdtt
{

namespace
{

const char* compound_crs(lsdtt_vertical_datum d)
{
    return d == lsdtt_vertical_datum::egm96 ? "EPSG:4326+5773" : "EPSG:4326+3855";
}

const char* datum_name(lsdtt_vertical_datum d)
{
    return d == lsdtt_vertical_datum::egm96 ? "EGM96" : d == lsdtt_vertical_datum::egm2008 ? "EGM2008" : "unspecified";
}

[[noreturn]] void missing_grid(lsdtt_vertical_datum from, lsdtt_vertical_datum to)
{
    std::vector<std::string> missing = missing_vertical_datum_grids(from, to);
    if (missing.empty()) missing = vertical_datum_grids(from, to);
    throw std::runtime_error(missing_vertical_datum_grids_message(from, to, missing));
}

} // namespace

std::vector<std::string> vertical_datum_grids(lsdtt_vertical_datum from, lsdtt_vertical_datum to)
{
    std::vector<std::string> g;
    if (from == to || from == lsdtt_vertical_datum::unspecified || to == lsdtt_vertical_datum::unspecified)
        return g;
    for (lsdtt_vertical_datum d : {from, to})
        g.push_back(d == lsdtt_vertical_datum::egm96 ? "us_nga_egm96_15.tif" : "us_nga_egm08_25.tif");
    return g;
}

std::vector<std::string> missing_vertical_datum_grids(lsdtt_vertical_datum from, lsdtt_vertical_datum to)
{
    std::vector<std::string> missing;
    const std::vector<std::string> needed = vertical_datum_grids(from, to);
    if (needed.empty()) return missing;
    const char* net = std::getenv("PROJ_NETWORK");
    if (net && (std::string(net) == "ON" || std::string(net) == "on" || std::string(net) == "1")) return missing;
    // Where lsdtt3 points PROJ: the bundled share/proj. Without one (a source
    // build), ask PROJ itself through the grid check.
    const std::string dir = pinned_proj_data_dir();
    if (!dir.empty())
    {
        for (const auto& g : needed)
        {
            std::error_code ec;
            if (!std::filesystem::is_regular_file(std::filesystem::path(dir) / g, ec)) missing.push_back(g);
        }
        return missing;
    }
    const proj_grid_check_result r =
        check_proj_transformation_grids(compound_crs(from), compound_crs(to));
    for (const auto& g : r.missing_grids)
        for (const auto& n : needed)
            if (g == n) missing.push_back(g);
    return missing;
}

std::string missing_vertical_datum_grids_message(lsdtt_vertical_datum from, lsdtt_vertical_datum to,
                                                 const std::vector<std::string>& missing)
{
    std::ostringstream s;
    s << "vertical datum conversion " << datum_name(from) << " -> " << datum_name(to)
      << " needs the PROJ geoid grid";
    if (missing.size() > 1) s << "s";
    for (std::size_t k = 0; k < missing.size(); ++k) s << (k ? ", " : " ") << missing[k];
    s << ", which " << (missing.size() > 1 ? "are" : "is") << " not installed.";
    for (const auto& g : missing)
        if (g == "us_nga_egm08_25.tif")
            s << " The EGM2008 grid (77 MB) is not bundled with lsdtt3 because of its size.";
    s << " Install with: python3 fetch_proj_grids.py fetch";
    for (const auto& g : missing) s << ' ' << g;
    s << " (tools/fetch_proj_grids.py in the lsdtt3 repository; see "
         "https://lsdtopotools.github.io/lsdtt3/what_is_lsdtt3/note_on_grids.html), or let PROJ read "
         "them from its CDN with the environment variable PROJ_NETWORK=ON.";
    return s.str();
}

lsdtt_vertical_datum parse_vertical_datum(const std::string& name)
{
    if (name == "unspecified") return lsdtt_vertical_datum::unspecified;
    if (name == "egm96") return lsdtt_vertical_datum::egm96;
    if (name == "egm2008") return lsdtt_vertical_datum::egm2008;
    throw std::invalid_argument("vertical datum must be unspecified, egm96 or egm2008 (got '" +
                                name + "')");
}

std::vector<float> vertical_datum_shift(
    const lsdtt_raster_meta& grid,
    lsdtt_vertical_datum     from,
    lsdtt_vertical_datum     to)
{
    if (from == lsdtt_vertical_datum::unspecified || to == lsdtt_vertical_datum::unspecified)
        throw std::invalid_argument("vertical_datum_shift: both datums must be specified");
    if (grid.crs_wkt.empty())
        throw std::invalid_argument("vertical_datum_shift: the grid has no CRS");

    const std::size_t n = static_cast<std::size_t>(grid.n_rows) * static_cast<std::size_t>(grid.n_cols);
    if (from == to) return std::vector<float>(n, 0.0f);

    gdal_init();
    OGRSpatialReference grid_srs, lonlat, source, target;
    if (grid_srs.importFromWkt(grid.crs_wkt.c_str()) != OGRERR_NONE)
        throw std::invalid_argument("vertical_datum_shift: the grid's CRS cannot be read");
    lonlat.importFromEPSG(4326);
    if (source.SetFromUserInput(compound_crs(from)) != OGRERR_NONE ||
        target.SetFromUserInput(compound_crs(to)) != OGRERR_NONE)
        missing_grid(from, to);
    for (OGRSpatialReference* s : {&grid_srs, &lonlat, &source, &target})
        s->SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    std::unique_ptr<OGRCoordinateTransformation> to_lonlat(
        OGRCreateCoordinateTransformation(&grid_srs, &lonlat));
    OGRCoordinateTransformationOptions options;
    options.SetBallparkAllowed(false);   // a ballpark transform would return a zero shift
    std::unique_ptr<OGRCoordinateTransformation> geoid(
        OGRCreateCoordinateTransformation(&source, &target, options));
    if (!to_lonlat || !geoid) missing_grid(from, to);

    const double cell_metres = grid_srs.IsGeographic() ? grid.cell_size * 111320.0 : grid.cell_size;
    const int step = std::max(1, static_cast<int>(std::floor(1000.0 / cell_metres)));
    const int lattice_rows = (grid.n_rows - 1) / step + 2;
    const int lattice_cols = (grid.n_cols - 1) / step + 2;
    std::vector<double> lattice(static_cast<std::size_t>(lattice_rows) * lattice_cols);
    for (int lr = 0; lr < lattice_rows; ++lr)
    {
        std::vector<double> x(lattice_cols), y(lattice_cols), z(lattice_cols, 0.0);
        std::vector<int> ok(lattice_cols, 0);
        for (int lc = 0; lc < lattice_cols; ++lc)
        {
            // Beyond the last row or column is fine: this is a smooth field.
            const lsdtt_xy p = cell_center_from_row_col(grid, lr * step, lc * step);
            x[lc] = p.x;
            y[lc] = p.y;
        }
        if (!to_lonlat->Transform(lattice_cols, x.data(), y.data(), nullptr, ok.data()) ||
            !std::all_of(ok.begin(), ok.end(), [](int v) { return v != 0; }))
            throw std::runtime_error("vertical_datum_shift: cannot transform the grid to longitude and latitude");
        std::fill(ok.begin(), ok.end(), 0);
        if (!geoid->Transform(lattice_cols, x.data(), y.data(), z.data(), ok.data()) ||
            !std::all_of(ok.begin(), ok.end(), [](int v) { return v != 0; }))
            missing_grid(from, to);
        for (int lc = 0; lc < lattice_cols; ++lc)
        {
            if (!std::isfinite(z[lc])) missing_grid(from, to);
            lattice[static_cast<std::size_t>(lr) * lattice_cols + lc] = z[lc];
        }
    }

    std::vector<float> out(n);
    for (int r = 0; r < grid.n_rows; ++r)
    {
        const int r0 = r / step;
        const double fr = static_cast<double>(r - r0 * step) / step;
        for (int c = 0; c < grid.n_cols; ++c)
        {
            const int c0 = c / step;
            const double fc = static_cast<double>(c - c0 * step) / step;
            const auto at = [&](int a, int b) {
                return lattice[static_cast<std::size_t>(a) * lattice_cols + b];
            };
            const double v = (1 - fr) * ((1 - fc) * at(r0, c0) + fc * at(r0, c0 + 1)) +
                             fr * ((1 - fc) * at(r0 + 1, c0) + fc * at(r0 + 1, c0 + 1));
            out[static_cast<std::size_t>(r) * grid.n_cols + c] = static_cast<float>(v);
        }
    }
    return out;
}

} // namespace lsdtt
