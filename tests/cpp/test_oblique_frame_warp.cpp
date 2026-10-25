// test_oblique_frame_warp.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Warping into an oblique frame and rotating back (warp_raster_to_frame,
// rotate_back_raster, rotation_frame_from_raster). Fixtures are GeoTIFFs
// written by this test into the work directory: a smooth 10 m DEM in UTM 12N
// whose origin is NOT on the 10 m lattice (580003.7, 4252601.3).
//
//   W1  the DEM warped into a frame loads in the frame CRS on the frame grid,
//       passes check_projected_metres, fills and routes.
//   W2  rotate-back with the source origin as target lands on the source DEM's
//       own pixels (origin offset a whole number of cells); values agree with
//       the source (measured, reported).
//   W3  rotate-back to the (0, 0) lattice (the fetch -tap case) lands on whole
//       multiples of the cell size.
//   W4  rotate_back_epsg naming another CRS: that CRS, (0, 0) lattice, warned.
//   W5  near keeps integer IDs intact (no new values).
//   W6  a raster at another cell size keeps its own cell size.
//   W7  the missing-frame-file fallback rebuilds the frame from the raster,
//       targets UTM 12N on the (0, 0) lattice and warns.
//   R*  refusals: a raster not in the frame CRS (naming it), an unknown
//       method, the fallback on a raster that is not in a frame.
//
// Usage: lsdtt-test-oblique-frame-warp <work_dir>
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include <cpl_conv.h>
#include <gdal_priv.h>
#include <ogr_spatialref.h>

#include "lsdtt_data_loader.hpp"
#include "lsdtt_gdal_init.hpp"
#include "lsdtt_oblique_frame.hpp"
#include "lsdtt_raster_warp.hpp"
#include "lsdtt_rotation_frame_file.hpp"
#include "../../src/cpu/lsdtt_flowinfo_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"

namespace {

int checks = 0;
int failures = 0;

void check(const std::string& name, bool ok)
{
    ++checks;
    if (!ok) ++failures;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name << '\n';
}

bool contains(const std::string& s, const std::string& part)
{
    return s.find(part) != std::string::npos;
}

std::string work;

std::string wkt_of_epsg(int code)
{
    OGRSpatialReference srs;
    srs.importFromEPSG(code);
    char* raw = nullptr;
    srs.exportToWkt(&raw);
    std::string wkt(raw ? raw : "");
    CPLFree(raw);
    return wkt;
}

template <typename T>
void write_tif(const std::string& path, const std::string& wkt, double x_min, double y_max,
               double cell, int cols, int rows, const std::vector<T>& data, GDALDataType type,
               double nodata)
{
    GDALDriver* drv = GetGDALDriverManager()->GetDriverByName("GTiff");
    GDALDataset* ds = drv->Create(path.c_str(), cols, rows, 1, type, nullptr);
    double gt[6] = {x_min, cell, 0.0, y_max, 0.0, -cell};
    ds->SetGeoTransform(gt);
    ds->SetProjection(wkt.c_str());
    GDALRasterBand* b = ds->GetRasterBand(1);
    b->SetNoDataValue(nodata);
    if (b->RasterIO(GF_Write, 0, 0, cols, rows, const_cast<T*>(data.data()), cols, rows, type, 0,
                    0) != CE_None)
        std::cerr << "write failed: " << path << '\n';
    GDALClose(ds);
}

lsdtt::lsdtt_raster_meta meta_of(const std::string& path)
{
    lsdtt::lsdtt_data_loader loader;
    loader.open_raster(path);
    return loader.get_meta();
}

lsdtt::lsdtt_raster_tile tile_of(const std::string& path)
{
    lsdtt::lsdtt_data_loader loader;
    loader.open_raster(path);
    return loader.read_full();
}

// Offset of a in cells from b, and its distance from a whole number.
double off_lattice(double a, double b, double cell)
{
    const double k = (a - b) / cell;
    return std::abs(k - std::round(k));
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "usage: lsdtt-test-oblique-frame-warp <work_dir>\n";
        return 2;
    }
    work = argv[1];
    std::filesystem::create_directories(work);
    lsdtt::gdal_init();

    // --- Fixture: smooth UTM 12N DEM on an unsnapped origin --------------------
    const std::string utm = wkt_of_epsg(32612);
    const double sx0 = 580003.7;
    const double sy1 = 4252601.3;
    const double cell = 10.0;
    const int scols = 300;
    const int srows = 240;
    const double sy0 = sy1 - srows * cell;
    std::vector<float> dem(static_cast<std::size_t>(scols) * srows);
    for (int r = 0; r < srows; ++r)
        for (int c = 0; c < scols; ++c)
        {
            const double dx = (c + 0.5) * cell;
            const double dy = (srows - r - 0.5) * cell;
            dem[static_cast<std::size_t>(r) * scols + c] = static_cast<float>(
                500.0 + 0.05 * dx + 0.03 * dy + 5.0 * std::sin(dx / 300.0) * std::cos(dy / 250.0));
        }
    const std::string src = work + "/source_utm.tif";
    write_tif(src, utm, sx0, sy1, cell, scols, srows, dem, GDT_Float32, -9999.0);

    // Frame across the middle of the DEM, azimuth about 35 deg.
    lsdtt::lsdtt_lonlat a, b;
    {
        OGRSpatialReference s, g;
        s.importFromEPSG(32612);
        g.SetWellKnownGeogCS("WGS84");
        s.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        g.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        OGRCoordinateTransformation* ct = OGRCreateCoordinateTransformation(&s, &g);
        const double cx = sx0 + 0.5 * scols * cell;
        const double cy = sy0 + 0.5 * srows * cell;
        double x1 = cx - 600.0, y1 = cy - 850.0, x2 = cx + 600.0, y2 = cy + 850.0;
        ct->Transform(1, &x1, &y1);
        ct->Transform(1, &x2, &y2);
        OGRCoordinateTransformation::DestroyCT(ct);
        a = {x1, y1};
        b = {x2, y2};
    }
    const lsdtt::lsdtt_oblique_frame frame =
        lsdtt::frame_from_points(a, b, 800.0, -1.0, cell);

    // --- W1: warp into the frame ---------------------------------------------
    const std::string rot = work + "/dem_oblique.tif";
    lsdtt::lsdtt_warp_options wopts;
    wopts.resampling_method = "bilinear";
    lsdtt::warp_raster_to_frame(src, frame, rot, wopts);
    {
        const lsdtt::lsdtt_raster_meta m = meta_of(rot);
        check("W1 frame raster in the frame CRS", lsdtt::crs_equal_semantically(m.crs_wkt,
                                                                                frame.to_wkt()));
        check("W1 frame grid (origin, cell, size)",
              std::abs(m.x_min - frame.x_min) < 1e-6 && std::abs(m.y_min - frame.y_min) < 1e-6 &&
              std::abs(m.cell_size - cell) < 1e-9 && m.n_cols == frame.n_cols() &&
              m.n_rows == frame.n_rows());
        bool metres_ok = true;
        try { m.check_projected_metres(); } catch (...) { metres_ok = false; }
        check("W1 passes check_projected_metres", metres_ok);

        lsdtt::lsdtt_data_loader loader;
        loader.open_raster(rot);
        const lsdtt::lsdtt_raster_cpu raw(loader);
        std::size_t valid = 0;
        for (auto v : raw.get_mask()) valid += v ? 1 : 0;
        const double frac = static_cast<double>(valid) / (m.n_cols * m.n_rows);
        const lsdtt::lsdtt_raster_cpu filled = lsdtt::fill_pits(raw);
        lsdtt::lsdtt_flowinfo_options fi;
        const lsdtt::lsdtt_flowinfo flow(filled.get_data_array(), filled.get_mask_array(),
                                         filled.get_meta(), fi);
        double max_area = 0.0;
        for (std::int64_t i = 0; i < static_cast<std::int64_t>(m.n_cols) * m.n_rows; ++i)
            max_area = std::max(max_area, flow.drainage_area(i));
        std::cout << "      frame raster " << m.n_cols << " x " << m.n_rows << ", valid "
                  << frac << ", max drainage area " << max_area << " m2\n";
        check("W1 fills and routes (valid cells, drainage area > 100 cells)",
              frac > 0.99 && max_area > 100 * cell * cell);
    }

    // Frame file with the source DEM's own lattice as target.
    lsdtt::lsdtt_rotation_target target;
    target.epsg = 32612;
    target.cell_size = cell;
    target.origin_x = sx0;
    target.origin_y = sy0;
    target.source_dem_fname = src;
    lsdtt::lsdtt_rotation_provenance prov;
    prov.created_by = "test";
    prov.rotated_dem_fname = rot;
    const std::string frame_path = work + "/dem_oblique_rotation_frame.txt";
    lsdtt::write_rotation_frame(frame_path, frame, target, prov);
    const lsdtt::lsdtt_rotation_frame_file ff = lsdtt::read_rotation_frame(frame_path);
    check("frame file reads", ff.ok);

    // --- W2: back onto the source pixels ---------------------------------------
    {
        const std::string out = work + "/dem_rotated_back.tif";
        const auto res = lsdtt::rotate_back_raster(rot, ff, "bilinear", 0, out);
        const lsdtt::lsdtt_raster_meta m = meta_of(out);
        const double ox = off_lattice(m.x_min, sx0, cell);
        const double oy = off_lattice(m.y_min, sy0, cell);
        std::cout << "      rotated back " << m.n_cols << " x " << m.n_rows << " at ("
                  << std::setprecision(12) << m.x_min << ", " << m.y_min
                  << "), off-lattice " << ox << ", " << oy << " cells\n";
        check("W2 on the source DEM's pixels",
              ox < 1e-6 && oy < 1e-6 && std::abs(m.cell_size - cell) < 1e-9 &&
              lsdtt::crs_equal_semantically(m.crs_wkt, utm) && res.warnings.empty());

        // Compare with the source where both are valid.
        const lsdtt::lsdtt_raster_tile t = tile_of(out);
        const int dc = static_cast<int>(std::lround((m.x_min - sx0) / cell));
        const int dr = static_cast<int>(std::lround((sy1 - m.y_max()) / cell));
        std::vector<double> diffs;
        for (int r = 0; r < m.n_rows; ++r)
            for (int c = 0; c < m.n_cols; ++c)
            {
                const int sr = r + dr, sc = c + dc;
                if (sr < 0 || sr >= srows || sc < 0 || sc >= scols) continue;
                const float v = t.data[static_cast<std::size_t>(r) * m.n_cols + c];
                if (v == -9999.0f || !std::isfinite(v)) continue;
                diffs.push_back(std::abs(v - dem[static_cast<std::size_t>(sr) * scols + sc]));
            }
        std::sort(diffs.begin(), diffs.end());
        const double p50 = diffs.empty() ? -1 : diffs[diffs.size() / 2];
        const double pmax = diffs.empty() ? -1 : diffs.back();
        std::cout << "      " << diffs.size() << " overlapping cells, |dz| median " << p50
                  << " m, max " << pmax << " m (bilinear both ways)\n";
        check("W2 values overlap the source (cells compared > 0)", diffs.size() > 1000);
    }

    // --- W3: the (0, 0) lattice --------------------------------------------------
    {
        lsdtt::lsdtt_rotation_frame_file tap = ff;
        tap.target.origin_x = 0.0;
        tap.target.origin_y = 0.0;
        const std::string out = work + "/dem_rotated_back_tap.tif";
        lsdtt::rotate_back_raster(rot, tap, "bilinear", 0, out);
        const lsdtt::lsdtt_raster_meta m = meta_of(out);
        check("W3 on whole multiples of the cell size",
              off_lattice(m.x_min, 0.0, cell) < 1e-6 && off_lattice(m.y_min, 0.0, cell) < 1e-6);
    }

    // --- W4: epsg override ---------------------------------------------------------
    {
        const std::string out = work + "/dem_rotated_back_nad83.tif";
        const auto res = lsdtt::rotate_back_raster(rot, ff, "bilinear", 26912, out);
        const lsdtt::lsdtt_raster_meta m = meta_of(out);
        check("W4 override: that CRS, (0, 0) lattice, warned",
              lsdtt::crs_equal_semantically(m.crs_wkt, wkt_of_epsg(26912)) &&
              off_lattice(m.x_min, 0.0, cell) < 1e-6 && off_lattice(m.y_min, 0.0, cell) < 1e-6 &&
              !res.warnings.empty());
    }

    // --- W5: near keeps IDs ----------------------------------------------------
    const int fcols = frame.n_cols();
    const int frows = frame.n_rows();
    {
        std::vector<std::int32_t> ids(static_cast<std::size_t>(fcols) * frows);
        std::set<std::int32_t> id_set;
        for (int r = 0; r < frows; ++r)
            for (int c = 0; c < fcols; ++c)
            {
                const std::int32_t id = 1 + (r / 13) * 1000 + (c / 17) * 7;
                ids[static_cast<std::size_t>(r) * fcols + c] = id;
                id_set.insert(id);
            }
        const std::string id_path = work + "/ids_oblique.tif";
        write_tif(id_path, frame.to_wkt(), frame.x_min, frame.y_max, cell, fcols, frows, ids,
                  GDT_Int32, -9999.0);
        const std::string out = work + "/ids_rotated_back.tif";
        lsdtt::rotate_back_raster(id_path, ff, "near", 0, out);
        const lsdtt::lsdtt_raster_tile t = tile_of(out);
        std::size_t foreign = 0;
        std::set<std::int32_t> seen;
        for (std::size_t i = 0; i < t.data.size(); ++i)
        {
            if (!t.mask.empty() && !t.mask[i]) continue;
            const float v = t.data[i];
            if (v == -9999.0f) continue;
            const auto id = static_cast<std::int32_t>(std::lround(v));
            if (static_cast<float>(id) != v || !id_set.count(id)) ++foreign;
            else seen.insert(id);
        }
        std::cout << "      " << seen.size() << " of " << id_set.size()
                  << " IDs present after near, " << foreign << " foreign values\n";
        check("W5 near: no new values, IDs kept", foreign == 0 && seen.size() == id_set.size());
    }

    // --- W6: own cell size -------------------------------------------------------
    {
        const double c2 = 20.0;
        const int c2cols = static_cast<int>(std::lround((frame.x_max - frame.x_min) / c2));
        const int c2rows = static_cast<int>(std::lround((frame.y_max - frame.y_min) / c2));
        std::vector<float> v(static_cast<std::size_t>(c2cols) * c2rows, 1.0f);
        const std::string p20 = work + "/coarse_oblique.tif";
        write_tif(p20, frame.to_wkt(), frame.x_min, frame.y_max, c2, c2cols, c2rows, v,
                  GDT_Float32, -9999.0);
        const std::string out = work + "/coarse_rotated_back.tif";
        const auto res = lsdtt::rotate_back_raster(p20, ff, "bilinear", 0, out);
        const lsdtt::lsdtt_raster_meta m = meta_of(out);
        check("W6 keeps its own 20 m cell, anchored at the target origin",
              std::abs(m.cell_size - c2) < 1e-9 && res.cell_size == c2 &&
              off_lattice(m.x_min, sx0, c2) < 1e-6 && off_lattice(m.y_min, sy0, c2) < 1e-6);
    }

    // --- W7: missing frame file ----------------------------------------------------
    {
        const lsdtt::lsdtt_rotation_frame_file fb = lsdtt::rotation_frame_from_raster(rot);
        bool warned = false;
        for (const auto& w : fb.warnings) warned = warned || contains(w, "cannot be recovered");
        check("W7 fallback: frame rebuilt, UTM 12N target, warned",
              fb.ok && fb.target.epsg == 32612 && warned &&
              fb.frame.definition == lsdtt::oblique_frame_definition::crs &&
              lsdtt::crs_equal_semantically(fb.frame.to_wkt(), frame.to_wkt()));
        if (!fb.ok) std::cout << "      " << fb.message << '\n';
        const std::string out = work + "/dem_rotated_back_fallback.tif";
        lsdtt::rotate_back_raster(rot, fb, "bilinear", 0, out);
        const lsdtt::lsdtt_raster_meta m = meta_of(out);
        check("W7 fallback rotate-back on the (0, 0) lattice",
              off_lattice(m.x_min, 0.0, cell) < 1e-6 && off_lattice(m.y_min, 0.0, cell) < 1e-6);

        const lsdtt::lsdtt_rotation_frame_file not_frame = lsdtt::rotation_frame_from_raster(src);
        check("R fallback refuses a raster not in a frame",
              !not_frame.ok && contains(not_frame.message, "not an oblique Mercator"));
    }

    // --- Refusals ----------------------------------------------------------------
    {
        bool refused = false;
        std::string what;
        try { lsdtt::rotate_back_raster(src, ff, "bilinear", 0, work + "/never.tif"); }
        catch (const std::invalid_argument& e) { refused = true; what = e.what(); }
        check("R raster not in the frame CRS refused, naming it",
              refused && contains(what, "UTM zone 12N") && contains(what, "not the frame CRS"));

        refused = false;
        try { lsdtt::rotate_back_raster(rot, ff, "nearest_neighbour", 0, work + "/never.tif"); }
        catch (const std::invalid_argument& e) { refused = contains(e.what(), "unknown"); }
        check("R unknown resampling method refused", refused);
    }

    std::cout << "\n" << checks - failures << " / " << checks << " checks passed\n";
    return failures == 0 ? 0 : 1;
}
