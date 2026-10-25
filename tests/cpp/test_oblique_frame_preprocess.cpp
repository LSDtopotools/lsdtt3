// test_oblique_frame_preprocess.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Rotating an existing DEM into an oblique frame, as lsdtt-raster-preprocessing
// does (frame_cell_size_for_raster, rotation_target_for_raster,
// frame_from_valid_data, rotate_raster_to_frame), and rotating its products
// back. Fixtures are written by this test into the work directory: a smooth
// 10 m DEM in UTM 12N with an origin NOT on the 10 m lattice (580003.7,
// 4252601.3), the same DEM in NAD83 / UTM 12N, a DEM whose valid cells are a
// rotated strip, a 1 arc-second geographic DEM, and GeoJSON feature files.
//
//   S*  cell size and rotate-back target from the source raster.
//   L*  line feature: the frame equals frame_from_points on the line's end
//       points; the rotated DEM and its frame file; DEM and ID products
//       rotated back land on the source pixels.
//   G*  polygon feature (a rotated rectangle): azimuth, centre and size.
//   F*  valid-data footprint: azimuth and length of the strip, frame_source
//       data_footprint read back without a rebuild warning.
//   D*  NAD83 source: the frame is on NAD83, the target is EPSG:26912 and the
//       round trip matches the WGS84 one.
//   R*  refusals: no valid cells, no CRS.
//
// Usage: lsdtt-test-oblique-frame-preprocess <work_dir>
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
#include <ogrsf_frmts.h>

#include "lsdtt_data_loader.hpp"
#include "lsdtt_gdal_init.hpp"
#include "lsdtt_oblique_frame.hpp"
#include "lsdtt_raster_warp.hpp"
#include "lsdtt_rotation_frame_file.hpp"

namespace {

int checks = 0;
int failures = 0;

void check(const std::string& name, bool ok)
{
    ++checks;
    if (!ok) ++failures;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name << '\n';
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
    if (!wkt.empty()) ds->SetProjection(wkt.c_str());
    GDALRasterBand* b = ds->GetRasterBand(1);
    b->SetNoDataValue(nodata);
    if (b->RasterIO(GF_Write, 0, 0, cols, rows, const_cast<T*>(data.data()), cols, rows, type, 0,
                    0) != CE_None)
        std::cerr << "write failed: " << path << '\n';
    GDALClose(ds);
}

// A GeoJSON layer in WGS84 lon/lat holding one line or one polygon.
void write_geojson(const std::string& path, const std::vector<lsdtt::lsdtt_lonlat>& pts,
                   bool polygon)
{
    std::filesystem::remove(path);
    GDALDriver* drv = GetGDALDriverManager()->GetDriverByName("GeoJSON");
    GDALDataset* ds = drv->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    OGRSpatialReference wgs;
    wgs.SetWellKnownGeogCS("WGS84");
    wgs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    OGRLayer* layer =
        ds->CreateLayer("feature", &wgs, polygon ? wkbPolygon : wkbLineString, nullptr);
    OGRFeature* f = OGRFeature::CreateFeature(layer->GetLayerDefn());
    if (polygon)
    {
        OGRPolygon poly;
        OGRLinearRing ring;
        for (const auto& p : pts) ring.addPoint(p.longitude, p.latitude);
        ring.closeRings();
        poly.addRing(&ring);
        f->SetGeometry(&poly);
    }
    else
    {
        OGRLineString line;
        for (const auto& p : pts) line.addPoint(p.longitude, p.latitude);
        f->SetGeometry(&line);
    }
    if (layer->CreateFeature(f) != OGRERR_NONE) std::cerr << "feature write failed: " << path << '\n';
    OGRFeature::DestroyFeature(f);
    GDALClose(ds);
}

// Projected points to lon/lat on the projected CRS's own geographic CRS.
std::vector<lsdtt::lsdtt_lonlat> to_lonlat(int epsg, const std::vector<double>& xy)
{
    OGRSpatialReference s;
    s.importFromEPSG(epsg);
    OGRSpatialReference* g = s.CloneGeogCS();
    s.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    g->SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    OGRCoordinateTransformation* ct = OGRCreateCoordinateTransformation(&s, g);
    std::vector<lsdtt::lsdtt_lonlat> out;
    for (std::size_t i = 0; i + 1 < xy.size(); i += 2)
    {
        double x = xy[i], y = xy[i + 1];
        ct->Transform(1, &x, &y);
        out.push_back({x, y});
    }
    OGRCoordinateTransformation::DestroyCT(ct);
    OGRSpatialReference::DestroySpatialReference(g);
    return out;
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

double off_lattice(double a, double b, double cell)
{
    const double k = (a - b) / cell;
    return std::abs(k - std::round(k));
}

double valid_fraction(const std::string& path)
{
    const lsdtt::lsdtt_raster_tile t = tile_of(path);
    std::size_t n = 0;
    for (auto v : t.mask) n += v ? 1 : 0;
    return t.mask.empty() ? 0.0 : static_cast<double>(n) / t.mask.size();
}

// Fold a difference of axis azimuths into [0, 90].
double axis_angle_diff(double a, double b)
{
    double d = std::fmod(std::abs(a - b), 180.0);
    return d > 90.0 ? 180.0 - d : d;
}

struct round_trip
{
    double off_x = 1.0;
    double off_y = 1.0;
    double median = -1.0;
    double max = -1.0;
    std::size_t n = 0;
    bool crs_ok = false;
};

// Rotate a frame raster back with the frame file and compare it with the
// source array on the source lattice.
round_trip rotate_back_and_compare(const std::string& rotated,
                                   const lsdtt::lsdtt_rotation_frame_file& ff,
                                   const std::string& out, const std::string& method,
                                   const std::vector<float>& source, double sx0, double sy1,
                                   double cell, int scols, int srows,
                                   const std::string& source_wkt)
{
    round_trip rt;
    lsdtt::rotate_back_raster(rotated, ff, method, 0, out);
    const lsdtt::lsdtt_raster_meta m = meta_of(out);
    const double sy0 = sy1 - srows * cell;
    rt.off_x = off_lattice(m.x_min, sx0, cell);
    rt.off_y = off_lattice(m.y_min, sy0, cell);
    rt.crs_ok = lsdtt::crs_equal_semantically(m.crs_wkt, source_wkt) &&
                std::abs(m.cell_size - cell) < 1e-9;
    const lsdtt::lsdtt_raster_tile t = tile_of(out);
    const int dc = static_cast<int>(std::lround((m.x_min - sx0) / cell));
    const int dr = static_cast<int>(std::lround((sy1 - m.y_max()) / cell));
    std::vector<double> diffs;
    for (int r = 0; r < m.n_rows; ++r)
        for (int c = 0; c < m.n_cols; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(r) * m.n_cols + c;
            const int sr = r + dr, sc = c + dc;
            if (sr < 0 || sr >= srows || sc < 0 || sc >= scols || !t.mask[i]) continue;
            const float s = source[static_cast<std::size_t>(sr) * scols + sc];
            if (s == -9999.0f) continue;
            diffs.push_back(std::abs(t.data[i] - s));
        }
    std::sort(diffs.begin(), diffs.end());
    rt.n = diffs.size();
    if (!diffs.empty())
    {
        rt.median = diffs[diffs.size() / 2];
        rt.max = diffs.back();
    }
    return rt;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "usage: lsdtt-test-oblique-frame-preprocess <work_dir>\n";
        return 2;
    }
    work = argv[1];
    std::filesystem::create_directories(work);
    lsdtt::gdal_init();

    // --- Fixtures ---------------------------------------------------------------
    const std::string utm = wkt_of_epsg(32612);
    const std::string nad83 = wkt_of_epsg(26912);
    const double sx0 = 580003.7;
    const double sy1 = 4252601.3;
    const double cell = 10.0;
    const int scols = 300;
    const int srows = 240;
    const double sy0 = sy1 - srows * cell;
    const double cx = sx0 + 0.5 * scols * cell;
    const double cy = sy0 + 0.5 * srows * cell;
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
    const std::string src_nad83 = work + "/source_nad83.tif";
    write_tif(src, utm, sx0, sy1, cell, scols, srows, dem, GDT_Float32, -9999.0);
    write_tif(src_nad83, nad83, sx0, sy1, cell, scols, srows, dem, GDT_Float32, -9999.0);

    // Line feature through the centre, grid azimuth about 35 deg, 2.07 km.
    const std::vector<double> line_xy = {cx - 600.0, cy - 850.0, cx - 100.0, cy - 100.0,
                                         cx + 600.0, cy + 850.0};
    const std::string line_file = work + "/line.geojson";
    write_geojson(line_file, to_lonlat(32612, line_xy), false);

    // --- S: cell size and target ------------------------------------------------
    const lsdtt::lsdtt_raster_meta src_meta = meta_of(src);
    check("S1 projected source keeps its own cell size",
          lsdtt::frame_cell_size_for_raster(src_meta, -1.0) == cell);
    check("S2 grid_spacing > 0 wins", lsdtt::frame_cell_size_for_raster(src_meta, 30.0) == 30.0);
    {
        // 1 arc-second geographic DEM near 38 N.
        const std::string geo = work + "/source_geo.tif";
        const double arcsec = 1.0 / 3600.0;
        std::vector<float> g(100 * 100, 100.0f);
        write_tif(geo, wkt_of_epsg(4326), -110.1, 38.4, arcsec, 100, 100, g, GDT_Float32,
                  -9999.0);
        const lsdtt::lsdtt_raster_meta gm = meta_of(geo);
        const double gcell = lsdtt::frame_cell_size_for_raster(gm, -1.0);
        std::cout << "      1 arc-second at 38.4 N: frame cell size " << gcell << " m\n";
        check("S3 geographic source: snapped analysis pixel size (25 m)", gcell == 25.0);

        const lsdtt::lsdtt_oblique_frame gf = lsdtt::frame_from_points(
            {-110.09, 38.39}, {-110.08, 38.38}, 500.0, -1.0, gcell);
        const lsdtt::lsdtt_rotation_target gt = lsdtt::rotation_target_for_raster(gm, geo, gf);
        check("S4 geographic source: target UTM 12N at the frame cell size on (0, 0)",
              gt.epsg == 32612 && gt.crs_wkt.empty() && gt.cell_size == gcell &&
              gt.origin_x == 0.0 && gt.origin_y == 0.0 && gt.source_dem_fname == "source_geo.tif");
    }

    // --- L: line feature ----------------------------------------------------------
    const std::string base_utm = lsdtt::base_geographic_wkt_of(utm);
    const lsdtt::lsdtt_frame_from_feature_result lr =
        lsdtt::frame_from_feature_file(line_file, "", 800.0, -1.0, cell, base_utm);
    check("L1 line feature reads", lr.ok);
    {
        const auto ends = to_lonlat(32612, {line_xy[0], line_xy[1], line_xy[4], line_xy[5]});
        const lsdtt::lsdtt_oblique_frame expect =
            lsdtt::frame_from_points(ends[0], ends[1], 800.0, -1.0, cell, base_utm);
        const auto& f = lr.frame;
        std::cout << std::setprecision(10) << "      line frame azimuth " << f.azimuth_degrees
                  << " deg, " << f.n_cols() << " x " << f.n_rows() << " cells\n";
        check("L2 line frame equals frame_from_points on its end points",
              f.definition == lsdtt::oblique_frame_definition::line &&
              std::abs(f.azimuth_degrees - expect.azimuth_degrees) < 1e-9 &&
              std::abs(f.centre.longitude - expect.centre.longitude) < 1e-9 &&
              std::abs(f.centre.latitude - expect.centre.latitude) < 1e-9 &&
              f.n_cols() == expect.n_cols() && f.n_rows() == expect.n_rows());
    }

    const std::string rot = work + "/line_oblique.tif";
    const std::string rot_frame = work + "/line_oblique_rotation_frame.txt";
    lsdtt::lsdtt_rotation_provenance prov;
    prov.created_by = "test";
    prov.feature_fname = "line.geojson";
    lsdtt::lsdtt_warp_options wopts;
    wopts.resampling_method = "bilinear";
    const lsdtt::lsdtt_rotate_to_frame_result rr =
        lsdtt::rotate_raster_to_frame(src, lr.frame, rot, rot_frame, prov, wopts);
    const lsdtt::lsdtt_rotation_frame_file ff = lsdtt::read_rotation_frame(rot_frame);
    {
        const lsdtt::lsdtt_raster_meta m = meta_of(rot);
        check("L3 rotated DEM on the frame grid",
              lsdtt::crs_equal_semantically(m.crs_wkt, lr.frame.to_wkt()) &&
              m.n_cols == lr.frame.n_cols() && m.n_rows == lr.frame.n_rows() &&
              std::abs(m.x_min - lr.frame.x_min) < 1e-6);
        check("L4 frame file reads, no warnings", ff.ok && ff.warnings.empty());
        check("L5 target is the source lattice (EPSG:32612, 10 m, source corner)",
              ff.target.epsg == 32612 && ff.target.cell_size == cell &&
              ff.target.origin_x == sx0 && ff.target.origin_y == sy0 &&
              ff.target.source_dem_fname == "source_utm.tif" &&
              rr.target.origin_x == sx0);
        check("L6 provenance file names only",
              ff.provenance.rotated_dem_fname == "line_oblique.tif" &&
              ff.provenance.feature_fname == "line.geojson" &&
              ff.frame.definition == lsdtt::oblique_frame_definition::line);
    }
    round_trip utm_rt;
    {
        utm_rt = rotate_back_and_compare(rot, ff, work + "/line_dem_rotated_back.tif", "bilinear",
                                         dem, sx0, sy1, cell, scols, srows, utm);
        std::cout << std::setprecision(6) << "      DEM round trip: " << utm_rt.n
                  << " cells, median |dz| " << utm_rt.median << " m, max " << utm_rt.max
                  << " m, off-lattice " << utm_rt.off_x << ", " << utm_rt.off_y << "\n";
        check("L7 DEM rotated back lands on the source pixels",
              utm_rt.crs_ok && utm_rt.off_x < 1e-6 && utm_rt.off_y < 1e-6);
        check("L8 DEM round trip: median |dz| < 0.05 m, max < 1 m",
              utm_rt.n > 1000 && utm_rt.median < 0.05 && utm_rt.max < 1.0);
    }
    {
        // An ID product made in the frame: 100 m x 100 m blocks.
        const lsdtt::lsdtt_raster_meta m = meta_of(rot);
        std::vector<std::int32_t> ids(static_cast<std::size_t>(m.n_cols) * m.n_rows);
        std::set<std::int32_t> made;
        for (int r = 0; r < m.n_rows; ++r)
            for (int c = 0; c < m.n_cols; ++c)
            {
                const std::int32_t id = 1 + (c / 10) + 1000 * (r / 10);
                ids[static_cast<std::size_t>(r) * m.n_cols + c] = id;
                made.insert(id);
            }
        const std::string idf = work + "/line_ids.tif";
        write_tif(idf, lr.frame.to_wkt(), m.x_min, m.y_max(), m.cell_size, m.n_cols, m.n_rows,
                  ids, GDT_Int32, -9999.0);
        const std::string out = work + "/line_ids_rotated_back.tif";
        lsdtt::rotate_back_raster(idf, ff, "near", 0, out);
        const lsdtt::lsdtt_raster_meta om = meta_of(out);
        const lsdtt::lsdtt_raster_tile t = tile_of(out);
        std::set<std::int32_t> seen;
        bool foreign = false;
        for (std::size_t i = 0; i < t.data.size(); ++i)
        {
            if (!t.mask[i]) continue;
            const auto v = static_cast<std::int32_t>(t.data[i]);
            seen.insert(v);
            if (!made.count(v)) foreign = true;
        }
        std::cout << "      ID round trip: " << seen.size() << " of " << made.size()
                  << " IDs, foreign values " << (foreign ? "yes" : "no") << "\n";
        check("L9 ID product (near) lands on the source pixels with no foreign values",
              off_lattice(om.x_min, sx0, cell) < 1e-6 && off_lattice(om.y_min, sy0, cell) < 1e-6 &&
              !foreign && seen.size() == made.size());
    }

    // --- G: polygon feature -----------------------------------------------------
    {
        // A 2000 m x 700 m rectangle, long side at grid azimuth 120 deg.
        const double th = 120.0 * M_PI / 180.0;
        const double ux = std::sin(th), uy = std::cos(th);  // long-side direction
        const double vx = uy, vy = -ux;                       // short-side direction
        std::vector<double> xy;
        for (const auto& s : {std::pair<double, double>{-1, -1}, {1, -1}, {1, 1}, {-1, 1}})
        {
            xy.push_back(cx + s.first * 1000.0 * ux + s.second * 350.0 * vx);
            xy.push_back(cy + s.first * 1000.0 * uy + s.second * 350.0 * vy);
        }
        const std::string poly_file = work + "/polygon.geojson";
        write_geojson(poly_file, to_lonlat(32612, xy), true);
        const auto pr = lsdtt::frame_from_feature_file(poly_file, "", 0.0, 0.0, cell, base_utm);
        check("G1 polygon feature reads", pr.ok);
        const auto mids = to_lonlat(32612, {cx - 1000.0 * ux, cy - 1000.0 * uy, cx + 1000.0 * ux,
                                            cy + 1000.0 * uy});
        const lsdtt::lsdtt_oblique_frame expect =
            lsdtt::frame_from_points(mids[0], mids[1], 700.0, 0.0, cell, base_utm);
        const auto& f = pr.frame;
        const double daz = axis_angle_diff(f.azimuth_degrees, expect.azimuth_degrees);
        const double dlon_m = (f.centre.longitude - expect.centre.longitude) * 111320.0 *
                              std::cos(f.centre.latitude * M_PI / 180.0);
        const double dlat_m = (f.centre.latitude - expect.centre.latitude) * 110574.0;
        std::cout << std::setprecision(8) << "      polygon frame azimuth " << f.azimuth_degrees
                  << " vs " << expect.azimuth_degrees << " deg (diff " << daz
                  << "), centre offset " << std::hypot(dlon_m, dlat_m) << " m, length "
                  << f.axis_length_metres << " m, width " << f.width_metres << " m\n";
        check("G2 polygon frame: azimuth within 0.05 deg, centre within 1 m",
              f.definition == lsdtt::oblique_frame_definition::polygon && daz < 0.05 &&
              std::hypot(dlon_m, dlat_m) < 1.0);
        check("G3 polygon frame: length and width within 2 m of 2000 x 700",
              std::abs(f.axis_length_metres - 2000.0) < 2.0 &&
              std::abs(f.width_metres - 700.0) < 2.0);
    }

    // --- F: valid-data footprint --------------------------------------------------
    {
        // Valid cells: a 2400 m x 600 m strip at grid azimuth 60 deg.
        const double th = 60.0 * M_PI / 180.0;
        const double ux = std::sin(th), uy = std::cos(th);
        std::vector<float> strip(dem.size(), -9999.0f);
        std::size_t n_valid = 0;
        for (int r = 0; r < srows; ++r)
            for (int c = 0; c < scols; ++c)
            {
                const double x = sx0 + (c + 0.5) * cell - cx;
                const double y = sy1 - (r + 0.5) * cell - cy;
                const double along = x * ux + y * uy;
                const double across = x * uy - y * ux;
                if (std::abs(along) <= 1200.0 && std::abs(across) <= 300.0)
                {
                    const std::size_t i = static_cast<std::size_t>(r) * scols + c;
                    strip[i] = dem[i];
                    ++n_valid;
                }
            }
        const std::string strip_src = work + "/source_strip.tif";
        write_tif(strip_src, utm, sx0, sy1, cell, scols, srows, strip, GDT_Float32, -9999.0);
        const auto fr = lsdtt::frame_from_valid_data(strip_src, 0.0, cell);
        check("F1 footprint frame builds", fr.ok);
        const auto mids = to_lonlat(32612, {cx - 1200.0 * ux, cy - 1200.0 * uy, cx + 1200.0 * ux,
                                            cy + 1200.0 * uy});
        const lsdtt::lsdtt_oblique_frame expect =
            lsdtt::frame_from_points(mids[0], mids[1], 600.0, 0.0, cell, base_utm);
        const double daz = axis_angle_diff(fr.frame.azimuth_degrees, expect.azimuth_degrees);
        std::cout << std::setprecision(8) << "      footprint azimuth "
                  << fr.frame.azimuth_degrees << " vs " << expect.azimuth_degrees
                  << " deg (diff " << daz << "), length " << fr.frame.axis_length_metres
                  << " m, width " << fr.frame.width_metres << " m\n";
        check("F2 footprint: data_footprint, azimuth within 0.5 deg, size within 2 cells",
              fr.frame.definition == lsdtt::oblique_frame_definition::data_footprint &&
              daz < 0.5 && std::abs(fr.frame.axis_length_metres - 2400.0) < 2 * cell &&
              std::abs(fr.frame.width_metres - 600.0) < 2 * cell);

        const std::string out = work + "/strip_oblique.tif";
        const std::string out_frame = work + "/strip_oblique_rotation_frame.txt";
        lsdtt::rotate_raster_to_frame(strip_src, fr.frame, out, out_frame, {}, wopts);
        const double before = static_cast<double>(n_valid) / dem.size();
        const double after = valid_fraction(out);
        std::cout << "      valid fraction: source grid " << before << ", footprint frame "
                  << after << "\n";
        check("F3 footprint frame holds more valid cells per cell than the source grid",
              after > 0.9 && after > before);
        const auto sff = lsdtt::read_rotation_frame(out_frame);
        check("F4 frame file: frame_source data_footprint, no rebuild warning",
              sff.ok && sff.warnings.empty() &&
              sff.frame.definition == lsdtt::oblique_frame_definition::data_footprint &&
              sff.frame.buffer_all_sides);
    }

    // --- D: NAD83 source ----------------------------------------------------------
    {
        const std::string base_nad = lsdtt::base_geographic_wkt_of(nad83);
        const auto nr = lsdtt::frame_from_feature_file(line_file, "", 800.0, -1.0, cell, base_nad);
        check("D1 line feature reads onto NAD83", nr.ok);
        OGRSpatialReference fs;
        fs.importFromWkt(nr.frame.to_wkt().c_str());
        const char* datum = fs.GetAttrValue("DATUM");
        const std::string dname = datum ? datum : "";
        std::cout << "      frame datum: " << dname << "\n";
        check("D2 frame on the source datum (NAD83)",
              dname.find("North") != std::string::npos && dname.find("1983") != std::string::npos);

        const std::string out = work + "/nad83_oblique.tif";
        const std::string out_frame = work + "/nad83_oblique_rotation_frame.txt";
        lsdtt::rotate_raster_to_frame(src_nad83, nr.frame, out, out_frame, {}, wopts);
        const auto nff = lsdtt::read_rotation_frame(out_frame);
        check("D3 target EPSG:26912 on the source corner",
              nff.ok && nff.target.epsg == 26912 && nff.target.origin_x == sx0 &&
              nff.target.origin_y == sy0);
        const round_trip nrt = rotate_back_and_compare(out, nff, work + "/nad83_rotated_back.tif",
                                                       "bilinear", dem, sx0, sy1, cell, scols,
                                                       srows, nad83);
        std::cout << "      NAD83 round trip: " << nrt.n << " cells, median |dz| " << nrt.median
                  << " m, max " << nrt.max << " m\n";
        check("D4 NAD83 round trip lands on the source pixels and matches the WGS84 one",
              nrt.crs_ok && nrt.off_x < 1e-6 && nrt.off_y < 1e-6 &&
              std::abs(nrt.median - utm_rt.median) < 0.01 && nrt.max < 1.0);
    }

    // --- R: refusals ----------------------------------------------------------------
    {
        const std::string empty = work + "/all_nodata.tif";
        std::vector<float> nd(100, -9999.0f);
        write_tif(empty, utm, sx0, sy1, cell, 10, 10, nd, GDT_Float32, -9999.0);
        const auto e = lsdtt::frame_from_valid_data(empty, 0.0, cell);
        check("R1 footprint of an all-nodata raster refused",
              !e.ok && e.message.find("no valid cells") != std::string::npos);

        const std::string nocrs = work + "/no_crs.tif";
        write_tif(nocrs, "", sx0, sy1, cell, scols, srows, dem, GDT_Float32, -9999.0);
        const auto n = lsdtt::frame_from_valid_data(nocrs, 0.0, cell);
        check("R2 footprint of a raster with no CRS refused",
              !n.ok && n.message.find("no CRS") != std::string::npos);
    }

    std::cout << "\n" << (checks - failures) << " / " << checks << " checks passed\n";
    return failures == 0 ? 0 : 1;
}
