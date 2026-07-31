// test_divide_migration.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// Unit driver for the divide-migration compute layer
// (src/cpu/lsdtt_divide_migration_cpu.cpp), plus the generator for the
// synthetic fixture used by the CLI end-to-end tests.
//
//   argv[1] (optional) - directory into which the CLI fixture (a three-valley
//                        DEM and the outlet point layer) is written. Without it
//                        only the in-memory unit checks run.
//
// The checks target the contract stated in lsdtt_divide_migration_cpu.hpp:
//
//   N1  the nearest divide point is the closest point on a SEGMENT, not the
//       closest vertex -- every fixture here puts the answer strictly between
//       two vertices, because a projection bug is invisible when it does not
//   N2  a head sitting on its divide yields a NaN gradient rather than an
//       infinity, and does not crash
//   N3  negative relief is a legitimate result and is passed through
//   N4  a head whose unit borders no divide keeps NaN fields and an empty
//       pair_id_nearest
//   N5  a nearest point on nodata leaves the elevation-derived fields NaN but
//       still reports distance_to_divide_m and pair_id_nearest
//   N6  the per-unit index: a head is NEVER matched to a divide between two
//       OTHER units, even when that divide is far closer in space. This is the
//       property the per-unit index design exists to guarantee
//   B1  apply_common_base_level_to_chi truncates below the base level, re-zeroes
//       chi at the root of EACH retained fragment, and leaves outlet_flat alone
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <ogr_spatialref.h>

#include "lsdtt_chi_extraction_cpu.hpp"
#include "lsdtt_data_loader.hpp"
#include "lsdtt_divide_migration_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_vector_utils.hpp"

namespace
{

int g_n_checks = 0;
int g_n_failed = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!ok && !detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << "\n";
}

bool close_to(double got, double expected, double tol)
{
    return std::isfinite(got) && std::fabs(got - expected) <= tol;
}

std::string num(double v)
{
    return std::to_string(v);
}

// -----------------------------------------------------------------------------
// A raster whose lower-left corner is at (0, 0) in projected units, so that
// test geometry can be written in plain metres.
// -----------------------------------------------------------------------------
lsdtt::lsdtt_raster_cpu make_uniform_elevation(
    int n_rows, int n_cols, double cell_size, float z)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = n_rows;
    meta.n_cols = n_cols;
    meta.cell_size = cell_size;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();
    meta.x_min = 0.0;
    meta.y_min = 0.0;

    std::vector<float> data(static_cast<std::size_t>(n_rows) * n_cols, z);
    std::vector<std::uint8_t> mask(data.size(), 1);
    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

lsdtt::lsdtt_divide_segment_record make_segment(
    const std::string& pair_id, int unit_a, int unit_b,
    std::vector<lsdtt::lsdtt_xy> points)
{
    lsdtt::lsdtt_divide_segment_record rec;
    rec.pair_id = pair_id;
    rec.unit_a = unit_a;
    rec.unit_b = unit_b;
    rec.segment_id = 1;
    rec.points = std::move(points);
    return rec;
}

lsdtt::lsdtt_channel_head_metric_record make_head(
    int head_id, int unit_id, double x, double y, double z)
{
    lsdtt::lsdtt_channel_head_metric_record rec;
    rec.head_id = head_id;
    rec.outlet_id = unit_id;
    rec.unit_id = unit_id;
    rec.head_x = x;
    rec.head_y = y;
    rec.head_elevation_m = z;
    return rec;
}

// =============================================================================
// N1, N2, N3 -- projection, zero distance, negative relief
// =============================================================================
void test_nearest_point_geometry()
{
    // Elevation is a uniform 40 m plate, so every sampled divide elevation is
    // known without reference to the grid geometry.
    const lsdtt::lsdtt_raster_cpu elevation =
        make_uniform_elevation(40, 40, 10.0, 40.0f);

    // One horizontal divide, y = 100, running from x = 0 to x = 300. The two
    // vertices are deliberately far from every head below.
    std::vector<lsdtt::lsdtt_divide_segment_record> segments;
    segments.push_back(make_segment("1|2", 1, 2,
        {lsdtt::lsdtt_xy{0.0, 100.0}, lsdtt::lsdtt_xy{300.0, 100.0}}));

    std::vector<lsdtt::lsdtt_channel_head_metric_record> heads;
    heads.push_back(make_head(1, 1, 150.0, 60.0, 10.0));   // N1: mid-segment
    heads.push_back(make_head(2, 1, 150.0, 100.0, 10.0));  // N2: on the divide
    heads.push_back(make_head(3, 2, 150.0, 60.0, 100.0));  // N3: head above divide

    const std::size_t n_matched =
        lsdtt::compute_nearest_divide_metrics(segments, elevation, heads);

    check("N0 every head matched", n_matched == 3,
          "n_matched=" + std::to_string(n_matched));

    // N1: the answer is the foot of the perpendicular at x = 150, 40 m away.
    // Nearest-VERTEX logic would return (0,100) or (300,100) at 156.5 m.
    check("N1 nearest point is mid-segment, not a vertex",
          close_to(heads[0].nearest_divide_x, 150.0, 1e-9) &&
          close_to(heads[0].nearest_divide_y, 100.0, 1e-9),
          "got (" + num(heads[0].nearest_divide_x) + ", " +
                    num(heads[0].nearest_divide_y) + ")");
    check("N1 distance is the perpendicular distance",
          close_to(heads[0].distance_to_divide_m, 40.0, 1e-9),
          num(heads[0].distance_to_divide_m));
    check("N1 pair_id is the owning segment", heads[0].pair_id_nearest == "1|2",
          heads[0].pair_id_nearest);
    check("N1 relief and gradient from the sampled divide elevation",
          close_to(heads[0].relief_m, 30.0, 1e-9) &&
          close_to(heads[0].gradient_to_divide, 0.75, 1e-9),
          "relief=" + num(heads[0].relief_m) +
          " gradient=" + num(heads[0].gradient_to_divide));

    // N2: zero distance must not divide by zero.
    check("N2 zero distance is reported as zero",
          close_to(heads[1].distance_to_divide_m, 0.0, 1e-12),
          num(heads[1].distance_to_divide_m));
    check("N2 zero distance yields a NaN gradient, not an infinity",
          std::isnan(heads[1].gradient_to_divide),
          num(heads[1].gradient_to_divide));
    check("N2 relief is still reported at zero distance",
          close_to(heads[1].relief_m, 30.0, 1e-9), num(heads[1].relief_m));

    // N3: a head above its divide is a real result, not an error to clamp.
    check("N3 negative relief is preserved, not clamped or absolutised",
          close_to(heads[2].relief_m, -60.0, 1e-9), num(heads[2].relief_m));
    check("N3 negative gradient is preserved",
          close_to(heads[2].gradient_to_divide, -1.5, 1e-9),
          num(heads[2].gradient_to_divide));
}

// =============================================================================
// N4 -- a unit that borders no divide
// =============================================================================
void test_unit_without_divide()
{
    const lsdtt::lsdtt_raster_cpu elevation =
        make_uniform_elevation(40, 40, 10.0, 40.0f);

    std::vector<lsdtt::lsdtt_divide_segment_record> segments;
    segments.push_back(make_segment("1|2", 1, 2,
        {lsdtt::lsdtt_xy{0.0, 100.0}, lsdtt::lsdtt_xy{300.0, 100.0}}));

    std::vector<lsdtt::lsdtt_channel_head_metric_record> heads;
    heads.push_back(make_head(1, 7, 150.0, 60.0, 10.0));   // unit 7 has no divide
    heads.push_back(make_head(2, 2, 150.0, 60.0, 10.0));

    const std::size_t n_matched =
        lsdtt::compute_nearest_divide_metrics(segments, elevation, heads);

    check("N4 only the head whose unit borders a divide is matched",
          n_matched == 1, "n_matched=" + std::to_string(n_matched));
    check("N4 unmatched head keeps an empty pair_id_nearest",
          heads[0].pair_id_nearest.empty(), heads[0].pair_id_nearest);
    check("N4 unmatched head keeps NaN nearest-divide fields",
          std::isnan(heads[0].nearest_divide_x) &&
          std::isnan(heads[0].nearest_divide_y) &&
          std::isnan(heads[0].divide_elevation_m) &&
          std::isnan(heads[0].distance_to_divide_m) &&
          std::isnan(heads[0].relief_m) &&
          std::isnan(heads[0].gradient_to_divide));
    check("N4 the neighbouring unit is unaffected",
          heads[1].pair_id_nearest == "1|2" &&
          close_to(heads[1].distance_to_divide_m, 40.0, 1e-9));
}

// =============================================================================
// N5 -- nearest point on nodata
// =============================================================================
void test_nearest_point_on_nodata()
{
    // 40 x 40 at 10 m from (0,0), so y_max = 400. The divide is put at y = 105
    // rather than on a cell boundary, so the cell containing the nearest point
    // (155, 105) is unambiguously col = 15, row = floor((400-105)/10) = 29.
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = 40;
    meta.n_cols = 40;
    meta.cell_size = 10.0;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();
    meta.x_min = 0.0;
    meta.y_min = 0.0;

    std::vector<float> data(1600, 40.0f);
    std::vector<std::uint8_t> mask(1600, 1);
    const std::size_t hole = static_cast<std::size_t>(29) * 40 + 15;
    data[hole] = -9999.0f;
    mask[hole] = 0;
    lsdtt::lsdtt_raster_cpu elevation(std::move(data), std::move(mask), meta);

    std::vector<lsdtt::lsdtt_divide_segment_record> segments;
    segments.push_back(make_segment("1|2", 1, 2,
        {lsdtt::lsdtt_xy{0.0, 105.0}, lsdtt::lsdtt_xy{300.0, 105.0}}));

    std::vector<lsdtt::lsdtt_channel_head_metric_record> heads;
    heads.push_back(make_head(1, 1, 155.0, 60.0, 10.0));

    const std::size_t n_matched =
        lsdtt::compute_nearest_divide_metrics(segments, elevation, heads);

    check("N5 a nodata nearest cell still counts as matched", n_matched == 1,
          "n_matched=" + std::to_string(n_matched));
    check("N5 geometry is still reported over nodata",
          heads[0].pair_id_nearest == "1|2" &&
          close_to(heads[0].distance_to_divide_m, 45.0, 1e-9) &&
          close_to(heads[0].nearest_divide_x, 155.0, 1e-9),
          "pair=" + heads[0].pair_id_nearest +
          " d=" + num(heads[0].distance_to_divide_m));
    check("N5 elevation-derived fields are NaN over nodata",
          std::isnan(heads[0].divide_elevation_m) &&
          std::isnan(heads[0].relief_m) &&
          std::isnan(heads[0].gradient_to_divide),
          "z=" + num(heads[0].divide_elevation_m) +
          " relief=" + num(heads[0].relief_m));
}

// =============================================================================
// N6 -- three units: never match across a divide the head's unit does not share
// =============================================================================
void test_three_unit_isolation()
{
    const lsdtt::lsdtt_raster_cpu elevation =
        make_uniform_elevation(60, 120, 10.0, 40.0f);

    // Units 1 | 2 | 3 side by side. The 1|2 divide runs up x = 100 and the 2|3
    // divide up x = 900. Unit 1 and unit 3 share no divide at all.
    std::vector<lsdtt::lsdtt_divide_segment_record> segments;
    segments.push_back(make_segment("1|2", 1, 2,
        {lsdtt::lsdtt_xy{100.0, 0.0}, lsdtt::lsdtt_xy{100.0, 500.0}}));
    segments.push_back(make_segment("2|3", 2, 3,
        {lsdtt::lsdtt_xy{900.0, 0.0}, lsdtt::lsdtt_xy{900.0, 500.0}}));

    std::vector<lsdtt::lsdtt_channel_head_metric_record> heads;
    // A unit-3 head parked 10 m from the 1|2 divide. The globally nearest
    // geometry is the WRONG answer: unit 3 does not border that divide, so the
    // 2|3 divide 790 m away is the only admissible match.
    heads.push_back(make_head(1, 3, 110.0, 250.0, 10.0));
    // The mirror case for unit 1 against the 2|3 divide.
    heads.push_back(make_head(2, 1, 890.0, 250.0, 10.0));
    // A unit-2 head, which borders both, must take the genuinely nearest one.
    heads.push_back(make_head(3, 2, 200.0, 250.0, 10.0));

    const std::size_t n_matched =
        lsdtt::compute_nearest_divide_metrics(segments, elevation, heads);
    check("N6 all three heads matched", n_matched == 3,
          "n_matched=" + std::to_string(n_matched));

    check("N6 a unit-3 head is never matched to the 1|2 divide",
          heads[0].pair_id_nearest == "2|3" &&
          close_to(heads[0].distance_to_divide_m, 790.0, 1e-9),
          "pair=" + heads[0].pair_id_nearest +
          " d=" + num(heads[0].distance_to_divide_m));
    check("N6 a unit-1 head is never matched to the 2|3 divide",
          heads[1].pair_id_nearest == "1|2" &&
          close_to(heads[1].distance_to_divide_m, 790.0, 1e-9),
          "pair=" + heads[1].pair_id_nearest +
          " d=" + num(heads[1].distance_to_divide_m));
    check("N6 a shared unit takes the genuinely nearest of its own divides",
          heads[2].pair_id_nearest == "1|2" &&
          close_to(heads[2].distance_to_divide_m, 100.0, 1e-9),
          "pair=" + heads[2].pair_id_nearest +
          " d=" + num(heads[2].distance_to_divide_m));
}

// =============================================================================
// B1 -- common base level applied to chi
// =============================================================================
void test_common_base_level()
{
    // A 1-column, 8-row grid at 10 m. Flat index i is row i, so elevation is
    // set per node directly.
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = 8;
    meta.n_cols = 1;
    meta.cell_size = 10.0;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();
    meta.x_min = 0.0;
    meta.y_min = 0.0;

    std::vector<float> data = {0.0f, 10.0f, 20.0f, 30.0f,
                               40.0f, 50.0f, 60.0f, 70.0f};
    std::vector<std::uint8_t> mask(8, 1);
    lsdtt::lsdtt_raster_cpu elevation(std::move(data), std::move(mask), meta);

    lsdtt::lsdtt_chi_collection_result chi;

    // Basin A: a single chain 5 -> 4 -> 3 -> 2 -> 1 -> 0, chi = node index.
    lsdtt::lsdtt_chi_basin_result basin_a;
    basin_a.outlet_flat = 0;
    basin_a.outlet_chi = 0.0f;
    basin_a.flats = {0, 1, 2, 3, 4, 5};
    basin_a.chi = {0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f};
    basin_a.receiver_flats = {-1, 0, 1, 2, 3, 4};
    chi.basins.push_back(basin_a);

    // Basin B: node 0 with TWO tributaries, so truncating just above node 0
    // leaves two disjoint retained fragments, each of which must be re-zeroed
    // at its own root. A single global offset would fail this.
    //   0(z=0) <- 4(z=40) <- 6(z=60)
    //   0(z=0) <- 5(z=50) <- 7(z=70)
    lsdtt::lsdtt_chi_basin_result basin_b;
    basin_b.outlet_flat = 0;
    basin_b.outlet_chi = 0.0f;
    basin_b.flats = {0, 4, 5, 6, 7};
    basin_b.chi = {0.0f, 1.0f, 2.0f, 3.0f, 5.0f};
    basin_b.receiver_flats = {-1, 0, 0, 4, 5};
    chi.basins.push_back(basin_b);

    const std::size_t n_removed =
        lsdtt::apply_common_base_level_to_chi(chi, elevation, 35.0);

    // Basin A keeps z >= 35, i.e. nodes 4 and 5 (three removed);
    // basin B keeps 4, 5, 6, 7 (one removed).
    check("B1 removed-node count is the sum over basins", n_removed == 5,
          "n_removed=" + std::to_string(n_removed));

    const lsdtt::lsdtt_chi_basin_result& a = chi.basins[0];
    check("B1 basin A retains exactly the nodes at or above base level",
          a.flats.size() == 2 && a.flats[0] == 4 && a.flats[1] == 5,
          "n=" + std::to_string(a.flats.size()));
    check("B1 basin A chi is re-zeroed at the lowest retained node",
          a.chi.size() == 2 && close_to(a.chi[0], 0.0, 1e-6) &&
          close_to(a.chi[1], 1.0, 1e-6));
    check("B1 basin A outlet_flat is left untouched", a.outlet_flat == 0,
          std::to_string(a.outlet_flat));
    check("B1 basin A outlet_chi is reset to zero",
          close_to(a.outlet_chi, 0.0, 1e-9), num(a.outlet_chi));

    const lsdtt::lsdtt_chi_basin_result& b = chi.basins[1];
    check("B1 basin B retains all four upslope nodes",
          b.flats.size() == 4 && b.flats[0] == 4 && b.flats[1] == 5 &&
          b.flats[2] == 6 && b.flats[3] == 7,
          "n=" + std::to_string(b.flats.size()));
    // Fragment roots are nodes 4 (chi 1) and 5 (chi 2); each subtree is offset
    // by its OWN root, never by a single collection-wide value.
    check("B1 chi is re-zeroed per retained fragment, not globally",
          b.chi.size() == 4 &&
          close_to(b.chi[0], 0.0, 1e-6) && close_to(b.chi[1], 0.0, 1e-6) &&
          close_to(b.chi[2], 2.0, 1e-6) && close_to(b.chi[3], 3.0, 1e-6),
          "chi=" + num(b.chi[0]) + "," + num(b.chi[1]) + "," +
                   num(b.chi[2]) + "," + num(b.chi[3]));
    check("B1 basin B receivers are kept in step with the retained nodes",
          b.receiver_flats.size() == 4 && b.receiver_flats[0] == 0 &&
          b.receiver_flats[1] == 0 && b.receiver_flats[2] == 4 &&
          b.receiver_flats[3] == 5);
    check("B1 basin B outlet_flat is left untouched", b.outlet_flat == 0,
          std::to_string(b.outlet_flat));
}

// =============================================================================
// CLI fixture: a deterministic three-valley DEM plus its outlet layer
//
// Elevation rises to the north and is corrugated across x with troughs at
// x = 2000, 6000 and 10000 m, so the grid holds three parallel south-draining
// valleys separated by two ridges. Units 1 and 3 are therefore NOT adjacent,
// which is exactly the arrangement N6 pins in memory.
//
// The outlet layer carries two cluster fields so that one fixture drives all
// three CLI cases:
//   cluster_id  -> 1, 1, 2  (two units: the middle divide is internal)
//   cluster_all -> 7, 7, 7  (one unit: must be rejected)
// =============================================================================
void write_cli_fixture(const std::string& out_dir)
{
    const int    n_rows = 80;
    const int    n_cols = 120;
    const double cell   = 100.0;

    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = n_rows;
    meta.n_cols = n_cols;
    meta.cell_size = cell;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();

    const double x0 = meta.x_min;
    const double y0 = meta.y_min;

    constexpr double k_pi = 3.14159265358979323846;
    // The troughs sit on cell CENTRES (col 19, 59, 99), not on cell edges. A
    // trough on an edge is a perfect elevation tie between two columns, which
    // splits each valley into two half-basins and destroys the adjacency the
    // divide extraction depends on.
    const double period = 4000.0;     // valley spacing (m)
    const double x_trough = 1950.0;   // easting of the first trough

    std::vector<float> data(static_cast<std::size_t>(n_rows) * n_cols, 0.0f);
    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const double x = (col + 0.5) * cell;               // local easting
            const double y = (n_rows - row - 0.5) * cell;      // local northing

            // Cross-valley corrugation: 0 at every trough, 120 m at every crest.
            const double cross =
                60.0 * (1.0 - std::cos(2.0 * k_pi * (x - x_trough) / period));
            // Down-valley fall to the south, plus a hillslope corrugation in y
            // that seeds tributaries (and therefore junctions) without
            // disturbing the trough, where cross is zero.
            const double down = 0.05 * y;
            const double trib = 1.0 + 0.12 * std::cos(2.0 * k_pi * y / 1600.0);

            data[static_cast<std::size_t>(row) * n_cols + col] =
                static_cast<float>(down + cross * trib);
        }
    }

    std::filesystem::create_directories(out_dir);
    const std::string dem_path = out_dir + "/divmig_valleys.tif";
    lsdtt::lsdtt_data_loader writer;
    writer.write_full(data, meta, dem_path,
                      lsdtt::raster_output_format::gtiff);
    std::cout << "Wrote CLI fixture DEM: " << dem_path << "\n";

    // One outlet at the mouth of each valley: the bottom row, on the trough.
    const int outlet_row = n_rows - 1;
    const std::vector<int> outlet_cols = {19, 59, 99};
    const std::vector<std::string> cluster_id  = {"1", "1", "2"};
    const std::vector<std::string> cluster_all = {"7", "7", "7"};

    lsdtt::lsdtt_vector_schema schema;
    schema.field_names = {"point_id", "cluster_id", "cluster_all"};
    schema.field_types = {lsdtt::lsdtt_vector_field_type::int64_type,
                          lsdtt::lsdtt_vector_field_type::string_type,
                          lsdtt::lsdtt_vector_field_type::string_type};

    const std::string outlet_path = out_dir + "/divmig_outlets.fgb";
    lsdtt::lsdtt_data_loader vwriter;
    vwriter.begin_vector_write(outlet_path, "outlets", schema, meta.crs_wkt,
                               lsdtt::vector_output_format::flatgeobuf, false);
    for (std::size_t i = 0; i < outlet_cols.size(); ++i)
    {
        const double x = x0 + (outlet_cols[i] + 0.5) * cell;
        const double y = y0 + (n_rows - outlet_row - 0.5) * cell;
        vwriter.write_point(x, y, {static_cast<long int>(i + 1),
                                   cluster_id[i], cluster_all[i]});
    }
    vwriter.finish_vector_write();
    std::cout << "Wrote CLI fixture outlets: " << outlet_path << "\n";

    // The same three outlets as a WGS84 CSV. Both forms are needed: the CLI
    // reads cluster fields from a CSV by column name and from a vector layer
    // through OGR, and only the CSV path is usable today (see the CLI defect
    // reported alongside these tests).
    OGRSpatialReference utm;
    OGRSpatialReference wgs84;
    if (utm.importFromWkt(meta.crs_wkt.c_str()) != OGRERR_NONE ||
        wgs84.importFromEPSG(4326) != OGRERR_NONE)
    {
        throw std::runtime_error("could not build the fixture CRS transform");
    }
    utm.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    wgs84.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    OGRCoordinateTransformation* ct =
        OGRCreateCoordinateTransformation(&utm, &wgs84);
    if (ct == nullptr)
    {
        throw std::runtime_error("could not create the fixture CRS transform");
    }

    const std::string csv_path = out_dir + "/divmig_outlets.csv";
    std::ofstream csv(csv_path);
    if (!csv) throw std::runtime_error("could not write " + csv_path);
    csv << "latitude,longitude,point_id,cluster_id,cluster_all\n";
    csv << std::setprecision(12);
    for (std::size_t i = 0; i < outlet_cols.size(); ++i)
    {
        double x = x0 + (outlet_cols[i] + 0.5) * cell;
        double y = y0 + (n_rows - outlet_row - 0.5) * cell;
        if (!ct->Transform(1, &x, &y))
        {
            throw std::runtime_error("fixture outlet transform failed");
        }
        csv << y << "," << x << "," << (i + 1) << ","
            << cluster_id[i] << "," << cluster_all[i] << "\n";
    }
    csv.close();
    OCTDestroyCoordinateTransformation(ct);
    std::cout << "Wrote CLI fixture outlet CSV: " << csv_path << "\n";
}

} // namespace

int main(int argc, char* argv[])
{
    std::cout << "=== lsdtt3 divide-migration unit tests ===\n";

    test_nearest_point_geometry();
    test_unit_without_divide();
    test_nearest_point_on_nodata();
    test_three_unit_isolation();
    test_common_base_level();

    if (argc > 1)
    {
        try
        {
            write_cli_fixture(argv[1]);
        }
        catch (const std::exception& e)
        {
            check("F1 CLI fixture written", false, e.what());
        }
    }

    std::cout << "\nChecks: " << g_n_checks
              << "   failures: " << g_n_failed << "\n";
    std::cout << (g_n_failed == 0 ? "RESULT: PASS\n" : "RESULT: FAIL\n");
    return g_n_failed == 0 ? 0 : 1;
}
