// test_lithology_channel_stats.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Tests for k_sn categorised by lithologic unit (lsdtt_lithology_channel_stats_cpu).
//
//   T1  contact distance on a two-unit grid (exact)
//   T3  node -> segment on a hand-built fit result (exact)
//   T7  summary and histogram bookkeeping (counts add up, rows present)
//   T10 known answer, n = 1: a steady-state landscape with K1 = 1e-5 west of a
//       vertical contact and K2 = 4e-5 east of it, U = 2e-4, map area, A0 = 1:
//       k_sn = (U/K)^(1/n) = 20 and 5, ratio 4
//   T11 the same with m = 0.9, n = 2: k_sn = sqrt(20) and sqrt(5), ratio 2
//
// The steady state integrates chi with the area of each cell (left end), the
// chi of lsdtt-chi-analysis is trapezoidal, so the fitted k_sn sits slightly
// off (U/K)^(1/n). Measured 2026-10-17 (120 x 160 at 30 m, seed 11): medians
// 1.14 % and 1.25 % high for both n, ratio within 0.108 %. Tolerances: the
// ratio at twice the measured error (0.25 %); the medians at 2 %, the cap
// (twice the measured 1.25 % would exceed it).
//
// Usage:
//   lsdtt-test-lithology-channel-stats <work_dir>
//   lsdtt-test-lithology-channel-stats --write-fixtures <dir> <tiny_dem.asc>
// -----------------------------------------------------------------------------

#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include "../../src/core/lsdtt_data_loader.hpp"
#include "../../src/core/lsdtt_gdal_init.hpp"
#include "../../src/cpu/lsdtt_chi_extraction_cpu.hpp"
#include "../../src/cpu/lsdtt_chi_network_analysis_cpu.hpp"
#include "../../src/cpu/lsdtt_flowinfo_cpu.hpp"
#include "../../src/cpu/lsdtt_lithology_channel_stats_cpu.hpp"
#include "../../src/cpu/lsdtt_lithology_map_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../../src/cpu/lsdtt_steady_state_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_synthetic_surfaces_cpu.hpp"

using namespace lsdtt;

namespace
{

int g_failures = 0;

void check(const std::string& what, bool ok, const std::string& detail = "")
{
    std::cout << (ok ? "PASS " : "FAIL ") << what;
    if (!detail.empty()) std::cout << " -- " << detail;
    std::cout << "\n";
    if (!ok) ++g_failures;
}

std::string fmt(double v)
{
    std::ostringstream s;
    s << std::setprecision(6) << v;
    return s.str();
}

lsdtt_raster_meta utm_grid(int nr, int nc, double dx)
{
    lsdtt_raster_meta m;
    m.n_rows = nr; m.n_cols = nc; m.cell_size = dx;
    m.x_min = 500000.0; m.y_min = 5000000.0; m.no_data_value = -9999.0;
    m.set_crs_from_epsg(32632);
    return m;
}

lsdtt_raster_cpu blank(const lsdtt_raster_meta& m)
{
    const std::size_t n = static_cast<std::size_t>(m.n_rows) * m.n_cols;
    return lsdtt_raster_cpu(std::vector<float>(n, 0.0f), std::vector<std::uint8_t>(n, 1), m);
}

// Units by column: code 1 for col < split, else 2.
lsdtt_lithology_units_on_grid column_units(const lsdtt_raster_meta& m, int split)
{
    lsdtt_lithology_units_on_grid u;
    u.code_of_pixel.resize(static_cast<std::size_t>(m.n_rows) * m.n_cols);
    for (int r = 0; r < m.n_rows; ++r)
        for (int c = 0; c < m.n_cols; ++c)
            u.code_of_pixel[static_cast<std::size_t>(r) * m.n_cols + c] = c < split ? 1 : 2;
    u.units = {{1, "west", {"west"}, 0}, {2, "east", {"east"}, 0}};
    return u;
}

// ---------------------------------------------------------------------------
// T1 contact distance
// ---------------------------------------------------------------------------
void test_contact_distance()
{
    std::cout << "-- T1 contact distance --\n";
    const lsdtt_raster_meta m = utm_grid(20, 20, 10.0);
    lsdtt_raster_cpu g = blank(m);
    std::vector<std::uint8_t> mask = g.get_mask();
    mask[static_cast<std::size_t>(5) * 20 + 2] = 0;   // a nodata pixel in the west
    g = lsdtt_raster_cpu(std::vector<float>(g.get_data()), std::move(mask), m);
    const auto u = column_units(m, 10);
    const std::vector<float> d = compute_lithology_contact_distance(g, u, 1);
    auto at = [&](int r, int c) { return d[static_cast<std::size_t>(r) * 20 + c]; };
    check("contact pixels (cols 9 and 10) are 0", at(3, 9) == 0.0f && at(3, 10) == 0.0f);
    check("col 0 is 9 dx = 90 m, col 19 is 90 m", at(0, 0) == 90.0f && at(19, 19) == 90.0f,
          fmt(at(0, 0)) + ", " + fmt(at(19, 19)));
    check("col 4 is 50 m", at(12, 4) == 50.0f, fmt(at(12, 4)));
    check("grid nodata is NaN", std::isnan(at(5, 2)));

    lsdtt_lithology_units_on_grid one = u;
    for (auto& c : one.code_of_pixel) c = 1;
    const auto none = compute_lithology_contact_distance(blank(m), one, 1);
    check("no contact: all NaN", std::isnan(none[0]) && std::isnan(none[399]));

    lsdtt_lithology_units_on_grid hole = one;
    hole.code_of_pixel[static_cast<std::size_t>(10) * 20 + 10] = -1;
    const auto h = compute_lithology_contact_distance(blank(m), hole, 1);
    check("unmapped counts as a unit (its neighbours are contacts)",
          h[static_cast<std::size_t>(10) * 20 + 11] == 0.0f && h[static_cast<std::size_t>(10) * 20 + 13] == 20.0f);
}

// ---------------------------------------------------------------------------
// T3 node -> segment
// ---------------------------------------------------------------------------
void test_node_segments()
{
    std::cout << "-- T3 node -> segment --\n";
    // Main stem 0 <- 1 <- 2 <- 3 <- 4 <- 5 (0 = outlet), tributary 7 <- 6 -> joins at 3
    // (6 drains to 3). Segments: main [5..3], [3..0]; tributary [7..3].
    lsdtt_chi_network_fit_result r;
    r.flats = {0, 1, 2, 3, 4, 5, 6, 7};
    r.receiver_flats = {k_outlet_receiver, 0, 1, 2, 3, 4, 3, 6};
    r.chi.assign(8, 0.0f);
    r.observed_elevation.assign(8, 0.0f);
    lsdtt_chi_network_fit_segment a, b, c;
    a.segment_id = 0; a.upstream_flat = 3; a.downstream_flat = 0; a.k_sn = 1.0f;
    b.segment_id = 1; b.upstream_flat = 5; b.downstream_flat = 3; b.k_sn = 2.0f;
    c.segment_id = 2; c.upstream_flat = 7; c.downstream_flat = 3; c.k_sn = 3.0f;
    r.segments = {a, b, c};
    const auto s = assign_chi_fit_nodes_to_segments(r);
    const std::vector<std::int32_t> want = {0, 0, 0, 0, 1, 1, 2, 2};
    std::string got;
    for (auto v : s) got += std::to_string(v) + " ";
    check("each node on one segment; confluence on the downstream segment; outlet on the last",
          s == want, got);
}

// ---------------------------------------------------------------------------
// The steady-state fixture
// ---------------------------------------------------------------------------
struct steady_fixture
{
    lsdtt_raster_cpu surface;      // the noisy surface the flow is routed on
    lsdtt_raster_cpu elevation;    // the steady state
    lsdtt_lithology_units_on_grid units;
};

const int k_nr = 120, k_nc = 160, k_split = 80;
const double k_dx = 30.0;
const float k_U = 2e-4f, k_K1 = 1e-5f, k_K2 = 4e-5f;

lsdtt_flowinfo d8_flow(const lsdtt_raster_cpu& dem)
{
    lsdtt_flowinfo_options fo;
    fo.boundary_conditions = "oooo";
    fo.n_threads = 1;
    fo.build_topological_order = true;
    fo.build_s_vector = true;
    return lsdtt_flowinfo(dem.get_data_array(), dem.get_mask_array(), dem.get_meta(), fo);
}

lsdtt_raster_cpu tilted_surface()
{
    lsdtt_raster_cpu z = generate_red_noise_surface(k_nr, k_nc, k_dx, 11, 0.0, 3.0, 6, 1);
    std::vector<float> d = z.get_data();
    for (int r = 0; r < k_nr; ++r)
        for (int c = 0; c < k_nc; ++c)
            d[static_cast<std::size_t>(r) * k_nc + c] += static_cast<float>(0.02 * k_dx * r);
    lsdtt_raster_meta m = z.get_meta();
    m.x_min = 500000.0; m.y_min = 5000000.0;
    m.set_crs_from_epsg(32632);
    return fill_pits(lsdtt_raster_cpu(std::move(d), std::vector<std::uint8_t>(z.get_mask()), m),
                     1e-4f, 1);
}

steady_fixture make_steady(float m_exp, float n_exp)
{
    steady_fixture f;
    f.surface = tilted_surface();
    const lsdtt_flowinfo flow = d8_flow(f.surface);
    lsdtt_fluvial_steady_state_options o;
    o.U = k_U;
    o.K = k_K1;
    o.m = m_exp;
    o.n = n_exp;
    o.divide_area_by_dx = false;   // map area, as chi with A0 = 1
    o.erodibility_field.resize(static_cast<std::size_t>(k_nr) * k_nc);
    for (int r = 0; r < k_nr; ++r)
        for (int c = 0; c < k_nc; ++c)
            o.erodibility_field[static_cast<std::size_t>(r) * k_nc + c] = c < k_split ? k_K1 : k_K2;
    f.elevation = compute_d8_fluvial_steady_state(flow, o).elevation;
    f.units = column_units(f.surface.get_meta(), k_split);
    return f;
}

lsdtt_chi_network_fit_collection_result fit_on(const lsdtt_raster_cpu& surface,
                                               const lsdtt_raster_cpu& elevation, double mn)
{
    const lsdtt_flowinfo flow = d8_flow(surface);
    const lsdtt_raster_result npx = flow.get_n_contributing_pixels_raster();
    std::vector<std::int64_t> outlets;
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(flow.full_grid_size()); ++i)
        if (flow.is_valid_cell(i) && flow.is_outlet(i) && npx.data[static_cast<std::size_t>(i)] >= 50.0f)
            outlets.push_back(i);
    lsdtt_chi_options co;
    co.m_over_n = mn;
    co.A_0 = 1.0;
    lsdtt_chi_collection_result chi = get_chi_from_outlets_zero_base_threshold(flow, outlets, 50, co);
    (void)reject_chi_basins_with_few_channel_pixels(chi, 30);
    lsdtt_chi_network_fit_collection_result fits;
    for (const auto& b : chi.basins) fits.basins.push_back(fit_chi_network(b, elevation, {}));
    return fits;
}

const lsdtt_chi_ksn_lithology_summary_row* row_for(const lsdtt_chi_ksn_lithology_summary& s, int code)
{
    for (const auto& r : s.rows)
        if (r.lith_code == code) return &r;
    return nullptr;
}

void test_known_answer(float m_exp, float n_exp, double tol_rel, double tol_ratio)
{
    const std::string tag = "n = " + fmt(n_exp) + ": ";
    std::cout << "-- T10/T11 steady state, m = " << m_exp << ", n = " << n_exp << " --\n";
    const steady_fixture f = make_steady(m_exp, n_exp);
    const double mn = static_cast<double>(m_exp) / n_exp;
    const auto fits = fit_on(f.surface, f.elevation, mn);
    const std::vector<float> cd = compute_lithology_contact_distance(f.surface, f.units, 1);
    const auto before = fits.basins.empty() ? lsdtt_chi_network_fit_result{} : fits.basins[0];
    const auto nodes = build_chi_lithology_node_table(fits, f.surface.get_meta(), f.units, cd);
    check(tag + "fits untouched by the tagging (observation only)",
          !fits.basins.empty() && fits.basins[0].flats == before.flats &&
              fits.basins[0].segments.size() == before.segments.size());
    const auto s = summarise_chi_ksn_by_lithology(nodes, f.units, mn, 1.0);

    const double want1 = std::pow(k_U / k_K1, 1.0 / n_exp), want2 = std::pow(k_U / k_K2, 1.0 / n_exp);
    const auto* w = row_for(s, 1);
    const auto* e = row_for(s, 2);
    if (!w || !e || w->n_nodes == 0 || e->n_nodes == 0)
    {
        check(tag + "both units have channel nodes", false);
        return;
    }
    const double err1 = std::fabs(w->k_sn_median / want1 - 1.0);
    const double err2 = std::fabs(e->k_sn_median / want2 - 1.0);
    const double ratio = w->k_sn_median / e->k_sn_median, want_ratio = want1 / want2;
    const double err_ratio = std::fabs(ratio / want_ratio - 1.0);
    std::cout << "   measured: west median " << w->k_sn_median << " (want " << want1 << ", "
              << w->n_nodes << " nodes, " << w->n_segments << " segments), east "
              << e->k_sn_median << " (want " << want2 << ", " << e->n_nodes << " nodes, "
              << e->n_segments << " segments), ratio " << ratio << "\n";
    check(tag + "west (K1) median k_sn = (U/K1)^(1/n)", err1 <= tol_rel,
          "rel err " + fmt(err1) + ", tolerance " + fmt(tol_rel));
    check(tag + "east (K2) median k_sn = (U/K2)^(1/n)", err2 <= tol_rel,
          "rel err " + fmt(err2) + ", tolerance " + fmt(tol_rel));
    check(tag + "ratio = (K2/K1)^(1/n)", err_ratio <= tol_ratio,
          fmt(ratio) + " vs " + fmt(want_ratio) + ", rel err " + fmt(err_ratio));

    // T7 bookkeeping on the same tables.
    const auto* all = row_for(s, -2);
    const auto* un = row_for(s, -1);
    check(tag + "all = west + east + unmapped (nodes, length)",
          all && un && all->n_nodes == w->n_nodes + e->n_nodes + un->n_nodes &&
              std::fabs(all->channel_length_m -
                        (w->channel_length_m + e->channel_length_m + un->channel_length_m)) < 1e-6 * all->channel_length_m);
    check(tag + "no unmapped nodes", un && un->n_nodes == 0 && s.n_nodes_unmapped == 0);
    check(tag + "m/n and A0 copied", s.m_over_n == mn && s.A_0 == 1.0);

    lsdtt_chi_ksn_histogram_options ho;
    ho.n_bins = 25;
    const auto h = make_chi_ksn_histogram_by_lithology(nodes, f.units, mn, 1.0, ho);
    std::map<int, std::size_t> sum, rows;
    for (const auto& r : h.rows) { sum[r.lith_code] += r.n_nodes; ++rows[r.lith_code]; }
    check(tag + "histogram: 27 rows (25 bins + below + above) per unit, unmapped and all",
          rows.size() == 4 && rows[1] == 27 && rows[2] == 27 && rows[-1] == 27 && rows[-2] == 27);
    check(tag + "histogram counts = summary node counts",
          sum[1] == w->n_nodes && sum[2] == e->n_nodes && sum[-2] == all->n_nodes);
    ho.log_bins = true;
    const auto hl = make_chi_ksn_histogram_by_lithology(nodes, f.units, mn, 1.0, ho);
    std::size_t lsum = 0;
    for (const auto& r : hl.rows) if (r.lith_code == -2) lsum += r.n_nodes;
    check(tag + "log histogram keeps every node", lsum == all->n_nodes);
}

// ---------------------------------------------------------------------------
// --write-fixtures: a steady-state DEM and its unit map for the CLI, and a
// two-unit code raster on tiny_dem.asc's grid for the observation-only check.
// ---------------------------------------------------------------------------
void write_polygon_units(const std::string& path, const lsdtt_raster_meta& m, int split)
{
    gdal_init();
    std::filesystem::remove(path);
    GDALDataset* ds = GetGDALDriverManager()->GetDriverByName("GPKG")->Create(
        path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    OGRSpatialReference srs;
    srs.importFromEPSG(32632);
    srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    OGRLayer* layer = ds->CreateLayer("units", &srs, wkbPolygon, nullptr);
    OGRFieldDefn fd("lith_primary", OFTString);
    (void)layer->CreateField(&fd);
    const double xs = m.x_min + split * m.cell_size, x1 = m.x_min + m.n_cols * m.cell_size;
    const double y0 = m.y_min, y1 = m.y_min + m.n_rows * m.cell_size;
    for (int k = 0; k < 2; ++k)
    {
        OGRFeature* f = OGRFeature::CreateFeature(layer->GetLayerDefn());
        f->SetField("lith_primary", k == 0 ? "sandstone" : "shale");
        const double a = k == 0 ? m.x_min : xs, b = k == 0 ? xs : x1;
        OGRLinearRing ring;
        ring.addPoint(a, y0); ring.addPoint(b, y0); ring.addPoint(b, y1); ring.addPoint(a, y1);
        ring.addPoint(a, y0);
        OGRPolygon p;
        p.addRing(&ring);
        f->SetGeometry(&p);
        if (layer->CreateFeature(f) != OGRERR_NONE) std::cerr << "could not write a unit polygon\n";
        OGRFeature::DestroyFeature(f);
    }
    GDALClose(ds);
}

int write_fixtures(const std::string& dir, const std::string& tiny_dem)
{
    std::filesystem::create_directories(dir);
    gdal_init();
    const steady_fixture f = make_steady(0.45f, 1.0f);
    lsdtt_data_loader w;
    w.write_full(f.elevation.get_data(), f.elevation.get_meta(), dir + "/lithks_steady_dem.tif",
                 raster_output_format::gtiff);
    write_polygon_units(dir + "/lithks_units.gpkg", f.elevation.get_meta(), k_split);

    lsdtt_data_loader l;
    l.open_raster(tiny_dem);
    const lsdtt_raster_meta tm = l.get_meta();
    std::vector<float> codes(static_cast<std::size_t>(tm.n_rows) * tm.n_cols);
    for (int r = 0; r < tm.n_rows; ++r)
        for (int c = 0; c < tm.n_cols; ++c)
            codes[static_cast<std::size_t>(r) * tm.n_cols + c] = c < tm.n_cols / 2 ? 1.0f : 2.0f;
    w.write_full(codes, tm, dir + "/lithks_tiny_codes.tif", raster_output_format::gtiff);
    // tiny_dem beside its code raster, so one read_path serves both.
    const std::filesystem::path src(tiny_dem);
    for (const char* ext : {".asc", ".prj"})
    {
        std::filesystem::path from = src;
        from.replace_extension(ext);
        if (std::filesystem::exists(from))
            std::filesystem::copy_file(from, std::filesystem::path(dir) / from.filename(),
                                       std::filesystem::copy_options::overwrite_existing);
    }
    std::cout << "fixtures written to " << dir << "\n";
    return 0;
}

} // namespace

int main(int argc, char* argv[])
{
    if (argc > 3 && std::string(argv[1]) == "--write-fixtures") return write_fixtures(argv[2], argv[3]);
    const std::string dir = argc > 1 ? argv[1] : ".";
    std::filesystem::create_directories(dir);
    gdal_init();
    std::cout << "== k_sn by lithologic unit ==\n";
    test_contact_distance();
    test_node_segments();
    test_known_answer(0.45f, 1.0f, 0.02, 0.0025);
    test_known_answer(0.9f, 2.0f, 0.02, 0.0025);
    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
