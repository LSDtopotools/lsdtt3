// test_erodibility_calibration.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Tests for the erodibility calibration from basin erosion rates: the parts
// built so far.
//
//   T-Y   York regression (fit_york_line, lsdtt_statistics):
//         Y1 Pearson's data with York's weights against the published values
//            (Cantrell 2008, ACP 8, 5477, Tables 1 and 2) and against an
//            independent solver (scipy.odr / ODRPACK);
//         Y2 an exact line; Y3 sigma_x = 0 equals weighted least squares;
//         Y4 degenerate input returns converged = false without throwing.
//   T-R   the sample reader and units (lsdtt_erodibility_calibration_cpu):
//         R1 unit factors; R2 the units column; R3 missing columns; R4 the
//         sample id fallback; R5 a UTM point layer reads the same pixels as
//         the CSV; plus bad_rate / missing_uncertainty and the reported area.
//   S-a   the shared K formula and its A0 cancellation; the K units string.
//   S-c   the shared property-table writer.
//   D     the two declare_* parameter groups.
//   C0    basin completeness at a grid edge under open boundaries.
//
//   F1    known answer: n and K from basins with per-basin uplift;
//   F2    known answer: K per lithologic unit (mono and unit_in_basin), and
//         its agreement with relative K at equal uplift.
//
// Usage: lsdtt-test-erodibility-calibration <work_dir>
//        lsdtt-test-erodibility-calibration --write-fixtures <dir>
//        lsdtt-test-erodibility-calibration --check-erodibility-calibration-fit ...
//        lsdtt-test-erodibility-calibration --check-erodibility-calibration-unit ...
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include "../../src/core/lsdtt_common_parameters.hpp"
#include "../../src/core/lsdtt_csv_table.hpp"
#include "../../src/core/lsdtt_data_loader.hpp"
#include "../../src/core/lsdtt_gdal_init.hpp"
#include "../../src/core/lsdtt_param_parser.hpp"
#include "../../src/cpu/lsdtt_basin_cpu.hpp"
#include "../../src/cpu/lsdtt_erodibility_calibration_cpu.hpp"
#include "../../src/cpu/lsdtt_flowinfo_cpu.hpp"
#include "../../src/cpu/lsdtt_lithology_relative_erodibility_cpu.hpp"
#include "../../src/cpu/lsdtt_chi_extraction_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_synthetic_surfaces_cpu.hpp"
#include "../../src/cpu/lsdtt_steady_state_cpu.hpp"
#include "../../src/utils/lsdtt_statistics.hpp"

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

template <typename F>
std::string thrown_message(F f)
{
    try { f(); }
    catch (const std::exception& e) { return e.what(); }
    return "";
}

bool contains(const std::string& s, const std::string& p) { return s.find(p) != std::string::npos; }
bool close(double a, double b, double tol) { return std::fabs(a - b) <= tol; }

std::string fmt(double v)
{
    std::ostringstream s;
    s << std::setprecision(15) << v;
    return s.str();
}

// ---------------------------------------------------------------------------
// York
// ---------------------------------------------------------------------------
void test_york()
{
    // Pearson's data with York's weights, transcribed from Cantrell (2008),
    // Table 1 (x, w_x, y, w_y); sigma = 1 / sqrt(w).
    const std::vector<double> x{0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4};
    const std::vector<double> wx{1000.0, 1000.0, 500.0, 800.0, 200.0, 80.0, 60.0, 20.0, 1.8, 1.0};
    const std::vector<double> y{5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5};
    const std::vector<double> wy{1.0, 1.8, 4.0, 8.0, 20.0, 20.0, 70.0, 70.0, 100.0, 500.0};
    std::vector<double> sx, sy;
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        sx.push_back(1.0 / std::sqrt(wx[i]));
        sy.push_back(1.0 / std::sqrt(wy[i]));
    }
    const lsdtt::lsdtt_york_fit f = lsdtt::fit_york_line(x, y, sx, sy);
    std::cout << "York, Pearson-York data: slope " << fmt(f.slope) << " intercept " << fmt(f.intercept)
              << " se " << fmt(f.se_slope) << " / " << fmt(f.se_intercept) << " scaled "
              << fmt(f.se_slope_scaled) << " / " << fmt(f.se_intercept_scaled) << " mswd "
              << fmt(f.mswd) << " corr " << fmt(f.corr_slope_intercept) << " iterations " << f.n_iter
              << "\n";
    // Cantrell (2008) Table 2, Williamson-York row, as printed (5 significant
    // figures): slope -0.48053 +- 0.0706, intercept 5.4799 +- 0.359. The
    // printed errors equal the analytic ones times sqrt(MSWD).
    check("T-Y1: slope and intercept equal the published values to their printed precision",
          f.converged && close(f.slope, -0.48053, 5e-6) && close(f.intercept, 5.4799, 5e-5));
    check("T-Y1: MSWD-scaled errors equal the published ones to their printed precision",
          close(f.se_slope_scaled, 0.0706, 5e-5) && close(f.se_intercept_scaled, 0.359, 5e-4));
    // scipy.odr (ODRPACK, scipy 1.18.0, sstol = partol = 1e-15) on the same
    // data: slope -0.480533392649, intercept 5.479910153327, unscaled errors
    // 0.0579850066 / 0.2949707257, reduced chi-square 1.4832941493. ODRPACK's
    // own convergence leaves about 1.5e-8 in the slope.
    check("T-Y1: agrees with an independent ODR solver",
          close(f.slope, -0.480533392649, 5e-8) && close(f.intercept, 5.479910153327, 2e-7) &&
          close(f.se_slope, 0.0579850066, 1e-8) && close(f.se_intercept, 0.2949707257, 1e-8) &&
          close(f.mswd, 1.4832941493, 1e-8));

    // Y2: an exact line, any sigmas.
    const std::vector<double> xe{0.0, 1.0, 2.5, 4.0, 7.0};
    std::vector<double> ye, sxe{0.1, 0.3, 0.05, 0.2, 0.4}, sye{0.2, 0.1, 0.3, 0.05, 0.5};
    for (double v : xe) ye.push_back(2.0 + 1.5 * v);
    const auto e = lsdtt::fit_york_line(xe, ye, sxe, sye);
    check("T-Y2: exact line y = 2 + 1.5 x: slope and intercept to 1e-12, MSWD 0",
          e.converged && close(e.slope, 1.5, 1e-12) && close(e.intercept, 2.0, 1e-12) && e.mswd < 1e-20,
          fmt(e.slope) + " " + fmt(e.intercept) + " " + fmt(e.mswd));

    // Y3: sigma_x = 0 is weighted least squares in y.
    const std::vector<double> zero(x.size(), 0.0);
    const auto w = lsdtt::fit_york_line(x, y, zero, sy);
    double S = 0, Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        const double wi = 1.0 / (sy[i] * sy[i]);
        S += wi; Sx += wi * x[i]; Sy += wi * y[i]; Sxx += wi * x[i] * x[i]; Sxy += wi * x[i] * y[i];
    }
    const double D = S * Sxx - Sx * Sx;
    const double b = (S * Sxy - Sx * Sy) / D, a = (Sxx * Sy - Sx * Sxy) / D;
    check("T-Y3: sigma_x = 0 equals weighted least squares to 1e-10",
          w.converged && close(w.slope, b, 1e-10) && close(w.intercept, a, 1e-10) &&
          close(w.se_slope, std::sqrt(S / D), 1e-10) && close(w.se_intercept, std::sqrt(Sxx / D), 1e-10),
          fmt(w.slope) + " vs " + fmt(b));

    // Y4: degenerate input.
    const auto two = lsdtt::fit_york_line({1.0, 2.0}, {1.0, 2.0}, {0.1, 0.1}, {0.1, 0.1});
    const auto flat = lsdtt::fit_york_line({3.0, 3.0, 3.0}, {1.0, 2.0, 4.0}, {0.1, 0.1, 0.1}, {0.1, 0.1, 0.1});
    check("T-Y4: fewer than 3 points or zero x range -> converged = false, NaN, no exception",
          !two.converged && std::isnan(two.slope) && !flat.converged && std::isnan(flat.intercept));
    check("T-Y4: sigma_y <= 0, sigma_x < 0 or unequal lengths throw",
          !thrown_message([] { (void)lsdtt::fit_york_line({1, 2, 3}, {1, 2, 3}, {0, 0, 0}, {1, 0, 1}); }).empty() &&
          !thrown_message([] { (void)lsdtt::fit_york_line({1, 2, 3}, {1, 2, 3}, {0, -1, 0}, {1, 1, 1}); }).empty() &&
          !thrown_message([] { (void)lsdtt::fit_york_line({1, 2, 3}, {1, 2}, {0, 0, 0}, {1, 1, 1}); }).empty());
}

// ---------------------------------------------------------------------------
// S-a and S-c
// ---------------------------------------------------------------------------
void test_seams(const std::string& work)
{
    using lsdtt::stream_power_erodibility_from_erosion_rate;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    check("S-a: K = E / (A0^m k^n)",
          close(stream_power_erodibility_from_erosion_rate(3e-5, 20.0, 0.45, 1.0, 1.0), 3e-5 / 20.0, 1e-20) &&
          close(stream_power_erodibility_from_erosion_rate(3e-5, 20.0, 0.9, 2.0, 4.0),
                3e-5 / (std::pow(4.0, 0.9) * 400.0), 1e-20));
    check("S-a: NaN for E <= 0, k <= 0 or k not finite",
          std::isnan(stream_power_erodibility_from_erosion_rate(0.0, 20.0, 0.45, 1.0, 1.0)) &&
          std::isnan(stream_power_erodibility_from_erosion_rate(1e-5, 0.0, 0.45, 1.0, 1.0)) &&
          std::isnan(stream_power_erodibility_from_erosion_rate(1e-5, -2.0, 0.45, 1.0, 1.0)) &&
          std::isnan(stream_power_erodibility_from_erosion_rate(1e-5, nan, 0.45, 1.0, 1.0)));
    // chi referenced to A0 scales k by A0^-theta, so K does not depend on A0.
    bool a0_ok = true;
    for (double n : {1.0, 1.5, 2.0})
    {
        const double theta = 0.45, m = theta * n, k1 = 37.0, E = 2e-5;
        const double K1 = stream_power_erodibility_from_erosion_rate(E, k1, m, n, 1.0);
        const double A0 = 1e6, kA0 = k1 / std::pow(A0, theta);
        const double KA0 = stream_power_erodibility_from_erosion_rate(E, kA0, m, n, A0);
        a0_ok = a0_ok && close(KA0 / K1, 1.0, 1e-10);
    }
    check("T-K4 (formula): A0 = 1 and A0 = 1e6 give the same K to 1e-10 relative", a0_ok);
    check("T-C5: K units for m = 0.45 and 0.9",
          lsdtt::lithology_erodibility_units_string(0.45) == "m^(0.1) yr^-1" &&
          lsdtt::lithology_erodibility_units_string(0.9) == "m^(-0.8) yr^-1");

    std::filesystem::create_directories(work);
    const std::string path = work + "/property_table.txt";
    std::vector<lsdtt::lsdtt_lithology_erodibility_table_row> units{
        {1, "sandstone", 2.5e-6, "from 7 basins"}, {2, "shale", nan, ""}, {3, "lime, stone", 1e-5, "a \"quoted\" note"}};
    const std::size_t n = lsdtt::write_lithology_erodibility_property_table(
        path, 0.45, 1.0, "test provenance", {"K units: m^(0.1) yr^-1"}, units, "no K");
    std::ifstream in(path);
    std::stringstream ss;
    ss << in.rdbuf();
    const std::string expected =
        "# lsdtt lithology property table v1\n"
        "# erodibility_convention: map_area\n"
        "# m_exponent: 0.45\n"
        "# n_exponent: 1\n"
        "# provenance: test provenance\n"
        "# K units: m^(0.1) yr^-1\n"
        "# not_estimated: 2 shale (no K)\n"
        "code,name,erodibility,reference\n"
        "1,sandstone,2.5e-06,from 7 basins\n"
        "3,\"lime, stone\",1e-05,\"a \"\"quoted\"\" note\"\n";
    check("S-c: the shared property-table writer, byte for byte", n == 2 && ss.str() == expected, ss.str());
    std::vector<lsdtt::lsdtt_lithology_erodibility_table_row> dup{{1, "Shale", 1e-6, ""}, {2, "shale", 2e-6, ""}};
    check("S-c: names equal ignoring case are refused",
          contains(thrown_message([&] {
                       (void)lsdtt::write_lithology_erodibility_property_table(path, 0.45, 1.0, "p", {}, dup, "");
                   }), "same name ignoring case"));
}

// ---------------------------------------------------------------------------
// Reader and units
// ---------------------------------------------------------------------------
lsdtt::lsdtt_raster_meta wv_meta()
{
    lsdtt::lsdtt_raster_meta m;
    m.n_rows = 100;
    m.n_cols = 100;
    m.cell_size = 30.0;
    m.x_min = 636000.0;
    m.y_min = 4303000.0;
    m.set_crs_from_epsg(32617);
    return m;
}

struct pt { double x, y; };

void utm_to_lonlat(const std::vector<pt>& in, std::vector<pt>& out)
{
    OGRSpatialReference utm, wgs;
    utm.importFromEPSG(32617);
    wgs.importFromEPSG(4326);
    utm.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    wgs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    std::unique_ptr<OGRCoordinateTransformation> t(OGRCreateCoordinateTransformation(&utm, &wgs));
    out = in;
    for (auto& p : out) t->Transform(1, &p.x, &p.y);
}

void write_text(const std::string& path, const std::string& text)
{
    std::ofstream(path) << text;
}

void test_reader(const std::string& work)
{
    // R1: the factors, exact.
    check("T-R1: mm_per_kyr = m_per_myr = 1e-6, cm_per_kyr 1e-5, mm_per_yr 1e-3, m_per_yr 1",
          lsdtt::erosion_rate_units_to_m_per_yr("mm_per_kyr") == 1e-6 &&
          lsdtt::erosion_rate_units_to_m_per_yr("m_per_myr") == 1e-6 &&
          lsdtt::erosion_rate_units_to_m_per_yr("cm_per_kyr") == 1e-5 &&
          lsdtt::erosion_rate_units_to_m_per_yr("mm_per_yr") == 1e-3 &&
          lsdtt::erosion_rate_units_to_m_per_yr("m_per_yr") == 1.0 &&
          lsdtt::erosion_rate_units_to_m_per_yr(" MM_PER_KYR ") == 1e-6 &&
          lsdtt::canonical_erosion_rate_units("m_per_myr") == "mm_per_kyr");
    const std::string unknown = thrown_message([] { (void)lsdtt::erosion_rate_units_to_m_per_yr("furlongs"); });
    check("T-R1: an unknown unit throws listing the accepted ones",
          contains(unknown, "erosion_rate_units must be one of mm_per_kyr, m_per_myr, cm_per_kyr, mm_per_yr, m_per_yr") &&
          contains(unknown, "'furlongs'"), unknown);
    check("T-R1: an empty unit throws",
          contains(thrown_message([] { (void)lsdtt::canonical_erosion_rate_units(""); }), "must be one of"));
    check("T-R1: a mass flux is refused with the density hint",
          contains(thrown_message([] { (void)lsdtt::canonical_erosion_rate_units("g_per_cm2_per_yr"); }),
                   "dividing by the rock density"));

    std::filesystem::create_directories(work);
    const lsdtt::lsdtt_raster_meta meta = wv_meta();
    // Pixel centres (row, col) = (10, 20), (50, 50), (90, 70), and a point off the grid.
    auto centre = [&](int r, int c) { return pt{meta.x_min + (c + 0.5) * 30.0, meta.y_min + (meta.n_rows - r - 0.5) * 30.0}; };
    const std::vector<pt> utm{centre(10, 20), centre(50, 50), centre(90, 70), {meta.x_min - 5000.0, meta.y_min}};
    std::vector<pt> ll;
    utm_to_lonlat(utm, ll);
    auto latlon = [&](std::size_t i) {
        std::ostringstream s;
        s << std::setprecision(12) << ll[i].y << "," << ll[i].x;
        return s.str();
    };

    // R2: the units column.
    const std::string good = work + "/samples_units.csv";
    write_text(good, "sample_id,latitude,longitude,erosion_rate,erosion_rate_uncertainty,erosion_rate_units,reported_area_km2\n"
                     "A1," + latlon(0) + ",12.88,2.62,mm_per_kyr,20.14\n"
                     "A2," + latlon(1) + ",0.5,,m_per_myr,\n"
                     "A3," + latlon(2) + ",-1,0.2,mm_per_kyr,3\n"
                     "A4," + latlon(3) + ",7,1,mm_per_kyr,4\n");
    lsdtt::lsdtt_erosion_rate_input_options o;
    o.path = good;
    o.units = "mm_per_kyr";
    const auto r = lsdtt::read_erosion_rate_samples(meta, o);
    check("T-R2: a matching units column (m_per_myr = mm_per_kyr) reads every row in order",
          r.samples.size() == 4 && r.has_units_column && r.units == "mm_per_kyr" &&
          r.samples[0].sample_id == "A1" && r.samples[3].sample_id == "A4" && r.samples[2].input_row == 3);
    check("T-R: E and sigma in m/yr, reported area in m2",
          r.samples.size() == 4 && r.samples[0].E_m_per_yr == 12.88 * 1e-6 &&
          r.samples[0].sigma_E_m_per_yr == 2.62 * 1e-6 && close(r.samples[0].reported_area_m2, 20.14e6, 1e-6) &&
          std::isnan(r.samples[1].reported_area_m2) && r.samples[0].rate_input == 12.88);
    check("T-R: statuses: missing_uncertainty, bad_rate, pending",
          r.samples.size() == 4 && r.samples[0].status.empty() && r.samples[1].status == "missing_uncertainty" &&
          r.samples[2].status == "bad_rate" && r.samples[3].status.empty());
    check("T-R: points land on the expected pixels; the off-grid one is kept as outside_raster_extent",
          r.samples.size() == 4 && r.samples[0].point.row == 10 && r.samples[0].point.col == 20 &&
          r.samples[1].point.row == 50 && r.samples[2].point.col == 70 &&
          r.samples[3].point.status == lsdtt::lsdtt_point_input_status::outside_raster_extent);

    const std::string bad = work + "/samples_units_conflict.csv";
    write_text(bad, "sample_id,latitude,longitude,erosion_rate,erosion_rate_uncertainty,erosion_rate_units\n"
                    "B1," + latlon(0) + ",1,0.1,m_per_yr\n"
                    "B2," + latlon(1) + ",1,0.1,mm_per_yr\n");
    lsdtt::lsdtt_erosion_rate_input_options ob = o;
    ob.path = bad;
    ob.units = "m_per_yr";
    const std::string conflict = thrown_message([&] { (void)lsdtt::read_erosion_rate_samples(meta, ob); });
    check("T-R2: a disagreeing units column throws naming the first row",
          contains(conflict, "row 2 (sample B2)") && contains(conflict, "disagrees with erosion_rate_units"), conflict);

    // R3: missing columns.
    const std::string nocol = work + "/samples_no_rate.csv";
    write_text(nocol, "sample_id,latitude,longitude,rate,sigma\nC1," + latlon(0) + ",1,0.1\n");
    lsdtt::lsdtt_erosion_rate_input_options oc = o;
    oc.path = nocol;
    const std::string missing_rate = thrown_message([&] { (void)lsdtt::read_erosion_rate_samples(meta, oc); });
    check("T-R3: a missing rate column throws listing the columns found",
          contains(missing_rate, "erosion_rate_column 'erosion_rate'") && contains(missing_rate, "rate, sigma"),
          missing_rate);
    oc.rate_column = "rate";
    const std::string missing_unc = thrown_message([&] { (void)lsdtt::read_erosion_rate_samples(meta, oc); });
    check("T-R3: a missing uncertainty column throws listing the columns found",
          contains(missing_unc, "erosion_rate_uncertainty_column 'erosion_rate_uncertainty'") &&
          contains(missing_unc, "sample_id, rate, sigma"), missing_unc);
    lsdtt::lsdtt_erosion_rate_input_options om = o;
    om.path = work + "/no_such_file.csv";
    check("missing file throws 'erosion_rate_fname does not exist'",
          contains(thrown_message([&] { (void)lsdtt::read_erosion_rate_samples(meta, om); }),
                   "erosion_rate_fname does not exist"));

    // R4: no sample id column.
    const std::string noid = work + "/samples_no_id.csv";
    write_text(noid, "latitude,longitude,erosion_rate,erosion_rate_uncertainty\n" + latlon(0) + ",1,0.1\n" +
                     latlon(1) + ",2,0.2\n");
    lsdtt::lsdtt_erosion_rate_input_options oi = o;
    oi.path = noid;
    const auto ri = lsdtt::read_erosion_rate_samples(meta, oi);
    check("T-R4: no sample id column -> row_<n>",
          ri.samples.size() == 2 && !ri.has_sample_id && ri.samples[0].sample_id == "row_1" &&
          ri.samples[1].sample_id == "row_2" && !ri.has_reported_area);

    // R5: a GeoPackage point layer in UTM, no latitude / longitude.
    const std::string gpkg = work + "/samples_utm.gpkg";
    {
        lsdtt::gdal_init();
        std::filesystem::remove(gpkg);
        GDALDataset* ds = GetGDALDriverManager()->GetDriverByName("GPKG")->Create(gpkg.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
        OGRSpatialReference srs;
        srs.importFromEPSG(32617);
        srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        OGRLayer* layer = ds->CreateLayer("samples", &srs, wkbPoint, nullptr);
        for (const char* n : {"sample_id", "erosion_rate_units"})
        {
            OGRFieldDefn fd(n, OFTString);
            (void)layer->CreateField(&fd);
        }
        for (const char* n : {"erosion_rate", "erosion_rate_uncertainty", "reported_area_km2"})
        {
            OGRFieldDefn fd(n, OFTReal);
            (void)layer->CreateField(&fd);
        }
        const char* ids[] = {"A1", "A2", "A3", "A4"};
        const double rate[] = {12.88, 0.5, -1.0, 7.0};
        for (std::size_t i = 0; i < 4; ++i)
        {
            OGRFeature* f = OGRFeature::CreateFeature(layer->GetLayerDefn());
            f->SetField("sample_id", ids[i]);
            f->SetField("erosion_rate_units", "mm_per_kyr");
            f->SetField("erosion_rate", rate[i]);
            if (i != 1) f->SetField("erosion_rate_uncertainty", 0.5);
            OGRPoint p(utm[i].x, utm[i].y);
            f->SetGeometry(&p);
            if (layer->CreateFeature(f) != OGRERR_NONE) std::cerr << "could not write a sample point\n";
            OGRFeature::DestroyFeature(f);
        }
        GDALClose(ds);
    }
    lsdtt::lsdtt_erosion_rate_input_options og = o;
    og.path = gpkg;
    const auto rg = lsdtt::read_erosion_rate_samples(meta, og);
    bool same = rg.samples.size() == r.samples.size();
    for (std::size_t i = 0; same && i < rg.samples.size(); ++i)
        same = rg.samples[i].point.row == r.samples[i].point.row &&
               rg.samples[i].point.col == r.samples[i].point.col &&
               rg.samples[i].point.flat == r.samples[i].point.flat &&
               rg.samples[i].point.status == r.samples[i].point.status &&
               rg.samples[i].sample_id == r.samples[i].sample_id && rg.samples[i].status == r.samples[i].status;
    check("T-R5: a UTM GeoPackage point layer reads the same pixels and statuses as the CSV", same);
}

// ---------------------------------------------------------------------------
// Parameter groups
// ---------------------------------------------------------------------------
void test_parameters()
{
    std::map<std::string, bool> b;
    std::map<std::string, int> i;
    std::map<std::string, float> f;
    std::map<std::string, std::string> s;
    lsdtt::lsdtt_help_registry h;
    lsdtt::declare_erosion_rate_input_parameters(f, s, h);
    lsdtt::declare_erodibility_calibration_parameters(b, i, f, s, h);
    check("declare: erosion-rate input keys and defaults",
          s.count("erosion_rate_fname") && s["erosion_rate_fname"].empty() && s["erosion_rate_units"].empty() &&
          s["erosion_rate_column"] == "erosion_rate" &&
          s["erosion_rate_uncertainty_column"] == "erosion_rate_uncertainty" &&
          s["erosion_rate_reported_area_column"] == "reported_area_km2" &&
          s["erosion_rate_snapping_method"] == "reported_area" &&
          f["erosion_rate_snap_max_distance_metres"] == 300.0f &&
          f["erosion_rate_area_mismatch_tolerance"] == 1.25f &&
          s["erosion_rate_area_mismatch_action"] == "exclude" && s["erosion_rate_nested_basins"] == "keep");
    check("declare: calibration keys and defaults",
          s["erosion_rate_ksn_statistic"] == "reach_median" && i["erosion_rate_min_reaches_per_basin"] == 5 &&
          i["erosion_rate_min_groups_for_n_fit"] == 8 && f["erosion_rate_min_ksn_range_ratio"] == 2.0f &&
          i["erosion_rate_bootstrap_replicates"] == 1000 && i["erosion_rate_bootstrap_random_seed"] == 5489 &&
          b["erosion_rate_fit_n"] && b["write_erodibility_calibration_basin_polygons"] &&
          s["erosion_rate_unit_k_methods"] == "mono,unit_in_basin");
    bool groups = true;
    for (const auto& k : h.key_order())
    {
        const auto& e = h.entries().at(k);
        groups = groups && e.size() == 7 &&
                 (e[4] == "Erosion-rate samples" || e[4] == "Erodibility calibration") &&
                 (k == "erosion_rate_fname" || e[6] == "erosion_rate_fname");
    }
    check("declare: two groups, every key but erosion_rate_fname under it; no n list declared",
          groups && h.key_order().size() == 29 && !s.count("stream_power_n_values") &&
          !s.count("relative_erodibility_n_values"), std::to_string(h.key_order().size()));
}

// ---------------------------------------------------------------------------
// C0: completeness at a grid edge
// ---------------------------------------------------------------------------
void test_edge_completeness()
{
    // A plane falling south, 12 x 9 at 10 m, routed "oooo". Column 4 drains
    // from the north edge (row 0) down to the south edge.
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = 12;
    meta.n_cols = 9;
    meta.cell_size = 10.0;
    meta.set_default_synthetic_georeference();
    std::vector<float> z(static_cast<std::size_t>(meta.n_rows * meta.n_cols));
    std::vector<std::uint8_t> mask(z.size(), 1);
    for (int r = 0; r < meta.n_rows; ++r)
        for (int c = 0; c < meta.n_cols; ++c)
            z[static_cast<std::size_t>(r * meta.n_cols + c)] = static_cast<float>(100 - r + 0.01 * std::abs(c - 4));
    // A nodata block in column 6 at rows 2-3 cuts a basin there.
    mask[static_cast<std::size_t>(2 * meta.n_cols + 6)] = 0;
    mask[static_cast<std::size_t>(3 * meta.n_cols + 6)] = 0;
    const lsdtt::lsdtt_raster_cpu dem(std::move(z), std::move(mask), meta);
    lsdtt::lsdtt_flowinfo_options fo;
    fo.boundary_conditions = "oooo";
    fo.build_topological_order = true;
    fo.build_s_vector = true;
    const lsdtt::lsdtt_flowinfo flow(dem.get_data_array(), dem.get_mask_array(), dem.get_meta(), fo);

    const lsdtt::lsdtt_basin_cpu edge(flow.flat_index(8, 4), flow);
    const lsdtt::lsdtt_basin_cpu cut(flow.flat_index(8, 6), flow);
    std::cout << "C0: basin of (8,4) has " << edge.get_number_of_basin_pixels() << " pixels, complete = "
              << edge.check_basin_completeness() << "; basin of (8,6) below nodata has "
              << cut.get_number_of_basin_pixels() << " pixels, complete = " << cut.check_basin_completeness()
              << "\n";
    check("C0: a basin reaching the open north edge is reported incomplete", !edge.check_basin_completeness());
    check("C0: a basin cut by a nodata block is reported incomplete", !cut.check_basin_completeness());
}

// ---------------------------------------------------------------------------
// Known answers F1 (n and K) and F2 (K per unit) on steady-state surfaces
//
// A red-noise surface (320 x 480 at 30 m, UTM 32N, sd 3 m, 6 blur passes)
// with a tent tilt of 0.02 dx per row towards the north and south edges
// (divide at mid-height), routed "oooo". The basins are the 10 largest
// interior, complete, non-overlapping sub-basins of 4000-10000 px at least
// 8 px from the grid edge (14 qualify), so no basin reaches the open edge.
// Channels start at 50 px. Basin b gets its own uplift U_b
// through an exactly equivalent erodibility field: K_eff = K_unit U_ref / U_b
// on its pixels, U_ref everywhere, so k_sn,b = (U_b / K_unit)^(1/n) exactly
// at steady state and E_b = U_b.
// ---------------------------------------------------------------------------
const int k_ka_nr = 320, k_ka_nc = 480, k_ka_split = 240;
const double k_ka_sd = 3.0;
const int k_ka_blur = 6;
const double k_ka_dx = 30.0, k_ka_Uref = 2e-4, k_ka_K = 1e-5, k_ka_K1 = 1e-5, k_ka_K2 = 4e-5;
const int k_ka_min_px = 4000, k_ka_max_px = 10000, k_ka_n_basins = 10, k_ka_channel_px = 50;
const int k_ka_nested_min_px = 400;
const int k_ka_margin = 8;
const double k_ka_tilt = 0.02;
const int k_ka_order[10] = {3, 7, 0, 9, 5, 1, 8, 2, 6, 4};   // U ranks, so U is not ordered by position

const std::string k_ka_bc = "oooo";

lsdtt::lsdtt_flowinfo ka_flow(const lsdtt::lsdtt_raster_cpu& dem)
{
    lsdtt::lsdtt_flowinfo_options fo;
    fo.boundary_conditions = k_ka_bc;
    fo.n_threads = 1;
    fo.build_topological_order = true;
    fo.build_s_vector = true;
    return lsdtt::lsdtt_flowinfo(dem.get_data_array(), dem.get_mask_array(), dem.get_meta(), fo);
}

lsdtt::lsdtt_raster_cpu ka_surface()
{
    lsdtt::lsdtt_raster_cpu z = lsdtt::generate_red_noise_surface(k_ka_nr, k_ka_nc, k_ka_dx, 11, 0.0, k_ka_sd, k_ka_blur, 1);
    std::vector<float> d = z.get_data();
    for (int r = 0; r < k_ka_nr; ++r)
        for (int c = 0; c < k_ka_nc; ++c)
            d[static_cast<std::size_t>(r) * k_ka_nc + c] +=
                static_cast<float>(k_ka_tilt * k_ka_dx * std::min(r, k_ka_nr - 1 - r));
    lsdtt::lsdtt_raster_meta m = z.get_meta();
    m.x_min = 500000.0;
    m.y_min = 5000000.0;
    m.set_crs_from_epsg(32632);
    return lsdtt::fill_pits(lsdtt::lsdtt_raster_cpu(std::move(d), std::vector<std::uint8_t>(z.get_mask()), m),
                            1e-4f, k_ka_bc, 1);
}

struct ka_fixture
{
    lsdtt::lsdtt_raster_cpu surface, elevation;
    std::vector<std::int64_t> outlets;     // the basins, flat order
    std::vector<double> U;                 // per basin
    std::size_t n_edge_outlets = 0;        // sub-basins that qualify (all of them)
    lsdtt::lsdtt_lithology_units_on_grid units;   // F2: west / east
};

// The basins: interior, complete, non-overlapping sub-basins of
// k_ka_min_px - k_ka_max_px px with every pixel at least k_ka_margin px from
// the grid edge, largest first (ties to the lower flat), at most k_ka_n_basins,
// in flat order. The margin keeps them complete when the CLI re-routes the
// steady state (its own D8 moves a few divide pixels; measured: the routed
// basins share 83-100 % of their pixels with the surface-routed ones).
std::vector<std::int64_t> ka_select_basins(const lsdtt::lsdtt_flowinfo& flow, std::size_t& n_qualifying)
{
    std::vector<std::int64_t> cand, chosen;
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(flow.full_grid_size()); ++i)
    {
        if (!flow.is_valid_cell(i) || flow.is_outlet(i)) continue;
        const std::uint64_t a = flow.s_vector_size(i);
        if (a >= static_cast<std::uint64_t>(k_ka_min_px) && a <= static_cast<std::uint64_t>(k_ka_max_px))
            cand.push_back(i);
    }
    std::stable_sort(cand.begin(), cand.end(), [&](std::int64_t a, std::int64_t b) {
        return flow.s_vector_size(a) > flow.s_vector_size(b);
    });
    std::vector<std::uint8_t> taken(flow.full_grid_size(), 0);
    n_qualifying = 0;
    for (const std::int64_t c : cand)
    {
        if (taken[static_cast<std::size_t>(c)]) continue;
        bool clash = false;
        for (const std::int64_t p : flow.upslope_nodes(c))
            if (taken[static_cast<std::size_t>(p)]) { clash = true; break; }
        if (clash || !lsdtt::lsdtt_basin_cpu(c, flow).check_basin_completeness()) continue;
        bool near_edge = false;
        for (const std::int64_t p : flow.upslope_nodes(c))
        {
            const int r = flow.row_from_flat(p), cc = flow.col_from_flat(p);
            if (r < k_ka_margin || cc < k_ka_margin || r >= k_ka_nr - k_ka_margin || cc >= k_ka_nc - k_ka_margin)
            { near_edge = true; break; }
        }
        if (near_edge) continue;
        ++n_qualifying;
        for (const std::int64_t p : flow.upslope_nodes(c)) taken[static_cast<std::size_t>(p)] = 1;
        if (static_cast<int>(chosen.size()) < k_ka_n_basins) chosen.push_back(c);
    }
    std::sort(chosen.begin(), chosen.end());
    return chosen;
}

// two_units: K1 west of the middle column, K2 east (F2); equal_U: every U_b
// = U_ref. The steady state is built on the routing of the noisy surface;
// f.surface is that surface.
ka_fixture make_ka_fixture(double n_exp, bool two_units, bool equal_U)
{
    ka_fixture f;
    f.surface = ka_surface();
    const lsdtt::lsdtt_flowinfo flow = ka_flow(f.surface);
    f.outlets = ka_select_basins(flow, f.n_edge_outlets);
    const std::size_t npix = static_cast<std::size_t>(k_ka_nr) * k_ka_nc;
    auto k_unit = [&](std::size_t p) {
        if (!two_units) return k_ka_K;
        return static_cast<int>(p % k_ka_nc) < k_ka_split ? k_ka_K1 : k_ka_K2;
    };
    const std::size_t nb = f.outlets.size();
    f.U.assign(nb, k_ka_Uref);
    for (std::size_t b = 0; b < nb; ++b)
        if (!equal_U) f.U[b] = 5e-5 * std::pow(16.0, k_ka_order[b % 10] / 9.0);
    std::vector<float> field(npix);
    for (std::size_t p = 0; p < npix; ++p) field[p] = static_cast<float>(k_unit(p));
    for (std::size_t b = 0; b < nb; ++b)
        for (const std::int64_t p : flow.upslope_nodes(f.outlets[b]))
            field[static_cast<std::size_t>(p)] =
                static_cast<float>(k_unit(static_cast<std::size_t>(p)) * k_ka_Uref / f.U[b]);
    lsdtt::lsdtt_fluvial_steady_state_options o;
    o.U = static_cast<float>(k_ka_Uref);
    o.K = static_cast<float>(two_units ? k_ka_K1 : k_ka_K);
    o.n = static_cast<float>(n_exp);
    o.m = static_cast<float>(0.45 * n_exp);
    o.divide_area_by_dx = false;
    o.erodibility_field = field;
    f.elevation = lsdtt::compute_d8_fluvial_steady_state(flow, o).elevation;
    if (two_units)
    {
        f.units.code_of_pixel.resize(npix);
        for (std::size_t p = 0; p < npix; ++p)
            f.units.code_of_pixel[p] = static_cast<int>(p % k_ka_nc) < k_ka_split ? 1 : 2;
        f.units.units = {{1, "west", {"west"}, 0}, {2, "east", {"east"}, 0}};
    }
    // Diagnostic: the same outlets on the routing of the steady state itself.
    return f;
}

struct ka_run
{
    std::vector<lsdtt::lsdtt_erosion_rate_sample> samples;
    lsdtt::lsdtt_erosion_rate_basins_result basins;
    lsdtt::lsdtt_chi_network_fit_collection_result fits;
    std::vector<lsdtt::lsdtt_erosion_rate_basin_ksn> ksn;
    lsdtt::lsdtt_erodibility_calibration_result cal;
    std::vector<std::size_t> sample_basin;   // basin of each sample row (-1 none)
};

lsdtt::lsdtt_erosion_rate_sample ka_sample(const lsdtt::lsdtt_flowinfo& flow, std::int64_t flat,
                                           double U, double area_factor, const std::string& id)
{
    lsdtt::lsdtt_erosion_rate_sample s;
    s.sample_id = id;
    const auto& m = flow.meta();
    if (flat >= 0)
    {
        s.point.status = lsdtt::lsdtt_point_input_status::success;
        s.point.row = flow.row_from_flat(flat);
        s.point.col = flow.col_from_flat(flat);
        s.point.flat = flat;
        s.point.x = m.x_min + (s.point.col + 0.5) * m.cell_size;
        s.point.y = m.y_min + (m.n_rows - s.point.row - 0.5) * m.cell_size;
        s.reported_area_m2 = area_factor * static_cast<double>(flow.s_vector_size(flat)) * m.cell_size * m.cell_size;
    }
    else
        s.point.status = lsdtt::lsdtt_point_input_status::outside_raster_extent;
    s.rate_input = U * 1e6;            // mm/kyr
    s.uncertainty_input = 0.1 * U * 1e6;
    s.E_m_per_yr = s.rate_input * 1e-6;
    s.sigma_E_m_per_yr = s.uncertainty_input * 1e-6;
    return s;
}

// The samples of the plan: the basin outlets, a nested sample in basin 3, a
// duplicate of basin 5, basin 8's outlet with twice its area, a point off the
// grid. Then the whole pipeline with snapping "none".
ka_run run_ka(const ka_fixture& f, const std::vector<double>& n_values, bool with_units,
              bool extras = true)
{
    ka_run run;
    const lsdtt::lsdtt_flowinfo flow = ka_flow(f.surface);
    for (std::size_t b = 0; b < f.outlets.size(); ++b)
    {
        run.samples.push_back(ka_sample(flow, f.outlets[b], f.U[b], 1.0, "B" + std::to_string(b)));
        run.sample_basin.push_back(b);
    }
    if (extras && f.outlets.size() > 8)
    {
        // Nested: the largest cell of basin 3 with 400 px <= A <= 0.6 A_outlet.
        const std::uint64_t a3 = flow.s_vector_size(f.outlets[3]);
        std::int64_t best = -1;
        std::uint64_t best_a = 0;
        for (const std::int64_t p : flow.upslope_nodes(f.outlets[3]))
        {
            const std::uint64_t a = flow.s_vector_size(p);
            if (a >= static_cast<std::uint64_t>(k_ka_nested_min_px) && 10 * a <= 6 * a3 && a > best_a)
            {
                best = p;
                best_a = a;
            }
        }
        run.samples.push_back(ka_sample(flow, best, f.U[3], 1.0, "N3"));
        run.sample_basin.push_back(3);
        run.samples.push_back(ka_sample(flow, f.outlets[5], f.U[5], 1.0, "D5"));
        run.sample_basin.push_back(5);
        run.samples.push_back(ka_sample(flow, f.outlets[8], f.U[8], 2.0, "A8"));
        run.sample_basin.push_back(8);
        run.samples.push_back(ka_sample(flow, -1, 1e-4, 1.0, "OUT"));
        run.sample_basin.push_back(static_cast<std::size_t>(-1));
    }
    lsdtt::lsdtt_erosion_rate_basin_options bo;
    bo.snapping_method = "none";
    bo.channel_threshold_pixels = k_ka_channel_px;
    run.basins = lsdtt::delineate_erosion_rate_basins(flow, run.samples, bo, with_units ? &f.units : nullptr);

    lsdtt::lsdtt_chi_options co;
    co.m_over_n = 0.45;
    co.A_0 = 1.0;
    const lsdtt::lsdtt_chi_collection_result chi =
        lsdtt::get_chi_from_outlets_zero_base_threshold(flow, run.basins.outlets, k_ka_channel_px, co);
    for (const auto& b : chi.basins) run.fits.basins.push_back(lsdtt::fit_chi_network(b, f.elevation, {}));

    lsdtt::lsdtt_erodibility_calibration_options opts;
    opts.n_values = n_values;
    std::mt19937_64 rng = lsdtt::erodibility_calibration_bootstrap_engine(5489, 0.45);
    run.ksn = lsdtt::compute_erosion_rate_basin_ksn(run.fits, f.surface.get_meta(), {}, opts, rng);
    lsdtt::lsdtt_chi_lithology_tables tables;
    if (with_units)
    {
        const std::vector<float> cd = lsdtt::compute_lithology_contact_distance(f.surface, f.units, 1);
        lsdtt::lsdtt_chi_lithology_options lo;
        tables = lsdtt::build_chi_lithology_tables(run.fits, f.surface.get_meta(), f.units, cd, {}, lo);
    }
    run.cal = lsdtt::calibrate_erodibility_from_erosion_rates(
        run.samples, run.basins, run.ksn, run.fits, with_units ? &tables : nullptr,
        with_units ? &f.units : nullptr, 0.45, 1.0, opts, rng);
    return run;
}

const lsdtt::lsdtt_erosion_rate_basin_ksn* ka_ksn_of(const ka_run& run, std::size_t sample)
{
    for (const auto& k : run.ksn)
        if (k.outlet_flat == run.basins.basins[sample].outlet_flat) return &k;
    return nullptr;
}

const lsdtt::lsdtt_erodibility_fit_row* ka_fit(const ka_run& run, const std::string& type, double n)
{
    for (const auto& f : run.cal.fits)
        if (f.fit_type == type && (type == "fitted_n" || f.n == n)) return &f;
    return nullptr;
}

void test_known_answer_f1(double n_true)
{
    const std::string tag = "F1 (n = " + fmt(n_true) + "): ";
    const ka_fixture f = make_ka_fixture(n_true, false, false);
    const double wrong_n = n_true == 1.0 ? 2.0 : 1.0;
    const ka_run run = run_ka(f, {n_true, wrong_n}, false);
    std::cout << tag << f.n_edge_outlets << " complete non-overlapping sub-basins of " << k_ka_min_px
              << "-" << k_ka_max_px << " px; using "
              << f.outlets.size() << "\n";

    // F1-a, F1-b per basin.
    double max_k_err = 0.0, max_K_err = 0.0;
    bool all_basins = true;
    for (std::size_t i = 0; i < run.samples.size(); ++i)
    {
        if (run.samples[i].status != "used") continue;
        const auto* k = ka_ksn_of(run, i);
        const std::size_t b = run.sample_basin[i];
        const double want = std::pow(f.U[b] / k_ka_K, 1.0 / n_true);
        const double err = std::fabs(k->k_sn / want - 1.0);
        max_k_err = std::max(max_k_err, err);
        std::cout << "   " << run.samples[i].sample_id << ": U " << f.U[b] << " k_sn " << k->k_sn
                  << " (want " << want << ", rel err " << err << ") [" << k->ci_low << ", " << k->ci_high
                  << "] reaches " << k->n_reaches_qualifying << "/" << k->n_reaches << "\n";
        for (const auto& r : run.cal.basin_K)
            if (r.sample == static_cast<int>(i) && r.n == n_true)
                max_K_err = std::max(max_K_err, std::fabs(r.K / k_ka_K - 1.0));
    }
    for (std::size_t b = 0; b < f.outlets.size(); ++b)
        all_basins = all_basins && run.samples[b].status == "used";
    const auto* fx = ka_fit(run, "fixed_n", n_true);
    const auto* fw = ka_fit(run, "fixed_n", wrong_n);
    const auto* fn = ka_fit(run, "fitted_n", 0.0);
    const double K_fixed_err = fx ? std::fabs(fx->K / k_ka_K - 1.0) : 1.0;
    const double n_err = fn && fn->status == "ok" ? std::fabs(fn->n - n_true) : 1.0;
    const double K_fit_err = fn && fn->status == "ok" ? std::fabs(fn->K / k_ka_K - 1.0) : 1.0;
    std::cout << std::setprecision(8) << tag << "MEASURED max |k_sn err| " << max_k_err
              << ", max |K_b err| " << max_K_err << ", fixed-n K " << (fx ? fx->K : 0.0) << " (err "
              << K_fixed_err << ", corr_lnK_lnk " << (fx ? fx->corr_lnK_lnk : 0.0) << ", dispersion "
              << (fx ? fx->ln_K_dispersion : 0.0) << ") [" << (fx ? fx->K_ci_low : 0.0) << ", "
              << (fx ? fx->K_ci_high : 0.0) << "]; wrong n = " << wrong_n << ": dispersion "
              << (fw ? fw->ln_K_dispersion : 0.0) << ", corr " << (fw ? fw->corr_lnK_lnk : 0.0)
              << "; fitted n " << (fn ? fn->n : 0.0) << " [" << (fn ? fn->n_ci_low : 0.0) << ", "
              << (fn ? fn->n_ci_high : 0.0) << "] (err " << n_err << ", York SE " << (fn ? fn->n_se_york : 0.0)
              << ", MSWD " << (fn ? fn->mswd : 0.0) << "), K " << (fn ? fn->K : 0.0) << " (err "
              << K_fit_err << ") [" << (fn ? fn->K_ci_low : 0.0) << ", " << (fn ? fn->K_ci_high : 0.0)
              << "], N_groups " << run.cal.n_groups << std::setprecision(6) << "\n";

    // Pins: 2 x the measured error, capped at the plan's caps (k_sn 2 %, K
    // n x 2 %, n 0.03 n, K at the fitted n 6 % / 12 %). Measured 2026-10-19:
    //   n = 1: k_sn 1.853 %, K_b 1.819 %, fixed-n K 1.484 % (corr_lnK_lnk
    //          0.495), n 0.00123, K at fitted n 1.852 %; wrong n = 2:
    //          dispersion 0.951, corr -0.99999765
    //   n = 2: k_sn 1.853 %, K_b 3.605 %, fixed-n K 2.946 % (corr 0.493),
    //          n 0.0049, K at fitted n 3.671 %; wrong n = 1: dispersion 0.478,
    //          corr 0.99996276
    // Every basin k_sn is HIGH by +1.1 to +1.9 %: the steady state integrates
    // chi with a left rectangle rule, chi extraction with the trapezoid rule.
    const bool n1 = n_true == 1.0;
    const double pin_k = 0.02, pin_Kb = n1 ? 0.02 : 0.04, pin_Kfixed = n1 ? 0.02 : 0.04;
    const double pin_n = n1 ? 0.0025 : 0.0099, pin_Kfit = n1 ? 0.037 : 0.074;
    check(tag + "F1-a every basin k_sn = (U_b/K)^(1/n) within " + fmt(pin_k), all_basins && max_k_err <= pin_k,
          fmt(max_k_err));
    check(tag + "F1-b every basin K = K_true within " + fmt(pin_Kb), max_K_err <= pin_Kb, fmt(max_K_err));
    check(tag + "F1-c fixed n = n_true: K within " + fmt(pin_Kfixed), K_fixed_err <= pin_Kfixed, fmt(K_fixed_err));
    check(tag + "F1-c wrong n: |corr_lnK_lnk| > 0.9", fw && std::fabs(fw->corr_lnK_lnk) > 0.9);
    if (n1)
        check(tag + "F1-c wrong n = 2 on the n = 1 surface: ln_K_dispersion > 0.5",
              fw && fw->ln_K_dispersion > 0.5);
    check(tag + "F1-d fitted n within " + fmt(pin_n) + ", interval brackets n",
          fn && fn->status == "ok" && n_err <= pin_n && fn->n_ci_low <= n_true && n_true <= fn->n_ci_high,
          fmt(n_err));
    check(tag + "F1-e K at the fitted n within " + fmt(pin_Kfit), K_fit_err <= pin_Kfit, fmt(K_fit_err));
    check(tag + "fit rows: one fixed_n per n and one fitted_n; MSWD and both York SEs written",
          run.cal.fits.size() == 3 && fn && std::isfinite(fn->mswd) && fn->n_se_york > 0.0 &&
              fn->n_se_york_mswd_scaled > 0.0);

    // Determinism: the same seed gives the same intervals bitwise.
    const ka_run again = run_ka(f, {n_true, wrong_n}, false);
    const auto* fn2 = ka_fit(again, "fitted_n", 0.0);
    check(tag + "same seed, same intervals (bitwise)",
          fn2 && std::memcmp(&fn2->n_ci_low, &fn->n_ci_low, sizeof(double)) == 0 &&
              std::memcmp(&fn2->K_ci_high, &fn->K_ci_high, sizeof(double)) == 0 &&
              again.cal.basin_K.size() == run.cal.basin_K.size() &&
              std::memcmp(&again.cal.basin_K[0].K_ci_low, &run.cal.basin_K[0].K_ci_low, sizeof(double)) == 0);

    // F1-f bookkeeping.
    const auto count = [&](const char* st) { return lsdtt::count_erosion_rate_samples_with_status(run.samples, st); };
    std::cout << tag << "statuses: used " << count("used") << ", area_mismatch " << count("area_mismatch")
              << ", outside_dem " << count("outside_dem") << ", too_few_reaches " << count("too_few_reaches")
              << "\n";
    check(tag + "F1-f status counts: 12 used, 1 area_mismatch, 1 outside_dem",
          count("used") == 12 && count("area_mismatch") == 1 && count("outside_dem") == 1);
    check(tag + "F1-f N_groups = 10", run.cal.n_groups == 10);
    const auto* kn = ka_ksn_of(run, 10);
    const double want3 = std::pow(f.U[3] / k_ka_K, 1.0 / n_true);
    check(tag + "F1-f the nested sample shares basin 3's E and group, k_sn within the cap",
          run.samples.size() > 10 && run.samples[10].E_m_per_yr == run.samples[3].E_m_per_yr &&
              run.basins.basins[10].nesting_group == run.basins.basins[3].nesting_group && kn &&
              std::fabs(kn->k_sn / want3 - 1.0) <= pin_k);
    check(tag + "F1-f the duplicate shares basin 5's outlet and group",
          run.basins.basins[11].outlet_flat == run.basins.basins[5].outlet_flat &&
              run.basins.basins[11].duplicates == std::vector<int>{5});
}

void test_known_answer_f2()
{
    const std::string tag = "F2: ";
    const ka_fixture f = make_ka_fixture(1.0, true, false);
    const ka_run run = run_ka(f, {1.0}, true, false);
    std::size_t n_mono_w = 0, n_mono_e = 0;
    for (std::size_t i = 0; i < run.samples.size(); ++i)
    {
        const auto& b = run.basins.basins[i];
        std::cout << "   " << run.samples[i].sample_id << " status " << run.samples[i].status << " dominant "
                  << b.dominant_code << " " << b.dominant_area_fraction << "\n";
        if (run.samples[i].status == "used" && b.dominant_area_fraction >= 0.8)
            (b.dominant_code == 1 ? n_mono_w : n_mono_e)++;
    }
    double errs[2][2] = {{1.0, 1.0}, {1.0, 1.0}};   // [method mono/uib][unit]
    for (const auto& r : run.cal.units)
    {
        const double want = r.lith_code == 1 ? k_ka_K1 : k_ka_K2;
        const int mi = r.method == "mono" ? 0 : 1;
        const double e = r.status == "ok" ? std::fabs(r.K / want - 1.0) : 1.0;
        errs[mi][r.lith_code - 1] = e;
        std::cout << "   " << r.method << " " << r.lith_name << ": K " << r.K << " [" << r.K_ci_low << ", "
                  << r.K_ci_high << "] (want " << want << ", rel err " << e << ") basins " << r.N_basins
                  << " groups " << r.N_groups << "\n";
    }
    const auto* fn = ka_fit(run, "fitted_n", 0.0);
    std::cout << std::setprecision(8) << tag << "MEASURED unit_in_basin err west " << errs[1][0] << " east "
              << errs[1][1] << "; mono err west " << errs[0][0] << " east " << errs[0][1] << "; mono basins "
              << n_mono_w << " / " << n_mono_e << "; bulk fitted n " << (fn ? fn->n : 0.0) << " ("
              << (fn ? fn->status : "") << ", MSWD " << (fn ? fn->mswd : 0.0) << ")" << std::setprecision(6)
              << "\n";
    // Pins: 2 x measured, capped at n x 2 %. Measured 2026-10-19:
    // unit_in_basin west 1.594 %, east 1.558 %; mono west 1.503 %, east 1.543 %;
    // mono basins 3 / 5; unit_in_basin entries 5 / 7; bulk fitted n 0.6779
    // (MSWD 36.7: the lithology confound, F2-e, not asserted). F2-c: 0.116 %
    // west, 0.030 % east.
    check(tag + "F2-a unit_in_basin recovers K1 and K2 within 0.02", errs[1][0] <= 0.02 && errs[1][1] <= 0.02);
    check(tag + "F2-b mono recovers K1 and K2 within 0.02", errs[0][0] <= 0.02 && errs[0][1] <= 0.02);
    std::size_t uib_w = 0, uib_e = 0;
    for (const auto& r : run.cal.units)
        if (r.method == "unit_in_basin") (r.lith_code == 1 ? uib_w : uib_e) = r.N_basins;
    check(tag + "F2-b mono basins 3 west / 5 east; unit_in_basin entries 5 / 7 (exact)",
          n_mono_w == 3 && n_mono_e == 5 && uib_w == 5 && uib_e == 7);

    // F2-c: with every U_b equal, unit_in_basin agrees with relative K's
    // absolute K (E = U) on the same fits, to the sum of their pins (2 % + 2 %).
    const ka_fixture fe = make_ka_fixture(1.0, true, true);
    const ka_run re = run_ka(fe, {1.0}, true, false);
    const std::vector<float> cd = lsdtt::compute_lithology_contact_distance(fe.surface, fe.units, 1);
    lsdtt::lsdtt_chi_lithology_options lo;
    const auto t = lsdtt::build_chi_lithology_tables(re.fits, fe.surface.get_meta(), fe.units, cd, {}, lo);
    std::mt19937_64 srng(5489);
    const auto sum = lsdtt::summarise_chi_ksn_by_lithology(t, fe.units, 0.45, 1.0, lo, srng);
    lsdtt::lsdtt_relative_erodibility_options ro;
    ro.n_values = {1.0};
    ro.erosion_rate_m_per_yr = k_ka_Uref;
    std::mt19937_64 rrng = lsdtt::lithology_relative_erodibility_bootstrap_engine(5489, 0.45);
    const auto rel = lsdtt::compute_lithology_relative_erodibility(t, sum, ro, rrng);
    double diff[3] = {1.0, 1.0, 1.0};
    for (const auto& u : re.cal.units)
        if (u.method == "unit_in_basin")
            for (const auto& r : rel.rows)
                if (r.lith_code == u.lith_code) diff[u.lith_code] = std::fabs(u.K / r.K_abs - 1.0);
    std::cout << std::setprecision(8) << tag << "MEASURED F2-c |K_unit_in_basin / K_relative - 1| west "
              << diff[1] << " east " << diff[2] << std::setprecision(6) << "\n";
    check(tag + "F2-c equal U: unit_in_basin K = relative-K absolute K within 0.04",
          diff[1] <= 0.04 && diff[2] <= 0.04);
}

// ---------------------------------------------------------------------------
// CLI fixtures (--write-fixtures <dir>)
//
// f1_dem.tif / f2_dem.tif: the F1 (n = 1) and F2 steady states, routed by the
// CLI as given (raster_is_filled: true, "oooo"). f1_samples.csv /
// f2_samples.csv: the samples of run_ka (outlets as WGS84 latitude/longitude
// of their pixel centres; E in mm/kyr with 10 % sigma; reported areas from
// routing the DEM exactly as the CLI does; the off-grid point 1 degree east).
// f1_samples_m_per_yr.csv: f1 in m/yr. f2_units.gpkg: west / east polygons
// (lith_primary). The CSV with a units column is checked by the reader.
// ---------------------------------------------------------------------------
std::vector<pt> ka_lonlat(const std::vector<pt>& xy)
{
    OGRSpatialReference utm, wgs;
    utm.importFromEPSG(32632);
    wgs.importFromEPSG(4326);
    utm.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    wgs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    std::unique_ptr<OGRCoordinateTransformation> t(OGRCreateCoordinateTransformation(&utm, &wgs));
    std::vector<pt> out = xy;
    for (auto& p : out) t->Transform(1, &p.x, &p.y);
    return out;
}

void write_ka_samples(const std::string& path, const ka_fixture& f, const lsdtt::lsdtt_flowinfo& cli_flow,
                      bool m_per_yr)
{
    const lsdtt::lsdtt_flowinfo flow = ka_flow(f.surface);
    struct row { std::string id; std::int64_t flat; double U, area_factor; };
    std::vector<row> rows;
    for (std::size_t b = 0; b < f.outlets.size(); ++b) rows.push_back({"B" + std::to_string(b), f.outlets[b], f.U[b], 1.0});
    // The extras as run_ka builds them.
    const std::uint64_t a3 = flow.s_vector_size(f.outlets[3]);
    std::int64_t best = -1;
    std::uint64_t best_a = 0;
    for (const std::int64_t p : flow.upslope_nodes(f.outlets[3]))
    {
        const std::uint64_t a = flow.s_vector_size(p);
        if (a >= static_cast<std::uint64_t>(k_ka_nested_min_px) && 10 * a <= 6 * a3 && a > best_a) { best = p; best_a = a; }
    }
    rows.push_back({"N3", best, f.U[3], 1.0});
    rows.push_back({"D5", f.outlets[5], f.U[5], 1.0});
    rows.push_back({"A8", f.outlets[8], f.U[8], 2.0});
    const auto& m = f.surface.get_meta();
    std::vector<pt> xy;
    for (const auto& r : rows)
        xy.push_back({m.x_min + (flow.col_from_flat(r.flat) + 0.5) * m.cell_size,
                      m.y_min + (m.n_rows - flow.row_from_flat(r.flat) - 0.5) * m.cell_size});
    xy.push_back({m.x_min + m.n_cols * m.cell_size + 80000.0, m.y_min + 100.0});
    const std::vector<pt> ll = ka_lonlat(xy);
    std::ofstream out(path);
    out << std::setprecision(15);
    out << "sample_id,latitude,longitude,erosion_rate,erosion_rate_uncertainty,erosion_rate_units,reported_area_km2\n";
    const double scale = m_per_yr ? 1.0 : 1e6;
    const std::string units = m_per_yr ? "m_per_yr" : "mm_per_kyr";
    for (std::size_t i = 0; i < rows.size(); ++i)
    {
        const double area_km2 = rows[i].area_factor * static_cast<double>(cli_flow.s_vector_size(rows[i].flat)) *
                                m.cell_size * m.cell_size / 1e6;
        out << rows[i].id << "," << ll[i].y << "," << ll[i].x << "," << rows[i].U * scale << ","
            << 0.1 * rows[i].U * scale << "," << units << "," << area_km2 << "\n";
    }
    out << "OUT," << ll.back().y << "," << ll.back().x << "," << 1e-4 * scale << "," << 1e-5 * scale << ","
        << units << ",5\n";
}

void write_ka_units(const std::string& path, const lsdtt::lsdtt_raster_meta& m)
{
    std::filesystem::remove(path);
    GDALDataset* ds = GetGDALDriverManager()->GetDriverByName("GPKG")->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    OGRSpatialReference srs;
    srs.importFromEPSG(32632);
    srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    OGRLayer* layer = ds->CreateLayer("units", &srs, wkbPolygon, nullptr);
    OGRFieldDefn fd("lith_primary", OFTString);
    (void)layer->CreateField(&fd);
    const double xs = m.x_min + k_ka_split * m.cell_size, x1 = m.x_min + m.n_cols * m.cell_size;
    const double y0 = m.y_min, y1 = m.y_min + m.n_rows * m.cell_size;
    for (int k = 0; k < 2; ++k)
    {
        OGRFeature* f = OGRFeature::CreateFeature(layer->GetLayerDefn());
        f->SetField("lith_primary", k == 0 ? "west" : "east");
        const double a = k == 0 ? m.x_min : xs, b = k == 0 ? xs : x1;
        OGRLinearRing ring;
        ring.addPoint(a, y0); ring.addPoint(b, y0); ring.addPoint(b, y1); ring.addPoint(a, y1); ring.addPoint(a, y0);
        OGRPolygon poly;
        poly.addRing(&ring);
        f->SetGeometry(&poly);
        if (layer->CreateFeature(f) != OGRERR_NONE) std::cerr << "could not write a unit polygon\n";
        OGRFeature::DestroyFeature(f);
    }
    GDALClose(ds);
}

int write_fixtures(const std::string& dir)
{
    std::filesystem::create_directories(dir);
    lsdtt::lsdtt_data_loader w;
    const ka_fixture f1 = make_ka_fixture(1.0, false, false);
    w.write_full(f1.elevation.get_data(), f1.elevation.get_meta(), dir + "/f1_dem.tif",
                 lsdtt::raster_output_format::gtiff);
    const lsdtt::lsdtt_flowinfo cli1 = ka_flow(f1.elevation);
    write_ka_samples(dir + "/f1_samples.csv", f1, cli1, false);
    write_ka_samples(dir + "/f1_samples_m_per_yr.csv", f1, cli1, true);
    const ka_fixture f2 = make_ka_fixture(1.0, true, false);
    w.write_full(f2.elevation.get_data(), f2.elevation.get_meta(), dir + "/f2_dem.tif",
                 lsdtt::raster_output_format::gtiff);
    const lsdtt::lsdtt_flowinfo cli2 = ka_flow(f2.elevation);
    write_ka_samples(dir + "/f2_samples.csv", f2, cli2, false);
    write_ka_units(dir + "/f2_units.gpkg", f2.surface.get_meta());
    std::cout << "wrote the erodibility calibration CLI fixtures to " << dir << "\n";
    return 0;
}

// --check-erodibility-calibration-fit <fit.csv> <fit_type> <n> <expected_n> <tol_n>
//     <expected_K> <tol_K>: a written _erodibility_calibration_fit_ CSV against
// expected values (n absolute, K relative). For fixed_n the row with that n is
// read and expected_n is ignored; status must be ok.
int check_fit(int argc, char* argv[])
{
    if (argc != 9)
    {
        std::cerr << "usage: --check-erodibility-calibration-fit <csv> <fixed_n|fitted_n> <n> <expected_n> "
                     "<tol_n> <expected_K> <tol_K>\n";
        return 2;
    }
    try
    {
        const std::string csv = argv[2], type = argv[3];
        const double n_sel = std::stod(argv[4]), want_n = std::stod(argv[5]), tol_n = std::stod(argv[6]);
        const double want_K = std::stod(argv[7]), tol_K = std::stod(argv[8]);
        const lsdtt::lsdtt_csv_table t = lsdtt::load_csv_table(csv);
        for (std::size_t i = 0; i < t.rows.size(); ++i)
        {
            const auto& r = t.rows[i];
            const int line = static_cast<int>(i) + 2;
            if (lsdtt::csv_row_string(r, "fit_type") != type) continue;
            if (type == "fixed_n" && lsdtt::csv_row_double_required(r, "n", csv, line) != n_sel) continue;
            const std::string status = lsdtt::csv_row_string(r, "status");
            const double n = lsdtt::csv_row_double_required(r, "n", csv, line);
            const double K = lsdtt::csv_row_double_required(r, "K", csv, line);
            const double errK = std::fabs(K / want_K - 1.0);
            const double errn = type == "fitted_n" ? std::fabs(n - want_n) : 0.0;
            std::cout << std::setprecision(10) << type << " (" << status << "): n " << n << " (abs err " << errn
                      << "), K " << K << " (rel err " << errK << "), interval [" << lsdtt::csv_row_string(r, "K_ci_low")
                      << ", " << lsdtt::csv_row_string(r, "K_ci_high") << "], n interval ["
                      << lsdtt::csv_row_string(r, "n_ci_low") << ", " << lsdtt::csv_row_string(r, "n_ci_high")
                      << "], MSWD " << lsdtt::csv_row_string(r, "mswd") << ", N_groups "
                      << lsdtt::csv_row_string(r, "N_groups") << "\n";
            const bool ok = status == "ok" && errK <= tol_K && errn <= tol_n;
            std::cout << (ok ? "PASS" : "FAIL") << "\n";
            return ok ? 0 : 1;
        }
        std::cerr << "FAIL: no " << type << " row in " << csv << "\n";
        return 1;
    }
    catch (const std::exception& e)
    {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}

// --check-erodibility-calibration-unit <by_lithology.csv> <method> <n> <unit> <expected_K> <tol_K>
int check_unit(int argc, char* argv[])
{
    if (argc != 8)
    {
        std::cerr << "usage: --check-erodibility-calibration-unit <csv> <method> <n> <unit> <expected_K> <tol_K>\n";
        return 2;
    }
    try
    {
        const std::string csv = argv[2], method = argv[3], unit = argv[5];
        const double n = std::stod(argv[4]), want = std::stod(argv[6]), tol = std::stod(argv[7]);
        const lsdtt::lsdtt_csv_table t = lsdtt::load_csv_table(csv);
        for (std::size_t i = 0; i < t.rows.size(); ++i)
        {
            const auto& r = t.rows[i];
            const int line = static_cast<int>(i) + 2;
            if (lsdtt::csv_row_string(r, "method") != method || lsdtt::csv_row_string(r, "lith_name") != unit ||
                lsdtt::csv_row_double_required(r, "n", csv, line) != n)
                continue;
            const double K = lsdtt::csv_row_double_required(r, "K", csv, line);
            const double err = std::fabs(K / want - 1.0);
            std::cout << std::setprecision(10) << method << " " << unit << " n = " << n << ": K " << K
                      << " (expected " << want << ", rel err " << err << "), basins "
                      << lsdtt::csv_row_string(r, "N_basins") << ", status " << lsdtt::csv_row_string(r, "status") << "\n";
            const bool ok = lsdtt::csv_row_string(r, "status") == "ok" && err <= tol;
            std::cout << (ok ? "PASS" : "FAIL") << "\n";
            return ok ? 0 : 1;
        }
        std::cerr << "FAIL: no " << method << " row for " << unit << " in " << csv << "\n";
        return 1;
    }
    catch (const std::exception& e)
    {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}

} // namespace

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr << "usage: lsdtt-test-erodibility-calibration <work_dir>\n";
        return 2;
    }
    lsdtt::gdal_init();
    const std::string mode = argv[1];
    if (mode == "--write-fixtures" && argc == 3) return write_fixtures(argv[2]);
    if (mode == "--check-erodibility-calibration-fit") return check_fit(argc, argv);
    if (mode == "--check-erodibility-calibration-unit") return check_unit(argc, argv);
    const std::string work = argv[1];
    std::cout << "lsdtt3 erodibility calibration tests (York, reader, seams, parameters)\n";
    test_york();
    test_seams(work);
    try { test_reader(work); }
    catch (const std::exception& e) { check("reader tests ran without an exception", false, e.what()); }
    test_parameters();
    test_edge_completeness();
    test_known_answer_f1(1.0);
    test_known_answer_f1(2.0);
    test_known_answer_f2();
    std::cout << g_n_checks - g_n_failed << " / " << g_n_checks << " checks passed\n"
              << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";
    return g_n_failed == 0 ? 0 : 1;
}
