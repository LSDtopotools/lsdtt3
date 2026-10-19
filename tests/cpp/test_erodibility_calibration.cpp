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
// Usage: lsdtt-test-erodibility-calibration <work_dir>
// -----------------------------------------------------------------------------

#include <cmath>
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
#include "../../src/core/lsdtt_gdal_init.hpp"
#include "../../src/core/lsdtt_param_parser.hpp"
#include "../../src/cpu/lsdtt_basin_cpu.hpp"
#include "../../src/cpu/lsdtt_erodibility_calibration_cpu.hpp"
#include "../../src/cpu/lsdtt_flowinfo_cpu.hpp"
#include "../../src/cpu/lsdtt_lithology_relative_erodibility_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"
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

} // namespace

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr << "usage: lsdtt-test-erodibility-calibration <work_dir>\n";
        return 2;
    }
    lsdtt::gdal_init();
    const std::string work = argv[1];
    std::cout << "lsdtt3 erodibility calibration tests (York, reader, seams, parameters)\n";
    test_york();
    test_seams(work);
    try { test_reader(work); }
    catch (const std::exception& e) { check("reader tests ran without an exception", false, e.what()); }
    test_parameters();
    test_edge_completeness();
    std::cout << g_n_checks - g_n_failed << " / " << g_n_checks << " checks passed\n"
              << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";
    return g_n_failed == 0 ? 0 : 1;
}
