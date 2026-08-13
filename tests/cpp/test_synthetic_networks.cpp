// test_synthetic_networks.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <ogr_spatialref.h>

#include "lsdtt_flow_base.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_meta.hpp"
#include "lsdtt_synthetic_network_cpu.hpp"

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

void check_case(const std::string& case_name,
                const std::string& name,
                bool               ok,
                const std::string& detail = "")
{
    check(case_name + " :: " + name, ok, detail);
}

struct synthetic_case
{
    std::string              name;
    lsdtt::lsdtt_ocn_options opts;
};

std::vector<std::uint8_t> build_island_mask(int n_rows, int n_cols)
{
    std::vector<std::uint8_t> mask(static_cast<std::size_t>(n_rows) *
                                       static_cast<std::size_t>(n_cols),
                                   static_cast<std::uint8_t>(0));
    for (int r = 1; r < n_rows - 1; ++r)
    {
        for (int c = 1; c < n_cols - 1; ++c)
        {
            mask[static_cast<std::size_t>(r) * static_cast<std::size_t>(n_cols) +
                 static_cast<std::size_t>(c)] = static_cast<std::uint8_t>(1);
        }
    }

    // Carve small coastal embayments for a simple island-like shape.
    mask[1 * n_cols + 1]                 = 0;
    mask[1 * n_cols + (n_cols - 2)]      = 0;
    mask[(n_rows - 2) * n_cols + 2]      = 0;
    mask[(n_rows - 2) * n_cols + n_cols - 3] = 0;
    return mask;
}

void run_case(const synthetic_case& sc)
{
    std::cout << "running case: " << sc.name << "\n";

    lsdtt::lsdtt_synthetic_network_cpu engine(sc.opts);

    const double e0 = engine.initial_energy();
    check_case(sc.name,
               "initial energy finite",
               std::isfinite(e0) && e0 > 0.0,
               "initial energy must be > 0");

    engine.run();

    const double e1 = engine.final_energy();
    check_case(sc.name,
               "final energy finite",
               std::isfinite(e1) && e1 > 0.0,
               "final energy must be > 0");
    if (sc.opts.temperature_scale == 0.0 && sc.opts.initial_temperature <= 0.0)
    {
        check_case(sc.name,
                   "energy does not increase",
                   e1 <= e0 + 1e-12 * std::abs(e0),
                   "anneal should not increase energy in this deterministic setup");
    }

    const auto dem = engine.get_elevation_raster();
    const auto area = engine.get_drainage_area_cells_raster();
    const std::size_t expected = static_cast<std::size_t>(sc.opts.n_rows) *
                                 static_cast<std::size_t>(sc.opts.n_cols);
    const int n_cols = sc.opts.n_cols;

    check_case(sc.name,
               "dem size matches grid",
               dem.data.size() == expected && dem.mask.size() == expected,
               "DEM vectors should be n_rows*n_cols");
    check_case(sc.name,
               "area size matches grid",
               area.data.size() == expected && area.mask.size() == expected,
               "area vectors should be n_rows*n_cols");

    const auto& areas = engine.areas();
    const auto& rec = engine.receivers();

    bool sea_receiver_minus_one = true;
    bool outlet_maps_to_self = true;
    bool non_outlet_land_not_to_sea = true;
    bool d8_receiver_ok = true;
    bool receiver_walk_terminates = true;
    bool area_sane = true;
    bool downstream_dem_strictly_decreases = true;
    bool outlet_dem_at_baselevel = true;

    std::size_t n_land = 0;
    double outlet_area_sum = 0.0;
    const double baselevel = static_cast<double>(sc.opts.baselevel_elevation);
    const double baselevel_tolerance =
        1e-9 * std::max(1.0, std::abs(baselevel));

    for (std::size_t i = 0; i < expected; ++i)
    {
        const auto flat = static_cast<std::int64_t>(i);
        if (!engine.is_land(flat))
        {
            if (rec[i] != -1) sea_receiver_minus_one = false;
            continue;
        }

        ++n_land;

        const double a = areas[i];
        const double nearest_int = std::nearbyint(a);
        if (!std::isfinite(a) || a < 1.0 || std::abs(a - nearest_int) > 1e-9)
            area_sane = false;

        if (engine.is_outlet(flat))
        {
            if (rec[i] != flat) outlet_maps_to_self = false;
            outlet_area_sum += a;

            if (dem.mask[i] != 0)
            {
                const double z_outlet = static_cast<double>(dem.data[i]);
                if (!std::isfinite(z_outlet) ||
                    std::abs(z_outlet - baselevel) > baselevel_tolerance)
                {
                    outlet_dem_at_baselevel = false;
                }
            }
            continue;
        }

        const std::int64_t r = rec[i];
        if (r < 0 || r >= static_cast<std::int64_t>(expected) || !engine.is_land(r))
        {
            non_outlet_land_not_to_sea = false;
            d8_receiver_ok = false;
            downstream_dem_strictly_decreases = false;
            continue;
        }

        if (dem.mask[i] != 0 && dem.mask[static_cast<std::size_t>(r)] != 0)
        {
            const double z_i = static_cast<double>(dem.data[i]);
            const double z_r = static_cast<double>(dem.data[static_cast<std::size_t>(r)]);
            if (!std::isfinite(z_i) || !std::isfinite(z_r) || !(z_i > z_r))
                downstream_dem_strictly_decreases = false;
        }

        const int row = static_cast<int>(flat) / n_cols;
        const int col = static_cast<int>(flat) % n_cols;
        const int r_row = static_cast<int>(r) / n_cols;
        const int r_col = static_cast<int>(r) % n_cols;
        const int dr = std::abs(r_row - row);
        const int dc = std::abs(r_col - col);
        if (dr > 1 || dc > 1 || (dr == 0 && dc == 0)) d8_receiver_ok = false;
    }

    for (std::size_t i = 0; i < expected; ++i)
    {
        const auto start = static_cast<std::int64_t>(i);
        if (!engine.is_land(start)) continue;

        std::int64_t cur = start;
        std::size_t  steps = 0;
        while (true)
        {
            const std::int64_t r = rec[static_cast<std::size_t>(cur)];
            if (r < 0 || r >= static_cast<std::int64_t>(expected) || !engine.is_land(r))
            {
                receiver_walk_terminates = false;
                break;
            }
            if (r == cur)
            {
                if (!engine.is_outlet(cur)) receiver_walk_terminates = false;
                break;
            }
            ++steps;
            if (steps > expected)
            {
                receiver_walk_terminates = false;
                break;
            }
            cur = r;
        }
        if (!receiver_walk_terminates) break;
    }

    const double area_sum_tolerance = 1e-9 * std::max<std::size_t>(n_land, 1);
    const bool outlet_area_partitions_land =
        std::abs(outlet_area_sum - static_cast<double>(n_land)) <= area_sum_tolerance;

    check_case(sc.name,
               "sea receiver is -1",
               sea_receiver_minus_one,
               "all sea pixels should have receiver -1");
    check_case(sc.name,
               "outlet maps to self",
               outlet_maps_to_self,
               "all outlets should route to themselves");
    check_case(sc.name,
               "non-outlet land does not drain to sea",
               non_outlet_land_not_to_sea,
               "non-outlet land receivers must be land cells");
    check_case(sc.name,
               "D8 receiver neighbor constraint",
               d8_receiver_ok,
               "non-outlet land receivers must be one-cell D8 neighbors");
    check_case(sc.name,
               "receiver forest walk terminates",
               receiver_walk_terminates,
               "land receiver walks should terminate at an outlet self-loop");
    check_case(sc.name,
               "land area basic sanity",
               area_sane,
               "land drainage area should be finite, >= 1, and near-integer");
    check_case(sc.name,
               "outlet area sums to land count",
               outlet_area_partitions_land,
               "sum(area[outlets]) should equal total number of land cells");
    check_case(sc.name,
               "DEM decreases strictly downstream",
               downstream_dem_strictly_decreases,
               "for non-outlet land pixels, DEM must be > DEM[receiver]");
    check_case(sc.name,
               "outlet DEM equals baselevel",
               outlet_dem_at_baselevel,
               "outlet DEM should equal baselevel_elevation");
}

void run_determinism_checks(const std::string& case_name,
                            lsdtt::lsdtt_ocn_options opts)
{
    std::cout << "running case: " << case_name << "\n";

    // Same-seed repeatability check.
    lsdtt::lsdtt_synthetic_network_cpu run1(opts);
    lsdtt::lsdtt_synthetic_network_cpu run2(opts);
    run1.run();
    run2.run();

    const auto& rec1 = run1.receivers();
    const auto& rec2 = run2.receivers();
    check_case(case_name,
               "same-seed receivers exactly equal",
               rec1 == rec2,
               "receiver arrays should match exactly for identical options and seed");

    const auto& a1 = run1.areas();
    const auto& a2 = run2.areas();
    bool areas_same = (a1.size() == a2.size());
    constexpr double area_tol = 1e-12;
    if (areas_same)
    {
        for (std::size_t i = 0; i < a1.size(); ++i)
        {
            const double lhs = static_cast<double>(a1[i]);
            const double rhs = static_cast<double>(a2[i]);
            const double tol = area_tol * std::max(1.0, std::max(std::abs(lhs), std::abs(rhs)));
            if (!std::isfinite(lhs) || !std::isfinite(rhs) || std::abs(lhs - rhs) > tol)
            {
                areas_same = false;
                break;
            }
        }
    }
    check_case(case_name,
               "same-seed areas equal (tiny tolerance)",
               areas_same,
               "drainage area arrays should be identical up to tiny floating-point tolerance");

    const auto dem1 = run1.get_elevation_raster();
    const auto dem2 = run2.get_elevation_raster();
    bool dem_same_valid = dem1.mask.size() == dem2.mask.size() &&
                          dem1.data.size() == dem2.data.size() &&
                          dem1.mask == dem2.mask;
    if (dem_same_valid)
    {
        for (std::size_t i = 0; i < dem1.data.size(); ++i)
        {
            if (dem1.mask[i] == 0) continue;
            if (dem1.data[i] != dem2.data[i])
            {
                dem_same_valid = false;
                break;
            }
        }
    }
    check_case(case_name,
               "same-seed DEM data equal on valid cells",
               dem_same_valid,
               "DEM masks should match and valid-cell DEM values should be exactly equal");

    // Different-seed diversity check.
    opts.seed = 1;
    lsdtt::lsdtt_synthetic_network_cpu seed1(opts);
    seed1.run();

    opts.seed = 2;
    lsdtt::lsdtt_synthetic_network_cpu seed2(opts);
    seed2.run();

    const auto& r1 = seed1.receivers();
    const auto& r2 = seed2.receivers();
    const std::size_t n = std::min(r1.size(), r2.size());

    std::size_t land_pixels = 0;
    std::size_t different_receivers = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
        const auto idx = static_cast<std::int64_t>(i);
        if (!seed1.is_land(idx)) continue;
        ++land_pixels;
        if (r1[i] != r2[i]) ++different_receivers;
    }

    const double diff_fraction =
        (land_pixels > 0)
            ? static_cast<double>(different_receivers) / static_cast<double>(land_pixels)
            : 0.0;
    const bool diff_seed_receivers_sufficiently_different = diff_fraction >= 0.01;

    std::ostringstream diff_detail;
    diff_detail << "different receivers=" << different_receivers
                << "/" << land_pixels
                << " (fraction=" << diff_fraction << ")";

    check_case(case_name,
               "different seeds change >=1% of land receivers",
               diff_seed_receivers_sufficiently_different,
               diff_detail.str());
}

// -----------------------------------------------------------------------------
// D8 round-trip routing fidelity
//
// The OCN is generated as a receiver graph and the DEM is then synthesised from
// it via z = baselevel + f*chi. Everything downstream of the generator assumes
// the DEM still CARRIES the generated network. The existing "DEM decreases
// strictly downstream" invariant is necessary but NOT sufficient for that: a
// pixel can be strictly higher than its OCN receiver while a different
// neighbour is lower still, in which case independent D8 routing picks a
// different receiver and the DEM no longer carries the OCN.
//
// This mirrors the CLI's "D8 round-trip fidelity (interior)" computation
// (cli/lsdtt-synthetic-networks.cpp) in process, so the number is asserted on
// rather than printed and discarded.
// -----------------------------------------------------------------------------
struct fidelity_measurement
{
    double    fidelity = -1.0;
    long long interior = 0;
    long long matches  = 0;
};

fidelity_measurement measure_d8_fidelity(lsdtt::lsdtt_synthetic_network_cpu& engine)
{
    fidelity_measurement out;

    lsdtt::lsdtt_raster_result elev = engine.get_elevation_raster();
    lsdtt::lsdtt_raster_cpu    dem(std::move(elev));

    lsdtt::lsdtt_flowinfo_options fi_opts;
    auto flow = lsdtt::make_flow_base(dem.get_data_array(),
                                      dem.get_mask_array(),
                                      dem.get_meta(),
                                      fi_opts);

    const auto& ocn_recv = engine.receivers();
    const int   nr       = engine.n_rows();
    const int   nc       = engine.n_cols();

    for (int r = 1; r < nr - 1; ++r)
    {
        for (int c = 1; c < nc - 1; ++c)
        {
            const std::int64_t f = static_cast<std::int64_t>(r) * nc + c;
            if (!engine.is_land(f))     continue;
            if (engine.is_outlet(f))    continue;
            if (!flow->is_valid_cell(f)) continue;

            ++out.interior;
            if (flow->receiver_of(f) == ocn_recv[static_cast<std::size_t>(f)])
                ++out.matches;
        }
    }

    if (out.interior > 0)
        out.fidelity = static_cast<double>(out.matches) /
                       static_cast<double>(out.interior);
    return out;
}

void run_routing_fidelity_checks(const std::string& case_name)
{
    std::cout << "running case: " << case_name << "\n";

    // MEASURED 2026-08-13, 48x48 / cell 50 m / gamma 0.5 / seed 42 and the
    // production default schedule (n_iter = -1 -> 40*n_rows*n_cols = 92160):
    //
    //   temperature_scale   fidelity   matches/interior
    //   0                   0.916352   1939/2116
    //   0.01                0.879962   1862/2116
    //   0.1                 0.856805   1813/2116
    //   0.3                 0.873819   1849/2116
    //   1                   0.845936   1790/2116
    //
    // Two things this records. First, fidelity is NOT 1 even for a fully
    // relaxed network: a greedy local minimum of sum(A^gamma) under single-flip
    // moves is not identically the steepest-descent graph of the chi surface,
    // so ~8% of interior pixels re-route under independent D8. Second, it is
    // strongly dependent on how far the anneal has actually run - the same
    // configuration truncated to n_iter = 20000 gives only 0.703 at
    // temperature_scale = 0. The floor below is pinned ~1.8% under the measured
    // relaxed value; the second check pins the existence of the temperature
    // limitation without pinning a fragile number.
    constexpr double fidelity_floor_relaxed = 0.90;  // measured 0.916352

    lsdtt::lsdtt_ocn_options opts;
    opts.n_rows            = 48;
    opts.n_cols            = 48;
    opts.cell_size         = 50.0;
    opts.gamma             = 0.5;
    opts.seed              = 42;
    opts.n_iter            = -1;   // OCNet default schedule, 40*n_rows*n_cols
    opts.outlet_mode       = lsdtt::ocn_outlet_mode::single_south;
    opts.temperature_scale = 0.0;

    lsdtt::lsdtt_synthetic_network_cpu cold(opts);
    cold.run();
    const fidelity_measurement f_cold = measure_d8_fidelity(cold);

    opts.temperature_scale = 1.0;
    lsdtt::lsdtt_synthetic_network_cpu hot(opts);
    hot.run();
    const fidelity_measurement f_hot = measure_d8_fidelity(hot);

    std::cout << std::fixed << std::setprecision(6)
              << "  fidelity(temperature_scale=0) = " << f_cold.fidelity
              << "  (" << f_cold.matches << "/" << f_cold.interior << ")\n"
              << "  fidelity(temperature_scale=1) = " << f_hot.fidelity
              << "  (" << f_hot.matches << "/" << f_hot.interior << ")\n"
              << std::defaultfloat;

    check_case(case_name,
               "interior sample is non-trivial",
               f_cold.interior > 1000 && f_hot.interior > 1000,
               "48x48 interior should supply >1000 comparable pixels");

    std::ostringstream cold_detail;
    cold_detail << "fidelity=" << f_cold.fidelity
                << " floor=" << fidelity_floor_relaxed
                << " (" << f_cold.matches << "/" << f_cold.interior << ")";
    check_case(case_name,
               "relaxed network D8 fidelity above pinned floor",
               f_cold.fidelity >= fidelity_floor_relaxed,
               cold_detail.str());

    std::ostringstream mono_detail;
    mono_detail << "fidelity(T=0)=" << f_cold.fidelity
                << " fidelity(T=1)=" << f_hot.fidelity;
    check_case(case_name,
               "fidelity degrades with annealing temperature",
               f_hot.fidelity < f_cold.fidelity,
               mono_detail.str());
}

// -----------------------------------------------------------------------------
// Georeferencing (lsdtt_raster_meta::set_synthetic_georeference and friends)
// -----------------------------------------------------------------------------
lsdtt::lsdtt_raster_meta make_bare_meta(int n_rows, int n_cols, double cell_size)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows        = n_rows;
    meta.n_cols        = n_cols;
    meta.cell_size     = cell_size;
    meta.no_data_value = -9999.0;
    return meta;
}

// Projects (lat, lon) from WGS84 into the supplied EPSG. Returns false if the
// transform could not be built (in which case the caller should skip).
bool project_point(double lat_deg, double lon_deg, int epsg,
                   double& easting, double& northing)
{
    OGRSpatialReference wgs84;
    OGRSpatialReference target;
    if (wgs84.importFromEPSG(4326) != OGRERR_NONE) return false;
    if (target.importFromEPSG(epsg) != OGRERR_NONE) return false;
    wgs84.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    target.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    OGRCoordinateTransformation* ct =
        OGRCreateCoordinateTransformation(&wgs84, &target);
    if (!ct) return false;
    double x = lon_deg;
    double y = lat_deg;
    const bool ok = ct->Transform(1, &x, &y) != 0;
    OGRCoordinateTransformation::DestroyCT(ct);
    if (!ok) return false;
    easting  = x;
    northing = y;
    return true;
}

bool unproject_point(double easting, double northing, int epsg,
                     double& lat_deg, double& lon_deg)
{
    OGRSpatialReference wgs84;
    OGRSpatialReference source;
    if (wgs84.importFromEPSG(4326) != OGRERR_NONE) return false;
    if (source.importFromEPSG(epsg) != OGRERR_NONE) return false;
    wgs84.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    source.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    OGRCoordinateTransformation* ct =
        OGRCreateCoordinateTransformation(&source, &wgs84);
    if (!ct) return false;
    double x = easting;
    double y = northing;
    const bool ok = ct->Transform(1, &x, &y) != 0;
    OGRCoordinateTransformation::DestroyCT(ct);
    if (!ok) return false;
    lon_deg = x;
    lat_deg = y;
    return true;
}

int epsg_of_wkt(const std::string& wkt)
{
    if (wkt.empty()) return 0;
    OGRSpatialReference srs;
    if (srs.importFromWkt(wkt.c_str()) != OGRERR_NONE) return 0;
    if (srs.AutoIdentifyEPSG() != OGRERR_NONE) return 0;
    const char* code = srs.GetAuthorityCode(nullptr);
    return code ? std::atoi(code) : 0;
}

void run_georeference_checks(const std::string& case_name)
{
    std::cout << "running case: " << case_name << "\n";

    // -- G1: default anchor regression pin -----------------------------------
    // Teanu Island / "Watopia" (11.633 S, 166.967 E) in EPSG:32758.
    // The expected x_min/y_min are derived independently here (project the
    // anchor with OGR, then offset by half the grid) so the pin catches a
    // change of anchor, CRS, or centring convention rather than restating the
    // implementation's arithmetic.
    {
        constexpr double watopia_lat = -11.633;
        constexpr double watopia_lon = 166.967;
        constexpr int    watopia_epsg = 32758;
        const int    n_rows    = 32;
        const int    n_cols    = 40;
        const double cell_size = 30.0;

        lsdtt::lsdtt_raster_meta meta = make_bare_meta(n_rows, n_cols, cell_size);
        meta.set_default_synthetic_georeference();

        double e = 0.0, n = 0.0;
        const bool projected = project_point(watopia_lat, watopia_lon,
                                             watopia_epsg, e, n);
        check_case(case_name,
                   "G1 anchor projects with OGR",
                   projected,
                   "could not build WGS84 -> EPSG:32758 transform");

        if (projected)
        {
            const double expect_x = e - 0.5 * n_cols * cell_size;
            const double expect_y = n - 0.5 * n_rows * cell_size;
            std::ostringstream d;
            d << std::setprecision(12)
              << "x_min=" << meta.x_min << " expected=" << expect_x
              << " y_min=" << meta.y_min << " expected=" << expect_y;
            check_case(case_name,
                       "G1 default georeference x_min/y_min pinned",
                       std::abs(meta.x_min - expect_x) < 1e-3 &&
                           std::abs(meta.y_min - expect_y) < 1e-3,
                       d.str());
        }

        std::ostringstream ed;
        ed << "epsg_code=" << meta.epsg_code
           << " wkt epsg=" << epsg_of_wkt(meta.crs_wkt);
        check_case(case_name,
                   "G1 default georeference CRS is EPSG:32758",
                   meta.epsg_code == 32758 &&
                       epsg_of_wkt(meta.crs_wkt) == 32758 &&
                       !meta.crs_wkt.empty(),
                   ed.str());
        check_case(case_name,
                   "G1 default georeference is projected metres",
                   !meta.is_geographic() && meta.is_in_metres(),
                   "default synthetic CRS must be projected and in metres");
    }

    // -- G2: epsg = 0 derives the WGS84 UTM zone from longitude --------------
    {
        struct utm_case { double lat; double lon; int expect_epsg; const char* what; };
        const std::vector<utm_case> utm_cases = {
            { 55.95,  -3.19, 32630, "Edinburgh (zone 30N)" },
            { -33.87, 151.21, 32756, "Sydney (zone 56S)" },
            {  0.0,   -180.0, 32601, "antimeridian west edge (zone 1N)" },
            { -0.5,    179.99, 32760, "antimeridian east edge (zone 60S)" },
            {  10.0,   -0.001, 32630, "just west of Greenwich (zone 30N)" },
            {  10.0,    0.001, 32631, "just east of Greenwich (zone 31N)" },
            { -10.0,    5.999, 32731, "just inside zone 31S" },
            { -10.0,    6.001, 32732, "just inside zone 32S" },
            {  0.0,     0.0,   32631, "equator on the zone boundary -> north" },
        };

        for (const auto& uc : utm_cases)
        {
            lsdtt::lsdtt_raster_meta meta = make_bare_meta(16, 16, 100.0);
            bool threw = false;
            try
            {
                meta.set_synthetic_georeference(uc.lat, uc.lon, 0);
            }
            catch (const std::exception&)
            {
                threw = true;
            }

            std::ostringstream d;
            d << uc.what << ": epsg_code=" << meta.epsg_code
              << " expected=" << uc.expect_epsg
              << (threw ? " (threw)" : "");
            check_case(case_name,
                       std::string("G2 UTM zone derived for ") + uc.what,
                       !threw && meta.epsg_code == uc.expect_epsg,
                       d.str());

            const bool hemisphere_ok =
                (uc.lat >= 0.0) ? (meta.epsg_code / 100 == 326)
                                : (meta.epsg_code / 100 == 327);
            check_case(case_name,
                       std::string("G2 hemisphere prefix for ") + uc.what,
                       !threw && hemisphere_ok,
                       d.str());
        }
    }

    // -- G3: polar latitudes with epsg = 0 -----------------------------------
    {
        const std::vector<double> polar_lats = { 84.5, -84.5, 89.9, -89.9 };
        for (double lat : polar_lats)
        {
            lsdtt::lsdtt_raster_meta meta = make_bare_meta(16, 16, 100.0);
            bool threw_runtime = false;
            try
            {
                meta.set_synthetic_georeference(lat, 10.0, 0);
            }
            catch (const std::runtime_error&) { threw_runtime = true; }
            catch (...) {}

            std::ostringstream d;
            d << "lat=" << lat << " did not throw std::runtime_error";
            check_case(case_name,
                       "G3 polar latitude with epsg=0 throws",
                       threw_runtime,
                       d.str());
        }

        // 84 exactly is still inside the UTM band and must not throw.
        {
            lsdtt::lsdtt_raster_meta meta = make_bare_meta(16, 16, 100.0);
            bool threw = false;
            try { meta.set_synthetic_georeference(84.0, 10.0, 0); }
            catch (const std::exception&) { threw = true; }
            check_case(case_name,
                       "G3 latitude exactly 84 with epsg=0 does not throw",
                       !threw && meta.epsg_code == 32632,
                       "|lat| > 84 is the documented cut, 84 itself is valid");
        }

        // An explicit polar EPSG at the same latitude must be accepted.
        {
            lsdtt::lsdtt_raster_meta meta = make_bare_meta(16, 16, 100.0);
            bool threw = false;
            try { meta.set_synthetic_georeference(89.0, 10.0, 3413); }
            catch (const std::exception&) { threw = true; }
            check_case(case_name,
                       "G3 explicit polar EPSG at 89N does not throw",
                       !threw && meta.epsg_code == 3413 && !meta.crs_wkt.empty(),
                       "explicit projection_epsg must bypass the UTM restriction");
        }
        {
            lsdtt::lsdtt_raster_meta meta = make_bare_meta(16, 16, 100.0);
            bool threw = false;
            try { meta.set_synthetic_georeference(-89.0, 10.0, 3976); }
            catch (const std::exception&) { threw = true; }
            check_case(case_name,
                       "G3 explicit polar EPSG at 89S does not throw",
                       !threw && meta.epsg_code == 3976 && !meta.crs_wkt.empty(),
                       "explicit projection_epsg must bypass the UTM restriction");
        }
    }

    // -- G4: non-finite anchors ----------------------------------------------
    {
        const double nan_v = std::nan("");
        const double inf_v = std::numeric_limits<double>::infinity();
        struct bad_case { double lat; double lon; const char* what; };
        const std::vector<bad_case> bad_cases = {
            { nan_v,  10.0,  "NaN latitude"  },
            { 10.0,   nan_v, "NaN longitude" },
            { inf_v,  10.0,  "+Inf latitude" },
            { 10.0,  -inf_v, "-Inf longitude" },
        };
        for (const auto& bc : bad_cases)
        {
            lsdtt::lsdtt_raster_meta meta = make_bare_meta(16, 16, 100.0);
            bool threw_runtime = false;
            try { meta.set_synthetic_georeference(bc.lat, bc.lon, 0); }
            catch (const std::runtime_error&) { threw_runtime = true; }
            catch (...) {}
            check_case(case_name,
                       std::string("G4 ") + bc.what + " throws",
                       threw_runtime,
                       "non-finite anchors must be rejected");

            // Also rejected when an explicit EPSG is supplied.
            lsdtt::lsdtt_raster_meta meta2 = make_bare_meta(16, 16, 100.0);
            bool threw_runtime2 = false;
            try { meta2.set_synthetic_georeference(bc.lat, bc.lon, 32630); }
            catch (const std::runtime_error&) { threw_runtime2 = true; }
            catch (...) {}
            check_case(case_name,
                       std::string("G4 ") + bc.what + " throws with explicit epsg",
                       threw_runtime2,
                       "the finiteness check must precede EPSG selection");
        }
    }

    // -- G5: the grid is genuinely centred on the anchor and round-trips -----
    {
        struct anchor_case { double lat; double lon; int epsg; const char* what; };
        const std::vector<anchor_case> anchors = {
            { 55.95,  -3.19,  32630, "Edinburgh explicit 32630" },
            { -33.87, 151.21, 0,     "Sydney derived UTM"       },
            { -11.633, 166.967, 32758, "Watopia explicit 32758" },
        };

        for (const auto& ac : anchors)
        {
            const int    n_rows    = 37;   // deliberately odd and unequal
            const int    n_cols    = 53;
            const double cell_size = 25.0;

            lsdtt::lsdtt_raster_meta meta =
                make_bare_meta(n_rows, n_cols, cell_size);
            meta.set_synthetic_georeference(ac.lat, ac.lon, ac.epsg);

            const int used_epsg = meta.epsg_code;
            double e = 0.0, n = 0.0;
            const bool projected = project_point(ac.lat, ac.lon, used_epsg, e, n);

            const double centre_x = meta.x_min + 0.5 * n_cols * cell_size;
            const double centre_y = meta.y_min + 0.5 * n_rows * cell_size;

            std::ostringstream cd;
            cd << std::setprecision(12) << ac.what
               << ": centre=(" << centre_x << "," << centre_y
               << ") projected anchor=(" << e << "," << n << ")";
            check_case(case_name,
                       std::string("G5 grid centre equals projected anchor: ") + ac.what,
                       projected && std::abs(centre_x - e) < 1e-6 &&
                           std::abs(centre_y - n) < 1e-6,
                       cd.str());

            double lat_back = 0.0, lon_back = 0.0;
            const bool back = unproject_point(centre_x, centre_y, used_epsg,
                                              lat_back, lon_back);
            std::ostringstream rd;
            rd << std::setprecision(12) << ac.what
               << ": round trip = (" << lat_back << "," << lon_back
               << ") requested (" << ac.lat << "," << ac.lon << ")";
            check_case(case_name,
                       std::string("G5 centre round-trips to the anchor: ") + ac.what,
                       back && std::abs(lat_back - ac.lat) < 1e-7 &&
                           std::abs(lon_back - ac.lon) < 1e-7,
                       rd.str());

            check_case(case_name,
                       std::string("G5 anchored CRS is projected metres: ") + ac.what,
                       !meta.is_geographic() && meta.is_in_metres(),
                       "anchored synthetic CRS must be projected and in metres");
        }
    }

    // -- G6: the OCN engine honours use_geographic_anchor --------------------
    {
        lsdtt::lsdtt_ocn_options opts;
        opts.n_rows      = 16;
        opts.n_cols      = 20;
        opts.cell_size   = 40.0;
        opts.n_iter      = 100;
        opts.seed        = 7;
        opts.gamma       = 0.5;
        opts.temperature_scale = 0.0;
        opts.outlet_mode = lsdtt::ocn_outlet_mode::single_south;

        // Default: Watopia anchor in EPSG:32758.
        {
            lsdtt::lsdtt_synthetic_network_cpu engine(opts);
            engine.run();
            const auto meta = engine.get_elevation_raster().meta;

            double e = 0.0, n = 0.0;
            const bool projected = project_point(-11.633, 166.967, 32758, e, n);
            const double expect_x = e - 0.5 * opts.n_cols * opts.cell_size;
            const double expect_y = n - 0.5 * opts.n_rows * opts.cell_size;

            std::ostringstream d;
            d << "epsg=" << meta.epsg_code << " x_min=" << std::setprecision(12)
              << meta.x_min << " expected " << expect_x;
            check_case(case_name,
                       "G6 OCN without geographic anchor uses the default anchor",
                       meta.epsg_code == 32758 && projected &&
                           std::abs(meta.x_min - expect_x) < 1e-3 &&
                           std::abs(meta.y_min - expect_y) < 1e-3,
                       d.str());
        }

        // Explicit anchor with an explicit EPSG.
        {
            lsdtt::lsdtt_ocn_options anchored = opts;
            anchored.use_geographic_anchor = true;
            anchored.centre_latitude_deg   = 55.95;
            anchored.centre_longitude_deg  = -3.19;
            anchored.projection_epsg       = 32630;

            lsdtt::lsdtt_synthetic_network_cpu engine(anchored);
            engine.run();
            const auto meta = engine.get_elevation_raster().meta;

            double e = 0.0, n = 0.0;
            const bool projected = project_point(55.95, -3.19, 32630, e, n);
            const double centre_x = meta.x_min + 0.5 * anchored.n_cols * anchored.cell_size;
            const double centre_y = meta.y_min + 0.5 * anchored.n_rows * anchored.cell_size;

            std::ostringstream d;
            d << "epsg=" << meta.epsg_code << " centre=(" << std::setprecision(12)
              << centre_x << "," << centre_y << ") anchor=(" << e << "," << n << ")";
            check_case(case_name,
                       "G6 OCN with geographic anchor carries the requested CRS",
                       meta.epsg_code == 32630 &&
                           epsg_of_wkt(meta.crs_wkt) == 32630 && projected &&
                           std::abs(centre_x - e) < 1e-6 &&
                           std::abs(centre_y - n) < 1e-6,
                       d.str());

            // The drainage-area raster must carry the same georeference.
            const auto area_meta = engine.get_drainage_area_cells_raster().meta;
            check_case(case_name,
                       "G6 OCN area raster shares the elevation georeference",
                       area_meta.epsg_code == meta.epsg_code &&
                           area_meta.x_min == meta.x_min &&
                           area_meta.y_min == meta.y_min &&
                           area_meta.crs_wkt == meta.crs_wkt,
                       "all OCN outputs must share one georeference");
        }

        // Explicit anchor with epsg = 0 -> derived UTM zone.
        {
            lsdtt::lsdtt_ocn_options anchored = opts;
            anchored.use_geographic_anchor = true;
            anchored.centre_latitude_deg   = -33.87;
            anchored.centre_longitude_deg  = 151.21;
            anchored.projection_epsg       = 0;

            lsdtt::lsdtt_synthetic_network_cpu engine(anchored);
            engine.run();
            const auto meta = engine.get_elevation_raster().meta;

            std::ostringstream d;
            d << "epsg=" << meta.epsg_code << " expected 32756";
            check_case(case_name,
                       "G6 OCN with epsg=0 derives the UTM zone from the anchor",
                       meta.epsg_code == 32756 &&
                           epsg_of_wkt(meta.crs_wkt) == 32756,
                       d.str());
        }
    }
}

} // namespace

int main(int argc, char** argv)
{
    const std::string out_dir = (argc > 1) ? argv[1] : ".";
    std::cout << "synthetic-networks scaffold output dir: " << out_dir << "\n";

    lsdtt::lsdtt_ocn_options base_opts;
    base_opts.n_rows = 8;
    base_opts.n_cols = 10;
    base_opts.cell_size = 50.0;
    base_opts.n_iter = 200;
    base_opts.seed = 42;
    base_opts.outlet_mode = lsdtt::ocn_outlet_mode::single_south;

    std::vector<synthetic_case> cases;

    synthetic_case case_a;
    case_a.name = "A rectangle 8x10 gamma=0.5 Tscale=0";
    case_a.opts = base_opts;
    case_a.opts.gamma = 0.5;
    case_a.opts.temperature_scale = 0.0;
    cases.push_back(case_a);

    synthetic_case case_b;
    case_b.name = "B rectangle 8x10 gamma=0.75 Tscale=1";
    case_b.opts = base_opts;
    case_b.opts.gamma = 0.75;
    case_b.opts.temperature_scale = 1.0;
    cases.push_back(case_b);

    synthetic_case case_c;
    case_c.name = "C masked island-like gamma=0.5 Tscale=0";
    case_c.opts = base_opts;
    case_c.opts.gamma = 0.5;
    case_c.opts.temperature_scale = 0.0;
    case_c.opts.domain_mask = build_island_mask(base_opts.n_rows, base_opts.n_cols);
    cases.push_back(case_c);

    for (const auto& sc : cases)
    {
        run_case(sc);
    }

    lsdtt::lsdtt_ocn_options deterministic_opts = base_opts;
    deterministic_opts.gamma = 0.5;
    deterministic_opts.temperature_scale = 0.0;
    run_determinism_checks("D determinism rectangle 8x10 gamma=0.5 Tscale=0",
                           deterministic_opts);

    run_routing_fidelity_checks("E D8 round-trip fidelity 48x48 gamma=0.5 seed=42");

    run_georeference_checks("G synthetic georeferencing");

    std::cout << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL")
              << "  (" << (g_n_checks - g_n_failed) << "/" << g_n_checks
              << " checks passed)\n";

    return (g_n_failed == 0) ? 0 : 1;
}
