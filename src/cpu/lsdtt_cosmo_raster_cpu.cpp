// =============================================================================
//  DO NOT CHANGE THE NUMBERS THIS CODE PRODUCES — ask Simon Mudd first.
//
//  This applies to the concentration and erosion-rate (inversion) calculations
//  in this file. The raster reading/writing and flow-routing helpers here are
//  ordinary code and can be changed normally.
//
//  Those cosmogenic calculations were checked by hand over several years. To
//  stop anyone changing their results by accident, the test suite keeps a stored
//  list of known-correct output values in
//  tests/expected/cosmo/numerics_sentinel.csv. The tests in
//  tests/cpp/test_cosmo_physics.cpp re-run these calculations and compare the
//  results against that stored list to 12 decimal places.
//
//  So if you edit those calculations and any output value changes — even from a
//  refactor or an algebra rewrite you believe is harmless — those tests will
//  fail. When that happens, do NOT edit the stored list to make the tests pass
//  again. A changed value has to be reviewed and approved by Simon Mudd before
//  it is committed.
// =============================================================================
#include "lsdtt_cosmo_raster_cpu.hpp"
#include "lsdtt_crn_particle.hpp"
#include "../utils/lsdtt_array2d.hpp"
#include "../core/lsdtt_raster_meta.hpp"
#include <ogr_spatialref.h>
#include <ogr_api.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt
{

// 
// Internal helpers
// 

static OGRCoordinateTransformation* make_to_wgs84(const lsdtt_raster_meta& meta)
{
    OGRSpatialReference src, dst;
    src.importFromWkt(meta.crs_wkt.c_str());
    dst.importFromEPSG(4326);
    dst.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    return OGRCreateCoordinateTransformation(&src, &dst);
}

// 
// Precomputed basin summary for O(4) per-evaluation inversion.
//
// With a spatially uniform erosion rate epsilon, the erosion-rate-weighted average
// concentration reduces to a plain mean (epsilon cancels). The per-pixel steady-state
// formula then factors as:
//
//   C_acc(epsilon) = P0 / n_pixels x Sigma_j Bsum[j] / (epsilon + Gamma[j]xlambda)
//
// where Bsum[j] = Sigma_pixels [ F_scaled[j] x W[j] ] is computed once in O(N).
//
// For a point sample (d_self == 0):  W[j] = Gamma[j] x exp(-d_snow/Gamma[j])
// For depth-integrated (d_self > 0): W[j] = Gamma[j]^2 x (exp(-d_top/G) - exp(-d_bot/G))
//                                             / d_self
//
// The analytic derivative is:
//   dC_acc/depsilon = -P0 / n_pixels x Sigma_j Bsum[j] / (epsilon + Gamma[j]xlambda)^2
// 
struct CrnBasinPrecomp
{
    double P0;
    double lambda;
    double Gamma[4];
    double Bsum[4];
    int    n_pixels;
};

static CrnBasinPrecomp precompute_basin_crn(
    const std::string&      nuclide,
    const std::string&      muon_scheme,
    const lsdtt_raster_cpu& prod_raster,
    const lsdtt_raster_cpu& topo_shield,
    const lsdtt_raster_cpu& snow_shield,
    const lsdtt_raster_cpu& self_shield,
    const lsdtt_flowinfo&   flow,
    std::int64_t            outlet_node)
{
    // Reference params for Gamma, P0, lambda (not modified by scale_F_values)
    lsdtt_crn_parameters ref;
    ref.reset_scaling(muon_scheme);

    CrnBasinPrecomp pre{};
    pre.P0 = nuclide == "Al26" ? ref.get_P0_26Al() :
             nuclide == "C14"  ? ref.get_P0_14C()  : ref.get_P0_10Be();
    pre.lambda = nuclide == "Al26" ? ref.get_lambda_26Al() :
                 nuclide == "C14"  ? ref.get_lambda_14C()  : ref.get_lambda_10Be();
    for (int j = 0; j < 4; ++j) { pre.Gamma[j] = ref.get_Gamma(j); pre.Bsum[j] = 0.0; }

    const auto& meta   = prod_raster.get_meta();
    int nc             = meta.n_cols;
    lsdtt_flat_index_slice upslope = flow.upslope_nodes(outlet_node);
    int n = static_cast<int>(upslope.size());

    const auto& p_data  = prod_raster.get_data_array();
    const auto& p_mask  = prod_raster.get_mask_array();
    const auto& ts_data = topo_shield.get_data_array();
    const auto& ts_mask = topo_shield.get_mask_array();
    const auto& snow_d  = snow_shield.get_data_array();
    const auto& snow_m  = snow_shield.get_mask_array();
    const auto& self_d  = self_shield.get_data_array();
    const auto& self_m  = self_shield.get_mask_array();

    double bsum0 = 0.0, bsum1 = 0.0, bsum2 = 0.0, bsum3 = 0.0;
    int n_pixels = 0;

#pragma omp parallel default(none) \
    shared(n, upslope, nc, p_data, p_mask, ts_data, ts_mask, \
           snow_d, snow_m, self_d, self_m, nuclide, muon_scheme, pre) \
    reduction(+:bsum0, bsum1, bsum2, bsum3, n_pixels)
    {
        lsdtt_crn_parameters params;

#pragma omp for schedule(static)
        for (int k = 0; k < n; ++k)
        {
            std::int64_t flat = upslope[k];
            int row = static_cast<int>(flat / nc);
            int col = static_cast<int>(flat % nc);

            if (!p_mask(row, col) || !ts_mask(row, col) ||
                !snow_m(row, col) || !self_m(row, col))
                continue;

            double shielding = static_cast<double>(p_data(row, col)) *
                               static_cast<double>(ts_data(row, col));
            double d_snow = static_cast<double>(snow_d(row, col));
            double d_self = static_cast<double>(self_d(row, col));
            double d_top  = d_snow;
            double d_bot  = d_snow + d_self;

            params.reset_scaling(muon_scheme);
            params.scale_F_values(shielding, nuclide);

            double F[4];
            for (int j = 0; j < 4; ++j)
                F[j] = nuclide == "Al26" ? params.get_F_26Al(j) :
                       nuclide == "C14"  ? params.get_F_14C(j)  : params.get_F_10Be(j);

            // Per-pathway depth weights
            double W[4];
            if (d_top == d_bot)
            {
                for (int j = 0; j < 4; ++j)
                    W[j] = pre.Gamma[j] * std::exp(-d_top / pre.Gamma[j]);
            }
            else
            {
                double inv_dself = 1.0 / (d_bot - d_top);
                for (int j = 0; j < 4; ++j)
                {
                    double G = pre.Gamma[j];
                    W[j] = G * G * (std::exp(-d_top / G) - std::exp(-d_bot / G)) * inv_dself;
                }
            }

            bsum0 += F[0] * W[0];
            bsum1 += F[1] * W[1];
            bsum2 += F[2] * W[2];
            bsum3 += F[3] * W[3];
            ++n_pixels;
        }
    } // end omp parallel

    pre.Bsum[0] = bsum0; pre.Bsum[1] = bsum1;
    pre.Bsum[2] = bsum2; pre.Bsum[3] = bsum3;
    pre.n_pixels = n_pixels;
    return pre;
}

// C_acc(epsilon) = P0/N x Sigma_j Bsum[j] / (epsilon + Gamma[j]xlambda)
static double eval_conc_analytic(const CrnBasinPrecomp& pre, double erate)
{
    if (pre.n_pixels == 0) return 0.0;
    double sum = 0.0;
    for (int j = 0; j < 4; ++j)
        sum += pre.Bsum[j] / (erate + pre.Gamma[j] * pre.lambda);
    return pre.P0 / pre.n_pixels * sum;
}

// dC_acc/depsilon = -P0/N x Sigma_j Bsum[j] / (epsilon + Gamma[j]xlambda)^2
static double eval_conc_deriv_analytic(const CrnBasinPrecomp& pre, double erate)
{
    if (pre.n_pixels == 0) return -1.0;  // non-zero to avoid divide-by-zero
    double sum = 0.0;
    for (int j = 0; j < 4; ++j)
    {
        double d = erate + pre.Gamma[j] * pre.lambda;
        sum += pre.Bsum[j] / (d * d);
    }
    return -pre.P0 / pre.n_pixels * sum;
}

// Initial erosion rate guess from average shielding (spallation-only approximation).
// epsilon ~= avg_shield x P0 x Gamma_0 / C_meas - Gamma_0 x lambda
static double erate_initial_guess(double measured_conc,
                                  const CrnBasinPrecomp& pre)
{
    if (pre.n_pixels == 0 || measured_conc <= 0.0) return 0.01;

    // Average Bsum[0]/Gamma[0] approximates avg_shield * P0 * Gamma[0] / Gamma[0] = avg_shield * P0
    // Simple spallation-only inversion of C ~= avg_scaled_prod x Gamma[0] / (epsilon + Gamma[0]xlambda)
    double avg_prod_term = pre.P0 * pre.Bsum[0] / pre.n_pixels;
    double guess = avg_prod_term / measured_conc - pre.Gamma[0] * pre.lambda;
    if (guess <= 1e-8) guess = 0.001;
    return guess;
}

// 
lsdtt_raster_cpu calculate_production_raster(const lsdtt_raster_cpu& dem,
                                              const CrnNcData&         nc,
                                              double                   age_yr)
{
    const auto& meta = dem.get_meta();
    const auto& data = dem.get_data_array();
    const auto& mask = dem.get_mask_array();
    int nr = meta.n_rows, nc_cols = meta.n_cols;
    float ndv = static_cast<float>(meta.no_data_value);

    lsdtt_array2d<float>   prod_data(nr, nc_cols, ndv);
    lsdtt_array2d<uint8_t> prod_mask(nr, nc_cols, 0);

    OGRCoordinateTransformation* ct = make_to_wgs84(meta);
    if (!ct)
        throw std::runtime_error(
            "calculate_production_raster: could not create coordinate transform. "
            "DEM must be in a projected CRS.");

    // Collect valid pixel indices and compute their projected coordinates.
    // Then batch-transform all at once (one GDAL call instead of one per pixel).
    std::vector<int>    valid_r, valid_c;
    std::vector<double> xs, ys;
    valid_r.reserve(static_cast<std::size_t>(nr * nc_cols / 2));
    valid_c.reserve(valid_r.capacity());
    xs.reserve(valid_r.capacity());
    ys.reserve(valid_r.capacity());

    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc_cols; ++c)
            if (mask(r, c))
            {
                valid_r.push_back(r);
                valid_c.push_back(c);
                xs.push_back(meta.x_min + (c + 0.5) * meta.cell_size);
                ys.push_back(meta.y_min + (meta.n_rows - r - 0.5) * meta.cell_size);
            }

    int n_valid = static_cast<int>(xs.size());
    if (n_valid == 0)
    {
        OCTDestroyCoordinateTransformation(ct);
        return lsdtt_raster_cpu(std::move(prod_data), std::move(prod_mask), meta);
    }

    // Single batch transform: projected -> WGS84. Not parallelised - GDAL transforms
    // are not thread-safe for the same object, and one call is already fast.
    ct->Transform(n_valid, xs.data(), ys.data());
    OCTDestroyCoordinateTransformation(ct);

    // xs now holds longitudes, ys holds latitudes.
    const bool use_era40  = nc.loaded;
    const bool use_paleomag = nc.loaded && (age_yr > 0.0);

    lsdtt_crn_parameters ref_params;

#pragma omp parallel for schedule(static) default(none) \
    shared(n_valid, valid_r, valid_c, xs, ys, data, prod_data, prod_mask, \
           ref_params, meta, nc, age_yr, use_era40, use_paleomag)
    for (int k = 0; k < n_valid; ++k)
    {
        int r = valid_r[k], c = valid_c[k];
        double lat  = ys[k];
        double lon  = xs[k];
        double elev = static_cast<double>(data(r, c));

        double pressure = use_era40
            ? lsdtt_crn_parameters::atmospheric_pressure_era40(lat, lon, elev, nc)
            : lsdtt_crn_parameters::atmospheric_pressure_standard(elev);

        double scaling = use_paleomag
            ? ref_params.stone2000sp_paleomag(lat, pressure, 1.0, age_yr, nc)
            : ref_params.stone2000sp(lat, pressure, 1.0);

        prod_data(r, c) = static_cast<float>(scaling);
        prod_mask(r, c) = 1;
    }

    return lsdtt_raster_cpu(std::move(prod_data), std::move(prod_mask), meta);
}

// 
lsdtt_raster_cpu calculate_crn_concentration_raster(
    const std::string&      nuclide,
    const std::string&      muon_scheme,
    const lsdtt_raster_cpu& erosion_raster,
    const lsdtt_raster_cpu& prod_raster,
    const lsdtt_raster_cpu& topo_shield,
    const lsdtt_raster_cpu& snow_shield,
    const lsdtt_raster_cpu& self_shield,
    const lsdtt_flowinfo&   flow,
    std::int64_t            outlet_node)
{
    const auto& meta = prod_raster.get_meta();
    int nr = meta.n_rows, nc = meta.n_cols;
    float ndv = static_cast<float>(meta.no_data_value);

    lsdtt_array2d<float>   out_data(nr, nc, ndv);
    lsdtt_array2d<uint8_t> out_mask(nr, nc, 0);

    lsdtt_flat_index_slice upslope = flow.upslope_nodes(outlet_node);
    int n = static_cast<int>(upslope.size());

    const auto& e_data  = erosion_raster.get_data_array();
    const auto& e_mask  = erosion_raster.get_mask_array();
    const auto& p_data  = prod_raster.get_data_array();
    const auto& p_mask  = prod_raster.get_mask_array();
    const auto& ts_data = topo_shield.get_data_array();
    const auto& ts_mask = topo_shield.get_mask_array();
    const auto& snow_d  = snow_shield.get_data_array();
    const auto& snow_m  = snow_shield.get_mask_array();
    const auto& self_d  = self_shield.get_data_array();
    const auto& self_m  = self_shield.get_mask_array();

    // Each pixel is independent. Each thread gets private params/particle.
#pragma omp parallel default(none) \
    shared(n, upslope, nc, e_data, e_mask, p_data, p_mask, ts_data, ts_mask, \
           snow_d, snow_m, self_d, self_m, out_data, out_mask, nuclide, muon_scheme)
    {
        lsdtt_crn_parameters params;
        lsdtt_crn_particle   particle;

#pragma omp for schedule(static)
        for (int k = 0; k < n; ++k)
        {
            std::int64_t flat = upslope[k];
            int row = static_cast<int>(flat / nc);
            int col = static_cast<int>(flat % nc);

            if (!e_mask(row, col) || !p_mask(row, col) ||
                !ts_mask(row, col) || !snow_m(row, col) || !self_m(row, col))
                continue;

            double erate         = static_cast<double>(e_data(row, col));
            double total_shielding = static_cast<double>(p_data(row, col)) *
                                     static_cast<double>(ts_data(row, col));
            double d_snow        = static_cast<double>(snow_d(row, col));
            double d_self        = static_cast<double>(self_d(row, col));
            double d_top         = d_snow;
            double d_bot         = d_snow + d_self;

            params.reset_scaling(muon_scheme);
            params.scale_F_values(total_shielding, nuclide);

            double conc = 0.0;
            if (nuclide == "Al26")
            {
                particle.update_26Al_ss(erate, params, d_top, d_bot);
                conc = particle.get_conc_26Al();
            }
            else if (nuclide == "C14")
            {
                particle.update_14C_ss(erate, params, d_top, d_bot);
                conc = particle.get_conc_14C();
            }
            else
            {
                particle.update_10Be_ss(erate, params, d_top, d_bot);
                conc = particle.get_conc_10Be();
            }

            out_data(row, col) = static_cast<float>(conc);
            out_mask(row, col) = 1;
        }
    } // end omp parallel

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}

// 
lsdtt_raster_cpu calculate_crn_concentration_raster_whole(
    const std::string&      nuclide,
    const std::string&      muon_scheme,
    const lsdtt_raster_cpu& erosion_raster,
    const lsdtt_raster_cpu& prod_raster,
    const lsdtt_raster_cpu& topo_shield,
    const lsdtt_raster_cpu& snow_shield,
    const lsdtt_raster_cpu& self_shield)
{
    const auto& meta = prod_raster.get_meta();
    int nr = meta.n_rows, nc = meta.n_cols;
    int n_total = nr * nc;
    float ndv = static_cast<float>(meta.no_data_value);

    lsdtt_array2d<float>   out_data(nr, nc, ndv);
    lsdtt_array2d<uint8_t> out_mask(nr, nc, 0);

    const auto& e_data  = erosion_raster.get_data_array();
    const auto& e_mask  = erosion_raster.get_mask_array();
    const auto& p_data  = prod_raster.get_data_array();
    const auto& p_mask  = prod_raster.get_mask_array();
    const auto& ts_data = topo_shield.get_data_array();
    const auto& ts_mask = topo_shield.get_mask_array();
    const auto& snow_d  = snow_shield.get_data_array();
    const auto& snow_m  = snow_shield.get_mask_array();
    const auto& self_d  = self_shield.get_data_array();
    const auto& self_m  = self_shield.get_mask_array();

#pragma omp parallel default(none) \
    shared(nr, nc, n_total, e_data, e_mask, p_data, p_mask, ts_data, ts_mask, \
           snow_d, snow_m, self_d, self_m, out_data, out_mask, nuclide, muon_scheme)
    {
        lsdtt_crn_parameters params;
        lsdtt_crn_particle   particle;

#pragma omp for schedule(static)
        for (int k = 0; k < n_total; ++k)
        {
            int row = k / nc;
            int col = k % nc;

            if (!e_mask(row, col) || !p_mask(row, col) ||
                !ts_mask(row, col) || !snow_m(row, col) || !self_m(row, col))
                continue;

            double erate           = static_cast<double>(e_data(row, col));
            double total_shielding = static_cast<double>(p_data(row, col)) *
                                     static_cast<double>(ts_data(row, col));
            double d_top = static_cast<double>(snow_d(row, col));
            double d_bot = d_top + static_cast<double>(self_d(row, col));

            params.reset_scaling(muon_scheme);
            params.scale_F_values(total_shielding, nuclide);

            double conc = 0.0;
            if (nuclide == "Al26")
            {
                particle.update_26Al_ss(erate, params, d_top, d_bot);
                conc = particle.get_conc_26Al();
            }
            else if (nuclide == "C14")
            {
                particle.update_14C_ss(erate, params, d_top, d_bot);
                conc = particle.get_conc_14C();
            }
            else
            {
                particle.update_10Be_ss(erate, params, d_top, d_bot);
                conc = particle.get_conc_10Be();
            }

            out_data(row, col) = static_cast<float>(conc);
            out_mask(row, col) = 1;
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}

// 
lsdtt_raster_cpu calculate_crn_concentration_raster_step_change(
    const std::string&      nuclide,
    const std::string&      muon_scheme,
    const lsdtt_raster_cpu& time_since_change,
    const lsdtt_raster_cpu& prod_raster,
    const lsdtt_raster_cpu& topo_shield,
    double                  eff_erate_old_g_cm2_yr,
    double                  eff_erate_new_g_cm2_yr)
{
    const auto& meta = prod_raster.get_meta();
    int nr = meta.n_rows, nc = meta.n_cols;
    int n_total = nr * nc;
    float ndv = static_cast<float>(meta.no_data_value);

    lsdtt_array2d<float>   out_data(nr, nc, ndv);
    lsdtt_array2d<uint8_t> out_mask(nr, nc, 0);

    const auto& p_data  = prod_raster.get_data_array();
    const auto& p_mask  = prod_raster.get_mask_array();
    const auto& ts_data = topo_shield.get_data_array();
    const auto& ts_mask = topo_shield.get_mask_array();
    const auto& t_data  = time_since_change.get_data_array();
    const auto& t_mask  = time_since_change.get_mask_array();

#pragma omp parallel default(none) \
    shared(nr, nc, n_total, p_data, p_mask, ts_data, ts_mask, t_data, t_mask, \
           out_data, out_mask, nuclide, muon_scheme, \
           eff_erate_old_g_cm2_yr, eff_erate_new_g_cm2_yr)
    {
        lsdtt_crn_parameters params;
        lsdtt_crn_particle   particle;

#pragma omp for schedule(static)
        for (int k = 0; k < n_total; ++k)
        {
            int row = k / nc;
            int col = k % nc;

            if (!p_mask(row, col) || !ts_mask(row, col) || !t_mask(row, col))
                continue;

            double total_shielding = static_cast<double>(p_data(row, col)) *
                                     static_cast<double>(ts_data(row, col));
            double t_step = static_cast<double>(t_data(row, col));

            params.reset_scaling(muon_scheme);
            params.scale_F_values(total_shielding, nuclide);

            double conc = 0.0;

            if (t_step > 0.0)
            {
                // Particle was at starting depth = eff_erate_new * t_step at t_step yr ago.
                // SS at eff_erate_old, then evolved for t_step yr at eff_erate_new.
                double d_start = eff_erate_new_g_cm2_yr * t_step;
                particle.set_effective_depth(d_start);
                particle.set_conc_10Be(0.0);
                particle.set_conc_26Al(0.0);
                particle.set_conc_14C(0.0);

                if (nuclide == "Al26")
                {
                    particle.update_26Al_ss(eff_erate_old_g_cm2_yr, params, d_start, d_start);
                    particle.update_conc(t_step, eff_erate_new_g_cm2_yr, params);
                    conc = particle.get_conc_26Al();
                }
                else if (nuclide == "C14")
                {
                    particle.update_14C_ss(eff_erate_old_g_cm2_yr, params, d_start, d_start);
                    particle.update_conc(t_step, eff_erate_new_g_cm2_yr, params);
                    conc = particle.get_conc_14C();
                }
                else
                {
                    particle.update_10Be_ss(eff_erate_old_g_cm2_yr, params, d_start, d_start);
                    particle.update_conc(t_step, eff_erate_new_g_cm2_yr, params);
                    conc = particle.get_conc_10Be();
                }
            }
            else
            {
                // No step change at this pixel: SS at new erosion rate (surface sample)
                if (nuclide == "Al26")
                {
                    particle.update_26Al_ss(eff_erate_new_g_cm2_yr, params, 0.0, 0.0);
                    conc = particle.get_conc_26Al();
                }
                else if (nuclide == "C14")
                {
                    particle.update_14C_ss(eff_erate_new_g_cm2_yr, params, 0.0, 0.0);
                    conc = particle.get_conc_14C();
                }
                else
                {
                    particle.update_10Be_ss(eff_erate_new_g_cm2_yr, params, 0.0, 0.0);
                    conc = particle.get_conc_10Be();
                }
            }

            out_data(row, col) = static_cast<float>(conc);
            out_mask(row, col) = 1;
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}

// 
lsdtt_raster_cpu calculate_accumulated_crn_raster(
    const lsdtt_raster_cpu& conc_raster,
    const lsdtt_raster_cpu& erosion_raster,
    const lsdtt_flowinfo&   flow)
{
    const auto& meta   = conc_raster.get_meta();
    int nr = meta.n_rows, nc_cols = meta.n_cols;
    int n_total = nr * nc_cols;
    float ndv = static_cast<float>(meta.no_data_value);

    const auto& c_data = conc_raster.get_data_array();
    const auto& c_mask = conc_raster.get_mask_array();
    const auto& e_data = erosion_raster.get_data_array();
    const auto& e_mask = erosion_raster.get_mask_array();

    // Build per-cell weight vectors for accumulation.
    // nodata cells contribute 0 to both sums.
    std::vector<double> atoms_flux(static_cast<std::size_t>(n_total), 0.0);
    std::vector<double> erosion_wt(static_cast<std::size_t>(n_total), 0.0);

    for (int k = 0; k < n_total; ++k)
    {
        int row = k / nc_cols;
        int col = k % nc_cols;
        if (!c_mask(row, col) || !e_mask(row, col)) continue;
        double erate = static_cast<double>(e_data(row, col));
        double conc  = static_cast<double>(c_data(row, col));
        atoms_flux[static_cast<std::size_t>(k)] = conc * erate;
        erosion_wt[static_cast<std::size_t>(k)] = erate;
    }

    std::vector<double> acc_atoms = flow.accumulate(atoms_flux);
    std::vector<double> acc_eros  = flow.accumulate(erosion_wt);

    lsdtt_array2d<float>   out_data(nr, nc_cols, ndv);
    lsdtt_array2d<uint8_t> out_mask(nr, nc_cols, 0);

    for (int k = 0; k < n_total; ++k)
    {
        int row = k / nc_cols;
        int col = k % nc_cols;
        double acc_e = acc_eros[static_cast<std::size_t>(k)];
        if (acc_e <= 0.0) continue;
        out_data(row, col) = static_cast<float>(acc_atoms[static_cast<std::size_t>(k)] / acc_e);
        out_mask(row, col) = 1;
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}

// 
double calculate_accumulated_crn_concentration(
    const lsdtt_raster_cpu& conc_raster,
    const lsdtt_raster_cpu& erosion_raster,
    const lsdtt_flowinfo&   flow,
    std::int64_t            outlet_node)
{
    const auto& meta   = conc_raster.get_meta();
    int nc             = meta.n_cols;
    const auto& c_data = conc_raster.get_data_array();
    const auto& c_mask = conc_raster.get_mask_array();
    const auto& e_data = erosion_raster.get_data_array();
    const auto& e_mask = erosion_raster.get_mask_array();

    lsdtt_flat_index_slice upslope = flow.upslope_nodes(outlet_node);

    double acc_mass = 0.0, acc_atoms = 0.0;
    for (std::int64_t flat : upslope)
    {
        int row = static_cast<int>(flat / nc);
        int col = static_cast<int>(flat % nc);
        if (!c_mask(row, col) || !e_mask(row, col)) continue;
        double erate = static_cast<double>(e_data(row, col));
        double conc  = static_cast<double>(c_data(row, col));
        acc_mass  += erate;
        acc_atoms += conc * erate;
    }

    if (acc_mass == 0.0) return 0.0;
    return acc_atoms / acc_mass;
}

// 
double calculate_eff_erate_from_conc(
    double                  measured_conc,
    const std::string&      nuclide,
    const std::string&      muon_scheme,
    const lsdtt_raster_cpu& prod_raster,
    const lsdtt_raster_cpu& topo_shield,
    const lsdtt_raster_cpu& snow_shield,
    const lsdtt_raster_cpu& self_shield,
    const lsdtt_flowinfo&   flow,
    std::int64_t            outlet_node,
    double&                 predicted_conc)
{
    // Guard: non-positive concentration has no physical meaning
    if (measured_conc <= 0.0)
    {
        std::cout << "calculate_eff_erate_from_conc: measured_conc <= 0 -- returning 0.\n";
        predicted_conc = 0.0;
        return 0.0;
    }

    //  Precompute basin summary - O(N_pixels) once 
    // All subsequent C(epsilon) and dC/depsilon evaluations cost O(4).
    CrnBasinPrecomp pre = precompute_basin_crn(
        nuclide, muon_scheme, prod_raster, topo_shield,
        snow_shield, self_shield, flow, outlet_node);

    if (pre.n_pixels == 0)
    {
        std::cout << "calculate_eff_erate_from_conc: no valid upslope pixels -- returning 0.\n";
        predicted_conc = 0.0;
        return 0.0;
    }

    //  Tolerances 
    const int    nr_max_steps  = 50;          // N-R iteration cap
    const double conc_tol_frac = 0.002;       // accept when |DeltaC|/C_meas < 0.2%
    const double erate_min     = 1e-8;        // floor to keep epsilon > 0

    const double conc_tol = conc_tol_frac * measured_conc;

    //  Initial guess 
    double eff_e = erate_initial_guess(measured_conc, pre);

    //  Newton-Raphson with analytic derivative 
    // C(epsilon) is monotonically decreasing, so the derivative is always negative
    // and the N-R step always moves in the correct direction.  The only failure
    // mode is overshooting to epsilon <= 0, which the step limiter handles.
    //
    // We also track the tightest bracket [lo_e, hi_e] where
    //   C(lo_e) > C_meas  and  C(hi_e) < C_meas.
    // If N-R converges cleanly we never need bisection.

    bool have_lo = false, have_hi = false;
    double lo_e = erate_min, hi_e = 1000.0;   // physical limits: ~380,000 mm/kyr
    double C_cur = 0.0;
    int n_steps = 0;

    for (; n_steps < nr_max_steps; ++n_steps)
    {
        C_cur = eval_conc_analytic(pre, eff_e);
        double fx    = C_cur - measured_conc;
        double deriv = eval_conc_deriv_analytic(pre, eff_e);

        // Update brackets
        if (fx > 0.0) { have_lo = true; if (eff_e < hi_e) lo_e = std::max(lo_e, eff_e); }
        if (fx < 0.0) { have_hi = true; if (eff_e > lo_e) hi_e = std::min(hi_e, eff_e); }

        if (std::fabs(fx) < conc_tol) break;    // converged

        if (deriv == 0.0) break;                 // flat - fall through to bisection

        double step = fx / deriv;               // deriv < 0, so step sign matches correction

        // Damp step: don't let epsilon drop below erate_min in one move
        if (eff_e - step < erate_min)
            step = eff_e - erate_min * 2.0;    // bring epsilon to 2xfloor

        eff_e -= step;
        if (eff_e < erate_min) eff_e = erate_min;
    }

    //  Check convergence; fall through to bisection if needed 
    C_cur = eval_conc_analytic(pre, eff_e);
    double fx_final = C_cur - measured_conc;
    if (std::fabs(fx_final) >= conc_tol)
    {
        // N-R did not converge to tolerance.  Ensure we have a valid bracket.
        if (!(have_lo && have_hi))
        {
            // Scan 40 log-spaced points across the full physical range to find both
            // sides of the root.  Covers epsilon in [1e-6, 1000] g/cm^2/yr.
            const int scan_n = 40;
            const double log_lo = std::log10(1e-6), log_hi = std::log10(1000.0);
            for (int s = 0; s <= scan_n && !(have_lo && have_hi); ++s)
            {
                double e_scan = std::pow(10.0, log_lo + s * (log_hi - log_lo) / scan_n);
                double C_scan = eval_conc_analytic(pre, e_scan);
                if (C_scan > measured_conc) { have_lo = true; lo_e = std::max(lo_e, e_scan); }
                if (C_scan < measured_conc) { have_hi = true; hi_e = std::min(hi_e, e_scan); }
            }
        }

        if (!have_lo || !have_hi)
        {
            // Solution not in [1e-6, 1000] - extreme concentration value.
            // Return the closest scan point and warn.
            std::cout << "calculate_eff_erate_from_conc: WARNING -- concentration "
                      << measured_conc << " is outside the bracketed range. "
                      << "Returning best estimate " << eff_e << " g/cm2/yr.\n";
            predicted_conc = C_cur;
            return eff_e;
        }

        // Ensure lo_e < hi_e (bracket is ordered)
        if (lo_e > hi_e) std::swap(lo_e, hi_e);

        //  Bisection 
        // lo_e: C(lo_e) > C_meas  (too low epsilon)
        // hi_e: C(hi_e) < C_meas  (too high epsilon)
        const int bisect_max = 200;
        for (int bs = 0; bs < bisect_max; ++bs)
        {
            double mid   = 0.5 * (lo_e + hi_e);
            double C_mid = eval_conc_analytic(pre, mid);
            double fx_mid = C_mid - measured_conc;

            if (C_mid > measured_conc) lo_e = mid;
            else                       hi_e = mid;

            // Converged in concentration or bracket collapsed
            if (std::fabs(fx_mid) < conc_tol || (hi_e - lo_e) < erate_min)
            {
                eff_e = mid;
                C_cur = C_mid;
                break;
            }
        }

        // If still not converged after bisection, use the midpoint of the bracket
        if (std::fabs(C_cur - measured_conc) >= conc_tol)
        {
            eff_e = 0.5 * (lo_e + hi_e);
            C_cur = eval_conc_analytic(pre, eff_e);
        }
    }

    //  Final convergence check 
    double rel_err = std::fabs((C_cur - measured_conc) / measured_conc);
    if (rel_err > 0.05)
    {
        std::cout << "calculate_eff_erate_from_conc: WARNING -- did not converge. "
                  << "Relative concentration error = " << rel_err * 100.0 << "%\n";
    }

    std::cout << "  n_steps=" << n_steps
              << "  erate=" << eff_e << " g/cm2/yr"
              << "  pred=" << C_cur << "  meas=" << measured_conc << "\n";

    predicted_conc = C_cur;
    return eff_e;
}

} // namespace lsdtt
