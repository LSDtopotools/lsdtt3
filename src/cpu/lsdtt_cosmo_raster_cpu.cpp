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

// ─────────────────────────────────────────────────────────────────────────────
// Internal helpers
// ─────────────────────────────────────────────────────────────────────────────

static OGRCoordinateTransformation* make_to_wgs84(const lsdtt_raster_meta& meta)
{
    OGRSpatialReference src, dst;
    src.importFromWkt(meta.crs_wkt.c_str());
    dst.importFromEPSG(4326);
    dst.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    return OGRCreateCoordinateTransformation(&src, &dst);
}

// ─────────────────────────────────────────────────────────────────────────────
// Precomputed basin summary for O(4) per-evaluation inversion.
//
// With a spatially uniform erosion rate ε, the erosion-rate-weighted average
// concentration reduces to a plain mean (ε cancels). The per-pixel steady-state
// formula then factors as:
//
//   C_acc(ε) = P0 / n_pixels × Σ_j Bsum[j] / (ε + Gamma[j]×λ)
//
// where Bsum[j] = Σ_pixels [ F_scaled[j] × W[j] ] is computed once in O(N).
//
// For a point sample (d_self == 0):  W[j] = Gamma[j] × exp(-d_snow/Gamma[j])
// For depth-integrated (d_self > 0): W[j] = Gamma[j]² × (exp(-d_top/G) - exp(-d_bot/G))
//                                             / d_self
//
// The analytic derivative is:
//   dC_acc/dε = -P0 / n_pixels × Σ_j Bsum[j] / (ε + Gamma[j]×λ)²
// ─────────────────────────────────────────────────────────────────────────────
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

// C_acc(ε) = P0/N × Σ_j Bsum[j] / (ε + Gamma[j]×λ)
static double eval_conc_analytic(const CrnBasinPrecomp& pre, double erate)
{
    if (pre.n_pixels == 0) return 0.0;
    double sum = 0.0;
    for (int j = 0; j < 4; ++j)
        sum += pre.Bsum[j] / (erate + pre.Gamma[j] * pre.lambda);
    return pre.P0 / pre.n_pixels * sum;
}

// dC_acc/dε = -P0/N × Σ_j Bsum[j] / (ε + Gamma[j]×λ)²
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
// ε ≈ avg_shield × P0 × Γ₀ / C_meas - Γ₀ × λ
static double erate_initial_guess(double measured_conc,
                                  const CrnBasinPrecomp& pre)
{
    if (pre.n_pixels == 0 || measured_conc <= 0.0) return 0.01;

    // Average Bsum[0]/Gamma[0] approximates avg_shield * P0 * Gamma[0] / Gamma[0] = avg_shield * P0
    // Simple spallation-only inversion of C ≈ avg_scaled_prod × Gamma[0] / (ε + Gamma[0]×λ)
    double avg_prod_term = pre.P0 * pre.Bsum[0] / pre.n_pixels;
    double guess = avg_prod_term / measured_conc - pre.Gamma[0] * pre.lambda;
    if (guess <= 1e-8) guess = 0.001;
    return guess;
}

// ─────────────────────────────────────────────────────────────────────────────
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

    // Single batch transform: projected → WGS84. Not parallelised — GDAL transforms
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

// ─────────────────────────────────────────────────────────────────────────────
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

// ─────────────────────────────────────────────────────────────────────────────
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

// ─────────────────────────────────────────────────────────────────────────────
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

    // ── Precompute basin summary — O(N_pixels) once ───────────────────────────
    // All subsequent C(ε) and dC/dε evaluations cost O(4).
    CrnBasinPrecomp pre = precompute_basin_crn(
        nuclide, muon_scheme, prod_raster, topo_shield,
        snow_shield, self_shield, flow, outlet_node);

    if (pre.n_pixels == 0)
    {
        std::cout << "calculate_eff_erate_from_conc: no valid upslope pixels -- returning 0.\n";
        predicted_conc = 0.0;
        return 0.0;
    }

    // ── Tolerances ────────────────────────────────────────────────────────────
    const int    nr_max_steps  = 50;          // N-R iteration cap
    const double conc_tol_frac = 0.002;       // accept when |ΔC|/C_meas < 0.2%
    const double erate_min     = 1e-8;        // floor to keep ε > 0

    const double conc_tol = conc_tol_frac * measured_conc;

    // ── Initial guess ─────────────────────────────────────────────────────────
    double eff_e = erate_initial_guess(measured_conc, pre);

    // ── Newton-Raphson with analytic derivative ───────────────────────────────
    // C(ε) is monotonically decreasing, so the derivative is always negative
    // and the N-R step always moves in the correct direction.  The only failure
    // mode is overshooting to ε ≤ 0, which the step limiter handles.
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

        if (deriv == 0.0) break;                 // flat — fall through to bisection

        double step = fx / deriv;               // deriv < 0, so step sign matches correction

        // Damp step: don't let ε drop below erate_min in one move
        if (eff_e - step < erate_min)
            step = eff_e - erate_min * 2.0;    // bring ε to 2×floor

        eff_e -= step;
        if (eff_e < erate_min) eff_e = erate_min;
    }

    // ── Check convergence; fall through to bisection if needed ────────────────
    C_cur = eval_conc_analytic(pre, eff_e);
    double fx_final = C_cur - measured_conc;
    if (std::fabs(fx_final) >= conc_tol)
    {
        // N-R did not converge to tolerance.  Ensure we have a valid bracket.
        if (!(have_lo && have_hi))
        {
            // Scan 40 log-spaced points across the full physical range to find both
            // sides of the root.  Covers ε in [1e-6, 1000] g/cm²/yr.
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
            // Solution not in [1e-6, 1000] — extreme concentration value.
            // Return the closest scan point and warn.
            std::cout << "calculate_eff_erate_from_conc: WARNING -- concentration "
                      << measured_conc << " is outside the bracketed range. "
                      << "Returning best estimate " << eff_e << " g/cm2/yr.\n";
            predicted_conc = C_cur;
            return eff_e;
        }

        // Ensure lo_e < hi_e (bracket is ordered)
        if (lo_e > hi_e) std::swap(lo_e, hi_e);

        // ── Bisection ─────────────────────────────────────────────────────────
        // lo_e: C(lo_e) > C_meas  (too low ε)
        // hi_e: C(hi_e) < C_meas  (too high ε)
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

    // ── Final convergence check ───────────────────────────────────────────────
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
