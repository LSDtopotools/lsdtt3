// test_cosmo_physics.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// Standalone validation driver for the cosmogenic-nuclide (CRN) physics core.
// Methodology, reference sources and findings: docs/testing_and_benchmarking.md.
//
// Covers: reference helpers, the single-basin synthetic fixture, CSV emitters,
// smoke tests, the parameter pins, the numerics sentinel, and the
// equation-fidelity battery E1-E9.
//
// *** COSMO NUMERICS WARNING ************************************************
// The analytical solutions exercised here (CAIRN steady-state / inversion;
// muon-based step change) encode years of manual verification by Simon Mudd.
// Any change to the numeric core (lsdtt_crn_parameters.*, lsdtt_crn_column.*,
// lsdtt_crn_particle.*, lsdtt_cosmo_raster_cpu.*) or to the reference values
// in tests/expected/cosmo/ requires his explicit review BEFORE merging.
// ***************************************************************************
//
// Reference equations (independent re-implementations, written from the papers
// and NOT calling the library, so library/reference agreement is meaningful):
//   - Mudd et al. (2016), "The CAIRN method", Earth Surf. Dynam. 4, 655-674:
//     Eq. (12) surface steady state; Eq. (14) depth-integrated steady state.
//   - Mudd & Roda-Boluda (draft), "Is it possible to detect a channel incision
//     wave using the 14C-10Be in-situ cosmogenic pair": Eq. (7) step-change
//     point solution; Eq. (8) its analytic depth integral.
// Further sources: Stone (2000); Granger & Smith (2000); Braucher et al.
// (2011); Borchers et al. (2016); Lupker et al. (2015); Vermeesch (2007);
// Balco et al. (2008).
//
// Usage: test-cosmo-physics [output_dir] [sentinel_csv_path]
//                           [--regenerate-sentinel-i-have-simons-approval]
//   output_dir (default ".") receives the plottable CSVs (plan Sec7).
//   sentinel_csv_path (default tests/expected/cosmo/numerics_sentinel.csv
//   relative to the CTest invocation) holds the Sec3b pinned values.
// Prints one PASS/FAIL line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_crn_parameters.hpp"
#include "lsdtt_cosmo_raster_cpu.hpp"
#include "lsdtt_exe_path.hpp"

namespace
{

// =============================================================================
// Tiny check harness
// =============================================================================
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

bool approx_rel(double a, double b, double rel_tol)
{
    const double scale = std::max(std::abs(a), std::abs(b));
    if (scale == 0.0) return true;
    return std::abs(a - b) <= rel_tol * scale;
}

std::string diff_str(double a, double b)
{
    std::ostringstream os;
    os << std::setprecision(17) << a << " vs " << b;
    return os.str();
}

// =============================================================================
// Consolidated frozen tolerances (plan Sec3b S6)
// =============================================================================
// Every frozen bound in this driver lives here, each carrying the measured
// value and the date it was frozen.  Loosening any of these to make a test
// pass is the SAME class of event as regenerating the numerics sentinel: it is
// gated and needs Simon Mudd's review (plan S6).  Do not relax a bound to
// silence a failure - a failure here is a result.
namespace tol
{
// Phases 0-3 (measured 2026-07-16):
constexpr double sentinel_rel = 1e-12;   // S1 numerics sentinel (pins the core)
constexpr double e7_analytic  = 5e-7;    // E7 analytic vs integrator; worst 9.25e-8
constexpr double e6_roundtrip = 3e-3;    // E6 forward/inverse round-trip; worst 9.93e-4
constexpr double tc2_cronus   = 10.0;    // T-C2 vs CRONUS v3 St, percent; worst 6.65
// Phase 5 geomagnetic (measured-then-frozen 2026-07-23; see plan Sec9b):
constexpr double g1_static    = 1e-12;   // G1 constant-VDM correction is exactly 1
constexpr double g4_golden    = 1e-9;    // G4 paleomag goldens; deterministic double core
constexpr double g6c_band     = 25.0;    // 6c(A) equatorial ON-OFF |delta|, percent
}

// =============================================================================
// Reference physics (plan Sec3) - independent of the library implementation
// =============================================================================

// One nuclide + scheme + per-pathway scaling, everything the reference
// equations need.  S[j] is the per-pathway scaling factor (the output of the
// Vermeesch Lambda_v redistribution); for an unscaled/unshielded surface case
// use S = {1,1,1,1}.
struct ref_params
{
    double P0     = 0.0;   // SLHL production, atoms g^-1 yr^-1
    double lambda = 0.0;   // decay constant, yr^-1
    double F[4]     = {0.0, 0.0, 0.0, 0.0};
    double Gamma[4] = {1.0, 1.0, 1.0, 1.0};   // g cm^-2
    double S[4]     = {1.0, 1.0, 1.0, 1.0};
};

// CAIRN Eq. (12): steady-state concentration at shielding depth d (g cm^-2)
// under steady effective erosion eps (g cm^-2 yr^-1).
double ref_conc_surface(const ref_params& p, double eps, double d)
{
    double c = 0.0;
    for (int j = 0; j < 4; ++j)
    {
        c += p.S[j] * p.F[j] * p.Gamma[j] * std::exp(-d / p.Gamma[j])
             / (eps + p.lambda * p.Gamma[j]);
    }
    return p.P0 * c;
}

// CAIRN Eq. (14): concentration depth-averaged between snow shielding depth
// d_s and d_s + d_t (eroded-layer thickness d_t), both g cm^-2.  This is the
// form the CAIRN inversion iterates on.
double ref_conc_depth_integrated(const ref_params& p, double eps,
                                 double d_s, double d_t)
{
    double c = 0.0;
    for (int j = 0; j < 4; ++j)
    {
        c += p.S[j] * p.F[j] * p.Gamma[j] * p.Gamma[j]
             * (std::exp(-d_s / p.Gamma[j])
                - std::exp(-(d_s + d_t) / p.Gamma[j]))
             / (eps + p.lambda * p.Gamma[j]);
    }
    return p.P0 * c / d_t;
}

// Knickpoint Eq. (7): concentration at shielding depth d, a time t after an
// instantaneous step change from eps_old to eps_new (both g cm^-2 yr^-1).
// Term 1: decayed inherited steady-state inventory (particle started at depth
// d + t*eps_new).  Term 2: accumulation at the new rate.
double ref_conc_step_change(const ref_params& p, double d, double t,
                            double eps_old, double eps_new)
{
    double c = 0.0;
    for (int j = 0; j < 4; ++j)
    {
        const double d0 = d + t * eps_new;   // initial depth
        const double inherited =
            std::exp(-p.lambda * t)
            * p.S[j] * p.F[j] * p.Gamma[j] * std::exp(-d0 / p.Gamma[j])
            / (p.lambda * p.Gamma[j] + eps_old);
        // Stable form of exp(-d0/G)*(exp(t*eps_n/G) - exp(-t*lambda)): the
        // printed Eq. (7) grouping underflows/overflows to 0*inf = NaN for
        // large t*eps_new/Gamma.  Algebraically identical (plan E3 identity).
        const double accumulated =
            p.S[j] * p.F[j] * p.Gamma[j]
            * (std::exp(-d / p.Gamma[j])
               - std::exp(-d0 / p.Gamma[j]) * std::exp(-t * p.lambda))
            / (p.lambda * p.Gamma[j] + eps_new);
        c += inherited + accumulated;
    }
    return p.P0 * c;
}

// Knickpoint Eq. (8): Eq. (7) integrated analytically between shielding
// depths d_t (top) and d_b (bottom), d_b > d_t, both g cm^-2.
double ref_conc_step_change_depth_integrated(const ref_params& p,
                                             double d_t, double d_b, double t,
                                             double eps_old, double eps_new)
{
    double c = 0.0;
    for (int j = 0; j < 4; ++j)
    {
        const double G  = p.Gamma[j];
        const double gl = G * p.lambda;
        const double term_ss =
            p.S[j] * G * G * p.F[j] / (gl + eps_new)
            * (std::exp(-d_t / G) - std::exp(-d_b / G));
        const double term_tr =
            p.S[j] * G * G * p.F[j] * (eps_new - eps_old)
            / ((gl + eps_old) * (gl + eps_new))
            * (std::exp(-(d_t + gl * t + t * eps_new) / G)
               - std::exp(-(d_b + gl * t + t * eps_new) / G));
        c += term_ss + term_tr;
    }
    return p.P0 * c / (d_b - d_t);
}

// Pulls a ref_params from the library's parameter object so the reference
// equations run with exactly the constants under test.
ref_params ref_params_from_scheme(const std::string& scheme,
                                  const std::string& nuclide)
{
    lsdtt::lsdtt_crn_parameters cp;
    cp.reset_scaling(scheme);
    ref_params p;
    for (int j = 0; j < 4; ++j)
    {
        p.Gamma[j] = cp.get_Gamma(j);
        if      (nuclide == "Be10") p.F[j] = cp.get_F_10Be(j);
        else if (nuclide == "Al26") p.F[j] = cp.get_F_26Al(j);
        else                        p.F[j] = cp.get_F_14C(j);
    }
    if      (nuclide == "Be10") { p.P0 = cp.get_P0_10Be(); p.lambda = cp.get_lambda_10Be(); }
    else if (nuclide == "Al26") { p.P0 = cp.get_P0_26Al(); p.lambda = cp.get_lambda_26Al(); }
    else                        { p.P0 = cp.get_P0_14C();  p.lambda = cp.get_lambda_14C();  }
    return p;
}

// =============================================================================
// Single-basin synthetic fixture (plan Sec4)
// =============================================================================
//
// A V-shaped valley draining to a single outlet at (row 0, centre column):
//   z(r,c) = down_slope * r * cell + cross_slope * |c - c_mid| * cell
// Every valid pixel has a monotone descending D8 path to the outlet, so the
// basin is single-outlet by construction (verified by a smoke test below).
// Projected CRS comes from set_default_synthetic_georeference() (Watopia /
// Teanu UTM 58S) so check_projected_metres passes.
struct single_basin_fixture
{
    lsdtt::lsdtt_raster_cpu dem;
    lsdtt::lsdtt_raster_meta meta;
    std::int64_t outlet_node = -1;   // flat index, row-major r * n_cols + c
    int outlet_row = 0;
    int outlet_col = 0;
};

single_basin_fixture make_single_basin(int n_rows, int n_cols,
                                       double cell_size = 30.0,
                                       double down_slope = 0.05,
                                       double cross_slope = 0.02)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows    = n_rows;
    meta.n_cols    = n_cols;
    meta.cell_size = cell_size;
    meta.set_default_synthetic_georeference();

    const int c_mid = n_cols / 2;
    lsdtt::lsdtt_array2d<float> z(n_rows, n_cols, 0.0f);
    lsdtt::lsdtt_array2d<std::uint8_t> mask(n_rows, n_cols, 1);
    for (int r = 0; r < n_rows; ++r)
    {
        for (int c = 0; c < n_cols; ++c)
        {
            z(r, c) = static_cast<float>(
                down_slope * r * cell_size
                + cross_slope * std::abs(c - c_mid) * cell_size);
        }
    }

    single_basin_fixture fx;
    fx.meta        = meta;
    fx.outlet_row  = 0;
    fx.outlet_col  = c_mid;
    fx.outlet_node = static_cast<std::int64_t>(0) * n_cols + c_mid;
    fx.dem = lsdtt::lsdtt_raster_cpu(std::move(z), std::move(mask), meta);
    return fx;
}

// Uniform-value raster on the fixture's grid (same mask/geometry as the DEM).
lsdtt::lsdtt_raster_cpu make_uniform_raster(const single_basin_fixture& fx,
                                            float value)
{
    const auto& meta = fx.meta;
    return lsdtt::lsdtt_raster_cpu(
        lsdtt::lsdtt_array2d<float>(meta.n_rows, meta.n_cols, value),
        lsdtt::lsdtt_array2d<std::uint8_t>(fx.dem.get_mask_array()),
        meta);
}

// Two-rate erosion raster: rows < split_row get value_upstream... note row 0
// is the OUTLET end, so "upstream" is the high-row (high-elevation) half.
// Used by the flux-weighting and step-change tests (plan E8).
lsdtt::lsdtt_raster_cpu make_two_rate_raster(const single_basin_fixture& fx,
                                             int split_row,
                                             float value_downstream,
                                             float value_upstream)
{
    const auto& meta = fx.meta;
    lsdtt::lsdtt_array2d<float> v(meta.n_rows, meta.n_cols, value_downstream);
    for (int r = split_row; r < meta.n_rows; ++r)
        for (int c = 0; c < meta.n_cols; ++c)
            v(r, c) = value_upstream;
    return lsdtt::lsdtt_raster_cpu(
        std::move(v),
        lsdtt::lsdtt_array2d<std::uint8_t>(fx.dem.get_mask_array()),
        meta);
}

// =============================================================================
// CSV emitters (plan Sec7) - long format, one row per evaluation
// =============================================================================
class csv_writer
{
public:
    csv_writer(const std::string& path, const std::string& header)
        : out_(path)
    {
        ok_ = out_.good();
        if (ok_) out_ << header << "\n";
        out_ << std::setprecision(17);
    }
    template <typename... Ts>
    void row(const Ts&... fields)
    {
        bool first = true;
        ((out_ << (first ? "" : ","), out_ << fields, first = false), ...);
        out_ << "\n";
    }
    bool good() const { return ok_ && out_.good(); }

private:
    std::ofstream out_;
    bool ok_ = false;
};

// =============================================================================
// Phase 1: S2 exact parameter pins (plan Sec3b)
// =============================================================================
//
// Every constant is compared with EXACT equality.  A changed constant is a
// deliberate act and must fail loudly; see the warning at the top of this
// file.  Literature anchors: Braucher = CAIRN Table 2; BraucherBorchers =
// knickpoint paper Table 1 (P0 = Borchers spallation / F0; 14C lambda updated
// to the IAEA 2024 5700 yr half-life); newCRONUS = CAIRN Table 4 (Marrero);
// Granger / Schaller have no paper table - frozen code values, COSMOCALC /
// Schaller et al. (2009, 2016) lineage.
void phase1_parameter_pins()
{
    struct scheme_pin
    {
        const char* scheme;
        double lambda_10Be, lambda_26Al, lambda_14C;
        double P0_10Be, P0_26Al, P0_14C;
        double Gamma[4];
        double F_10Be[4], F_26Al[4], F_14C[4];
    };
    const scheme_pin pins[] = {
        {"Braucher",
         500e-9, 980e-9, 121e-6,
         4.30, 31.10, 15.21,
         {160.0, 1500.0, 1500.0, 4320.0},
         {0.9887, 0.0027, 0.0, 0.0086},
         {0.9699, 0.0275, 0.0, 0.0026},
         {0.83, 0.15, 0.0, 0.02}},
        {"BraucherBorchers",
         500e-9, 980e-9, 121.58e-6,
         4.061, 28.851, 15.533,
         {160.0, 1500.0, 1500.0, 4320.0},
         {0.9874, 0.0030, 0.0, 0.0096},
         {0.9681, 0.0291, 0.0, 0.0028},
         {0.788, 0.212, 0.0, 0.0}},
        {"Granger",
         500e-9, 980e-9, 121e-6,
         5.11, 30.31, 5.86,
         {160.0, 738.6, 2688.0, 4360.0},
         {0.9724, 0.0186, 0.004, 0.005},
         {0.9655, 0.0233, 0.005, 0.0062},
         {0.83, 0.0691, 0.0809, 0.02}},
        {"Schaller",
         500e-9, 980e-9, 121e-6,
         4.01, 28.84, 12.1,
         {160.0, 1500.0, 1500.0, 4320.0},
         {0.9884, 0.0042, 0.0, 0.0074},
         {0.9697, 0.0270, 0.0, 0.0033},
         {0.83, 0.15, 0.0, 0.02}},
        {"newCRONUS",
         500e-9, 980e-9, 121e-6,
         4.075213, 27.9, 12.2,
         {160.0, 1459.76761923, 11039.2402217, 4320.0},
         {0.98374, 0.0137188126531, 0.00252519164093, 0.0},
         {0.9699, 0.0275, 0.0, 0.0026},
         {0.83, 0.15, 0.0, 0.02}},
    };

    for (const scheme_pin& pin : pins)
    {
        lsdtt::lsdtt_crn_parameters cp;
        cp.reset_scaling(pin.scheme);
        bool ok = cp.get_lambda_10Be() == pin.lambda_10Be
               && cp.get_lambda_26Al() == pin.lambda_26Al
               && cp.get_lambda_14C()  == pin.lambda_14C
               && cp.get_P0_10Be()     == pin.P0_10Be
               && cp.get_P0_26Al()     == pin.P0_26Al
               && cp.get_P0_14C()      == pin.P0_14C;
        for (int j = 0; j < 4; ++j)
        {
            ok = ok && cp.get_Gamma(j)  == pin.Gamma[j]
                    && cp.get_F_10Be(j) == pin.F_10Be[j]
                    && cp.get_F_26Al(j) == pin.F_26Al[j]
                    && cp.get_F_14C(j)  == pin.F_14C[j];
        }
        check(std::string("pin: scheme constants exact [") + pin.scheme + "]",
              ok, "a physics constant changed - see Sec3b warning");
    }

    // Documented F-sum facts (do NOT "fix" these): BraucherBorchers F sums to
    // exactly 1 for all three nuclides (knickpoint Table 1); Braucher 26Al
    // sums to 0.9699+0.0275+0.0026 = 1.0 in code (NOTE: CAIRN Table 2 prints
    // 0.00275/0.0026 summing to 0.97525 - the code uses 0.0275).
    {
        lsdtt::lsdtt_crn_parameters cp;
        cp.reset_scaling("BraucherBorchers");
        double s10 = 0, s26 = 0, s14 = 0;
        for (int j = 0; j < 4; ++j)
        {
            s10 += cp.get_F_10Be(j);
            s26 += cp.get_F_26Al(j);
            s14 += cp.get_F_14C(j);
        }
        check("pin: BraucherBorchers F sums == 1 (all nuclides)",
              std::abs(s10 - 1.0) < 1e-12 && std::abs(s26 - 1.0) < 1e-12
                  && std::abs(s14 - 1.0) < 1e-12,
              "sums: " + std::to_string(s10) + ", " + std::to_string(s26)
                  + ", " + std::to_string(s14));
    }
}

// E9: AMS standardisation factors - exact pins against the CRONUS conversion
// table (make_al_be_consts / riversand params).  Be10 reference standard is
// 07KNSTD; Al26 reference is KNSTD.
void phase1_standardisation_pins()
{
    lsdtt::lsdtt_crn_parameters cp;

    const std::pair<const char*, double> be10_pins[] = {
        {"07KNSTD", 1.0},      {"KNSTD", 0.9042},   {"NIST_Certified", 1.0425},
        {"LLNL31000", 0.8761}, {"LLNL10000", 0.9042}, {"LLNL3000", 0.8644},
        {"LLNL1000", 0.9313},  {"LLNL300", 0.8562},  {"NIST_30000", 0.9313},
        {"NIST_30200", 0.9251},{"NIST_30300", 0.9221},{"NIST_30600", 0.9130},
        {"NIST_27900", 1.0},   {"S555", 0.9124},     {"S2007", 0.9124},
        {"BEST433", 0.9124},   {"BEST433N", 1.0},    {"S555N", 1.0},
        {"S2007N", 1.0},
    };
    bool be_ok = true;
    std::string be_bad;
    for (const auto& [name, factor] : be10_pins)
    {
        if (cp.be10_standardisation_factor(name) != factor)
        {
            be_ok = false;
            be_bad = name;
        }
    }
    check("pin: Be10 standardisation factors exact (19 standards)", be_ok,
          "first mismatch: " + be_bad);

    const std::pair<const char*, double> al26_pins[] = {
        {"KNSTD", 1.0},  {"ZAL94", 0.9134}, {"SMAL11", 1.021}, {"0", 1.0},
        {"ZAL94N", 1.0}, {"ASTER", 1.021},  {"Z92-0222", 1.0},
    };
    bool al_ok = true;
    std::string al_bad;
    for (const auto& [name, factor] : al26_pins)
    {
        if (cp.al26_standardisation_factor(name) != factor)
        {
            al_ok = false;
            al_bad = name;
        }
    }
    check("pin: Al26 standardisation factors exact (7 standards)", al_ok,
          "first mismatch: " + al_bad);

    // apply_standardisation multiplies concentrations by the factor; empty
    // standard name means "already on the reference standard" (factor 1).
    const std::vector<double> conc = {100000.0, 200000.0};
    const std::vector<std::string> names = {"KNSTD", ""};
    const std::vector<double> out = cp.apply_standardisation(conc, names, "Be10");
    check("E9: apply_standardisation scales correctly",
          out.size() == 2 && approx_rel(out[0], 100000.0 * 0.9042, 1e-15)
              && out[1] == 200000.0,
          out.size() == 2 ? diff_str(out[0], 100000.0 * 0.9042) : "bad size");
}

// =============================================================================
// Phase 1: equation-fidelity tests E1-E5 (plan Sec3)
// =============================================================================

// E1 - depth-integration limits: Eq. (14) -> Eq. (12) as d_t -> 0 and
// Eq. (8) -> Eq. (7) as the slab shrinks, with first-order convergence.
void phase1_E1_depth_integration_limits()
{
    const ref_params p = ref_params_from_scheme("BraucherBorchers", "Be10");
    const double eps = 0.01;

    const double c_surface = ref_conc_surface(p, eps, 0.0);
    double prev_err = -1.0;
    bool converges = true;
    for (double d_t : {1.0, 0.1, 0.01})
    {
        const double err =
            std::abs(ref_conc_depth_integrated(p, eps, 0.0, d_t) - c_surface)
            / c_surface;
        if (prev_err > 0.0 && err > 0.2 * prev_err) converges = false;
        prev_err = err;
    }
    check("E1: Eq14 -> Eq12 first-order in d_t", converges,
          "final rel err = " + std::to_string(prev_err));

    const double t = 20000.0, eps_o = 0.005, eps_n = 0.05, d_mid = 30.0;
    const double c_point = ref_conc_step_change(p, d_mid, t, eps_o, eps_n);
    prev_err = -1.0;
    converges = true;
    for (double half : {1.0, 0.1, 0.01})
    {
        const double err =
            std::abs(ref_conc_step_change_depth_integrated(
                         p, d_mid - half, d_mid + half, t, eps_o, eps_n)
                     - c_point) / c_point;
        if (prev_err > 0.0 && err > 0.2 * prev_err) converges = false;
        prev_err = err;
    }
    check("E1: Eq8 -> Eq7 as slab shrinks", converges,
          "final rel err = " + std::to_string(prev_err));
}

// E2 - step-change end members, LIBRARY implementation
// (calculate_crn_concentration_raster_step_change).  Raster storage is float,
// so library-to-library comparisons use a 1e-5 relative tolerance.
void phase1_E2_step_change_library(const single_basin_fixture& fx)
{
    const double eps_o = 0.005, eps_n = 0.05;
    const lsdtt::lsdtt_raster_cpu prod  = make_uniform_raster(fx, 0.85f);
    const lsdtt::lsdtt_raster_cpu topo  = make_uniform_raster(fx, 1.0f);
    const lsdtt::lsdtt_raster_cpu zero  = make_uniform_raster(fx, 0.0f);
    const std::size_t idx = static_cast<std::size_t>(fx.outlet_node);

    auto ss_conc = [&](double eps) {
        return lsdtt::calculate_crn_concentration_raster_whole(
                   "Be10", "BraucherBorchers",
                   make_uniform_raster(fx, static_cast<float>(eps)),
                   prod, topo, zero, zero).get_data()[idx];
    };
    auto step_conc = [&](double t, double eo, double en) {
        return lsdtt::calculate_crn_concentration_raster_step_change(
                   "Be10", "BraucherBorchers",
                   make_uniform_raster(fx, static_cast<float>(t)),
                   prod, topo, eo, en).get_data()[idx];
    };

    const double c_ss_old = ss_conc(eps_o);
    const double c_ss_new = ss_conc(eps_n);

    // Documented behaviour: time_since_change <= 0 returns SS at eps_new.
    check("E2: library step-change t=0 == SS(eps_new) [documented]",
          approx_rel(step_conc(0.0, eps_o, eps_n), c_ss_new, 1e-5),
          diff_str(step_conc(0.0, eps_o, eps_n), c_ss_new));

    // t -> infinity relaxes to SS at eps_new.
    check("E2: library step-change t->inf -> SS(eps_new)",
          approx_rel(step_conc(5.0e7, eps_o, eps_n), c_ss_new, 1e-4),
          diff_str(step_conc(5.0e7, eps_o, eps_n), c_ss_new));

    // eps_old == eps_new: SS at that rate for any t.
    check("E2: library step-change eps_o==eps_n == SS",
          approx_rel(step_conc(23456.0, eps_n, eps_n), c_ss_new, 1e-4),
          diff_str(step_conc(23456.0, eps_n, eps_n), c_ss_new));

    // Trajectory: monotone from SS(old) toward SS(new), never outside the
    // [SS(new), SS(old)] envelope (eps_n > eps_o so concentration falls).
    bool inside = true, monotone = true;
    double c_prev = c_ss_old;
    const double slack = 1e-4;   // float-storage slack on the envelope
    for (double t : {1000.0, 5000.0, 20000.0, 80000.0, 320000.0})
    {
        const double c = step_conc(t, eps_o, eps_n);
        if (c > c_ss_old * (1.0 + slack) || c < c_ss_new * (1.0 - slack))
            inside = false;
        if (c > c_prev * (1.0 + slack)) monotone = false;
        c_prev = c;
    }
    check("E2: library step-change trajectory inside SS envelope", inside);
    check("E2: library step-change trajectory monotone", monotone);
}

// E3 - the exponential identity in Eq. (7): both algebraic groupings agree
// on the domain where the printed form does not overflow.
void phase1_E3_exponential_identity()
{
    const ref_params p = ref_params_from_scheme("BraucherBorchers", "Be10");
    bool ok = true;
    std::string bad;
    for (double t : {0.0, 1000.0, 50000.0, 1.0e6})
    {
        for (double eps_n : {0.001, 0.05})
        {
            for (int j = 0; j < 4; ++j)
            {
                const double G = p.Gamma[j];
                if (t * eps_n / G > 600.0) continue;   // overflow domain
                const double d = 15.0, d0 = d + t * eps_n;
                const double printed =
                    std::exp(-d0 / G)
                    * (std::exp(t * eps_n / G) - std::exp(-t * p.lambda));
                const double stable =
                    std::exp(-d / G)
                    - std::exp(-d0 / G) * std::exp(-t * p.lambda);
                if (!approx_rel(printed, stable, 1e-12))
                {
                    ok = false;
                    bad = diff_str(printed, stable);
                }
            }
        }
    }
    check("E3: Eq7 exponential identity (printed vs stable grouping)", ok, bad);
}

// E4 - the Vermeesch virtual-attenuation-length redistribution
// (scale_F_values / the Lambda_v Newton solve, CAIRN Eqs. 15-16).
void phase1_E4_lambda_v_solve()
{
    for (const char* scheme : {"BraucherBorchers", "Braucher", "newCRONUS"})
    {
        for (double s_tot : {0.35, 0.7, 0.95})
        {
            lsdtt::lsdtt_crn_parameters cp;
            cp.reset_scaling(scheme);
            cp.scale_F_values(s_tot, "Be10");
            double sum = 0.0;
            for (int j = 0; j < 4; ++j) sum += cp.get_F_10Be(j);
            check(std::string("E4: sum(F_scaled) == S_tot [") + scheme + ", "
                      + std::to_string(s_tot) + "]",
                  approx_rel(sum, s_tot, 1e-6), diff_str(sum, s_tot));
        }
    }

    // S_tot = 1 with an F-sum-1 scheme leaves F unchanged (Lambda_v = 0).
    {
        lsdtt::lsdtt_crn_parameters cp;
        cp.reset_scaling("BraucherBorchers");
        double f_before[4];
        for (int j = 0; j < 4; ++j) f_before[j] = cp.get_F_10Be(j);
        cp.scale_F_values(1.0, "Be10");
        bool unchanged = true;
        for (int j = 0; j < 4; ++j)
            if (!approx_rel(cp.get_F_10Be(j), f_before[j], 1e-9))
                unchanged = false;
        check("E4: S_tot=1 leaves F unchanged (Lambda_v=0)", unchanged);
    }

    // Attenuation ordering: for S_tot < 1 the effective depth x is positive,
    // so the pathway with the smallest Gamma (spallation) is attenuated most:
    // F0_scaled/F0 < F1_scaled/F1 <= F3_scaled/F3 (Gamma 160 < 1500 < 4320).
    {
        lsdtt::lsdtt_crn_parameters cp;
        cp.reset_scaling("BraucherBorchers");
        const double f0 = cp.get_F_10Be(0), f1 = cp.get_F_10Be(1),
                     f3 = cp.get_F_10Be(3);
        cp.scale_F_values(0.5, "Be10");
        const double r0 = cp.get_F_10Be(0) / f0;
        const double r1 = cp.get_F_10Be(1) / f1;
        const double r3 = cp.get_F_10Be(3) / f3;
        check("E4: spallation attenuated most (r0 < r1 < r3)",
              r0 < r1 && r1 < r3,
              "ratios: " + std::to_string(r0) + ", " + std::to_string(r1)
                  + ", " + std::to_string(r3));
    }

    // S_tot > 1 (production scaling at altitude): the Vermeesch solve must
    // amplify with a NEGATIVE effective depth.  Sum identity holds, and the
    // ordering reverses: spallation (smallest Gamma) amplified most.  This
    // regime is what the >= 1 clamp silently skipped (findings log Sec9b,
    // fix approved 2026-07-16); values span sea-level-ish to the Stone
    // factor at 60 deg / 3000 m (~11).
    for (double s_tot : {1.8, 4.5, 11.0})
    {
        lsdtt::lsdtt_crn_parameters cp;
        cp.reset_scaling("BraucherBorchers");
        const double f0 = cp.get_F_10Be(0), f1 = cp.get_F_10Be(1),
                     f3 = cp.get_F_10Be(3);
        cp.scale_F_values(s_tot, "Be10");
        double sum = 0.0;
        for (int j = 0; j < 4; ++j) sum += cp.get_F_10Be(j);
        check("E4: sum(F_scaled) == S_tot for S_tot>1 ["
                  + std::to_string(s_tot) + "]",
              approx_rel(sum, s_tot, 1e-6), diff_str(sum, s_tot));

        const double r0 = cp.get_F_10Be(0) / f0;
        const double r1 = cp.get_F_10Be(1) / f1;
        const double r3 = cp.get_F_10Be(3) / f3;
        check("E4: spallation amplified most for S_tot>1 ["
                  + std::to_string(s_tot) + "]",
              r0 > r1 && r1 > r3 && r3 > 1.0,
              "ratios: " + std::to_string(r0) + ", " + std::to_string(r1)
                  + ", " + std::to_string(r3));
    }
}

// E5 - linearity and asymptotic limits of the reference steady state.
void phase1_E5_linearity_and_limits()
{
    ref_params p = ref_params_from_scheme("BraucherBorchers", "Be10");
    const double eps = 0.01;

    // Linear in P0 (exact).
    {
        ref_params p2 = p;
        p2.P0 *= 2.0;
        check("E5: C linear in P0",
              approx_rel(ref_conc_surface(p2, eps, 0.0),
                         2.0 * ref_conc_surface(p, eps, 0.0), 1e-15));
    }

    // Saturation: eps -> 0 gives C -> (P0/lambda) * sum(S_j F_j).
    {
        double f_sum = 0.0;
        for (int j = 0; j < 4; ++j) f_sum += p.S[j] * p.F[j];
        const double c_sat = p.P0 / p.lambda * f_sum;
        check("E5: eps->0 saturation limit",
              approx_rel(ref_conc_surface(p, 0.0, 0.0), c_sat, 1e-12),
              diff_str(ref_conc_surface(p, 0.0, 0.0), c_sat));
    }

    // Fast erosion: eps >> Gamma0*lambda gives C*eps -> P0 * sum(S F Gamma).
    {
        double sfg = 0.0;
        for (int j = 0; j < 4; ++j) sfg += p.S[j] * p.F[j] * p.Gamma[j];
        const double eps_fast = 100.0;   // g cm^-2 yr^-1, >> 160*5e-7
        check("E5: fast-erosion 1/eps limit",
              approx_rel(ref_conc_surface(p, eps_fast, 0.0) * eps_fast,
                         p.P0 * sfg, 1e-4),
              diff_str(ref_conc_surface(p, eps_fast, 0.0) * eps_fast,
                       p.P0 * sfg));
    }
}

// =============================================================================
// Phase 2: E6-E8 (plan Sec3)
// =============================================================================

// Builds ref_params whose F values have been scaled by the library's own
// Lambda_v redistribution for total shielding s_tot, so the reference and the
// particle path see identical inputs and any disagreement is in the
// concentration math itself.
ref_params ref_params_scaled(const std::string& scheme,
                             const std::string& nuclide, double s_tot)
{
    lsdtt::lsdtt_crn_parameters cp;
    cp.reset_scaling(scheme);
    cp.scale_F_values(s_tot, nuclide);
    ref_params p = ref_params_from_scheme(scheme, nuclide);
    for (int j = 0; j < 4; ++j)
    {
        p.F[j] = (nuclide == "Be10") ? cp.get_F_10Be(j)
                 : (nuclide == "Al26") ? cp.get_F_26Al(j)
                                       : cp.get_F_14C(j);
        p.S[j] = 1.0;
    }
    return p;
}

// E7 - analytic reference (double precision, this file) vs the library's
// particle integrator (compute_ss via calculate_crn_concentration_raster_whole),
// point and depth-integrated, shielded and unshielded, across an eps sweep.
// Library rasters store float, so the agreement floor is ~1e-7 relative.
// Measured worst case on first run (2026-07-16): 9.25e-8, i.e. the two paths
// are algebraically identical and the difference is pure float storage.
// Bound frozen at 5e-7 per plan open question 2 (measure first, then freeze).
void phase2_E7_analytic_vs_integrator(const single_basin_fixture& fx,
                                      const std::string& out_dir)
{
    const std::size_t idx = static_cast<std::size_t>(fx.outlet_node);
    csv_writer csv(out_dir + "/analytic_vs_integrator.csv",
                   "nuclide,scheme,erate_gcm2yr,s_tot,snow_gcm2,self_gcm2,"
                   "conc_reference,conc_library,rel_diff");

    double worst = 0.0;
    std::string worst_case;
    for (const char* scheme : {"BraucherBorchers", "newCRONUS"})
    {
        for (const char* nuc : {"Be10", "Al26", "C14"})
        {
            for (double s_tot : {1.0, 0.85})
            {
                for (double snow : {0.0, 10.0})
                {
                    for (double self : {0.0, 20.0})
                    {
                        for (double eps : {0.0005, 0.005, 0.05})
                        {
                            const auto conc =
                                lsdtt::calculate_crn_concentration_raster_whole(
                                    nuc, scheme,
                                    make_uniform_raster(fx, static_cast<float>(eps)),
                                    make_uniform_raster(fx, static_cast<float>(s_tot)),
                                    make_uniform_raster(fx, 1.0f),
                                    make_uniform_raster(fx, static_cast<float>(snow)),
                                    make_uniform_raster(fx, static_cast<float>(self)));
                            const double c_lib = conc.get_data()[idx];

                            const ref_params p =
                                ref_params_scaled(scheme, nuc, s_tot);
                            const double c_ref =
                                (self == 0.0)
                                    ? ref_conc_surface(p, eps, snow)
                                    : ref_conc_depth_integrated(p, eps, snow,
                                                                self);

                            const double rd =
                                std::abs(c_lib - c_ref)
                                / std::max(std::abs(c_ref), 1e-30);
                            csv.row(nuc, scheme, eps, s_tot, snow, self,
                                    c_ref, c_lib, rd);
                            if (rd > worst)
                            {
                                worst = rd;
                                std::ostringstream os;
                                os << nuc << "/" << scheme << " eps=" << eps
                                   << " s=" << s_tot << " snow=" << snow
                                   << " self=" << self;
                                worst_case = os.str();
                            }
                        }
                    }
                }
            }
        }
    }
    std::cout << "      E7 measured worst rel diff = "
              << std::setprecision(3) << worst << "  (" << worst_case << ")\n"
              << std::setprecision(17);
    check("E7: analytic vs particle integrator agree (< 5e-7, float floor)",
          worst < tol::e7_analytic, "worst = " + std::to_string(worst));
    check("csv: analytic_vs_integrator.csv written", csv.good());
}

// E6 - forward/inverse round-trip on the fixture basin with the REAL Stone
// production raster (per-pixel variation), for a sweep of imposed erosion
// rates, nuclides/schemes, with and without shielding.  Measured worst
// recovery error on first run (2026-07-16): 9.93e-4 - that is the inversion
// solver's own convergence tolerance, not a physics error (E7 shows the
// forward paths agree to 1e-7).  Bound frozen at 3e-3 for headroom.
void phase2_E6_round_trip(const single_basin_fixture& fx,
                          const lsdtt::lsdtt_flowinfo& flow)
{
    const lsdtt::lsdtt_raster_cpu prod =
        lsdtt::calculate_production_raster(fx.dem);

    struct rt_case
    {
        const char* nuclide;
        const char* scheme;
        float topo, snow, self;
    };
    const rt_case cases[] = {
        {"Be10", "BraucherBorchers", 1.0f, 0.0f, 0.0f},
        {"Be10", "BraucherBorchers", 0.85f, 8.0f, 12.0f},
        {"C14",  "BraucherBorchers", 1.0f, 0.0f, 0.0f},
        {"Al26", "newCRONUS",        1.0f, 0.0f, 0.0f},
    };

    double worst = 0.0;
    for (const rt_case& rc : cases)
    {
        const lsdtt::lsdtt_raster_cpu topo = make_uniform_raster(fx, rc.topo);
        const lsdtt::lsdtt_raster_cpu snow = make_uniform_raster(fx, rc.snow);
        const lsdtt::lsdtt_raster_cpu self = make_uniform_raster(fx, rc.self);

        for (double eps : {0.001, 0.01, 0.08})
        {
            const lsdtt::lsdtt_raster_cpu ero =
                make_uniform_raster(fx, static_cast<float>(eps));
            const auto conc = lsdtt::calculate_crn_concentration_raster(
                rc.nuclide, rc.scheme, ero, prod, topo, snow, self,
                flow, fx.outlet_node);
            const double c_meas = lsdtt::calculate_accumulated_crn_concentration(
                conc, ero, flow, fx.outlet_node);

            double c_pred = 0.0;
            const double eps_out = lsdtt::calculate_eff_erate_from_conc(
                c_meas, rc.nuclide, rc.scheme, prod, topo, snow, self,
                flow, fx.outlet_node, c_pred);

            const double rd = std::abs(eps_out - eps) / eps;
            worst = std::max(worst, rd);
            check(std::string("E6: round-trip recovers eps [") + rc.nuclide
                      + "/" + rc.scheme + " eps=" + std::to_string(eps)
                      + (rc.topo != 1.0f ? " shielded" : "") + "]",
                  rd < tol::e6_roundtrip,
                  "recovered " + std::to_string(eps_out) + ", rel err "
                      + std::to_string(rd));
        }
    }
    std::cout << "      E6 measured worst round-trip rel err = "
              << std::setprecision(3) << worst << std::setprecision(17)
              << "\n";

    // Edge cases: the inversion must not diverge, hang, or return non-finite
    // values for pathological concentrations.
    {
        const lsdtt::lsdtt_raster_cpu topo = make_uniform_raster(fx, 1.0f);
        const lsdtt::lsdtt_raster_cpu zero = make_uniform_raster(fx, 0.0f);

        // Above saturation: higher concentration than eps -> 0 can produce.
        const lsdtt::lsdtt_raster_cpu ero_tiny = make_uniform_raster(fx, 1e-8f);
        const auto conc_sat = lsdtt::calculate_crn_concentration_raster(
            "Be10", "BraucherBorchers", ero_tiny, prod, topo, zero, zero,
            flow, fx.outlet_node);
        const double c_sat = lsdtt::calculate_accumulated_crn_concentration(
            conc_sat, ero_tiny, flow, fx.outlet_node);
        bool ok = true;
        double eps_out = 0.0, c_pred = 0.0;
        try
        {
            eps_out = lsdtt::calculate_eff_erate_from_conc(
                2.0 * c_sat, "Be10", "BraucherBorchers", prod, topo, zero,
                zero, flow, fx.outlet_node, c_pred);
        }
        catch (const std::exception&) { ok = false; }
        check("E6: above-saturation concentration handled (finite, no throw)",
              ok && std::isfinite(eps_out) && eps_out >= 0.0,
              "eps_out = " + std::to_string(eps_out));

        // Very high implied erosion rate (tiny concentration).
        ok = true;
        try
        {
            eps_out = lsdtt::calculate_eff_erate_from_conc(
                50.0, "Be10", "BraucherBorchers", prod, topo, zero, zero,
                flow, fx.outlet_node, c_pred);
        }
        catch (const std::exception&) { ok = false; }
        check("E6: tiny concentration (huge erate) handled",
              ok && std::isfinite(eps_out) && eps_out > 0.0,
              "eps_out = " + std::to_string(eps_out));
    }
}

// E8 - averaging identities of the accumulation (prediction) path:
// flux-weighted mean sum(C*eps)/sum(eps), which reduces to the plain mean for
// uniform eps and must differ from it measurably on a two-rate basin.
void phase2_E8_averaging_identities(const single_basin_fixture& fx,
                                    const lsdtt::lsdtt_flowinfo& flow)
{
    const lsdtt::lsdtt_raster_cpu prod = make_uniform_raster(fx, 0.85f);
    const lsdtt::lsdtt_raster_cpu topo = make_uniform_raster(fx, 1.0f);
    const lsdtt::lsdtt_raster_cpu zero = make_uniform_raster(fx, 0.0f);
    const auto& mask = fx.dem.get_mask();

    // Uniform eps: accumulated concentration == plain mean over the basin.
    // (The conc raster varies per pixel here?  With uniform prod it does not,
    // so instead use the REAL production raster to get spatial variation and
    // make the mean identity non-trivial.)
    const lsdtt::lsdtt_raster_cpu prod_real =
        lsdtt::calculate_production_raster(fx.dem);
    {
        const double eps = 0.01;
        const lsdtt::lsdtt_raster_cpu ero =
            make_uniform_raster(fx, static_cast<float>(eps));
        const auto conc = lsdtt::calculate_crn_concentration_raster(
            "Be10", "BraucherBorchers", ero, prod_real, topo, zero, zero,
            flow, fx.outlet_node);
        const double acc = lsdtt::calculate_accumulated_crn_concentration(
            conc, ero, flow, fx.outlet_node);

        const auto& cd = conc.get_data();
        double sum = 0.0;
        std::int64_t n = 0;
        for (std::size_t i = 0; i < cd.size(); ++i)
        {
            if (mask[i]) { sum += cd[i]; ++n; }
        }
        const double plain_mean = sum / static_cast<double>(n);
        check("E8: uniform eps -> accumulation == plain mean",
              approx_rel(acc, plain_mean, 1e-6), diff_str(acc, plain_mean));
    }

    // Two-rate basin: accumulated == sum(C*eps)/sum(eps) exactly, and it
    // must differ measurably from the plain mean (flux weighting matters).
    {
        const lsdtt::lsdtt_raster_cpu ero =
            make_two_rate_raster(fx, fx.meta.n_rows / 2, 0.005f, 0.05f);
        const auto conc = lsdtt::calculate_crn_concentration_raster(
            "Be10", "BraucherBorchers", ero, prod_real, topo, zero, zero,
            flow, fx.outlet_node);
        const double acc = lsdtt::calculate_accumulated_crn_concentration(
            conc, ero, flow, fx.outlet_node);

        const auto& cd = conc.get_data();
        const auto& ed = ero.get_data();
        double sum_ce = 0.0, sum_e = 0.0, sum_c = 0.0;
        std::int64_t n = 0;
        for (std::size_t i = 0; i < cd.size(); ++i)
        {
            if (!mask[i]) continue;
            sum_ce += static_cast<double>(cd[i]) * ed[i];
            sum_e  += ed[i];
            sum_c  += cd[i];
            ++n;
        }
        const double flux_weighted = sum_ce / sum_e;
        const double plain_mean    = sum_c / static_cast<double>(n);

        check("E8: two-rate basin -> accumulation == flux-weighted mean",
              approx_rel(acc, flux_weighted, 1e-6),
              diff_str(acc, flux_weighted));
        check("E8: flux weighting differs from plain mean on two-rate basin",
              std::abs(flux_weighted - plain_mean) / plain_mean > 1e-3,
              diff_str(flux_weighted, plain_mean));

        // Consistency: the whole-raster accumulation at the outlet equals the
        // scalar accumulated concentration.
        const auto acc_raster = lsdtt::calculate_accumulated_crn_raster(
            conc, ero, flow);
        const double acc_at_outlet =
            acc_raster.get_data()[static_cast<std::size_t>(fx.outlet_node)];
        check("E8: accumulated raster outlet == scalar accumulation",
              approx_rel(acc_at_outlet, acc, 1e-5),
              diff_str(acc_at_outlet, acc));
    }
}

// =============================================================================
// Phase 3: CRONUS comparison fixtures and tests T-C1..T-C4 (plan Sec5)
// =============================================================================

// T-C3 - Stone (2000) scaling factors pinned against an INDEPENDENT
// evaluation of Stone's Table 1 spallation coefficients (computed outside
// this codebase, 2026-07-16; also what CRONUS v2.3 stone2000.m implements).
// Latitudes >= 60 use the 60-degree row (Stone 2000; verified to match).
void phase3_TC3_stone_pins()
{
    struct stone_pin { double lat, p_mbar, expected; };
    const stone_pin pins[] = {
        {0.0,  1013.25, 0.5872729071401146},
        {0.0,  700.0,   4.5358221953964959},
        {30.0, 1013.25, 0.83002237518144995},
        {30.0, 700.0,   6.8069850492071069},
        {60.0, 1013.25, 1.0021331925734103},
        {60.0, 700.0,   11.011637548444043},
        {70.0, 1013.25, 1.0021331925734103},   // >=60 clamps to the 60 row
        {70.0, 700.0,   11.011637548444043},
    };
    lsdtt::lsdtt_crn_parameters cp;
    bool ok = true;
    std::string bad;
    for (const stone_pin& s : pins)
    {
        const double got = cp.stone2000sp(s.lat, s.p_mbar, 1.0);
        if (!approx_rel(got, s.expected, 1e-10))
        {
            ok = false;
            bad = "lat=" + std::to_string(s.lat) + " P="
                  + std::to_string(s.p_mbar) + ": " + diff_str(got, s.expected);
        }
    }
    check("T-C3: stone2000sp matches independent Stone (2000) Table 1", ok, bad);
}

// T-C1 - forward concentration curves at the CAIRN Figs. 8-9 hypothetical
// site (70 N, 1007 hPa, Lal/Stone, no shielding), eps swept across the
// decay-dominated to erosion-dominated regimes.  Plottable analogue of
// CAIRN Fig. 9a; CRONUS forward concentrations can be merged into the same
// CSV externally for the Fig. 9b percent-difference panel.
void phase3_TC1_site_curve(const std::string& out_dir)
{
    lsdtt::lsdtt_crn_parameters cp;
    const double s_site = cp.stone2000sp(70.0, 1007.0, 1.0);

    csv_writer csv(out_dir + "/conc_vs_erate_70N.csv",
                   "nuclide,scheme,lat,pressure_hPa,stone_factor,"
                   "erate_gcm2yr,conc_atoms_g,source");
    bool ok = true;
    for (const char* scheme : {"Braucher", "BraucherBorchers", "newCRONUS"})
    {
        const ref_params p = ref_params_scaled(scheme, "Be10", s_site);
        double c_prev = -1.0;
        for (double eps = 0.0005; eps <= 0.2000001; eps += 0.0005)
        {
            const double c = ref_conc_surface(p, eps, 0.0);
            if (!std::isfinite(c) || c <= 0.0) ok = false;
            if (c_prev > 0.0 && c >= c_prev) ok = false;   // must fall with eps
            csv.row("Be10", scheme, 70.0, 1007.0, s_site, eps, c, "lsdtt");
            c_prev = c;
        }
    }
    check("T-C1: 70N/1007hPa forward curves finite and monotone",
          ok && csv.good());
}

// The CRONUS comparison fixture: synthetic point samples spanning latitude,
// elevation, and erosion rate.  Concentrations are OUR forward model (Stone
// scaling at the standard-atmosphere pressure, BraucherBorchers, surface,
// no shielding), reported on the reference standards (07KNSTD / KNSTD), so
// submitting them to CRONUS and comparing inverted rates against the imposed
// eps measures the total inter-calculator difference (plan Sec5).
struct cronus_fixture_sample
{
    std::string name;
    double lat, lon, elev_m, pressure_hPa, stone_factor;
    double erate_gcm2yr;
    double conc_be10, conc_al26;
};

std::vector<cronus_fixture_sample> make_cronus_fixture_samples()
{
    std::vector<cronus_fixture_sample> out;
    const double lats[]  = {0.0, 20.0, 40.0, 60.0, 70.0};
    const double elevs[] = {250.0, 1500.0, 3000.0};
    // Erosion rates cycle log-spaced across the CAIRN-relevant range.
    const double erates[] = {3e-4, 1e-3, 3e-3, 1e-2, 3e-2};

    lsdtt::lsdtt_crn_parameters cp;
    int k = 0;
    for (double lat : lats)
    {
        for (double elev : elevs)
        {
            cronus_fixture_sample s;
            s.name         = "LSDTT-" + std::to_string(k + 1);
            s.lat          = lat;
            s.lon          = 0.0;
            s.elev_m       = elev;
            s.pressure_hPa =
                lsdtt::lsdtt_crn_parameters::atmospheric_pressure_standard(elev);
            s.stone_factor = cp.stone2000sp(lat, s.pressure_hPa, 1.0);
            s.erate_gcm2yr = erates[k % 5];

            const ref_params p10 =
                ref_params_scaled("BraucherBorchers", "Be10", s.stone_factor);
            const ref_params p26 =
                ref_params_scaled("BraucherBorchers", "Al26", s.stone_factor);
            s.conc_be10 = ref_conc_surface(p10, s.erate_gcm2yr, 0.0);
            s.conc_al26 = ref_conc_surface(p26, s.erate_gcm2yr, 0.0);
            out.push_back(s);
            ++k;
        }
    }
    return out;
}

void phase3_write_cronus_fixture(const std::string& out_dir)
{
    const auto samples = make_cronus_fixture_samples();
    csv_writer csv(out_dir + "/cronus_fixture_samples.csv",
                   "name,lat,lon,elev_m,pressure_hPa,stone_factor,"
                   "erate_gcm2yr,conc_Be10_07KNSTD,conc_Al26_KNSTD");
    for (const auto& s : samples)
    {
        csv.row(s.name, s.lat, s.lon, s.elev_m, s.pressure_hPa,
                s.stone_factor, s.erate_gcm2yr, s.conc_be10, s.conc_al26);
    }
    check("T-C2: cronus_fixture_samples.csv written ("
              + std::to_string(samples.size()) + " samples)",
          csv.good());
}

// T-C2 / T-C4 - comparison against the FROZEN CRONUS v3 reference
// (tests/expected/cosmo/cronus_v3_reference.csv, fetched once by
// tests/scripts/fetch_cronus_reference.py).  Skips with a NOTE when the
// reference has not been fetched yet; never touches the network.
void phase3_TC2_cronus_comparison(const std::string& expected_dir,
                                  const std::string& out_dir)
{
    const std::string ref_path = expected_dir + "/cronus_v3_reference.csv";
    std::ifstream in(ref_path);
    if (!in.good())
    {
        std::cout << "NOTE  T-C2 skipped: no frozen CRONUS reference at "
                  << ref_path << "\n"
                  << "      (fetch once with tests/scripts/fetch_cronus_reference.py)\n";
        check("T-C2: CRONUS comparison (SKIPPED - reference not fetched)", true);
        return;
    }

    // reference rows: name,nuclide,scheme,erate_gcm2yr,erate_mMyr,unc_int_mMyr
    std::map<std::string, double> cronus;   // key: name|nuclide|scheme
    std::string line;
    while (std::getline(in, line))
    {
        if (line.empty() || line[0] == '#' || line.rfind("name,", 0) == 0)
            continue;
        std::istringstream ls(line);
        std::string name, nuc, scheme, field;
        std::getline(ls, name, ',');
        std::getline(ls, nuc, ',');
        std::getline(ls, scheme, ',');
        std::getline(ls, field, ',');
        cronus[name + "|" + nuc + "|" + scheme] = std::stod(field);
    }

    const auto samples = make_cronus_fixture_samples();
    csv_writer csv(out_dir + "/cronus_comparison.csv",
                   "name,nuclide,scheme,lat,elev_m,erate_lsdtt_gcm2yr,"
                   "erate_cronus_gcm2yr,pct_diff");
    double worst_st = 0.0;
    int n_compared = 0;
    for (const auto& s : samples)
    {
        for (const char* scheme : {"St", "Lm", "LSDn"})
        {
            const auto it = cronus.find(s.name + "|Be-10|" + scheme);
            if (it == cronus.end()) continue;
            const double pct =
                100.0 * (it->second - s.erate_gcm2yr) / s.erate_gcm2yr;
            csv.row(s.name, "Be-10", scheme, s.lat, s.elev_m,
                    s.erate_gcm2yr, it->second, pct);
            if (std::string(scheme) == "St")
            {
                worst_st = std::max(worst_st, std::abs(pct));
                ++n_compared;
            }
        }
    }
    std::cout << "      T-C2 measured worst |pct diff| vs CRONUS v3 St = "
              << std::setprecision(4) << worst_st << " % over " << n_compared
              << " samples\n" << std::setprecision(17);
    // Band: our BraucherBorchers muons + Borchers-lineage P0 differ from
    // CRONUS v3 St (Balco 2017 model-1A muons, 07KNSTD-normalised rates) by
    // a bounded, production-dependent amount.  Measured on the frozen
    // reference (fetched 2026-07-16, post S_tot>1 fix): +1.82 % to +6.65 %,
    // mean +3.61 % (CRONUS higher, growing with erosion rate as the muon
    // share rises) - consistent with CAIRN Sec7's 2.5-5 % for
    // Borchers-consistent parameters.  Band frozen at 10 %.
    // For the record: Lm ran mean +14.6 % (max +45 %) and LSDn mean +11.8 %
    // vs our time-INdependent rates - that spread is the size of the
    // geomagnetic correction itself (baseline for the Sec6 G-tests).
    check("T-C2: vs CRONUS v3 St within band (< 10 %)",
          n_compared > 0 && worst_st < tol::tc2_cronus,
          "worst = " + std::to_string(worst_st) + " % (n=" + std::to_string(n_compared) + ")");
    check("csv: cronus_comparison.csv written", csv.good());
}

// =============================================================================
// Phase 1: S1/S3/S4 numerics sentinel (plan Sec3b)
// =============================================================================
//
// Full-double-precision pinned values for the library's numeric core.  The
// comparison tolerance (1e-12 relative) is deliberately far stricter than
// physical accuracy: the sentinel's job is to trip on ANY behavioural change.
// DO NOT loosen it or regenerate the file to make a failure go away - that is
// exactly the event the sentinel exists to gate (Sec3b S6).

constexpr double k_sentinel_rel_tol = tol::sentinel_rel;

struct sentinel_entry
{
    std::string quantity, nuclide, scheme;
    double p1 = 0.0, p2 = 0.0, p3 = 0.0;
    double value = 0.0;
};

std::string sentinel_key(const sentinel_entry& e)
{
    std::ostringstream os;
    os << std::setprecision(17) << e.quantity << "|" << e.nuclide << "|"
       << e.scheme << "|" << e.p1 << "|" << e.p2 << "|" << e.p3;
    return os.str();
}

// Computes every sentinel quantity live from the library.
std::vector<sentinel_entry> compute_sentinel_values(
    const single_basin_fixture& fx, const lsdtt::lsdtt_flowinfo& flow)
{
    std::vector<sentinel_entry> out;
    const char* nuclides[] = {"Be10", "Al26", "C14"};
    const char* schemes[]  = {"Braucher", "BraucherBorchers", "Granger",
                              "Schaller", "newCRONUS"};

    // stone2000sp scaling factors (scheme-independent).
    {
        lsdtt::lsdtt_crn_parameters cp;
        for (double lat : {0.0, 30.0, 60.0, 70.0})
        {
            for (double p_mbar : {1013.25, 700.0})
            {
                out.push_back({"stone2000sp", "-", "-", lat, p_mbar, 0.0,
                               cp.stone2000sp(lat, p_mbar, 1.0)});
            }
        }
    }

    // scale_F_values (Lambda_v Newton solve) results.
    for (const char* scheme : schemes)
    {
        for (const char* nuc : nuclides)
        {
            for (double s_tot : {0.35, 0.85})
            {
                lsdtt::lsdtt_crn_parameters cp;
                cp.reset_scaling(scheme);
                cp.scale_F_values(s_tot, nuc);
                for (int j = 0; j < 4; ++j)
                {
                    const double f =
                        (std::string(nuc) == "Be10") ? cp.get_F_10Be(j)
                        : (std::string(nuc) == "Al26") ? cp.get_F_26Al(j)
                                                       : cp.get_F_14C(j);
                    out.push_back({"scaled_F", nuc, scheme,
                                   s_tot, static_cast<double>(j), 0.0, f});
                }
            }
        }
    }

    // Whole-raster steady-state concentration (centre pixel, uniform inputs:
    // prod scaling 0.85, topo shielding 1, no snow/self).
    {
        const lsdtt::lsdtt_raster_cpu prod = make_uniform_raster(fx, 0.85f);
        const lsdtt::lsdtt_raster_cpu topo = make_uniform_raster(fx, 1.0f);
        const lsdtt::lsdtt_raster_cpu zero = make_uniform_raster(fx, 0.0f);
        const std::size_t idx = static_cast<std::size_t>(fx.outlet_node);
        for (const char* scheme : schemes)
        {
            for (const char* nuc : nuclides)
            {
                for (double eps : {0.001, 0.02})
                {
                    const auto conc =
                        lsdtt::calculate_crn_concentration_raster_whole(
                            nuc, scheme,
                            make_uniform_raster(fx, static_cast<float>(eps)),
                            prod, topo, zero, zero);
                    out.push_back({"conc_whole", nuc, scheme, eps, 0.0, 0.0,
                                   conc.get_data()[idx]});
                }
            }
        }
    }

    // Step-change concentration (BraucherBorchers, eps 0.005 -> 0.05).
    {
        const lsdtt::lsdtt_raster_cpu prod = make_uniform_raster(fx, 0.85f);
        const lsdtt::lsdtt_raster_cpu topo = make_uniform_raster(fx, 1.0f);
        const std::size_t idx = static_cast<std::size_t>(fx.outlet_node);
        const double eps_o = 0.005, eps_n = 0.05;
        for (const char* nuc : nuclides)
        {
            for (double t : {0.0, 5000.0, 50000.0})
            {
                const auto conc =
                    lsdtt::calculate_crn_concentration_raster_step_change(
                        nuc, "BraucherBorchers",
                        make_uniform_raster(fx, static_cast<float>(t)),
                        prod, topo, eps_o, eps_n);
                out.push_back({"conc_step", nuc, "BraucherBorchers",
                               t, eps_o, eps_n, conc.get_data()[idx]});
            }
        }
    }

    // Basin inversion on the fixture (production from the real Stone raster,
    // so this also pins calculate_production_raster along the way).
    {
        const lsdtt::lsdtt_raster_cpu prod =
            lsdtt::calculate_production_raster(fx.dem);
        const lsdtt::lsdtt_raster_cpu topo = make_uniform_raster(fx, 1.0f);
        const lsdtt::lsdtt_raster_cpu zero = make_uniform_raster(fx, 0.0f);
        out.push_back({"production_outlet", "-", "-", 0.0, 0.0, 0.0,
                       prod.get_data()[static_cast<std::size_t>(fx.outlet_node)]});
        for (double meas : {100000.0, 20000.0})
        {
            double pred = 0.0;
            const double erate = lsdtt::calculate_eff_erate_from_conc(
                meas, "Be10", "BraucherBorchers", prod, topo, zero, zero,
                flow, fx.outlet_node, pred);
            out.push_back({"erate_inversion", "Be10", "BraucherBorchers",
                           meas, 0.0, 0.0, erate});
            out.push_back({"erate_inversion_predconc", "Be10",
                           "BraucherBorchers", meas, 0.0, 0.0, pred});
        }
    }

    return out;
}

void write_sentinel_file(const std::string& path,
                         const std::vector<sentinel_entry>& entries)
{
    std::ofstream out(path);
    out << std::setprecision(17);
    out << "# COSMO NUMERICS SENTINEL - see docs/testing_and_benchmarking.md\n"
        << "# These values encode years of manual verification of the CAIRN\n"
        << "# and step-change analytical solutions.  DO NOT regenerate this\n"
        << "# file to silence a test failure.  Any change requires review\n"
        << "# and sign-off by Simon Mudd BEFORE merging.\n"
        << "quantity,nuclide,scheme,p1,p2,p3,value\n";
    for (const auto& e : entries)
    {
        out << e.quantity << "," << e.nuclide << "," << e.scheme << ","
            << e.p1 << "," << e.p2 << "," << e.p3 << "," << e.value << "\n";
    }
}

bool read_sentinel_file(const std::string& path,
                        std::map<std::string, double>& values)
{
    std::ifstream in(path);
    if (!in.good()) return false;
    std::string line;
    while (std::getline(in, line))
    {
        if (line.empty() || line[0] == '#' || line.rfind("quantity,", 0) == 0)
            continue;
        std::istringstream ls(line);
        sentinel_entry e;
        std::string field;
        std::getline(ls, e.quantity, ',');
        std::getline(ls, e.nuclide, ',');
        std::getline(ls, e.scheme, ',');
        std::getline(ls, field, ','); e.p1 = std::stod(field);
        std::getline(ls, field, ','); e.p2 = std::stod(field);
        std::getline(ls, field, ','); e.p3 = std::stod(field);
        std::getline(ls, field, ','); e.value = std::stod(field);
        values[sentinel_key(e)] = e.value;
    }
    return true;
}

void phase1_sentinel(const single_basin_fixture& fx,
                     const lsdtt::lsdtt_flowinfo& flow,
                     const std::string& sentinel_path, bool regenerate)
{
    const std::vector<sentinel_entry> live = compute_sentinel_values(fx, flow);

    if (regenerate)
    {
        write_sentinel_file(sentinel_path, live);
        std::cout << "NOTE  sentinel REGENERATED at " << sentinel_path << " ("
                  << live.size() << " values).  This file must be reviewed\n"
                  << "      and signed off by Simon Mudd before committing.\n";
        check("sentinel: regenerated (review required before commit)", true);
        return;
    }

    std::map<std::string, double> pinned;
    if (!read_sentinel_file(sentinel_path, pinned))
    {
        check("sentinel: reference file present", false,
              "missing " + sentinel_path
                  + " - bootstrap with --regenerate-sentinel-i-have-simons-approval");
        return;
    }

    int n_bad = 0, n_missing = 0;
    std::ostringstream detail;
    for (const auto& e : live)
    {
        const auto it = pinned.find(sentinel_key(e));
        if (it == pinned.end())
        {
            ++n_missing;
            continue;
        }
        if (!approx_rel(e.value, it->second, k_sentinel_rel_tol))
        {
            ++n_bad;
            if (n_bad <= 5)
            {
                detail << "\n        " << sentinel_key(e) << ": "
                       << std::setprecision(17) << e.value << " vs pinned "
                       << it->second;
            }
        }
    }

    const bool ok = (n_bad == 0 && n_missing == 0
                     && live.size() == pinned.size());
    if (!ok)
    {
        std::cout
            << "**********************************************************\n"
            << "* COSMO NUMERICS CHANGED.  These values encode years of  *\n"
            << "* manual verification of the CAIRN and step-change       *\n"
            << "* analytical solutions.  Do NOT regenerate the sentinel  *\n"
            << "* file.  Any change requires review and sign-off by      *\n"
            << "* Simon Mudd before merging.                             *\n"
            << "**********************************************************\n"
            << "  mismatched=" << n_bad << " missing=" << n_missing
            << " live=" << live.size() << " pinned=" << pinned.size()
            << detail.str() << "\n";
    }
    check("sentinel: all " + std::to_string(live.size())
              + " pinned numerics unchanged (rel tol 1e-12)",
          ok);
}

// =============================================================================
// Phase 4: bundled-data auto-discovery + production-vs-depth emitter (plan Sec9)
// =============================================================================

// The cosmo CLIs locate the shipped crn_parameters.nc relative to the running
// executable (src/core/lsdtt_exe_path.hpp -> default_crn_nc_path()).  The test
// binary lives in build/tests/, so the build-tree layout <exe>/../share/lsdtt3
// resolves to build/share/lsdtt3/crn_parameters.nc.  This guards the discovery
// path the CLIs rely on when the user gives no crn_nc_file.
lsdtt::CrnNcData phase4_auto_discovery()
{
    const std::string nc_path = lsdtt::default_crn_nc_path();
    const bool found = !nc_path.empty();
    check("Phase4: bundled crn_parameters.nc auto-discovered via exe path",
          found, "default_crn_nc_path() returned empty (share/lsdtt3 not staged?)");

    lsdtt::CrnNcData nc;
    if (found)
    {
        nc = lsdtt::lsdtt_crn_parameters::load_nc_data(nc_path);
        check("Phase4: discovered crn_parameters.nc loads (ERA40 + VDM present)",
              nc.loaded,
              "load_nc_data reported loaded=false for " + nc_path);
    }
    return nc;
}

// production_vs_depth.csv (plan Sec7 / Fig. 8 analogue): per-pathway and total
// production P_j(z) = P0 * F_j * exp(-z/Gamma_j) over depth for each scheme and
// nuclide.  Pure library values (no external oracle - the Balco full-integration
// reference of plan Sec6d D1 is out of scope here).  Asserts the physical
// invariants: total production is positive, strictly decreasing with depth, and
// equals P0*sum(F) at the surface.
void phase4_production_vs_depth(const std::string& out_dir)
{
    csv_writer csv(out_dir + "/production_vs_depth.csv",
                   "nuclide,scheme,depth_gcm2,P_spallation,P_slowmuon,"
                   "P_fastmuon,P_veryfastmuon,P_total,source");
    const char* schemes[]  = {"Braucher", "BraucherBorchers", "newCRONUS"};
    const char* nuclides[] = {"Be10", "Al26", "C14"};
    const double depths[]  = {0.0, 50.0, 100.0, 200.0, 500.0, 1000.0,
                              2000.0, 3000.0, 5000.0};
    bool ok = true;
    std::string bad;
    for (const char* scheme : schemes)
    {
        for (const char* nuc : nuclides)
        {
            const ref_params p = ref_params_from_scheme(scheme, nuc);
            double f_sum = 0.0;
            for (int j = 0; j < 4; ++j) f_sum += p.F[j];

            double prev_total = -1.0;
            for (double z : depths)
            {
                double pj[4], total = 0.0;
                for (int j = 0; j < 4; ++j)
                {
                    pj[j] = p.P0 * p.F[j] * std::exp(-z / p.Gamma[j]);
                    total += pj[j];
                }
                csv.row(nuc, scheme, z, pj[0], pj[1], pj[2], pj[3], total,
                        "lsdtt");
                if (!std::isfinite(total) || total <= 0.0) { ok = false; bad = "non-positive total"; }
                if (prev_total > 0.0 && total >= prev_total) { ok = false; bad = "not decreasing with depth"; }
                if (z == 0.0 && !approx_rel(total, p.P0 * f_sum, 1e-12))
                { ok = false; bad = "surface != P0*sum(F)"; }
                prev_total = total;
            }
        }
    }
    check("Phase4: production-vs-depth positive, monotone, surface==P0*sum(F)",
          ok, bad);
    check("csv: production_vs_depth.csv written", csv.good());
}

// =============================================================================
// Phase 5: geomagnetic / time-varying production tests G1-G4 + Sec6c(A)
// =============================================================================
//
// Exercises stone2000sp_paleomag (the Muscheler/Knudsen VDM time-integrated
// Stone/Lal correction) and the ERA40 atmosphere - the ONLY physics in the
// core that had zero coverage.  The core is LOCKED: these tests call it
// read-only, build synthetic CrnNcData in-code (no external file), and pin our
// trajectory against the committed data/crn_parameters.nc (Knudsen VDM record).
// They never assert exact agreement with CRONUS's Lm/LSDn (which use a
// different VDM record); see plan Sec6c.

// Synthetic CrnNcData with a CONSTANT VDM record (G1 constant-field limit).
// loaded=true is enough for stone2000sp_paleomag: it only reads muscheler_*.
// A constant VDM makes the time-average equal the modern reference, so the
// correction is exactly 1 for any age/lat/pressure - this validates the
// time-integration loop itself, independent of any real record.
lsdtt::CrnNcData make_constant_vdm_nc(double vdm_value, int n = 50,
                                      double dt_yr = 1000.0)
{
    lsdtt::CrnNcData nc;
    for (int i = 0; i < n; ++i)
    {
        nc.muscheler_time.push_back(static_cast<double>(i) * dt_yr);
        nc.muscheler_vdm.push_back(static_cast<float>(vdm_value));
    }
    nc.loaded = true;
    return nc;
}

// Synthetic record: modern (index 0) VDM, then a uniform "past" level for t>0.
// Lets a test dial the sign/size of the correction deterministically.
lsdtt::CrnNcData make_two_level_vdm_nc(double modern_vdm, double past_vdm,
                                       int n = 50, double dt_yr = 1000.0)
{
    lsdtt::CrnNcData nc;
    for (int i = 0; i < n; ++i)
    {
        nc.muscheler_time.push_back(static_cast<double>(i) * dt_yr);
        nc.muscheler_vdm.push_back(
            static_cast<float>(i == 0 ? modern_vdm : past_vdm));
    }
    nc.loaded = true;
    return nc;
}

// G1 - static-field / degenerate limits.
void phase5_G1_static_limits(const lsdtt::CrnNcData& nc_real)
{
    lsdtt::lsdtt_crn_parameters cp;

    // G1(a) constant-VDM record reproduces the time-independent Stone factor to
    // machine precision for ANY age/lat/pressure (correction == 1 exactly).
    {
        const lsdtt::CrnNcData nc = make_constant_vdm_nc(7.5);
        bool ok = true;
        std::string bad;
        for (double lat : {0.0, 30.0, 60.0})
            for (double P : {1013.25, 700.0})
                for (double age : {1000.0, 20000.0, 49000.0})
                {
                    const double on  = cp.stone2000sp_paleomag(lat, P, 1.0, age, nc);
                    const double off = cp.stone2000sp(lat, P, 1.0);
                    if (!approx_rel(on, off, tol::g1_static))
                    { ok = false; bad = diff_str(on, off); }
                }
        check("G1: constant-VDM record == time-independent Stone (any age)",
              ok, bad);
    }

    // G1(b) degenerate fallbacks return the time-independent result EXACTLY:
    // not loaded, age <= 0, and fewer than two record points within the age.
    {
        const lsdtt::CrnNcData empty;                 // loaded == false
        const lsdtt::CrnNcData nc = make_constant_vdm_nc(7.5, 50, 1000.0);
        const double lat = 0.0, P = 1013.25;
        const double off = cp.stone2000sp(lat, P, 1.0);
        const bool ok =
            cp.stone2000sp_paleomag(lat, P, 1.0, 20000.0, empty) == off   // not loaded
         && cp.stone2000sp_paleomag(lat, P, 1.0, 0.0,     nc)    == off   // age == 0
         && cp.stone2000sp_paleomag(lat, P, 1.0, -5.0,    nc)    == off   // age < 0
         && cp.stone2000sp_paleomag(lat, P, 1.0, 500.0,   nc)    == off;  // < 2 pts <= age
        check("G1: degenerate fallbacks return time-independent Stone (exact)",
              ok);
    }

    // G1(a') real record at age 0 is exactly the time-independent Stone factor
    // (the geomag-OFF path the CLIs use when sample_age_yr == 0).
    if (nc_real.loaded)
    {
        const double lat = 0.0, P = 1013.25;
        check("G1: real nc at age=0 == time-independent Stone (exact)",
              cp.stone2000sp_paleomag(lat, P, 1.0, 0.0, nc_real)
                  == cp.stone2000sp(lat, P, 1.0));
    }
}

// G2 - physics invariants (synthetic low-/high-field records).
void phase5_G2_invariants()
{
    lsdtt::lsdtt_crn_parameters cp;

    // G2(a) HIGH geomagnetic latitude: above the cutoff-rigidity knee the Stone
    // factor is on its polar plateau, so production is insensitive to VDM - a
    // large field swing barely changes the time-averaged scaling at 70 deg.
    {
        const double lat = 70.0, P = 1007.0, age = 30000.0;
        const double base    = cp.stone2000sp(lat, P, 1.0);
        const double sf_low  = cp.stone2000sp_paleomag(lat, P, 1.0, age,
                                   make_constant_vdm_nc(3.0));
        const double sf_high = cp.stone2000sp_paleomag(lat, P, 1.0, age,
                                   make_two_level_vdm_nc(8.0, 12.0));
        check("G2a: high-lat production insensitive to VDM (<0.5%)",
              std::abs(sf_low - base) / base < 5e-3
                  && std::abs(sf_high - base) / base < 5e-3,
              "low=" + diff_str(sf_low, base) + " high=" + diff_str(sf_high, base));
    }

    // G2(b) LOW latitude sign: at the equator higher past field -> higher cutoff
    // rigidity -> LOWER production (correction < 1); lower past field -> HIGHER
    // production (correction > 1).
    {
        const double lat = 0.0, P = 1013.25, age = 30000.0;
        const double base = cp.stone2000sp(lat, P, 1.0);
        const double sf_lowfield  = cp.stone2000sp_paleomag(lat, P, 1.0, age,
                                        make_two_level_vdm_nc(8.0, 4.0));
        const double sf_highfield = cp.stone2000sp_paleomag(lat, P, 1.0, age,
                                        make_two_level_vdm_nc(8.0, 12.0));
        check("G2b: equator low past-field -> more production (corr>1)",
              sf_lowfield > base, diff_str(sf_lowfield, base));
        check("G2b: equator high past-field -> less production (corr<1)",
              sf_highfield < base, diff_str(sf_highfield, base));
    }

    // G2(c) monotonicity: at the equator a deeper low-field excursion gives a
    // larger correction (production rises as the past field falls).
    {
        const double lat = 0.0, P = 1013.25, age = 30000.0;
        double prev = -1.0;
        bool mono = true;
        for (double past : {8.0, 6.0, 4.0, 2.0})   // modern 8, past field falling
        {
            const double sf = cp.stone2000sp_paleomag(lat, P, 1.0, age,
                                  make_two_level_vdm_nc(8.0, past));
            if (prev > 0.0 && sf < prev) mono = false;
            prev = sf;
        }
        check("G2c: equatorial correction monotone in past-field depth", mono);
    }

    // G2(d) ERA40 interpolation on a synthetic UNIFORM grid: bilinear returns
    // the constant exactly, pressure at z=0 equals the SLP, and pressure falls
    // with elevation; the difference from the standard atmosphere is bounded.
    {
        lsdtt::CrnNcData nc;
        nc.era40_lat = {90.0, -90.0};    // descending
        nc.era40_lon = {0.0, 360.0};     // ascending
        const float Pconst = 1013.0f, Tconst = 288.0f;
        nc.era40_meanP = {Pconst, Pconst, Pconst, Pconst};
        nc.era40_meanT = {Tconst, Tconst, Tconst, Tconst};
        nc.lifton_lr   = {-6.1517e-3, 0, 0, 0, 0, 0, 0};   // constant lapse rate
        nc.muscheler_time = {0.0, 1000.0};
        nc.muscheler_vdm  = {7.5f, 7.5f};
        nc.loaded = true;

        const double p0 = lsdtt::lsdtt_crn_parameters::atmospheric_pressure_era40(
            0.0, 45.0, 0.0, nc);
        check("G2d: ERA40 uniform-grid pressure at z=0 == SLP",
              approx_rel(p0, Pconst, 1e-6), diff_str(p0, Pconst));

        const double p1000 = lsdtt::lsdtt_crn_parameters::atmospheric_pressure_era40(
            0.0, 45.0, 1000.0, nc);
        const double p_std = lsdtt::lsdtt_crn_parameters::atmospheric_pressure_standard(1000.0);
        check("G2d: ERA40 pressure falls with elevation and is bounded",
              p1000 < p0 && std::abs(p1000 - p_std) < 50.0,
              "era40=" + std::to_string(p1000) + " std=" + std::to_string(p_std));
    }
}

// G4 - golden trajectories against the committed Knudsen VDM record.  These pin
// OUR time-integrated Stone factor at a few (lat, pressure, age) triplets;
// measured 2026-07-23 against data/crn_parameters.nc and frozen at tol::g4_golden.
// They are NOT a CRONUS comparison (different VDM record - plan Sec6c).
void phase5_G4_golden_trajectories(const lsdtt::CrnNcData& nc_real,
                                   const std::string& out_dir)
{
    if (!nc_real.loaded)
    {
        check("G4: golden trajectories (SKIPPED - real nc not loaded)", true);
        return;
    }
    lsdtt::lsdtt_crn_parameters cp;
    struct golden { double lat, P, age, expected; };
    const golden g[] = {
        {0.0,  1013.25, 20000.0, 0.60474384759202182},
        {0.0,  700.0,   20000.0, 4.738527144189737},
        {30.0, 1013.25, 20000.0, 0.82691194517522892},
        {60.0, 1013.25, 20000.0, 1.0026856559878936},
        {0.0,  1013.25, 50000.0, 0.70258020553273326},
        {70.0, 1007.0,  20000.0, 1.0612665967394139},
    };
    csv_writer csv(out_dir + "/geomag_scaling.csv",
                   "lat,pressure_hPa,t_yr,sf_paleomag,sf_modern,ratio,source");
    bool ok = true;
    std::string bad;
    for (const golden& x : g)
    {
        const double got = cp.stone2000sp_paleomag(x.lat, x.P, 1.0, x.age, nc_real);
        const double off = cp.stone2000sp(x.lat, x.P, 1.0);
        csv.row(x.lat, x.P, x.age, got, off, got / off, "lsdtt_knudsen");
        if (!approx_rel(got, x.expected, tol::g4_golden))
        { ok = false; bad = "lat=" + std::to_string(x.lat) + " age="
                            + std::to_string(x.age) + ": " + diff_str(got, x.expected); }
    }
    check("G4: paleomag goldens match frozen Knudsen trajectory (1e-9)", ok, bad);
    check("csv: geomag_scaling.csv written", csv.good());
}

// Sec6c(A) - equatorial geomag ON/OFF at the LSDTT-1 fixture (lat 0, 250 m),
// where the VDM/cutoff-rigidity effect is strongest.  Reports the ON-OFF
// production-scaling delta across integration ages and asserts its documented
// SIGN behaviour and a frozen magnitude band.  The Knudsen record has HIGHER-
// than-modern VDM in the last few kyr (vdm[1..] > vdm[0]) and lower-than-modern
// field deeper in time, so the delta reverses sign from short to long ages.
//
// The frozen CRONUS reference's Lm-vs-St spread at LSDTT-1 is read for CONTEXT
// only: it is the size of a time-dependent correction, but its sign need not
// match ours - CRONUS uses a different (Pavon-family) VDM record, so absolute
// agreement is not expected and is NOT asserted (plan Sec6c).
void phase5_G6c_equator_on_off(const lsdtt::CrnNcData& nc_real,
                               const std::string& expected_dir,
                               const std::string& out_dir)
{
    if (!nc_real.loaded)
    {
        check("6cA: equatorial ON/OFF (SKIPPED - real nc not loaded)", true);
        return;
    }
    lsdtt::lsdtt_crn_parameters cp;
    const double lat = 0.0;
    const double P =
        lsdtt::lsdtt_crn_parameters::atmospheric_pressure_standard(250.0);
    const double off = cp.stone2000sp(lat, P, 1.0);

    csv_writer csv(out_dir + "/geomag_on_off_LSDTT1.csv",
                   "lat,elev_m,t_yr,sf_off,sf_on,delta_pct,source");
    double worst_abs = 0.0;
    for (double age : {5000.0, 10000.0, 20000.0, 49500.0})
    {
        const double on  = cp.stone2000sp_paleomag(lat, P, 1.0, age, nc_real);
        const double dpct = 100.0 * (on - off) / off;
        csv.row(lat, 250.0, age, off, on, dpct, "lsdtt_knudsen");
        worst_abs = std::max(worst_abs, std::abs(dpct));
    }

    const double d_short =
        100.0 * (cp.stone2000sp_paleomag(lat, P, 1.0, 5000.0, nc_real) - off) / off;
    const double d_long =
        100.0 * (cp.stone2000sp_paleomag(lat, P, 1.0, 49500.0, nc_real) - off) / off;

    check("6cA: short-age ON-OFF delta negative (recent VDM above modern)",
          d_short < 0.0, "delta = " + std::to_string(d_short) + " %");
    check("6cA: record-spanning ON-OFF delta positive (deep field below modern)",
          d_long > 0.0, "delta = " + std::to_string(d_long) + " %");
    check("6cA: equatorial ON-OFF |delta| within frozen band (<=25%)",
          worst_abs <= tol::g6c_band,
          "worst = " + std::to_string(worst_abs) + " %");

    // CONTEXT (reported, not asserted): CRONUS St/Lm spread at LSDTT-1.
    const std::string ref_path = expected_dir + "/cronus_v3_reference.csv";
    std::ifstream in(ref_path);
    double st = 0.0, lm = 0.0;
    if (in.good())
    {
        std::string line;
        while (std::getline(in, line))
        {
            if (line.rfind("LSDTT-1,Be-10,St,", 0) == 0)
                st = std::stod(line.substr(line.find(",St,") + 4));
            else if (line.rfind("LSDTT-1,Be-10,Lm,", 0) == 0)
                lm = std::stod(line.substr(line.find(",Lm,") + 4));
        }
    }
    if (st > 0.0 && lm > 0.0)
    {
        std::cout << "      6cA CONTEXT: CRONUS LSDTT-1 Be-10 St=" << st
                  << " Lm=" << lm << " (Lm/St=" << lm / st << ", +"
                  << std::setprecision(3) << 100.0 * (lm - st) / st << " %); "
                  << "our equatorial ON-OFF production delta spans "
                  << std::setprecision(4) << d_short << " to " << d_long
                  << " % across 5-49.5 kyr (Knudsen VDM; opposite net sign to "
                  << "CRONUS is the record difference, plan Sec6c).\n"
                  << std::setprecision(17);
    }
    check("csv: geomag_on_off_LSDTT1.csv written", csv.good());
}

// =============================================================================
// Phase 0 smoke tests
// =============================================================================

// Reference-helper identities: quick versions of plan E1/E2 that prove the
// reference implementations are self-consistent.  (The full battery with
// sweeps, all schemes, and tight documented tolerances is Phase 1.)
void smoke_reference_identities()
{
    const ref_params p = ref_params_from_scheme("BraucherBorchers", "Be10");
    const double eps_o = 0.01, eps_n = 0.05;   // g cm^-2 yr^-1
    const double d = 15.0, t = 25000.0;

    // Eq. (7) at t = 0 equals Eq. (12) at eps_old.
    check("ref: step-change t=0 == steady(eps_old)",
          approx_rel(ref_conc_step_change(p, d, 0.0, eps_o, eps_n),
                     ref_conc_surface(p, eps_o, d), 1e-14),
          diff_str(ref_conc_step_change(p, d, 0.0, eps_o, eps_n),
                   ref_conc_surface(p, eps_o, d)));

    // Eq. (7) with eps_old == eps_new equals Eq. (12) at any t.
    check("ref: step-change eps_o==eps_n == steady",
          approx_rel(ref_conc_step_change(p, d, t, eps_n, eps_n),
                     ref_conc_surface(p, eps_n, d), 1e-14),
          diff_str(ref_conc_step_change(p, d, t, eps_n, eps_n),
                   ref_conc_surface(p, eps_n, d)));

    // Eq. (7) at large t approaches Eq. (12) at eps_new.
    const double t_inf = 5.0e7;
    check("ref: step-change t->inf -> steady(eps_new)",
          approx_rel(ref_conc_step_change(p, d, t_inf, eps_o, eps_n),
                     ref_conc_surface(p, eps_n, d), 1e-9),
          diff_str(ref_conc_step_change(p, d, t_inf, eps_o, eps_n),
                   ref_conc_surface(p, eps_n, d)));

    // Eq. (14) -> Eq. (12) as d_t -> 0 (d_s = 0 surface case).
    check("ref: depth-integrated d_t->0 -> surface",
          approx_rel(ref_conc_depth_integrated(p, eps_o, 0.0, 1e-6),
                     ref_conc_surface(p, eps_o, 0.0), 1e-6),
          diff_str(ref_conc_depth_integrated(p, eps_o, 0.0, 1e-6),
                   ref_conc_surface(p, eps_o, 0.0)));

    // Eq. (8) at t = 0 equals the depth-integrated steady state at eps_old
    // (CAIRN Eq. 14 with d_s = d_t_top, d_t = thickness).
    const double d_top = 5.0, d_bot = 55.0;
    check("ref: Eq8 t=0 == depth-integrated steady(eps_old)",
          approx_rel(
              ref_conc_step_change_depth_integrated(p, d_top, d_bot, 0.0,
                                                    eps_o, eps_n),
              ref_conc_depth_integrated(p, eps_o, d_top, d_bot - d_top),
              1e-14),
          diff_str(
              ref_conc_step_change_depth_integrated(p, d_top, d_bot, 0.0,
                                                    eps_o, eps_n),
              ref_conc_depth_integrated(p, eps_o, d_top, d_bot - d_top)));

    // Eq. (8) -> Eq. (7) as the slab shrinks.
    const double d_mid = 20.0, half = 1e-5;
    check("ref: Eq8 thin slab -> Eq7",
          approx_rel(
              ref_conc_step_change_depth_integrated(p, d_mid - half,
                                                    d_mid + half, t,
                                                    eps_o, eps_n),
              ref_conc_step_change(p, d_mid, t, eps_o, eps_n), 1e-8),
          diff_str(
              ref_conc_step_change_depth_integrated(p, d_mid - half,
                                                    d_mid + half, t,
                                                    eps_o, eps_n),
              ref_conc_step_change(p, d_mid, t, eps_o, eps_n)));
}

// Fixture sanity: single outlet that captures every valid pixel.
void smoke_fixture(const single_basin_fixture& fx,
                   const lsdtt::lsdtt_flowinfo& flow)
{
    const auto& meta = fx.meta;
    const std::int64_t n_valid =
        static_cast<std::int64_t>(meta.n_rows) * meta.n_cols;

    const auto npix_res = flow.get_n_contributing_pixels_raster();
    const double outlet_npix =
        npix_res.data[static_cast<std::size_t>(fx.outlet_node)];
    check("fixture: all pixels drain to the single outlet",
          static_cast<std::int64_t>(outlet_npix + 0.5) == n_valid,
          "outlet contributing pixels = "
              + std::to_string(static_cast<std::int64_t>(outlet_npix))
              + ", expected " + std::to_string(n_valid));
}

// Library smoke: the full forward pipeline runs on the fixture and behaves
// physically (finite, positive, monotone-decreasing in erosion rate); the
// erosion-weighted accumulation of a uniform concentration field is that
// constant.  Also emits the first Sec7 CSV skeletons.
void smoke_library_pipeline(const single_basin_fixture& fx,
                            const lsdtt::lsdtt_flowinfo& flow,
                            const std::string& out_dir)
{
    const lsdtt::lsdtt_raster_cpu prod =
        lsdtt::calculate_production_raster(fx.dem);
    const lsdtt::lsdtt_raster_cpu shield_topo = make_uniform_raster(fx, 1.0f);
    const lsdtt::lsdtt_raster_cpu shield_zero = make_uniform_raster(fx, 0.0f);

    const double prod_outlet =
        prod.get_data()[static_cast<std::size_t>(fx.outlet_node)];
    check("library: production raster finite and positive at outlet",
          std::isfinite(prod_outlet) && prod_outlet > 0.0,
          "prod = " + std::to_string(prod_outlet));

    // Whole-raster steady state at two erosion rates: C must fall as eps rises.
    csv_writer conc_csv(out_dir + "/conc_vs_erate.csv",
                        "nuclide,scheme,fixture,erate_gcm2yr,conc_atoms_g,source");
    double c_prev = -1.0;
    bool monotone = true, all_finite = true;
    const std::vector<double> eps_sweep = {0.001, 0.01, 0.1};
    for (double eps : eps_sweep)
    {
        const lsdtt::lsdtt_raster_cpu ero =
            make_uniform_raster(fx, static_cast<float>(eps));
        const lsdtt::lsdtt_raster_cpu conc =
            lsdtt::calculate_crn_concentration_raster_whole(
                "Be10", "BraucherBorchers", ero, prod,
                shield_topo, shield_zero, shield_zero);
        const double c_outlet =
            conc.get_data()[static_cast<std::size_t>(fx.outlet_node)];
        if (!std::isfinite(c_outlet) || c_outlet <= 0.0) all_finite = false;
        if (c_prev > 0.0 && c_outlet >= c_prev) monotone = false;
        conc_csv.row("Be10", "BraucherBorchers", "single_basin",
                     eps, c_outlet, "lsdtt");
        c_prev = c_outlet;
    }
    check("library: whole-raster concentration finite and positive", all_finite);
    check("library: concentration monotone decreasing in erosion rate", monotone);
    check("csv: conc_vs_erate.csv written", conc_csv.good());

    // Accumulation of a uniform concentration field is that constant,
    // independent of the erosion weighting (weighted mean of a constant).
    {
        const float c0 = 12345.0f;
        const lsdtt::lsdtt_raster_cpu conc_u = make_uniform_raster(fx, c0);
        const lsdtt::lsdtt_raster_cpu ero =
            make_two_rate_raster(fx, fx.meta.n_rows / 2, 0.01f, 0.05f);
        const double acc = lsdtt::calculate_accumulated_crn_concentration(
            conc_u, ero, flow, fx.outlet_node);
        check("library: accumulation of uniform concentration == constant",
              approx_rel(acc, static_cast<double>(c0), 1e-6),
              diff_str(acc, c0));
    }

    // Step-change relaxation CSV skeleton (knickpoint Fig. 1 analogue):
    // reference-equation curves for 10Be and 14C at the outlet's production.
    {
        csv_writer relax_csv(out_dir + "/step_change_relaxation.csv",
                             "nuclide,scheme,t_yr,conc_atoms_g,conc_over_ss_new,source");
        const double eps_o = 0.005, eps_n = 0.05;
        bool good = true;
        for (const std::string& nuc : {std::string("Be10"), std::string("C14")})
        {
            ref_params p = ref_params_from_scheme("BraucherBorchers", nuc);
            for (int j = 0; j < 4; ++j) p.S[j] = 1.0;
            const double c_ss_new = ref_conc_surface(p, eps_n, 0.0);
            for (double t = 0.0; t <= 2.0e5; t += 5000.0)
            {
                const double c = ref_conc_step_change(p, 0.0, t, eps_o, eps_n);
                relax_csv.row(nuc, "BraucherBorchers", t, c, c / c_ss_new,
                              "reference");
                if (!std::isfinite(c)) good = false;
            }
        }
        check("csv: step_change_relaxation.csv written",
              good && relax_csv.good());
    }
}

} // namespace

// =============================================================================
int main(int argc, char** argv)
{
    std::string out_dir = ".";
    std::string sentinel_path = "numerics_sentinel.csv";
    bool regenerate_sentinel = false;
    int positional = 0;
    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if (arg == "--regenerate-sentinel-i-have-simons-approval")
            regenerate_sentinel = true;
        else if (positional == 0) { out_dir = arg; ++positional; }
        else if (positional == 1) { sentinel_path = arg; ++positional; }
    }

    std::cout << "=========================================================\n"
              << "|| lsdtt3 cosmogenic physics test driver               ||\n"
              << "|| Method: docs/testing_and_benchmarking.md            ||\n"
              << "|| Fixtures, parameter pins, numerics sentinel,        ||\n"
              << "|| and equation fidelity E1-E9                         ||\n"
              << "=========================================================\n";

    try
    {
        smoke_reference_identities();

        single_basin_fixture fx = make_single_basin(40, 31);
        lsdtt::lsdtt_flowinfo flow(fx.dem.get_data_array(),
                                   fx.dem.get_mask_array(),
                                   fx.meta);

        smoke_fixture(fx, flow);
        smoke_library_pipeline(fx, flow, out_dir);

        phase1_parameter_pins();
        phase1_standardisation_pins();
        phase1_E1_depth_integration_limits();
        phase1_E2_step_change_library(fx);
        phase1_E3_exponential_identity();
        phase1_E4_lambda_v_solve();
        phase1_E5_linearity_and_limits();

        phase2_E7_analytic_vs_integrator(fx, out_dir);
        phase2_E6_round_trip(fx, flow);
        phase2_E8_averaging_identities(fx, flow);

        // expected/ dir (frozen references) sits next to the sentinel file.
        std::string expected_dir = ".";
        const auto slash = sentinel_path.find_last_of("/\\");
        if (slash != std::string::npos)
            expected_dir = sentinel_path.substr(0, slash);
        phase3_TC3_stone_pins();
        phase3_TC1_site_curve(out_dir);
        phase3_write_cronus_fixture(out_dir);
        phase3_TC2_cronus_comparison(expected_dir, out_dir);

        // Phase 4: bundled-data auto-discovery + production-vs-depth emitter.
        const lsdtt::CrnNcData crn_nc = phase4_auto_discovery();
        phase4_production_vs_depth(out_dir);

        // Phase 5: geomagnetic / time-varying production (G1-G4 + Sec6c(A)).
        phase5_G1_static_limits(crn_nc);
        phase5_G2_invariants();
        phase5_G4_golden_trajectories(crn_nc, out_dir);
        phase5_G6c_equator_on_off(crn_nc, expected_dir, out_dir);

        phase1_sentinel(fx, flow, sentinel_path, regenerate_sentinel);
    }
    catch (const std::exception& e)
    {
        std::cout << "FAIL  unhandled exception: " << e.what() << "\n";
        ++g_n_checks;
        ++g_n_failed;
    }

    std::cout << "\n" << (g_n_checks - g_n_failed) << "/" << g_n_checks
              << " checks passed\n";
    std::cout << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";
    return g_n_failed == 0 ? 0 : 1;
}
