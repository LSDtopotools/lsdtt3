// test_spectral_significant_2d.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// 2D significant frequencies and background-seed reproducibility of
// full_spectral_analysis (lsdtt_raster_spectral_cpu).
//
// Modes, one per CTest test:
//
//   --threshold   spectral_normalised_psd_threshold against the closed form
//                 -ln(1 - p), and its rejection of p outside (0, 1).
//   --bearings    Known-answer sinusoids. A plane wave whose frequency falls
//                 exactly on a spectrum cell puts its peak in that cell (the
//                 Hann window spreads power to the neighbours but keeps the
//                 maximum in place), so the strongest significant pair must
//                 sit at the wave's bearing and that bearing + 180, at its
//                 frequency. Cases: wave vector north (0/180), east (90/270),
//                 north-east (45/225), and north-east on a 64 x 128 grid where
//                 dfy = 2 dfx. On that grid the legacy index-offset bearing
//                 would give atan2(8, 4) = 63.43 deg, not 45.
//   --seed        The same background_seed gives bit-identical background
//                 spectra at 1 and 4 threads; a different seed does not.
//
// Row 0 of a raster is its northern edge, so northing increases as the row
// index decreases. A wave exp(2 pi i (kx col - ky row) / N) therefore has a
// wave vector pointing kx east and ky north.
//
// Prints one PASS/FAIL line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_spectral_cpu.hpp"

namespace
{

int g_n_checks = 0;
int g_n_failed = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;
    std::cout << (ok ? "PASS: " : "FAIL: ") << name;
    if (!detail.empty()) std::cout << "  [" << detail << "]";
    std::cout << "\n";
}

std::string num(double v)
{
    std::ostringstream os;
    os << std::setprecision(10) << v;
    return os.str();
}

int finish()
{
    std::cout << "\n" << g_n_checks - g_n_failed << " / " << g_n_checks
              << " checks passed\nRESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL")
              << "\n";
    return g_n_failed == 0 ? 0 : 1;
}

constexpr double k_pi = 3.14159265358979323846;

// A plane wave of amplitude 1 with kx cycles across the columns and ky cycles
// across the rows (north positive), plus Gaussian noise of standard deviation
// noise_sd from a fixed seed. 1 m cells, synthetic UTM georeference.
lsdtt::lsdtt_raster_cpu make_wave_dem(int n_rows, int n_cols, int kx, int ky,
                                      double noise_sd)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows        = n_rows;
    meta.n_cols        = n_cols;
    meta.cell_size     = 1.0;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();

    std::mt19937 rng(12345);
    std::normal_distribution<double> noise(0.0, noise_sd);

    std::vector<float>   data(static_cast<std::size_t>(n_rows) * n_cols);
    std::vector<uint8_t> mask(data.size(), 1);
    for (int r = 0; r < n_rows; ++r)
    {
        for (int c = 0; c < n_cols; ++c)
        {
            const double phase = 2.0 * k_pi *
                (static_cast<double>(kx) * c / n_cols -
                 static_cast<double>(ky) * r / n_rows);
            data[static_cast<std::size_t>(r) * n_cols + c] =
                static_cast<float>(100.0 + std::sin(phase) + noise(rng));
        }
    }
    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

double angular_difference_deg(double a, double b)
{
    double d = std::fmod(std::fabs(a - b), 360.0);
    return d > 180.0 ? 360.0 - d : d;
}

// -----------------------------------------------------------------------------
// --threshold
// -----------------------------------------------------------------------------
void run_threshold_checks()
{
    const std::vector<double> levels = {0.5, 0.95, 0.99, 0.9999};
    for (double p : levels)
    {
        const double got  = lsdtt::spectral_normalised_psd_threshold(p);
        const double want = -std::log(1.0 - p);
        check("threshold(" + num(p) + ") = -ln(1-p)",
              std::fabs(got - want) < 1e-12 * want,
              "got " + num(got) + " want " + num(want));
    }
    // The legacy constant was 0.5 * 18.42185 = 9.210925.
    check("threshold(0.9999) within 1e-3 of legacy 9.210925",
          std::fabs(lsdtt::spectral_normalised_psd_threshold(0.9999) - 9.210925) < 1e-3);

    const std::vector<double> bad = {0.0, 1.0, -0.1, 1.5,
                                     std::numeric_limits<double>::quiet_NaN()};
    for (double p : bad)
    {
        bool threw = false;
        try { lsdtt::spectral_normalised_psd_threshold(p); }
        catch (const std::invalid_argument&) { threw = true; }
        check("threshold(" + num(p) + ") throws invalid_argument", threw);
    }
}

// -----------------------------------------------------------------------------
// --bearings
// -----------------------------------------------------------------------------
void run_bearing_case(const std::string& label, int n_rows, int n_cols,
                      int kx, int ky, double expected_bearing_deg)
{
    const lsdtt::lsdtt_raster_cpu dem = make_wave_dem(n_rows, n_cols, kx, ky, 0.2);

    lsdtt::lsdtt_spectral_analysis_options opts;
    opts.log_bin_width           = 0.05;
    opts.n_background_iterations = 5;
    opts.background_seed         = 7;
    opts.n_threads               = 1;
    const lsdtt::lsdtt_spectral_analysis_result res =
        lsdtt::full_spectral_analysis(dem, opts);

    const lsdtt::lsdtt_spectral_significant_2d sig =
        lsdtt::find_significant_frequencies_2d(
            res, lsdtt::spectral_normalised_psd_threshold(0.9999));

    check(label + ": some cells significant", !sig.frequency.empty(),
          num(static_cast<double>(sig.frequency.size())) + " cells");
    if (sig.frequency.empty()) return;

    // The two strongest cells are the symmetric pair of the wave's peak.
    std::size_t i1 = 0;
    for (std::size_t i = 1; i < sig.normalised_psd.size(); ++i)
        if (sig.normalised_psd[i] > sig.normalised_psd[i1]) i1 = i;
    std::size_t i2 = (i1 == 0) ? 1 : 0;
    for (std::size_t i = 0; i < sig.normalised_psd.size(); ++i)
        if (i != i1 && sig.normalised_psd[i] > sig.normalised_psd[i2]) i2 = i;

    const double f_want = std::hypot(static_cast<double>(kx) / n_cols,
                                     static_cast<double>(ky) / n_rows);
    const double b_a = std::min(sig.bearing_deg[i1], sig.bearing_deg[i2]);
    const double b_b = std::max(sig.bearing_deg[i1], sig.bearing_deg[i2]);
    const double e_a = expected_bearing_deg;
    const double e_b = expected_bearing_deg + 180.0;

    check(label + ": peak pair at bearing " + num(e_a) + " and " + num(e_b),
          angular_difference_deg(b_a, e_a) < 1e-6 &&
          angular_difference_deg(b_b, e_b) < 1e-6,
          "got " + num(b_a) + ", " + num(b_b));
    check(label + ": peak pair at frequency " + num(f_want) + " cycles/m",
          std::fabs(sig.frequency[i1] - f_want) < 1e-9 &&
          std::fabs(sig.frequency[i2] - f_want) < 1e-9,
          "got " + num(sig.frequency[i1]) + ", " + num(sig.frequency[i2]));
    check(label + ": wavelength = 1 / frequency",
          std::fabs(sig.wavelength[i1] * sig.frequency[i1] - 1.0) < 1e-12);
    check(label + ": bearing_rad and bearing_deg agree, in [0, 2 pi)",
          std::fabs(sig.bearing_rad[i1] * 180.0 / k_pi - sig.bearing_deg[i1]) < 1e-9 &&
          sig.bearing_rad[i1] >= 0.0 && sig.bearing_rad[i1] < 2.0 * k_pi);

    // Point symmetry: every significant cell has a partner at bearing + 180.
    // Cells on the most negative row/column of the shifted spectrum (north
    // component +Nyquist_y, or east component -Nyquist_x) are their own
    // wrap-around partners and are skipped.
    constexpr double nyquist = 0.5;  // cycles/m, for the 1 m cells used here
    std::size_t n_unpaired = 0;
    for (std::size_t i = 0; i < sig.bearing_deg.size(); ++i)
    {
        const double f_east  = sig.frequency[i] * std::sin(sig.bearing_rad[i]);
        const double f_north = sig.frequency[i] * std::cos(sig.bearing_rad[i]);
        if (std::fabs(f_north - nyquist) < 1e-9 || std::fabs(f_east + nyquist) < 1e-9)
            continue;
        bool found = false;
        for (std::size_t j = 0; j < sig.bearing_deg.size() && !found; ++j)
            found = angular_difference_deg(sig.bearing_deg[j],
                                           sig.bearing_deg[i] + 180.0) < 1e-6 &&
                    std::fabs(sig.frequency[j] - sig.frequency[i]) < 1e-9;
        if (!found) ++n_unpaired;
    }
    check(label + ": every significant interior cell has its +180 partner",
          n_unpaired == 0, num(static_cast<double>(n_unpaired)) + " unpaired");
}

void run_bearing_checks()
{
    run_bearing_case("north 128x128",     128, 128, 0, 8, 0.0);
    run_bearing_case("east 128x128",      128, 128, 8, 0, 90.0);
    run_bearing_case("north-east 128x128", 128, 128, 8, 8, 45.0);
    // 64 rows x 128 cols: dfx = 1/128, dfy = 1/64. kx = 8, ky = 4 gives equal
    // physical components (1/16 cycles/m each), so the bearing is 45.
    run_bearing_case("north-east 64x128", 64, 128, 8, 4, 45.0);
}

// -----------------------------------------------------------------------------
// --seed
// -----------------------------------------------------------------------------
lsdtt::lsdtt_spectral_analysis_result run_seeded(const lsdtt::lsdtt_raster_cpu& dem,
                                                 int seed, int n_threads)
{
    lsdtt::lsdtt_spectral_analysis_options opts;
    opts.log_bin_width           = 0.05;
    opts.n_background_iterations = 5;
    opts.background_seed         = seed;
    opts.n_threads               = n_threads;
    return lsdtt::full_spectral_analysis(dem, opts);
}

bool same_2d(const lsdtt::lsdtt_array2d_f& a, const lsdtt::lsdtt_array2d_f& b)
{
    if (a.n_rows() != b.n_rows() || a.n_cols() != b.n_cols()) return false;
    for (int r = 0; r < a.n_rows(); ++r)
        for (int c = 0; c < a.n_cols(); ++c)
            if (a(r, c) != b(r, c)) return false;
    return true;
}

void run_seed_checks()
{
    const lsdtt::lsdtt_raster_cpu dem = make_wave_dem(64, 64, 4, 2, 0.5);

    const auto a1 = run_seeded(dem, 42, 1);
    const auto a4 = run_seeded(dem, 42, 4);
    const auto b1 = run_seeded(dem, 43, 1);
    const auto d1 = run_seeded(dem, -1, 1);
    const auto d4 = run_seeded(dem, -1, 4);
    const auto e1 = run_seeded(dem, 1000, 1);

    check("seed 42: background_psd_2d identical at 1 and 4 threads",
          same_2d(a1.background_psd_2d, a4.background_psd_2d));
    check("seed 42: binned normalised_psd identical at 1 and 4 threads",
          a1.binned.normalised_psd == a4.binned.normalised_psd);
    check("seed 42 vs 43: background_psd_2d differs",
          !same_2d(a1.background_psd_2d, b1.background_psd_2d));
    check("default seed (-1): identical at 1 and 4 threads",
          same_2d(d1.background_psd_2d, d4.background_psd_2d));
    check("default seed (-1) is the same as seed 1000",
          same_2d(d1.background_psd_2d, e1.background_psd_2d));
    check("observed psd_2d does not depend on the seed",
          same_2d(a1.psd_2d, b1.psd_2d));

    const double t = lsdtt::spectral_normalised_psd_threshold(0.95);
    const auto s1 = lsdtt::find_significant_frequencies_2d(a1, t);
    const auto s4 = lsdtt::find_significant_frequencies_2d(a4, t);
    check("seed 42: significant 2D cells identical at 1 and 4 threads",
          s1.bearing_rad == s4.bearing_rad && s1.frequency == s4.frequency &&
          s1.normalised_psd == s4.normalised_psd);
}

} // namespace

int main(int argc, char* argv[])
{
    try
    {
        const std::string mode = argc == 2 ? argv[1] : "";
        if      (mode == "--threshold") run_threshold_checks();
        else if (mode == "--bearings")  run_bearing_checks();
        else if (mode == "--seed")      run_seed_checks();
        else
        {
            std::cerr << "usage: lsdtt-test-spectral-significant-2d "
                         "--threshold | --bearings | --seed\n";
            return 2;
        }
        return finish();
    }
    catch (const std::exception& e)
    {
        std::cerr << "EXCEPTION: " << e.what() << "\n";
        return 1;
    }
}
