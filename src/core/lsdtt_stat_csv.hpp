// lsdtt_stat_csv.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// CSV writers for statistical outputs (QQ, CDF, distribution fits).
//
// These functions are shared across multiple CLIs and avoid repeating boilerplate
// in each CLI file.
//
// All functions:
//   - take a lsdtt_data_loader& for writing (the sole GDAL/file-I/O path)
//   - take a file path as a std::string
//   - take the statistical result struct to serialise
//   - return void; throw on I/O error propagated from the loader
// -----------------------------------------------------------------------------

#pragma once

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include "lsdtt_data_loader.hpp"
#include "../utils/lsdtt_distributions.hpp"
#include "../utils/lsdtt_histogram.hpp"

namespace lsdtt {

// -----------------------------------------------------------------------------
// csv_number
//
// Formats a double for a CSV field. This is the single spelling of numeric
// output for every CSV lsdtt3 writes by hand; it replaces four private copies
// that had grown in cli/ (named number, format_double, fmt and precise_string)
// which between them used three precisions and three different non-finite
// spellings.
//
// Format: std::ostringstream in its default (defaultfloat) float format with
// precision(precision). That is exactly what all four copies did, so for
// finite values the output is byte-identical to what they produced at the same
// precision.
//
// NON-FINITE CONTRACT -- relied on by downstream readers:
//   NaN        -> "nan"
//   +infinity  -> "inf"
//   -infinity  -> "-inf"
// All lowercase. pandas.read_csv, numpy.loadtxt and R's read.csv all parse
// those three spellings natively, so a CSV written through this function reads
// back with its non-finite values intact and needs no na_values= argument.
// This also removes the platform-dependent "nan" / "-nan" that the two copies
// with no non-finite handling emitted by accident, since the sign of a NaN is
// not meaningful.
//
// PRECISION IS A PER-CALL-SITE PROPERTY, NOT A STYLE CHOICE. The default of 17
// is the round-trip precision for IEEE-754 double and is what new code should
// use. Existing call sites that were written at another precision must keep
// passing it explicitly (the swath-profile writers pass 9, the basin-extraction
// writers pass 10) so that their numeric output stays byte-identical to what
// they emitted before the unification. Do NOT "tidy" those call sites onto the
// default: doing so silently changes the contents of published CSVs.
// -----------------------------------------------------------------------------
inline std::string csv_number(double v, int precision = 17)
{
    if (std::isnan(v))     return "nan";
    if (!std::isfinite(v)) return (v > 0.0) ? "inf" : "-inf";

    std::ostringstream ss;
    ss.precision(precision);
    ss << v;
    return ss.str();
}

// -----------------------------------------------------------------------------
// write_qq_csv
//
// Writes a QQ-plot data table with three columns:
//   standard_normal_variate - z = Phi^-^1(p_k)
//   data_quantile           - empirical quantile of the sample at p_k
//   normal_reference        - IQR-fitted normal reference line at p_k
// -----------------------------------------------------------------------------
void write_qq_csv(
    lsdtt_data_loader&    loader,
    const std::string&    out_path,
    const lsdtt_qq_data&  qq
);

// -----------------------------------------------------------------------------
// write_smoothed_cdf_csv
//
// Writes the empirical CDF alongside the best-fitting parametric CDF.
// Columns:
//   bin_centre               - x value at the bin centre
//   cdf                      - empirical cumulative probability
//   fitted_cdf_<family>      - parametric CDF of the best valid fit at x
//                              (column header includes the family name)
//   cdf_minus_fitted         - empirical CDF minus fitted CDF (divergence)
//
// The best fit is the valid entry in `fits` with the lowest KS statistic.
// If no valid fit exists the third and fourth columns are written as NaN.
// -----------------------------------------------------------------------------
void write_smoothed_cdf_csv(
    lsdtt_data_loader&                         loader,
    const std::string&                         out_path,
    const lsdtt_smoothed_cdf&                  scdf,
    const std::vector<lsdtt_distribution_fit>& fits
);

// -----------------------------------------------------------------------------
// write_dist_fits_csv
//
// Writes one row per fitted distribution family with columns:
//   family, param1, param2, ks_statistic, max_divergence_x,
//   n_values, ks_critical_value_05, valid
//
// ks_critical_value_05 = 1.36 / sqrt(n_values) is the KS critical value at
// the 5 % significance level for a known-parameter test.  Because parameters
// are estimated from the same sample (MLE), the actual critical value is
// somewhat smaller (Lilliefors correction), so this column is a conservative
// upper bound.
// -----------------------------------------------------------------------------
void write_dist_fits_csv(
    lsdtt_data_loader&                         loader,
    const std::string&                         out_path,
    const std::vector<lsdtt_distribution_fit>& fits
);

} // namespace lsdtt
