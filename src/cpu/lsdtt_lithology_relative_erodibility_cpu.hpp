// lsdtt_lithology_relative_erodibility_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Relative erodibility of lithologic units from their channel steepness.
//
// With the stream power law E = K A^m S^n and uniform erosion at steady state
// (E = U everywhere), the channel steepness k_sn = S A^theta (theta = m/n) of a
// unit is k_sn,i = (E / K_i)^(1/n), so
//
//     K_i / K_ref = (k_sn,ref / k_sn,i)^n
//
// for any reference unit, and with a known erosion rate E (m/yr)
//
//     K_i = E / (A0^m k_i^n),   k_i = dz/dchi with chi referenced to A0,
//
// which is E / k_sn,i^n when A0 = 1 m^2. k_i is the unit's REACH-MEDIAN k_sn:
// the median of k_sn_ols over its qualifying lithologic reaches
// (lsdtt_lithology_channel_stats_cpu), the same value as reach_k_sn_median of
// _chi_ksn_by_lithology_mn_*.csv. m = (m/n) n.
//
// Units of K (map-area convention, A in m^2, S dimensionless, E in m/yr):
// K = E / (A^m S^n) has units m^(1 - 2m) yr^-1, whatever n is (S carries no
// unit). Written as "m^(<1 - 2m>) yr^-1", e.g. m^(0.1) yr^-1 for m = 0.45.
//
// Reference unit. Named, or by default the unit with the most qualifying
// reaches; ties to the larger channel length (channel_length_m of the k_sn
// summary), then the lower code. Unmapped (-1), all (-2), mixed (-3) and arms
// disagree (-4) are never a reference and get no ratio row.
//
// Interval: the ratio is resampled, never formed from the ends of two
// intervals. Each bootstrap replicate resamples the reaches of every unit and
// of the reference, takes each unit's replicate median k*, and forms
// (k*_ref / k*_i)^n and E / (A0^m k*_i^n) in that replicate; the 95 % interval
// is the 2.5 and 97.5 % quantiles of those replicate values
// (compute_quantile_interpolated_in_place). Every n of the list uses the same
// replicate medians.
//   * cluster reach: within a replicate each unit's qualifying reaches are
//     drawn with replacement, as many as it has, independently per unit (the
//     reaches of two units are disjoint sets, so there is nothing to draw
//     jointly). The reference row uses the reference's own replicate median on
//     both sides, so its ratio is 1 with the interval [1, 1].
//   * cluster basin: within a replicate the basins holding a qualifying reach
//     of any lithologic unit are drawn with replacement, as many as there are,
//     and every unit pools its reaches from the drawn basins: a basin draw
//     brings the reaches of the reference and of unit i together. A replicate
//     in which a unit (or the reference) drew no reach has no value for that
//     unit and is left out; n_bootstrap counts the replicates used for the
//     ratio (the K interval uses those in which the unit itself drew reaches).
// A unit with fewer than 5 qualifying reaches (low_sample, as in the k_sn
// summary) gets its point values and NaN intervals, and draws nothing. When
// the reference is low_sample every ratio interval is NaN. A replicate median
// <= 0 gives no ratio or K in that replicate.
//
// Random stream: lithology_relative_erodibility_bootstrap_engine(seed, m/n) is
// one std::mt19937_64 per m/n, seeded from std::seed_seq{seed, round(m/n *
// 10000), 1}: per m/n like the k_sn summary's engine, so an interval does not
// depend on the other m/n values of a run, but a different stream, so the k_sn
// summary's intervals are unchanged by switching this on. Draws go replicate
// by replicate; inside a replicate: cluster reach, the reference then the
// other units in ascending code; cluster basin, one basin draw.
//
// Registry of the per-m/n bootstrap streams of the k_sn / erodibility
// analyses (a new analysis takes the next free tag, so switching it on never
// moves another analysis's intervals):
//   * k_sn by lithology summary:  seed_seq{seed, round(m/n * 1e4)}
//   * relative erodibility:       seed_seq{seed, round(m/n * 1e4), 1}
//   * erodibility calibration:    seed_seq{seed, round(m/n * 1e4), 2}
//     (from per-basin erosion rates; reserved)
//
// The K formula (stream_power_erodibility_from_erosion_rate) and the property
// table writer (write_lithology_erodibility_property_table) are shared with the
// erodibility calibration from per-basin erosion rates; the relative-K writer
// is a thin wrapper over the shared one.
//
// Draft property table: when E is given, a table in the format
// load_lithology_property_table reads (lsdtt-lithology-erodibility):
//   # lsdtt lithology property table v1
//   # erodibility_convention: map_area
//   # m_exponent: <m>
//   # n_exponent: <n>
//   # provenance: <DEM, map, m/n, A0, E, reference, "assumes uniform erosion
//                  and steady state">
//   code,name,erodibility,reference
// with one row per unit that has a finite K > 0 (the `reference` column
// describes where the value came from). Units with no K are listed in a
// `# not_estimated:` comment line; lsdtt-lithology-erodibility refuses a map
// with a code missing from the table, so rows for them must be added by hand.
//
// Interpretation conditions (stated, not tested here): uniform erosion rate,
// steady state, k_sn of a reach near a contact mixing both units, transient
// knickzones, and the number of reaches per unit.
//
// Nothing here prints.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include "lsdtt_lithology_channel_stats_cpu.hpp"
#include "lsdtt_lithology_map_cpu.hpp"

namespace lsdtt
{

struct lsdtt_relative_erodibility_reference
{
    std::int32_t code = -1;
    std::string  name;
    std::string  reason;   // why it was chosen, for the run report
};

// The reference unit from the reach-era k_sn summary (rows with a unit code
// >= 0 only). requested_name "" = the default rule; otherwise the unit with
// exactly that name. Throws std::invalid_argument for a name that is no unit
// (listing the unit names), for unmapped / all / mixed / arms disagree, and
// when the chosen unit has no qualifying reach.
[[nodiscard]] lsdtt_relative_erodibility_reference choose_relative_erodibility_reference_unit(
    const lsdtt_chi_ksn_lithology_summary& summary, const std::string& requested_name);

// Checks a requested reference name against the units of a run before any
// analysis (the CLI's early FATAL). Throws as above, without the
// qualifying-reach check. "" is always accepted.
void check_relative_erodibility_reference_name(const lsdtt_lithology_units_on_grid& units,
                                               const std::string& requested_name);

struct lsdtt_relative_erodibility_options
{
    std::string reference_unit;                 // "" = default rule
    std::vector<double> n_values{1.0};          // each finite and > 0
    double erosion_rate_m_per_yr = 0.0;         // <= 0 = not given (no K, no table)
    int bootstrap_replicates = 1000;            // 0 = no interval
    lsdtt_lithology_bootstrap_cluster bootstrap_cluster = lsdtt_lithology_bootstrap_cluster::reach;
};

struct lsdtt_relative_erodibility_row
{
    static constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    std::int32_t lith_code = -1;
    std::string  lith_name;
    double n = 1.0, m = 0.0;                    // m = (m/n) n
    std::size_t n_reaches_qualifying = 0;
    double reach_k_sn_median = nan;
    double k_ratio = nan, k_ratio_ci_low = nan, k_ratio_ci_high = nan;   // K_i / K_ref
    int    n_bootstrap = 0;                     // replicates behind the ratio interval
    bool   low_sample = false;                  // fewer than 5 qualifying reaches
    bool   is_reference = false;
    double K_abs = nan, K_abs_ci_low = nan, K_abs_ci_high = nan;   // with E only
};

struct lsdtt_relative_erodibility_result
{
    double m_over_n = 0.0, A_0 = 1.0;
    lsdtt_relative_erodibility_reference reference;
    lsdtt_lithology_bootstrap_cluster bootstrap_cluster = lsdtt_lithology_bootstrap_cluster::reach;
    double erosion_rate_m_per_yr = 0.0;         // 0 = not given
    std::vector<double> n_values;
    std::vector<lsdtt_relative_erodibility_row> rows;   // n in list order, then units ascending
};

// The stream for one m/n (see the header comment).
[[nodiscard]] std::mt19937_64 lithology_relative_erodibility_bootstrap_engine(std::int64_t seed,
                                                                              double m_over_n);

// Relative (and, with E, absolute) erodibility per lithologic unit and n.
// `summary` is the reach-era summary of `tables` (its m/n, A0, reach medians,
// counts and channel lengths are used). Throws std::invalid_argument for a bad
// n, a bad reference (choose_relative_erodibility_reference_unit) or a negative
// erosion rate.
[[nodiscard]] lsdtt_relative_erodibility_result compute_lithology_relative_erodibility(
    const lsdtt_chi_lithology_tables& tables, const lsdtt_chi_ksn_lithology_summary& summary,
    const lsdtt_relative_erodibility_options& opts, std::mt19937_64& rng);

// "m^(0.1) yr^-1" for m = 0.45: the units of map-area K = E / (A^m S^n).
[[nodiscard]] std::string lithology_erodibility_units_string(double m);

// Map-area stream-power erodibility from an erosion rate and a channel
// steepness: K = E / (A0^m k_sn^n), with E in m/yr and k_sn = dz/dchi for chi
// referenced to A0 (m^2). NaN when E <= 0, k_sn <= 0 or k_sn is not finite.
[[nodiscard]] double stream_power_erodibility_from_erosion_rate(double E_m_per_yr, double k_sn,
                                                                double m, double n, double A0);

// n in a file name: "1", "2", "1p5".
[[nodiscard]] std::string relative_erodibility_n_suffix(double n);

// _chi_relative_erodibility_by_lithology_mn_*.csv: one row per unit and n.
// Columns lith_code, lith_name, m_over_n, A_0, n, m, reference_unit,
// n_reaches_qualifying, reach_k_sn_median, k_ratio, k_ratio_ci_low,
// k_ratio_ci_high, n_bootstrap, bootstrap_cluster, low_sample, then with E
// erosion_rate_m_per_yr, K_abs, K_abs_ci_low, K_abs_ci_high, K_units.
void write_lithology_relative_erodibility_csv(const std::string& path,
                                              const lsdtt_relative_erodibility_result& r);

// One unit of a lithology property table. K that is not finite and > 0, or a
// code outside 0..65534, puts the unit on the "# not_estimated:" line instead.
struct lsdtt_lithology_erodibility_table_row
{
    std::int32_t code = -1;
    std::string  name;
    double       K = std::numeric_limits<double>::quiet_NaN();
    std::string  reference;   // the reference column: where the value came from
};

// The shared property-table writer (the format load_lithology_property_table
// reads). Writes the version, convention, m_exponent and n_exponent lines,
// "# provenance: <provenance>", one "# <line>" per extra comment line, the
// "# not_estimated: <code> <name>; ... (<not_estimated_note>)" line when a
// unit has no K, then the code,name,erodibility,reference rows in the given
// order. Throws std::invalid_argument when two written units' names are equal
// ignoring case (the reader refuses that), std::runtime_error when the file
// cannot be written. Returns the number of unit rows written.
std::size_t write_lithology_erodibility_property_table(
    const std::string& path, double m, double n, const std::string& provenance,
    const std::vector<std::string>& extra_comment_lines,
    const std::vector<lsdtt_lithology_erodibility_table_row>& units,
    const std::string& not_estimated_note);

struct lsdtt_relative_erodibility_provenance
{
    std::string dem;              // file name of the DEM
    std::string lithology_map;    // file name of the map
    std::string field;            // lithology_field_name (vectors)
    std::string crosswalk;        // "" = none
};

// The draft property table for one n (needs E > 0; throws otherwise, and when
// two units' names are equal ignoring case, which the table reader refuses).
// Returns the number of unit rows written.
std::size_t write_lithology_relative_erodibility_property_table(
    const std::string& path, const lsdtt_relative_erodibility_result& r, double n,
    const lsdtt_relative_erodibility_provenance& provenance);

} // namespace lsdtt
