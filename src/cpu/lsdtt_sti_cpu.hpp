// lsdtt_sti_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt {

std::string normalize_sti_series_code(const std::string& raw_code);
std::string normalize_sti_land_use(const std::string& raw_land_use);

enum class sti_term_mode
{
    ksd,
    spr
};

enum class sti_spr_transform
{
    ln_1_minus_spr,
    ln_spr
};

struct sti_lnkd_options
{
    double depth_cap_cm = 100.0;
    double min_valid_depth_cm = 30.0;
    double ks_unit_scale_to_m_per_day = 1.0;
    double eps = 1.0e-12;
};

struct sti_lnkd_table
{
    std::unordered_map<std::string, double> ln_kd_by_key;
    std::unordered_map<std::string, int> n_profiles_by_key;
};

struct sti_lnkd_build_result
{
    sti_lnkd_table table;
    std::size_t rows_read = 0;
    std::size_t rows_used = 0;
    std::size_t profiles_total = 0;
    std::size_t profiles_valid = 0;
    std::size_t pairs_built = 0;
};

struct sti_spr_table
{
    std::unordered_map<std::string, double> spr_fraction_by_series;
};

struct sti_spr_parse_result
{
    sti_spr_table table;
    std::size_t rows_read = 0;
    std::size_t rows_used = 0;
    std::size_t clamped_rows = 0;
};

struct sti_polygon_string_rasterise_result
{
    lsdtt_raster_cpu raster_codes;
    std::vector<std::pair<std::uint32_t, std::string>> lookup;
};

struct sti_compute_options
{
    sti_term_mode mode = sti_term_mode::ksd;
    sti_spr_transform spr_transform = sti_spr_transform::ln_1_minus_spr;
    double eps = 1.0e-12;
};

struct sti_compute_stats
{
    std::size_t valid_twi_pixels = 0;
    std::size_t valid_output_pixels = 0;
    std::size_t missing_series_code_pixels = 0;
    std::size_t missing_land_use_pixels = 0;
    std::size_t missing_lookup_pixels = 0;
    std::size_t fallback_pixels = 0;
    std::size_t clamped_pixels = 0;
};

struct sti_compute_result
{
    lsdtt_raster_cpu sti_raster;
    lsdtt_raster_cpu term_raster;
    lsdtt_raster_cpu valid_mask_raster;
    sti_compute_stats stats;
    std::unordered_map<std::string, std::string> source_by_series_landuse;
};

struct sti_lookup_row
{
    std::string series_code;
    std::string land_use;
    int n_profiles = 0;
    bool has_ln_kd = false;
    double ln_kd_median = 0.0;
    bool has_spr_prc = false;
    double spr_prc = 0.0;
    bool has_spr_term = false;
    double spr_term = 0.0;
    std::string source;
};

// ---------------------------------------------------------------------------
// Look up ln(Ks*D) for a (series code, land use) pair, falling back to related
// land uses when the exact pair is absent.
//
// Exposed because the distributed-hydrology work needs exactly this lookup to
// build a transmissivity field: `ln_kd` IS ln(T0), so a land-cover scenario
// changes transmissivity by RE-KEYING this table to the measured profile for
// the new land use -- not by applying a multiplier. Writing a second lookup
// with its own fallback order would be a silent way for the two to disagree.
//
// NOTE the units trap: ln_kd is ln of m.cm/day, not m^2/day. Convert with
// lsdtt_topmodel_cpu.hpp's sti_ln_kd_to_ln_t0_m2_per_day/hour() before using it
// as TOPMODEL's ln_t0.
//
// `out_source` is "direct" for an exact match, "fallback_other_landuse" when a
// related land use supplied the value, and "missing" when nothing did.
// Returns false when nothing was found.
// ---------------------------------------------------------------------------
bool lookup_ln_kd_with_fallback(
    const sti_lnkd_table& table,
    const std::string& series_code,
    const std::string& land_use,
    double& out_ln_kd,
    std::string& out_source,
    int& out_profile_count);

sti_lnkd_build_result build_sti_lnkd_table_from_csv(
    const std::string& csv_path,
    const sti_lnkd_options& options);

sti_spr_parse_result load_sti_spr_table_from_csv(
    const std::string& csv_path,
    double eps = 1.0e-12);

bool sti_check_xlsx_available(const std::string& xlsx_path);

sti_polygon_string_rasterise_result rasterise_polygon_string_field(
    const lsdtt_raster_cpu& reference_raster,
    const std::string& vector_path,
    const std::string& layer_name,
    const std::string& field_name,
    int oversample = 8);

sti_compute_result compute_sti_from_rasters(
    const lsdtt_raster_cpu& twi_raster,
    const sti_polygon_string_rasterise_result& series_code_raster,
    const sti_polygon_string_rasterise_result& land_use_raster,
    const sti_lnkd_table& lnkd_table,
    const sti_spr_table& spr_table,
    const sti_compute_options& options);

std::vector<sti_lookup_row> build_sti_lookup_rows(
    const std::unordered_map<std::string, std::string>& source_by_series_landuse,
    const sti_lnkd_table& lnkd_table,
    const sti_spr_table& spr_table,
    sti_spr_transform spr_transform,
    double eps = 1.0e-12);

void write_sti_lookup_csv(
    const std::vector<sti_lookup_row>& rows,
    const std::string& out_csv_path);

} // namespace lsdtt
