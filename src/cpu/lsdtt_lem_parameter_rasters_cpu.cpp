// lsdtt_lem_parameter_rasters_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_lem_parameter_rasters_cpu.hpp"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "lsdtt_lem_input_grid_cpu.hpp"

namespace lsdtt
{

namespace
{

std::string num(double v)
{
    std::ostringstream s;
    s << std::setprecision(6) << v;
    return s.str();
}

bool differs(float a, float b) { return std::fabs(static_cast<double>(a) - b) > 1e-6; }

std::string pixel(std::size_t i, int n_cols)
{
    const std::size_t c = static_cast<std::size_t>(n_cols);
    return "(" + std::to_string(i / c) + ", " + std::to_string(i % c) + ")";
}

std::string k_must_be_positive(const std::string& label, std::size_t i, int n_cols, float k)
{
    std::ostringstream s;
    s << label << ": K must be > 0; pixel " << pixel(i, n_cols) << " has K = " << k;
    return s.str();
}

} // namespace

lsdtt_erodibility_calibration_check check_lem_erodibility_calibration(
    const lsdtt_erodibility_calibration&    stated,
    const lsdtt_erodibility_run_convention& run,
    const std::string&                      source_label)
{
    lsdtt_erodibility_calibration_check out;
    out.run_units = describe_erodibility_units(run.m_exponent, run.divide_area_by_dx);

    if (stated.m_exponent < 0.0f)
        out.warnings.push_back(source_label + ": the m exponent K was calibrated for is not "
                               "stated (m_exponent); assuming the run's m = " +
                               num(run.m_exponent));
    else if (differs(stated.m_exponent, run.m_exponent))
        throw std::invalid_argument(source_label + ": K was calibrated for m_exponent = " +
                                    num(stated.m_exponent) + " but the run uses m = " +
                                    num(run.m_exponent) + "; K is tied to m");

    if (stated.n_exponent < 0.0f)
        out.warnings.push_back(source_label + ": the n exponent K was calibrated for is not "
                               "stated (n_exponent); assuming the run's n = " +
                               num(run.n_exponent));
    else if (differs(stated.n_exponent, run.n_exponent))
        throw std::invalid_argument(source_label + ": K was calibrated for n_exponent = " +
                                    num(stated.n_exponent) + " but the run uses n = " +
                                    num(run.n_exponent) + "; K is tied to n");

    const double m = static_cast<double>(run.m_exponent);

    if (!stated.area_convention)
    {
        if (run.divide_area_by_dx)
            out.warnings.push_back(
                source_label + ": the area convention of K is not stated (area_convention); "
                "taken as specific area at this grid's cell size. If this K is in map-area "
                "units it is too small by a factor dx^m = " +
                num(std::pow(run.cell_size, m)));
        else
            out.warnings.push_back(source_label + ": the area convention of K is not stated "
                                   "(area_convention); taken as map area, as the run uses");
        return out;
    }

    const erodibility_area_convention conv = *stated.area_convention;
    if (conv == erodibility_area_convention::map_area)
    {
        if (run.divide_area_by_dx)
        {
            out.conversion_factor = erodibility_area_conversion_factor(
                run.m_exponent, conv, -1.0, true, run.cell_size);
            out.notes.push_back("converted map-area K to specific-area K at cell size " +
                                num(run.cell_size) + ": factor dx^m = " +
                                num(out.conversion_factor));
        }
        return out;
    }

    // specific_area stated
    if (!(stated.cell_size > 0.0))
    {
        if (!run.divide_area_by_dx)
            throw std::invalid_argument(
                source_label + ": K is in specific-area units but its calibration cell size "
                "is not stated (calibration_cell_size), so it cannot be converted to the "
                "run's map-area convention");
        out.warnings.push_back(source_label + ": specific-area K with no calibration cell "
                               "size stated (calibration_cell_size); assumed calibrated at "
                               "this grid's cell size " + num(run.cell_size));
        return out;
    }

    if (!run.divide_area_by_dx)
    {
        out.conversion_factor = erodibility_area_conversion_factor(
            run.m_exponent, conv, stated.cell_size, false, run.cell_size);
        out.notes.push_back("converted specific-area K at calibration cell size " +
                            num(stated.cell_size) + " to map-area K: factor dx^(-m) = " +
                            num(out.conversion_factor));
        return out;
    }

    if (std::fabs(stated.cell_size - run.cell_size) > 1e-6 * run.cell_size)
    {
        out.conversion_factor = erodibility_area_conversion_factor(
            run.m_exponent, conv, stated.cell_size, true, run.cell_size);
        out.notes.push_back("converted specific-area K from calibration cell size " +
                            num(stated.cell_size) + " to " + num(run.cell_size) +
                            ": factor (" + num(run.cell_size) + "/" + num(stated.cell_size) +
                            ")^m = " + num(out.conversion_factor));
    }
    return out;
}

// -----------------------------------------------------------------------------
// Erodibility field values
// -----------------------------------------------------------------------------
lsdtt_lem_erodibility_summary summarise_lem_erodibility_field(
    const std::vector<float>& field, const std::vector<std::uint8_t>& dem_mask)
{
    lsdtt_lem_erodibility_summary out;
    double sum = 0.0;
    bool first = true;
    const std::size_t n = std::min(field.size(), dem_mask.size());
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!dem_mask[i]) continue;
        const float k = field[i];
        if (first) { out.min_value = out.max_value = k; first = false; }
        out.min_value = std::min(out.min_value, k);
        out.max_value = std::max(out.max_value, k);
        sum += k;
        ++out.n_valid;
    }
    if (out.n_valid > 0) out.mean_value = sum / static_cast<double>(out.n_valid);
    return out;
}

lsdtt_lem_erodibility_field validate_lem_erodibility_values(
    const std::vector<float>&                  values,
    const std::vector<std::uint8_t>&           values_valid,
    const lsdtt_raster_cpu&                    dem,
    const lsdtt_lem_erodibility_field_options& opts)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::size_t n = static_cast<std::size_t>(meta.n_rows) * meta.n_cols;
    const std::string& label = opts.source_label;
    if (values.size() != n)
        throw std::invalid_argument(label + " has " + std::to_string(values.size()) +
                                    " values but the DEM has " + std::to_string(n) + " pixels");
    if (!values_valid.empty() && values_valid.size() != n)
        throw std::invalid_argument(label + ": validity mask has " +
                                    std::to_string(values_valid.size()) +
                                    " entries but the DEM has " + std::to_string(n) + " pixels");
    if (!(opts.reference_erodibility > 0.0f) || !std::isfinite(opts.reference_erodibility))
        throw std::invalid_argument("erodibility must be > 0 (got " +
                                    num(opts.reference_erodibility) + ")");

    lsdtt_erodibility_run_convention run;
    run.m_exponent = opts.model_m_exponent;
    run.n_exponent = opts.model_n_exponent;
    run.divide_area_by_dx = opts.divide_area_by_dx;
    run.cell_size = meta.cell_size;
    const lsdtt_erodibility_calibration_check cal =
        check_lem_erodibility_calibration(opts.stated, run, label);

    lsdtt_lem_erodibility_field out;
    out.conversion_factor = cal.conversion_factor;
    out.units = cal.run_units;
    out.notes = cal.notes;
    out.warnings = cal.warnings;
    out.values.assign(n, 0.0f);

    const std::vector<std::uint8_t>& dem_mask = dem.get_mask();
    std::size_t n_missing = 0, first_missing = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!dem_mask[i]) continue;
        if (!values_valid.empty() && !values_valid[i])
        {
            if (n_missing == 0) first_missing = i;
            ++n_missing;
        }
    }
    if (n_missing > 0)
        throw std::invalid_argument(label + " has nodata on " + std::to_string(n_missing) +
                                    " DEM-valid pixels; first at " +
                                    pixel(first_missing, meta.n_cols));

    const bool convert = cal.conversion_factor != 1.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!dem_mask[i]) continue;
        const float k_in = values[i];
        const float k = convert ? static_cast<float>(static_cast<double>(k_in) * cal.conversion_factor)
                                : k_in;
        if (!(k > 0.0f) || !std::isfinite(k))
            throw std::invalid_argument(k_must_be_positive(label, i, meta.n_cols, k_in));
        out.values[i] = k;
    }

    const lsdtt_lem_erodibility_summary sm = summarise_lem_erodibility_field(out.values, dem_mask);
    out.n_valid = sm.n_valid;
    out.min_value = sm.min_value;
    out.max_value = sm.max_value;
    out.mean_value = sm.mean_value;

    // A field whose typical value is orders of magnitude from the scalar
    // reference is usually a units mistake (map- vs specific-area, or m).
    if (out.n_valid > 0)
    {
        std::vector<float> v;
        v.reserve(out.n_valid);
        for (std::size_t i = 0; i < n; ++i)
            if (dem_mask[i]) v.push_back(out.values[i]);
        std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
        const double ratio = static_cast<double>(v[v.size() / 2]) / opts.reference_erodibility;
        if (ratio > 100.0 || ratio < 0.01)
            out.warnings.push_back(label + ": median K " + num(v[v.size() / 2]) +
                                   " differs from the scalar erodibility " +
                                   num(opts.reference_erodibility) +
                                   " by more than 100x; check its units (" + out.units + ")");
    }
    return out;
}

lsdtt_lem_erodibility_field validate_lem_erodibility_field(
    const lsdtt_raster_cpu&                    erodibility_raster,
    const lsdtt_raster_cpu&                    dem,
    const lsdtt_lem_erodibility_field_options& opts)
{
    const std::vector<std::string> grid_warnings =
        check_lem_field_grid_alignment(erodibility_raster.get_meta(), dem.get_meta(),
                                       opts.source_label);
    lsdtt_lem_erodibility_field out = validate_lem_erodibility_values(
        erodibility_raster.get_data(), erodibility_raster.get_mask(), dem, opts);
    out.warnings.insert(out.warnings.begin(), grid_warnings.begin(), grid_warnings.end());
    return out;
}

lsdtt_lem_erodibility_field load_lem_erodibility_field(
    const std::string&                              path,
    const lsdtt_raster_cpu&                         dem,
    const lsdtt_lem_input_raster_placement_options& placement,
    const lsdtt_lem_erodibility_field_options&      opts)
{
    const lsdtt_lem_prepared_input_raster prepared =
        prepare_lem_input_raster(path, dem.get_meta(), placement, opts.source_label);
    lsdtt_lem_erodibility_field out = validate_lem_erodibility_field(prepared.raster, dem, opts);
    out.notes.insert(out.notes.begin(), prepared.notes.begin(), prepared.notes.end());
    out.warnings.insert(out.warnings.begin(), prepared.warnings.begin(), prepared.warnings.end());
    return out;
}

void check_lem_erodibility_field(const std::vector<float>& field, const lsdtt_raster_cpu& dem)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::size_t n = static_cast<std::size_t>(meta.n_rows) * meta.n_cols;
    if (field.size() != n)
        throw std::invalid_argument("erodibility_field has " + std::to_string(field.size()) +
                                    " values but the DEM has " + std::to_string(n) + " pixels");
    const std::vector<std::uint8_t>& mask = dem.get_mask();
    for (std::size_t i = 0; i < n; ++i)
        if (mask[i] && (!(field[i] > 0.0f) || !std::isfinite(field[i])))
            throw std::invalid_argument(
                k_must_be_positive("erodibility_field", i, meta.n_cols, field[i]));
}

void scale_lem_erodibility_field(const std::vector<float>& base, float reference_K,
                                 float applied_K, std::vector<float>& out)
{
    if (!(reference_K > 0.0f) || !std::isfinite(reference_K))
        throw std::invalid_argument("scale_lem_erodibility_field: the reference K must be > 0 "
                                    "(got " + num(reference_K) + ")");
    if (!(applied_K > 0.0f) || !std::isfinite(applied_K))
        throw std::invalid_argument("scale_lem_erodibility_field: the applied K must be > 0 "
                                    "(got " + num(applied_K) + ")");
    out.resize(base.size());
    const double a = static_cast<double>(applied_K);
    const double r = static_cast<double>(reference_K);
    for (std::size_t i = 0; i < base.size(); ++i)
    {
        if (base[i] == 0.0f) { out[i] = 0.0f; continue; }
        const double v = (static_cast<double>(base[i]) * a) / r;
        if (v > static_cast<double>(FLT_MAX) || static_cast<float>(v) == 0.0f)
            throw std::range_error("scale_lem_erodibility_field: pixel " + std::to_string(i) +
                                   " K = " + num(base[i]) + " scaled by " + num(applied_K) +
                                   "/" + num(reference_K) + " = " + num(v) +
                                   " is outside the float range");
        out[i] = static_cast<float>(v);
    }
}

} // namespace lsdtt
