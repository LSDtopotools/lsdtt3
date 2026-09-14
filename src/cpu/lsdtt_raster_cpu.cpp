// lsdtt_raster_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_raster_cpu.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt {

// -----------------------------------------------------------------------------
// Constructors
// -----------------------------------------------------------------------------

lsdtt_raster_cpu::lsdtt_raster_cpu(const lsdtt_data_loader& loader)
{
    lsdtt_raster_tile tile = loader.read_full();

    set_from_vectors(
        std::move(tile.data),
        std::move(tile.mask),
        tile.meta
    );
}

lsdtt_raster_cpu::lsdtt_raster_cpu(const lsdtt_raster_tile& tile)
    : meta_(tile.meta)
{
    data_ = lsdtt_array2d_f(
        tile.meta.n_rows,
        tile.meta.n_cols,
        lsdtt_array2d_f::from_ptr,
        tile.data.data()
    );

    mask_ = lsdtt_array2d_u8(
        tile.meta.n_rows,
        tile.meta.n_cols,
        lsdtt_array2d_u8::from_ptr,
        tile.mask.data()
    );
}

lsdtt_raster_cpu::lsdtt_raster_cpu(
    lsdtt_array2d_f&&        data,
    lsdtt_array2d_u8&&       mask,
    const lsdtt_raster_meta& meta
)
    : data_(std::move(data))
    , mask_(std::move(mask))
    , meta_(meta)
{
}

lsdtt_raster_cpu::lsdtt_raster_cpu(
    std::vector<float>&&     data,
    std::vector<uint8_t>&&   mask,
    const lsdtt_raster_meta& meta
)
    : data_(meta.n_rows, meta.n_cols, std::move(data))
    , mask_(meta.n_rows, meta.n_cols, std::move(mask))
    , meta_(meta)
{
}

lsdtt_raster_cpu::lsdtt_raster_cpu(lsdtt_raster_result&& result)
    : data_(result.meta.n_rows, result.meta.n_cols, std::move(result.data))
    , mask_(result.meta.n_rows, result.meta.n_cols, std::move(result.mask))
    , meta_(result.meta)
{
}

// -----------------------------------------------------------------------------
// Loader-to-backend transfer
// -----------------------------------------------------------------------------

void lsdtt_raster_cpu::set_from_vectors(
    std::vector<float>&&        data,
    std::vector<std::uint8_t>&& mask,
    const lsdtt_raster_meta&    meta
)
{
    meta_ = meta;

    data_ = lsdtt_array2d_f(
        meta.n_rows,
        meta.n_cols,
        std::move(data)
    );

    mask_ = lsdtt_array2d_u8(
        meta.n_rows,
        meta.n_cols,
        std::move(mask)
    );
}

// -----------------------------------------------------------------------------
// Private helpers
// -----------------------------------------------------------------------------

void lsdtt_raster_cpu::validate_storage_sizes(const char* function_name) const
{
    if (data_.size() != mask_.size())
    {
        throw std::runtime_error(
            std::string("lsdtt_raster_cpu::") + function_name +
            "(): data and mask sizes do not match"
        );
    }
}

std::int64_t lsdtt_raster_cpu::checked_cell_count(const char* function_name) const
{
    validate_storage_sizes(function_name);

    if (data_.size() >
        static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
    {
        throw std::runtime_error(
            std::string("lsdtt_raster_cpu::") + function_name +
            "(): raster contains more cells than fit in int64_t"
        );
    }

    return static_cast<std::int64_t>(data_.size());
}

// -----------------------------------------------------------------------------
// Whole-raster operations over valid pixels (mask != 0)
// -----------------------------------------------------------------------------

std::int64_t lsdtt_raster_cpu::count_valid_pixels() const
{
    const std::int64_t n_cells = checked_cell_count("count_valid_pixels");
    const std::uint8_t* mask_ptr = mask_.data();

    std::int64_t n_valid = 0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:n_valid)
#endif
    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        if (mask_ptr[static_cast<std::size_t>(i)] != 0)
        {
            ++n_valid;
        }
    }

    return n_valid;
}

std::int64_t lsdtt_raster_cpu::assert_mask()
{
    const std::int64_t n_cells = checked_cell_count("assert_mask");

    const float*        data_ptr = data_.data();
    std::uint8_t*       mask_ptr = mask_.data();
    const float         nodata   = static_cast<float>(meta_.no_data_value);
    const bool          nodata_is_nan = std::isnan(static_cast<double>(nodata));

    std::int64_t n_changed = 0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:n_changed)
#endif
    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        const std::size_t idx = static_cast<std::size_t>(i);
        const float z = data_ptr[idx];

        const bool is_nodata =
            nodata_is_nan
                ? std::isnan(static_cast<double>(z))
                : (z == nodata || std::isnan(static_cast<double>(z)));

        if (is_nodata && mask_ptr[idx] != 0)
        {
            mask_ptr[idx] = 0;
            ++n_changed;
        }
    }

    return n_changed;
}


lsdtt_raster_cpu make_uniform_raster(const lsdtt_raster_cpu& reference,
                                     float                   value,
                                     MaskedFill              outside)
{
    const lsdtt_raster_meta& meta = reference.get_meta();
    const float fill = (outside == MaskedFill::Value)
                     ? value
                     : static_cast<float>(meta.no_data_value);

    lsdtt_array2d<float>   d(meta.n_rows, meta.n_cols, fill);
    lsdtt_array2d<uint8_t> m(meta.n_rows, meta.n_cols, 0);

    const auto& ref_mask = reference.get_mask_array();
    for (int r = 0; r < meta.n_rows; ++r)
        for (int c = 0; c < meta.n_cols; ++c)
            if (ref_mask(r, c)) { d(r, c) = value; m(r, c) = 1; }

    return lsdtt_raster_cpu(std::move(d), std::move(m), meta);
}

} // namespace lsdtt

// =============================================================================
// Vertical units
//
// See the header for why none of this is ever applied automatically.
// =============================================================================
namespace lsdtt
{

vertical_unit parse_vertical_unit(const std::string& name)
{
    // Trim the ENDS only, and lower-case. Interior whitespace is left in place
    // so that "metre s" stays a typo rather than collapsing into "metres" --
    // this parser is the last line of defence against a silent 3.28x error.
    std::size_t b = 0, e = name.size();
    while (b < e && std::isspace(static_cast<unsigned char>(name[b]))) ++b;
    while (e > b && std::isspace(static_cast<unsigned char>(name[e - 1]))) --e;

    std::string k;
    k.reserve(e - b);
    for (std::size_t i = b; i < e; ++i)
        k += static_cast<char>(std::tolower(static_cast<unsigned char>(name[i])));

    if (k == "metres" || k == "meters" || k == "metre" || k == "meter" || k == "m")
        return vertical_unit::metres;
    if (k == "feet" || k == "foot" || k == "ft")
        return vertical_unit::feet;
    if (k == "us_survey_feet" || k == "us_survey_foot" || k == "us_feet" ||
        k == "usft" || k == "us_ft")
        return vertical_unit::us_survey_feet;

    // Deliberately fatal. Falling back to metres on an unrecognised string
    // would turn a typo into a silent 3.28x error in every derived quantity.
    throw std::runtime_error(
        "unrecognised vertical_units '" + name +
        "' - expected one of: metres, feet, us_survey_feet");
}

namespace
{
    double metres_per(vertical_unit u)
    {
        switch (u)
        {
            case vertical_unit::feet:           return 0.3048;
            case vertical_unit::us_survey_feet: return 1200.0 / 3937.0;
            case vertical_unit::metres:         break;
        }
        return 1.0;
    }
}

lsdtt_raster_cpu convert_vertical_units_to_metres(
    const lsdtt_raster_cpu& dem,
    vertical_unit           from)
{
    lsdtt_raster_cpu out(dem);
    const double f = metres_per(from);
    if (f == 1.0) return out;

    // Only valid pixels are scaled: nodata is a sentinel, not an elevation, and
    // scaling it would move it off the value the metadata declares.
    std::vector<float>&         d = out.get_data();
    const std::vector<uint8_t>& m = out.get_mask();
    for (std::size_t i = 0; i < d.size(); ++i)
        if (m[i]) d[i] = static_cast<float>(static_cast<double>(d[i]) * f);

    return out;
}

vertical_plausibility_report check_vertical_plausibility(
    const lsdtt_raster_cpu& dem,
    double                  median_slope_limit_deg)
{
    vertical_plausibility_report rep;

    const int    nr = dem.n_rows();
    const int    nc = dem.n_cols();
    const double cs = dem.cell_size();
    if (nr < 3 || nc < 3 || cs <= 0.0) return rep;

    const std::vector<float>&   d = dem.get_data();
    const std::vector<uint8_t>& m = dem.get_mask();

    std::vector<float> slopes;
    slopes.reserve(static_cast<std::size_t>(nr) * static_cast<std::size_t>(nc) / 4);

    float zmin =  std::numeric_limits<float>::max();
    float zmax = -std::numeric_limits<float>::max();

    // Interior pixels with all four cardinal neighbours valid: a centred
    // difference at the edge of a nodata block measures the block boundary,
    // not the terrain, and those are exactly the pixels that would bias this.
    for (int r = 1; r < nr - 1; ++r)
    {
        for (int c = 1; c < nc - 1; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(r) * nc + c;
            if (!m[i]) continue;

            if (d[i] < zmin) zmin = d[i];
            if (d[i] > zmax) zmax = d[i];

            const std::size_t n = i - nc, s = i + nc, w = i - 1, e = i + 1;
            if (!m[n] || !m[s] || !m[w] || !m[e]) continue;

            const double dzdx = (static_cast<double>(d[e]) - d[w]) / (2.0 * cs);
            const double dzdy = (static_cast<double>(d[s]) - d[n]) / (2.0 * cs);
            slopes.push_back(static_cast<float>(std::hypot(dzdx, dzdy)));
        }
    }

    if (slopes.empty()) return rep;

    const auto pct = [&slopes](double p) -> double
    {
        const std::size_t k = static_cast<std::size_t>(
            p / 100.0 * static_cast<double>(slopes.size() - 1));
        std::nth_element(slopes.begin(), slopes.begin() + k, slopes.end());
        return static_cast<double>(slopes[k]);
    };

    const auto deg = [](double g) { return std::atan(g) * 180.0 / M_PI; };

    rep.n_sampled                = slopes.size();
    const double med             = pct(50.0);
    rep.median_slope_deg         = deg(med);
    rep.p90_slope_deg            = deg(pct(90.0));
    rep.p99_slope_deg            = deg(pct(99.0));
    rep.median_slope_deg_if_feet = deg(med * 0.3048);
    rep.relief_m                 = static_cast<double>(zmax - zmin);
    rep.implausible              = rep.median_slope_deg > median_slope_limit_deg;

    return rep;
}

vertical_plausibility_report warn_if_vertically_implausible(
    const lsdtt_raster_cpu& dem,
    std::ostream&           out,
    const std::string&      vertical_units_parameter)
{
    const vertical_plausibility_report rep = check_vertical_plausibility(dem);
    if (!rep.implausible) return rep;

    out << "\n"
        << "===============================================================\n"
        << "WARNING: this DEM's slopes are not physically plausible.\n"
        << "===============================================================\n"
        << "  median slope   " << rep.median_slope_deg << " degrees\n"
        << "  90th pct       " << rep.p90_slope_deg    << " degrees\n"
        << "  99th pct       " << rep.p99_slope_deg    << " degrees\n"
        << "  relief         " << rep.relief_m
        << " (in whatever unit the elevations are in)\n\n"
        << "Soil and regolith fail at roughly 30-35 degrees, so a landscape\n"
        << "whose MEDIAN pixel is steeper than that cannot be soil-mantled\n"
        << "terrain. Isolated cliffs do not move a median.\n\n"
        << "THE USUAL CAUSE is elevations in FEET on a metre grid, which no\n"
        << "raster format records and no reader can detect: US state lidar is\n"
        << "often delivered as NAVD88 feet on a UTM metre grid. If that is the\n"
        << "case here, the median slope is really "
        << rep.median_slope_deg_if_feet << " degrees, and EVERY\n"
        << "slope, curvature, relief and hillslope gradient in this run is\n"
        << "too large by a factor of 3.28.\n\n"
        << "lsdtt3 WILL NOT CONVERT THIS FOR YOU. The vertical unit cannot be\n"
        << "inferred from the data -- a genuinely alpine catchment in metres\n"
        << "and a gentle one in feet overlap in every statistic -- so guessing\n"
        << "would corrupt correct data. Check the vertical datum in the data\n"
        << "provider's metadata. If it is feet, set:\n\n"
        << "    " << vertical_units_parameter << ": feet\n\n"
        << "and if it is genuinely metres, ignore this warning.\n"
        << "===============================================================\n\n";

    return rep;
}

} // namespace lsdtt
