// lsdtt_lem_input_grid_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_lem_input_grid_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "../core/lsdtt_raster_warp.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"

namespace lsdtt
{

namespace
{

// Grid coordinates round-trip exactly through GeoTIFF, so this only absorbs
// print-rounding in hand-built metas; it is far below any real misalignment.
constexpr double k_grid_tol_m = 1e-6;

// The 1 % distance rule (A11, Q6).
constexpr double k_max_scale_error = 0.01;

std::string num(double v)
{
    std::ostringstream s;
    s << std::setprecision(12) << v;
    return s.str();
}

std::string describe_grid(const lsdtt_raster_meta& m)
{
    std::ostringstream s;
    s << std::setprecision(12) << m.n_rows << "x" << m.n_cols << ", cell " << m.cell_size
      << ", origin (" << m.x_min << ", " << m.y_min << "), CRS "
      << (m.crs_wkt.empty() ? std::string("none") : std::string("set"));
    return s.str();
}

std::string describe_extent(const lsdtt_lem_clip_extent& e)
{
    return "[" + num(e.x_min) + ", " + num(e.x_max) + "] x [" + num(e.y_min) + ", " +
           num(e.y_max) + "]";
}

lsdtt_lem_clip_extent extent_of(const lsdtt_raster_meta& m)
{
    lsdtt_lem_clip_extent e;
    e.x_min = m.x_min;
    e.y_min = m.y_min;
    e.x_max = m.x_min + m.n_cols * m.cell_size;
    e.y_max = m.y_min + m.n_rows * m.cell_size;
    return e;
}

bool is_web_mercator(const lsdtt_raster_meta& m)
{
    // EPSG:3857 and its deprecated EPSG alias 3785. Non-EPSG aliases (900913,
    // ESRI 102100 / 102113) describe the same CRS and so compare equal to it.
    for (int code : {3857, 3785})
    {
        lsdtt_raster_meta w;
        try { w.set_crs_from_epsg(code); }
        catch (const std::exception&) { continue; }
        if (crs_equal_semantically(m.crs_wkt, w.crs_wkt)) return true;
    }
    return false;
}

// A clip source must be a projected raster in metres that preserves distance.
void check_clip_source_crs(const lsdtt_raster_meta& m, const std::string& label)
{
    if (m.crs_wkt.empty())
        throw std::invalid_argument(label + " raster has no CRS; a clip point in the source CRS "
                                    "needs one (use relabel_whole_raster for a raster with no "
                                    "georeference)");
    if (m.is_geographic())
        throw std::invalid_argument(label + " raster is in a geographic (degree) CRS; the clip "
                                    "coordinates and cell size need a projected CRS in metres. "
                                    "Reproject it first (e.g. lsdtt-raster-preprocessing)");
    if (!m.is_in_metres())
        throw std::invalid_argument(label + " raster CRS unit is '" + m.linear_unit_name() +
                                    "', not metres; reproject it to a metric CRS first");
    check_lem_projection_distortion(m, label);
}

lsdtt_raster_cpu relabel_to(std::vector<float>&& data, std::vector<std::uint8_t>&& mask,
                            const lsdtt_raster_meta& model_grid)
{
    // Pixels flagged nodata take the model grid's nodata value.
    const float nd = static_cast<float>(model_grid.no_data_value);
    for (std::size_t i = 0; i < data.size(); ++i)
        if (!mask[i]) data[i] = nd;
    return lsdtt_raster_cpu(std::move(data), std::move(mask), model_grid);
}

void check_resampling(const lsdtt_lem_input_raster_placement_options& o, const std::string& label)
{
    const std::string& r = o.resampling_method;
    const bool ok = o.categorical ? (r == "near" || r == "mode")
                                  : (r == "near" || r == "bilinear" || r == "average" || r == "mode");
    if (!ok)
        throw std::invalid_argument(
            label + ": input_raster_resampling_method '" + r + "' is not allowed for " +
            (o.categorical ? "a categorical raster (use near or mode)"
                           : "a continuous raster (use near, bilinear, average or mode)"));
}

} // namespace

// -----------------------------------------------------------------------------
// Parsing
// -----------------------------------------------------------------------------
lem_model_grid_source parse_lem_model_grid_source(const std::string& name)
{
    if (name == "auto") return lem_model_grid_source::automatic;
    if (name == "parameters") return lem_model_grid_source::parameters;
    if (name == "dem") return lem_model_grid_source::dem;
    if (name == "erodibility_raster") return lem_model_grid_source::erodibility_raster;
    if (name == "lithology_raster") return lem_model_grid_source::lithology_raster;
    throw std::invalid_argument("unknown model_grid_source '" + name +
                                "': expected auto, parameters, dem, erodibility_raster or "
                                "lithology_raster");
}

lem_input_raster_placement parse_lem_input_raster_placement(const std::string& name)
{
    if (name == "as_georeferenced") return lem_input_raster_placement::as_georeferenced;
    if (name == "warp_to_model_grid") return lem_input_raster_placement::warp_to_model_grid;
    if (name == "relabel_whole_raster") return lem_input_raster_placement::relabel_whole_raster;
    if (name == "relabel_centre_clip") return lem_input_raster_placement::relabel_centre_clip;
    if (name == "relabel_lower_left_clip") return lem_input_raster_placement::relabel_lower_left_clip;
    throw std::invalid_argument("unknown input_raster_placement '" + name +
                                "': expected as_georeferenced, warp_to_model_grid, "
                                "relabel_whole_raster, relabel_centre_clip or "
                                "relabel_lower_left_clip");
}

std::string lem_input_raster_placement_name(lem_input_raster_placement p)
{
    switch (p)
    {
    case lem_input_raster_placement::as_georeferenced: return "as_georeferenced";
    case lem_input_raster_placement::warp_to_model_grid: return "warp_to_model_grid";
    case lem_input_raster_placement::relabel_whole_raster: return "relabel_whole_raster";
    case lem_input_raster_placement::relabel_centre_clip: return "relabel_centre_clip";
    case lem_input_raster_placement::relabel_lower_left_clip: return "relabel_lower_left_clip";
    }
    return "unknown";
}

// -----------------------------------------------------------------------------
// Source grids
// -----------------------------------------------------------------------------
lsdtt_lem_source_grid make_lem_source_grid(const lsdtt_raster_meta& meta, int n_bands,
                                           const std::string& label)
{
    if (meta.n_rows <= 0 || meta.n_cols <= 0)
        throw std::invalid_argument(label + " raster grid is empty (" +
                                    std::to_string(meta.n_rows) + "x" +
                                    std::to_string(meta.n_cols) + ")");
    if (!(meta.cell_size > 0.0))
        throw std::invalid_argument(label + " raster cell size must be > 0 (got " +
                                    num(meta.cell_size) + ")");
    lsdtt_lem_source_grid g;
    g.meta = meta;
    g.n_bands = std::max(1, n_bands);
    return g;
}

lsdtt_lem_source_grid describe_lem_source_grid(const lsdtt_data_loader& loader,
                                               const std::string& label)
{
    return make_lem_source_grid(loader.get_meta(), loader.band_count(), label);
}

// -----------------------------------------------------------------------------
// Strict alignment
// -----------------------------------------------------------------------------
std::vector<std::string> check_lem_field_grid_alignment(const lsdtt_raster_meta& field_meta,
                                                        const lsdtt_raster_meta& model_grid,
                                                        const std::string&       field_label)
{
    const bool geometry_ok = same_grid(field_meta, model_grid, k_grid_tol_m);
    const bool both_crs = !field_meta.crs_wkt.empty() && !model_grid.crs_wkt.empty();
    const bool crs_ok = !both_crs || crs_equal_semantically(field_meta.crs_wkt, model_grid.crs_wkt);

    if (!geometry_ok || !crs_ok)
    {
        std::ostringstream msg;
        msg << field_label << " raster grid (" << describe_grid(field_meta)
            << ") does not match the model grid (" << describe_grid(model_grid) << ")";
        if (geometry_ok && !crs_ok) msg << ": the CRS differs";
        msg << ". Use input_raster_placement (warp_to_model_grid or a relabel_* option), "
               "or align it with lsdtt-raster-calculator.";
        throw std::invalid_argument(msg.str());
    }

    std::vector<std::string> warnings;
    if (!both_crs)
        warnings.push_back(field_label + " raster: no CRS on " +
                           (field_meta.crs_wkt.empty() ? std::string("the raster")
                                                       : std::string("the model grid")) +
                           "; accepted because its rows, cols, cell size and origin match the "
                           "model grid");
    return warnings;
}

// -----------------------------------------------------------------------------
// Projection distortion
// -----------------------------------------------------------------------------
void check_lem_projection_distortion(const lsdtt_raster_meta& meta, const std::string& label)
{
    if (meta.crs_wkt.empty() || meta.is_geographic()) return;
    const double err = meta.max_linear_scale_error();
    if (is_web_mercator(meta))
        throw std::invalid_argument(
            label + " CRS is Web Mercator (EPSG:3857), which does not preserve distances (scale "
            "factor " + num(1.0 + err) + " here); reproject the raster to a conformal metric "
            "CRS such as the local UTM zone");
    if (err > k_max_scale_error)
    {
        std::ostringstream msg;
        msg << label << " CRS distorts distances by " << std::setprecision(3) << 100.0 * err
            << "% at this location; the LEM needs distances within 1 % of true. Reproject the "
               "raster, for example to the local UTM zone";
        throw std::invalid_argument(msg.str());
    }
}

// -----------------------------------------------------------------------------
// Model grids
// -----------------------------------------------------------------------------
lsdtt_raster_meta make_lem_synthetic_grid(int n_rows, int n_cols, double cell_size,
                                          const std::string& synthetic_crs_wkt,
                                          std::vector<std::string>* warnings)
{
    if (n_rows <= 0 || n_cols <= 0 || !(cell_size > 0.0))
        throw std::invalid_argument("synthetic model grid needs n_rows, n_cols > 0 and a cell "
                                    "size > 0 (got " + std::to_string(n_rows) + "x" +
                                    std::to_string(n_cols) + ", cell " + num(cell_size) + ")");
    lsdtt_raster_meta m;
    m.n_rows = n_rows;
    m.n_cols = n_cols;
    m.cell_size = cell_size;
    m.x_min = 0.0;
    m.y_min = 0.0;
    m.no_data_value = -9999.0;
    m.crs_wkt = synthetic_crs_wkt;
    if (synthetic_crs_wkt.empty())
    {
        // Watopia (UTM zone 58S), centred on the grid: dims must be set first.
        m.set_default_synthetic_georeference();
        if (m.crs_wkt.empty() && warnings)
            warnings->push_back("could not build the default synthetic georeference (Watopia, "
                                "UTM zone 58S); the synthetic surface has no CRS");
    }
    return m;
}

lsdtt_raster_meta make_lem_model_grid_from_raster(const lsdtt_lem_source_grid& source,
                                                  const std::string& label)
{
    const lsdtt_raster_meta& s = source.meta;
    if (s.n_rows <= 0 || s.n_cols <= 0)
        throw std::invalid_argument(label + " raster grid is empty");
    if (s.crs_wkt.empty())
        throw std::invalid_argument(label + " raster has no CRS, so it cannot place the model on "
                                    "a real landscape; use input_raster_placement: "
                                    "relabel_whole_raster to use its dimensions on the "
                                    "synthetic georeference");
    if (s.is_geographic())
        throw std::invalid_argument(label + " raster is in a geographic (degree) CRS; the LEM grid "
                                    "must be projected in metres. Reproject it first (e.g. "
                                    "lsdtt-raster-preprocessing)");
    if (!s.is_in_metres())
        throw std::invalid_argument(label + " raster CRS unit is '" + s.linear_unit_name() +
                                    "', not metres; the LEM grid must be in metres");
    check_lem_projection_distortion(s, label);
    lsdtt_raster_meta m = s;
    m.no_data_value = -9999.0;
    return m;
}

lsdtt_raster_meta make_lem_model_grid_from_raster_dimensions(const lsdtt_lem_source_grid& source,
                                                             const std::string& synthetic_crs_wkt,
                                                             std::vector<std::string>* warnings)
{
    return make_lem_synthetic_grid(source.meta.n_rows, source.meta.n_cols, source.meta.cell_size,
                                   synthetic_crs_wkt, warnings);
}

std::size_t apply_lem_source_footprint(lsdtt_raster_cpu& dem,
                                       const std::vector<std::uint8_t>& source_valid,
                                       const std::string& label)
{
    std::uint8_t* mask = dem.get_mask_array().data();
    float* z = dem.get_data_array().data();
    const std::size_t n = static_cast<std::size_t>(dem.get_meta().n_rows) * dem.get_meta().n_cols;
    if (source_valid.size() != n)
        throw std::invalid_argument(label + " raster has " + std::to_string(source_valid.size()) +
                                    " pixels but the model grid has " + std::to_string(n));
    const float nd = static_cast<float>(dem.get_meta().no_data_value);
    std::size_t removed = 0, remaining = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (mask[i] && !source_valid[i])
        {
            mask[i] = 0;
            z[i] = nd;
            ++removed;
        }
        if (mask[i]) ++remaining;
    }
    if (remaining == 0)
        throw std::invalid_argument(label + " raster has no valid pixels; the model domain would "
                                    "be empty");
    return removed;
}

// -----------------------------------------------------------------------------
// Clip windows
// -----------------------------------------------------------------------------
lsdtt_lem_clip_extent compute_lem_clip_extent(lem_input_raster_placement placement,
                                              double easting, double northing, int n_rows,
                                              int n_cols, double cell_size)
{
    if (!std::isfinite(easting) || !std::isfinite(northing))
        throw std::invalid_argument("a clip placement needs input_raster_clip_easting_metres and "
                                    "input_raster_clip_northing_metres");
    lsdtt_lem_clip_extent e;
    const double w = n_cols * cell_size, h = n_rows * cell_size;
    if (placement == lem_input_raster_placement::relabel_centre_clip)
    {
        e.x_min = easting - 0.5 * w;
        e.y_min = northing - 0.5 * h;
    }
    else if (placement == lem_input_raster_placement::relabel_lower_left_clip)
    {
        e.x_min = easting;
        e.y_min = northing;
    }
    else
        throw std::invalid_argument("compute_lem_clip_extent: placement must be "
                                    "relabel_centre_clip or relabel_lower_left_clip");
    e.x_max = e.x_min + w;
    e.y_max = e.y_min + h;
    return e;
}

void check_lem_clip_extent_within(const lsdtt_lem_clip_extent& window,
                                  const lsdtt_lem_clip_extent& source_extent,
                                  const std::string& label)
{
    const double tol = k_grid_tol_m;
    if (window.x_min < source_extent.x_min - tol || window.x_max > source_extent.x_max + tol ||
        window.y_min < source_extent.y_min - tol || window.y_max > source_extent.y_max + tol)
        throw std::invalid_argument(label + " clip window " + describe_extent(window) +
                                    " runs off the source extent " +
                                    describe_extent(source_extent));
}

lsdtt_lem_clip_window compute_lem_clip_window(const lsdtt_lem_source_grid& source,
                                              lem_input_raster_placement placement,
                                              double easting, double northing, int n_rows,
                                              int n_cols, double cell_size,
                                              const std::string& label)
{
    if (placement != lem_input_raster_placement::relabel_centre_clip &&
        placement != lem_input_raster_placement::relabel_lower_left_clip)
        throw std::invalid_argument("compute_lem_clip_window: placement must be "
                                    "relabel_centre_clip or relabel_lower_left_clip");
    if (!std::isfinite(easting) || !std::isfinite(northing))
        throw std::invalid_argument(label + ": a clip placement needs "
                                    "input_raster_clip_easting_metres and "
                                    "input_raster_clip_northing_metres");
    const lsdtt_raster_meta& s = source.meta;
    check_clip_source_crs(s, label);

    lsdtt_lem_clip_window w;
    w.n_rows = n_rows;
    w.n_cols = n_cols;
    w.source_cell_size = s.cell_size;
    w.target_cell_size = cell_size;
    const lsdtt_lem_clip_extent src_ext = extent_of(s);

    if (std::fabs(s.cell_size / cell_size - 1.0) <= 1e-9)
    {
        const double dx = s.cell_size;
        // Pixels counted from the south-west corner; a point on a grid line
        // belongs to the pixel to its north / east (hence floor, with a
        // rounding guard for coordinates that sit exactly on a line).
        const long ic = static_cast<long>(std::floor((easting - s.x_min) / dx + 1e-9));
        const long jb = static_cast<long>(std::floor((northing - s.y_min) / dx + 1e-9));
        const long row_of_jb = static_cast<long>(s.n_rows) - 1 - jb;   // top-row indexing
        long col0, row0;
        if (placement == lem_input_raster_placement::relabel_centre_clip)
        {
            col0 = ic - n_cols / 2;
            row0 = row_of_jb - n_rows / 2;
        }
        else
        {
            col0 = ic;
            row0 = row_of_jb - (n_rows - 1);
        }
        w.pixel_aligned = true;
        w.extent.x_min = s.x_min + col0 * dx;
        w.extent.x_max = w.extent.x_min + n_cols * dx;
        w.extent.y_max = s.y_min + (static_cast<double>(s.n_rows) - row0) * dx;
        w.extent.y_min = w.extent.y_max - n_rows * dx;
        if (col0 < 0 || row0 < 0 || col0 + n_cols > s.n_cols || row0 + n_rows > s.n_rows)
            throw std::invalid_argument(label + " clip window " + describe_extent(w.extent) +
                                        " runs off the source extent " + describe_extent(src_ext));
        w.col0 = static_cast<int>(col0);
        w.row0 = static_cast<int>(row0);
        return w;
    }

    w.extent = compute_lem_clip_extent(placement, easting, northing, n_rows, n_cols, cell_size);
    check_lem_clip_extent_within(w.extent, src_ext, label);
    return w;
}

// -----------------------------------------------------------------------------
// prepare_lem_input_raster
// -----------------------------------------------------------------------------
lsdtt_lem_prepared_input_raster prepare_lem_input_raster(
    const std::string& path, const lsdtt_raster_meta& model_grid,
    const lsdtt_lem_input_raster_placement_options& opts, const std::string& label)
{
    // The method matters for a warp, and for a continuous raster resampled in a
    // clip; a categorical clip always uses near, whatever was asked for.
    const bool clip = opts.placement == lem_input_raster_placement::relabel_centre_clip ||
                      opts.placement == lem_input_raster_placement::relabel_lower_left_clip;
    if (opts.placement == lem_input_raster_placement::warp_to_model_grid ||
        (clip && !opts.categorical))
        check_resampling(opts, label);
    lsdtt_data_loader loader;
    loader.open_raster(path);
    const lsdtt_lem_source_grid src = describe_lem_source_grid(loader, label);
    const lsdtt_raster_meta& s = src.meta;

    lsdtt_lem_prepared_input_raster out;
    const lem_input_raster_placement p = opts.placement;
    const bool relabel = p == lem_input_raster_placement::relabel_whole_raster ||
                         p == lem_input_raster_placement::relabel_centre_clip ||
                         p == lem_input_raster_placement::relabel_lower_left_clip;
    if (relabel && !opts.model_grid_is_synthetic)
        throw std::invalid_argument(label + ": input_raster_placement " +
                                    lem_input_raster_placement_name(p) +
                                    " discards the raster's georeference, which is only allowed "
                                    "on a synthetic model grid; a loaded DEM's georeference is "
                                    "kept (use as_georeferenced or warp_to_model_grid)");

    switch (p)
    {
    case lem_input_raster_placement::as_georeferenced:
    {
        out.warnings = check_lem_field_grid_alignment(s, model_grid, label);
        lsdtt_raster_tile t = loader.read_full();
        out.raster = relabel_to(std::move(t.data), std::move(t.mask), model_grid);
        return out;
    }
    case lem_input_raster_placement::warp_to_model_grid:
    {
        if (s.crs_wkt.empty() || model_grid.crs_wkt.empty())
            throw std::invalid_argument(label + ": warp_to_model_grid needs a CRS on both the "
                                        "raster and the model grid");
        lsdtt_warp_options wo;
        wo.resampling_method = opts.categorical && opts.resampling_method != "mode"
                                   ? std::string("near") : opts.resampling_method;
        wo.error_if_no_overlap = true;
        wo.n_threads = opts.n_threads;
        lsdtt_raster_tile t;
        try { t = warp_raster_to_grid_in_memory(path, model_grid, wo); }
        catch (const std::runtime_error& e)
        {
            throw std::invalid_argument(label + ": could not warp onto the model grid: " +
                                        e.what());
        }
        out.raster = relabel_to(std::move(t.data), std::move(t.mask), model_grid);
        out.warped = true;
        out.resampled = true;
        out.notes.push_back(label + ": warped onto the model grid (" + describe_grid(model_grid) +
                            ") with " + wo.resampling_method + " resampling" +
                            (crs_equal_semantically(s.crs_wkt, model_grid.crs_wkt)
                                 ? std::string("")
                                 : std::string(", changing CRS")));
        return out;
    }
    case lem_input_raster_placement::relabel_whole_raster:
    {
        if (s.n_rows != model_grid.n_rows || s.n_cols != model_grid.n_cols ||
            std::fabs(s.cell_size / model_grid.cell_size - 1.0) > 1e-9)
            throw std::invalid_argument(
                label + ": relabel_whole_raster needs the raster's rows, cols and cell size to "
                "equal the model grid's (" + std::to_string(s.n_rows) + "x" +
                std::to_string(s.n_cols) + " at " + num(s.cell_size) + " vs " +
                std::to_string(model_grid.n_rows) + "x" + std::to_string(model_grid.n_cols) +
                " at " + num(model_grid.cell_size) + "). It never rescales; use a clip "
                "placement to take a window at the model's cell size");
        lsdtt_raster_tile t = loader.read_full();
        out.raster = relabel_to(std::move(t.data), std::move(t.mask), model_grid);
        out.relabelled = true;
        out.notes.push_back(label + ": relabelled onto the model grid (values copied unchanged); "
                            "the raster's georeference (origin (" + num(s.x_min) + ", " +
                            num(s.y_min) + "), CRS " + s.crs_name() + ") is discarded");
        return out;
    }
    case lem_input_raster_placement::relabel_centre_clip:
    case lem_input_raster_placement::relabel_lower_left_clip:
    {
        const lsdtt_lem_clip_window w = compute_lem_clip_window(
            src, p, opts.clip_easting, opts.clip_northing, model_grid.n_rows, model_grid.n_cols,
            model_grid.cell_size, label);
        if (w.pixel_aligned)
        {
            lsdtt_tile_window tw;
            tw.x_off = tw.read_x_off = w.col0;
            tw.y_off = tw.read_y_off = w.row0;
            tw.tile_width = tw.read_width = w.n_cols;
            tw.tile_height = tw.read_height = w.n_rows;
            tw.halo = 0;
            lsdtt_raster_tile t = loader.read_tile(tw);
            out.raster = relabel_to(std::move(t.data), std::move(t.mask), model_grid);
        }
        else
        {
            lsdtt_raster_meta target;
            target.n_rows = model_grid.n_rows;
            target.n_cols = model_grid.n_cols;
            target.cell_size = model_grid.cell_size;
            target.x_min = w.extent.x_min;
            target.y_min = w.extent.y_min;
            target.no_data_value = s.no_data_value;
            target.crs_wkt = s.crs_wkt;
            lsdtt_warp_options wo;
            wo.resampling_method = opts.categorical ? std::string("near") : opts.resampling_method;
            wo.n_threads = opts.n_threads;
            lsdtt_raster_tile t = warp_raster_to_grid_in_memory(path, target, wo);
            out.raster = relabel_to(std::move(t.data), std::move(t.mask), model_grid);
            out.resampled = true;
            out.notes.push_back(label + ": resampled from cell size " + num(s.cell_size) +
                                " to " + num(model_grid.cell_size) + " with " +
                                wo.resampling_method + " resampling");
        }
        out.relabelled = true;
        out.notes.push_back(label + ": clip window " + describe_extent(w.extent) + " in " +
                            s.crs_name() + " relabelled onto the model grid; the raster's "
                            "georeference is discarded");
        return out;
    }
    }
    throw std::logic_error("prepare_lem_input_raster: unhandled placement");
}

} // namespace lsdtt
