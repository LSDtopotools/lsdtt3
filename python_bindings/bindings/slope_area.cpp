// slope_area.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// lsdtt-slope-area's routines (src/cpu/lsdtt_slope_area_cpu): the same C++
// calls and parameter names as the CLI (plan_slope_area_gradient_and_low_
// concavity.md, step 7). The CLI's routing-surface conditioning (carve, fill)
// is done in Python with the bound conditioning routines; the conditioned
// pixels are then passed as `exclude`, as the CLI does.
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cpu/lsdtt_slope_area_cpu.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

void require_same_grid(const lsdtt::lsdtt_raster_meta& a, const lsdtt::lsdtt_raster_meta& b, const char* what)
{
    if (a.n_rows != b.n_rows || a.n_cols != b.n_cols || a.x_min != b.x_min || a.y_min != b.y_min ||
        a.cell_size != b.cell_size)
        throw std::invalid_argument(std::string(what) + " must share the grid");
}

void require_flow_grid(const FlowInfo& flow, const Raster& r, const char* what)
{
    require_same_grid(r.native().get_meta(), flow.native().meta(), what);
}

// exclude: None, or a Raster whose nonzero valid cells are excluded (e.g. the
// result of conditioned_pixel_mask).
std::vector<std::uint8_t> exclude_mask(nb::object exclude, const lsdtt::lsdtt_raster_meta& meta)
{
    std::vector<std::uint8_t> out;
    if (exclude.is_none()) return out;
    const Raster& r = nb::cast<const Raster&>(exclude);
    require_same_grid(r.native().get_meta(), meta, "exclude");
    const auto& d = r.native().get_data();
    const auto& m = r.native().get_mask();
    out.assign(d.size(), 0);
    for (std::size_t i = 0; i < d.size(); ++i) out[i] = (m[i] && std::isfinite(d[i]) && d[i] != 0.0f) ? 1 : 0;
    return out;
}

Raster conditioned_pixel_mask(const Raster& input_dem, const Raster& routing_dem)
{
    require_same_grid(input_dem.native().get_meta(), routing_dem.native().get_meta(), "routing_dem");
    const auto codes = lsdtt::conditioned_pixel_mask(input_dem.native(), routing_dem.native());
    const auto& meta = input_dem.native().get_meta();
    std::vector<float> data(codes.begin(), codes.end());
    std::vector<std::uint8_t> mask(codes.size(), 1);
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta));
}

nb::object slope_area_bins(const Raster& drainage_area, const Raster& slope, int n_log_bins, nb::object exclude)
{
    require_same_grid(drainage_area.native().get_meta(), slope.native().get_meta(), "slope");
    if (n_log_bins < 1) throw std::invalid_argument("slope_area_bins: n_log_bins must be at least 1");
    const auto ex = exclude_mask(exclude, drainage_area.native().get_meta());
    const auto b = lsdtt::slope_area_log_bins(drainage_area.native(), slope.native(), n_log_bins,
                                              ex.empty() ? nullptr : &ex);
    std::vector<std::int64_t> counts(b.counts.begin(), b.counts.end());
    nb::dict d;
    d["drainage_area"] = b.geometric_centre;
    d["d05_slope"] = b.quantiles[0];
    d["d25_slope"] = b.quantiles[1];
    d["d50_slope"] = b.quantiles[2];
    d["d75_slope"] = b.quantiles[3];
    d["d95_slope"] = b.quantiles[4];
    d["number_of_pixels"] = counts;
    return nb::module_::import_("pandas").attr("DataFrame")(d);
}

Raster d8_flow_path_gradient(const FlowInfo& flow, const Raster& routing_dem, const Raster& drainage_area,
                             double window_metres, nb::object exclude, int n_threads)
{
    require_flow_grid(flow, routing_dem, "routing_dem");
    require_flow_grid(flow, drainage_area, "drainage_area");
    const auto ex = exclude_mask(exclude, flow.native().meta());
    lsdtt::lsdtt_d8_flow_path_gradient_options o;
    o.window_metres = window_metres;
    o.exclude = ex.empty() ? nullptr : &ex;
    o.n_threads = n_threads;
    return Raster(lsdtt::d8_flow_path_gradient(flow.native(), routing_dem.native(), drainage_area.native(), o));
}

lsdtt::lsdtt_low_concavity_result run_low_concavity(const FlowInfo& flow, const Raster& drainage_area,
                                                    const Raster& gradient, std::int64_t threshold_contributing_pixels,
                                                    int window_pixels, int min_pixels, double max_abs_concavity,
                                                    const std::vector<std::uint8_t>& ex, int n_threads)
{
    require_flow_grid(flow, drainage_area, "drainage_area");
    require_flow_grid(flow, gradient, "gradient");
    if (threshold_contributing_pixels < 1)
        throw std::invalid_argument("low_concavity: threshold_contributing_pixels must be at least 1");
    lsdtt::lsdtt_low_concavity_options lo;
    lo.threshold_contributing_pixels = static_cast<std::uint64_t>(threshold_contributing_pixels);
    lo.window_pixels = window_pixels;
    lo.min_pixels = min_pixels;
    lo.max_abs_concavity = max_abs_concavity;
    lo.exclude = ex.empty() ? nullptr : &ex;
    lo.n_threads = n_threads;
    return lsdtt::find_low_concavity(flow.native(), drainage_area.native(), gradient.native(), lo);
}

nb::object low_concavity(const FlowInfo& flow, const Raster& drainage_area, const Raster& gradient,
                         std::int64_t threshold_contributing_pixels, int window_pixels, int min_pixels,
                         double max_abs_concavity, nb::object exclude, int n_threads)
{
    const auto ex = exclude_mask(exclude, flow.native().meta());
    const auto r = run_low_concavity(flow, drainage_area, gradient, threshold_contributing_pixels, window_pixels,
                                     min_pixels, max_abs_concavity, ex, n_threads);
    const auto& meta = flow.native().meta();
    std::vector<double> x(r.flat.size()), y(r.flat.size());
    std::vector<std::int64_t> flag(r.flag.begin(), r.flag.end());
    for (std::size_t k = 0; k < r.flat.size(); ++k)
    {
        const std::int64_t row = r.flat[k] / meta.n_cols, col = r.flat[k] % meta.n_cols;
        x[k] = meta.x_min + (static_cast<double>(col) + 0.5) * meta.cell_size;
        y[k] = meta.y_min + (static_cast<double>(meta.n_rows - row) - 0.5) * meta.cell_size;
    }
    nb::dict d;
    d["flat"] = r.flat;
    d["x"] = x;
    d["y"] = y;
    d["drainage_area_m2"] = r.area_m2;
    d["slope"] = r.slope;
    d["theta"] = r.theta;
    d["theta_se"] = r.theta_se;
    d["r2"] = r.r2;
    d["n_pixels"] = r.n_pixels;
    d["log10_area_range"] = r.log10_area_range;
    d["n_junctions_in_window"] = r.n_junctions;
    d["low_concavity"] = flag;
    return nb::module_::import_("pandas").attr("DataFrame")(d);
}

nb::object low_concavity_segments(const FlowInfo& flow, const Raster& drainage_area, const Raster& gradient,
                                  std::int64_t threshold_contributing_pixels, int window_pixels, int min_pixels,
                                  double max_abs_concavity, int min_run_pixels, nb::object exclude,
                                  nb::object include, int n_threads)
{
    const auto ex = exclude_mask(exclude, flow.native().meta());
    const auto r = run_low_concavity(flow, drainage_area, gradient, threshold_contributing_pixels, window_pixels,
                                     min_pixels, max_abs_concavity, ex, n_threads);
    const auto inc = exclude_mask(include, flow.native().meta());     // nonzero = included
    const auto segs = lsdtt::build_low_concavity_segments(flow.native(), r, min_run_pixels,
                                                          inc.empty() ? nullptr : &inc);
    const auto& meta = flow.native().meta();
    auto xy = [&](std::int64_t f) {
        const std::int64_t row = f / meta.n_cols, col = f % meta.n_cols;
        return nb::make_tuple(meta.x_min + (static_cast<double>(col) + 0.5) * meta.cell_size,
                              meta.y_min + (static_cast<double>(meta.n_rows - row) - 0.5) * meta.cell_size);
    };
    nb::object line_string = nb::module_::import_("shapely.geometry").attr("LineString");
    nb::list geometry;
    std::vector<std::int64_t> id, link, flag, npx, order, up, down;
    std::vector<double> length, mth, mse, msl, amin, amax;
    for (const auto& s : segs)
    {
        nb::list coords;
        for (const auto f : s.flats) coords.append(xy(f));
        if (s.end_flat >= 0) coords.append(xy(s.end_flat));
        if (nb::len(coords) < 2) coords.append(xy(s.flats.front()));
        geometry.append(line_string(coords));
        id.push_back(s.segment_id);
        link.push_back(s.link_id);
        flag.push_back(s.flag);
        npx.push_back(static_cast<std::int64_t>(s.flats.size()));
        length.push_back(s.length_m);
        mth.push_back(s.median_theta);
        mse.push_back(s.median_theta_se);
        msl.push_back(s.median_slope);
        amin.push_back(s.min_area_m2);
        amax.push_back(s.max_area_m2);
        order.push_back(s.stream_order);
        up.push_back(s.flats.front());
        down.push_back(s.flats.back());
    }
    nb::dict d;
    d["segment_id"] = id;
    d["link_id"] = link;
    d["low_concavity"] = flag;
    d["n_pixels"] = npx;
    d["length_m"] = length;
    d["median_theta"] = mth;
    d["median_theta_se"] = mse;
    d["median_slope"] = msl;
    d["min_area_m2"] = amin;
    d["max_area_m2"] = amax;
    d["stream_order"] = order;
    d["upstream_flat"] = up;
    d["downstream_flat"] = down;
    nb::object ctor = nb::module_::import_("geopandas").attr("GeoDataFrame");
    return meta.crs_wkt.empty() ? ctor(d, nb::arg("geometry") = geometry)
                                : ctor(d, nb::arg("geometry") = geometry, nb::arg("crs") = meta.crs_wkt);
}

} // namespace

void bind_slope_area(nb::module_& module, nb::class_<FlowInfo>& flow)
{
    module.def("conditioned_pixel_mask", &conditioned_pixel_mask, nb::arg("input_dem"), nb::arg("routing_dem"),
        nb::sig("def conditioned_pixel_mask(input_dem: Raster, routing_dem: Raster) -> Raster"), R"doc(
Mark the pixels a DEM conditioning changed: 0 untouched, 1 filled, 2 carved.

lsdtt-slope-area's conditioned-pixel mask. A pixel is filled where
``routing_dem`` (the carved and / or filled surface flow is routed on) is
above ``input_dem``, and carved where it is below. lsdtt-slope-area leaves
both out of its measurements (``exclude_carved_and_filled_pixels``); pass the
result as ``exclude`` to ``slope_area_bins``, ``FlowInfo.d8_flow_path_gradient``
and ``FlowInfo.low_concavity`` to do the same.

Parameters
----------
input_dem : Raster
    The DEM before any carve or fill.
routing_dem : Raster
    The conditioned DEM, on the same grid.

Returns
-------
Raster
    Codes 0, 1, 2 (float), every pixel valid.

Examples
--------
>>> routed = lsdtt3.fill_pits(lsdtt3.breach_least_cost(dem))
>>> conditioned = lsdtt3.conditioned_pixel_mask(dem, routed)
)doc");

    module.def("slope_area_bins", &slope_area_bins, nb::arg("drainage_area"), nb::arg("slope"), nb::kw_only(),
        nb::arg("n_log_bins") = 20, nb::arg("exclude") = nb::none(),
        nb::sig("def slope_area_bins(drainage_area: Raster, slope: Raster, *, n_log_bins: int = 20, "
                "exclude: Raster | None = None) -> pandas.DataFrame"), R"doc(
Bin slope by drainage area in logarithmic bins: the lsdtt-slope-area lambda CSV.

The same routine as lsdtt-slope-area. Pixels valid in both rasters with area
> 0 are binned in ``n_log_bins`` bins of equal width in ln(area); empty bins
are dropped. ``drainage_area`` in the result is the geometric bin centre.

Parameters
----------
drainage_area : Raster
    Drainage area (m^2), e.g. from ``FlowGraph`` (D-infinity) or ``FlowInfo``.
slope : Raster
    Gradient (m/m) on the same grid, e.g. polyfit slope or
    ``FlowInfo.d8_flow_path_gradient``.
n_log_bins : int, default=20, keyword-only
    Number of bins (>= 1).
exclude : Raster or None, default=None, keyword-only
    Pixels with a nonzero value are left out, e.g. ``conditioned_pixel_mask``.

Returns
-------
pandas.DataFrame
    ``drainage_area, d05_slope, d25_slope, d50_slope, d75_slope, d95_slope,
    number_of_pixels``, one row per non-empty bin.
)doc");

    flow.def("d8_flow_path_gradient", &d8_flow_path_gradient, nb::arg("routing_dem"), nb::arg("drainage_area"),
        nb::kw_only(), nb::arg("window_metres") = 150.0, nb::arg("exclude") = nb::none(), nb::arg("n_threads") = 0,
        nb::sig("def d8_flow_path_gradient(self, routing_dem: Raster, drainage_area: Raster, *, "
                "window_metres: float = 150.0, exclude: Raster | None = None, n_threads: int = 0) -> Raster"), R"doc(
Gradient at every pixel along its D8 flow path, over a window centred on it.

lsdtt-slope-area ``slope_algorithm_choice: d8_flow_path``, same routine. The
window of total flow length ``window_metres`` runs downstream along this
FlowInfo's receivers and upstream along the D8 donor of largest
``drainage_area``. The gradient is the least-squares slope of elevation
against flow distance; the valid cells must span half the window less one
diagonal step, so a one-sided window (divide, outlet) can qualify.

Parameters
----------
routing_dem : Raster
    The DEM this FlowInfo was built on.
drainage_area : Raster
    Drainage area (m^2) used to pick the upstream path (D-infinity or D8).
window_metres : float, default=150.0, keyword-only
    Total window length (m), at least four cell sizes.
exclude : Raster or None, default=None, keyword-only
    Nonzero pixels (e.g. carved / filled) are walked through but not fitted,
    and their own gradient is nodata.
n_threads : int, default=0, keyword-only
    OpenMP threads (0 = default); the result does not depend on it.

Returns
-------
Raster
    Gradient (m/m); nodata where the fit is not usable.
)doc");

    flow.def("low_concavity", &low_concavity, nb::arg("drainage_area"), nb::arg("gradient"), nb::kw_only(),
        nb::arg("threshold_contributing_pixels") = 1000, nb::arg("window_pixels") = 10, nb::arg("min_pixels") = 5,
        nb::arg("max_abs_concavity") = 0.15, nb::arg("exclude") = nb::none(), nb::arg("n_threads") = 0,
        nb::sig("def low_concavity(self, drainage_area: Raster, gradient: Raster, *, "
                "threshold_contributing_pixels: int = 1000, window_pixels: int = 10, min_pixels: int = 5, "
                "max_abs_concavity: float = 0.15, exclude: Raster | None = None, n_threads: int = 0) "
                "-> pandas.DataFrame"), R"doc(
Fit log S against log A downstream of every channel pixel: the low-concavity finder.

lsdtt-slope-area ``run_low_concavity_detection``, same routine and names
(``low_concavity_window_pixels`` -> ``window_pixels`` etc.). At each channel
pixel (``threshold_contributing_pixels``) the pixel and the next ones down
the D8 receivers (``window_pixels`` in all; through junctions) are fitted as
log10 S = a + b log10 A; theta = -b. A pixel is low concavity if
|theta| < ``max_abs_concavity``.

Parameters
----------
drainage_area : Raster
    Drainage area (m^2).
gradient : Raster
    The D8 flow-path gradient (``d8_flow_path_gradient``), as the CLI uses.
threshold_contributing_pixels : int, default=1000, keyword-only
    Channel pixels have at least this many D8 contributing pixels.
window_pixels : int, default=10, keyword-only
    Pixels in the downstream window (>= 2).
min_pixels : int, default=5, keyword-only
    Fewest usable (area, gradient) pairs for a valid fit (3 to window_pixels).
max_abs_concavity : float, default=0.15, keyword-only
    Threshold on |theta|.
exclude : Raster or None, default=None, keyword-only
    Nonzero pixels (carved / filled) are left out of the fits; an excluded
    channel pixel gets no fit.
n_threads : int, default=0, keyword-only
    OpenMP threads.

Returns
-------
pandas.DataFrame
    One row per channel pixel, the columns of the CLI's point table:
    ``flat, x, y, drainage_area_m2, slope, theta, theta_se, r2, n_pixels,
    log10_area_range, n_junctions_in_window, low_concavity`` (1, 0, or -1 for
    no valid fit). ``theta_se`` is the standard error of theta.
)doc");

    flow.def("low_concavity_segments", &low_concavity_segments, nb::arg("drainage_area"), nb::arg("gradient"),
        nb::kw_only(), nb::arg("threshold_contributing_pixels") = 1000, nb::arg("window_pixels") = 10,
        nb::arg("min_pixels") = 5, nb::arg("max_abs_concavity") = 0.15, nb::arg("min_run_pixels") = 1,
        nb::arg("exclude") = nb::none(), nb::arg("include") = nb::none(), nb::arg("n_threads") = 0,
        nb::sig("def low_concavity_segments(self, drainage_area: Raster, gradient: Raster, *, "
                "threshold_contributing_pixels: int = 1000, window_pixels: int = 10, min_pixels: int = 5, "
                "max_abs_concavity: float = 0.15, min_run_pixels: int = 1, exclude: Raster | None = None, "
                "include: Raster | None = None, n_threads: int = 0) -> geopandas.GeoDataFrame"), R"doc(
The channel network as segments marked low concavity or not: the CLI's segments file.

Runs ``low_concavity`` with the same arguments, then cuts the channel network
into links (source or junction to the next junction or outlet) and splits
each where the flag changes. A flagged run shorter than ``min_run_pixels`` is
relabelled 0 and merged. Pixels outside ``include`` (nonzero = inside, e.g. a
basin mask) are left out. Each line ends at the first pixel below its last.

Parameters
----------
drainage_area, gradient, threshold_contributing_pixels, window_pixels, min_pixels, max_abs_concavity, exclude, n_threads
    As ``low_concavity``.
min_run_pixels : int, default=1, keyword-only
    Shortest flagged run kept as a low-concavity segment (>= 1).
include : Raster or None, default=None, keyword-only
    Nonzero pixels are included; None includes all.

Returns
-------
geopandas.GeoDataFrame
    LineStrings with ``segment_id, link_id, low_concavity, n_pixels,
    length_m, median_theta, median_theta_se, median_slope, min_area_m2,
    max_area_m2, stream_order`` (Strahler), ``upstream_flat,
    downstream_flat``.
)doc");
}

} // namespace lsdtt3_python
