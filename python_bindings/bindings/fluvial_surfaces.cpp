// Bindings for the lsdtt-fluvial-surfaces building blocks: the Clubb et al.
// (2017) threshold classifier, the shared threshold picker, level clustering
// and the comparison of surface classes between methods and against mapped
// surfaces. Each wraps the routine the CLI calls, with the CLI's parameter
// names; the method presets themselves are composed in the CLI.

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <core/lsdtt_gdal_init.hpp>
#include <cpu/lsdtt_floodplain_thresholds_cpu.hpp>
#include <cpu/lsdtt_proximity_cpu.hpp>
#include <cpu/lsdtt_surface_level_clustering_cpu.hpp>
#include <cpu/lsdtt_surface_method_comparison_cpu.hpp>
#include <cpu/lsdtt_surface_training_cpu.hpp>
#include <cpu/lsdtt_vector_rasterise_cpu.hpp>
#include <utils/lsdtt_threshold_selection.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

void require_same_grid(const Raster& first, const Raster& second, const char* message)
{
    const auto& a = first.native().get_meta();
    const auto& b = second.native().get_meta();
    if (a.n_rows != b.n_rows || a.n_cols != b.n_cols || a.x_min != b.x_min ||
        a.y_min != b.y_min || a.cell_size != b.cell_size || a.crs_wkt != b.crs_wkt)
        throw std::invalid_argument(message);
}

// Finite values of a Raster's valid cells, or of any array-like.
std::vector<double> sample_values(nb::object values)
{
    std::vector<double> out;
    if (nb::isinstance<Raster>(values)) {
        const auto& r = nb::cast<const Raster&>(values).native();
        for (std::size_t i = 0; i < r.get_data().size(); ++i)
            if (r.get_mask()[i] && std::isfinite(r.get_data()[i])) out.push_back(r.get_data()[i]);
        return out;
    }
    const auto all = nb::cast<std::vector<double>>(
        nb::module_::import_("numpy").attr("asarray")(values, nb::arg("dtype") = "float64")
            .attr("ravel")());
    for (const double v : all) if (std::isfinite(v)) out.push_back(v);
    return out;
}

std::vector<std::uint8_t> nonzero_cells(const Raster& raster)
{
    const auto& r = raster.native();
    std::vector<std::uint8_t> out(r.get_data().size(), 0);
    for (std::size_t i = 0; i < out.size(); ++i)
        out[i] = (r.get_mask()[i] && r.get_data()[i] > 0.5f) ? 1 : 0;
    return out;
}

// A surface-class raster (0 none, 1 floodplain, 2 terrace) as the comparison
// routines take it; invalid cells are 0.
lsdtt::method_class_map class_map(const std::string& name, const Raster& raster)
{
    const auto& r = raster.native();
    lsdtt::method_class_map m;
    m.name = name;
    m.cls.assign(r.get_data().size(), 0);
    for (std::size_t i = 0; i < m.cls.size(); ++i) {
        if (!r.get_mask()[i]) continue;
        const long v = std::lround(r.get_data()[i]);
        if (v < 0 || v > 2)
            throw std::invalid_argument("surface class values must be 0 (none), 1 (floodplain) or 2 (terrace)");
        m.cls[i] = static_cast<std::uint8_t>(v);
    }
    return m;
}

nb::dict select_threshold(nb::object values, const std::string& threshold_algorithm_choice,
                          double fixed_value, double qq_lower_percentile,
                          double qq_upper_percentile, double qq_threshold_fraction)
{
    const std::vector<double> sample = sample_values(values);
    lsdtt::threshold_selection_request request;
    request.fixed_value = fixed_value;
    request.qq_lower_percentile = qq_lower_percentile;
    request.qq_upper_percentile = qq_upper_percentile;
    request.qq_threshold_fraction = qq_threshold_fraction;
    lsdtt::threshold_selection_result result;
    if (threshold_algorithm_choice == "qq_upper_departure") {
        result = lsdtt::select_qq_upper_departure(sample, request);
    } else {
        request.method = lsdtt::parse_threshold_selection_choice(threshold_algorithm_choice);
        result = lsdtt::select_threshold(sample, request);
    }
    nb::dict out;
    out["threshold"] = result.threshold;
    out["found"] = result.found;
    out["n_values"] = sample.size();
    out["best_family"] = result.best_family;
    out["threshold_percentile"] = result.threshold_percentile;
    out["max_cdf_divergence"] = result.max_cdf_divergence;
    return out;
}

class FloodplainThresholdResult {
public:
    explicit FloodplainThresholdResult(lsdtt::floodplain_threshold_result result)
        : result_(std::move(result)) {}
    const lsdtt::floodplain_threshold_result& native() const noexcept { return result_; }
private:
    lsdtt::floodplain_threshold_result result_;
};

FloodplainThresholdResult classify_floodplain_terraces_by_thresholds(
    const Raster& admitted, const Raster& hand, const Raster& channel_mask,
    nb::object gradient, nb::object terrace_field, double max_floodplain_hand_metres,
    double max_floodplain_gradient, bool seed_channel_pixels, bool require_channel_contact,
    int min_band_component_pixels, bool fill_floodplain_holes, bool find_terraces,
    double min_terrace_height_above_floodplain_metres, double max_terrace_height_metres)
{
    require_same_grid(admitted, hand, "admitted and hand must share a grid");
    require_same_grid(admitted, channel_mask, "admitted and channel_mask must share a grid");
    const Raster* gradient_raster = nullptr;
    if (!gradient.is_none()) {
        gradient_raster = &nb::cast<const Raster&>(gradient);
        require_same_grid(admitted, *gradient_raster, "gradient must share the admitted grid");
    } else if (max_floodplain_gradient > 0.0) {
        throw std::invalid_argument("max_floodplain_gradient > 0 needs a gradient raster");
    }
    const Raster* terrace_raster = nullptr;
    if (!terrace_field.is_none()) {
        terrace_raster = &nb::cast<const Raster&>(terrace_field);
        require_same_grid(admitted, *terrace_raster, "terrace_field must share the admitted grid");
    }
    lsdtt::floodplain_threshold_request request;
    request.max_floodplain_hand_metres = max_floodplain_hand_metres;
    request.max_floodplain_gradient = max_floodplain_gradient;
    request.seed_channel_pixels = seed_channel_pixels;
    request.require_channel_contact = require_channel_contact;
    request.min_component_pixels = min_band_component_pixels;
    request.fill_holes = fill_floodplain_holes;
    request.find_terraces = find_terraces;
    request.min_terrace_height_above_floodplain_metres = min_terrace_height_above_floodplain_metres;
    request.max_terrace_height_metres = max_terrace_height_metres;
    // Without a gradient test the engine reads no gradient; pass HAND.
    return FloodplainThresholdResult(lsdtt::classify_floodplain_terraces_by_thresholds(
        admitted.native(), hand.native(),
        gradient_raster ? gradient_raster->native() : hand.native(),
        nonzero_cells(channel_mask),
        terrace_raster ? &terrace_raster->native() : nullptr, request));
}

nb::dict cluster_surface_levels(nb::object values, double required_gap_metres, int n_null,
                                std::uint64_t seed)
{
    const lsdtt::level_clustering_result result = lsdtt::test_level_clustering(
        sample_values(values), required_gap_metres, n_null, seed);
    nb::list rows;
    for (const auto& c : result.clusters) {
        nb::dict row;
        row["lower_metres"] = c.lo;
        row["upper_metres"] = c.hi;
        row["centre_metres"] = c.centre;
        row["n_values"] = c.n;
        rows.append(std::move(row));
    }
    nb::dict out;
    out["clusters"] = nb::module_::import_("pandas").attr("DataFrame")(rows);
    out["n_values"] = result.n_values;
    out["empty_fraction"] = result.empty_fraction;
    out["null_mean_empty_fraction"] = result.null_mean_empty_fraction;
    out["p_value"] = result.p_value;
    out["n_null"] = result.n_null;
    return out;
}

std::vector<lsdtt::method_class_map> class_maps(nb::dict surface_classes,
                                                std::vector<std::uint8_t>& valid,
                                                lsdtt::lsdtt_raster_meta* grid = nullptr)
{
    std::vector<lsdtt::method_class_map> maps;
    const Raster* first = nullptr;
    for (auto [key, value] : surface_classes) {
        const Raster& raster = nb::cast<const Raster&>(value);
        if (first == nullptr) {
            first = &raster;
            valid.assign(raster.native().get_mask().begin(), raster.native().get_mask().end());
        } else {
            require_same_grid(*first, raster, "every surface-class raster must share one grid");
            for (std::size_t i = 0; i < valid.size(); ++i)
                if (!raster.native().get_mask()[i]) valid[i] = 0;
        }
        maps.push_back(class_map(nb::cast<std::string>(key), raster));
    }
    if (maps.size() < 2) throw std::invalid_argument("surface_classes needs at least two methods");
    if (grid != nullptr) *grid = first->native().get_meta();
    return maps;
}

nb::object compare_surface_classes(nb::dict surface_classes)
{
    std::vector<std::uint8_t> valid;
    const auto maps = class_maps(surface_classes, valid);
    nb::list rows;
    for (const auto& c : lsdtt::compare_surface_classes(maps, valid)) {
        nb::dict row;
        row["method_a"] = maps[c.a].name;
        row["method_b"] = maps[c.b].name;
        row["floodplain_iou"] = c.floodplain_iou;
        row["terrace_iou"] = c.terrace_iou;
        row["classified_agreement"] = c.classified_agreement;
        nb::list confusion;
        for (const auto& line : c.confusion) {
            nb::list l;
            for (const auto v : line) l.append(v);
            confusion.append(std::move(l));
        }
        row["confusion"] = std::move(confusion);
        rows.append(std::move(row));
    }
    return nb::module_::import_("pandas").attr("DataFrame")(rows);
}

Raster surface_class_agreement(nb::dict surface_classes, const std::string& surface_class)
{
    std::uint8_t cls = 0;
    if (surface_class == "floodplain") cls = 1;
    else if (surface_class == "terrace") cls = 2;
    else throw std::invalid_argument("surface_class must be floodplain or terrace");
    std::vector<std::uint8_t> valid;
    lsdtt::lsdtt_raster_meta meta;
    const auto maps = class_maps(surface_classes, valid, &meta);
    const auto bits = lsdtt::surface_class_agreement_bits(maps, valid, cls);
    std::vector<float> data(bits.begin(), bits.end());
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(data), std::move(valid), meta));
}

nb::dict score_surface_classes(const Raster& surface_class,
                               const std::string& mapped_surfaces_fname,
                               const std::string& mapped_surfaces_layer_name,
                               const std::string& mapped_surfaces_class_field_name,
                               const std::string& mapped_surfaces_extent_fname)
{
    // The mapped classes and surveyed extent, read as lsdtt-fluvial-surfaces
    // reads them: polygons by cell centre, the extent by majority coverage.
    lsdtt::gdal_init();
    const auto& grid = surface_class.native();
    const lsdtt::training_set mapped_set = lsdtt::read_training_surfaces(
        grid, mapped_surfaces_fname, mapped_surfaces_layer_name,
        mapped_surfaces_class_field_name);
    const std::size_t n = grid.get_data().size();
    std::vector<std::uint8_t> mapped(n, 0), scored(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
        scored[i] = grid.get_mask()[i] ? 1 : 0;
        const std::int32_t f = mapped_set.polygon_of_pixel[i];
        if (f >= 0)
            mapped[i] = static_cast<std::uint8_t>(
                mapped_set.polygons[static_cast<std::size_t>(f)].cls + 1);
    }
    if (!mapped_surfaces_extent_fname.empty()) {
        const lsdtt::lsdtt_raster_cpu extent = lsdtt::rasterise_polygons(
            grid, mapped_surfaces_extent_fname, lsdtt::polygon_rasterise_rule::majority_coverage);
        for (std::size_t i = 0; i < n; ++i)
            if (!(extent.get_data()[i] > 0.5f)) scored[i] = 0;
    }
    const lsdtt::mapped_surface_score score =
        lsdtt::score_against_mapped(class_map("method", surface_class), mapped, scored);
    nb::dict out;
    out["n_scored"] = score.n_scored;
    out["floodplain_precision"] = score.precision[1];
    out["floodplain_recall"] = score.recall[1];
    out["floodplain_iou"] = score.iou[1];
    out["terrace_precision"] = score.precision[2];
    out["terrace_recall"] = score.recall[2];
    out["terrace_iou"] = score.iou[2];
    nb::list confusion;
    for (const auto& line : score.confusion) {
        nb::list l;
        for (const auto v : line) l.append(v);
        confusion.append(std::move(l));
    }
    out["confusion"] = std::move(confusion);
    return out;
}

Raster relief_above_nearest_baseline(const Raster& dem, const Raster& baseline_mask,
                                     double swath_half_width_metres, int n_threads)
{
    require_same_grid(dem, baseline_mask, "baseline_mask must share the DEM grid");
    const auto cells = nonzero_cells(baseline_mask);
    std::vector<std::int64_t> baseline;
    for (std::size_t i = 0; i < cells.size(); ++i)
        if (cells[i]) baseline.push_back(static_cast<std::int64_t>(i));
    if (baseline.empty()) throw std::invalid_argument("baseline_mask has no cells > 0.5");
    return Raster(lsdtt::relief_above_nearest_baseline(
        dem.native(), baseline, swath_half_width_metres, n_threads).relief);
}

} // namespace

void bind_fluvial_surfaces(nb::module_& m)
{
    nb::class_<FloodplainThresholdResult>(m, "FloodplainThresholdResult", R"doc(
Floodplain and terrace masks from :func:`classify_floodplain_terraces_by_thresholds`.

``band_index`` is 0 on floodplain, 1 on terrace and -1 elsewhere, the
encoding of the ladder's band index. The counts report what the patch
filters removed.
)doc")
        .def_prop_ro("floodplain_mask", [](const FloodplainThresholdResult& r) { return Raster(r.native().floodplain_mask); })
        .def_prop_ro("terrace_mask", [](const FloodplainThresholdResult& r) { return Raster(r.native().terrace_mask); })
        .def_prop_ro("band_index", [](const FloodplainThresholdResult& r) { return Raster(r.native().band_index); })
        .def_prop_ro("n_floodplain_candidates", [](const FloodplainThresholdResult& r) { return r.native().n_floodplain_candidates; })
        .def_prop_ro("n_floodplain_components_kept", [](const FloodplainThresholdResult& r) { return r.native().n_floodplain_components_kept; })
        .def_prop_ro("n_floodplain_dropped_small", [](const FloodplainThresholdResult& r) { return r.native().n_floodplain_dropped_small; })
        .def_prop_ro("n_floodplain_dropped_no_channel", [](const FloodplainThresholdResult& r) { return r.native().n_floodplain_dropped_no_channel; })
        .def_prop_ro("n_floodplain_pixels", [](const FloodplainThresholdResult& r) { return r.native().n_floodplain_pixels; })
        .def_prop_ro("n_terrace_components", [](const FloodplainThresholdResult& r) { return r.native().n_terrace_components; })
        .def_prop_ro("n_terrace_pixels", [](const FloodplainThresholdResult& r) { return r.native().n_terrace_pixels; });

    m.def("select_threshold", &select_threshold, nb::arg("values"),
        nb::arg("threshold_algorithm_choice") = "qq_gaussian", nb::arg("fixed_value") = 0.0,
        nb::arg("qq_lower_percentile") = 25.0, nb::arg("qq_upper_percentile") = 75.0,
        nb::arg("qq_threshold_fraction") = 0.01, R"doc(
Pick a threshold from a sample, as lsdtt-fluvial-surfaces picks its floodplain bounds.

The routine behind ``floodplain_threshold_algorithm_choice`` (and the
``qq_upper_departure`` floodplain second pass).

Parameters
----------
values : Raster or array-like of float
    The sample; for a Raster, its valid cells. Non-finite values are dropped.
threshold_algorithm_choice : {"fixed", "qq_gaussian", "cdf_divergence", "qq_upper_departure"}, default="qq_gaussian"
    ``qq_gaussian`` is Clubb et al. (2017): the lowest value at which the
    quantile-quantile plot rejoins a Gaussian drawn through two percentiles
    and stays on it. ``cdf_divergence`` fits a family of distributions and
    takes the value where the empirical CDF most exceeds the best fit.
    ``qq_upper_departure`` mirrors ``qq_gaussian`` for the upper tail.
    ``fixed`` returns ``fixed_value``.
fixed_value : float, default=0.0
    The threshold for ``fixed``.
qq_lower_percentile, qq_upper_percentile : float, default=25.0, 75.0
    Percentiles the Gaussian reference line is drawn through.
qq_threshold_fraction : float, default=0.01
    Fraction of the data range within which a quantile counts as on the line.

Returns
-------
dict
    ``threshold``, ``found`` (False when no threshold could be chosen, and
    ``threshold`` is then meaningless), ``n_values``, and for
    ``cdf_divergence`` ``best_family``, ``threshold_percentile`` and
    ``max_cdf_divergence``.
)doc");

    m.def("classify_floodplain_terraces_by_thresholds",
        &classify_floodplain_terraces_by_thresholds,
        nb::arg("admitted"), nb::arg("hand"), nb::arg("channel_mask"),
        nb::arg("gradient") = nb::none(), nb::arg("terrace_field") = nb::none(),
        nb::arg("max_floodplain_hand_metres") = 0.0, nb::arg("max_floodplain_gradient") = 0.0,
        nb::arg("seed_channel_pixels") = true, nb::arg("require_channel_contact") = true,
        nb::arg("min_band_component_pixels") = 0, nb::arg("fill_floodplain_holes") = false,
        nb::arg("find_terraces") = false,
        nb::arg("min_terrace_height_above_floodplain_metres") = 0.0,
        nb::arg("max_terrace_height_metres") = 0.0, R"doc(
Floodplains and terraces from a HAND bound and a gradient bound (Clubb et al. 2017).

The classifier of lsdtt-fluvial-surfaces' ``clubb2017`` method and of
lsdtt-hillslope-channel-coupling's floodplain rule. Floodplain: admitted
cells with HAND below ``max_floodplain_hand_metres`` and gradient below
``max_floodplain_gradient`` (when > 0), plus the channel cells when
``seed_channel_pixels``; 8-connected components smaller than
``min_band_component_pixels`` are dropped, and with
``require_channel_contact`` only components containing a channel cell are
kept. Terraces (``find_terraces``): admitted, non-floodplain cells of
gentle gradient whose ``terrace_field`` lies more than
``min_terrace_height_above_floodplain_metres`` above the floodplain bound
and at most ``max_terrace_height_metres`` (when > 0). Pick the bounds from
the data with :func:`select_threshold`.

Parameters
----------
admitted : Raster
    Cells > 0.5 may be floodplain or terrace; its mask is the output domain.
hand : Raster
    Height above the channel (m), co-registered.
channel_mask : Raster
    Channel cells (> 0.5), for seeding and the contact test.
gradient : Raster or None, default=None
    Gradient (m/m); required when ``max_floodplain_gradient`` > 0.
terrace_field : Raster or None, default=None
    The terrace reference (HAND, or :func:`relief_above_nearest_baseline`);
    required with ``find_terraces``.
max_floodplain_hand_metres : float, default=0.0
    Floodplain HAND bound (m).
max_floodplain_gradient : float, default=0.0
    Floodplain gradient bound (m/m); 0 applies none.
seed_channel_pixels, require_channel_contact : bool, default=True
    Add the channel cells to the floodplain; keep only components that
    contain one.
min_band_component_pixels : int, default=0
    Smallest component kept (pixels); 0 or 1 keeps all.
fill_floodplain_holes : bool, default=False
    Fill holes in the floodplain components.
find_terraces : bool, default=False
    Also classify terraces.
min_terrace_height_above_floodplain_metres : float, default=0.0
    Terraces start this far (m, >= 0) above the floodplain bound.
max_terrace_height_metres : float, default=0.0
    Terrace ceiling (m) of ``terrace_field``; 0 means none.

Returns
-------
FloodplainThresholdResult
)doc");

    m.def("cluster_surface_levels", &cluster_surface_levels, nb::arg("values"),
        nb::arg("required_gap_metres"), nb::arg("n_null") = 999, nb::arg("seed") = 1, R"doc(
Group surface heights into levels and test whether they cluster.

The test behind lsdtt-fluvial-surfaces' ``write_patch_level_clusters_csv``:
sorted heights are split wherever consecutive values are at least
``required_gap_metres`` apart. The statistic is the empty fraction, the
share of the height range taken by those gaps, compared with ``n_null``
samples of as many heights drawn uniformly over the same range:
``p = (1 + #{null >= observed}) / (1 + n_null)``. A small p says the
heights sit on levels more than random heights would.

Parameters
----------
values : Raster or array-like of float
    Heights (m), for example the ``median_field`` of
    ``FlatSurfaceClassification.patches``.
required_gap_metres : float
    Gap (m, > 0) that separates two levels, typically kappa x E_v.
n_null : int, default=999
    Number of uniform null samples.
seed : int, default=1
    Random seed of the null samples.

Returns
-------
dict
    ``clusters`` (DataFrame: lower_metres, upper_metres, centre_metres (the
    median), n_values, ascending), ``n_values``, ``empty_fraction``,
    ``null_mean_empty_fraction``, ``p_value`` (1 with fewer than 3 values)
    and ``n_null``.
)doc");

    m.def("compare_surface_classes", &compare_surface_classes, nb::arg("surface_classes"), R"doc(
Agreement between the surface classes of several methods, pair by pair.

The comparison lsdtt-fluvial-surfaces writes when it runs several methods.

Parameters
----------
surface_classes : dict of str to Raster
    Method name to a surface-class raster (0 none, 1 floodplain, 2
    terrace; ``FlatSurfaceClassification.surface_class``), all on one grid.
    Cells invalid in any raster are not counted.

Returns
-------
pandas.DataFrame
    One row per pair (method_a before method_b, in dict order):
    ``floodplain_iou``, ``terrace_iou`` (NaN when neither method has the
    class), ``classified_agreement`` (of the cells either method classed,
    the fraction both gave the same class) and ``confusion`` (3 x 3 counts,
    ``[class_a][class_b]``).
)doc");

    m.def("surface_class_agreement", &surface_class_agreement, nb::arg("surface_classes"),
        nb::arg("surface_class") = "floodplain", R"doc(
Which methods gave each cell a class, as a bit field.

Parameters
----------
surface_classes : dict of str to Raster
    As in :func:`compare_surface_classes` (at most 8 methods).
surface_class : {"floodplain", "terrace"}, default="floodplain"
    The class to count.

Returns
-------
Raster
    Bit k set where method k (dict order) gave ``surface_class``; 0 elsewhere.
)doc");

    m.def("score_surface_classes", &score_surface_classes, nb::arg("surface_class"),
        nb::arg("mapped_surfaces_fname"), nb::arg("mapped_surfaces_layer_name") = "",
        nb::arg("mapped_surfaces_class_field_name") = "class",
        nb::arg("mapped_surfaces_extent_fname") = "", R"doc(
Score a surface classification against mapped floodplains and terraces.

As lsdtt-fluvial-surfaces' ``mapped_surfaces_*`` parameters: the polygons
are burned by cell centre, and only cells inside the surveyed extent (and
valid in ``surface_class``) count.

Parameters
----------
surface_class : Raster
    0 none, 1 floodplain, 2 terrace.
mapped_surfaces_fname : str
    Vector file of mapped polygons classed ``floodplain``, ``terrace`` or
    ``terraceN``.
mapped_surfaces_layer_name : str, default=""
    Its layer; empty reads the first.
mapped_surfaces_class_field_name : str, default="class"
    Its class field.
mapped_surfaces_extent_fname : str, default=""
    Polygons of the surveyed extent (majority coverage); empty counts the
    whole raster as surveyed.

Returns
-------
dict
    ``n_scored``; ``floodplain_precision``, ``floodplain_recall``,
    ``floodplain_iou`` and the same for terrace (NaN on an empty
    denominator); ``confusion`` (3 x 3 counts, ``[mapped][method]``).
)doc");

    m.def("relief_above_nearest_baseline", &relief_above_nearest_baseline, nb::arg("dem"),
        nb::arg("baseline_mask"), nb::arg("swath_half_width_metres") = 0.0,
        nb::arg("n_threads") = 0, R"doc(
Height of each cell above its nearest baseline cell (swath baseline relief).

lsdtt-fluvial-surfaces' ``terrace_reference_field_algorithm_choice
swath_baseline`` with the channel as the baseline.

Parameters
----------
dem : Raster
    Elevation (m).
baseline_mask : Raster
    Baseline cells (> 0.5), for example the channel mask.
swath_half_width_metres : float, default=0.0
    Largest distance (m) from the baseline at which relief is computed; 0
    means no limit.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
Raster
    Elevation minus the elevation of the nearest baseline cell (m).
)doc");
}

} // namespace lsdtt3_python
