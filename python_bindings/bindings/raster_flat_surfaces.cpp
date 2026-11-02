#include <cmath>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <core/lsdtt_common_parameters.hpp>
#include <core/lsdtt_gdal_init.hpp>
#include <cpu/lsdtt_flat_surface_channel_classify_cpu.hpp>
#include <cpu/lsdtt_flat_surface_finder.hpp>
#include <cpu/lsdtt_raster_surface_metrics_cpu.hpp>
#include <cpu/lsdtt_surface_training_cpu.hpp>
#include <cpu/lsdtt_valley_cross_sections_cpu.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

lsdtt::iqr_fit_method parse_iqr_method(const std::string& method)
{
    if (method == "fftw_highpass") return lsdtt::iqr_fit_method::fftw_highpass;
    if (method == "wavelet_quadratic") return lsdtt::iqr_fit_method::wavelet_quadratic;
    if (method == "moving_window_mean") return lsdtt::iqr_fit_method::moving_window_mean;
    throw std::invalid_argument(
        "iqr_method must be fftw_highpass, wavelet_quadratic, or moving_window_mean"
    );
}

lsdtt::window_shape parse_window_shape(const std::string& shape)
{
    if (shape == "square") return lsdtt::window_shape::square;
    if (shape == "circular") return lsdtt::window_shape::circular;
    throw std::invalid_argument("window_shape must be square or circular");
}

// The ladder choices, with lsdtt-fluvial-surfaces' names and messages.
lsdtt::ladder_break_choice parse_ladder_break(const std::string& name)
{
    if (name == "antimode") return lsdtt::ladder_break_choice::antimode;
    if (name == "jenks") return lsdtt::ladder_break_choice::jenks;
    if (name == "gap") return lsdtt::ladder_break_choice::gap;
    throw std::invalid_argument("unknown ladder_break_algorithm_choice '" + name +
                                "'. Valid: antimode, jenks, gap.");
}

lsdtt::gap_break_position parse_gap_break_position(const std::string& name)
{
    if (name == "lower") return lsdtt::gap_break_position::lower;
    if (name == "midpoint") return lsdtt::gap_break_position::midpoint;
    if (name == "upper") return lsdtt::gap_break_position::upper;
    throw std::invalid_argument("unknown ladder_gap_break_position '" + name +
                                "'. Valid: lower, midpoint, upper.");
}

lsdtt::gap_rank_choice parse_gap_rank(const std::string& name)
{
    if (name == "lowest") return lsdtt::gap_rank_choice::lowest;
    if (name == "width") return lsdtt::gap_rank_choice::width;
    throw std::invalid_argument("unknown ladder_gap_rank '" + name +
                                "'. Valid: lowest, width.");
}

const char* ladder_status_name(lsdtt::ladder_status status)
{
    switch (status) {
        case lsdtt::ladder_status::resolved:   return "resolved";
        case lsdtt::ladder_status::unresolved: return "unresolved";
        case lsdtt::ladder_status::forced:     return "forced";
        case lsdtt::ladder_status::degenerate: return "degenerate";
        case lsdtt::ladder_status::supplied:   return "supplied";
    }
    return "unknown";
}

void require_matching_grids(const Raster& first, const Raster& second,
                            const char* message)
{
    const auto& a = first.native().get_meta();
    const auto& b = second.native().get_meta();
    if (a.n_rows != b.n_rows || a.n_cols != b.n_cols || a.x_min != b.x_min ||
        a.y_min != b.y_min || a.cell_size != b.cell_size || a.crs_wkt != b.crs_wkt)
    {
        throw std::invalid_argument(message);
    }
}

class FlatSurfaceResult {
public:
    explicit FlatSurfaceResult(lsdtt::flat_surface_result result)
        : result_(std::move(result)) {}
    const lsdtt::flat_surface_result& native() const noexcept { return result_; }
private:
    lsdtt::flat_surface_result result_;
};

class FlatSurfaceClassification {
public:
    explicit FlatSurfaceClassification(lsdtt::flat_surface_classify_result result)
        : result_(std::move(result)) {}
    const lsdtt::flat_surface_classify_result& native() const noexcept { return result_; }
private:
    lsdtt::flat_surface_classify_result result_;
};

FlatSurfaceResult find_flat_surfaces(
    const Raster& dem, bool use_madk2, double madk2_window_radius_metres,
    double madk2_threshold_percentile, bool use_residual_iqr,
    const std::string& iqr_method, double iqr_fit_lambda_metres,
    double iqr_detrend_radius_metres, double iqr_window_radius_metres,
    const std::string& iqr_window_shape, double iqr_threshold_percentile,
    bool use_nvd, double nvd_polyfit_lambda, double nvd_window_radius_metres,
    const std::string& nvd_window_shape, double nvd_threshold_percentile,
    bool require_all_metrics_pass, int directional_close_line_length_pixels,
    int directional_close_n_angles, int min_component_pixels,
    double relief_window_radius_metres, double max_relief_metres,
    bool remove_risers, int n_threads)
{
    lsdtt::flat_surface_request request;
    request.use_madk2 = use_madk2;
    request.madk2_window_radius_metres = madk2_window_radius_metres;
    request.madk2_threshold_percentile = madk2_threshold_percentile;
    request.use_residual_iqr = use_residual_iqr;
    request.iqr_method = parse_iqr_method(iqr_method);
    request.iqr_fit_lambda_metres = iqr_fit_lambda_metres;
    request.iqr_detrend_radius_metres = iqr_detrend_radius_metres;
    request.iqr_window_radius_metres = iqr_window_radius_metres;
    request.iqr_window_shape = parse_window_shape(iqr_window_shape);
    request.iqr_threshold_percentile = iqr_threshold_percentile;
    request.use_nvd = use_nvd;
    request.nvd_polyfit_lambda = nvd_polyfit_lambda;
    request.nvd_window_radius_metres = nvd_window_radius_metres;
    request.nvd_window_shape = parse_window_shape(nvd_window_shape);
    request.nvd_threshold_percentile = nvd_threshold_percentile;
    request.require_all_metrics_pass = require_all_metrics_pass;
    request.directional_close_line_length_pixels = directional_close_line_length_pixels;
    request.directional_close_n_angles = directional_close_n_angles;
    request.min_component_pixels = min_component_pixels;
    // A radius r cells is a (2r + 1)-cell side; 0 is 2 cells, as in
    // lsdtt-flat-surface-finder.
    request.relief_window_pixels = 2 * lsdtt::relief_window_radius_pixels(
        relief_window_radius_metres, dem.native().get_meta().cell_size, 2) + 1;
    request.relief_threshold_metres = max_relief_metres;
    request.remove_risers = remove_risers;
    request.n_threads = n_threads;
    return FlatSurfaceResult(lsdtt::find_flat_surfaces(dem.native(), request));
}

std::vector<std::int64_t> donor_vector(nb::object donor_flat, std::size_t n_pixels)
{
    if (donor_flat.is_none()) return {};
    nb::object array = nb::module_::import_("numpy").attr("asarray")(
        donor_flat, nb::arg("dtype") = "int64"
    );
    if (nb::cast<int>(array.attr("size")) != static_cast<int>(n_pixels))
        throw std::invalid_argument("donor_flat must have one value per raster pixel");
    return nb::cast<std::vector<std::int64_t>>(array.attr("ravel")());
}

FlatSurfaceClassification classify_flat_surfaces(
    const Raster& labelled_flat_mask, const Raster& height_above_channel,
    const Raster& channel_source_ids, double max_floodplain_hand_metres,
    double min_terrace_hand_metres, double channel_adjacency_radius_metres,
    double valley_corridor_half_width_metres, bool enable_leak_split,
    bool enable_hand_banding, nb::object donor_flat,
    nb::object channel_flow_distance, int n_threads)
{
    require_matching_grids(labelled_flat_mask, height_above_channel,
        "labelled_flat_mask and height_above_channel must share a grid");
    require_matching_grids(labelled_flat_mask, channel_source_ids,
        "labelled_flat_mask and channel_source_ids must share a grid");

    const auto& source = channel_source_ids.native();
    const std::size_t n_pixels = source.get_data().size();
    std::vector<std::uint8_t> channel_mask(n_pixels, 0);
    std::vector<std::int32_t> source_ids(n_pixels, -1);
    for (std::size_t i = 0; i < n_pixels; ++i) {
        if (source.get_mask()[i] == 0 || source.get_data()[i] < 0.0f) continue;
        channel_mask[i] = 1;
        source_ids[i] = static_cast<std::int32_t>(source.get_data()[i]);
    }

    const std::vector<std::int64_t> donors = donor_vector(donor_flat, n_pixels);
    const Raster* flow_distance = nullptr;
    if (!channel_flow_distance.is_none()) {
        const Raster& raster = nb::cast<const Raster&>(channel_flow_distance);
        require_matching_grids(labelled_flat_mask, raster,
            "channel_flow_distance must share the labelled-flat grid");
        flow_distance = &raster;
    }
    if (enable_hand_banding && (donors.empty() || flow_distance == nullptr)) {
        throw std::invalid_argument(
            "reach-relative HAND banding requires donor_flat and channel_flow_distance"
        );
    }

    lsdtt::flat_surface_classify_request request;
    request.floodplain_max_height_metres = max_floodplain_hand_metres;
    request.terrace_min_height_metres = min_terrace_hand_metres;
    request.channel_adjacency_radius_metres = channel_adjacency_radius_metres;
    request.terrace_max_distance_to_channel_metres = valley_corridor_half_width_metres;
    request.enable_leak_split = enable_leak_split;
    request.enable_hand_banding = enable_hand_banding;
    request.n_threads = n_threads;
    return FlatSurfaceClassification(lsdtt::classify_flat_surfaces_by_channel(
        labelled_flat_mask.native(), height_above_channel.native(), channel_mask,
        source_ids, request, donors.empty() ? nullptr : &donors,
        flow_distance ? &flow_distance->native() : nullptr
    ));
}

FlatSurfaceClassification classify_flat_surfaces_with_scope(
    const Raster& labelled_flat_mask, const Raster& dem, const ChannelScope& scope,
    double max_floodplain_hand_metres, double min_terrace_hand_metres,
    double channel_adjacency_radius_metres, double valley_corridor_half_width_metres,
    bool enable_leak_split, bool enable_hand_banding, int n_threads,
    const std::string& ladder_break_algorithm_choice, int ladder_level_count,
    int ladder_max_level_count, double vertical_accuracy_metres,
    const std::string& vertical_accuracy_algorithm_choice,
    double vertical_accuracy_gentle_metres, double vertical_accuracy_steep_metres,
    double vertical_accuracy_break_gradient,
    double vertical_accuracy_gradient_window_radius_metres,
    double ladder_resolvability_kappa, int ladder_k_bootstrap, int ladder_k_bootstrap_seed,
    double ladder_reach_window_metres, int ladder_min_window_pixels, int ladder_hist_bins,
    double ladder_min_field_metres, double ladder_max_valley_fraction,
    const std::string& ladder_gap_break_position, const std::string& ladder_gap_rank,
    double ladder_gap_trough_max_ratio, double ladder_gap_empty_fraction,
    double ladder_gap_min_width_metres, int ladder_gap_hist_bins,
    double ladder_max_separation_ratio, double ladder_min_break_metres,
    double ladder_max_break_metres, int min_band_component_pixels,
    const std::string& ladder_unresolved_window_algorithm_choice,
    bool clip_corridor_per_pixel, double max_band_spread_metres,
    const std::string& segmentation_algorithm_choice, nb::object riser_mask,
    bool break_flats_with_channels,
    const std::map<std::int32_t, std::vector<double>>& supplied_levels,
    const std::string& training_surfaces_fname,
    const std::string& training_surfaces_layer_name,
    const std::string& training_surfaces_class_field_name,
    double training_spread_factor, double training_max_relative_spread,
    int min_training_cells_per_class)
{
    require_matching_grids(labelled_flat_mask, dem,
        "labelled_flat_mask and dem must share a grid");
    const auto& flow = scope.flow();
    const auto& meta = flow.meta();
    const auto& dem_meta = dem.native().get_meta();
    if (meta.n_rows != dem_meta.n_rows || meta.n_cols != dem_meta.n_cols ||
        meta.x_min != dem_meta.x_min || meta.y_min != dem_meta.y_min ||
        meta.cell_size != dem_meta.cell_size || meta.crs_wkt != dem_meta.crs_wkt)
        throw std::invalid_argument("dem must share the ChannelScope flow grid");
    if (segmentation_algorithm_choice != "within_bands" &&
        segmentation_algorithm_choice != "whole_mask")
        throw std::invalid_argument("unknown segmentation_algorithm_choice '" +
            segmentation_algorithm_choice + "'. Valid: within_bands, whole_mask.");
    if (vertical_accuracy_algorithm_choice != "fixed" &&
        vertical_accuracy_algorithm_choice != "slope_dependent")
        throw std::invalid_argument("unknown vertical_accuracy_algorithm_choice '" +
            vertical_accuracy_algorithm_choice + "'. Valid: fixed, slope_dependent.");
    if (ladder_unresolved_window_algorithm_choice != "nearest_resolved" &&
        ladder_unresolved_window_algorithm_choice != "single_band")
        throw std::invalid_argument("unknown ladder_unresolved_window_algorithm_choice '" +
            ladder_unresolved_window_algorithm_choice +
            "'. Valid: nearest_resolved, single_band.");

    const auto& scope_data = scope.native();
    auto hand = flow.get_height_above_channel_mask_donor_raster(
        dem.native().get_data(), dem.native().get_mask(), scope_data.channel_mask);
    lsdtt::lsdtt_raster_cpu hand_raster(
        std::move(hand.height.data), std::move(hand.height.mask), hand.height.meta);

    lsdtt::flat_surface_classify_request request;
    request.floodplain_max_height_metres = max_floodplain_hand_metres;
    request.terrace_min_height_metres = min_terrace_hand_metres;
    request.channel_adjacency_radius_metres = channel_adjacency_radius_metres;
    request.terrace_max_distance_to_channel_metres = valley_corridor_half_width_metres;
    request.enable_leak_split = enable_leak_split;
    request.enable_hand_banding = enable_hand_banding;
    request.n_threads = n_threads;
    request.terrace_corridor_per_pixel = clip_corridor_per_pixel;
    request.terrace_max_hand_spread_metres = max_band_spread_metres;
    request.label_within_bands = true;
    request.ladder_break_algorithm = parse_ladder_break(ladder_break_algorithm_choice);
    request.ladder_level_count = ladder_level_count;
    request.ladder_max_level_count = ladder_max_level_count;
    request.vertical_accuracy_metres = vertical_accuracy_metres;
    request.k_bootstrap = ladder_k_bootstrap;
    request.k_bootstrap_seed = ladder_k_bootstrap_seed;
    const bool ev_slope_dependent = vertical_accuracy_algorithm_choice == "slope_dependent";
    request.vertical_accuracy_mode = ev_slope_dependent
        ? lsdtt::vertical_accuracy_choice::slope_dependent
        : lsdtt::vertical_accuracy_choice::fixed;
    request.vertical_accuracy_gentle_metres = vertical_accuracy_gentle_metres;
    request.vertical_accuracy_steep_metres = vertical_accuracy_steep_metres;
    request.vertical_accuracy_break_gradient = vertical_accuracy_break_gradient;
    request.resolvability_kappa = ladder_resolvability_kappa;
    request.band_reach_window_metres = ladder_reach_window_metres;
    request.band_min_window_pixels = ladder_min_window_pixels;
    request.band_hist_bins = ladder_hist_bins;
    request.band_min_hand_metres = ladder_min_field_metres;
    request.band_max_valley_fraction = ladder_max_valley_fraction;
    request.ladder_gap_break_at = parse_gap_break_position(ladder_gap_break_position);
    request.ladder_gap_rank = parse_gap_rank(ladder_gap_rank);
    request.ladder_gap_trough_max_ratio = ladder_gap_trough_max_ratio;
    request.ladder_gap_empty_fraction = ladder_gap_empty_fraction;
    request.ladder_gap_min_width_metres = ladder_gap_min_width_metres;
    request.ladder_gap_hist_bins = ladder_gap_hist_bins;
    request.ladder_max_separation_ratio = ladder_max_separation_ratio;
    request.ladder_min_break_metres = ladder_min_break_metres;
    request.ladder_max_break_metres = ladder_max_break_metres;
    request.band_min_component_pixels = min_band_component_pixels;
    request.ladder_unresolved_borrow_breaks =
        ladder_unresolved_window_algorithm_choice == "nearest_resolved";
    if (!supplied_levels.empty())
        request.ladder_supplied_breaks = lsdtt::breaks_from_levels(supplied_levels);

    // Risers as barriers for whole_mask segmentation (cells > 0.5).
    lsdtt::lsdtt_raster_cpu riser_barrier;
    request.segment_whole_mask = segmentation_algorithm_choice == "whole_mask";
    if (!riser_mask.is_none()) {
        if (!request.segment_whole_mask)
            throw std::invalid_argument(
                "riser_mask is used only with segmentation_algorithm_choice whole_mask");
        const Raster& risers = nb::cast<const Raster&>(riser_mask);
        require_matching_grids(dem, risers, "riser_mask must share the DEM grid");
        const std::size_t n = risers.native().get_data().size();
        std::vector<float> rd(n, 0.0f);
        for (std::size_t i = 0; i < n; ++i)
            if (risers.native().get_mask()[i] && risers.native().get_data()[i] > 0.5f) rd[i] = 1.0f;
        riser_barrier = lsdtt::lsdtt_raster_cpu(std::move(rd), std::vector<std::uint8_t>(n, 1), dem_meta);
        request.whole_mask_barrier = &riser_barrier;
    }

    // Training polygons, read against HAND as lsdtt-fluvial-surfaces does.
    lsdtt::training_set training;
    if (!training_surfaces_fname.empty()) {
        lsdtt::gdal_init();
        training = lsdtt::read_training_surfaces(
            hand_raster, training_surfaces_fname, training_surfaces_layer_name,
            training_surfaces_class_field_name);
        request.ladder_training = &training;
    }
    request.ladder_training_spread_factor = training_spread_factor;
    request.ladder_training_relative_spread = training_max_relative_spread;
    request.ladder_training_min_cells_per_class = min_training_cells_per_class;

    // break_flats_with_channels: channel cells leave the admitted set.
    lsdtt::lsdtt_raster_cpu labels = labelled_flat_mask.native();
    if (break_flats_with_channels) {
        std::vector<float>& ld = labels.get_data();
        for (std::size_t i = 0; i < ld.size(); ++i)
            if (scope_data.channel_mask[i] && ld[i] > 0.5f) ld[i] = 0.0f;
    }

    // Slope-dependent E_v reads the window gradient, measured as the CLI
    // does: a polyfit over twice the radius (0: five cells).
    lsdtt::lsdtt_raster_cpu ev_gradient;
    if (ev_slope_dependent) {
        double lambda = 2.0 * vertical_accuracy_gradient_window_radius_metres;
        if (!(lambda > 0.0)) lambda = 5.0 * dem_meta.cell_size;
        const lsdtt::poly_fit_coeffs fit = lsdtt::direct_polyfit(dem.native(), lambda, n_threads);
        lsdtt::polyfit_metric_request mreq;
        mreq.metrics = {lsdtt::polyfit_metric::slope};
        mreq.n_threads = n_threads;
        lsdtt::polyfit_metric_result mres = lsdtt::compute_polyfit_metrics(fit, mreq);
        ev_gradient = mres.at(lsdtt::polyfit_metric::slope);
    }

    auto flow_distance_result = flow.get_distance_from_outlet_raster();
    lsdtt::lsdtt_raster_cpu flow_distance(
        std::move(flow_distance_result.data), std::move(flow_distance_result.mask),
        flow_distance_result.meta);
    // The ladder windows along the channel through the HAND donors.
    return FlatSurfaceClassification(lsdtt::classify_flat_surfaces_by_channel(
        labels, hand_raster, scope_data.channel_mask,
        scope_data.channel_source_id, request, &hand.donor_flat, &flow_distance,
        ev_slope_dependent ? &ev_gradient : nullptr));
}

nb::object class_table(const FlatSurfaceClassification& classification)
{
    nb::list rows;
    for (const auto& record : classification.native().labels) {
        nb::dict row;
        row["label"] = record.label;
        row["surface_class"] = record.cls == lsdtt::surface_class::floodplain
            ? "floodplain" : record.cls == lsdtt::surface_class::terrace ? "terrace" : "none";
        row["source_id"] = record.source_id;
        row["n_pixels"] = record.n_pixels;
        row["median_height"] = record.median_height;
        row["mean_height"] = record.mean_height;
        row["min_distance_to_channel_metres"] = record.min_distance_to_channel_metres;
        row["adjacent_to_channel"] = record.adjacent_to_channel;
        rows.append(std::move(row));
    }
    return nb::module_::import_("pandas").attr("DataFrame")(rows);
}

nb::object ladder_window_table(const FlatSurfaceClassification& classification)
{
    nb::list rows;
    for (const auto& w : classification.native().ladder_windows) {
        nb::dict row;
        row["source_id"] = w.source_id;
        row["window_centre_metres"] = w.window_centre_m;
        row["n_pixels"] = w.n_pixels;
        row["k_selected"] = w.k_selected;
        row["status"] = ladder_status_name(w.status);
        row["used_global_breaks"] = w.used_global_breaks;
        row["breaks_metres"] = w.breaks_m;
        row["band_centres_metres"] = w.band_centres_m;
        row["required_gap_metres"] = w.required_gap_m;
        row["vertical_accuracy_metres"] = w.vertical_accuracy_m;
        row["median_gradient"] = w.median_gradient;
        row["n_trained_breaks"] = w.n_trained_breaks;
        row["training_source"] = w.training_source;
        rows.append(std::move(row));
    }
    return nb::module_::import_("pandas").attr("DataFrame")(rows);
}

nb::object patch_table(const FlatSurfaceClassification& classification)
{
    nb::list rows;
    for (const auto& p : classification.native().patches) {
        nb::dict row;
        row["label"] = p.label;
        row["band"] = p.band;
        row["n_pixels"] = p.n_pixels;
        row["area_m2"] = p.area_m2;
        row["median_field"] = p.median_field;
        rows.append(std::move(row));
    }
    return nb::module_::import_("pandas").attr("DataFrame")(rows);
}

nb::object training_polygon_table(const FlatSurfaceClassification& classification)
{
    nb::list rows;
    for (const auto& t : classification.native().ladder_training_polygons) {
        nb::dict row;
        row["label"] = t.label;
        row["surface_class"] = t.cls == 0 ? "floodplain" : "terrace";
        row["n_cells"] = t.n_cells;
        row["n_with_field"] = t.n_with_field;
        row["median"] = t.median;
        row["p10"] = t.p10;
        row["p25"] = t.p25;
        row["p75"] = t.p75;
        row["p90"] = t.p90;
        row["flagged"] = t.flagged;
        rows.append(std::move(row));
    }
    return nb::module_::import_("pandas").attr("DataFrame")(rows);
}

} // namespace

void bind_raster_flat_surfaces(nb::module_& m)
{
    nb::class_<FlatSurfaceResult>(m, "FlatSurfaceResult", R"doc(
Results from texture- and relief-based flat-surface detection.
)doc")
        .def_prop_ro("madk2_roughness", [](const FlatSurfaceResult& r) { return Raster(r.native().madk2_roughness); })
        .def_prop_ro("iqr_residuals", [](const FlatSurfaceResult& r) { return Raster(r.native().iqr_residuals); })
        .def_prop_ro("iqr_texture", [](const FlatSurfaceResult& r) { return Raster(r.native().iqr_texture); })
        .def_prop_ro("nvd_s3", [](const FlatSurfaceResult& r) { return Raster(r.native().nvd_s3); })
        .def_prop_ro("madk2_flat_mask", [](const FlatSurfaceResult& r) { return Raster(r.native().madk2_flat_mask); })
        .def_prop_ro("iqr_flat_mask", [](const FlatSurfaceResult& r) { return Raster(r.native().iqr_flat_mask); })
        .def_prop_ro("nvd_flat_mask", [](const FlatSurfaceResult& r) { return Raster(r.native().nvd_flat_mask); })
        .def_prop_ro("relief_flat_mask", [](const FlatSurfaceResult& r) { return Raster(r.native().relief_flat_mask); })
        .def_prop_ro("combined_flat_mask", [](const FlatSurfaceResult& r) { return Raster(r.native().combined_flat_mask); })
        .def_prop_ro("healed_flat_mask", [](const FlatSurfaceResult& r) { return Raster(r.native().healed_flat_mask); })
        .def_prop_ro("labelled_flat_mask", [](const FlatSurfaceResult& r) { return Raster(r.native().labelled_flat_mask); })
        .def_prop_ro("relief", [](const FlatSurfaceResult& r) { return Raster(r.native().relief); });
    nb::class_<FlatSurfaceClassification>(m, "FlatSurfaceClassification", R"doc(
Floodplain and terrace classification of labelled flats.

Made by :func:`classify_flat_surfaces` or
:func:`classify_flat_surfaces_with_scope`. ``table`` has one row per label:
``label``, ``surface_class`` (``"floodplain"``, ``"terrace"`` or
``"none"``), ``source_id``, ``n_pixels``, ``median_height``,
``mean_height`` (HAND, m), ``min_distance_to_channel_metres``,
``adjacent_to_channel``.
)doc")
        .def_prop_ro("surface_class", [](const FlatSurfaceClassification& r) { return Raster(r.native().surface_class_raster); })
        .def_prop_ro("floodplain_mask", [](const FlatSurfaceClassification& r) { return Raster(r.native().floodplain_mask); })
        .def_prop_ro("hand_band", [](const FlatSurfaceClassification& r) { return Raster(r.native().hand_band_raster); })
        .def_prop_ro("n_floodplain_labels", [](const FlatSurfaceClassification& r) { return r.native().n_floodplain_labels; })
        .def_prop_ro("n_terrace_labels", [](const FlatSurfaceClassification& r) { return r.native().n_terrace_labels; })
        .def_prop_ro("global_breaks_metres", [](const FlatSurfaceClassification& r) { return r.native().global_breaks_m; })
        .def_prop_ro("banded_labels", [](const FlatSurfaceClassification& r) { return Raster(r.native().banded_labels_raster); },
            "Patch labels after band-then-label (1..N, 0 = background); the ``patches`` labels.")
        .def_prop_ro("ladder_window_index", [](const FlatSurfaceClassification& r) { return Raster(r.native().ladder_window_index); },
            "Index into ``ladder_windows`` of the reach window each pixel was banded in (-1 = global).")
        .def_prop_ro("n_ladder_windows_trained", [](const FlatSurfaceClassification& r) { return r.native().ladder_n_windows_trained; })
        .def_prop_ro("table", &class_table, nb::sig("def table(self) -> pandas.DataFrame"),
            "One row per flat label; columns as in the class documentation.")
        .def_prop_ro("ladder_windows", &ladder_window_table, nb::sig("def ladder_windows(self) -> pandas.DataFrame"),
            "One row per ladder reach window: source_id, window_centre_metres, n_pixels, k_selected, status (resolved, unresolved, forced, degenerate, supplied), used_global_breaks, breaks_metres, band_centres_metres, required_gap_metres (kappa x E_v), vertical_accuracy_metres, median_gradient, n_trained_breaks, training_source.")
        .def_prop_ro("patches", &patch_table, nb::sig("def patches(self) -> pandas.DataFrame"),
            "One row per patch of ``banded_labels``: label, band, n_pixels, area_m2, median_field (HAND, m).")
        .def_prop_ro("training_polygons", &training_polygon_table, nb::sig("def training_polygons(self) -> pandas.DataFrame"),
            "One row per training polygon (empty without training_surfaces_fname): label, surface_class, n_cells, n_with_field, median, p10, p25, p75, p90 (HAND, m), flagged (spread too wide to be one level).");

    m.def("find_flat_surfaces", &find_flat_surfaces, nb::arg("dem"),
        nb::arg("use_madk2") = true, nb::arg("madk2_window_radius_metres") = 100.0,
        nb::arg("madk2_threshold_percentile") = 25.0,
        nb::arg("use_residual_iqr") = true, nb::arg("iqr_method") = "fftw_highpass",
        nb::arg("iqr_fit_lambda_metres") = 500.0,
        nb::arg("iqr_detrend_radius_metres") = 500.0,
        nb::arg("iqr_window_radius_metres") = 50.0,
        nb::arg("iqr_window_shape") = "circular",
        nb::arg("iqr_threshold_percentile") = 25.0,
        nb::arg("use_nvd") = false, nb::arg("nvd_polyfit_lambda") = 10.0,
        nb::arg("nvd_window_radius_metres") = 100.0,
        nb::arg("nvd_window_shape") = "circular",
        nb::arg("nvd_threshold_percentile") = 25.0,
        nb::arg("require_all_metrics_pass") = true,
        nb::arg("directional_close_line_length_pixels") = 0,
        nb::arg("directional_close_n_angles") = 18,
        nb::arg("min_component_pixels") = 100,
        nb::arg("relief_window_radius_metres") = 0.0,
        nb::arg("max_relief_metres") = 1.0,
        nb::arg("remove_risers") = false, nb::arg("n_threads") = 0,
        R"doc(Detect and label flat surfaces from texture and local relief.

Each enabled roughness metric is computed in a moving window and thresholded
at a percentile of its own distribution (low roughness = flat); the metric
masks are combined, local relief must also be below
``max_relief_metres``, and connected flats are labelled.

Parameters
----------
dem : Raster
    Elevation raster (m).
use_madk2 : bool, default=True
    Use MADk2 roughness (median absolute deviation of second differences;
    insensitive to uniform slope).
madk2_window_radius_metres : float, default=100.0
    MADk2 window radius (m).
madk2_threshold_percentile : float, default=25.0
    Pixels below this percentile of MADk2 count as flat.
use_residual_iqr : bool, default=True
    Use the interquartile range of residuals to a smooth reference surface.
iqr_method : {"fftw_highpass", "wavelet_quadratic", "moving_window_mean"}, default="fftw_highpass"
    How the reference surface is removed (the first two need FFTW).
iqr_fit_lambda_metres : float, default=500.0
    Reference-surface wavelength (m) for ``fftw_highpass`` and
    ``wavelet_quadratic``.
iqr_detrend_radius_metres : float, default=500.0
    Moving-window radius (m) for ``moving_window_mean``.
iqr_window_radius_metres : float, default=50.0
    Radius (m) of the window in which the IQR is computed.
iqr_window_shape : {"circular", "square"}, default="circular"
    Shape of that window.
iqr_threshold_percentile : float, default=25.0
    Pixels below this percentile of the IQR count as flat.
use_nvd : bool, default=False
    Use normal-vector dispersion (needs FFTW).
nvd_polyfit_lambda : float, default=10.0
    Wavelength (m) of the surface fit giving the normals.
nvd_window_radius_metres : float, default=100.0
    NVD window radius (m).
nvd_window_shape : {"circular", "square"}, default="circular"
    NVD window shape.
nvd_threshold_percentile : float, default=25.0
    Pixels below this percentile of NVD count as flat.
require_all_metrics_pass : bool, default=True
    True: a pixel must pass every enabled metric (AND); False: any (OR).
directional_close_line_length_pixels : int, default=0
    Length (pixels) of the line used to close small gaps in narrow non-flat
    corridors (channels, scarps) before labelling; below 3 disables it.
directional_close_n_angles : int, default=18
    Number of line orientations tried (18 = 10-degree steps).
min_component_pixels : int, default=100
    Smallest labelled flat (pixels).
relief_window_radius_metres : float, default=0.0
    Radius (m) of the square local-relief window, converted to the nearest
    whole number of cells (at least one); a radius of r cells is a
    (2r + 1)-cell side. 0 means 2 cells (a 5 x 5 window). Replaces the
    deprecated ``relief_window_side_pixels`` (a full side in cells).
max_relief_metres : float, default=1.0
    Pixels with local relief (max - min, m) at or above this are not flat.
    Replaces the deprecated ``relief_threshold_metres``.
remove_risers : bool, default=False
    Fit a plane to each flat and remove pixels with high residual relief
    (risers between terrace levels), then relabel.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
FlatSurfaceResult
    Roughness rasters, per-metric and combined masks, and labelled flats.
)doc");
    m.def("classify_flat_surfaces", &classify_flat_surfaces,
        nb::arg("labelled_flat_mask"), nb::arg("height_above_channel"),
        nb::arg("channel_source_ids"), nb::arg("max_floodplain_hand_metres") = 1.0,
        nb::arg("min_terrace_hand_metres") = 0.0,
        nb::arg("channel_adjacency_radius_metres") = 0.0,
        nb::arg("valley_corridor_half_width_metres") = 0.0,
        nb::arg("enable_leak_split") = true, nb::arg("enable_hand_banding") = false,
        nb::arg("donor_flat") = nb::none(), nb::arg("channel_flow_distance") = nb::none(),
        nb::arg("n_threads") = 0, R"doc(
Classify labelled flats as floodplain or terrace from HAND and channel IDs.

``channel_source_ids`` is a co-registered Raster containing non-negative
source IDs on in-scope channel pixels and negative values elsewhere.

Parameters
----------
labelled_flat_mask : Raster
    Flat labels (``FlatSurfaceResult.labelled_flat_mask``).
height_above_channel : Raster
    HAND (m) on the same grid (``FlowInfo.height_above_channel``).
channel_source_ids : Raster
    Source ID of each in-scope channel pixel, negative elsewhere.
max_floodplain_hand_metres : float, default=1.0
    A channel-adjacent label is floodplain when its median HAND (m) is at
    most this. Replaces the deprecated ``floodplain_max_height_metres``.
min_terrace_hand_metres : float, default=0.0
    A label linked to a channel is terrace when its median HAND (m) is at
    least this (and it is not floodplain). Replaces the deprecated
    ``terrace_min_height_metres``.
channel_adjacency_radius_metres : float, default=0.0
    A label is channel-adjacent when a pixel lies within this distance (m)
    of an in-scope channel; 0 means touching (8-neighbour).
valley_corridor_half_width_metres : float, default=0.0
    A terrace must have a pixel within this distance (m) of a channel,
    excluding upland flats; 0 means no limit (lsdtt-fluvial-surfaces
    defaults to 5000). Replaces the deprecated
    ``terrace_max_distance_to_channel_metres``.
enable_leak_split : bool, default=True
    Split a label whose HAND distribution is bimodal (a floodplain leaking
    into a terrace across a broken riser) at the antimode.
enable_hand_banding : bool, default=False
    Split labels into reach-relative HAND bands along the channel before
    classifying.
donor_flat : array-like of int or None, default=None
    HAND banding only: flat index of the channel pixel each pixel drains to
    (``FlowInfo.height_above_channel_with_donor()["channel_donor"]``).
channel_flow_distance : Raster or None, default=None
    HAND banding only: flow distance (m) used to place donors along the
    channel.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
FlatSurfaceClassification
)doc");
    m.def("classify_flat_surfaces_with_scope", &classify_flat_surfaces_with_scope,
        nb::arg("labelled_flat_mask"), nb::arg("dem"), nb::arg("scope"),
        nb::arg("max_floodplain_hand_metres") = 1.0,
        nb::arg("min_terrace_hand_metres") = 0.0,
        nb::arg("channel_adjacency_radius_metres") = 0.0,
        nb::arg("valley_corridor_half_width_metres") = 0.0,
        nb::arg("enable_leak_split") = true, nb::arg("enable_hand_banding") = false,
        nb::arg("n_threads") = 0,
        nb::arg("ladder_break_algorithm_choice") = "antimode",
        nb::arg("ladder_level_count") = 0, nb::arg("ladder_max_level_count") = 12,
        nb::arg("vertical_accuracy_metres") = 0.0,
        nb::arg("vertical_accuracy_algorithm_choice") = "fixed",
        nb::arg("vertical_accuracy_gentle_metres") = 0.0,
        nb::arg("vertical_accuracy_steep_metres") = 0.0,
        nb::arg("vertical_accuracy_break_gradient") = 0.20,
        nb::arg("vertical_accuracy_gradient_window_radius_metres") = 0.0,
        nb::arg("ladder_resolvability_kappa") = 3.0,
        nb::arg("ladder_k_bootstrap") = 0, nb::arg("ladder_k_bootstrap_seed") = 1,
        nb::arg("ladder_reach_window_metres") = 0.0,
        nb::arg("ladder_min_window_pixels") = 1000, nb::arg("ladder_hist_bins") = 64,
        nb::arg("ladder_min_field_metres") = 0.05,
        nb::arg("ladder_max_valley_fraction") = 0.6,
        nb::arg("ladder_gap_break_position") = "midpoint",
        nb::arg("ladder_gap_rank") = "lowest",
        nb::arg("ladder_gap_trough_max_ratio") = 0.6,
        nb::arg("ladder_gap_empty_fraction") = 0.02,
        nb::arg("ladder_gap_min_width_metres") = 0.5,
        nb::arg("ladder_gap_hist_bins") = 80,
        nb::arg("ladder_max_separation_ratio") = 0.0,
        nb::arg("ladder_min_break_metres") = 0.0,
        nb::arg("ladder_max_break_metres") = 0.0,
        nb::arg("min_band_component_pixels") = 1,
        nb::arg("ladder_unresolved_window_algorithm_choice") = "single_band",
        nb::arg("clip_corridor_per_pixel") = true,
        nb::arg("max_band_spread_metres") = 50.0,
        nb::arg("segmentation_algorithm_choice") = "within_bands",
        nb::arg("riser_mask") = nb::none(),
        nb::arg("break_flats_with_channels") = false,
        nb::arg("supplied_levels") = std::map<std::int32_t, std::vector<double>>{},
        nb::arg("training_surfaces_fname") = "",
        nb::arg("training_surfaces_layer_name") = "",
        nb::arg("training_surfaces_class_field_name") = "class",
        nb::arg("training_spread_factor") = 3.0,
        nb::arg("training_max_relative_spread") = 0.25,
        nb::arg("min_training_cells_per_class") = 20, R"doc(
Classify flat surfaces using a ``ChannelScope``: the HAND ladder of
``lsdtt-fluvial-surfaces``.

HAND, its nearest-channel donors, flow distance, channel mask, and source IDs
are derived from the scope's retained D8 flow topology. Pixels are banded by
HAND in reach windows along each channel (the reference-field ladder), each
band is labelled into patches, and patches are classed floodplain or
terrace. The ladder keywords have the names and meanings of the
lsdtt-fluvial-surfaces parameters; see that CLI's reference for the method.
Two defaults differ from the CLI so earlier results are unchanged:
``min_band_component_pixels`` (CLI 50) and
``ladder_unresolved_window_algorithm_choice`` (CLI ``nearest_resolved``).

Parameters
----------
labelled_flat_mask : Raster
    Flat labels (``FlatSurfaceResult.labelled_flat_mask``).
dem : Raster
    Elevation (m) on the scope's flow grid.
scope : ChannelScope
    In-scope channels (``FlowInfo.channel_scope``).
max_floodplain_hand_metres : float, default=1.0
    A channel-adjacent label is floodplain when its median HAND (m) is at
    most this. Replaces the deprecated ``floodplain_max_height_metres``.
min_terrace_hand_metres : float, default=0.0
    A label linked to a channel is terrace when its median HAND (m) is at
    least this (and it is not floodplain). Replaces the deprecated
    ``terrace_min_height_metres``.
channel_adjacency_radius_metres : float, default=0.0
    A label is channel-adjacent when a pixel lies within this distance (m)
    of an in-scope channel; 0 means touching (8-neighbour).
valley_corridor_half_width_metres : float, default=0.0
    A terrace must have a pixel within this distance (m) of a channel,
    excluding upland flats; 0 means no limit (lsdtt-fluvial-surfaces
    defaults to 5000). Replaces the deprecated
    ``terrace_max_distance_to_channel_metres``.
enable_leak_split : bool, default=True
    Split a label whose HAND distribution is bimodal (a floodplain leaking
    into a terrace across a broken riser) at the antimode.
enable_hand_banding : bool, default=False
    Split labels into reach-relative HAND bands along the channel before
    classifying.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.
ladder_break_algorithm_choice : {"antimode", "jenks", "gap"}, default="antimode"
    How each window's HAND population is split into levels.
ladder_level_count : int, default=0
    Force this many levels; 0 picks the largest resolvable number.
ladder_max_level_count : int, default=12
    Most levels tried.
vertical_accuracy_metres : float, default=0.0
    DEM vertical accuracy E_v (m); levels closer than
    ``ladder_resolvability_kappa`` x E_v are not resolved. 0 disables it.
vertical_accuracy_algorithm_choice : {"fixed", "slope_dependent"}, default="fixed"
    ``slope_dependent`` uses ``vertical_accuracy_gentle_metres`` below and
    ``vertical_accuracy_steep_metres`` above ``vertical_accuracy_break_gradient``
    (m/m), the window gradient measured by a polyfit of radius
    ``vertical_accuracy_gradient_window_radius_metres`` (0: 2.5 cells).
vertical_accuracy_gentle_metres, vertical_accuracy_steep_metres : float, default=0.0
    E_v (m) on gentle and steep ground, for ``slope_dependent``.
vertical_accuracy_break_gradient : float, default=0.20
    Gradient (m/m) between gentle and steep.
vertical_accuracy_gradient_window_radius_metres : float, default=0.0
    Polyfit radius (m) of that gradient; 0 means 2.5 cells.
ladder_resolvability_kappa : float, default=3.0
    Multiple of E_v two levels must be apart.
ladder_k_bootstrap, ladder_k_bootstrap_seed : int, default=0, 1
    Bootstrap resamples of each window's level count (0 = none), and seed.
ladder_reach_window_metres : float, default=0.0
    Reach window length (m) along the channel; 0 is max(2000 m, span / 8).
ladder_min_window_pixels : int, default=1000
    Smallest window population; smaller windows use the global breaks.
ladder_hist_bins : int, default=64
    Histogram bins of the antimode finder.
ladder_min_field_metres : float, default=0.05
    HAND (m) below which pixels are not banded.
ladder_max_valley_fraction : float, default=0.6
    Antimode: a dip is a break when its smoothed density is below this
    fraction of the smaller flanking peak.
ladder_gap_break_position : {"lower", "midpoint", "upper"}, default="midpoint"
    Where in an empty gap the gap finder puts the break.
ladder_gap_rank : {"lowest", "width"}, default="lowest"
    Which gaps the gap finder uses first.
ladder_gap_trough_max_ratio : float, default=0.6
    Gap finder: a local minimum is a candidate break when its floor is below
    this fraction of the smaller flanking peak; 0 means empty spans only.
ladder_gap_empty_fraction : float, default=0.02
    Gap finder: a bin below this fraction of the tallest bin is empty.
ladder_gap_min_width_metres : float, default=0.5
    Gap finder: narrowest gap (m).
ladder_gap_hist_bins : int, default=80
    Gap finder histogram bins.
ladder_max_separation_ratio : float, default=0.0
    Reject a partition whose worst break sits where the density is above
    this fraction of the populations either side; 0 disables it.
ladder_min_break_metres, ladder_max_break_metres : float, default=0.0
    Breaks outside this HAND range (m) are dropped; 0 disables each bound.
min_band_component_pixels : int, default=1
    Smallest patch kept within a band (lsdtt-fluvial-surfaces default 50).
ladder_unresolved_window_algorithm_choice : {"single_band", "nearest_resolved"}, default="single_band"
    A window with no resolvable levels gets one band, or borrows the
    nearest resolved window's breaks (the CLI default).
clip_corridor_per_pixel : bool, default=True
    Apply the valley corridor per pixel as well as per patch, so a patch
    touching the corridor edge does not carry its class across an
    interfluve.
max_band_spread_metres : float, default=50.0
    A patch whose robust HAND spread still exceeds this (m) cannot take the
    terrace class; 0 disables it.
segmentation_algorithm_choice : {"within_bands", "whole_mask"}, default="within_bands"
    Label within each band, or segment the whole admitted mask (with
    ``riser_mask`` as barriers).
riser_mask : Raster or None, default=None
    ``whole_mask`` only: riser cells (> 0.5) that patches may not cross,
    for example ``find_risers(...).centreline``.
break_flats_with_channels : bool, default=False
    Remove channel cells from the flats before classifying.
supplied_levels : dict of int to list of float, default={}
    Level heights (HAND, m) per channel source id, banded at the midpoints
    between levels (lsdtt-fluvial-surfaces ``supplied_levels``; for
    example ``ValleyCrossSections.levels()`` grouped by source).
training_surfaces_fname : str, default=""
    Vector file of floodplain / terraceN training polygons that set the
    HAND breaks per window (lsdtt-fluvial-surfaces training).
training_surfaces_layer_name : str, default=""
    Its layer; empty reads the first.
training_surfaces_class_field_name : str, default="class"
    Its class field (``floodplain``, ``terrace`` or ``terraceN``).
training_spread_factor, training_max_relative_spread : float, default=3.0, 0.25
    A training polygon whose HAND p75 - p25 exceeds both
    ``training_spread_factor`` x the global required gap and
    ``training_max_relative_spread`` x its median is flagged as not one
    level.
min_training_cells_per_class : int, default=20
    Fewest core cells a training population needs in a window.

Returns
-------
FlatSurfaceClassification
)doc");
}

} // namespace lsdtt3_python
