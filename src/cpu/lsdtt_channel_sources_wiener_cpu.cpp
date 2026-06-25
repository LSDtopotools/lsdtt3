#include "lsdtt_channel_sources_wiener_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "lsdtt_raster_binary_mask_cpu.hpp"
#include "lsdtt_raster_skeleton_graph_cpu.hpp"
#include "lsdtt_raster_spectral_cpu.hpp"
#include "lsdtt_raster_surface_metrics_cpu.hpp"
#include "../utils/lsdtt_statistics.hpp"

namespace lsdtt {
namespace {

double compute_qq_threshold(
    const lsdtt_qq_data&       qq,
    const std::vector<double>& positive_values)
{
    if (qq.data_quantile.empty() || qq.normal_reference.empty())
    {
        return compute_percentile(positive_values, 90.0);
    }

    double best_deviation = -std::numeric_limits<double>::infinity();
    double best_threshold = std::numeric_limits<double>::quiet_NaN();

    const std::size_t n = std::min(qq.standard_normal_variate.size(), qq.data_quantile.size());
    for (std::size_t i = 0; i < n; ++i)
    {
        const double z  = qq.standard_normal_variate[i];
        const double dq = qq.data_quantile[i];
        const double rq = qq.normal_reference[i];

        if (!std::isfinite(z) || !std::isfinite(dq) || !std::isfinite(rq)) continue;
        if (z <= 0.0 || dq <= 0.0) continue;

        const double deviation = dq - rq;
        if (deviation > best_deviation)
        {
            best_deviation = deviation;
            best_threshold = dq;
        }
    }

    if (!std::isfinite(best_threshold) || best_threshold <= 0.0)
    {
        return compute_percentile(positive_values, 90.0);
    }

    return best_threshold;
}

lsdtt_raster_cpu make_binary_channel_mask(
    const lsdtt_raster_cpu& tangential_curvature,
    double                  threshold)
{
    const auto& meta = tangential_curvature.get_meta();
    const auto& src  = tangential_curvature.get_data();
    const auto& msk  = tangential_curvature.get_mask();

    std::vector<float> out_data(src.size(), 0.0f);
    std::vector<std::uint8_t> out_mask = msk;

    for (std::size_t i = 0; i < src.size(); ++i)
    {
        if (msk[i] == 0) continue;
        if (!std::isfinite(src[i]))
        {
            out_mask[i] = 0;
            continue;
        }

        out_data[i] = (src[i] > 0.0f && static_cast<double>(src[i]) >= threshold)
            ? 1.0f
            : 0.0f;
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}

void prune_channel_mask_by_area_m2(
    lsdtt_raster_cpu&       channel_mask,
    const lsdtt_raster_cpu& drainage_area,
    double                  min_area_m2)
{
    auto& channel_data = channel_mask.get_data();
    auto& channel_mask_valid = channel_mask.get_mask();

    const auto& area_data = drainage_area.get_data();
    const auto& area_mask = drainage_area.get_mask();

    for (std::size_t i = 0; i < channel_data.size(); ++i)
    {
        if (channel_mask_valid[i] == 0) continue;
        if (channel_data[i] <= 0.5f) continue;
        if (area_mask[i] == 0)
        {
            channel_data[i] = 0.0f;
            continue;
        }

        const double area_m2 = static_cast<double>(area_data[i]);
        if (!std::isfinite(area_m2) || area_m2 < min_area_m2)
        {
            channel_data[i] = 0.0f;
        }
    }
}

} // namespace

wiener_source_extraction_result extract_channel_sources_wiener(
    const lsdtt_raster_cpu&               prepared_dem,
    const lsdtt_flow_base&                flow,
    const wiener_source_extraction_params& params)
{
    if (params.pruning_drainage_area_m2 <= 0.0 ||
        !std::isfinite(params.pruning_drainage_area_m2))
    {
        throw std::invalid_argument(
            "extract_channel_sources_wiener: wiener_pruning_drainage_area_m2 must be finite and > 0");
    }

    if (params.surface_fitting_radius_m <= 0.0 ||
        !std::isfinite(params.surface_fitting_radius_m))
    {
        throw std::invalid_argument(
            "extract_channel_sources_wiener: wiener_surface_fitting_radius_m must be finite and > 0");
    }

    if (params.connected_components_threshold_pixels < 1)
    {
        throw std::invalid_argument(
            "extract_channel_sources_wiener: wiener_connected_components_threshold_pixels must be >= 1");
    }

    prepared_dem.get_meta().check_projected_metres();

    const lsdtt_raster_cpu filtered_dem = fftw2D_wiener_filter(prepared_dem);

    const poly_fit_coeffs coeffs = direct_polyfit(
        filtered_dem,
        params.surface_fitting_radius_m
    );

    polyfit_metric_request metric_request;
    metric_request.metrics = {polyfit_metric::tangential_curvature};
    const polyfit_metric_result metrics =
        compute_polyfit_metrics(coeffs, metric_request);

    const lsdtt_raster_cpu tangential_curvature =
        metrics.at(polyfit_metric::tangential_curvature);

    std::vector<double> positive_curvature_values;
    positive_curvature_values.reserve(tangential_curvature.full_grid_size() / 4);

    const auto& tc_data = tangential_curvature.get_data();
    const auto& tc_mask = tangential_curvature.get_mask();
    for (std::size_t i = 0; i < tc_data.size(); ++i)
    {
        if (tc_mask[i] == 0) continue;
        const double v = static_cast<double>(tc_data[i]);
        if (std::isfinite(v) && v > 0.0) positive_curvature_values.push_back(v);
    }

    if (positive_curvature_values.size() < 10)
    {
        throw std::runtime_error(
            "extract_channel_sources_wiener: insufficient positive tangential curvature pixels for QQ thresholding");
    }

    const std::size_t qq_points = std::min<std::size_t>(512, positive_curvature_values.size());
    const lsdtt_qq_data qq_data = compute_qq(positive_curvature_values, qq_points);
    const double qq_threshold = compute_qq_threshold(qq_data, positive_curvature_values);

    lsdtt_raster_cpu channel_mask =
        make_binary_channel_mask(tangential_curvature, qq_threshold);

    const lsdtt_raster_cpu drainage_area_m2(flow.get_drainage_area_raster());
    prune_channel_mask_by_area_m2(
        channel_mask,
        drainage_area_m2,
        params.pruning_drainage_area_m2
    );

    channel_mask = filter_binary_components_by_size(
        channel_mask,
        params.connected_components_threshold_pixels,
        pixel_connectivity::d8
    );

    lsdtt_raster_cpu skeleton = thin_zhang_suen(channel_mask);
    const skeleton_node_result nodes = classify_skeleton_nodes(skeleton);

    if (nodes.endpoints.empty())
    {
        wiener_source_extraction_result out;
        if (params.write_filtered_dem) out.diagnostics.filtered_dem = filtered_dem;
        if (params.write_tangential_curvature)
            out.diagnostics.tangential_curvature = tangential_curvature;
        if (params.write_qq_table) out.diagnostics.qq_data = qq_data;
        return out;
    }

    const lsdtt_raster_cpu component_labels =
        label_binary_components_by_size(skeleton, 1, pixel_connectivity::d8);
    const auto& labels_data = component_labels.get_data();
    const auto& labels_mask = component_labels.get_mask();

    struct endpoint_candidate
    {
        std::int64_t flat = -1;
        int          component = 0;
        double       drainage_area = 0.0;
    };

    std::vector<endpoint_candidate> candidates;
    candidates.reserve(nodes.endpoints.size());

    std::unordered_map<int, std::vector<std::size_t>> endpoints_by_component;

    for (const skeleton_node& endpoint : nodes.endpoints)
    {
        const std::int64_t flat = flow.flat_index(endpoint.row, endpoint.col);
        if (!flow.is_valid_cell(flat)) continue;

        const std::size_t idx = static_cast<std::size_t>(flat);
        if (labels_mask[idx] == 0) continue;
        const int component = static_cast<int>(std::lround(labels_data[idx]));
        if (component <= 0) continue;

        const double area = static_cast<double>(drainage_area_m2.get_data()[idx]);

        const std::size_t cidx = candidates.size();
        candidates.push_back(endpoint_candidate{flat, component, area});
        endpoints_by_component[component].push_back(cidx);
    }

    std::vector<bool> keep(candidates.size(), true);
    for (const auto& kv : endpoints_by_component)
    {
        const std::vector<std::size_t>& comp = kv.second;
        if (comp.size() <= 1) continue;

        std::size_t downstream_idx = comp.front();
        double downstream_area = candidates[downstream_idx].drainage_area;

        for (std::size_t i = 1; i < comp.size(); ++i)
        {
            const std::size_t cidx = comp[i];
            const double area = candidates[cidx].drainage_area;
            if (area > downstream_area)
            {
                downstream_area = area;
                downstream_idx = cidx;
            }
        }

        keep[downstream_idx] = false;
    }

    const auto& skel_data = skeleton.get_data();
    const auto& skel_mask = skeleton.get_mask();

    std::unordered_set<std::int64_t> unique_sources;
    unique_sources.reserve(candidates.size());

    for (std::size_t i = 0; i < candidates.size(); ++i)
    {
        if (!keep[i]) continue;

        const std::int64_t src = candidates[i].flat;
        if (!flow.is_valid_cell(src)) continue;

        const std::int64_t rec = flow.receiver_of_unchecked(src);
        bool is_single_pixel_channel = false;

        if (rec == k_outlet_receiver || rec == k_nodata_receiver)
        {
            is_single_pixel_channel = true;
        }
        else
        {
            const std::size_t rec_idx = static_cast<std::size_t>(rec);
            if (skel_mask[rec_idx] == 0 || skel_data[rec_idx] <= 0.5f)
            {
                is_single_pixel_channel = true;
            }
        }

        if (!is_single_pixel_channel)
        {
            unique_sources.insert(src);
        }
    }

    wiener_source_extraction_result out;
    out.source_flats.assign(unique_sources.begin(), unique_sources.end());
    std::sort(out.source_flats.begin(), out.source_flats.end());

    if (params.write_filtered_dem) out.diagnostics.filtered_dem = filtered_dem;
    if (params.write_tangential_curvature)
    {
        out.diagnostics.tangential_curvature = tangential_curvature;
    }
    if (params.write_qq_table) out.diagnostics.qq_data = qq_data;

    return out;
}

} // namespace lsdtt
