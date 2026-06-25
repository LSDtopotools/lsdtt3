#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include "../core/lsdtt_flow_base.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "../utils/lsdtt_histogram.hpp"

namespace lsdtt {

struct wiener_source_extraction_params
{
    double pruning_drainage_area_m2 = 1000.0;
    double surface_fitting_radius_m = 6.0;
    int    connected_components_threshold_pixels = 100;

    bool write_filtered_dem = false;
    bool write_tangential_curvature = false;
    bool write_qq_table = false;
};

struct wiener_source_extraction_diagnostics
{
    std::optional<lsdtt_raster_cpu> filtered_dem;
    std::optional<lsdtt_raster_cpu> tangential_curvature;
    lsdtt_qq_data                   qq_data;
};

struct wiener_source_extraction_result
{
    std::vector<std::int64_t>            source_flats;
    wiener_source_extraction_diagnostics diagnostics;
};

[[nodiscard]] wiener_source_extraction_result extract_channel_sources_wiener(
    const lsdtt_raster_cpu&               prepared_dem,
    const lsdtt_flow_base&                flow,
    const wiener_source_extraction_params& params
);

} // namespace lsdtt
