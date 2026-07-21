import os

from lsdtt_qgis.naming import expected_output_path, extension_for_format

# Output suffixes written by lsdtt-raster-preprocessing, keyed by the
# boolean parameter-file key that enables each one.  Kept in sync with
# cli/lsdtt-raster-preprocessing.cpp.
RASTER_PREPROCESSING_RASTER_SUFFIXES = {
    "write_fill_raster": "_fill.tif",
    "write_fill_diff_raster": "_fill_diff.tif",
    "write_breached_raster": "_breached.tif",
    "write_carve_diff_raster": "_carve_diff.tif",
    "write_hillshade_raster": "_hillshade.tif",
    "write_multidirectional_hillshade_raster": "_hillshade_multidirectional.tif",
    "write_filled_nodata_raster": "_nodata_filled.tif",
    "write_clipped_rectangle_raster": "_clipped_rect.tif",
    "apply_binary_mask": "_binary_masked.tif",
    "clip_raster_to_polygon": "_polygon_clipped.tif",
}


def expected_raster_preprocessing_outputs(output_folder, output_prefix, enabled_keys):
    return {
        key: expected_output_path(output_folder, output_prefix, suffix)
        for key, suffix in RASTER_PREPROCESSING_RASTER_SUFFIXES.items()
        if key in enabled_keys
    }


def test_expected_outputs_include_only_enabled_layers(tmp_path):
    outputs = expected_raster_preprocessing_outputs(
        str(tmp_path), "demo", {"write_fill_raster", "write_hillshade_raster"}
    )
    assert outputs == {
        "write_fill_raster": os.path.join(str(tmp_path), "demo_fill.tif"),
        "write_hillshade_raster": os.path.join(str(tmp_path), "demo_hillshade.tif"),
    }


def test_expanded_outputs_resolve_expected_suffixes(tmp_path):
    outputs = expected_raster_preprocessing_outputs(
        str(tmp_path),
        "demo",
        set(RASTER_PREPROCESSING_RASTER_SUFFIXES),
    )
    assert outputs["write_multidirectional_hillshade_raster"].endswith(
        "demo_hillshade_multidirectional.tif"
    )
    assert outputs["write_filled_nodata_raster"].endswith("demo_nodata_filled.tif")
    assert outputs["clip_raster_to_polygon"].endswith("demo_polygon_clipped.tif")
    assert len(outputs) == len(RASTER_PREPROCESSING_RASTER_SUFFIXES)


# Suffixes written by lsdtt-channel-extraction, keyed by enabling flag.
# Kept in sync with cli/lsdtt-channel-extraction.cpp.
CHANNEL_EXTRACTION_SUFFIXES = {
    "write_channel_network_lines": "_channel_network_lines",
    "write_channel_network_points": "_channel_network_points",
    "write_channel_network_junctions": "_channel_network_junctions",
    "write_flow_lines": "_flow_lines",
    "write_arcmap_d8_flowdirection_raster": "_arcmap_d8_flowdirection.tif",
    "write_burned_stream_profiles": "_burned_stream_profiles",
    "wiener_write_filtered_dem": "_wiener_filtered_dem.tif",
    "wiener_write_tangential_curvature": "_wiener_tangential_curvature.tif",
    "wiener_write_qq_table": "_wiener_qq.csv",
}

# Suffixes written by lsdtt-basin-extraction.  Basin outlines switch between
# polygon and line suffixes on write_basin_outlines_as_polygon; the basins
# raster is _junction_basins.tif (NOT _basins.tif).
# Kept in sync with cli/lsdtt-basin-extraction.cpp.
BASIN_EXTRACTION_SUFFIXES = {
    "write_basin_outlines(polygon)": "_junction_basin_polygons",
    "write_basin_outlines(lines)": "_junction_basin_outlines",
    "write_basin_outlet_points": "_basin_outlet_points.fgb",
    "write_basins_raster": "_junction_basins.tif",
    "write_basins_datamap": "_junction_basins_datamap.csv",
    "write_basin_outlines_to_csv": "_basin_outlines.csv",
    "write_selected_basin_outlets_csv": "_selected_basin_outlets.csv",
    "select_basins_from_line_crossings": "_line_crossing_basins.csv",
    "write_side_tributary_main_channel": "_side_tributary_main_channel",
    "write_side_tributary_outlets": "_side_tributary_outlets",
    "write_side_tributary_basin_outlines": "_side_tributary_basin_outlines",
}


def test_basin_outline_suffixes_depend_on_polygon_flag():
    polygon = expected_output_path(
        "/out", "demo", BASIN_EXTRACTION_SUFFIXES["write_basin_outlines(polygon)"] + ".fgb"
    )
    lines = expected_output_path(
        "/out", "demo", BASIN_EXTRACTION_SUFFIXES["write_basin_outlines(lines)"] + ".fgb"
    )
    assert polygon.endswith("demo_junction_basin_polygons.fgb")
    assert lines.endswith("demo_junction_basin_outlines.fgb")


def test_basins_raster_uses_junction_basins_suffix():
    path = expected_output_path(
        "/out", "demo", BASIN_EXTRACTION_SUFFIXES["write_basins_raster"]
    )
    assert path.endswith("demo_junction_basins.tif")


def test_channel_extraction_flow_outputs():
    assert expected_output_path(
        "/out", "demo",
        CHANNEL_EXTRACTION_SUFFIXES["write_flow_lines"]
        + "." + extension_for_format("flatgeobuf"),
    ).endswith("demo_flow_lines.fgb")
    assert expected_output_path(
        "/out", "demo",
        CHANNEL_EXTRACTION_SUFFIXES["write_arcmap_d8_flowdirection_raster"],
    ).endswith("demo_arcmap_d8_flowdirection.tif")
    assert expected_output_path(
        "/out", "demo",
        CHANNEL_EXTRACTION_SUFFIXES["write_burned_stream_profiles"]
        + "." + extension_for_format("csv"),
    ).endswith("demo_burned_stream_profiles.csv")


# Suffixes written by lsdtt-chi-analysis, keyed by enabling flag.
# Kept in sync with cli/lsdtt-chi-analysis.cpp and
# src/cpu/lsdtt_chi_data_maps_cpu.cpp.
CHI_ANALYSIS_SUFFIXES = {
    "fit_chi_networks": "_chi_disorder_sweep.csv",
    "print_chi_data_maps": "_chi_data_maps",
    "write_chi_data_map_lines": "_chi_data_map_lines",
    "run_inversion(best)": "_chi_ksn_inversion_best.csv",
    "run_inversion(trials)": "_chi_ksn_inversion_trials.csv",
}


def test_chi_analysis_suffixes():
    assert expected_output_path(
        "/out", "demo", CHI_ANALYSIS_SUFFIXES["fit_chi_networks"]
    ).endswith("demo_chi_disorder_sweep.csv")
    assert expected_output_path(
        "/out", "demo",
        CHI_ANALYSIS_SUFFIXES["print_chi_data_maps"]
        + "." + extension_for_format("csv"),
    ).endswith("demo_chi_data_maps.csv")
    assert expected_output_path(
        "/out", "demo",
        CHI_ANALYSIS_SUFFIXES["write_chi_data_map_lines"]
        + "." + extension_for_format("flatgeobuf"),
    ).endswith("demo_chi_data_map_lines.fgb")
    assert expected_output_path(
        "/out", "demo", CHI_ANALYSIS_SUFFIXES["run_inversion(best)"]
    ).endswith("demo_chi_ksn_inversion_best.csv")


def test_vector_output_suffixes():
    for fmt, ext in [("flatgeobuf", "fgb"), ("gpkg", "gpkg")]:
        assert expected_output_path(
            "/out", "demo", "_depressions.{}".format(extension_for_format(fmt))
        ).endswith("demo_depressions." + ext)
        assert expected_output_path(
            "/out", "demo", "_clipped_vector.{}".format(extension_for_format(fmt))
        ).endswith("demo_clipped_vector." + ext)
        assert expected_output_path(
            "/out", "demo", "_sampled_points.{}".format(extension_for_format(fmt))
        ).endswith("demo_sampled_points." + ext)


# Suffixes written by lsdtt-basin-stats, keyed by enabling flag ("always" for
# the statistics table, which is written on every run).
# Kept in sync with cli/lsdtt-basin-stats.cpp.
BASIN_STATS_SUFFIXES = {
    "always": "_basin_stats.csv",
    "write_basin_polygons": "_basin_polygons",
    "calculate_junction_area_ratios": "_junction_area_ratios.csv",
    "calculate_junction_area_ratios(summary)": "_junction_area_ratio_summary.csv",
    "write_hillshade_raster": "_hillshade.tif",
}


def test_basin_stats_suffixes():
    assert expected_output_path(
        "/out", "demo", BASIN_STATS_SUFFIXES["always"]
    ).endswith("demo_basin_stats.csv")
    assert expected_output_path(
        "/out", "demo",
        BASIN_STATS_SUFFIXES["write_basin_polygons"]
        + "." + extension_for_format("flatgeobuf"),
    ).endswith("demo_basin_polygons.fgb")
    assert expected_output_path(
        "/out", "demo", BASIN_STATS_SUFFIXES["calculate_junction_area_ratios"]
    ).endswith("demo_junction_area_ratios.csv")


# Suffixes written by lsdtt-channel-network-analysis, keyed by enabling flag.
# The bearing/angle/channel layers take the chosen vector extension; the
# summary tables are always CSV.
# Kept in sync with cli/lsdtt-channel-network-analysis.cpp.
CHANNEL_NETWORK_ANALYSIS_SUFFIXES = {
    "write_channel_link_bearings": "_channel_link_bearings",
    "write_junction_angles": "_junction_angles",
    "write_strahler_channels": "_strahler_channels",
    "write_tokunaga_channels": "_tokunaga_channels",
    "write_tokunaga_channel_data": "_tokunaga_channel_data.csv",
    "write_tokunaga_summary_data": "_tokunaga_summary.csv",
    "write_junction_area_ratio_data": "_junction_area_ratios.csv",
    "write_junction_area_ratio_summary": "_junction_area_ratio_summary.csv",
    "write_junction_centre_of_mass_bearings": "_junction_com_bearings.csv",
    "write_channel_link_bearing_vectors": "_channel_link_bearing_vectors",
    "write_junction_angle_bearing_vectors": "_junction_angle_bearing_vectors",
}


def test_channel_network_analysis_layer_suffixes_follow_format():
    for fmt, ext in [("flatgeobuf", "fgb"), ("gpkg", "gpkg"), ("csv", "csv")]:
        assert expected_output_path(
            "/out", "demo",
            CHANNEL_NETWORK_ANALYSIS_SUFFIXES["write_channel_link_bearings"]
            + "." + extension_for_format(fmt),
        ).endswith("demo_channel_link_bearings." + ext)
        assert expected_output_path(
            "/out", "demo",
            CHANNEL_NETWORK_ANALYSIS_SUFFIXES["write_tokunaga_channels"]
            + "." + extension_for_format(fmt),
        ).endswith("demo_tokunaga_channels." + ext)


def test_channel_network_analysis_summary_tables_are_csv():
    for key in (
        "write_tokunaga_summary_data",
        "write_junction_area_ratio_data",
        "write_junction_centre_of_mass_bearings",
    ):
        assert CHANNEL_NETWORK_ANALYSIS_SUFFIXES[key].endswith(".csv")
    assert expected_output_path(
        "/out", "demo",
        CHANNEL_NETWORK_ANALYSIS_SUFFIXES["write_junction_centre_of_mass_bearings"],
    ).endswith("demo_junction_com_bearings.csv")


# Suffixes written by lsdtt-surface-metrics. Polyfit-derived rasters carry the
# fit method in the suffix; the wrapper pins the method to "direct".
# Kept in sync with cli/lsdtt-surface-metrics.cpp.
SURFACE_METRICS_SUFFIXES = {
    "write_hillshade_raster": "_hillshade.tif",
    "write_local_relief_raster": "_local_relief.tif",
    "write_shape_classification": "_direct_shape_class.tif",
    "write_shape_curvature_rasters": "_direct_shape_gaussian_curvature.tif",
    "write_madk2_rasters": "_madk2_ra.tif",
    "write_cloth_envelope_raster": "_cloth_envelope.tif",
    "write_cloth_incision_raster": "_cloth_incision.tif",
    "write_cloth_crest_mask_raster": "_cloth_crest_mask.tif",
}


def polyfit_metric_suffix(metric, method="direct"):
    """Mirror of polyfit_metric_suffix() in cli/lsdtt-surface-metrics.cpp."""
    return "_{}_polyfit_{}.tif".format(method, metric)


def test_surface_metrics_suffixes():
    assert expected_output_path(
        "/out", "demo", SURFACE_METRICS_SUFFIXES["write_local_relief_raster"]
    ).endswith("demo_local_relief.tif")
    assert expected_output_path(
        "/out", "demo", SURFACE_METRICS_SUFFIXES["write_shape_classification"]
    ).endswith("demo_direct_shape_class.tif")
    assert expected_output_path(
        "/out", "demo", SURFACE_METRICS_SUFFIXES["write_cloth_incision_raster"]
    ).endswith("demo_cloth_incision.tif")


def test_polyfit_metric_suffix_includes_method():
    assert expected_output_path(
        "/out", "demo", polyfit_metric_suffix("plan_curvature")
    ).endswith("demo_direct_polyfit_plan_curvature.tif")
    assert polyfit_metric_suffix("slope_degrees", "wavelet") == (
        "_wavelet_polyfit_slope_degrees.tif"
    )
