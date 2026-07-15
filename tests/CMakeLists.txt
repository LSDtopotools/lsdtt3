find_package(Python3 COMPONENTS Interpreter REQUIRED)

set(LSDTT_TEST_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cli_smoke.py")
set(LSDTT_TEST_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data")
set(LSDTT_TEST_PARAM_DIR "${CMAKE_CURRENT_SOURCE_DIR}/params")

function(add_lsdtt_cli_smoke_test TEST_NAME TARGET_NAME PARAM_TEMPLATE EXPECTED_OUTPUT)
    set(work_dir "${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}")

    add_test(
        NAME ${TEST_NAME}
        COMMAND
            ${Python3_EXECUTABLE} "${LSDTT_TEST_SCRIPT}"
            --exe "$<TARGET_FILE:${TARGET_NAME}>"
            --param-template "${LSDTT_TEST_PARAM_DIR}/${PARAM_TEMPLATE}"
            --read-path "${LSDTT_TEST_DATA_DIR}"
            --work-dir "${work_dir}"
            --expect "${EXPECTED_OUTPUT}"
            ${ARGN}
    )

    set(test_labels "cli")
    if(TEST_NAME MATCHES "^rigorous_")
        list(APPEND test_labels "rigorous")
    elseif(TEST_NAME MATCHES "^failure_")
        list(APPEND test_labels "failure")
    else()
        list(APPEND test_labels "smoke")
    endif()

    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "${test_labels}")
endfunction()

add_lsdtt_cli_smoke_test(
    smoke_raster_preprocessing
    lsdtt-raster-preprocessing
    raster_preprocessing_smoke.param.in
    tiny_fill.tif
)

add_lsdtt_cli_smoke_test(
    smoke_raster_preprocessing_tif
    lsdtt-raster-preprocessing
    raster_preprocessing_tif_smoke.param.in
    tiny30_fill.tif
)

add_lsdtt_cli_smoke_test(
    smoke_channel_extraction
    lsdtt-channel-extraction
    channel_extraction_smoke.param.in
    tiny_channel_network_lines.fgb
)

add_lsdtt_cli_smoke_test(
    rigorous_channel_lines_feature_count
    lsdtt-channel-extraction
    channel_extraction_smoke.param.in
    tiny_channel_network_lines.fgb
    --expect-vector-feature-min tiny_channel_network_lines.fgb:1
)

add_lsdtt_cli_smoke_test(
    smoke_basin_extraction
    lsdtt-basin-extraction
    basin_extraction_smoke.param.in
    tiny_junction_basin_outlines.fgb
)

add_lsdtt_cli_smoke_test(
    rigorous_basin_outlines_feature_count
    lsdtt-basin-extraction
    basin_extraction_smoke.param.in
    tiny_junction_basin_outlines.fgb
    --expect-vector-feature-min tiny_junction_basin_outlines.fgb:1
)

add_lsdtt_cli_smoke_test(
    smoke_nodata_fill
    lsdtt-raster-preprocessing
    raster_preprocessing_nodata_smoke.param.in
    tiny_nodata_nodata_filled.tif
    --expect-valid-percent-100 tiny_nodata_nodata_filled.tif
)

add_lsdtt_cli_smoke_test(
    smoke_nodata_edge_unfilled
    lsdtt-raster-preprocessing
    raster_preprocessing_nodata_edge_smoke.param.in
    tiny_edge_nodata_filled.tif
    --expect-valid-percent tiny_edge_nodata_filled.tif:96
)

add_lsdtt_cli_smoke_test(
    smoke_drainage_analysis
    lsdtt-drainage-analysis
    drainage_analysis_smoke.param.in
    tiny_d8_drainage_area.tif
)

add_lsdtt_cli_smoke_test(
    smoke_chi_analysis
    lsdtt-chi-analysis
    chi_analysis_smoke.param.in
    tiny_chi_fit_nodes_mn_0p500.csv
)

add_lsdtt_cli_smoke_test(
    smoke_raster_calculator_difference
    lsdtt-raster-calculator
    raster_calculator_difference.param.in
    tiny_difference.tif
)

add_lsdtt_cli_smoke_test(
    failure_channel_extraction_invalid_threshold
    lsdtt-channel-extraction
    channel_extraction_invalid_threshold.param.in
    tiny_invalid_channel_network_lines.fgb
    --expect-failure
    --stderr-contains "FATAL: source_pixels_threshold must be > 0."
)

add_lsdtt_cli_smoke_test(
    failure_raster_calculator_missing_second
    lsdtt-raster-calculator
    raster_calculator_missing_second.param.in
    tiny_missing_second_difference.tif
    --expect-failure
    --stderr-contains "FATAL: second_raster_fname must be set."
)

add_lsdtt_cli_smoke_test(
    rigorous_tif_fill_raster_size
    lsdtt-raster-preprocessing
    raster_preprocessing_tif_smoke.param.in
    tiny30_fill.tif
    --expect-raster-size tiny30_fill.tif:32x39
)

add_lsdtt_cli_smoke_test(
    rigorous_raster_calculator_difference_checksum
    lsdtt-raster-calculator
    raster_calculator_difference.param.in
    tiny_difference.tif
    --expect-raster-checksum tiny_difference.tif:65511
)

add_lsdtt_cli_smoke_test(
    rigorous_drainage_d8_checksum
    lsdtt-drainage-analysis
    drainage_analysis_smoke.param.in
    tiny_d8_drainage_area.tif
    --expect-raster-checksum tiny_d8_drainage_area.tif:112
)

add_lsdtt_cli_smoke_test(
    rigorous_chi_nodes_sha256
    lsdtt-chi-analysis
    chi_analysis_smoke.param.in
    tiny_chi_fit_nodes_mn_0p500.csv
    --expect-sha256 tiny_chi_fit_nodes_mn_0p500.csv:5ef85f68c63b5bd16d45051be2a2a645fd48c73a98123024724e58987902adab
)
