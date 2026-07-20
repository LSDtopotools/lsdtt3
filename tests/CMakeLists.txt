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

# -----------------------------------------------------------------------------
# Flow-routing benchmark vs RichDEM (docs/testing_and_benchmarking.md)
#
# OFF by default: these tests need network access for the fixtures and an
# optional RichDEM install, neither of which may turn CI red.  When enabled but
# unavailable, the driver exits 77 and CTest records a skip.
#
# RichDEM comes from conda-forge via pixi -- NOT pip, whose richdem 0.3.4 pins
# numpy<2 and will downgrade numpy in whatever environment it lands in:
#   pixi run fetch-benchmark-dem
#   pixi run -e benchmark flow-benchmark --exe build/lsdtt-raster-preprocessing
# -----------------------------------------------------------------------------
option(LSDTT_ENABLE_BENCHMARKS "Enable RichDEM flow-routing benchmark tests" OFF)

if(LSDTT_ENABLE_BENCHMARKS)
    # Launch through `pixi run`, not the environment's python directly.  pixi
    # exports PROJ_DATA / GDAL_DATA when it activates the environment; invoking
    # .pixi/envs/benchmark/bin/python straight from CTest skips activation and
    # GDAL then fails to open the PROJ database at write time.
    find_program(LSDTT_PIXI_EXECUTABLE pixi
                 HINTS "$ENV{HOME}/.pixi/bin" "$ENV{PIXI_HOME}/bin")

    if(NOT LSDTT_PIXI_EXECUTABLE)
        message(WARNING
            "LSDTT_ENABLE_BENCHMARKS is ON but pixi was not found. "
            "Install pixi (https://pixi.sh) or set LSDTT_PIXI_EXECUTABLE. "
            "Benchmark tests will not be registered.")
    else()
        set(LSDTT_BENCHMARK_LAUNCHER
            "${LSDTT_PIXI_EXECUTABLE}" run
            --manifest-path "${CMAKE_SOURCE_DIR}/pixi.toml"
            -e benchmark python)

    # accum_d8 is included because it is asserted against TopoToolbox, whose D8
    # uses the same steepest-descent rule as lsdtt3. The RichDEM row for the
    # same comparison is recorded but not asserted -- see REPORT_ONLY in the
    # driver and docs/testing_and_benchmarking.md.
    foreach(bench_fixture alpine_small lakeland_small)
        foreach(bench_case fill fill_epsilon breach accum_d8)
            add_test(
                NAME benchmark_${bench_case}_${bench_fixture}
                COMMAND ${LSDTT_BENCHMARK_LAUNCHER}
                        "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_flow_benchmark.py"
                        --exe "$<TARGET_FILE:lsdtt-raster-preprocessing>"
                        --fixture ${bench_fixture}
                        --comparison ${bench_case}
                        --work-dir "${CMAKE_CURRENT_BINARY_DIR}/flow_benchmark"
                        --check
            )
            set_tests_properties(benchmark_${bench_case}_${bench_fixture}
                PROPERTIES LABELS "benchmark;flow" SKIP_RETURN_CODE 77)
        endforeach()
    endforeach()
    endif()
endif()

# -----------------------------------------------------------------------------
# C++ unit tests
# -----------------------------------------------------------------------------

add_executable(lsdtt-test-param-parser-help cpp/test_param_parser_help.cpp)
target_link_libraries(lsdtt-test-param-parser-help PRIVATE lsdtt3_core)

add_test(NAME unit_param_parser_help COMMAND lsdtt-test-param-parser-help)
set_tests_properties(unit_param_parser_help PROPERTIES LABELS "unit")

# Cosmogenic physics driver (docs/cosmo_cli_testing_plan.md).  Writes its
# plottable CSVs into the test working directory.
add_executable(lsdtt-test-cosmo-physics cpp/test_cosmo_physics.cpp)
target_link_libraries(lsdtt-test-cosmo-physics PRIVATE lsdtt3_cpu)

# Flow routing / depression handling property tests
# (docs/testing_and_benchmarking.md).  Synthetic fixtures only: no
# network, no RichDEM, so this runs in the default CTest suite.
add_executable(lsdtt-test-flow-properties cpp/test_flow_properties.cpp)
target_link_libraries(lsdtt-test-flow-properties PRIVATE lsdtt3_cpu)

add_test(NAME unit_flow_properties
         COMMAND lsdtt-test-flow-properties
                 ${CMAKE_CURRENT_BINARY_DIR}/flow_properties_output)
set_tests_properties(unit_flow_properties PROPERTIES LABELS "unit;flow")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/flow_properties_output)

# -----------------------------------------------------------------------------
# lsdtt-fetch-raster CLI (docs/testing_and_benchmarking.md)
#
# These run offline: they exercise the argument handling and validation that
# happen before any network access is attempted.
# -----------------------------------------------------------------------------
add_lsdtt_cli_smoke_test(
    failure_fetch_raster_missing_bbox
    lsdtt-fetch-raster
    fetch_raster_missing_bbox.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: no area given"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_partial_bbox
    lsdtt-fetch-raster
    fetch_raster_partial_bbox.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: north must be set"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_bad_source
    lsdtt-fetch-raster
    fetch_raster_bad_source.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: unknown dem_source"
)

# -h must work without a parameter file, and without touching the network.
add_test(NAME smoke_fetch_raster_help
         COMMAND lsdtt-fetch-raster -h)
set_tests_properties(smoke_fetch_raster_help PROPERTIES
    LABELS "smoke;cli;fetch"
    PASS_REGULAR_EXPRESSION "Bounding box")

# The CLI must never grow a credential parameter: parameter files get committed
# and print_parameters() echoes every parameter at the end of a run. This guards
# the rule stated in docs/testing_and_benchmarking.md and in the CLI's own header.
add_test(NAME rigorous_fetch_raster_no_credential_parameters
         COMMAND lsdtt-fetch-raster -h)
set_tests_properties(rigorous_fetch_raster_no_credential_parameters PROPERTIES
    LABELS "rigorous;cli;fetch"
    FAIL_REGULAR_EXPRESSION "api_key;API_Key;password;secret;token;credential")

add_lsdtt_cli_smoke_test(
    failure_fill_vincent_needs_zero_slope
    lsdtt-raster-preprocessing
    fill_vincent_needs_zero_slope.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "cannot impose a gradient across flats"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_bounds_conflict
    lsdtt-fetch-raster
    fetch_raster_bounds_conflict.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: give either bounds_fname or west/south/east/north"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_bounds_three_points
    lsdtt-fetch-raster
    fetch_raster_bounds_three_points.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "exactly two"
)

if(LSDTT_ENABLE_NETWORK_TESTS)
    add_lsdtt_cli_smoke_test(
        network_fetch_raster_bounds_cli
        lsdtt-fetch-raster
        fetch_raster_bounds_network.param.in
        frompoly_DEM.tif
        --expect-raster-size frompoly_DEM.tif:128x128
    )
    set_tests_properties(network_fetch_raster_bounds_cli
        PROPERTIES LABELS "network;cli;fetch")

    add_lsdtt_cli_smoke_test(
        network_fetch_raster_cli
        lsdtt-fetch-raster
        fetch_raster_network.param.in
        ortler_test_DEM.tif
        --expect-raster-size ortler_test_DEM.tif:128x128
    )
    set_tests_properties(network_fetch_raster_cli PROPERTIES LABELS "network;cli;fetch")
endif()

# Flow-routing performance (docs/testing_and_benchmarking.md).
#
# Reported, never asserted: wall-clock thresholds in CI are flaky on shared
# runners. OFF by default because the sweep takes minutes; the value is a
# record comparable across commits, each row carrying a git SHA.
#
# Timings are IN-PROCESS. The lsdtt_seconds column of benchmark_results.csv is
# whole-process CLI wall time and is NOT comparable with RichDEM's in-routine
# numbers -- that is precisely what this driver exists to fix.
option(LSDTT_ENABLE_PERF_TESTS "Enable flow-routing performance benchmarks" OFF)

add_executable(lsdtt-test-flow-performance cpp/test_flow_performance.cpp)
target_link_libraries(lsdtt-test-flow-performance PRIVATE lsdtt3_cpu)

if(LSDTT_ENABLE_PERF_TESTS)
    add_test(NAME performance_flow_routines
             COMMAND lsdtt-test-flow-performance
                     ${CMAKE_CURRENT_BINARY_DIR}/perf_output --quick)
    set_tests_properties(performance_flow_routines
        PROPERTIES LABELS "performance;flow" TIMEOUT 1800)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/perf_output)
endif()

# Credential handling (docs/testing_and_benchmarking.md).
#
# lsdtt_secret is the only place in lsdtt3 that touches secrets, so it gets its
# own driver. Offline and deterministic; always in the default suite.
add_executable(lsdtt-test-secret cpp/test_secret.cpp)
target_link_libraries(lsdtt-test-secret PRIVATE lsdtt3_core)
target_include_directories(lsdtt-test-secret PRIVATE ${CMAKE_SOURCE_DIR}/src/core)

add_test(NAME unit_secret
         COMMAND lsdtt-test-secret ${CMAKE_CURRENT_BINARY_DIR}/secret_test_output)
set_tests_properties(unit_secret PROPERTIES LABELS "unit;fetch")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/secret_test_output)

# Native DEM fetching (docs/testing_and_benchmarking.md).
#
# The offline half needs no network and runs in the default suite: it covers UTM
# zone derivation, tile naming, option validation, and -- most importantly --
# that an unobtainable source degrades to "skipped" rather than failing.
add_executable(lsdtt-test-raster-fetch cpp/test_raster_fetch.cpp)
target_link_libraries(lsdtt-test-raster-fetch PRIVATE lsdtt3_core)
target_include_directories(lsdtt-test-raster-fetch PRIVATE ${CMAKE_SOURCE_DIR}/src/core)

add_test(NAME unit_raster_fetch
         COMMAND lsdtt-test-raster-fetch --offline --data-dir=${CMAKE_CURRENT_SOURCE_DIR}/data/bounds)
set_tests_properties(unit_raster_fetch PROPERTIES LABELS "unit;fetch")

# The network half actually downloads. Off by default; skips (77) rather than
# failing when the source is unreachable.
option(LSDTT_ENABLE_NETWORK_TESTS "Enable tests that download data" OFF)

if(LSDTT_ENABLE_NETWORK_TESTS)
    add_test(NAME network_raster_fetch
             COMMAND lsdtt-test-raster-fetch --network
                     ${CMAKE_CURRENT_BINARY_DIR}/fetch_test_output)
    set_tests_properties(network_raster_fetch
        PROPERTIES LABELS "network;fetch" SKIP_RETURN_CODE 77)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fetch_test_output)
endif()

add_test(NAME unit_cosmo_physics
         COMMAND lsdtt-test-cosmo-physics
                 ${CMAKE_CURRENT_BINARY_DIR}/cosmo_test_output
                 ${CMAKE_CURRENT_SOURCE_DIR}/expected/cosmo/numerics_sentinel.csv)
set_tests_properties(unit_cosmo_physics PROPERTIES LABELS "unit;cosmo")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cosmo_test_output)
