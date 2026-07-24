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

    # Surface-metrics SPEED benchmark vs TopoToolbox (surface-metrics plan,
    # Phases 3-4).  TopoToolbox is a PERFORMANCE reference here, not a
    # correctness one -- correctness comes from the Tier A analytic oracle in
    # lsdtt-test-surface-metrics, which needs no network.  So --check enforces
    # ONLY the rank-correlation guard: it fails when the two tools turn out not
    # to be computing the same kind of quantity, which would make the timing
    # meaningless.  No wall-clock threshold is ever asserted.
    #
    # The _medium fixtures, not the _small ones the flow benchmark uses: at
    # 256x256 the polyfit is ~5 ms and every timing is noise.
    foreach(surf_fixture alpine_medium lakeland_medium)
        foreach(surf_case slope curvature_all4 hillshade)
            add_test(
                NAME benchmark_surface_${surf_case}_${surf_fixture}
                COMMAND ${LSDTT_BENCHMARK_LAUNCHER}
                        "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_surface_benchmark.py"
                        --exe "$<TARGET_FILE:lsdtt-surface-metrics>"
                        --fixture ${surf_fixture}
                        --comparison ${surf_case}
                        --repeats 5
                        --skip-studies
                        --work-dir "${CMAKE_CURRENT_BINARY_DIR}/surface_benchmark"
                        --check
            )
            set_tests_properties(benchmark_surface_${surf_case}_${surf_fixture}
                PROPERTIES LABELS "benchmark;surface" SKIP_RETURN_CODE 77
                           TIMEOUT 1800)
        endforeach()
    endforeach()
    endif()
endif()

# -----------------------------------------------------------------------------
# C++ unit tests
# -----------------------------------------------------------------------------

# Surface-metrics validation: Tier A analytic oracle (global quadratic and
# periodic-wave surfaces with exact Taylor coefficients, sphere, cylinder,
# plane) plus Tier C property tests (algebraic identities, SMAP classification,
# rotation invariance, thread determinism, MADk2 invariants).  Synthetic
# fixtures only: no network and no committed DEM, so it runs in the default
# suite.
#
# NOTE: two MADk2 checks are KNOWN FAILURES (the anisotropy/DMC estimator treats
# axial directions as vectors).  They are reported as KNOWN-FAIL and do not fail
# the run; flip k_madk2_axial_is_fixed in the driver once the estimator is
# corrected.
add_executable(lsdtt-test-surface-metrics cpp/test_surface_metrics.cpp)
target_link_libraries(lsdtt-test-surface-metrics PRIVATE lsdtt3_cpu)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/surface_metrics_output)
add_test(NAME unit_surface_metrics
         COMMAND lsdtt-test-surface-metrics
                 ${CMAKE_CURRENT_BINARY_DIR}/surface_metrics_output)
set_tests_properties(unit_surface_metrics PROPERTIES LABELS "unit;surface")

add_executable(lsdtt-test-param-parser-help cpp/test_param_parser_help.cpp)
target_link_libraries(lsdtt-test-param-parser-help PRIVATE lsdtt3_core)

add_test(NAME unit_param_parser_help COMMAND lsdtt-test-param-parser-help)
set_tests_properties(unit_param_parser_help PROPERTIES LABELS "unit")

# Cosmogenic physics driver (docs/testing_and_benchmarking.md).  Writes its
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

# lsdtt-surface-metrics.  Uses tiny_dem.asc (UTM 30N) rather than tiny30_dem.tif,
# which is EPSG:4326 and is rejected by the CLI's projected-CRS check - that
# rejection is itself covered by failure_surface_metrics_geographic_crs below.
add_lsdtt_cli_smoke_test(
    smoke_surface_metrics
    lsdtt-surface-metrics
    surface_metrics_smoke.param.in
    tiny_hillshade.tif
    --expect tiny_local_relief.tif
    --expect tiny_direct_polyfit_slope.tif
    --expect tiny_direct_polyfit_aspect_degrees.tif
    --expect tiny_direct_polyfit_mean_curvature.tif
    --expect tiny_direct_polyfit_maximum_curvature.tif
    --expect tiny_direct_polyfit_minimum_curvature.tif
    --expect tiny_direct_shape_class.tif
    --expect tiny_direct_shape_gaussian_curvature.tif
    --expect tiny_direct_shape_mean_curvature.tif
)

add_lsdtt_cli_smoke_test(
    smoke_surface_metrics_madk2
    lsdtt-surface-metrics
    surface_metrics_madk2_smoke.param.in
    tiny_madk2_or.tif
    --expect tiny_madk2_dmc.tif
    --expect tiny_madk2_ra.tif
)

add_lsdtt_cli_smoke_test(
    failure_surface_metrics_invalid_metric
    lsdtt-surface-metrics
    surface_metrics_invalid_metric.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "Unknown polyfit metric: not_a_metric"
)

add_lsdtt_cli_smoke_test(
    failure_surface_metrics_geographic_crs
    lsdtt-surface-metrics
    surface_metrics_geographic_crs.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "geographic coordinate system"
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

    # Surface-metrics timing sweep: every comparison, plus the
    # marginal-cost-per-extra-metric and window-size (O(n)) studies.
    # Reported, never asserted -- no --check, so a slow shared runner cannot
    # turn it red.  Needs pixi + topotoolbox + a fetched fixture; skips (77)
    # when any of the three is missing.
    if(LSDTT_ENABLE_BENCHMARKS AND LSDTT_PIXI_EXECUTABLE)
        add_test(
            NAME performance_surface_metrics
            COMMAND ${LSDTT_BENCHMARK_LAUNCHER}
                    "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_surface_benchmark.py"
                    --exe "$<TARGET_FILE:lsdtt-surface-metrics>"
                    --fixture alpine_medium
                    --comparison all
                    --repeats 5
                    --work-dir "${CMAKE_CURRENT_BINARY_DIR}/surface_benchmark"
        )
        set_tests_properties(performance_surface_metrics
            PROPERTIES LABELS "performance;surface" SKIP_RETURN_CODE 77
                       TIMEOUT 3600)
    endif()
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

# -----------------------------------------------------------------------------
# Cosmo CLI integration tests (plan Phase 4; docs/testing_and_benchmarking.md)
#
# Offline: they exercise parameter-file parsing, help rendering, and the
# projected-CRS guard on the three cosmo CLIs. The help tests also drive the
# bundled-crn_parameters.nc auto-discovery path (src/core/lsdtt_exe_path.hpp)
# indirectly; the direct auto-discovery check lives in unit_cosmo_physics.
# -----------------------------------------------------------------------------

# Robustness: a geographic (non-projected) DEM must be a clean FATAL, not a
# crash or a silent wrong-units run. tiny30_dem.tif is EPSG:4326.
add_lsdtt_cli_smoke_test(
    failure_cosmo_tool_geographic_crs
    lsdtt-cosmo-tool
    cosmo_tool_geographic_crs.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "DEM must be in a projected CRS"
)
set_tests_properties(failure_cosmo_tool_geographic_crs
    PROPERTIES LABELS "failure;cli;cosmo")

# Help must render every parameter without a parameter file and without
# touching the DEM or the network. Guards that the parameter declarations for
# each cosmo CLI stay well-formed.
add_test(NAME smoke_cosmo_tool_help COMMAND lsdtt-cosmo-tool -h)
set_tests_properties(smoke_cosmo_tool_help PROPERTIES
    LABELS "smoke;cli;cosmo"
    PASS_REGULAR_EXPRESSION "sample_age_yr")

add_test(NAME smoke_cosmo_prediction_help COMMAND lsdtt-cosmo-prediction -h)
set_tests_properties(smoke_cosmo_prediction_help PROPERTIES
    LABELS "smoke;cli;cosmo"
    PASS_REGULAR_EXPRESSION "sample_age_yr")

add_test(NAME smoke_cosmo_column_tool_help COMMAND lsdtt-cosmo-column-tool -h)
set_tests_properties(smoke_cosmo_column_tool_help PROPERTIES
    LABELS "smoke;cli;cosmo"
    PASS_REGULAR_EXPRESSION "nuclide")

# -----------------------------------------------------------------------------
# Real-landscape cosmogenic validation tier (plan Phase 6 / §5b;
# tests/expected/cosmo/published/README.md)
#
# "download -> run -> delete": run_cosmo_validation.py fetches a DEM KEYLESS via
# lsdtt-fetch-raster (Copernicus GLO-30 / Planetary Computer -- no API key),
# clips it, runs lsdtt-cosmo-tool, compares to the OCTOPUS-harvested published
# rates in the tracked sample table, then deletes the temporary working directory
# INCLUDING ON FAILURE (cleanup trap). Without --check a run FREEZES the summary
# (append); the ctest test below uses --check to REGRESS a fresh run against the
# frozen summary WITHOUT modifying the tracked file.
# -----------------------------------------------------------------------------

# Harness self-test: OFFLINE. Exercises the cleanup-trap contract (the trap fires
# and removes the workdir on a simulated failed run), the manifest/sample-table
# schema guards, and the summary-append round-trip -- with a mock workdir, no
# network, no DEM. This is the cheap tracked guard that the harness logic stays
# correct; it runs in the default unit;cosmo suite.
add_test(NAME unit_cosmo_validation_harness
         COMMAND ${Python3_EXECUTABLE}
                 "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cosmo_validation.py"
                 --self-test)
set_tests_properties(unit_cosmo_validation_harness PROPERTIES LABELS "unit;cosmo")

# Real download -> run -> delete for the pilot dataset, in --check (regression)
# mode: it reproduces the fresh inverted rates and compares them to the FROZEN
# tracked summary without modifying it. NEEDS NETWORK and time; deliberately kept
# OUT of the default unit;cosmo CI suite (label cosmo_realworld, invoked manually
# via `ctest -L cosmo_realworld`). Skips (exit 77) when the network/DEM source is
# unavailable OR the manifest/sample table is still a template, so it can never
# turn CI red.
add_test(NAME cosmo_realworld_dethier2014
         COMMAND ${Python3_EXECUTABLE}
                 "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cosmo_validation.py"
                 dethier2014 --check
                 --exe "$<TARGET_FILE:lsdtt-cosmo-tool>"
                 --fetch-exe "$<TARGET_FILE:lsdtt-fetch-raster>")
set_tests_properties(cosmo_realworld_dethier2014 PROPERTIES
    LABELS "cosmo_realworld" SKIP_RETURN_CODE 77 TIMEOUT 1800)
