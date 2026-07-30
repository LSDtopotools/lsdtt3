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
    smoke_hydraulic_routing
    lsdtt-hydraulic-routing
    hydraulic_routing_smoke.param.in
    tiny_hydraulic_water_depth.tif
    --expect tiny_hydraulic_water_surface.tif
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
    failure_hydraulic_routing_invalid_mode
    lsdtt-hydraulic-routing
    hydraulic_routing_invalid_mode.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: hydraulic_run_mode must be 'stationary' or 'transient'."
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
                PROPERTIES LABELS "benchmark;surface_metrics" SKIP_RETURN_CODE 77
                           TIMEOUT 1800)
        endforeach()
    endforeach()
    endif()

    # Hydraulic DAGGER parity benchmark (plan Phases 3-4).
    # In --check mode it enforces mask agreement and pinned parity tolerances.
    # Prefer launching through pixi's benchmark env when available so GDAL/PROJ
    # and Python deps match the documented benchmark path.
    set(LSDTT_DAGGER_ROOT_DEFAULT "${CMAKE_SOURCE_DIR}/../DAGGER")
    set(LSDTT_DAGGER_REF_DEFAULT
        "${CMAKE_CURRENT_SOURCE_DIR}/scripts/dagger_graphflood_reference.py")

    if(LSDTT_PIXI_EXECUTABLE)
        set(LSDTT_HYDRAULIC_BENCHMARK_LAUNCHER ${LSDTT_BENCHMARK_LAUNCHER})
    else()
        set(LSDTT_HYDRAULIC_BENCHMARK_LAUNCHER ${Python3_EXECUTABLE})
    endif()

    foreach(hyd_fixture tiny_dem.asc alpine_small)
        foreach(hyd_scenario stationary_fixed transient_dynamic)
            add_test(
                NAME benchmark_hydraulic_dagger_${hyd_scenario}_${hyd_fixture}
                COMMAND ${LSDTT_HYDRAULIC_BENCHMARK_LAUNCHER}
                        "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_hydraulic_dagger_benchmark.py"
                        --exe "$<TARGET_FILE:lsdtt-hydraulic-routing>"
                        --fetch-exe "$<TARGET_FILE:lsdtt-fetch-raster>"
                        --dagger-root "${LSDTT_DAGGER_ROOT_DEFAULT}"
                        --dagger-ref "${LSDTT_DAGGER_REF_DEFAULT}"
                        --fixture ${hyd_fixture}
                        --scenario ${hyd_scenario}
                        --work-dir "${CMAKE_CURRENT_BINARY_DIR}/hydraulic_benchmark/${hyd_scenario}_${hyd_fixture}"
                        --write-json
                        --check
            )
            set_tests_properties(
                benchmark_hydraulic_dagger_${hyd_scenario}_${hyd_fixture}
                PROPERTIES LABELS "benchmark;hydraulic" SKIP_RETURN_CODE 77
                           TIMEOUT 1800)
        endforeach()
    endforeach()
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
set_tests_properties(unit_surface_metrics PROPERTIES LABELS "unit;surface_metrics")

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

add_executable(lsdtt-test-hydraulic-routing cpp/test_hydraulic_routing.cpp)
target_link_libraries(lsdtt-test-hydraulic-routing PRIVATE lsdtt3_cpu)

add_test(NAME unit_hydraulic_routing
         COMMAND lsdtt-test-hydraulic-routing)
set_tests_properties(unit_hydraulic_routing PROPERTIES LABELS "unit;hydraulic")

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

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_bad_bounds_mode
    lsdtt-fetch-raster
    fetch_raster_bad_bounds_mode.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: bounds_mode must be one of"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_outlets_missing_file
    lsdtt-fetch-raster
    fetch_raster_outlets_missing_file.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: bounds_mode=outlets requires outlet_points_fname"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_outlets_conflict_area
    lsdtt-fetch-raster
    fetch_raster_outlets_conflict_area.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: bounds_mode=outlets forbids west/south/east/north and bounds_fname"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_bad_prefetch_source
    lsdtt-fetch-raster
    fetch_raster_bad_prefetch_source.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: unknown prefetch_source"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_bad_auth_fallback_source
    lsdtt-fetch-raster
    fetch_raster_bad_auth_fallback_source.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: auth_fallback_source must be keyless"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_outlets_missing_dataset
    lsdtt-fetch-raster
    fetch_raster_outlets_not_implemented.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: could not open outlet points CSV"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_ot_missing_key_strict
    lsdtt-fetch-raster
    fetch_raster_ot_missing_key_strict.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: no credential found"
)

add_lsdtt_cli_smoke_test(
    failure_fetch_raster_ot_missing_key_fallback
    lsdtt-fetch-raster
    fetch_raster_ot_missing_key_fallback.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: lsdtt::fetch_dem: west must be less than east"
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

    add_lsdtt_cli_smoke_test(
        network_fetch_raster_outlets_cli
        lsdtt-fetch-raster
        fetch_raster_outlets_network.param.in
        outlets_net_group_0001_DEM.tif
        --expect outlets_net_group_0002_DEM.tif
        --expect outlets_net_dem_groups.csv
        --expect outlets_net_outlet_bounds_status.csv
        --expect-file-contains "outlets_net_dem_groups.csv:group_id"
        --expect-file-contains "outlets_net_dem_groups.csv:group_0002_DEM.tif"
    )
    set_tests_properties(network_fetch_raster_outlets_cli PROPERTIES LABELS "network;cli;fetch")
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
            PROPERTIES LABELS "performance;surface_metrics" SKIP_RETURN_CODE 77
                       TIMEOUT 3600)
    endif()
    # Per-pixel LSDn production cost (LSDn plan Sec3.8). Four cases: Stone vs
    # LSDn at age 0; LSDn age 0 vs 20 ka vs 1 Ma (1 Ma pulls in the long record);
    # the time-averaged surface vs the direct per-pixel loop, both driven
    # explicitly, with the measured crossover pixel count; and the PRE-Phase-8
    # stone2000sp_paleomag baseline. Emits lsdn_performance.csv.
    #
    # REPORTED, NEVER ASSERTED -- it can only fail by not running. That the
    # surface and the direct path AGREE is asserted in L10 (unit_cosmo_physics),
    # which is where a correctness claim belongs.
    add_test(NAME perf_cosmo_lsdn_production
             COMMAND lsdtt-test-cosmo-physics
                     ${CMAKE_CURRENT_BINARY_DIR}/perf_output --performance)
    set_tests_properties(perf_cosmo_lsdn_production
        PROPERTIES LABELS "performance;cosmo" TIMEOUT 1800)

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

add_executable(lsdtt-test-dem-grouping cpp/test_dem_grouping.cpp)
target_link_libraries(lsdtt-test-dem-grouping PRIVATE lsdtt3_core)
target_include_directories(lsdtt-test-dem-grouping PRIVATE ${CMAKE_SOURCE_DIR}/src/core)

add_test(NAME unit_dem_grouping
         COMMAND lsdtt-test-dem-grouping)
set_tests_properties(unit_dem_grouping PROPERTIES LABELS "unit;fetch")

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

# Arguments: output dir, Stone sentinel, LSDn sentinel (a SEPARATE file with its
# own gated regeneration flag; plan Sec3b), and the tests/data fixture dir (the
# pre-LSDn crn_parameters.nc used by L5).
add_test(NAME unit_cosmo_physics
         COMMAND lsdtt-test-cosmo-physics
                 ${CMAKE_CURRENT_BINARY_DIR}/cosmo_test_output
                 ${CMAKE_CURRENT_SOURCE_DIR}/expected/cosmo/numerics_sentinel.csv
                 ${CMAKE_CURRENT_SOURCE_DIR}/expected/cosmo/numerics_sentinel_lsdn.csv
                 ${CMAKE_CURRENT_SOURCE_DIR}/data)
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

# --- LSDn scaling scheme through the CLI (LSDn plan L8) ----------------------
#
# tiny_dem.asc is a 5x5 tilted plane in EPSG:32630 (UTM 30N) with its origin at
# easting/northing 0, which places it near the equator at about 7.489 deg W;
# cosmo_tiny_points.csv holds one sample at the WGS84 position of the (4.5, 0.5)
# cell centre, i.e. the low corner the plane drains to. threshold_contributing_
# pixels is dropped to 2 so a 25-pixel DEM has a channel to snap to.
#
# The run must SUCCEED, write _PROD.tif with every pixel valid and positive, and
# report the scheme in its run header. The header assertion is the point: a
# silent fall back to Stone would otherwise pass every other check here.
add_lsdtt_cli_smoke_test(
    smoke_cosmo_tool_lsdn
    lsdtt-cosmo-tool
    cosmo_tool_lsdn_smoke.param.in
    cosmo_lsdn_PROD.tif
    --expect cosmo_lsdn_CN_ERATE.csv
    --expect-raster-min cosmo_lsdn_PROD.tif:0.0001
    --stdout-contains "Scaling scheme: LSDn"
    --stdout-contains "LSDn scaling (Be10)"
)
set_tests_properties(smoke_cosmo_tool_lsdn
    PROPERTIES LABELS "smoke;cli;cosmo")

# self_consistent_age: the sample age and the erosion rate are solved together,
# t = Gamma_sp/eps.  The assertions are chosen so the test cannot pass if the
# iteration silently degenerates into the fixed-age path:
#   - "self-consistent iter 2" proves it went round MORE THAN ONCE.  A one-shot
#     loop, or one that converged on the first comparison, would not print it.
#   - the CSV must carry the three diagnostic columns, so the converged age is
#     recorded per sample rather than only narrated on stdout.
# Verified to fail on perturbation: forcing the loop to a single pass drops the
# "iter 2" line, and dropping the CSV columns fails the file-contains checks.
add_lsdtt_cli_smoke_test(
    smoke_cosmo_tool_self_consistent_age
    lsdtt-cosmo-tool
    cosmo_tool_self_consistent_age.param.in
    cosmo_sca_CN_ERATE.csv
    --stdout-contains "self-consistent iter 1"
    --stdout-contains "self-consistent iter 2"
    --stdout-contains "self-consistent: t = Gamma_sp/eps per sample"
    --expect-file-contains "cosmo_sca_CN_ERATE.csv:self_consistent_t_yr"
    --expect-file-contains "cosmo_sca_CN_ERATE.csv:self_consistent_n_iter"
    --expect-file-contains "cosmo_sca_CN_ERATE.csv:self_consistent_converged"
)
set_tests_properties(smoke_cosmo_tool_self_consistent_age
    PROPERTIES LABELS "smoke;cli;cosmo")

# A sample that runs out of iterations must be REPORTED, not quietly accepted.
# max_iter 1 cannot converge from the modern field (the first step moves the rate
# ~7 %), so this pins the warning text and the converged=0 flag in the CSV.  The
# run still SUCCEEDS -- non-convergence is a caveat on a number, not a failure.
add_lsdtt_cli_smoke_test(
    smoke_cosmo_tool_self_consistent_no_converge
    lsdtt-cosmo-tool
    cosmo_tool_self_consistent_no_converge.param.in
    cosmo_sca_nc_CN_ERATE.csv
    --stderr-contains "did not converge in 1 self_consistent_age iterations"
    --stderr-contains "Treat it as unconverged, not as a result"
    --expect-file-contains "cosmo_sca_nc_CN_ERATE.csv:,1,0"
)
set_tests_properties(smoke_cosmo_tool_self_consistent_no_converge
    PROPERTIES LABELS "smoke;cli;cosmo")

# A non-positive tolerance would make the loop exit on its first comparison and
# report a "converged" age that is nothing of the kind, so it is FATAL up front.
add_lsdtt_cli_smoke_test(
    failure_cosmo_tool_bad_self_consistent_tol
    lsdtt-cosmo-tool
    cosmo_tool_bad_self_consistent_tol.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "self_consistent_age_tol must be > 0"
)
set_tests_properties(failure_cosmo_tool_bad_self_consistent_tol
    PROPERTIES LABELS "failure;cli;cosmo")

# The forward -> inverse CLI BRIDGE (cosmo CLI test plan Phase 4), plus knob
# coverage. lsdtt-cosmo-prediction forward-models a concentration for a known
# erosion rate and writes the points CSV it advertises as ingestible;
# lsdtt-cosmo-tool inverts that CSV; the imposed rate must come back. This is
# the only test that drives the two CLIs against each other, so it is the only
# one that can catch plumbing faults -- a column renamed on one side, a knob
# honoured asymmetrically, a shielding term applied twice -- all of which
# produce plausible numbers and pass every single-CLI test.
#
# Offline and fast (a 5x5 fixture DEM), so it runs in the default suite.
# Verified to fail on perturbation: ignoring the standardisation factor, swapping
# the +/-1 sigma bounds, and dropping topographic shielding on the inverse side
# each turn it red.
add_test(NAME rigorous_cosmo_cli_bridge
         COMMAND ${Python3_EXECUTABLE}
                 "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cosmo_bridge_test.py"
                 --pred-exe "$<TARGET_FILE:lsdtt-cosmo-prediction>"
                 --tool-exe "$<TARGET_FILE:lsdtt-cosmo-tool>"
                 --data-dir "${LSDTT_TEST_DATA_DIR}")
set_tests_properties(rigorous_cosmo_cli_bridge PROPERTIES
    LABELS "rigorous;cli;cosmo" SKIP_RETURN_CODE 77 TIMEOUT 900)

# Multi-isotope step change: Be-10, Al-26 and C-14 through ONE scenario, plus
# the two ratio rasters. The isotope RATIOS are what carry a transient signal,
# and a ratio is only meaningful if both members saw an identical scenario --
# hence one particle evolution per pixel, not three separate runs.
#
# Bitwise equality of each nuclide against the single-nuclide kernel is asserted
# in the unit tier (M1); this pins the CLI plumbing: all five rasters appear,
# every pixel valid, and the ratio rasters carry positive, physically sane
# values. Verified to fail on perturbation.
add_lsdtt_cli_smoke_test(
    smoke_cosmo_prediction_multi_isotope
    lsdtt-cosmo-prediction
    cosmo_prediction_multi_isotope.param.in
    mi_Be10_step_crn_conc.tif
    --expect mi_Al26_step_crn_conc.tif
    --expect mi_C14_step_crn_conc.tif
    --expect mi_R2610_step.tif
    --expect mi_R1410_step.tif
    --expect-valid-percent-100 mi_R2610_step.tif
    --expect-valid-percent-100 mi_R1410_step.tif
    --expect-raster-min mi_R2610_step.tif:1.0
    --expect-raster-min mi_R1410_step.tif:0.01
    --stdout-contains "multi-isotope step-change rasters (Be10, Al26, C14)"
    --stdout-contains "Al26/Be10 ratio raster"
    --stdout-contains "C14/Be10 ratio raster"
)
set_tests_properties(smoke_cosmo_prediction_multi_isotope
    PROPERTIES LABELS "smoke;cli;cosmo")

# Multi-isotope step change under LSDn scaling: smoke that the prediction CLI
# reaches the LSDn path and writes the same 3 nuclide + 2 ratio rasters.
add_lsdtt_cli_smoke_test(
    smoke_cosmo_prediction_multi_isotope_lsdn
    lsdtt-cosmo-prediction
    cosmo_prediction_multi_isotope_lsdn.param.in
    mi_lsdn_Be10_step_crn_conc.tif
    --expect mi_lsdn_Al26_step_crn_conc.tif
    --expect mi_lsdn_C14_step_crn_conc.tif
    --expect mi_lsdn_R2610_step.tif
    --expect mi_lsdn_R1410_step.tif
    --expect-valid-percent-100 mi_lsdn_R2610_step.tif
    --expect-valid-percent-100 mi_lsdn_R1410_step.tif
    --expect-raster-min mi_lsdn_R2610_step.tif:1.0
    --expect-raster-min mi_lsdn_R1410_step.tif:0.01
    --stdout-contains "Scaling scheme: LSDn"
    --stdout-contains "LSDn scaling (Be10)"
)
set_tests_properties(smoke_cosmo_prediction_multi_isotope_lsdn
    PROPERTIES LABELS "smoke;cli;cosmo")

# Multi-isotope ARBITRARY erosion history: the general case of the step-change
# mode. Bitwise equality with the step-change kernel for a two-row history is
# asserted in the unit tier (H1); this pins the CLI plumbing -- the history file
# is read and validated, all five rasters appear, and the ratios are sane.
add_lsdtt_cli_smoke_test(
    smoke_cosmo_prediction_multi_transient
    lsdtt-cosmo-prediction
    cosmo_prediction_multi_transient.param.in
    tr_Be10_step_crn_conc.tif
    --expect tr_Al26_step_crn_conc.tif
    --expect tr_C14_step_crn_conc.tif
    --expect tr_R2610_step.tif
    --expect tr_R1410_step.tif
    --expect-valid-percent-100 tr_R2610_step.tif
    --expect-raster-min tr_R2610_step.tif:1.0
    --stdout-contains "multi-isotope TRANSIENT rasters (Be10, Al26, C14)"
    --stdout-contains "every pixel walks the whole history"
)
set_tests_properties(smoke_cosmo_prediction_multi_transient
    PROPERTIES LABELS "smoke;cli;cosmo")

# An unknown scheme is FATAL, and the message names both allowed values.
add_lsdtt_cli_smoke_test(
    failure_cosmo_tool_bad_scaling_scheme
    lsdtt-cosmo-tool
    cosmo_tool_bad_scaling_scheme.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "unknown scaling_scheme"
    --stderr-contains "stone"
    --stderr-contains "lsdn"
)
set_tests_properties(failure_cosmo_tool_bad_scaling_scheme
    PROPERTIES LABELS "failure;cli;cosmo")

# scaling_scheme: lsdn against an nc with NO LSDn grids must fail hard and name
# the nc it searched -- never silently degrade to Stone. crn_parameters_pre_lsdn.nc
# is the tracked pre-LSDn nc, which is exactly that situation. Library-level
# counterpart: L5(e) in unit_cosmo_physics.
add_lsdtt_cli_smoke_test(
    failure_cosmo_tool_lsdn_without_grids
    lsdtt-cosmo-tool
    cosmo_tool_lsdn_no_grids.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "requires the LSDn scaling grids"
    --stderr-contains "crn_parameters_pre_lsdn.nc"
)
set_tests_properties(failure_cosmo_tool_lsdn_without_grids
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

# Functional coverage for lsdtt-cosmo-column-tool, which until 2026-07-27 had
# only the -h test below. It is the natural place to compare lsdtt3 against
# POINT calculators (CRONUS v3, CRONUScalc): no DEM, no flow routing, no basin
# averaging, so nothing can mask a scaling or muon discrepancy.
#
# The checks are LIMIT IDENTITIES of the column model, needing no external
# reference: step_change at t=0 must equal steady state at the old rate and at
# t->inf the new rate; a constant-rate full_transient must relax to steady
# state; monotonicity in depth and rate; production_scaling exactly linear in P0
# while total_shielding is not (Vermeesch redistribution); C-14 relaxing faster
# than Be-10. Verified to fail on perturbation.
add_test(NAME rigorous_cosmo_column_tool
         COMMAND ${Python3_EXECUTABLE}
                 "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cosmo_column_test.py"
                 --exe "$<TARGET_FILE:lsdtt-cosmo-column-tool>")
set_tests_properties(rigorous_cosmo_column_tool PROPERTIES
    LABELS "rigorous;cli;cosmo" SKIP_RETURN_CODE 77 TIMEOUT 600)

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

# Riversand-comparison harness self-test (Phase 5 carried-forward item §6c(B)).
# OFFLINE scaffold check only: manifest/sample preflight + cleanup-trap contract.
add_test(NAME unit_cosmo_riversand_harness
         COMMAND ${Python3_EXECUTABLE}
                 "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cosmo_riversand_comparison.py"
                 --self-test)
set_tests_properties(unit_cosmo_riversand_harness PROPERTIES LABELS "unit;cosmo")

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

# Additional real-landscape datasets (harvested with tests/scripts/harvest_octopus.py),
# chosen to CONTRAST with the moderate-relief Front Range pilot across terrain and
# elevation: clapp2002 = Yuma Wash arid LOW-relief / LOW-elevation (OCTOPUS S014);
# dibiase2010 = San Gabriel STEEP / HIGH-relief (OCTOPUS S055). Each is fully
# manifest-driven (bbox, target UTM, run knobs) and frozen from its own accepted
# networked run. Same --check regression mode, same skip-77 / out-of-CI wiring.
add_test(NAME cosmo_realworld_clapp2002
         COMMAND ${Python3_EXECUTABLE}
                 "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cosmo_validation.py"
                 clapp2002 --check
                 --exe "$<TARGET_FILE:lsdtt-cosmo-tool>"
                 --fetch-exe "$<TARGET_FILE:lsdtt-fetch-raster>")
set_tests_properties(cosmo_realworld_clapp2002 PROPERTIES
    LABELS "cosmo_realworld" SKIP_RETURN_CODE 77 TIMEOUT 1800)

add_test(NAME cosmo_realworld_dibiase2010
         COMMAND ${Python3_EXECUTABLE}
                 "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cosmo_validation.py"
                 dibiase2010 --check
                 --exe "$<TARGET_FILE:lsdtt-cosmo-tool>"
                 --fetch-exe "$<TARGET_FILE:lsdtt-fetch-raster>")
set_tests_properties(cosmo_realworld_dibiase2010 PROPERTIES
    LABELS "cosmo_realworld" SKIP_RETURN_CODE 77 TIMEOUT 1800)

# LSDn tier of the same three datasets (lsdn plan Phase 6 deliverable (iii)): the
# END-TO-END PER-PIXEL LSDn run, (a) per-pixel LSDn basin mean vs (c) OCTOPUS EBEL.
# Each <name>_lsdn manifest REUSES its St counterpart's bbox, run knobs and SAMPLE
# TABLE -- there is exactly one copy of every measured concentration in the repo --
# and changes only scaling_scheme (lsdn) and the published target (EBE_* -> EBEL,
# via published_erate_override reading the tracked octopus_lsdn_point_reference.csv).
# The frozen band is a MONITORING BAND, NOT AN ORACLE: it mixes the kernel residual
# (pinned DEM-free by L12), the per-pixel-vs-lumped difference (pinned by L13, and a
# result rather than an error) and the DEM/delineation difference, and it runs at the
# CLI default sample_age_yr = 0 while EBEL is time-averaged over Gamma/eps. Attribute
# it with L12/L13 (see each summary header); never tune anything to shrink it.
# Same --check regression mode, same skip-77 / out-of-CI wiring as the St tier: these
# compare the fresh lsdtt3 rates against the frozen rows and must NOT rewrite them.
add_test(NAME cosmo_realworld_dethier2014_lsdn
         COMMAND ${Python3_EXECUTABLE}
                 "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cosmo_validation.py"
                 dethier2014_lsdn --check
                 --exe "$<TARGET_FILE:lsdtt-cosmo-tool>"
                 --fetch-exe "$<TARGET_FILE:lsdtt-fetch-raster>")
set_tests_properties(cosmo_realworld_dethier2014_lsdn PROPERTIES
    LABELS "cosmo_realworld" SKIP_RETURN_CODE 77 TIMEOUT 1800)

add_test(NAME cosmo_realworld_clapp2002_lsdn
         COMMAND ${Python3_EXECUTABLE}
                 "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cosmo_validation.py"
                 clapp2002_lsdn --check
                 --exe "$<TARGET_FILE:lsdtt-cosmo-tool>"
                 --fetch-exe "$<TARGET_FILE:lsdtt-fetch-raster>")
set_tests_properties(cosmo_realworld_clapp2002_lsdn PROPERTIES
    LABELS "cosmo_realworld" SKIP_RETURN_CODE 77 TIMEOUT 1800)

add_test(NAME cosmo_realworld_dibiase2010_lsdn
         COMMAND ${Python3_EXECUTABLE}
                 "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cosmo_validation.py"
                 dibiase2010_lsdn --check
                 --exe "$<TARGET_FILE:lsdtt-cosmo-tool>"
                 --fetch-exe "$<TARGET_FILE:lsdtt-fetch-raster>")
set_tests_properties(cosmo_realworld_dibiase2010_lsdn PROPERTIES
    LABELS "cosmo_realworld" SKIP_RETURN_CODE 77 TIMEOUT 1800)

# §6e production-vs-BEPROD isolation (decisive test for the ~+4% Dethier residual).
# At fixed concentration inferred erosion is proportional to basin-mean production;
# this fetches the keyless COP30 DEM, runs lsdtt-cosmo-tool MATCHED to OCTOPUS
# (Braucher / NCEP / paleomag-off) with the reporting-only basin-production CSV, and
# compares lsdtt's basin-mean Stone-2000 production scaling to OCTOPUS BEPROD/BETOTS
# and its basin-mean elevation to ELEV_AVE (harvested live from the OCTOPUS WFS).
# NEEDS NETWORK; kept OUT of the default unit;cosmo suite (label cosmo_realworld,
# `ctest -L cosmo_realworld`). Skips (exit 77) when the network is unavailable.
# Runs in check-only mode: it compares against the tracked frozen artefact and must
# NOT rewrite it. Re-freezing is a deliberate act -- run the script with --freeze by
# hand and review the diff.
add_test(NAME cosmo_realworld_production_isolation
         COMMAND ${Python3_EXECUTABLE}
                 "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cosmo_production_isolation.py"
                 --exe "$<TARGET_FILE:lsdtt-cosmo-tool>"
                 --fetch-exe "$<TARGET_FILE:lsdtt-fetch-raster>")
set_tests_properties(cosmo_realworld_production_isolation PROPERTIES
    LABELS "cosmo_realworld" SKIP_RETURN_CODE 77 TIMEOUT 1800)

# L13 GENERATION STEP -- the PER-PIXEL side of per-pixel-vs-lumped LSDn (lsdn plan
# Phase 6 deliverable (ii); cosmo plan §6f row L13). Fetches each of the three
# tracked datasets' keyless COP30 DEMs and records lsdtt-cosmo-tool's basin-mean
# production scaling (Stone, and LSDn at three sample ages) over the exact upslope
# node set the inversion uses.
#
# Deliberately a TWO-STEP design: this expensive networked step is frozen into
# tests/expected/cosmo/lsdn_pixel_vs_lumped_perpixel.csv, and the CHEAP comparison
# against the lumped LSDn kernel is then regression-checked on every build inside
# unit_cosmo_physics (default unit;cosmo suite), which is where the assertions and
# the lsdn_pixel_vs_lumped.csv product live.
#
# NEEDS NETWORK; kept OUT of the default suite (label cosmo_realworld). Skips
# (exit 77) when the DEM source is unreachable. Runs check-only: it compares
# against the tracked frozen artefact and must NOT rewrite it. Re-freezing is a
# deliberate act -- run the script with --freeze by hand and review the diff.
add_test(NAME cosmo_realworld_lsdn_pixel_vs_lumped
         COMMAND ${Python3_EXECUTABLE}
                 "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cosmo_lsdn_pixel_vs_lumped.py"
                 all
                 --exe "$<TARGET_FILE:lsdtt-cosmo-tool>"
                 --fetch-exe "$<TARGET_FILE:lsdtt-fetch-raster>")
set_tests_properties(cosmo_realworld_lsdn_pixel_vs_lumped PROPERTIES
    LABELS "cosmo_realworld" SKIP_RETURN_CODE 77 TIMEOUT 3600)
