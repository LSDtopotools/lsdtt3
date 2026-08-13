find_package(Python3 COMPONENTS Interpreter REQUIRED)

set(LSDTT_TEST_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cli_smoke.py")
set(LSDTT_TEST_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data")
set(LSDTT_TEST_PARAM_DIR "${CMAKE_CURRENT_SOURCE_DIR}/params")

# Written by lsdtt-test-divide-migration (the unit_divide_migration CTest
# fixture) and read by the divide-migration CLI tests.
set(LSDTT_DIVIDE_MIGRATION_FIXTURE_DIR
    "${CMAKE_CURRENT_BINARY_DIR}/divide_migration_fixture")

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
    smoke_synthetic_networks_rect
    lsdtt-synthetic-networks
    synthetic_networks_smoke.param.in
    ocnsmoke_ocn_dem.tif
    --expect ocnsmoke_ocn_drainage_cells.tif
    --expect ocnsmoke_ocn_summary.csv
    --expect-raster-size ocnsmoke_ocn_dem.tif:32x32
)

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
    tiny_and_tiny_dem_plus1_difference.tif
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
    --stderr-contains "FATAL: second_raster_read_fname must be set."
)

add_lsdtt_cli_smoke_test(
    failure_divide_migration_missing_outlets
    lsdtt-divide-migration
    divide_migration_missing_outlets.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: basin_outlet_file must be set."
)

add_lsdtt_cli_smoke_test(
    failure_divide_migration_bad_integer_cluster_mode
    lsdtt-divide-migration
    divide_migration_bad_integer_cluster_mode.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: outlet cluster value 'alpha' is not an integer"
)

# --- divide-migration end to end (plan §11.1I) -------------------------------
#
# These run the whole pipeline on the synthetic three-valley DEM written by
# unit_divide_migration, so they need the fixture directory rather than
# tests/data as their read path. run_cli_smoke.py stores --read-path (it does
# not append), so the value passed here through ARGN overrides the default the
# function supplies.
#
# The assertions are on CONTENT, not exit status. The divide-segment count is
# the topology: three parallel valleys give exactly two divides (1|2 and 2|3)
# and never a 1|3, so a run that matched heads across non-adjacent units, or
# that ignored the cluster field, cannot produce these numbers.
add_lsdtt_cli_smoke_test(
    smoke_divide_migration_no_clusters
    lsdtt-divide-migration
    divide_migration_smoke.param.in
    divmig_nocluster_divides.fgb
    --read-path "${LSDTT_DIVIDE_MIGRATION_FIXTURE_DIR}"
    --expect divmig_nocluster_channel_heads.fgb
    --expect divmig_nocluster_divide_migration_metrics.csv
    --expect-vector-feature-min divmig_nocluster_divides.fgb:2
    --expect-vector-feature-min divmig_nocluster_channel_heads.fgb:69
    --expect-file-contains "divmig_nocluster_divide_migration_metrics.csv:head_id,outlet_id,unit_id,pair_id_nearest,head_x,head_y,head_elevation_m,chi,nearest_divide_x,nearest_divide_y,divide_elevation_m,distance_to_divide_m,relief_m,gradient_to_divide"
    --stdout-contains "mapped units     : 3"
    --stdout-contains "selected basins  : 3"
    --stdout-contains "divide segments  : 2"
    --stdout-contains "channel heads    : 69"
    --stdout-contains "heads with divide: 69"
    --stdout-contains "heads without    : 0"
)
set_tests_properties(smoke_divide_migration_no_clusters PROPERTIES
    LABELS "smoke;cli;divide_migration"
    FIXTURES_REQUIRED divide_migration_fixture)

# Grouping outlets 1 and 2 into one unit makes the 1-2 divide internal, so the
# two-divide topology above must collapse to exactly one segment.
add_lsdtt_cli_smoke_test(
    smoke_divide_migration_clusters
    lsdtt-divide-migration
    divide_migration_clusters_smoke.param.in
    divmig_cluster_divides.fgb
    --read-path "${LSDTT_DIVIDE_MIGRATION_FIXTURE_DIR}"
    --expect divmig_cluster_channel_heads.fgb
    --expect divmig_cluster_divide_migration_metrics.csv
    --expect-vector-feature-min divmig_cluster_divides.fgb:1
    --stdout-contains "cluster field    : cluster_id"
    --stdout-contains "cluster groups   : 2"
    --stdout-contains "selected basins  : 3"
    --stdout-contains "divide segments  : 1"
    --stdout-contains "heads with divide: 69"
)
set_tests_properties(smoke_divide_migration_clusters PROPERTIES
    LABELS "smoke;cli;divide_migration"
    FIXTURES_REQUIRED divide_migration_fixture)

# The same clustered run, but the cluster labels arrive in a VECTOR outlet file
# instead of a CSV. Regression test: read_vector_layer("") resolves through
# GetLayerByName(), which never matches an empty string, so this path used to
# fail outright for every vector outlet file. Asserting the same topology as the
# CSV run is the point -- the container must not change the answer.
add_lsdtt_cli_smoke_test(
    smoke_divide_migration_clusters_vector
    lsdtt-divide-migration
    divide_migration_clusters_vector.param.in
    divmig_cluster_vec_divides.fgb
    --read-path "${LSDTT_DIVIDE_MIGRATION_FIXTURE_DIR}"
    --expect divmig_cluster_vec_channel_heads.fgb
    --expect divmig_cluster_vec_divide_migration_metrics.csv
    --expect-vector-feature-min divmig_cluster_vec_divides.fgb:1
    --stdout-contains "cluster field    : cluster_id"
    --stdout-contains "cluster groups   : 2"
    --stdout-contains "selected basins  : 3"
    --stdout-contains "divide segments  : 1"
    --stdout-contains "heads with divide: 69"
)
set_tests_properties(smoke_divide_migration_clusters_vector PROPERTIES
    LABELS "smoke;cli;divide_migration"
    FIXTURES_REQUIRED divide_migration_fixture)

# A clustered run whose outlets all collapse into ONE unit has no unit pair to
# compare. Counting snapped flats rather than units let this through and simply
# produced an empty divide set; the guard now rejects it and names both counts.
add_lsdtt_cli_smoke_test(
    failure_divide_migration_single_unit
    lsdtt-divide-migration
    divide_migration_single_unit_cluster.param.in
    unused_expected_output
    --read-path "${LSDTT_DIVIDE_MIGRATION_FIXTURE_DIR}"
    --expect-failure
    --stderr-contains "FATAL: fewer than two valid snapped outlets in two distinct units remain for unit mapping (outlets=3, units=1)."
    --stderr-contains "Divide extraction needs at least two units to compare."
)
set_tests_properties(failure_divide_migration_single_unit PROPERTIES
    LABELS "failure;cli;divide_migration"
    FIXTURES_REQUIRED divide_migration_fixture)

# --- divide-fragment merging end to end (plan §11.1E) -------------------------
#
# A second fixture DEM: one E-W ridge with a single valley draining north and
# two draining south, so the northern basin faces BOTH southern basins. The
# three-valley fixture above cannot test merging -- every unit pair there is
# adjacent along exactly one stretch of divide.
#
# The pair of tests is the assertion. Unclustered, all three basins are pairwise
# adjacent and the run emits THREE segments. Clustering the two southern outlets
# into one unit leaves that unit bordering the northern unit along TWO fragments
# that meet on the main ridge, and merging must turn those into ONE segment. A
# run that did not merge would report 2 here, and a fixture where the second
# fragment did not exist would fail the unclustered control.
add_lsdtt_cli_smoke_test(
    smoke_divide_migration_merge_unmerged
    lsdtt-divide-migration
    divide_migration_merge_unmerged.param.in
    divmig_merge_nc_divides.fgb
    --read-path "${LSDTT_DIVIDE_MIGRATION_FIXTURE_DIR}"
    --expect-vector-feature-min divmig_merge_nc_divides.fgb:3
    --stdout-contains "mapped units     : 3"
    --stdout-contains "selected basins  : 3"
    --stdout-contains "divide segments  : 3"
)
set_tests_properties(smoke_divide_migration_merge_unmerged PROPERTIES
    LABELS "smoke;cli;divide_migration"
    FIXTURES_REQUIRED divide_migration_fixture)

add_lsdtt_cli_smoke_test(
    smoke_divide_migration_merge_fragments
    lsdtt-divide-migration
    divide_migration_merge_fragments.param.in
    divmig_merge_divides.fgb
    --read-path "${LSDTT_DIVIDE_MIGRATION_FIXTURE_DIR}"
    --expect divmig_merge_divide_migration_metrics.csv
    --expect-vector-feature-min divmig_merge_divides.fgb:1
    --stdout-contains "cluster field    : cluster_id"
    --stdout-contains "cluster groups   : 2"
    --stdout-contains "selected basins  : 3"
    --stdout-contains "divide segments  : 1"
    --stdout-contains "heads without    : 0"
)
set_tests_properties(smoke_divide_migration_merge_fragments PROPERTIES
    LABELS "smoke;cli;divide_migration"
    FIXTURES_REQUIRED divide_migration_fixture)

# --- pair-summary statistics end to end (plan §11.2 L/M) ----------------------
#
# The same ridge fixture and clustering as the merge test above, so the topology
# is known: exactly two units bordering along one merged divide. The pair
# summary must therefore be exactly ONE row, on pair 1|2, with heads on both
# sides -- which is what makes the header assertion worth anything.
#
# The header is pinned in full because it is the published schema (§6.3 plus the
# per-metric valid-value counts this work added). The row assertion pins the
# leading pair_id/unit_a/unit_b, and the trailing uncertainty_method proves the
# parameter was read rather than defaulted: the CLI is run with bootstrap, not
# the std_err default.
add_lsdtt_cli_smoke_test(
    smoke_divide_migration_pair_summary
    lsdtt-divide-migration
    divide_migration_pair_summary.param.in
    divmig_pairsum_divide_pair_summary.csv
    --read-path "${LSDTT_DIVIDE_MIGRATION_FIXTURE_DIR}"
    --expect divmig_pairsum_divides.fgb
    --expect divmig_pairsum_divide_migration_metrics.csv
    --expect-file-contains "divmig_pairsum_divide_pair_summary.csv:pair_id,unit_a,unit_b,n_heads_a,n_heads_b,mean_head_elevation_a,mean_head_elevation_b,delta_head_elevation,n_head_elevation_a,n_head_elevation_b,mean_gradient_a,mean_gradient_b,delta_gradient,n_gradient_a,n_gradient_b,mean_relief_a,mean_relief_b,delta_relief,n_relief_a,n_relief_b,mean_chi_a,mean_chi_b,delta_chi,n_chi_a,n_chi_b,uncertainty_method,uncertainty_head_elevation,uncertainty_gradient,uncertainty_relief,uncertainty_chi"
    --expect-file-contains "divmig_pairsum_divide_pair_summary.csv:1|2,1,2,"
    --expect-file-contains "divmig_pairsum_divide_pair_summary.csv:,bootstrap,"
    --stdout-contains "divide segments  : 1"
    --stdout-contains "pair summaries   : 1"
    --stdout-contains "uncertainty      : bootstrap"
    --stdout-contains "Pair summaries report side A minus side B"
    --stdout-contains "heads without    : 0"
)
set_tests_properties(smoke_divide_migration_pair_summary PROPERTIES
    LABELS "smoke;cli;divide_migration"
    FIXTURES_REQUIRED divide_migration_fixture)

# --- lsdtt-lem: the in-LEM coupled hillslope path -----------------------------
#
# Why this exists (docs/testing_and_benchmarking.md §8.4): the semi-implicit fast
# path applied the DAMPED Picard omega to its single IMEX solve, under-diffusing
# every step. The standalone driver (unit_hillslope_implicit) could not see it --
# its steady-state tests all finish on full Picard -- while coupled LEM steady
# states came out about 40 % HIGH. A fast path that is only taken in regimes the
# unit tests do not finish in is a fast path the unit tests cannot validate.
#
# This run spends essentially all of its steps in that fast path: U*L/D = 0.2 is
# a quarter of Sc = 0.8, well under the max-face-ratio < 0.5 semi-implicit
# threshold. It therefore exercises exactly the path the unit driver structurally
# cannot reach.
#
# The assertion is the CONVERGED STEADY STATE, two-sided:
#   peak 4.9429 m -- the recorded in-LEM strip value (§8.3), against the analytic
#                    Roering peak of ~4.85 and the explicit model's 4.91
#   mean 3.2719 m -- pins the profile SHAPE, not just its crest
# Verified to fail on perturbation: reinstating the omega bug (step_omega =
# cur_omega for a lone semi-implicit solve) moves the peak to 6.862 m (+38.8 %)
# and the mean to 4.516 m, both far outside the bands. Starting from a FLAT
# initial profile is what stops "the solver did nothing" from landing inside the
# band -- inaction leaves the strip at 0 m.
#
# ~0.5 s, so it stays in the default suite.
add_lsdtt_cli_smoke_test(
    smoke_lem_strip_fixture
    lsdtt-raster-synthetic-tests
    lem_strip_fixture.param.in
    lemstrip_hillslope_strip.tif
    --expect lemstrip_hillslope_strip_channel_mask.tif
    --expect-raster-size lemstrip_hillslope_strip.tif:101x11
)
set_tests_properties(smoke_lem_strip_fixture PROPERTIES
    LABELS "smoke;cli;hillslope"
    FIXTURES_SETUP lem_strip_fixture)

add_lsdtt_cli_smoke_test(
    smoke_lem_roering_implicit_strip
    lsdtt-lem
    lem_roering_implicit_strip.param.in
    lemstrip_run_modelled_topography.tif
    --read-path "${CMAKE_CURRENT_BINARY_DIR}/smoke_lem_strip_fixture/out"
    --expect lemstrip_run_model_elevation_stats.csv
    --expect-raster-stat lemstrip_run_modelled_topography.tif:max:4.90:4.99
    --expect-raster-stat lemstrip_run_modelled_topography.tif:mean:3.24:3.31
    --expect-raster-stat lemstrip_run_modelled_topography.tif:min:0.0:0.0
    --expect-valid-percent-100 lemstrip_run_modelled_topography.tif
    --stdout-contains "Hillslope transp. : on  (model=roering_implicit"
    --stdout-contains "Loading channel mask:"
)
set_tests_properties(smoke_lem_roering_implicit_strip PROPERTIES
    LABELS "smoke;cli;hillslope"
    FIXTURES_REQUIRED lem_strip_fixture)

add_lsdtt_cli_smoke_test(
    smoke_lem_equilibration_cycle
    lsdtt-lem
    lem_equilibration_cycle_smoke.param.in
    tiny_eqcycle_steady_state_elevation.tif
    --expect tiny_eqcycle_snap_convergence.csv
    --expect-file-contains "tiny_eqcycle_snap_convergence.csv:max_abs_dz_m"
    --expect-file-contains "tiny_eqcycle_snap_convergence.csv:divide_step_median_m"
    --expect-file-contains "tiny_eqcycle_snap_convergence.csv:divide_gradient_median"
    --expect-file-contains "tiny_eqcycle_snap_convergence.csv:divide_gradient_p90"
    --expect-file-contains "tiny_eqcycle_snap_convergence.csv:divide_gradient_target_used"
    --expect-file-contains "tiny_eqcycle_snap_convergence.csv:divide_gradient_target_margin"
    --expect-file-contains "tiny_eqcycle_snap_convergence.csv:divide_gradient_target_met"
    --stdout-contains "Path           : landscape_equilibration_cycle"
)

add_lsdtt_cli_smoke_test(
    smoke_lem_equilibration_divide_step
    lsdtt-lem
    lem_equilibration_divide_step_smoke.param.in
    tiny_eqdiv_steady_state_elevation.tif
    --expect tiny_eqdiv_divide_step.tif
    --expect tiny_eqdiv_snap_convergence.csv
    --expect-raster-size tiny_eqdiv_divide_step.tif:5x5
    --expect-raster-stat tiny_eqdiv_divide_step.tif:max:0.0001:1000000.0
    --expect-file-contains "tiny_eqdiv_snap_convergence.csv:divide_step_median_m"
)

add_lsdtt_cli_smoke_test(
    smoke_lem_equilibration_divide_step_nodiag
    lsdtt-lem
    lem_equilibration_divide_step_nodiag_smoke.param.in
    tiny_eqdiv_nodiag_steady_state_elevation.tif
    --expect tiny_eqdiv_nodiag_divide_step.tif
    --expect tiny_eqdiv_nodiag_snap_convergence.csv
    --expect-file-contains "tiny_eqdiv_nodiag_snap_convergence.csv:divide_step_median_m"
)

add_lsdtt_cli_smoke_test(
    smoke_lem_equilibration_divide_target
    lsdtt-lem
    lem_equilibration_divide_target_smoke.param.in
    tiny_eqtarget_steady_state_elevation.tif
    --expect tiny_eqtarget_snap_convergence.csv
    --stdout-contains "Divide grad med:"
    --stdout-contains "met"
    --expect-file-contains "tiny_eqtarget_snap_convergence.csv:divide_gradient_median"
    --expect-file-contains "tiny_eqtarget_snap_convergence.csv:divide_gradient_p90"
    --expect-file-contains "tiny_eqtarget_snap_convergence.csv:divide_gradient_target_used"
    --expect-file-contains "tiny_eqtarget_snap_convergence.csv:divide_gradient_target_margin"
    --expect-file-contains "tiny_eqtarget_snap_convergence.csv:divide_gradient_target_met"
)

add_lsdtt_cli_smoke_test(
    smoke_lem_equilibration_divide_target_nodiag
    lsdtt-lem
    lem_equilibration_divide_target_nodiag_smoke.param.in
    tiny_eqtarget_nodiag_steady_state_elevation.tif
    --expect tiny_eqtarget_nodiag_snap_convergence.csv
    --stdout-contains "Divide grad med:"
    --stdout-contains "met"
    --expect-file-contains "tiny_eqtarget_nodiag_snap_convergence.csv:divide_gradient_target_met"
)

add_lsdtt_cli_smoke_test(
    smoke_lem_equilibration_divide_target_notmet
    lsdtt-lem
    lem_equilibration_divide_target_notmet_smoke.param.in
    tiny_eqtarget_notmet_steady_state_elevation.tif
    --expect tiny_eqtarget_notmet_snap_convergence.csv
    --stdout-contains "Divide grad med:"
    --stdout-contains "not met"
    --expect-file-contains "tiny_eqtarget_notmet_snap_convergence.csv:divide_gradient_median"
    --expect-file-contains "tiny_eqtarget_notmet_snap_convergence.csv:divide_gradient_p90"
    --expect-file-contains "tiny_eqtarget_notmet_snap_convergence.csv:divide_gradient_target_used"
    --expect-file-contains "tiny_eqtarget_notmet_snap_convergence.csv:divide_gradient_target_margin"
    --expect-file-contains "tiny_eqtarget_notmet_snap_convergence.csv:divide_gradient_target_met"
)

add_lsdtt_cli_smoke_test(
    smoke_lem_equilibration_convergence
    lsdtt-lem
    lem_equilibration_convergence_smoke.param.in
    tiny_eqconv_steady_state_elevation.tif
    --expect tiny_eqconv_snap_convergence.csv
    --stdout-contains "Converged at cycle"
)

add_lsdtt_cli_smoke_test(
    smoke_lem_equilibration_fixed_cycles
    lsdtt-lem
    lem_equilibration_fixed_cycles_smoke.param.in
    tiny_eqfixed_steady_state_elevation.tif
    --expect tiny_eqfixed_snap_convergence.csv
    --expect-file-contains "tiny_eqfixed_snap_convergence.csv:cycle"
    --stdout-contains "Cycles completed: 4"
)

add_lsdtt_cli_smoke_test(
    failure_lem_equilibration_bad_divide_steps
    lsdtt-lem
    lem_equilibration_bad_divide_steps.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "divide_merge_search_steps must be > 0"
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
    tiny_and_tiny_dem_plus1_difference.tif
    --expect-raster-checksum tiny_and_tiny_dem_plus1_difference.tif:65511
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

add_executable(lsdtt-test-cosmo-ratio cpp/test_cosmo_ratio.cpp)
target_link_libraries(lsdtt-test-cosmo-ratio PRIVATE lsdtt3_cpu)

add_executable(lsdtt-test-cosmo-units cpp/test_cosmo_units.cpp)
target_link_libraries(lsdtt-test-cosmo-units PRIVATE lsdtt3_cpu)

# Analytical step-change rasters: the Royden & Perron / Mitchell & Yanites
# transient solution in chi, plus the imposed critical-slope hillslopes.
# Synthetic ridge-and-valley fixture, no network, no external data.
add_executable(lsdtt-test-step-change cpp/test_step_change.cpp)
target_link_libraries(lsdtt-test-step-change PRIVATE lsdtt3_cpu)

# Landslide size distribution (plan Phase 0): the imposed magnitude-frequency
# distribution and the area-volume scaling, with no landscape involved.
# Synthetic only -- no DEM, no network.
add_executable(lsdtt-test-landslide-size cpp/test_landslide_size.cpp)
target_link_libraries(lsdtt-test-landslide-size PRIVATE lsdtt3_cpu)

# Reproducing Niemi et al. (2005): the ensemble of basin samples, the
# apparent-versus-true erosion rate, and the area dependence of the scatter.
add_executable(lsdtt-test-niemi-ensemble cpp/test_niemi_ensemble.cpp)
target_link_libraries(lsdtt-test-niemi-ensemble PRIVATE lsdtt3_cpu)

# Quartz fertility as a mixing weight, and the ice mask that shares its
# machinery. Synthetic rasters only.
add_executable(lsdtt-test-quartz-fertility cpp/test_quartz_fertility.cpp)
target_link_libraries(lsdtt-test-quartz-fertility PRIVATE lsdtt3_cpu)

# Landslide scars superimposed on a step change, plus the density convention
# (real metres vs mass depth, and the per-material scar density). Synthetic
# rasters only -- no DEM, no flow routing.
add_executable(lsdtt-test-landslide-step-change cpp/test_landslide_step_change.cpp)
target_link_libraries(lsdtt-test-landslide-step-change PRIVATE lsdtt3_cpu)

# Landslide placement (plan Phase 1): eligible hosts and the s_vector footprint.
# Synthetic ridge-and-valley fixture, no network.
add_executable(lsdtt-test-landslide-placement cpp/test_landslide_placement.cpp)
target_link_libraries(lsdtt-test-landslide-placement PRIVATE lsdtt3_cpu)

# Landslide event (plan Phase 2): budget -> N -> placed scars, with mass balance.
add_executable(lsdtt-test-landslide-event cpp/test_landslide_event.cpp)
target_link_libraries(lsdtt-test-landslide-event PRIVATE lsdtt3_cpu)

# Landslide -> CRN coupling (plan Phase 4, T16). Composes existing crn_column
# calls only; no new cosmogenic physics.
add_executable(lsdtt-test-landslide-crn cpp/test_landslide_crn.cpp)
target_link_libraries(lsdtt-test-landslide-crn PRIVATE lsdtt3_cpu)

# Landslide series (plan Phases 3 and 5): multi-event integration, plus the
# Niemi rate-closure and Yanites depth-regime reproductions.
add_executable(lsdtt-test-landslide-series cpp/test_landslide_series.cpp)
target_link_libraries(lsdtt-test-landslide-series PRIVATE lsdtt3_cpu)

# Flow routing / depression handling property tests
# (docs/testing_and_benchmarking.md).  Synthetic fixtures only: no
# network, no RichDEM, so this runs in the default CTest suite.
add_executable(lsdtt-test-flow-properties cpp/test_flow_properties.cpp)
target_link_libraries(lsdtt-test-flow-properties PRIVATE lsdtt3_cpu)

add_executable(lsdtt-test-synthetic-networks cpp/test_synthetic_networks.cpp)
target_link_libraries(lsdtt-test-synthetic-networks PRIVATE lsdtt3_cpu)

add_test(NAME unit_synthetic_networks
         COMMAND lsdtt-test-synthetic-networks
                 ${CMAKE_CURRENT_BINARY_DIR}/synthetic_networks_output)
set_tests_properties(unit_synthetic_networks PROPERTIES LABELS "unit;synthetic_networks")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/synthetic_networks_output)

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

# Divide-migration compute layer (docs/plan_docs/divide_migration_cli_plan.md
# §11.1I).  Covers the routines in src/cpu/lsdtt_divide_migration_cpu.cpp:
# nearest-divide geometry (mid-segment projection, zero-distance guard, negative
# relief, nodata sampling, and the PER-UNIT index that stops a head being matched
# to a divide between two other units), the common-base-level chi rebasing,
# fragment merging, the pair summaries (which are driven by the DIVIDES, so a
# pair with no matched heads still gets a row), and the eight-connected basin
# adjacency pass that decides which pairs are worth a ridgeline extraction.
#
# The same binary writes the synthetic three-valley DEM and outlet layer used by
# the divide-migration CLI end-to-end tests, so the fixture and the in-memory
# geometry come from one generator.  Test domains are generated in C++, never in
# Python (docs/testing_and_benchmarking.md §8).
add_executable(lsdtt-test-divide-migration cpp/test_divide_migration.cpp)
target_link_libraries(lsdtt-test-divide-migration PRIVATE lsdtt3_cpu)

add_executable(lsdtt-test-divide-diagnostics cpp/test_divide_diagnostics.cpp)
target_link_libraries(lsdtt-test-divide-diagnostics PRIVATE lsdtt3_cpu)

add_executable(lsdtt-test-landscape-equilibration cpp/test_landscape_equilibration.cpp)
target_link_libraries(lsdtt-test-landscape-equilibration PRIVATE lsdtt3_cpu)

add_test(NAME unit_divide_diagnostics
         COMMAND lsdtt-test-divide-diagnostics)
set_tests_properties(unit_divide_diagnostics PROPERTIES
    LABELS "unit;divide_diagnostics")

add_test(NAME unit_divide_migration
         COMMAND lsdtt-test-divide-migration
                 "${LSDTT_DIVIDE_MIGRATION_FIXTURE_DIR}")
set_tests_properties(unit_divide_migration PROPERTIES
    LABELS "unit;divide_migration"
    FIXTURES_SETUP divide_migration_fixture)

add_test(NAME unit_landscape_equilibration
         COMMAND lsdtt-test-landscape-equilibration)
set_tests_properties(unit_landscape_equilibration PROPERTIES
    LABELS "unit;landscape_equilibration")

# Implicit (backward-Euler / TR-BDF2) Roering nonlinear hillslope solver
# (docs/testing_and_benchmarking.md).  Rectangular-domain validation against the
# analytic Roering steady profile, the SOR steady-state solver, and the
# conservation / Sc-enforcement / timestep-independence properties.  Synthetic
# surfaces only: no network, no committed DEM, so it runs in the default suite.
# It writes its plottable profile / grid-refinement CSVs into the test working
# directory (consumed by lsdviztools.lsdtransientvalidation).
#
# The Phase 4 linear-solve TIMING sweep in the same driver is NOT part of this
# test -- it is reported, never asserted, and is registered separately as
# performance_hillslope_implicit below.
add_executable(lsdtt-test-hillslope-implicit cpp/test_hillslope_implicit.cpp)
target_link_libraries(lsdtt-test-hillslope-implicit PRIVATE lsdtt3_cpu)

add_test(NAME unit_hillslope_implicit
         COMMAND lsdtt-test-hillslope-implicit
                 ${CMAKE_CURRENT_BINARY_DIR}/hillslope_implicit_output)
set_tests_properties(unit_hillslope_implicit PROPERTIES LABELS "unit;hillslope")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/hillslope_implicit_output)

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
    failure_fetch_raster_bad_nested_stop_level
    lsdtt-fetch-raster
    fetch_raster_bad_nested_stop_level.param.in
    unused_expected_output
    --expect-failure
    --stderr-contains "FATAL: basin_nested_stop_level must be one of: coarse, 90m, 30m"
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

    # Implicit hillslope solver: per-solve linear-algebra cost and grid scaling
    # (128^2 / 256^2 / 512^2) for each Eigen backend, from the SAME driver as
    # unit_hillslope_implicit via --performance.
    #
    # REPORTED, NEVER ASSERTED -- it can only fail by not running. Kept out of the
    # default suite because the 512^2 LDLT solve alone is seconds of wall time and
    # the timings are worthless (and flaky) on a loaded shared runner.
    add_test(NAME performance_hillslope_implicit
             COMMAND lsdtt-test-hillslope-implicit
                     ${CMAKE_CURRENT_BINARY_DIR}/perf_output --performance)
    set_tests_properties(performance_hillslope_implicit
        PROPERTIES LABELS "performance;hillslope" SKIP_RETURN_CODE 77 TIMEOUT 1800)

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

add_executable(lsdtt-test-coastal-mask cpp/test_coastal_mask.cpp)
target_link_libraries(lsdtt-test-coastal-mask PRIVATE lsdtt3_core)
target_include_directories(lsdtt-test-coastal-mask PRIVATE ${CMAKE_SOURCE_DIR}/src/core)

add_test(NAME unit_coastal_mask
         COMMAND lsdtt-test-coastal-mask)
set_tests_properties(unit_coastal_mask PROPERTIES LABELS "unit;fetch")

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
# Step-change rasters. Checks the analytic algebra (ksn, celerity, knickpoint
# position, the A_0 dependence) and then the rasters themselves: time-since-change
# equals the step age at the outlet and falls to zero at the knickpoint, erosion
# rates split old/new across it, the chi profile is piecewise linear, and n < 1
# produces a rarefaction fan. Emits step_change_profiles.csv.
add_test(NAME unit_step_change
         COMMAND lsdtt-test-step-change
                 ${CMAKE_CURRENT_BINARY_DIR}/step_change_test_output)
set_tests_properties(unit_step_change PROPERTIES LABELS "unit;step_change")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/step_change_test_output)

# Landslide size distribution T1-T5b. Emits landslide_truncation.csv, which
# tabulates how much of the distribution the A_max cutoff removes by number and
# by volume (plan Sec 2.5).
add_test(NAME unit_landslide_size
         COMMAND lsdtt-test-landslide-size
                 ${CMAKE_CURRENT_BINARY_DIR}/landslide_test_output)
set_tests_properties(unit_landslide_size PROPERTIES LABELS "unit;landslide")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/landslide_test_output)

# Niemi ensemble N1-N6. The CONTROL (N1) is the one that matters most: with no
# landslides the inversion must return the imposed rate, or nothing else here
# means anything. Emits niemi_basins.csv, niemi_area_bins.csv, niemi_ensemble.csv.
# Its OWN output directory: sharing landslide_test_output with the other
# landslide tests lets them race on the CSVs under `ctest -j`.
add_test(NAME unit_niemi_ensemble
         COMMAND lsdtt-test-niemi-ensemble
                 ${CMAKE_CURRENT_BINARY_DIR}/niemi_test_output)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/niemi_test_output)
set_tests_properties(unit_niemi_ensemble PROPERTIES
    LABELS "unit;cosmo;landslide" TIMEOUT 1800)

# Quartz fertility Q1-Q6 and the ice mask I1-I3. Pins the case that makes
# fertility a correctness issue rather than a refinement: a quartz-free unit
# eroding ten times faster than its neighbour must contribute nothing.
add_test(NAME unit_quartz_fertility COMMAND lsdtt-test-quartz-fertility)
set_tests_properties(unit_quartz_fertility PROPERTIES LABELS "unit;cosmo")

# Landslide scar coupling and densities. Emits landslide_scar_recovery.csv,
# which shows how a scar's concentration recovers with age at three depths and
# what carrying it at the wrong material density does to the answer.
add_test(NAME unit_landslide_step_change
         COMMAND lsdtt-test-landslide-step-change
                 ${CMAKE_CURRENT_BINARY_DIR}/landslide_test_output)
set_tests_properties(unit_landslide_step_change PROPERTIES
    LABELS "unit;landslide;step_change")

# Landslide placement T6-T10d. Emits landslide_host_fill.csv, which records how
# completely landslides of each size fill their host basin (plan Sec 3.2).
add_test(NAME unit_landslide_placement
         COMMAND lsdtt-test-landslide-placement
                 ${CMAKE_CURRENT_BINARY_DIR}/landslide_test_output)
set_tests_properties(unit_landslide_placement PROPERTIES LABELS "unit;landslide")

# Landslide event T11-T13b. Emits landslide_budget_vs_amax.csv.
add_test(NAME unit_landslide_event
         COMMAND lsdtt-test-landslide-event
                 ${CMAKE_CURRENT_BINARY_DIR}/landslide_test_output)
set_tests_properties(unit_landslide_event PROPERTIES LABELS "unit;landslide")

# Landslide CRN response T16. Emits landslide_crn_response.csv.
add_test(NAME unit_landslide_crn
         COMMAND lsdtt-test-landslide-crn
                 ${CMAKE_CURRENT_BINARY_DIR}/landslide_test_output)
set_tests_properties(unit_landslide_crn PROPERTIES LABELS "unit;landslide;cosmo")

# Landslide series: Phase 3 multi-event + Phase 5 T14/T15. Emits
# landslide_niemi_closure.csv and landslide_yanites_regimes.csv.
add_test(NAME unit_landslide_series
         COMMAND lsdtt-test-landslide-series
                 ${CMAKE_CURRENT_BINARY_DIR}/landslide_test_output)
set_tests_properties(unit_landslide_series PROPERTIES LABELS "unit;landslide" TIMEOUT 900)

# The isotope-ratio output policy. One division per pixel, and every failure
# mode is silent: an inf written into a raster a basin average then consumes, a
# nodata numerator divided and emitted as a number, a mask disagreeing with its
# own data. R8 pins the 2026-08-05 promotion from a lambda in
# lsdtt-cosmo-prediction as BITWISE faithful on the input shape it was written
# for, so the refactor is proved at unit level and not only through the CLI
# oracle.
# The metre <-> g/cm2 conversions. Four-line functions, and worth pinning
# because the density split of 2026-08-04 exists precisely because one of the
# eight hand-written `x * rho * 0.1` sites chose the wrong density. U7 asserts
# the promoted functions reproduce those hand-written expressions BITWISE --
# `0.1 * rho * m` and `m * rho * 0.1` are not the same double, and the refactor
# oracle compares file hashes.
add_test(NAME unit_cosmo_units COMMAND lsdtt-test-cosmo-units)
set_tests_properties(unit_cosmo_units PROPERTIES LABELS "unit;cosmo")

add_test(NAME unit_cosmo_ratio COMMAND lsdtt-test-cosmo-ratio)
set_tests_properties(unit_cosmo_ratio PROPERTIES LABELS "unit;cosmo")

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

# SA2 rigorous check: in multi-isotope step-change + accumulation mode, sampled
# ratios must be computed as ratios of accumulated basin means (not as means of
# per-pixel ratios).
add_test(NAME rigorous_cosmo_multi_accum_ratio
         COMMAND ${Python3_EXECUTABLE}
                 "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cosmo_multi_accum_ratio_check.py"
                 --pred-exe "$<TARGET_FILE:lsdtt-cosmo-prediction>"
                 --data-dir "${LSDTT_TEST_DATA_DIR}")
set_tests_properties(rigorous_cosmo_multi_accum_ratio PROPERTIES
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

# SA2 end-to-end path: one multi-isotope step-change run produces accumulated
# basin means, accumulated ratio rasters and sampled point ratios without
# rerunning for each isotope.
add_lsdtt_cli_smoke_test(
    smoke_cosmo_prediction_multi_isotope_accum_sampling
    lsdtt-cosmo-prediction
    cosmo_prediction_multi_isotope_accum_sampling.param.in
    mi_acc_sampled_crn_conc.csv
    --expect mi_acc_Be10_crn_acc.tif
    --expect mi_acc_Al26_crn_acc.tif
    --expect mi_acc_C14_crn_acc.tif
    --expect mi_acc_R2610_acc.tif
    --expect mi_acc_R1410_acc.tif
    --expect mi_acc_sampled_crn_conc_Be10.csv
    --expect mi_acc_sampled_crn_conc_Al26.csv
    --expect mi_acc_sampled_crn_conc_C14.csv
    --expect-file-contains "mi_acc_sampled_crn_conc.csv:Be10_CONC"
    --expect-file-contains "mi_acc_sampled_crn_conc.csv:R2610_ACC"
    --expect-file-contains "mi_acc_sampled_crn_conc.csv:R1410_ACC"
    --expect-file-contains "mi_acc_sampled_crn_conc_Be10.csv:Be10_CONC"
    --expect-file-contains "mi_acc_sampled_crn_conc_Al26.csv:Al26_CONC"
    --expect-file-contains "mi_acc_sampled_crn_conc_C14.csv:C14_CONC"
    --stdout-contains "Accumulating multi-isotope concentration downstream"
)
set_tests_properties(smoke_cosmo_prediction_multi_isotope_accum_sampling
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
# lsdtt-cosmo-landslides (plan Phase 6).
add_test(NAME smoke_cosmo_landslides_help COMMAND lsdtt-cosmo-landslides -h)
set_tests_properties(smoke_cosmo_landslides_help PROPERTIES
    LABELS "smoke;cli;landslide"
    PASS_REGULAR_EXPRESSION "long_term_erosion_rate_mm_kyr")

# End-to-end on the tiny fixture: an event is placed and the inventory written.
add_lsdtt_cli_smoke_test(
    smoke_cosmo_landslides_single_event
    lsdtt-cosmo-landslides
    cosmo_landslides_single_event.param.in
    lssmoke_LS_INVENTORY.csv
    --stdout-contains "landslide event:"
    --expect-file-contains "lssmoke_LS_INVENTORY.csv:n_pixels"
)

# R2c: CRN concentrations in lsdtt-cosmo-landslides can now build a real site
# (per-pixel production + shielding) instead of only the scalar convenience
# path.
add_lsdtt_cli_smoke_test(
    smoke_cosmo_landslides_crn_site_mode
    lsdtt-cosmo-landslides
    cosmo_landslides_crn_site_mode.param.in
    lscrn_LS_CONC.tif
    --stdout-contains "CRN site mode: per-pixel"
)

add_lsdtt_cli_smoke_test(
    smoke_cosmo_landslides_sample_ensemble
    lsdtt-cosmo-landslides
    cosmo_landslides_sample_ensemble.param.in
    lsens_LS_SAMPLE_EFFECT.csv
    --expect-file-contains "lsens_LS_SAMPLE_EFFECT.csv:ratio_to_control_median"
)

# A geographic-CRS DEM must exit cleanly with a diagnosis, NOT abort on an
# uncaught throw. The deeper check throws; the CLI catches it.
add_lsdtt_cli_smoke_test(
    failure_cosmo_landslides_geographic_crs
    lsdtt-cosmo-landslides
    cosmo_landslides_geographic.param.in
    lsgeo_LS_INVENTORY.csv
    --expect-failure
    --stderr-contains "projected"
)

# ---------------------------------------------------------------------------
# Cosmogenic parameter validation in lsdtt-cosmo-landslides (plan R2b).
#
# Until 2026-08-05 a bad muon_scheme or nuclide was caught inside the CRN block,
# downgraded to a stderr warning, and the tool exited **0** having written no
# concentration raster. To a scripted workflow that is a success that silently
# produced nothing -- a worse failure mode than a crash, because nothing
# notices. All three cases below were measured on the shipped binary before the
# fix; the first two exited 0.
add_lsdtt_cli_smoke_test(
    failure_cosmo_landslides_bad_muon_scheme
    lsdtt-cosmo-landslides
    cosmo_landslides_bad_muon_scheme.param.in
    lsbadmuon_LS_CONC.tif
    --expect-failure
    --stderr-contains "unknown muon_scheme"
)

add_lsdtt_cli_smoke_test(
    failure_cosmo_landslides_bad_scaling_scheme
    lsdtt-cosmo-landslides
    cosmo_landslides_bad_scaling_scheme.param.in
    lsbadscale_LS_CONC.tif
    --expect-failure
    --stderr-contains "unknown scaling_scheme"
)

# VALIDITY IS NOT CAPABILITY. is_valid_muon_scheme() accepts the CRONUS v3
# emulator schemes, but they are pressure- and nuclide-dependent and this tool
# reaches the column through reset_scaling(), which throws for them. A shared
# validator that only checked the name would accept something that fails later.
# The message must also point at the tool that DOES support them.
add_lsdtt_cli_smoke_test(
    failure_cosmo_landslides_emulator_scheme
    lsdtt-cosmo-landslides
    cosmo_landslides_emulator_scheme.param.in
    lsemul_LS_CONC.tif
    --expect-failure
    --stderr-contains "lsdtt-cosmo-prediction"
)

# The other side of the guard, and the reason it is conditional: muon_scheme is
# INERT with no cosmogenic output requested, so a stale line must not fail the
# run. Without this test the natural "fix" is to validate unconditionally, which
# would break every pure-geomorphology parameter file carrying a cosmo line.
add_lsdtt_cli_smoke_test(
    smoke_cosmo_landslides_geomorph_only_ignores_muon_scheme
    lsdtt-cosmo-landslides
    cosmo_landslides_geomorph_only.param.in
    lsgeomorph_LS_INVENTORY.csv
    --expect lsgeomorph_LS_TFAIL.tif
)

add_test(NAME smoke_cosmo_tool_help COMMAND lsdtt-cosmo-tool -h)
set_tests_properties(smoke_cosmo_tool_help PROPERTIES
    LABELS "smoke;cli;cosmo"
    PASS_REGULAR_EXPRESSION "sample_age_yr")

add_test(NAME smoke_divide_migration_help COMMAND lsdtt-divide-migration -h)
set_tests_properties(smoke_divide_migration_help PROPERTIES
    LABELS "smoke;cli"
    PASS_REGULAR_EXPRESSION "basin_outlet_file")

# Generate the step change from the analytical fluvial solution and feed the
# generated time-since-change field straight into the concentration kernel.
# Checks that the derived K is reported (it is back-calculated, so the whole
# configuration hangs off it) and that the tag raster is written.
add_lsdtt_cli_smoke_test(
    smoke_cosmo_prediction_step_change
    lsdtt-cosmo-prediction
    cosmo_prediction_step_change.param.in
    scgen_Be10_step_crn_conc.tif
    --stdout-contains "K (back-calculated)"
    --stdout-contains "wave form"
    --expect scgen_SC_STEPTIME.tif
    --expect scgen_SC_TAG.tif
)

# The generated time-since-change field is consumed in memory. This asserts
# that route agrees with writing it to disk and reading it back in a second
# run -- the workflow the legacy tool forced. See the script header for why the
# comparison is a checksum plus a 1e-6 statistics tolerance rather than bitwise.
add_test(NAME regression_step_change_disk_roundtrip
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_step_change_roundtrip.py
                 --exe $<TARGET_FILE:lsdtt-cosmo-prediction>
                 --data-dir ${LSDTT_TEST_DATA_DIR}
                 --work-dir ${CMAKE_CURRENT_BINARY_DIR}/step_change_roundtrip)
set_tests_properties(regression_step_change_disk_roundtrip PROPERTIES
    LABELS "cli;cosmo;step_change")

# The scar coupling end to end: lsdtt-cosmo-landslides writes the failure-time and
# excavation rasters, lsdtt-cosmo-prediction consumes them. Asserts that
# switching the coupling on leaves unfailed ground BITWISE unchanged -- the
# property that stops a landslide run silently perturbing the whole raster.
add_test(NAME regression_landslide_scar_coupling
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_landslide_scar_coupling.py
                 --landslides-exe $<TARGET_FILE:lsdtt-cosmo-landslides>
                 --prediction-exe $<TARGET_FILE:lsdtt-cosmo-prediction>
                 --data-dir ${LSDTT_TEST_DATA_DIR}
                 --work-dir ${CMAKE_CURRENT_BINARY_DIR}/landslide_scar_coupling)
set_tests_properties(regression_landslide_scar_coupling PROPERTIES
    LABELS "cli;cosmo;landslide;step_change")

# SC path smoke/regression: the combined multi-isotope step-change + landslide
# coupling route that was previously blocked by a fatal guard.
add_test(NAME regression_cosmo_multi_landslide_coupling
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/scripts/run_cosmo_multi_landslide_coupling.py
                 --landslides-exe $<TARGET_FILE:lsdtt-cosmo-landslides>
                 --prediction-exe $<TARGET_FILE:lsdtt-cosmo-prediction>
                 --data-dir ${LSDTT_TEST_DATA_DIR}
                 --work-dir ${CMAKE_CURRENT_BINARY_DIR}/cosmo_multi_landslide_coupling)
set_tests_properties(regression_cosmo_multi_landslide_coupling PROPERTIES
    LABELS "cli;cosmo;landslide;step_change")

# ---------------------------------------------------------------------------
# The REFACTOR ORACLE (plan R0a).
#
# tests/data/cosmo_oracle.csv holds a bitwise hash of every raster and CSV the
# four cosmogenic CLIs produce across a matrix of scaling scheme, muon scheme,
# shielding path and prediction mode. It exists because the production/pressure/
# shielding setup is about to be extracted from the CLIs into a shared
# component, and THE ORACLE FOR A REFACTOR IS THE PRE-REFACTOR BINARY -- which
# stops existing the moment the refactor lands.
#
# Every subsequent change is one of two kinds and must declare which:
#   behaviour-PRESERVING  this test must stay green. That is the proof.
#   behaviour-CHANGING    re-freeze deliberately (--freeze) and say in the
#                         commit message which artefacts moved and why.
#
# Do NOT re-freeze to make a red test pass.
#
# Verified to fail: a 1-part-in-1e6 perturbation of the production raster moves
# 16 of the 22 cheap-tier artefacts, including lsdtt-cosmo-tool's inverted
# erosion-rate CSVs -- so the oracle reaches through the whole forward/inverse
# chain, not just the raster it perturbs.
add_test(NAME regression_cosmo_oracle
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/scripts/freeze_cosmo_oracle.py
                 --exe lsdtt-cosmo-prediction=$<TARGET_FILE:lsdtt-cosmo-prediction>
                 --exe lsdtt-cosmo-tool=$<TARGET_FILE:lsdtt-cosmo-tool>
                 --exe lsdtt-cosmo-column-tool=$<TARGET_FILE:lsdtt-cosmo-column-tool>
                 --exe lsdtt-cosmo-landslides=$<TARGET_FILE:lsdtt-cosmo-landslides>
                 --read-path ${LSDTT_TEST_DATA_DIR}
                 --oracle ${LSDTT_TEST_DATA_DIR}/cosmo_oracle.csv
                 --work-dir ${CMAKE_CURRENT_BINARY_DIR}/cosmo_oracle
                 --tier cheap)
set_tests_properties(regression_cosmo_oracle PROPERTIES
    LABELS "cli;cosmo;regression")

# The LSDn and CRONUS-emulator rows. Separated because they need the LSDn
# scaling grids in crn_parameters.nc and carry the most setup state -- which is
# exactly why they matter to the extraction, and why they are worth the label
# rather than being dropped.
add_test(NAME regression_cosmo_oracle_heavy
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/scripts/freeze_cosmo_oracle.py
                 --exe lsdtt-cosmo-prediction=$<TARGET_FILE:lsdtt-cosmo-prediction>
                 --exe lsdtt-cosmo-tool=$<TARGET_FILE:lsdtt-cosmo-tool>
                 --exe lsdtt-cosmo-column-tool=$<TARGET_FILE:lsdtt-cosmo-column-tool>
                 --exe lsdtt-cosmo-landslides=$<TARGET_FILE:lsdtt-cosmo-landslides>
                 --read-path ${LSDTT_TEST_DATA_DIR}
                 --oracle ${LSDTT_TEST_DATA_DIR}/cosmo_oracle.csv
                 --work-dir ${CMAKE_CURRENT_BINARY_DIR}/cosmo_oracle_heavy
                 --tier all)
set_tests_properties(regression_cosmo_oracle_heavy PROPERTIES
    LABELS "cli;cosmo;regression;heavy")

# Two prediction modes that write the SAME filenames must not both run: the
# transient goes first and the step change overwrote all five of its outputs,
# so the run reported success having produced half its results. Found while
# promoting the ratio helper (plan SA1); guarded rather than renamed, because
# changing shipped output filenames is a separate decision.
add_lsdtt_cli_smoke_test(
    failure_cosmo_prediction_both_multi_modes
    lsdtt-cosmo-prediction
    cosmo_prediction_both_multi_modes.param.in
    bothmulti_Be10_step_crn_conc.tif
    --expect-failure
    --stderr-contains "silently overwrite"
)

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

# =============================================================================
# REDUNDANCY GUARD: no hand-written metre <-> g/cm2 conversion (plan R0b)
#
# `x * rho * 0.1` was open-coded at nine sites across the library and two CLIs,
# each independently choosing WHICH density to use. That is not a style
# complaint: the density split of 2026-08-04 exists because one of those choices
# was wrong, and a soil scar carried at bedrock density is 1.77x too deep in
# mass -- a 2.05x error in spallation production at the scar base.
#
# lsdtt_crn_units.hpp now owns the arithmetic, with a REQUIRED density argument.
# This guard stops the pattern coming back one call site at a time, which is how
# it accumulated in the first place. It is a CONFIGURE ERROR, matching the
# "unclassified test" guard below: a lint nobody runs is a lint that does not
# exist.
#
# EXEMPT, deliberately:
#   * src/cpu/lsdtt_crn_column.{hpp,cpp} and lsdtt_crn_particle.hpp -- LOCKED
#     NUMERICS. to_eff(), remove_mass() and the four depth_m recoveries live
#     inside the validated core and must not be touched. They are excluded by
#     NOT being in the glob below; that is deliberate, not an oversight.
#     See feedback_cosmo_numerics_locked.
#   * lsdtt_crn_units.hpp itself, which is where the arithmetic now lives.
# =============================================================================

set(_lsdtt_unit_guard_files "")
file(GLOB_RECURSE _lsdtt_conv_candidates
     "${CMAKE_SOURCE_DIR}/src/cpu/lsdtt_landslide_*.cpp"
     "${CMAKE_SOURCE_DIR}/src/cpu/lsdtt_cosmo_*.cpp"
     "${CMAKE_SOURCE_DIR}/src/cpu/lsdtt_quartz_*.cpp"
     "${CMAKE_SOURCE_DIR}/src/cpu/lsdtt_step_change_*.cpp"
     "${CMAKE_SOURCE_DIR}/src/cpu/lsdtt_niemi_*.cpp"
     "${CMAKE_SOURCE_DIR}/src/cpu/lsdtt_landslide_*.hpp"
     "${CMAKE_SOURCE_DIR}/src/cpu/lsdtt_cosmo_*.hpp"
     "${CMAKE_SOURCE_DIR}/src/cpu/lsdtt_quartz_*.hpp"
     "${CMAKE_SOURCE_DIR}/src/cpu/lsdtt_step_change_*.hpp"
     "${CMAKE_SOURCE_DIR}/src/cpu/lsdtt_niemi_*.hpp"
     "${CMAKE_SOURCE_DIR}/cli/lsdtt-cosmo-*.cpp")

foreach(_f IN LISTS _lsdtt_conv_candidates)
    file(STRINGS "${_f}" _hits REGEX "(density|rho)[a-z_0-9]* *\\* *0\\.1[^0-9]")
    foreach(_h IN LISTS _hits)
        # Skip comments and help/doc strings -- they describe the conversion,
        # which is exactly what we want them to keep doing.
        string(STRIP "${_h}" _stripped)
        if(NOT _stripped MATCHES "^//" AND NOT _stripped MATCHES "^\\*" AND
           NOT _stripped MATCHES "^\"")
            list(APPEND _lsdtt_unit_guard_files "${_f}: ${_stripped}")
        endif()
    endforeach()
endforeach()

if(_lsdtt_unit_guard_files)
    string(REPLACE ";" "\n    " _guard_report "${_lsdtt_unit_guard_files}")
    message(FATAL_ERROR
        "Hand-written metre <-> g/cm2 conversion found:\n    ${_guard_report}\n\n"
        "Use lsdtt_crn_units.hpp (metres_to_g_cm2 / m_yr_to_eff_erate / "
        "g_cm2_to_metres).\nThe density argument is REQUIRED there, which is the "
        "whole point: the bug this\nprevents was never the arithmetic, it was "
        "using 2650 for a soil scar.")
endif()

# =============================================================================
# COMPARTMENTALISATION: domain and cost labels
#
# Running all 100+ tests after every change is wasteful: most of them cannot be
# affected by a given edit, and a handful are computationally heavy enough to
# dominate the wall clock (the implicit hillslope solver saturates every core
# for ~9 minutes on its own; the real-world cosmogenic tiers take longer still).
#
# So every test carries a DOMAIN label saying which part of the codebase it
# exercises, and the expensive ones carry a COST label so they can be excluded.
# The intended workflow:
#
#   ctest -L cosmo -LE heavy      # I changed cosmogenic code
#   ctest -L network -LE heavy    # I changed flow routing / chi / basins
#   ctest -LE heavy               # everything cheap; the normal inner loop
#   ctest                         # everything; before a commit
#
# The domain is derived from the TEST NAME rather than declared at each of the
# hundred-odd call sites, so there is one table to read and one place to change.
# A test whose name matches no rule is a CONFIGURE ERROR, not a silently
# unlabelled test -- that is what stops the scheme rotting as tests are added.
#
# DOMAIN RULES, first match wins. Add a rule here when you add a subsystem.
# =============================================================================

set(LSDTT_DOMAIN_RULES
    "synthetic_networks|_ocn_" "synthetic_networks"
    "niemi"                  "cosmo"
    "quartz|fertility"       "quartz"
    "cosmo"                  "cosmo"
    "landslide"              "landslide"
    "step_change"            "step_change"
    "landscape_equilibration" "landscape_equilibration"
    "divide_diagnostics"     "divide_diagnostics"
    "divide_migration"       "divide_migration"
    "fetch"                  "fetch"
    "hydraulic"              "hydraulic"
    "surface_metrics"        "surface_metrics"
    "hillslope"              "hillslope"
    "lem"                    "lem"
    "chi|channel|drainage|basin|flow_properties"          "network"
    "raster|nodata|fill|tif|coastal|dem_grouping"         "raster"
    "param|secret"           "core"
)

# COST: tests that are expensive enough to keep out of the inner loop. Measured
# 2026-08-04 on an 8-core machine; the figure is wall-clock when the test has
# the machine to itself.
#
#   unit_hillslope_implicit         ~9 min, saturates all 8 cores
#   smoke_lem_roering_implicit_strip ~5 min, likewise
#   cosmo_realworld_*                67 s to 30 min, and most need the network
set(LSDTT_HEAVY_TESTS
    unit_hillslope_implicit
    smoke_lem_roering_implicit_strip
    cosmo_realworld_clapp2002
    cosmo_realworld_clapp2002_lsdn
    cosmo_realworld_dethier2014
    cosmo_realworld_dethier2014_lsdn
    cosmo_realworld_dibiase2010
    cosmo_realworld_dibiase2010_lsdn
    cosmo_realworld_lsdn_pixel_vs_lumped
    cosmo_realworld_production_isolation
)

# Tests that use every core they can get. Under `ctest -j` these must not share
# the machine: doing so does not merely slow them, it starves them past their
# timeout, which is how unit_hillslope_implicit came to "fail" at 1500 s while
# taking 9 minutes when run alone.
set(LSDTT_SERIAL_TESTS
    unit_hillslope_implicit
    smoke_lem_roering_implicit_strip
)

get_property(_lsdtt_tests DIRECTORY PROPERTY TESTS)
set(_lsdtt_unclassified "")

foreach(_test IN LISTS _lsdtt_tests)
    set(_domains "")
    list(LENGTH LSDTT_DOMAIN_RULES _n_rules)
    math(EXPR _last "${_n_rules} - 1")
    foreach(_i RANGE 0 ${_last} 2)
        if(_domains STREQUAL "")
            list(GET LSDTT_DOMAIN_RULES ${_i} _pattern)
            math(EXPR _j "${_i} + 1")
            list(GET LSDTT_DOMAIN_RULES ${_j} _domain)
            if(_test MATCHES "${_pattern}")
                set(_domains "${_domain}")
            endif()
        endif()
    endforeach()

    if(_domains STREQUAL "")
        list(APPEND _lsdtt_unclassified "${_test}")
        continue()
    endif()

    # The implicit solvers are hillslope physics wherever they are driven from,
    # so an LEM test that runs one belongs to both compartments.
    if(_test MATCHES "roering|implicit" AND NOT _domains STREQUAL "hillslope")
        list(APPEND _domains "hillslope")
    endif()

    # lsdtt-cosmo-landslides is BOTH a cosmogenic tool and the landslide
    # module's CLI. "cosmo" precedes "landslide" in the rules above and first
    # match wins, so without this its tests would silently leave the landslide
    # compartment at the 2026-08-05 rename and `ctest -L landslide` would stop
    # seeing the CLI for the module it is named after.
    if(_test MATCHES "landslide" AND NOT _domains STREQUAL "landslide")
        list(APPEND _domains "landslide")
    endif()

    if("${_test}" IN_LIST LSDTT_HEAVY_TESTS)
        list(APPEND _domains "heavy")
    endif()

    get_test_property(${_test} LABELS _existing)
    if(_existing)
        list(APPEND _domains ${_existing})
    endif()
    list(REMOVE_DUPLICATES _domains)
    set_tests_properties(${_test} PROPERTIES LABELS "${_domains}")

    if("${_test}" IN_LIST LSDTT_SERIAL_TESTS)
        set_tests_properties(${_test} PROPERTIES RUN_SERIAL TRUE)
        get_test_property(${_test} TIMEOUT _has_timeout)
        if(NOT _has_timeout)
            set_tests_properties(${_test} PROPERTIES TIMEOUT 3600)
        endif()
    endif()
endforeach()

if(NOT _lsdtt_unclassified STREQUAL "")
    string(REPLACE ";" "\n    " _pretty "${_lsdtt_unclassified}")
    message(FATAL_ERROR
        "These tests match no domain rule and would be unreachable by a "
        "compartment run:\n    ${_pretty}\n"
        "Add a rule to LSDTT_DOMAIN_RULES above, or rename the test so it "
        "carries its subsystem in the name. Do not leave a test unlabelled -- "
        "it will silently stop being run by anyone working on that subsystem.")
endif()
