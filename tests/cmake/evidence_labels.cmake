# evidence_labels.cmake
#
# Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
#
# Copyright (C) 2026 Simon M. Mudd
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# EVIDENCE AND LEVEL LABELS
#
# Included by tests/CMakeLists.txt once every test is registered. Every test
# gets two more kinds of label, beside its domain and its tier / flag / cost
# labels:
#
#   EVIDENCE, one or more: how the test knows the right answer (its oracle).
#   The sections of docs/tex_technical/testing_regime.tex, one label each:
#
#     exact_solution    output compared with an answer worked out analytically
#                       or by hand from the fixture's construction
#     property          a statement every correct output satisfies: no pits,
#                       conservation, weights sum to one, monotonicity, range,
#                       no nodata or NaN, a write/read round trip
#     metamorphic       two lsdtt3 runs on related inputs or settings whose
#                       outputs must be related: rotation, thread count, two
#                       methods or interfaces that must agree, a re-run from
#                       the ingested parameters
#     differential      comparison with an implementation lsdtt3 does not own
#                       (RichDEM, TopoToolbox, Whitebox, CRONUS-Earth, HYPE,
#                       DAGGER, OpenCV, Barnes' code, pyproj, scipy.odr,
#                       published tables), including frozen outputs of those
#                       programs and a brute-force reference written inside
#                       the test independently of the code under test
#     characterisation  lsdtt3's own earlier output stored and compared:
#                       checksums, sha256, frozen files, pinned measured
#                       numbers, measured bands, and, for a refactor, the
#                       replaced code kept verbatim as the reference
#     smoke             the program runs to completion and writes its outputs
#     negative          invalid input or configuration is refused with the
#                       expected message or exit code
#     performance       timings reported, never asserted (no evidence of
#                       correctness; kept out of the default suite)
#     conformance       help text, documentation, version strings, naming
#                       rules, build currency or the test tooling conform to a
#                       stated rule; no computed output is judged
#
#   LEVEL, exactly one:
#
#     unit              lsdtt3 routines called in process (a C++ driver), or a
#                       script that inspects repository files or test tooling
#     system            an output judged was written by a shipped program run
#                       end to end: a CLI (directly, through a harness, or as
#                       a fixture this test compares), the project's configure,
#                       or the built Python package
#
# Select ANCHORED, because -L is a regex: `ctest -L '^metamorphic$'`,
# `ctest -L '^unit$'`. `ctest -L '^smoke$'` is every test with the smoke oracle,
# which includes most CLI tests that also pin checksums.
#
# WHERE THE LABELS COME FROM, unioned for evidence:
#
#   1. registration: lsdtt_mark_harness_test() reads a run_cli_smoke.py test's
#      assertion flags (negative or smoke; characterisation for a checksum or
#      sha256; property for an all-valid raster or a no-NaN file);
#      lsdtt_add_lem_k_compare() records exact_solution or metamorphic.
#   2. LSDTT_RERUN_FROM_INGESTED_TESTS: metamorphic (the run is repeated from
#      its ingested parameters and must write the same bytes).
#   3. the name rules below (failure_*, fixture_*).
#   4. the hand table below, one lsdtt_evidence() row per test that needs more
#      than the rules give. A row ADDS evidence and SETS the level.
#
# Level: the table row, else the level recorded at registration, else the name
# rule, else `unit` for a native (non-Python) unit_* test.
#
# This pass is the ONLY writer of these labels. A test that arrives here already
# carrying one (or a retired kind label: rigorous, regression, failure) stops
# configure, as does a test with no evidence, a level other than exactly one, a
# table row naming no test, or a label outside the two lists.
# =============================================================================

set(LSDTT_EVIDENCE_LABELS
    exact_solution property metamorphic differential characterisation
    smoke negative performance conformance)
set(LSDTT_LEVEL_LABELS unit system)
# Kind labels retired on 2026-11-04 in favour of the two axes above.
set(LSDTT_RETIRED_KIND_LABELS rigorous regression failure)

# Name rules, first match wins: pattern, evidence, level.
set(LSDTT_EVIDENCE_NAME_RULES
    # Every failure_* test runs a CLI through the harness and expects it to
    # stop with a given message or exit code.
    "^failure_"  "negative"  "system"
    # Fixture writers: they pass when the fixture is written. A harness fixture
    # (a CLI run) is recorded as system at registration; the rest are C++
    # generators, `cmake -E copy` or Python generators.
    "^fixture_"  "smoke"     "unit"
)

set(_lsdtt_ev_rows "")
set(_lsdtt_ev_errors "")

# lsdtt_evidence(<test> <unit|system> <evidence>...)
#   One hand-table row. A test listed twice, an unknown level or an unknown
#   evidence label is a configure error.
function(lsdtt_evidence TEST_NAME LEVEL)
    if("${TEST_NAME}" IN_LIST _lsdtt_ev_rows)
        set(_lsdtt_ev_errors "${_lsdtt_ev_errors};${TEST_NAME}: listed twice in the evidence table" PARENT_SCOPE)
        return()
    endif()
    list(LENGTH LEVEL _n_levels)
    if(NOT _n_levels EQUAL 1 OR NOT "${LEVEL}" IN_LIST LSDTT_LEVEL_LABELS)
        string(REPLACE ";" "," _given "${LEVEL}")
        set(_lsdtt_ev_errors "${_lsdtt_ev_errors};${TEST_NAME}: level '${_given}' is not exactly one of unit or system" PARENT_SCOPE)
    endif()
    if(NOT ARGN)
        set(_lsdtt_ev_errors "${_lsdtt_ev_errors};${TEST_NAME}: evidence table row names no evidence" PARENT_SCOPE)
    endif()
    foreach(_ev IN LISTS ARGN)
        if(NOT "${_ev}" IN_LIST LSDTT_EVIDENCE_LABELS)
            set(_lsdtt_ev_errors "${_lsdtt_ev_errors};${TEST_NAME}: '${_ev}' is not an evidence label" PARENT_SCOPE)
        endif()
    endforeach()
    set(_lsdtt_ev_rows "${_lsdtt_ev_rows};${TEST_NAME}" PARENT_SCOPE)
    set(_lsdtt_ev_level_${TEST_NAME} "${LEVEL}" PARENT_SCOPE)
    set(_lsdtt_ev_table_${TEST_NAME} "${ARGN}" PARENT_SCOPE)
endfunction()

# =============================================================================
# HAND TABLE
#
# Reasons are in the comments beside each group; the test's own registration
# comment and driver header are the authority. Labels only: no row changes what
# a test runs. Rows for the locked numerics (cosmogenic CRN/CAIRN core, the
# hillslope solvers) are marked; their drivers and pins are not touched here.
# =============================================================================

# -----------------------------------------------------------------------------
# Help and documentation checks (`lsdtt-<cli> -h`). A parameter name in the
# PASS_REGULAR_EXPRESSION says the help renders (smoke); documentation prose
# or a help default says the text conforms (conformance).
# -----------------------------------------------------------------------------
foreach(_t IN ITEMS
        synthetic_networks distributed_hydrology internal_drainage fetch_raster
        fetch_lithology fetch_erosion_rates cosmo_landslides cosmo_tool
        divide_migration cosmo_prediction cosmo_column_tool flat_surface_finder
        valley_centreline fluvial_surfaces)
    lsdtt_evidence(smoke_${_t}_help system smoke)
endforeach()
foreach(_t IN ITEMS drainage_analysis fluvial_surfaces valley_centreline
                    stream_burning feature_finder)
    lsdtt_evidence(smoke_${_t}_threshold_help system smoke)
endforeach()
lsdtt_evidence(smoke_feature_finder_fill_help system smoke)
foreach(_t IN ITEMS dem_source dem_source_polar dem_source_wcs
                    dem_source_hrdem_swiss dem_source_grid_keyed dem_source_gsi
                    linear_feature)
    lsdtt_evidence(smoke_fetch_raster_help_${_t} system conformance)
endforeach()
lsdtt_evidence(smoke_raster_preprocessing_help_rotation system conformance)
foreach(_t IN ITEMS
        surface_metrics_n_threads surface_metrics_frangi_c
        synthetic_surfaces_n_rows synthetic_surfaces_n_cols
        spectral_analysis_n_threads spectral_analysis_valley_spacing_area_exponent
        valley_centreline_write_output_format
        feature_finder_big_channel_threshold_contributing_pixels
        mountain_front_basin_minimum_contributing_pixels
        lem_incision_threshold_specific_area lem_channel_mask_threshold_specific_area
        drainage_analysis_connectivity_area_threshold_m2
        cosmo_prediction_builtin_io)
    lsdtt_evidence(regression_${_t}_help system conformance)
endforeach()
# FAIL_REGULAR_EXPRESSION on -h: no credential is a parameter.
lsdtt_evidence(rigorous_fetch_raster_no_credential_parameters system conformance)

# -----------------------------------------------------------------------------
# Harness system tests that check more than the flags show. The harness gives
# smoke (plus characterisation / property from its flags) and the re-run list
# gives metamorphic; these rows add the rest.
# -----------------------------------------------------------------------------
# lsdtt-synthetic-networks
lsdtt_evidence(smoke_synthetic_networks_island system smoke exact_solution)        # 890/1600 land pixels by construction
lsdtt_evidence(smoke_synthetic_networks_polygon_area system smoke exact_solution)  # 14x14 = ceil(3000/250)+2
lsdtt_evidence(smoke_synthetic_networks_equilibrate system smoke property characterisation)  # "Converged at cycle 5" measured
lsdtt_evidence(rigorous_synthetic_networks_dem_band system smoke characterisation) # measured bands, "not an oracle"

# lsdtt-raster-preprocessing
lsdtt_evidence(smoke_raster_preprocessing_replace_pixels system smoke exact_solution)
lsdtt_evidence(smoke_raster_preprocessing_breach_least_cost system smoke exact_solution)  # hand-built 7x7 bowl
lsdtt_evidence(rigorous_raster_preprocessing_breach_least_cost system smoke characterisation exact_solution)
lsdtt_evidence(rigorous_raster_preprocessing_breach_lc_whitebox_edge system smoke characterisation exact_solution)
foreach(_t IN ITEMS line polygon footprint back)
    # Frame sizes and origins pinned; the geometry is asserted in
    # unit_oblique_frame_preprocess.
    lsdtt_evidence(smoke_raster_preprocessing_rotate_${_t} system smoke characterisation)
endforeach()
lsdtt_evidence(smoke_basin_extraction_outlets_off_network system smoke characterisation)

# Interior nodata and min_elevation: counts follow from the fixture.
lsdtt_evidence(smoke_interior_nodata_fill_drainage_analysis system smoke property exact_solution)
lsdtt_evidence(smoke_interior_nodata_sink_drainage_analysis system smoke exact_solution)  # 96 % = 24/25
lsdtt_evidence(smoke_interior_nodata_edge_not_filled_drainage_analysis system smoke exact_solution)
lsdtt_evidence(smoke_interior_nodata_fill_preprocessing system smoke property exact_solution)
lsdtt_evidence(smoke_nodata_min_elevation_negative_preprocessing system smoke exact_solution)  # 7 cells, min -79
lsdtt_evidence(smoke_nodata_edge_unfilled system smoke exact_solution)

# True-north bearings on planes dipping to a known azimuth.
lsdtt_evidence(smoke_surface_metrics_true_north_frame system smoke exact_solution)  # 200
lsdtt_evidence(smoke_surface_metrics_grid_north_frame system smoke exact_solution)  # 200 + 39.87
lsdtt_evidence(smoke_surface_metrics_true_north_utm system smoke exact_solution)
lsdtt_evidence(smoke_feature_finder_true_north_structure_tensor system smoke exact_solution)  # 110
lsdtt_evidence(smoke_ridgeline_true_north system smoke exact_solution)

# lsdtt-drainage-analysis and TWI/STI
lsdtt_evidence(rigorous_drainage_analysis_canonical_threshold system smoke property)  # HAND >= 0, 0 on channels
lsdtt_evidence(rigorous_drainage_analysis_legacy_threshold system smoke property)
lsdtt_evidence(rigorous_drainage_sti_deprecated_names system smoke characterisation metamorphic)  # canonical pins
lsdtt_evidence(rigorous_drainage_sti_reprojected_layers system smoke characterisation metamorphic) # UTM sibling pins

# lsdtt-distributed-hydrology
lsdtt_evidence(smoke_distributed_hydrology_run system smoke characterisation property metamorphic)  # TWI = drainage-analysis pin
lsdtt_evidence(smoke_distributed_hydrology_hype_outlets system smoke property exact_solution)       # 540 = 9 x 60
lsdtt_evidence(smoke_distributed_hydrology_hype_snap_square system smoke property characterisation)
lsdtt_evidence(smoke_distributed_hydrology_hype_snap_downstream system smoke property characterisation)
lsdtt_evidence(smoke_distributed_hydrology_hype_snap_largest_area system smoke property exact_solution)
lsdtt_evidence(smoke_distributed_hydrology_hype_class_holes system smoke property exact_solution)
lsdtt_evidence(smoke_distributed_hydrology_hype_short_run system smoke exact_solution)  # 30 rows for 30 steps

# Channels, sections, widths, centrelines
lsdtt_evidence(smoke_stream_burning_lcp_burn system smoke characterisation)
lsdtt_evidence(smoke_cross_sections system smoke characterisation metamorphic)
lsdtt_evidence(smoke_bankfull_width system smoke property characterisation metamorphic)
lsdtt_evidence(rigorous_cross_sections_legacy_points_key system smoke characterisation)
lsdtt_evidence(rigorous_bankfull_width_legacy_points_key system smoke characterisation)
lsdtt_evidence(smoke_marsh_topography system smoke characterisation metamorphic)
lsdtt_evidence(smoke_water_centreline system smoke characterisation metamorphic)
lsdtt_evidence(smoke_advection_diffusion_resume system smoke characterisation metamorphic)  # resumed = uninterrupted hash
lsdtt_evidence(smoke_slope_area system smoke characterisation metamorphic)
lsdtt_evidence(smoke_slope_area_d8 system smoke characterisation)
lsdtt_evidence(smoke_slope_area_no_flowinfo system smoke characterisation metamorphic)    # same row, flowinfo off
lsdtt_evidence(smoke_slope_area_clip system smoke characterisation)
lsdtt_evidence(smoke_outlet_matching_identity system smoke property characterisation metamorphic)
lsdtt_evidence(smoke_outlet_matching_resampled system smoke characterisation)
lsdtt_evidence(smoke_outlet_matching_edge_only system smoke characterisation)
lsdtt_evidence(smoke_water_line_repair system smoke exact_solution metamorphic)  # one fixture case per rule
foreach(_t IN ITEMS mutual_best_on mutual_best_off wide_road_gap separate_roads)
    lsdtt_evidence(smoke_water_line_repair_${_t} system smoke exact_solution)
endforeach()
lsdtt_evidence(smoke_swath_profile_plane system smoke exact_solution metamorphic)  # analytic answers

# Mountain front, feature finder, channel network analysis, raster calculator
lsdtt_evidence(rigorous_mountain_front_legacy_keys system smoke characterisation metamorphic)
lsdtt_evidence(rigorous_feature_finder_risers_staircase system smoke property exact_solution)
lsdtt_evidence(rigorous_feature_finder_risers_mapped system smoke property characterisation)
lsdtt_evidence(rigorous_feature_finder_risers_mapped_30m system smoke characterisation)
lsdtt_evidence(rigorous_feature_finder_flow_discrimination_canonical system smoke characterisation property)
lsdtt_evidence(rigorous_feature_finder_flow_discrimination_legacy_keys system smoke characterisation property metamorphic)
lsdtt_evidence(smoke_channel_network_analysis_area_scaling system smoke characterisation metamorphic)
lsdtt_evidence(rigorous_channel_network_analysis_bearing_statistics_reference system smoke exact_solution)  # W1 = 90 deg
lsdtt_evidence(rigorous_channel_network_analysis_bearing_statistics_region system smoke exact_solution)     # W1 = 0
# The fixture's shift is (4, -3, 0.5) m; the test pins the recovered dx=4.01.
lsdtt_evidence(smoke_raster_calculator_coregistration system smoke characterisation)

# lsdtt-divide-migration: divide counts follow from the three-valley topology.
lsdtt_evidence(smoke_divide_migration_no_clusters system smoke exact_solution characterisation)
lsdtt_evidence(rigorous_divide_migration_legacy_threshold system smoke characterisation)
lsdtt_evidence(smoke_divide_migration_clusters system smoke exact_solution characterisation metamorphic)
lsdtt_evidence(smoke_divide_migration_clusters_vector system smoke exact_solution metamorphic)
foreach(_t IN ITEMS merge_unmerged merge_fragments pair_summary)
    lsdtt_evidence(smoke_divide_migration_${_t} system smoke exact_solution)
endforeach()

# lsdtt-lem (smoke_lem_roering_implicit_strip drives the locked hillslope solver)
lsdtt_evidence(smoke_lem_erodibility_map_area_conversion system smoke exact_solution)  # 10^0.45
lsdtt_evidence(smoke_lem_grid_from_erodibility_raster_footprint system smoke exact_solution)  # 75 % footprint
lsdtt_evidence(smoke_lem_roering_implicit_strip system smoke property characterisation)
lsdtt_evidence(smoke_lem_hillslope_on_loaded_dem system smoke property)
lsdtt_evidence(smoke_lem_hillslope_on_loaded_dem_critical_slope system smoke property)
lsdtt_evidence(smoke_lem_equilibration_divide_step system smoke property)

# lsdtt-internal-drainage (crater fixture)
lsdtt_evidence(smoke_internal_drainage_both system smoke exact_solution metamorphic)
lsdtt_evidence(smoke_internal_drainage_route system smoke exact_solution)
lsdtt_evidence(smoke_internal_drainage_moved_route system smoke characterisation)

# lsdtt-ridgeline-analysis, lsdtt-hillslope-channel-coupling
lsdtt_evidence(smoke_ridgeline_paired_flowpaths_sl_on system smoke characterisation)
lsdtt_evidence(rigorous_ridgeline_legacy_threshold system smoke characterisation)
lsdtt_evidence(smoke_ridgeline_paired_flowpaths_sl_off system smoke characterisation)
lsdtt_evidence(smoke_ridgeline_paired_flowpaths_all_basins system smoke characterisation)
lsdtt_evidence(rigorous_hillslope_channel_coupling_floodplain_legacy_keys system smoke characterisation metamorphic)
lsdtt_evidence(rigorous_hillslope_channel_coupling_legacy_threshold system smoke characterisation)

# lsdtt-chi-analysis and lithology
lsdtt_evidence(smoke_chi_analysis_erodibility_calibration system smoke characterisation)
lsdtt_evidence(unit_channel_gradient unit exact_solution)  # floor gradient 0.01 exact
lsdtt_evidence(smoke_channel_extraction_slope_points system smoke exact_solution)
lsdtt_evidence(smoke_channel_extraction_point_attributes system smoke exact_solution)
lsdtt_evidence(smoke_channel_extraction_point_attributes_fgb system smoke)
lsdtt_evidence(unit_channel_nearest_point_snap unit exact_solution property)  # vs brute force
lsdtt_evidence(unit_chi_knickpoints unit exact_solution)  # hand-built fit, known k_sn / areas
lsdtt_evidence(smoke_chi_tributaries_knickpoints system smoke characterisation)
lsdtt_evidence(smoke_chi_lithology_ksn system smoke exact_solution)  # 9600 px per unit, k_sn 20 / 5 bins
lsdtt_evidence(rigorous_chi_lithology_observation_only system smoke characterisation metamorphic)
lsdtt_evidence(smoke_chi_lithology_crosswalk system smoke characterisation)
lsdtt_evidence(smoke_chi_lithology_relative_erodibility_k_raster system smoke exact_solution)
lsdtt_evidence(smoke_channel_network_analysis_lithology system smoke exact_solution)  # 24000 - 364

# Flat resolution
lsdtt_evidence(rigorous_fill_flat_resolution_valley_float_step system smoke exact_solution)  # 48 x 2^-17 m
lsdtt_evidence(rigorous_fill_flat_resolution_valley_gradient system smoke exact_solution)
lsdtt_evidence(rigorous_fill_flat_resolution_knob_rounds system smoke characterisation)
foreach(_t IN ITEMS
        smoke_chi_analysis_flat_resolution_gradient
        smoke_chi_analysis_flat_resolution_filled_float_step
        smoke_channel_extraction_flat_resolution_gradient
        smoke_drainage_analysis_flat_resolution_float_step
        smoke_basin_extraction_flat_resolution_gradient
        smoke_basin_extraction_flat_resolution_filled_float_step
        smoke_basin_stats_flat_resolution_float_step
        smoke_basin_surface_metrics_flat_resolution_gradient
        smoke_channel_network_analysis_flat_resolution_float_step
        smoke_cross_sections_flat_resolution_gradient
        smoke_bankfull_width_flat_resolution_gradient_warning
        smoke_bankfull_width_flat_resolution_float_step_warning
        smoke_distributed_hydrology_flat_resolution_hype_gradient
        smoke_distributed_hydrology_flat_resolution_topmodel_float_step
        smoke_fluvial_surfaces_flat_resolution_float_step
        smoke_ridgeline_flat_resolution_gradient
        smoke_valley_centreline_flat_resolution_float_step
        rigorous_drainage_analysis_flat_resolution_real_terrain_float_step
        rigorous_drainage_analysis_flat_resolution_real_terrain_gradient)
    # Report lines ("Flat regions: 1, flat cells: 198 ...", round counts) pinned.
    lsdtt_evidence(${_t} system smoke characterisation)
endforeach()

# Spectral analysis, fetch, cosmo-landslides, channel extraction
lsdtt_evidence(smoke_spectral_analysis_valley_spacing system smoke characterisation metamorphic)
lsdtt_evidence(smoke_spectral_analysis_significant_2d system smoke characterisation)
lsdtt_evidence(network_fetch_raster_linear_feature_cli system smoke property characterisation)
lsdtt_evidence(smoke_fetch_lithology_offline system smoke characterisation metamorphic)
lsdtt_evidence(smoke_fetch_river_lines_local system smoke exact_solution)
lsdtt_evidence(smoke_fetch_river_lines_to_stream_burning system metamorphic characterisation)
lsdtt_evidence(unit_stream_burning_source_seed unit exact_solution)
lsdtt_evidence(unit_stream_burning_loop_labels unit exact_solution)
lsdtt_evidence(unit_stream_burning_corridor unit exact_solution property metamorphic)
lsdtt_evidence(smoke_stream_burning_corridor system smoke property characterisation)
lsdtt_evidence(smoke_stream_burning_corridor_x2 system smoke property)
lsdtt_evidence(smoke_stream_burning_corridor_min_elevation_below_land system smoke metamorphic)
lsdtt_evidence(smoke_stream_burning_corridor_min_elevation_sea system smoke property)
lsdtt_evidence(smoke_stream_burning_corridor_nodata_pieces system smoke property)
lsdtt_evidence(unit_water_line_repair_point_snapper unit exact_solution)
lsdtt_evidence(smoke_fetch_erosion_rates_offline system smoke exact_solution property metamorphic)  # = committed response
lsdtt_evidence(smoke_fetch_erosion_rates_offline_lsdn system smoke exact_solution)
lsdtt_evidence(smoke_cosmo_landslides_sample_ensemble system smoke characterisation)
lsdtt_evidence(smoke_channel_extraction_perimeter_sources_tributary system smoke exact_solution)
lsdtt_evidence(smoke_channel_extraction_perimeter_sources_two_network system smoke exact_solution)

# lsdtt-fluvial-surfaces, lsdtt-feature-finder risers, lsdtt-valley-centreline,
# lsdtt-flat-surface-finder
lsdtt_evidence(smoke_fluvial_surfaces_generators system smoke exact_solution)  # ids as configured
lsdtt_evidence(smoke_feature_finder_risers_connect system smoke characterisation)
lsdtt_evidence(smoke_fluvial_surfaces_band_then_label system smoke property)  # band value set
lsdtt_evidence(smoke_fluvial_surfaces system smoke property metamorphic)
lsdtt_evidence(smoke_fluvial_surfaces_fitted_surface system smoke property)
lsdtt_evidence(rigorous_fluvial_surfaces_clubb2017 system smoke characterisation property)
lsdtt_evidence(rigorous_fluvial_surfaces_legacy_keys system smoke characterisation metamorphic)
lsdtt_evidence(rigorous_fluvial_surfaces_clubb2017_floodplain_only system smoke characterisation metamorphic)
lsdtt_evidence(rigorous_fluvial_surfaces_clubb2017_second_pass system smoke characterisation metamorphic)
lsdtt_evidence(rigorous_fluvial_surfaces_relative_elevation system smoke characterisation)
lsdtt_evidence(rigorous_valley_centreline_relative_elevation_whole_dem system smoke characterisation)
lsdtt_evidence(rigorous_valley_centreline_terrace_profile system smoke characterisation)
lsdtt_evidence(rigorous_fluvial_surfaces_clubb2017_low_relief_auto system smoke characterisation metamorphic)
lsdtt_evidence(rigorous_fluvial_surfaces_patch_polygons system smoke exact_solution)  # 4 levels x 8
lsdtt_evidence(rigorous_fluvial_surfaces_methods_list system smoke characterisation metamorphic)
lsdtt_evidence(rigorous_fluvial_surfaces_fill_admitted_holes system smoke characterisation)
lsdtt_evidence(rigorous_fluvial_surfaces_training system smoke characterisation)
lsdtt_evidence(rigorous_fluvial_surfaces_mapped system smoke characterisation)
lsdtt_evidence(rigorous_fluvial_surfaces_mapped_gpkg system smoke characterisation metamorphic)   # = FlatGeobuf sha
lsdtt_evidence(rigorous_fluvial_surfaces_mapped_4326 system smoke characterisation metamorphic)
lsdtt_evidence(rigorous_fluvial_surfaces_mapped_extent system smoke exact_solution)  # 30000 of 60000
lsdtt_evidence(rigorous_fluvial_surfaces_mapped_extent_multi system smoke characterisation)
lsdtt_evidence(rigorous_fluvial_surfaces_clubb2017_order1 system smoke characterisation)
lsdtt_evidence(rigorous_fluvial_surfaces_levels_hand system smoke characterisation property)
lsdtt_evidence(rigorous_flat_surface_finder_relief_qq system smoke characterisation metamorphic)  # = mountain-front QQ
lsdtt_evidence(rigorous_mountain_front_relief_qq system smoke characterisation metamorphic)

# -----------------------------------------------------------------------------
# System tests outside the harness
# -----------------------------------------------------------------------------
# cmake -E compare_files / compare_output_sets.cmake between two CLI runs.
foreach(_t IN ITEMS
        rigorous_drainage_analysis_threshold_equivalence
        rigorous_chi_erodibility_calibration_units_invariant
        rigorous_chi_erodibility_calibration_observation_only
        rigorous_chi_lithology_bootstrap_per_mn
        rigorous_chi_lithology_relative_erodibility_observation_only
        rigorous_chi_lithology_relative_erodibility_per_mn
        rigorous_channel_network_analysis_lithology_observation_only)
    lsdtt_evidence(${_t} system metamorphic)
endforeach()

# lsdtt-lem rasters: lsdtt_add_lem_k_compare records exact_solution or
# metamorphic. The two lithology routes below reach the own-grid expected
# raster from a different input, so they are both.
lsdtt_evidence(regression_lem_lithology_erodibility_code_raster system exact_solution metamorphic)
lsdtt_evidence(regression_lem_lithology_erodibility_crosswalk system exact_solution metamorphic)
lsdtt_evidence(regression_lem_synthetic_default_crs system exact_solution)   # EPSG 32758, grid centre
lsdtt_evidence(regression_lem_synthetic_supplied_crs system exact_solution)  # EPSG 32630

# Known answers recovered from CLI output by a C++ checker.
lsdtt_evidence(rigorous_chi_erodibility_calibration_n system exact_solution)         # n = 1
lsdtt_evidence(rigorous_chi_erodibility_calibration_fixed_n_K system exact_solution) # K = 1e-5
lsdtt_evidence(rigorous_chi_erodibility_calibration_lithology_units_west system exact_solution)
lsdtt_evidence(rigorous_chi_erodibility_calibration_lithology_units_east system exact_solution)
lsdtt_evidence(rigorous_chi_lithology_ksn_ratio system exact_solution)                # ratio 4
lsdtt_evidence(rigorous_chi_lithology_relative_erodibility_ratio system exact_solution)
lsdtt_evidence(rigorous_chi_lithology_relative_erodibility_reference system exact_solution)

# Cosmogenic CLI scripts (locked CRN/CAIRN core: labels only)
lsdtt_evidence(rigorous_cosmo_cli_bridge system metamorphic exact_solution)  # forward then inverse returns the rate
lsdtt_evidence(rigorous_cosmo_multi_accum_ratio system property)              # ratio = Al26/Be10 quotient
lsdtt_evidence(rigorous_cosmo_column_tool system exact_solution property metamorphic)
lsdtt_evidence(regression_step_change_disk_roundtrip system metamorphic)      # in-memory = written and read
lsdtt_evidence(regression_landslide_scar_coupling system metamorphic)         # unfailed ground unchanged
lsdtt_evidence(regression_cosmo_multi_landslide_coupling system smoke property)
lsdtt_evidence(regression_cosmo_oracle system characterisation)               # file hashes, cosmo_oracle.csv
lsdtt_evidence(regression_cosmo_oracle_heavy system characterisation)
# Real-world tier. run_cosmo_validation.py --check asserts only that the
# inverted rate is within 2 % of the frozen lsdtt3 summary; the published
# OCTOPUS rate and the difference from it are written, not asserted. Likewise
# the production-isolation run compares with its own frozen artefact.
foreach(_t IN ITEMS dethier2014 clapp2002 dibiase2010
                    dethier2014_lsdn clapp2002_lsdn dibiase2010_lsdn
                    production_isolation lsdn_pixel_vs_lumped)
    lsdtt_evidence(cosmo_realworld_${_t} system characterisation)
endforeach()

# lsdtt-fluvial-surfaces / valley-centreline / flat-surface-finder scripts
lsdtt_evidence(rigorous_fluvial_surfaces_label_band_identity system property)  # joins consistent across outputs
lsdtt_evidence(rigorous_fluvial_surfaces_clubb2017_joins system property)
lsdtt_evidence(rigorous_fluvial_surfaces_empty_result system property)
lsdtt_evidence(rigorous_flat_surface_finder_admission_modes system property metamorphic)
lsdtt_evidence(rigorous_fluvial_surfaces_cross_section_levels_workflow system exact_solution property)  # 8/16/24/32 m
lsdtt_evidence(rigorous_fluvial_surfaces_whole_mask_segmentation system exact_solution metamorphic)
lsdtt_evidence(rigorous_valley_centreline_legacy_threshold system smoke)   # legacy key runs and warns
lsdtt_evidence(rigorous_valley_centreline_legacy_points_key system smoke)
lsdtt_evidence(smoke_fluvial_surfaces_pipeline system smoke)
lsdtt_evidence(smoke_fluvial_surfaces_hand_input system metamorphic)        # supplied HAND = computed, by checksum
lsdtt_evidence(smoke_multipart_cli_sweep system smoke metamorphic)
lsdtt_evidence(smoke_flat_resolution_real_terrain_routing system smoke property)

# network_* CLI runs: the harness flags give smoke and the size/minimum checks.
# network_raster_fetch is a C++ driver calling the fetch routines in process.
lsdtt_evidence(network_raster_fetch unit smoke property metamorphic characterisation)

# Configure, banners, conventions, packaging
lsdtt_evidence(configure_contract_native_default system conformance)
lsdtt_evidence(configure_contract_geo_without_python system negative)
lsdtt_evidence(configure_contract_geo_missing_modules system negative)
lsdtt_evidence(configure_contract_gdal_programs_missing system negative)
lsdtt_evidence(unit_synthetic_surfaces_help_defaults_exact system conformance)
lsdtt_evidence(unit_synthetic_cli_banners system conformance)
lsdtt_evidence(unit_test_env_isolation system conformance)   # runs a CLI --version
lsdtt_evidence(unit_cli_parameter_conventions system conformance)  # -h json of every CLI
lsdtt_evidence(unit_cli_reference_is_current system conformance)
# Builds the bindings and runs python_bindings/tests; their own evidence is
# that of the pytest suite, which is not labelled here.
lsdtt_evidence(python_bindings_smoke system smoke)

# -----------------------------------------------------------------------------
# Unit tests: flow, depressions, internal drainage, channels
# -----------------------------------------------------------------------------
lsdtt_evidence(unit_flow_properties unit property exact_solution metamorphic)  # P1-P11
lsdtt_evidence(unit_flow_d8_identity unit characterisation property metamorphic)  # pins from 7e81a39
lsdtt_evidence(unit_flat_resolution unit differential exact_solution property metamorphic negative characterisation)  # frozen RichDEM raises
lsdtt_evidence(unit_flat_resolution_real_terrain unit characterisation property metamorphic)  # pinned rounds; routes
lsdtt_evidence(unit_breach_least_cost unit differential property negative)  # frozen Whitebox 2.0.6
lsdtt_evidence(unit_breach_least_cost_real_terrain unit differential property)
lsdtt_evidence(unit_flow_snapping unit exact_solution negative metamorphic)
lsdtt_evidence(unit_mountain_front_outlet_pruning unit exact_solution negative)
# Asserts a complexity exponent from timings: a property, not performance.
lsdtt_evidence(unit_depression_hierarchy unit property exact_solution)
lsdtt_evidence(unit_depression_hierarchy_oracle unit differential)  # Barnes' upstream code
lsdtt_evidence(unit_internal_drainage unit exact_solution property negative)
lsdtt_evidence(unit_internal_drainage_detect unit exact_solution metamorphic)
lsdtt_evidence(unit_drainage_injection_weights unit characterisation)  # removed loops kept verbatim
lsdtt_evidence(unit_boundary_conditions unit characterisation)
lsdtt_evidence(unit_junction_graph_orders unit characterisation)
foreach(_t IN ITEMS circular summaries comparisons junction_angles trellis grouping)
    lsdtt_evidence(unit_channel_bearing_statistics_${_t} unit exact_solution)
endforeach()
lsdtt_evidence(unit_channel_bearing_statistics_isotropic_null unit property)
lsdtt_evidence(unit_channel_profile_attributes unit exact_solution property)  # P2 one cell per cell; P3/P4
lsdtt_evidence(unit_link_length_periodic unit exact_solution)
lsdtt_evidence(unit_channel_threshold_units unit exact_solution)
lsdtt_evidence(unit_drainage_area_exceedance_ols unit characterisation)  # old ols_fit vs new
lsdtt_evidence(unit_lithology_channel_network_stats unit exact_solution)
lsdtt_evidence(unit_chi_lithology_channel_stats unit exact_solution property metamorphic)
lsdtt_evidence(unit_chi_erodibility_calibration unit differential exact_solution metamorphic)  # Cantrell, scipy.odr
lsdtt_evidence(unit_divide_migration unit exact_solution)
# Brute-force pair references written in the driver; K invalid throws;
# idempotent on a ramp; thread counts; no nodata in the mask.
lsdtt_evidence(unit_divide_diagnostics unit differential property metamorphic negative)
lsdtt_evidence(unit_paired_flowpaths unit exact_solution)
lsdtt_evidence(unit_hilltop_flow_routing unit exact_solution)
lsdtt_evidence(unit_hilltop_flow_routing_integration unit property)
lsdtt_evidence(unit_streamlines unit exact_solution)

# -----------------------------------------------------------------------------
# Unit tests: distributed hydrology and HYPE (HYPE 5.36.0 oracle files are
# frozen outputs of another program: differential)
# -----------------------------------------------------------------------------
lsdtt_evidence(unit_topmodel_column unit exact_solution property)
lsdtt_evidence(unit_distributed_hydrology_forcing unit exact_solution negative)
lsdtt_evidence(unit_hype_parameters unit exact_solution negative property)
lsdtt_evidence(unit_hype_column unit differential property)
lsdtt_evidence(unit_hype_routing unit differential exact_solution property characterisation)
lsdtt_evidence(unit_hype_model unit differential property metamorphic negative)
lsdtt_evidence(unit_hype_raster_setup unit exact_solution negative property)
lsdtt_evidence(unit_hype_thread_invariance unit metamorphic)
lsdtt_evidence(unit_subbasin_partition unit exact_solution)
lsdtt_evidence(unit_channel_extent unit exact_solution metamorphic)
lsdtt_evidence(unit_hydro_routing unit property metamorphic)
lsdtt_evidence(unit_canopy_interception unit exact_solution property)
lsdtt_evidence(unit_hydro_roughness unit exact_solution property metamorphic)
lsdtt_evidence(unit_landcover_scenario unit exact_solution)
# test_efficiency asserts the recurrence is faster than pow-and-factorial: a
# property of the method, not a reported timing.
lsdtt_evidence(unit_infiltration_excess unit exact_solution property metamorphic negative)
lsdtt_evidence(unit_topmodel_dynamic unit exact_solution property)
lsdtt_evidence(unit_hydro_velocity_field unit exact_solution property metamorphic)
lsdtt_evidence(unit_hydro_velocity_manning_formula unit characterisation)
lsdtt_evidence(unit_wetness_index unit exact_solution)  # ln(100); reference float evaluation
lsdtt_evidence(unit_distributed_hydrology_argument_checks unit characterisation negative)
lsdtt_evidence(unit_giuh_discharge unit property metamorphic)

# -----------------------------------------------------------------------------
# Unit tests: surface metrics, image operations, spectral, raster utilities
# -----------------------------------------------------------------------------
lsdtt_evidence(unit_surface_metrics unit exact_solution property metamorphic characterisation)
lsdtt_evidence(unit_surface_metrics_polyfit_metric_table unit property)
lsdtt_evidence(unit_image_filters unit differential exact_solution property)
lsdtt_evidence(unit_image_ops unit differential characterisation property metamorphic)  # OpenCV 4.13 hashes
lsdtt_evidence(unit_spectral_analysis_significant_2d_threshold unit exact_solution negative)
foreach(_t IN ITEMS align bearings edges)
    lsdtt_evidence(unit_spectral_analysis_significant_2d_${_t} unit exact_solution)
endforeach()
lsdtt_evidence(unit_spectral_analysis_significant_2d_seed unit metamorphic)
lsdtt_evidence(unit_vertical_units unit exact_solution property negative)
# Consolidation characterisation: the removed copies kept verbatim.
foreach(_t IN ITEMS raster_cell_predicates raster_sampling grid_neighbours
                    connected_components grid_geodesic raster_central_difference
                    order_statistics)
    lsdtt_evidence(unit_${_t} unit characterisation)
endforeach()
lsdtt_evidence(unit_raster_vectorise_polygons unit exact_solution)
lsdtt_evidence(unit_raster_resampling unit exact_solution property)
lsdtt_evidence(unit_raster_hierarchy unit exact_solution)  # dyadic layout and activation cases
lsdtt_evidence(unit_raster_jenks_classification unit exact_solution property negative)
lsdtt_evidence(unit_raster_active_masks unit exact_solution negative)
lsdtt_evidence(unit_raster_transform unit exact_solution)
lsdtt_evidence(unit_multipart_geometry unit exact_solution)
lsdtt_evidence(unit_oblique_frame unit differential metamorphic)  # Vincenty in the test, pyproj
lsdtt_evidence(unit_oblique_frame_feature unit metamorphic negative)
lsdtt_evidence(unit_oblique_frame_file unit property negative)
lsdtt_evidence(unit_oblique_frame_warp unit property)
lsdtt_evidence(unit_oblique_frame_preprocess unit metamorphic property)
lsdtt_evidence(unit_true_north unit differential exact_solution property)
lsdtt_evidence(unit_crs_semantic_matching unit exact_solution)
lsdtt_evidence(unit_vector_layer_transform unit characterisation)
lsdtt_evidence(unit_gdal_init_idempotent unit property)
# The Edinburgh coordinates are MEASURED and pinned (test_gdal_init.cpp).
lsdtt_evidence(unit_gdal_init_pinned unit characterisation)
lsdtt_evidence(unit_gdal_init_fallthrough unit property)
lsdtt_evidence(unit_gdal_format_support unit property)  # every required driver round-trips
lsdtt_evidence(unit_raster_loader_geotransform_guard unit negative)
lsdtt_evidence(unit_gdal_init_boundary unit conformance)  # headers confined to src/

# -----------------------------------------------------------------------------
# Unit tests: LEM, hillslope solvers (LOCKED numerics: labels only), synthetic
# networks
# -----------------------------------------------------------------------------
lsdtt_evidence(unit_lem_erodibility unit metamorphic exact_solution negative property)
lsdtt_evidence(rigorous_lem_erodibility_analytic unit exact_solution)
lsdtt_evidence(unit_lem_timestep unit property)
lsdtt_evidence(unit_lem_channel_coupling unit property)
lsdtt_evidence(unit_lem_forcing_history unit exact_solution)
lsdtt_evidence(unit_lem_lithology_erodibility unit exact_solution negative)
lsdtt_evidence(unit_thread_count_invariance unit metamorphic property)
lsdtt_evidence(unit_hillslope_implicit unit exact_solution metamorphic property)        # locked
lsdtt_evidence(unit_hillslope_newton unit metamorphic exact_solution property)          # locked
lsdtt_evidence(unit_hillslope_steady_state unit exact_solution property metamorphic)    # locked
# Analytic axisymmetric equilibrium; Jacobian against finite differences
# computed in the test; rows sum to zero.
lsdtt_evidence(unit_hillslope_perron_qimp unit differential exact_solution property)    # locked
lsdtt_evidence(unit_hillslope_ren_modular_implicit unit differential exact_solution property negative)  # locked; published RMSE, analytic axisymmetric, row sum, supercritical rejected
# 1-D ramp profile z(receiver)+dx*S_c; refusals; masked cells unchanged;
# m/n = 0 equals the fixed-gradient path.
lsdtt_evidence(unit_hillslope_critical_slope unit exact_solution property metamorphic negative)  # locked
lsdtt_evidence(unit_hillslope_transport_bc unit metamorphic property)                   # locked
lsdtt_evidence(unit_hillslope_thread_count_invariance unit metamorphic)                 # locked
lsdtt_evidence(unit_hillslope_length_statistic unit exact_solution)
lsdtt_evidence(unit_hillslope_length_search unit property exact_solution)
lsdtt_evidence(unit_hillslope_peclet unit property exact_solution)
# Writes reference fields; asserts convergence, translation invariance of the
# periodic case, the uplift volume and slope-ratio bounds (no stored field is
# compared).
lsdtt_evidence(unit_cellular_hillslope_references unit metamorphic exact_solution property)  # locked
lsdtt_evidence(unit_cellular_hillslope_constitutive unit property exact_solution)
lsdtt_evidence(unit_cellular_hillslope_solver unit metamorphic property)                # locked
lsdtt_evidence(unit_cellular_hillslope_common_limit unit metamorphic)                  # locked
lsdtt_evidence(unit_hillslope_particle_io unit property)  # write/read round trip
lsdtt_evidence(unit_landscape_equilibration unit property)
# C1/C2 inequalities; C3 outputs unmoved by a threshold; C6 conventions agree.
lsdtt_evidence(unit_chi_channel_head unit property metamorphic)
lsdtt_evidence(unit_synthetic_networks unit property differential metamorphic)
lsdtt_evidence(unit_synthetic_networks_sweep unit smoke)  # table reported, never asserted

# -----------------------------------------------------------------------------
# Unit tests: cosmogenic (LOCKED CRN/CAIRN core: labels only) and landslides
# -----------------------------------------------------------------------------
lsdtt_evidence(unit_cosmo_physics unit exact_solution property metamorphic differential characterisation)  # locked
lsdtt_evidence(unit_cosmo_units unit characterisation exact_solution)   # locked-adjacent
lsdtt_evidence(unit_cosmo_ratio unit characterisation property)
lsdtt_evidence(unit_cosmo_lsdn_grid_lookup unit characterisation)       # old bilinear_regular verbatim
lsdtt_evidence(unit_cosmo_raster_axis_order unit metamorphic)
lsdtt_evidence(unit_quartz_fertility unit exact_solution)
lsdtt_evidence(unit_step_change unit exact_solution property)
lsdtt_evidence(unit_landslide_size unit exact_solution)
lsdtt_evidence(unit_niemi_ensemble unit exact_solution property)  # N1 control; N3-N6 trends
lsdtt_evidence(unit_landslide_step_change unit exact_solution property)
lsdtt_evidence(unit_landslide_placement unit exact_solution property)
lsdtt_evidence(unit_landslide_event unit exact_solution property)
lsdtt_evidence(unit_landslide_crn unit property)
lsdtt_evidence(unit_landslide_series unit exact_solution property)
lsdtt_evidence(unit_landslide_rng unit characterisation)
lsdtt_evidence(unit_splitmix_copies unit characterisation)
lsdtt_evidence(unit_cosmo_validation_harness unit conformance)   # harness self-test, no CLI
lsdtt_evidence(unit_cosmo_riversand_harness unit conformance)

# -----------------------------------------------------------------------------
# Unit tests: fluvial surfaces, features, fetch, core
# -----------------------------------------------------------------------------
lsdtt_evidence(unit_feature_finder_risers unit exact_solution property)
lsdtt_evidence(unit_fluvial_surface_method_comparison unit exact_solution property)
lsdtt_evidence(unit_valley_width_network unit exact_solution)
lsdtt_evidence(unit_valley_relative_elevation unit exact_solution)
lsdtt_evidence(unit_perimeter_sources unit exact_solution metamorphic)
lsdtt_evidence(unit_fluvial_surface_training unit exact_solution)
lsdtt_evidence(unit_fluvial_surfaces_fixture unit exact_solution property)  # fixture self-checks, staircase
lsdtt_evidence(unit_fluvial_surface_ladder unit property characterisation)  # T5 legacy cut bit for bit
# Cone gradient recovered within 10 %; fewer components than HAND; refusal.
lsdtt_evidence(unit_fluvial_surface_fitted_field unit exact_solution metamorphic negative)
lsdtt_evidence(unit_fluvial_surface_perona_malik unit property metamorphic negative)
lsdtt_evidence(unit_fluvial_surface_bootstrap_rbar unit property exact_solution)
lsdtt_evidence(unit_fluvial_surface_slope_dependent_ev unit property)
lsdtt_evidence(unit_fluvial_surface_channel_order_swath unit differential)  # brute force in the test
lsdtt_evidence(unit_fluvial_surface_floodplain_thresholds unit exact_solution negative)
lsdtt_evidence(unit_fluvial_surface_cross_section_levels unit exact_solution)
lsdtt_evidence(unit_fluvial_surface_coherence_diffusion unit property metamorphic)
lsdtt_evidence(unit_fluvial_surfaces_trough_upstream unit property)  # > 200 cells excluded
lsdtt_evidence(unit_fluvial_surfaces_band_then_label unit exact_solution property)  # one label per known tread
lsdtt_evidence(unit_fluvial_surface_threshold_selection unit exact_solution)
lsdtt_evidence(unit_basin_outline_metrics unit exact_solution)
lsdtt_evidence(unit_raster_fetch unit exact_solution negative)
lsdtt_evidence(unit_dem_grouping unit exact_solution)
lsdtt_evidence(unit_fetch_lithology unit exact_solution)
lsdtt_evidence(unit_fetch_erosion_rates unit exact_solution)
lsdtt_evidence(unit_fetch_river_lines_prepare unit exact_solution)
lsdtt_evidence(unit_fetch_river_lines_adapters unit exact_solution differential metamorphic negative)
lsdtt_evidence(unit_coastal_mask unit exact_solution)
# Security properties of credential handling; resolution order by hand.
lsdtt_evidence(unit_secret unit exact_solution property)
foreach(_t IN ITEMS cli_helper_string_utils cli_helper_formats csv_table_quoting
                    param_parser_help common_parameters param_parser_strict
                    run_record)
    lsdtt_evidence(unit_${_t} unit exact_solution)  # hand-written expected strings
endforeach()
lsdtt_evidence(unit_csv_table_field_parsers unit exact_solution negative)
lsdtt_evidence(unit_point_input unit exact_solution metamorphic negative)
lsdtt_evidence(unit_point_attribute_output unit exact_solution metamorphic negative)
# Removed and misspelt keys are refused (negative); a deprecated alias or an
# unrelated key is accepted with a warning (conformance, rule R1 of the
# pytest-marker plan); the exact messages and carried values are pinned by
# hand (exact_solution).
lsdtt_evidence(unit_param_parser_deprecation unit exact_solution negative conformance)
lsdtt_evidence(unit_output_file unit characterisation)  # removed pattern vs helper
lsdtt_evidence(unit_hydraulic_routing unit property metamorphic)
foreach(_t IN ITEMS wheel_cli_inventory version_consistency build_is_current
                    test_inventory_audit docs_test_commands
                    test_inventory_audit_selftest test_report_selftest)
    lsdtt_evidence(unit_${_t} unit conformance)
endforeach()
# The pytest markers of python_bindings/tests, read with ast; the self-test
# runs the checker on snippets whose verdicts are written by hand.
lsdtt_evidence(unit_python_test_markers unit conformance)
lsdtt_evidence(unit_python_test_markers_selftest unit conformance exact_solution)
# The five binding-free pytest files: their markers are conformance (7),
# exact_solution (3) and negative (1), all unit.
lsdtt_evidence(unit_python_binding_free_tests unit conformance exact_solution negative)
# The generated test inventory must match the registered tests (a docs check).
lsdtt_evidence(unit_test_inventory_is_current unit conformance)

# -----------------------------------------------------------------------------
# Benchmarks (LSDTT_ENABLE_BENCHMARKS) and performance (LSDTT_ENABLE_PERF_TESTS)
# -----------------------------------------------------------------------------
foreach(_fx IN ITEMS alpine_small lakeland_small)
    foreach(_case IN ITEMS fill fill_epsilon breach breach_least_cost accum_d8 flats_float_step)
        # --check against RichDEM / TopoToolbox / Whitebox; timings reported.
        lsdtt_evidence(benchmark_${_case}_${_fx} system differential performance)
    endforeach()
endforeach()
lsdtt_evidence(benchmark_breach_least_cost_fuzz system differential)  # 200 DEMs bitwise vs Whitebox
foreach(_fx IN ITEMS alpine_medium lakeland_medium)
    foreach(_case IN ITEMS slope curvature_all4 hillshade)
        # A timing benchmark; --check is only a rank-correlation guard
        # against TopoToolbox.
        lsdtt_evidence(benchmark_surface_${_case}_${_fx} system performance differential)
    endforeach()
endforeach()
lsdtt_evidence(benchmark_chi_value_topotoolbox system differential)
lsdtt_evidence(benchmark_flat_resolution_richdem system differential property)
foreach(_fx IN ITEMS tiny_dem.asc alpine_small)
    foreach(_sc IN ITEMS stationary_fixed transient_dynamic)
        lsdtt_evidence(benchmark_hydraulic_dagger_${_sc}_${_fx} system differential)
    endforeach()
endforeach()
lsdtt_evidence(performance_flow_routines unit performance)
lsdtt_evidence(performance_hillslope_implicit unit performance)
lsdtt_evidence(perf_cosmo_lsdn_production unit performance)
lsdtt_evidence(performance_surface_metrics system performance)
lsdtt_evidence(performance_hype_model unit performance property)  # fails on water balance only

# =============================================================================
# FINAL PASS AND GUARD
# =============================================================================
get_property(_lsdtt_ev_tests DIRECTORY PROPERTY TESTS)
get_property(_lsdtt_ev_py GLOBAL PROPERTY LSDTT_RUNTIME_python_TESTS)
get_property(_lsdtt_ev_off GLOBAL PROPERTY LSDTT_UNREGISTERED_TESTS)

# A table row for a test that is not registered is a typo, unless an option
# switched the test off: recorded by lsdtt_mark_test_unregistered(), a
# benchmark or performance test (as in the DEPENDENCY TIERS lists), or a
# downloading test while LSDTT_ENABLE_NETWORK_TESTS is OFF.
foreach(_row IN LISTS _lsdtt_ev_rows)
    if(_row STREQUAL "")
        continue()
    endif()
    set(_switched_off FALSE)
    if("${_row}" IN_LIST _lsdtt_ev_off OR _row MATCHES "^(benchmark_|performance_|perf_)")
        set(_switched_off TRUE)
    elseif(NOT LSDTT_ENABLE_NETWORK_TESTS AND _row MATCHES "${LSDTT_NEEDS_NETWORK_PATTERN}")
        set(_switched_off TRUE)
    endif()
    if(NOT TEST "${_row}" AND NOT _switched_off)
        list(APPEND _lsdtt_ev_errors "${_row}: in the evidence table but not a registered test")
    endif()
endforeach()

set(_lsdtt_ev_handwritten "")
set(_lsdtt_ev_no_evidence "")
set(_lsdtt_ev_no_level "")
set(_lsdtt_ev_n_multi 0)
foreach(_lab IN LISTS LSDTT_EVIDENCE_LABELS LSDTT_LEVEL_LABELS)
    set(_lsdtt_ev_count_${_lab} 0)
endforeach()

foreach(_test IN LISTS _lsdtt_ev_tests)
    get_test_property(${_test} LABELS _existing)
    if(NOT _existing)
        set(_existing "")
    endif()
    foreach(_lab IN LISTS _existing)
        if("${_lab}" IN_LIST LSDTT_EVIDENCE_LABELS OR "${_lab}" IN_LIST LSDTT_LEVEL_LABELS
           OR "${_lab}" IN_LIST LSDTT_RETIRED_KIND_LABELS)
            list(APPEND _lsdtt_ev_handwritten "${_test} (${_lab})")
        endif()
    endforeach()

    # Evidence: registration, re-run list, name rule, table.
    get_property(_ev GLOBAL PROPERTY LSDTT_EVIDENCE_${_test})
    get_property(_level GLOBAL PROPERTY LSDTT_LEVEL_${_test})
    if("${_test}" IN_LIST LSDTT_RERUN_FROM_INGESTED_TESTS)
        list(APPEND _ev metamorphic)
    endif()
    set(_rule_level "")
    list(LENGTH LSDTT_EVIDENCE_NAME_RULES _n)
    math(EXPR _last "${_n} - 1")
    foreach(_i RANGE 0 ${_last} 3)
        math(EXPR _j "${_i} + 1")
        math(EXPR _k "${_i} + 2")
        list(GET LSDTT_EVIDENCE_NAME_RULES ${_i} _pattern)
        if(_test MATCHES "${_pattern}")
            list(GET LSDTT_EVIDENCE_NAME_RULES ${_j} _rule_ev)
            list(GET LSDTT_EVIDENCE_NAME_RULES ${_k} _rule_level)
            list(APPEND _ev ${_rule_ev})
            break()
        endif()
    endforeach()
    if(DEFINED _lsdtt_ev_table_${_test})
        list(APPEND _ev ${_lsdtt_ev_table_${_test}})
        set(_level "${_lsdtt_ev_level_${_test}}")
    endif()
    if("${_level}" STREQUAL "")
        if(NOT "${_rule_level}" STREQUAL "")
            set(_level "${_rule_level}")
        elseif(_test MATCHES "^unit_" AND NOT "${_test}" IN_LIST _lsdtt_ev_py)
            set(_level unit)
        endif()
    endif()

    if(_ev)
        list(REMOVE_DUPLICATES _ev)
    endif()
    foreach(_lab IN LISTS _ev)
        if(NOT "${_lab}" IN_LIST LSDTT_EVIDENCE_LABELS)
            list(APPEND _lsdtt_ev_errors "${_test}: '${_lab}' is not an evidence label")
        endif()
    endforeach()
    list(LENGTH _ev _n_ev)
    if(_n_ev EQUAL 0)
        list(APPEND _lsdtt_ev_no_evidence "${_test}")
    elseif(_n_ev GREATER 1)
        math(EXPR _lsdtt_ev_n_multi "${_lsdtt_ev_n_multi} + 1")
    endif()
    list(LENGTH _level _n_level)
    if(NOT _n_level EQUAL 1 OR NOT "${_level}" IN_LIST LSDTT_LEVEL_LABELS)
        string(REPLACE ";" "," _given "${_level}")
        list(APPEND _lsdtt_ev_no_level "${_test} ('${_given}')")
        continue()
    endif()

    set_property(TEST ${_test} APPEND PROPERTY LABELS ${_ev} ${_level})
    foreach(_lab IN LISTS _ev _level)
        if(DEFINED _lsdtt_ev_count_${_lab})
            math(EXPR _lsdtt_ev_count_${_lab} "${_lsdtt_ev_count_${_lab}} + 1")
        endif()
    endforeach()
endforeach()

set(_lsdtt_ev_report "")
if(_lsdtt_ev_handwritten)
    string(REPLACE ";" "\n    " _p "${_lsdtt_ev_handwritten}")
    string(APPEND _lsdtt_ev_report
        "These tests set an evidence, level or retired kind label at "
        "registration. Remove it there: tests/cmake/evidence_labels.cmake is "
        "the only writer of those labels.\n    ${_p}\n")
endif()
if(_lsdtt_ev_no_evidence)
    string(REPLACE ";" "\n    " _p "${_lsdtt_ev_no_evidence}")
    string(APPEND _lsdtt_ev_report
        "These tests have no evidence label. Add an "
        "lsdtt_evidence() row saying how each knows its answer:\n    ${_p}\n")
endif()
if(_lsdtt_ev_no_level)
    string(REPLACE ";" "\n    " _p "${_lsdtt_ev_no_level}")
    string(APPEND _lsdtt_ev_report
        "These tests do not have exactly one level label (unit or system). "
        "Give each an lsdtt_evidence() row with its level:\n    ${_p}\n")
endif()
list(REMOVE_ITEM _lsdtt_ev_errors "")
if(_lsdtt_ev_errors)
    list(REMOVE_DUPLICATES _lsdtt_ev_errors)
endif()
if(_lsdtt_ev_errors)
    string(REPLACE ";" "\n    " _p "${_lsdtt_ev_errors}")
    string(APPEND _lsdtt_ev_report "Evidence table errors:\n    ${_p}\n")
endif()
if(NOT _lsdtt_ev_report STREQUAL "")
    message(FATAL_ERROR "${_lsdtt_ev_report}"
        "The evidence and level labels are defined in "
        "tests/cmake/evidence_labels.cmake.")
endif()

set(_lsdtt_ev_summary "")
foreach(_lab IN LISTS LSDTT_EVIDENCE_LABELS)
    string(APPEND _lsdtt_ev_summary "${_lab} ${_lsdtt_ev_count_${_lab}}, ")
endforeach()
message(STATUS "Test evidence : ${_lsdtt_ev_summary}"
               "more than one ${_lsdtt_ev_n_multi}")
message(STATUS "Test level    : unit ${_lsdtt_ev_count_unit}, system ${_lsdtt_ev_count_system}")
