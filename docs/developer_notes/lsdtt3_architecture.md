# lsdtt3 Architecture Document

A guide to the architecture of lsdtt3. 
Use for context when interacting with AI  during code development. 
Can also be used as a developer guide.

---

## Project identity

- **Package name:** lsdtt3
- **Repository:** https://github.com/LSDtopotools/lsdtt3
- **Copyright:** Simon M. Mudd
- **Licence:** AGPL-3.0-or-later
- **C++ standard:** C++17
- **Build system:** CMake 3.15+

### Code text policy

- In source code, build scripts, and machine-read config files, use ASCII characters only.
- Replace typographic punctuation and Unicode symbols with ASCII equivalents (for example: `-`, `->`, `<=`, `>=`, `x`).
- Keep Unicode only in user-facing documentation where needed.

---

## Directory structure

```
lsdtt3/
├── apptainer_def/
│   └── apptainer_setup.md
├── .github/workflows/                     — release.yml (Linux bundle on v* tags), python-wheels.yml (manual), docs.yml
├── cli/                                   — command-line executables (lsdtt-* dashed names)
│   ├── lsdtt-advection-diffusion.cpp
│   ├── lsdtt-bankfull-width.cpp
│   ├── lsdtt-basin-extraction.cpp
│   ├── lsdtt-basin-stats.cpp
│   ├── lsdtt-basin-surface-metrics.cpp
│   ├── lsdtt-channel-extraction.cpp
│   ├── lsdtt-channel-network-analysis.cpp
│   ├── lsdtt-chi-analysis.cpp             *(requires Eigen3)*
│   ├── lsdtt-cosmo-column-tool.cpp
│   ├── lsdtt-cosmo-landslides.cpp
│   ├── lsdtt-cosmo-prediction.cpp
│   ├── lsdtt-cosmo-tool.cpp
│   ├── lsdtt-cross-sections.cpp
│   ├── lsdtt-distributed-hydrology.cpp
│   ├── lsdtt-divide-migration.cpp
│   ├── lsdtt-drainage-analysis.cpp
│   ├── lsdtt-feature-finder.cpp
│   ├── lsdtt-fetch-erosion-rates.cpp
│   ├── lsdtt-fetch-lithology.cpp
│   ├── lsdtt-fetch-raster.cpp
│   ├── lsdtt-flat-surface-finder.cpp
│   ├── lsdtt-fluvial-surfaces.cpp
│   ├── lsdtt-graphflood.cpp
│   ├── lsdtt-hillslope-channel-coupling.cpp
│   ├── lsdtt-internal-drainage.cpp
│   ├── lsdtt-lem.cpp
│   ├── lsdtt-linear-stats.cpp
│   ├── lsdtt-lithology-erodibility.cpp
│   ├── lsdtt-marsh-topography.cpp
│   ├── lsdtt-mountain-front.cpp
│   ├── lsdtt-outlet-matching.cpp
│   ├── lsdtt-raster-calculator.cpp
│   ├── lsdtt-raster-preprocessing.cpp
│   ├── lsdtt-ridgeline-analysis.cpp
│   ├── lsdtt-sediment-routing.cpp
│   ├── lsdtt-slope-area.cpp
│   ├── lsdtt-spectral-analysis.cpp        *(requires FFTW3)*
│   ├── lsdtt-stream-burning.cpp
│   ├── lsdtt-surface-metrics.cpp
│   ├── lsdtt-swath-profile.cpp
│   ├── lsdtt-synthetic-networks.cpp
│   ├── lsdtt-synthetic-surfaces.cpp
│   ├── lsdtt-valley-centreline.cpp
│   ├── lsdtt-water-centreline.cpp
│   └── lsdtt-water-line-repair.cpp
├── conda-recipe/                           — conda-build recipe (local builds; conda-forge once public)
├── data/                                   — runtime data shipped with the tools (crn_parameters.nc)
├── debian/                                 — Debian package definition (see debian_packaging.md)
├── docs/                                   — MkDocs source, organised by purpose
│   ├── index.md                            — documentation site root
│   ├── what_is_lsdtt3/                     — overview, install, data, routing, workflows
│   ├── design_and_algorithms/              — scientific design and algorithms
│   ├── cli_reference/                      — one parameter reference per CLI
│   ├── developer_notes/                    — how to work on the code: architecture, building, packaging, testing
│   │   ├── lsdtt3_architecture.md          — this file
│   │   ├── developer_onboarding.md         — start here as a new contributor
│   │   ├── building_from_source.md         — build hub + reference Linux build
│   │   ├── building_on_windows.md          — WSL2 and native (QGIS backend)
│   │   ├── building_on_macos.md            — native and Multipass
│   │   ├── building_the_docs.md
│   │   ├── building_python_wheels.md       — wheel builds, local and on GitHub (Linux-only or all platforms)
│   │   ├── python_bindings_architecture.md — bindings layers, API contract files, adding a binding
│   │   ├── publishing_to_pypi.md
│   │   ├── making_release_packages.md      — one version number, what is built where
│   │   ├── qgis_plugin_development.md
│   │   ├── debian_packaging.md
│   │   ├── beta_testing_notes.md
│   │   ├── git_notes.md
│   │   ├── backend_abstraction_rules.md, basin_object_architecture.md, gpu_compute.md
│   │   ├── lsdtt_utils_guide.md, using_lsdtt_flowinfo.md, cli_parameter_conventions.md
│   │   └── testing_and_benchmarking.md
│   ├── python/                             — Python package guides
│   ├── python_api/                         — generated Python API reference (do not edit)
│   ├── qgis.md                             — the QGIS plugin, for users
│   └── tex_technical/                      — LaTeX method write-ups (untracked)
├── mkdocs.yml                              — documentation site configuration
├── pixi.toml                               — Python environments for tests, docs, packaging and benchmarks
├── python_bindings/                        — the lsdtt3 Python package (nanobind)
│   ├── bindings/                           — C++ adapters building the _core extension
│   ├── src/lsdtt3/                         — public package, wrappers, console-script dispatch
│   ├── api/                                — public_api.json, decisions.json: reviewed API contracts
│   ├── tools/, tests/
│   └── pyproject.toml, CMakeLists.txt      — a standalone CMake project compiling ../src
├── python_helpers/                         — stand-alone download, data-preparation and plotting scripts
├── qgis_plugin/                            — QGIS Processing plugin (lsdtt_qgis/) and package_plugin.py
├── src/
│   ├── core/
│   │   ├── lsdtt_types.hpp                  — all lightweight shared primitive types (see below)
│   │   ├── lsdtt_raster_meta.hpp            — metadata struct + CRS validation
│   │   ├── lsdtt_raster_base.hpp/.cpp       — abstract raster dispatch interface
│   │   ├── lsdtt_raster_result.hpp          — plain data-transfer struct for rasters returned by core-layer engines
│   │   ├── lsdtt_gdal_init.hpp/.cpp         — one-time GDAL/PROJ init and data pinning
│   │   ├── lsdtt_pixel_type.hpp/.cpp        — GDAL-free sample-type enum
│   │   ├── lsdtt_pixel_type_gdal.hpp        — the GDAL bridge (.cpp only)
│   │   ├── lsdtt_data_loader.hpp/.cpp       — GDAL/OGR I/O, RAII, COG output
│   │   ├── lsdtt_backend_factory.hpp/.cpp   — make_raster_engine(), gpu_available()
│   │   ├── lsdtt_param_parser.hpp/.cpp      — parameter file parser
│   │   ├── lsdtt_common_parameters.hpp/.cpp — shared parameter declarations reused across CLIs
│   │   ├── lsdtt_csv_table.hpp/.cpp         — minimal header-keyed CSV reader for CLI input tables
│   │   ├── lsdtt_secret.hpp/.cpp            — credential resolution and masking (API keys never in parameter files)
│   │   ├── lsdtt_exe_path.hpp/.cpp          — locating the running executable and its bundled data files
│   │   ├── lsdtt_basin_base.hpp             — abstract basin interface
│   │   ├── lsdtt_junction_network_base.hpp  — abstract channel network interface + lsdtt_junction_network_options
│   │   ├── lsdtt_flow_base.hpp              — abstract D8 flow interface + lsdtt_flowinfo_options
│   │   ├── lsdtt_hydraulic_base.hpp         — abstract graph-based hydraulic routing interface + lsdtt_hydraulic_options
│   │   ├── lsdtt_raster_meta.hpp/.cpp       — raster metadata; CRS methods in the .cpp
│   │   ├── lsdtt_crs_suggest.hpp/.cpp       — CRS inspection and projection recommendation
│   │   ├── lsdtt_raster_reproject.hpp/.cpp  — CRS reprojection
│   │   ├── lsdtt_raster_warp.hpp/.cpp       — GDAL warp helpers
│   │   ├── lsdtt_point_input.hpp/.cpp       — point/coordinate input handling
│   │   ├── lsdtt_outlet_bounds_discovery.hpp/.cpp— WGS84 outlet points to fetch bounds
│   │   ├── lsdtt_outlet_clusters.hpp/.cpp   — grouping outlets into fetch clusters
│   │   ├── lsdtt_dem_grouping.hpp/.cpp      — grouping bounding boxes for higher-resolution re-fetching
│   │   ├── lsdtt_coastal_mask.hpp/.cpp      — stops bulk tile downloads over the ocean (with dem_grouping)
│   │   ├── lsdtt_raster_fetch.hpp/.cpp      — native DEM fetching: WGS84 box to analysis-ready raster
│   │   ├── lsdtt_http_fetch.hpp/.cpp        — HTTP GET/POST through GDAL's curl (timeouts, retries)
│   │   ├── lsdtt_macrostrat_fetch.hpp/.cpp  — Macrostrat map units for a box (vector tiles + API)
│   │   ├── lsdtt_octopus_fetch.hpp/.cpp     — OCTOPUS basin erosion rates for a box (keyless WFS)
│   │   ├── lsdtt_stat_csv.hpp/.cpp          — CSV writers for QQ, CDF, and distribution-fit tables
│   │   ├── lsdtt_distributed_hydrology_forcing.hpp/.cpp— per-timestep CSV forcing readers for lsdtt-distributed-hydrology
│   │   └── lsdtt_vector_utils.hpp/.cpp      — GDAL/OGR-backed vector helpers
│   ├── cpu/                                   — 154 CPU-specific analysis modules (flow routing, surface metrics, channel analysis, LEM, cosmogenic nuclides, etc.)
│   │   ├── lsdtt_raster_cpu.hpp/.cpp                         — concrete CPU raster class
│   │   ├── lsdtt_flow_base.cpp                               — factory impl for make_flow_base()
│   │   ├── lsdtt_flowinfo_cpu.hpp/.cpp                       — D8 CPU flow routing
│   │   ├── lsdtt_flowgraph_cpu.hpp/.cpp                      — D-inf/MFD DAG routing (CSR layout, Kahn BFS)
│   │   ├── lsdtt_depression_hierarchy_cpu.hpp/.cpp           — depression hierarchy construction
│   │   ├── lsdtt_fill_spill_merge_cpu.hpp/.cpp               — fill-spill-merge depression routing
│   │   ├── lsdtt_flow_snapping_cpu.hpp/.cpp                  — snapping points to channel network
│   │   ├── lsdtt_junction_network_cpu.hpp/.cpp               — CPU channel junction network
│   │   ├── lsdtt_junction_network_base.cpp                   — factory impl for junction_network_base
│   │   ├── lsdtt_basin_cpu.hpp/.cpp                          — CPU basin objects
│   │   ├── lsdtt_junction_basin_cpu.hpp/.cpp                 — junction-based basin extraction
│   │   ├── lsdtt_junction_basin_manager_cpu.hpp/.cpp         — basin selection and management
│   │   ├── lsdtt_basin_outline_cpu.hpp/.cpp                  — basin outline/polygon extraction
│   │   ├── lsdtt_basin_outline_metrics_cpu.hpp/.cpp              — outline curvature flips and slope-area scale metrics on basin outlines
│   │   ├── lsdtt_basin_stats_cpu.hpp/.cpp                        — per-basin statistics for a set of outlet points
│   │   ├── lsdtt_zero_order_basin_manager_cpu.hpp/.cpp           — basins whose outlets are channel-head source pixels
│   │   ├── lsdtt_outlet_resolution_cpu.hpp/.cpp                  — turning a request for basin outlets into flat indices
│   │   ├── lsdtt_outlet_basin_verification_cpu.hpp/.cpp          — checking outlet basins fit inside a fetched DEM
│   │   ├── lsdtt_raster_algebra.hpp/.cpp                     — element-wise raster algebra
│   │   ├── lsdtt_raster_surface_metrics_cpu.hpp/.cpp         — hillshade, moving window, polyfit metrics, surface-shape (SMAP) classification
│   │   ├── lsdtt_surface_texture_cpu.hpp/.cpp                — MADk2, moving_window_iqr, NVD roughness metrics
│   │   ├── lsdtt_raster_surface_fitting_cpu.hpp/.cpp         — plane and polynomial surface fitting to masks
│   │   ├── lsdtt_raster_flow_preprocessing_cpu.hpp/.cpp      — fill (Priority-Flood, reconstruction), Lindsay 2016 breach, nodata
│   │   ├── lsdtt_raster_breach_least_cost_cpu.hpp/.cpp       — least-cost breaching (Lindsay 2020; bitwise port of Whitebox 2.0.6)
│   │   ├── lsdtt_raster_flat_resolution_cpu.hpp/.cpp         — Barnes 2014b flat resolution (mask, rounds, final pass)
│   │   ├── lsdtt_raster_nodata_filling_cpu.hpp/.cpp          — interior void filling
│   │   ├── lsdtt_fill_flat_resolution_reporting_cpu.hpp/.cpp     — fill_flat_resolution_choice checks and reporting shared by the CLIs
│   │   ├── lsdtt_raster_resampling_cpu.hpp/.cpp                  — GDAL-free resampling of in-memory rasters
│   │   ├── lsdtt_raster_hierarchy_cpu.hpp/.cpp                   — in-memory resolution hierarchy (restrict / prolong)
│   │   ├── lsdtt_raster_sampling_cpu.hpp/.cpp                    — sampling a raster at a map coordinate
│   │   ├── lsdtt_raster_cell_predicates_cpu.hpp/.cpp             — per-cell predicates and index helpers for flat raster data
│   │   ├── lsdtt_raster_central_difference_cpu.hpp/.cpp          — four-cardinal centred-difference gradient at one cell
│   │   ├── lsdtt_usgs_elevation.hpp/.cpp                         — settling a DEM's vertical unit against USGS 3DEP
│   │   ├── lsdtt_raster_spectral_cpu.hpp/.cpp                — wavelet quadratic fit; spectral filters; fftw2D_masked_surface_recovery
│   │   ├── lsdtt_raster_image_ops_cpu.hpp/.cpp               — Otsu, adaptive threshold, Canny, probabilistic Hough (native; replaced OpenCV)
│   │   ├── lsdtt_raster_feature_finder_cpu.hpp/.cpp          — feature-specific raster workflows
│   │   ├── lsdtt_raster_linear_features_cpu.hpp/.cpp         — curvilinear primitives; path opening, feature_type classification
│   │   ├── lsdtt_raster_binary_mask_cpu.hpp/.cpp             — binary mask morphology (dilation, erosion, opening, closing, labelling)
│   │   ├── lsdtt_raster_binary_mask_vectorise_cpu.hpp/.cpp   — binary mask to vector polygon conversion
│   │   ├── lsdtt_raster_skeleton_graph_cpu.hpp/.cpp          — 1-pixel skeleton graph analysis, pruning, reconnection
│   │   ├── lsdtt_raster_vectorise_cpu.hpp/.cpp               — raster to vector conversion
│   │   ├── lsdtt_raster_mosaic_cpu.hpp/.cpp                  — raster mosaicking
│   │   ├── lsdtt_raster_contours_cpu.hpp/.cpp                — contour extraction and vectorisation
│   │   ├── lsdtt_raster_coregistration_cpu.hpp/.cpp          — raster coregistration and alignment
│   │   ├── lsdtt_vector_rasterise_cpu.hpp/.cpp               — vector to raster rasterisation
│   │   ├── lsdtt_vector_ops_cpu.hpp/.cpp                     — vector geometry operations and repairs
│   │   ├── lsdtt_vector_network_d8_clean_cpu.hpp/.cpp        — D8 network vector cleaning and repair
│   │   ├── lsdtt_vector_line_repair_cpu.hpp/.cpp             — polyline repair and topology fixing
│   │   ├── lsdtt_vector_to_junction_network.hpp/.cpp         — convert vector networks to junction networks
│   │   ├── lsdtt_burn_vectormap_cpu.hpp/.cpp                 — burn rasterized vector into DEM
│   │   ├── lsdtt_connectivity_cpu.hpp/.cpp                   — raster connectivity analysis
│   │   ├── lsdtt_proximity_cpu.hpp/.cpp                      — proximity / distance-to-feature rasters
│   │   ├── lsdtt_swath_cpu.hpp/.cpp                          — swath profile extraction
│   │   ├── lsdtt_cross_sections_cpu.hpp/.cpp                 — cross-section profile extraction
│   │   ├── lsdtt_mountain_front_cpu.hpp/.cpp                 — mountain front identification from DEM
│   │   ├── lsdtt_streamlines_cpu.hpp/.cpp                    — kinematic streamline tracing (Lea 1992)
│   │   ├── lsdtt_lem_cpu.hpp/.cpp                            — LEM step orchestrator (fluvial + hillslope + uplift)
│   │   ├── lsdtt_lem_timestep_cpu.hpp/.cpp                       — timestep control for the coupled LEM
│   │   ├── lsdtt_lem_parameter_rasters_cpu.hpp/.cpp          — per-pixel LEM parameters: K calibration, conversion, validation, forcing scaling
│   │   ├── lsdtt_lem_input_grid_cpu.hpp/.cpp                 — LEM model grid and input-raster placement (clip / relabel / warp)
│   │   ├── lsdtt_lithology_map_cpu.hpp/.cpp                  — lithology maps to erodibility: property table, crosswalk, K raster
│   │   ├── lsdtt_lithology_channel_stats_cpu.hpp/.cpp            — channel statistics (k_sn) by lithologic unit
│   │   ├── lsdtt_lithology_relative_erodibility_cpu.hpp/.cpp     — relative erodibility of lithologic units from channel steepness
│   │   ├── lsdtt_erodibility_calibration_cpu.hpp/.cpp            — stream-power K (and n) from basin-averaged erosion rates
│   │   ├── lsdtt_stream_power_cpu.hpp/.cpp                   — implicit stream-power incision (D8 + D-inf)
│   │   ├── lsdtt_steady_state_cpu.hpp/.cpp                   — analytical chi-based fluvial steady state (shared core compute_fluvial_steady_state_on_tree, reused by the OCN generator)
│   │   ├── lsdtt_hillslope_transport_cpu.hpp/.cpp             — transient hillslope transport (cidre, roering, linear, adi)
│   │   ├── lsdtt_hillslope_implicit_cpu.hpp/.cpp              — implicit transient Roering solver (Eigen; backward-Euler/TR-BDF2 + damped Picard, channels as Dirichlet sinks)
│   │   ├── lsdtt_hillslope_steady_state_cpu.hpp/.cpp          — steady-state hillslope solver (SOR + pseudo-transient Cidre)
│   │   ├── lsdtt_hillslope_newton_cpu.hpp/.cpp                — steady/transient Newton Roering solver
│   │   ├── lsdtt_hillslope_cellular_cpu.hpp/.cpp              — soft-critical face transport, nonlocal routing, sparse particles
│   │   ├── lsdtt_hillslope_operator_cpu.hpp/.cpp                 — THE shared discrete hillslope operator every solver uses
│   │   ├── lsdtt_hillslope_length_search_cpu.hpp/.cpp            — channel threshold for a target hillslope length
│   │   ├── lsdtt_hillslope_peclet_cpu.hpp/.cpp                   — K and D from a target relief and hillslope length
│   │   ├── lsdtt_perron_qimp_cpu.hpp/.cpp                        — independent reference implementation of Perron (2011) Q-imp
│   │   ├── lsdtt_ren_modular_implicit_cpu.hpp/.cpp               — independent port of the Ren et al. (2026) modular implicit method
│   │   ├── lsdtt_hillslope_particle_io.hpp/.cpp               — particle CSV restart/output
│   │   ├── lsdtt_landscape_drape_cpu.hpp/.cpp                — Verlet cloth drape: upper envelope surface and incision depth
│   │   ├── lsdtt_landscape_equilibration_cpu.hpp/.cpp        — divide equilibration and landscape relaxation
│   │   ├── lsdtt_flat_surface_finder.hpp/.cpp                — flat surface identification (MADk2, residual IQR, NVD)
│   │   ├── lsdtt_flat_surface_plane_analysis.hpp/.cpp        — plane analysis on detected flat surfaces
│   │   ├── lsdtt_flat_surface_channel_classify_cpu.hpp/.cpp  — classify flat surfaces against a channel; shared surface/patch/ladder CSV writers
│   │   ├── lsdtt_reference_field_ladder_cpu.hpp/.cpp        — reference-field ladder: band a height field, label components WITHIN each band
│   │   ├── lsdtt_hydraulic_graph_cpu.hpp/.cpp                — CPU hydraulic routing (lsdtt_hydraulic_base impl, Gailleton et al. 2024)
│   │   ├── lsdtt_graphflood_optimised_cpu.hpp/.cpp               — native coarse-to-fine coordinator for the v2 GraphFlood solver
│   │   ├── lsdtt_manning_formula_cpu.hpp/.cpp                    — Manning's equation for mean flow velocity
│   │   ├── lsdtt_manning_hydraulics_cpu.hpp/.cpp                 — Manning's-equation bankfull width estimation
│   │   ├── lsdtt_hydro_roughness_cpu.hpp/.cpp                    — the roughness (n) field
│   │   ├── lsdtt_hydro_routing_cpu.hpp/.cpp                      — routing ladder: translation, storage, kinematic wave
│   │   ├── lsdtt_hydro_velocity_field_cpu.hpp/.cpp               — distributed velocity, travel time, width-function IUH
│   │   ├── lsdtt_wetness_index_cpu.hpp/.cpp                      — topographic wetness index
│   │   ├── lsdtt_internal_drainage_cpu.hpp/.cpp                  — closed-depression detection and the review round-trip
│   │   ├── lsdtt_subbasin_partition_cpu.hpp/.cpp                 — incremental sub-basin partition and per-sub-basin set-up
│   │   ├── lsdtt_topmodel_cpu.hpp/.cpp                           — TOPMODEL column physics, index mode
│   │   ├── lsdtt_topmodel_dynamic_cpu.hpp/.cpp                   — TOPMODEL with a prognostic saturated zone
│   │   ├── lsdtt_canopy_interception_cpu.hpp/.cpp                — canopy interception store
│   │   ├── lsdtt_infiltration_excess_cpu.hpp/.cpp                — Green-Ampt infiltration excess
│   │   ├── lsdtt_landcover_scenario_cpu.hpp/.cpp                 — land-cover scenarios: class raster to derived fields
│   │   ├── lsdtt_hype_parameters.hpp/.cpp                        — HYPE parameter tables and sub-basin set-up struct
│   │   ├── lsdtt_hype_column_cpu.hpp/.cpp                        — HYPE land column for one (sub-basin, class) unit
│   │   ├── lsdtt_hype_routing_cpu.hpp/.cpp                       — HYPE river routing
│   │   ├── lsdtt_hype_raster_setup_cpu.hpp/.cpp                  — HYPE set-up from rasters (sub-basins, GeoData, class fractions)
│   │   ├── lsdtt_hype_model_cpu.hpp/.cpp                         — the HYPE model of lsdtt-distributed-hydrology
│   │   ├── lsdtt_distributed_hydrology_argument_checks_cpu.hpp/.cpp— argument checks shared by the distributed-hydrology kernels
│   │   ├── lsdtt_chi_extraction_cpu.hpp/.cpp                 — chi coordinate extraction (full raster, outlet-seeded, inherited-base)
│   │   ├── lsdtt_chi_network_analysis_cpu.hpp/.cpp           — piecewise-linear fitting of chi channel networks
│   │   ├── lsdtt_chi_disorder_cpu.hpp/.cpp                   — chi disorder metric, jackknife, bootstrap, m/n sweep
│   │   ├── lsdtt_chi_data_maps_cpu.hpp/.cpp                  — chi data-map output (per-node points, basin outlines, polylines)
│   │   ├── lsdtt_chi_ksn_inversion_cpu.hpp/.cpp              — spatial k_sn/u* inversion with m/n model selection *(requires Eigen3)*
│   │   ├── lsdtt_channel_network_analysis_cpu.hpp/.cpp       — channel link bearings, junction angles
│   │   ├── lsdtt_channel_network_grouping_cpu.hpp/.cpp           — grouping keys for link and junction records (Strahler, basin, region, lithology)
│   │   ├── lsdtt_channel_bearing_statistics_cpu.hpp/.cpp         — grouped bearing and junction-angle distributions and comparisons
│   │   ├── lsdtt_channel_extent_cpu.hpp/.cpp                     — channel extent: where the channels are
│   │   ├── lsdtt_channel_sources_wiener_cpu.hpp/.cpp         — channel source detection via Wiener filtering
│   │   ├── lsdtt_channel_scope_cpu.hpp/.cpp                  — channel scope resolution; shared CLI plumbing (point loading, scope diagnostics, network-layer writer)
│   │   ├── lsdtt_strahler_tokunaga_analysis_cpu.hpp/.cpp     — Strahler channel aggregation, Tokunaga self-similarity, junction area-ratio distribution
│   │   ├── lsdtt_burn_streams_cpu.hpp/.cpp                   — burn polyline vector into DEM with monotone profile
│   │   ├── lsdtt_ridgeline_cpu.hpp/.cpp                      — ridgeline extraction and analysis
│   │   ├── lsdtt_valley_centreline_cpu.hpp/.cpp              — valley centreline extraction and tracing
│   │   ├── lsdtt_valley_cross_sections_cpu.hpp/.cpp          — valley cross-section analysis
│   │   ├── lsdtt_valley_spacing_cpu.hpp/.cpp                     — characteristic valley spacing (spectral and two other measures)
│   │   ├── lsdtt_water_area_centreline_cpu.hpp/.cpp          — water/floodplain centreline extraction
│   │   ├── lsdtt_marsh_platform_cpu.hpp/.cpp                 — marsh and tidal platform analysis
│   │   ├── lsdtt_tidal_channel_cpu.hpp/.cpp                  — tidal channel analysis and metrics
│   │   ├── lsdtt_tidal_skeleton_metrics_cpu.hpp/.cpp         — tidal network skeleton analysis
│   │   ├── lsdtt_terrace_edge_tracer_cpu.hpp/.cpp            — terrace edge and scarp tracing
│   │   ├── lsdtt_sediment_transport_cpu.hpp/.cpp             — sediment transport and routing
│   │   ├── lsdtt_giuh_cpu.hpp/.cpp                           — geomorphologic unit hydrograph
│   │   ├── lsdtt_side_tributary_basins_cpu.hpp/.cpp          — side tributary basin extraction
│   │   ├── lsdtt_slope_break_finder_cpu.hpp/.cpp             — slope break and terrain transition detection
│   │   ├── lsdtt_labelled_raster_stats_cpu.hpp/.cpp          — statistics on labelled raster regions
│   │   ├── lsdtt_least_cost_path.hpp/.cpp                    — least-cost path computation
│   │   ├── lsdtt_cosmo_raster_cpu.hpp/.cpp                   — raster CRN: production, SS/step-change concentration, flow accumulation
│   │   ├── lsdtt_cosmo_site_cpu.hpp/.cpp                     — single-site CRN calculations
│   │   ├── lsdtt_crn_parameters.hpp/.cpp                     — Stone 2000 scaling, muon schemes, ERA40/palaeomagnetic correction
│   │   ├── lsdtt_crn_column.hpp/.cpp                         — 1-D transient CRN column simulation
│   │   ├── lsdtt_crn_particle.hpp                            — header-only: SS and transient nuclide concentration particle
│   │   ├── lsdtt_crn_diagnostics.hpp/.cpp                    — CRN diagnostic utilities
│   │   ├── lsdtt_crn_units.hpp/.cpp                          — CRN unit conversion and constants
│   │   ├── lsdtt_cronus_v3_input.hpp/.cpp                    — CRONUS v3 input file preparation
│   │   ├── lsdtt_step_change_cpu.hpp/.cpp                    — step-change erosion scenario modelling
│   │   ├── lsdtt_sti_cpu.hpp/.cpp                            — sediment trace isotope calculations
│   │   ├── lsdtt_niemi_ensemble_cpu.hpp/.cpp                 — Niemi ensemble equilibrium models
│   │   ├── lsdtt_landslide_crn_cpu.hpp/.cpp                  — cosmogenic nuclide analysis for landslide samples
│   │   ├── lsdtt_landslide_event_cpu.hpp/.cpp                — landslide event detection and dating
│   │   ├── lsdtt_landslide_placement_cpu.hpp/.cpp            — landslide placement and spatial distribution
│   │   ├── lsdtt_landslide_size_cpu.hpp/.cpp                 — landslide size distribution analysis
│   │   ├── lsdtt_divide_migration_cpu.hpp/.cpp               — divide migration analysis and tracking
│   │   ├── lsdtt_divide_diagnostics_cpu.hpp/.cpp             — divide stability and diagnostic metrics
│   │   ├── lsdtt_critical_slope_cpu.hpp/.cpp                 — critical slope calculations
│   │   ├── lsdtt_paired_flowpaths_cpu.hpp/.cpp               — paired flowpath analysis
│   │   ├── lsdtt_hilltop_flow_routing_cpu.hpp/.cpp           — hilltop flow routing and spacing
│   │   ├── lsdtt_quartz_fertility_cpu.hpp/.cpp               — quartz fertility calculations for nuclides
│   │   ├── lsdtt_raster_synthetic_surfaces_cpu.hpp/.cpp      — test surface generators including mountain range geometries
│   │   ├── lsdtt_synthetic_network_cpu.hpp/.cpp              — OCN simulated-annealing network generator (masks + multi-outlet; tunable Tokunaga exponents)
│   │   └── lsdtt_synthetic_domain_cpu.hpp/.cpp               — island domain masks (EFD, red-noise, polygon import) for the OCN generator
│   ├── gpu/
│   │   └── lsdtt_surface_gpu.hpp/.cu        — selective-dispatch GPU surface-fitting kernels (LSDTT_GPU_KERNELS)
│   └── utils/
│       ├── lsdtt_array2d.hpp                — header-only 2D array
│       ├── lsdtt_raster_geometry.hpp        — header-only raster geometry helpers
│       ├── lsdtt_polygon_geometry.hpp       — header-only: polygon area/perimeter/moments, convex hull, min-area rect, Feret, basin shape metrics
│       ├── lsdtt_linear_algebra.hpp         — header-only: sym3x3/sym2x2 eigenvalues, 3×3 inversion, 2D geometry helpers
│       ├── lsdtt_thread_utils.hpp           — header-only: resolve_thread_count() (OpenMP / hardware fallback)
│       ├── lsdtt_random.hpp                 — header-only: splitmix64 generator, bounded(), shuffle (identical on every platform)
│       ├── lsdtt_image_filters.hpp/.cpp     — 8-bit/float image primitives: box mean, Gaussian, Sobel, Otsu, Canny, probabilistic Hough, line drawing
│       ├── lsdtt_binary_array2d_tracing.hpp/.cpp — binary array tracing utilities
│       ├── lsdtt_distributions.hpp/.cpp     — parametric distribution fitting: normal, exponential, Weibull, Pareto (MLE + KS)
│       ├── lsdtt_statistics.hpp/.cpp        — descriptive statistics (median, MAD, stddev, etc.)
│       ├── lsdtt_histogram.hpp/.cpp         — histogram, PDF, CDF
│       ├── lsdtt_statistical_tests.hpp/.cpp — nonparametric tests, KS, circular statistics
│       ├── lsdtt_splines.hpp/.cpp           — 1D spline interpolation
│       ├── lsdtt_regression.hpp/.cpp        — linear and polynomial regression helpers
│       ├── lsdtt_natural_breaks.hpp/.cpp    — Jenks natural-break thresholds
│       ├── lsdtt_signal1d.hpp/.cpp          — 1-D profile signal processing (moving average, LOWESS, derivative)
│       ├── lsdtt_boundary_conditions.hpp    — header-only: parsing the four-character N E S W boundary-condition string
│       ├── lsdtt_connected_components.hpp   — header-only: connected-component labelling on a row-major grid
│       ├── lsdtt_grid_geodesic.hpp          — header-only: Dijkstra on an 8-connected grid
│       ├── lsdtt_grid_neighbours.hpp        — header-only: neighbour offset tables for row-major rasters
│       ├── lsdtt_junction_graph_orders.hpp  — header-only: traversal orders on a compact junction graph
│       ├── lsdtt_numeric_utils.hpp          — header-only: inclusive numeric sweeps
│       ├── lsdtt_output_file.hpp            — header-only: open_output_file(): open for writing or throw with a caller's message
│       ├── lsdtt_string_utils.hpp           — header-only: string helpers shared by CLIs and library
│       └── lsdtt_statistics_detail.hpp      — internal helpers shared by the statistics implementation
├── tests/                                  — CTest suite: cpp/ drivers, scripts/, data/, cmake/
├── tools/                                  — release and packaging scripts, set_version.py, generate_cli_reference.py
├── workflows/                              — Snakemake layer for multi-stage pipelines
├── CMakeLists.txt                          — project(lsdtt3 VERSION ...): the one version number
├── LICENSE                                 — AGPL-3.0-or-later
└── README.md
```

---

## CMakeLists.txt summary

Three main static libraries:

- **`lsdtt3_core`** — abstract interfaces, factories, GDAL/OGR I/O, shared types
- **`lsdtt3_utils`** — generic mathematical / statistical utilities on valid flat data
- **`lsdtt3_cpu`** — concrete CPU implementations, raster operations, flow routing

`lsdtt3_cpu` depends on `lsdtt3_core` and `lsdtt3_utils`.
All CLI executables link against `lsdtt3_cpu` and therefore also gain access to
the utility-layer numerical functions.



**`CMakeLists.txt` is the authoritative list.** The inventories below group the
modules by subsystem so you can find your way around; they are not maintained
line-by-line against the build. To see exactly what is compiled, read
`LSDTT3_CORE_SOURCES`, `LSDTT3_UTILS_SOURCES` and `LSDTT3_CPU_SOURCES` in
`CMakeLists.txt`. At the time of writing that is 25 core, 10 utils and 148 CPU
translation units, plus `lsdtt_chi_ksn_inversion_cpu` when Eigen3 is found.


**`lsdtt3_core` sources** (`src/core/`) — GDAL-facing I/O, metadata,
parameter parsing, projection and vector plumbing:

`lsdtt_backend_factory`, `lsdtt_coastal_mask`, `lsdtt_common_parameters`, `lsdtt_crs_suggest`, `lsdtt_csv_table`, `lsdtt_data_loader`, `lsdtt_dem_grouping`, `lsdtt_gdal_init`, `lsdtt_http_fetch`, `lsdtt_macrostrat_fetch`, `lsdtt_octopus_fetch`, `lsdtt_outlet_bounds_discovery`, `lsdtt_outlet_clusters`, `lsdtt_param_parser`, `lsdtt_pixel_type`, `lsdtt_point_input`, `lsdtt_raster_base`, `lsdtt_raster_fetch`, `lsdtt_raster_meta`, `lsdtt_raster_reproject`, `lsdtt_raster_warp`, `lsdtt_secret`, `lsdtt_stat_csv`, `lsdtt_vector_utils`, `lsdtt_distributed_hydrology_forcing`


**`lsdtt3_utils` sources** (`src/utils/`) — backend-neutral numerics with no
raster or GDAL dependency (see [`lsdtt_utils_guide.md`](lsdtt_utils_guide.md)):

`lsdtt_binary_array2d_tracing`, `lsdtt_distributions`, `lsdtt_histogram`, `lsdtt_image_filters`, `lsdtt_natural_breaks`, `lsdtt_regression`, `lsdtt_signal1d`, `lsdtt_splines`, `lsdtt_statistical_tests`, `lsdtt_statistics`


**`lsdtt3_cpu` sources** (`src/cpu/`), by subsystem:


*Raster core and algebra*: `lsdtt_raster_cpu`, `lsdtt_raster_algebra`, `lsdtt_raster_hierarchy_cpu`, `lsdtt_raster_resampling_cpu`, `lsdtt_raster_mosaic_cpu`, `lsdtt_raster_coregistration_cpu`, `lsdtt_raster_nodata_filling_cpu`, `lsdtt_raster_image_ops_cpu`, `lsdtt_labelled_raster_stats_cpu`, `lsdtt_raster_sampling_cpu`


*Surface metrics and texture*: `lsdtt_raster_surface_metrics_cpu`, `lsdtt_raster_surface_fitting_cpu`, `lsdtt_surface_texture_cpu`, `lsdtt_raster_spectral_cpu`, `lsdtt_critical_slope_cpu`, `lsdtt_slope_break_finder_cpu`


*Flow routing and preprocessing*: `lsdtt_flow_base`, `lsdtt_flowinfo_cpu`, `lsdtt_flowgraph_cpu`, `lsdtt_raster_flow_preprocessing_cpu`, `lsdtt_raster_breach_least_cost_cpu`, `lsdtt_raster_flat_resolution_cpu`, `lsdtt_fill_flat_resolution_reporting_cpu`, `lsdtt_flow_snapping_cpu`, `lsdtt_depression_hierarchy_cpu`, `lsdtt_fill_spill_merge_cpu`, `lsdtt_internal_drainage_cpu`, `lsdtt_burn_streams_cpu`, `lsdtt_burn_vectormap_cpu`, `lsdtt_paired_flowpaths_cpu`, `lsdtt_least_cost_path`, `lsdtt_connectivity_cpu`, `lsdtt_wetness_index_cpu`


*Channel networks and chi*: `lsdtt_junction_network_base`, `lsdtt_junction_network_cpu`, `lsdtt_channel_network_analysis_cpu`, `lsdtt_channel_network_grouping_cpu`, `lsdtt_channel_bearing_statistics_cpu`, `lsdtt_channel_extent_cpu`, `lsdtt_channel_scope_cpu`, `lsdtt_channel_sources_wiener_cpu`, `lsdtt_chi_extraction_cpu`, `lsdtt_chi_network_analysis_cpu`, `lsdtt_chi_disorder_cpu`, `lsdtt_chi_data_maps_cpu`, `lsdtt_strahler_tokunaga_analysis_cpu`, `lsdtt_vector_to_junction_network`, `lsdtt_streamlines_cpu`


*Basins*: `lsdtt_basin_cpu`, `lsdtt_basin_outline_cpu`, `lsdtt_basin_outline_metrics_cpu`, `lsdtt_basin_stats_cpu`, `lsdtt_junction_basin_cpu`, `lsdtt_junction_basin_manager_cpu`, `lsdtt_zero_order_basin_manager_cpu`, `lsdtt_side_tributary_basins_cpu`, `lsdtt_outlet_basin_verification_cpu`, `lsdtt_outlet_resolution_cpu`


*Hillslopes and divides*: `lsdtt_hillslope_transport_cpu`, `lsdtt_hillslope_steady_state_cpu`, `lsdtt_hillslope_implicit_cpu`, `lsdtt_hillslope_newton_cpu`, `lsdtt_hillslope_cellular_cpu`, `lsdtt_hillslope_peclet_cpu`, `lsdtt_hillslope_length_search_cpu`, `lsdtt_hillslope_particle_io`, `lsdtt_hilltop_flow_routing_cpu`, `lsdtt_ridgeline_cpu`, `lsdtt_divide_diagnostics_cpu`, `lsdtt_divide_migration_cpu`, `lsdtt_perron_qimp_cpu`, `lsdtt_ren_modular_implicit_cpu`


*Landscape evolution*: `lsdtt_lem_cpu`, `lsdtt_lem_parameter_rasters_cpu`, `lsdtt_lem_input_grid_cpu`, `lsdtt_lithology_map_cpu`, `lsdtt_lem_timestep_cpu`, `lsdtt_stream_power_cpu`, `lsdtt_steady_state_cpu`, `lsdtt_sediment_transport_cpu`, `lsdtt_landscape_equilibration_cpu`, `lsdtt_landscape_drape_cpu`, `lsdtt_step_change_cpu`, `lsdtt_landcover_scenario_cpu`


*Hydrology and hydraulics*: `lsdtt_topmodel_cpu`, `lsdtt_topmodel_dynamic_cpu`, `lsdtt_hydro_routing_cpu`, `lsdtt_hydro_velocity_field_cpu`, `lsdtt_hydro_roughness_cpu`, `lsdtt_hydraulic_graph_cpu`, `lsdtt_infiltration_excess_cpu`, `lsdtt_canopy_interception_cpu`, `lsdtt_giuh_cpu`, `lsdtt_graphflood_optimised_cpu`, `lsdtt_sti_cpu`, `lsdtt_manning_hydraulics_cpu`, `lsdtt_subbasin_partition_cpu`, `lsdtt_hype_parameters`, `lsdtt_hype_column_cpu`, `lsdtt_hype_routing_cpu`, `lsdtt_hype_raster_setup_cpu`, `lsdtt_hype_model_cpu`


*Lithology and erodibility*: `lsdtt_lithology_channel_stats_cpu`, `lsdtt_lithology_relative_erodibility_cpu`, `lsdtt_erodibility_calibration_cpu` (the lithology-map reader, `lsdtt_lithology_map_cpu`, is listed under landscape evolution)


*Cosmogenic nuclides*: `lsdtt_cosmo_raster_cpu`, `lsdtt_cosmo_site_cpu`, `lsdtt_crn_parameters`, `lsdtt_crn_column`, `lsdtt_crn_diagnostics`, `lsdtt_cronus_v3_input`, `lsdtt_quartz_fertility_cpu`, `lsdtt_niemi_ensemble_cpu`


*Landslides*: `lsdtt_landslide_crn_cpu`, `lsdtt_landslide_event_cpu`, `lsdtt_landslide_placement_cpu`, `lsdtt_landslide_size_cpu`


*Valleys, floodplains and terraces*: `lsdtt_valley_centreline_cpu`, `lsdtt_valley_cross_sections_cpu`, `lsdtt_valley_spacing_cpu`, `lsdtt_cross_sections_cpu`, `lsdtt_terrace_edge_tracer_cpu`, `lsdtt_flat_surface_finder`, `lsdtt_flat_surface_plane_analysis`, `lsdtt_flat_surface_channel_classify_cpu`, `lsdtt_water_area_centreline_cpu`, `lsdtt_swath_cpu`, `lsdtt_proximity_cpu`


*Coastal and tidal*: `lsdtt_marsh_platform_cpu`, `lsdtt_tidal_channel_cpu`, `lsdtt_tidal_skeleton_metrics_cpu`


*Feature detection and image analysis*: `lsdtt_raster_feature_finder_cpu`, `lsdtt_raster_linear_features_cpu`, `lsdtt_raster_skeleton_graph_cpu`, `lsdtt_raster_binary_mask_cpu`, `lsdtt_raster_binary_mask_vectorise_cpu`, `lsdtt_raster_contours_cpu`, `lsdtt_mountain_front_cpu`


*Vector I/O and conversion*: `lsdtt_vector_ops_cpu`, `lsdtt_vector_rasterise_cpu`, `lsdtt_raster_vectorise_cpu`, `lsdtt_vector_line_repair_cpu`, `lsdtt_vector_network_d8_clean_cpu`


*Synthetic surfaces and reference fields*: `lsdtt_raster_synthetic_surfaces_cpu`, `lsdtt_synthetic_domain_cpu`, `lsdtt_synthetic_network_cpu`, `lsdtt_reference_field_ladder_cpu`


*Data acquisition*: `lsdtt_usgs_elevation`


`lsdtt_chi_ksn_inversion_cpu` is compiled only when Eigen3 is available.

`lsdtt3_cpu` depends on `lsdtt3_core` and `lsdtt3_utils`. All CLI executables
link against `lsdtt3_cpu` and `lsdtt3_core`. OpenMP, FFTW, NetCDF, and Eigen3
link to `lsdtt3_cpu`.


**Dependencies:**
- GDAL — required
- OpenMP — optional
- FFTW3 — optional (enables `lsdtt-spectral-analysis`)
- Eigen3 — optional (enables `lsdtt-chi-analysis` and `lsdtt_chi_ksn_inversion_cpu`)
- CUDA — optional


**FFTW:** Uses double API (`fftw_*`), NOT threaded FFTW, NOT float (`fftwf_*`).
Compile definition: `LSDTT_HAS_FFTW=1`. Stub in `lsdtt_raster_spectral_cpu.cpp`
throws at runtime if FFTW not available.

**`docs/developer_notes/lsdtt_utils_guide.md`**  
Developer guide to the `src/utils/` layer. Explains what belongs in utils,
what should remain in `core` / `cpu` / `gpu`, the valid-data contract for
utility functions, and when to use the main statistical tests and distance
measures.

**`docs/developer_notes/python_bindings_architecture.md`**
explains the in-memory Python backend, public API contract, nanobind adapters,
manifest/check tooling, bundled CLI drivers, and the safe binding change
workflow.


---

## Naming conventions

- All snake_case, `lsdtt_` prefix on all classes and files
- No camelCase anywhere
- Opening braces always on their own line
- All parameters use snake_case with no spaces

**CLI parameter names** are governed by
[`cli_parameter_conventions.md`](cli_parameter_conventions.md), which is the
authoritative standard: read it before naming or renaming any parameter. The
rules most often needed are:

- Output raster switches: `write_<name>_raster` (`_rasters` where a set is written)
- Output switches are `write_<output>`, never `print_` for something that writes a file
- One raster-format parameter per CLI, `write_output_format`, and one vector-format
  parameter, `write_vector_format`
- Algorithm selection parameters: `*_algorithm_choice`
- Units spelled out: `_metres`, `_pixels`, `_degrees` (never `_m`, `_px`, `_deg`)
- Bounds: `min_` / `max_` prefixes

The mechanically checkable rules are enforced by the CTest case
`unit_cli_parameter_conventions` (`tests/scripts/check_cli_parameter_conventions.py`),
which reads every built CLI's `-h json`. Its baseline allow-list is now empty,
so it fails on every violation. An extra format selector that genuinely governs
a different output domain (for example diagnostic frames rather than analysis
rasters) must be added to `REVIEWED_FORMAT_OVERRIDES` in that script, with a
comment saying why. See the "Conformance testing" section of the conventions
document.

**CLI executables use dashes, not underscores:**
- Object/library: `lsdtt_flowinfo_cpu`, `lsdtt_raster_cpu`, `lsdtt_param_parser`
- CLI executable: `lsdtt-channel-extraction`, `lsdtt-raster-preprocessing`

**Output filename convention:**
`<write_fname>_<descriptor>.tif`. If `write_fname` not set, defaults to stem
of `read_fname`. E.g. input `my_DEM.tif` → `my_DEM_n_contributing_pixels.tif`.

**Projection suffix convention** (applied by `lsdtt_raster_reproject`):
- `_utm`     — UTM zone of raster centre
- `_polar`   — NSIDC polar stereographic (|lat| > 75°)
- `_lambert` — custom Lambert Azimuthal Equal Area (wide-extent rasters, opt-in)

**File suffix convention:**
CPU-specific analysis files use `_cpu` suffix: `lsdtt_raster_spectral_cpu.hpp`.
The CPU raster class (`lsdtt_raster_cpu`) and CPU flow routing (`lsdtt_flowinfo_cpu`)
do not use the suffix because they are the primary named objects, not variants.

---

## Architecture design principles

**`src/core/`** — abstract interfaces, factories, shared types, I/O.
Nothing backend-specific. No CPU arrays. No CUDA.

**`src/cpu/`** — everything that runs on host CPU: the CPU raster class,
all CPU analysis free functions, and the CPU flow routing object.

**`src/gpu/`** — GPU implementations. Currently holds only
`lsdtt_surface_gpu.{hpp,cu}`, which has not been compiled or parity-tested on
real hardware. There is no `lsdtt_raster_gpu`; `lsdtt_backend_factory.cpp`
includes such a header under `#ifdef LSDTT_GPU_ENABLED`, a macro CMake never
defines. See [`gpu_compute.md`](gpu_compute.md) for what is actually built.

**`src/utils/`** — header-only utilities with no backend dependency.

**Abstract interface pattern:**
```
src/core/lsdtt_raster_base    ← abstract interface for raster backends
src/cpu/lsdtt_raster_cpu      ← CPU implementation
src/gpu/lsdtt_raster_gpu      ← GPU implementation (not written)

src/core/lsdtt_flow_base      ← abstract interface for flow routing
src/cpu/lsdtt_flowinfo_cpu    ← CPU implementation
src/gpu/lsdtt_flowinfo_gpu    ← GPU implementation (not written)
```

**Factory pattern:**
```
lsdtt_backend_factory  →  lsdtt_raster_cpu (the gpu branch is not built)
make_flow_base()       →  lsdtt_flowinfo_cpu (lsdtt_flowinfo_gpu not written)
```

**CLIs hold abstract interfaces via unique_ptr — never concrete types:**
```cpp
// Flow routing: include lsdtt_flow_base.hpp, never lsdtt_flowinfo_cpu.hpp
auto flow = lsdtt::make_flow_base(dem_data, dem_mask, meta, opts);
// returns unique_ptr<lsdtt_flow_base>
```

Any CLI written against `lsdtt_flow_base` works without modification when the
GPU flow backend is added — only the factory changes.


### Backend abstraction rules

`lsdtt3` uses **two backend-abstraction patterns**. See
`docs/developer_notes/backend_abstraction_rules.md` for the full design notes, rationale, and examples.

1. **Selective dispatch for broad compute substrates**  
   `lsdtt_raster_cpu` is a special case and remains a concrete CPU-side working
   object. Most raster analyses stay as CPU free functions taking
   `const lsdtt_raster_cpu&` and returning `lsdtt_raster_cpu`. Only the small
   subset of heavy kernels that are realistic GPU targets should use dispatch
   wrappers. Dispatch declarations live in `src/core/`; implementations live in
   `src/cpu/`; optional GPU kernels live in `src/gpu/`.

2. **Narrow abstract bases for persistent stateful engines**  
   Flow-like engines should use abstract bases in `src/core/` plus factories
   implemented in `src/cpu/`, so CLIs and workflows hold
   `std::unique_ptr<..._base>` rather than concrete CPU classes. This applies to
   `lsdtt_flow_base`, should also apply to `lsdtt_junction_network_base`, and is
   the intended long-term pattern for the basin-manager and particle-tracer
   layers. These interfaces must remain narrow, stable, and behaviour-oriented.

In short: **rasters stay concrete and selectively dispatched; stateful engines
use narrow abstract bases and factories.**


### GPU compute kernels (selective dispatch)

The selective-dispatch pattern (#1 above) is realised for GPU offload of heavy
per-pixel free functions. The first ported kernel is `direct_polyfit`; more
surface-fitting / roughness kernels follow the same template. Full design notes
are in `docs/developer_notes/gpu_compute.md`; the essentials:

**Two independent GPU axes — do not conflate them:**

| Axis | Macro | What it is | Status |
|---|---|---|---|
| GPU **compute kernels** | `LSDTT_GPU_KERNELS` | dispatch kernels for heavy free functions (`direct_polyfit`, future roughness metrics) | being built |
| GPU **raster-engine backend** | `LSDTT_GPU_ENABLED` | full `lsdtt_raster_base` GPU impl returned by `make_raster_engine()` | future, not built |

Both are gated by the single CMake option `-DLSDTT_GPU=ON`, but currently only
the kernel library (`lsdtt3_gpu`) is built and only `LSDTT_GPU_KERNELS` is
defined.

**Global toggle (a preference, not a mode):**
`set_compute_backend()` / `get_compute_backend()` live in
`src/core/lsdtt_backend_factory.{hpp,cpp}` as pure host state (no CUDA in core).
Each dispatched function resolves it as:

```
if (preference != cpu) AND (a GPU kernel exists for THIS function)
   AND (a device is present)   ->  run the GPU kernel
else                           ->  run the CPU implementation
```

Functions without a GPU port ignore the toggle and run on CPU; file I/O
(`lsdtt_data_loader`) never consults it; a failed GPU launch falls back to CPU
automatically, so correctness never depends on the GPU path.

**Dependency direction and seam placement:**

```
src/gpu/lsdtt_surface_gpu.hpp    POD launcher interface (raw pointers + dims,
                                 no CUDA types, no lsdtt_raster_cpu, no GDAL)
src/gpu/lsdtt_surface_gpu.cu     CUDA kernels + launchers   (lib lsdtt3_gpu)
src/cpu/..._surface_metrics_cpu.cpp
                                 direct_polyfit(): #ifdef LSDTT_GPU_KERNELS
                                 dispatch block at top, else CPU body
```

The dispatch decision is made in `lsdtt3_cpu` (which links `lsdtt3_gpu`), keeping
`src/core/` CUDA-free. The kernel interface is plain POD pointers, so the
one-way dependency is `lsdtt3_cpu → lsdtt3_gpu` with no reverse coupling.

**Parity rule:** GPU kernels accumulate OLS / reduction moments in **double**
and reuse the CPU's exact algebra (e.g. the symmetric-window decoupled normal
equations). A tolerance-based CPU-vs-GPU comparison test must pass on hardware
before a GPU kernel is trusted. The `.cu` is design-ahead until then; the default
`-DLSDTT_GPU=OFF` build compiles no CUDA.

**Porting another kernel:** add a POD launcher (`.hpp`) + kernel (`.cu`), add an
`#ifdef LSDTT_GPU_KERNELS` dispatch block at the top of the CPU free function,
add a parity test. No CMake change (the kernel joins the existing `lsdtt3_gpu`
sources).


### Returning rasters from core-layer engine objects

**Problem:** abstract base classes in `src/core/` sometimes need to return a
raster (e.g. a flow-direction raster, a channel-mask raster). They cannot
include `lsdtt_raster_cpu` because that would create a forbidden `core → cpu`
dependency.

**Solution: `lsdtt_raster_result`** — a plain data-transfer struct in
`src/core/lsdtt_raster_result.hpp`:

```cpp
struct lsdtt_raster_result
{
    std::vector<float>   data;
    std::vector<uint8_t> mask;
    lsdtt_raster_meta    meta;
};
```

Base-class virtual methods that produce rasters return `lsdtt_raster_result`.
The CPU implementation fills it from its internal storage and returns it.

**Two patterns for callers (CLIs, workflows):**

```cpp
// 1. Write directly to disk — no lsdtt_raster_cpu needed
lsdtt_raster_result r = engine->get_some_raster();
loader.write_full(r, out_path, format);

// 2. Promote to lsdtt_raster_cpu for further analysis (move — O(1))
lsdtt_raster_result r = engine->get_some_raster();
lsdtt_raster_cpu raster(std::move(r));
auto hillshade = lsdtt::compute_hillshade(raster);
```

**Choosing between `lsdtt_raster_result` and a `write_*(path, format)` method:**

- If the CLI only ever writes the raster to disk and never analyses it,
  a `write_*(path, format)` virtual method is simpler — the implementation
  creates the data internally and calls `loader.write_full()` directly.
  `print_arcmap_d8_flowdirection` on `lsdtt_flow_base` is an example of this.

- If the CLI may want to analyse the raster (pass it to an algebra function,
  compute derivatives, etc.), the method should return `lsdtt_raster_result`
  so the caller can promote it to `lsdtt_raster_cpu`.

**Convention:** pixels with `mask == 0` must already carry `meta.no_data_value`
in the `data` array. `lsdtt_data_loader::write_full` does not consult the mask.


---

## The GDAL boundary

GDAL, OGR and CPL headers appear **only in `.cpp` files under `src/`**. They are
not reachable from any lsdtt3 header, and not from `cli/` at all.

`unit_gdal_init_boundary` (`tests/scripts/check_gdal_boundary.py`) enforces this.
It must be a source check rather than a runtime one: a violation compiles
perfectly well and simply undoes the arrangement silently.

**Why it matters.** `lsdtt_raster_meta.hpp` is included by essentially every
compute object, so when it included `<ogr_spatialref.h>` the OGR headers reached
several hundred translation units. That cost compile time, and it meant any file
could reach for GDAL without declaring that it had.

**Three ways a header can refer to a GDAL type without including one**, all
demonstrated in `lsdtt_data_loader.hpp`:

1. **An lsdtt-native equivalent** — `lsdtt_pixel_type` in place of
   `GDALDataType`.
2. **A forward declaration** — `class GDALDataset;`. A `unique_ptr` needs the
   complete type only where the deleter runs, so defining
   `gdal_dataset_deleter::operator()` in the `.cpp` is sufficient. No pimpl is
   required as long as the destructor and move operations are also out of line.
3. **An opaque pointer** — `write_record_to_layer(..., void* ogr_layer)`.

**If a CLI needs GDAL, a library function is missing.** Add it to `src/core/`
and call that instead. `lsdtt-fetch-raster` used to build
`OGRCoordinateTransformation` objects by hand; it now uses
`lsdtt_coordinate_transform`, which already existed.

---

## lsdtt_gdal_init

File: `src/core/lsdtt_gdal_init.hpp/.cpp`.

`lsdtt::gdal_init()` is the **only** caller of `GDALAllRegister()` in the tree.
It is `std::call_once`-guarded, so it is cheap to call defensively at the top of
any function that touches GDAL — and every such function does, instead of
calling `GDALAllRegister()` itself.

Beyond registering drivers it pins the runtime data directories:

- Probes `<exe>/../share/proj`, `<exe>/share/proj`, `<exe>/../../share/proj`
  for `proj.db`, and pins the first hit with `OSRSetPROJSearchPaths()`.
- Does the same for `<exe>/../share/gdal` (sentinel `gdalvrt.xsd`) via
  `CPLSetConfigOption("GDAL_DATA", …)`.
- **A bundled copy always wins**, including over `PROJ_DATA` / `PROJ_LIB` /
  `GDAL_DATA`. Only when nothing is bundled next to the binary does it fall
  through to GDAL's own search — which is the correct behaviour for a
  conda-forge install, where sharing one GDAL is the point.

**Why this exists.** A mismatched `proj.db` does not always fail loudly. It can
resolve a datum differently and return coordinates that are plausible and wrong.
`-v` therefore reports the GDAL and PROJ versions and which data directories are
in force; that is the first thing to ask for in a projection bug report.

`set_bundled_data_root()` overrides the probe root before first init, for hosts
where `executable_path()` is not the right anchor — the Python extension module
being the case that needs it.

The executable-path resolution itself is in `src/core/lsdtt_exe_path.hpp`,
shared with the cosmo tools' `crn_parameters.nc` lookup.

---

## lsdtt_raster_meta

Plain struct in `src/core/lsdtt_raster_meta.hpp`.
Requires: `<cmath>`, `<stdexcept>`, `<string>`, `lsdtt_pixel_type.hpp`.

**No GDAL header.** Every compute object in the tree carries an
`lsdtt_raster_meta`, so a GDAL include here would reach most of the codebase.
The CRS methods below are declared in the header and defined in
`lsdtt_raster_meta.cpp`, which is where `<ogr_spatialref.h>` lives.

Key fields: `n_rows`, `n_cols`, `x_min`, `y_min`, `cell_size`, `no_data_value`,
`crs_wkt`, `epsg_code`, `data_type`.

`data_type` is an `lsdtt_pixel_type` (`src/core/lsdtt_pixel_type.hpp`), not a
`GDALDataType`: `unknown`, `uint8`, `int8`, `uint16`, `int16`, `uint32`,
`int32`, `float32`, `float64`. The GDAL mapping lives in
`lsdtt_pixel_type_gdal.hpp`, which only `.cpp` files at the GDAL boundary may
include. `tiff_predictor_for()` also lives with the enum, since the right
TIFF predictor is a property of the sample type.

**CRS / georeferencing methods:**
- `set_crs_from_epsg(epsg)` — sets `crs_wkt` + `epsg_code` from an EPSG code (via OGR)
- `set_default_synthetic_georeference()` — gives a generated raster (no input DEM)
  a sensible default: the grid is centred on Teanu Island, Solomon Islands
  (11.633°S, 166.967°E — the real-world location of Zwift's "Watopia"), in UTM
  zone 58S (EPSG:32758), away from the antimeridian and poles so coordinates
  round-trip cleanly to/from lat/lon. Used by the synthetic-surface generators,
  `lsdtt_synthetic_network_cpu`, and `lsdtt-lem` when no CRS is supplied.
- `is_geographic()` — true if lat/lon CRS
- `linear_unit_name()` — returns e.g. "metre", "foot"
- `is_in_metres()` — true if linear unit is metres (within 1e-6)
- `check_projected_metres()` — throws `std::runtime_error` with messages:
  - *"You are using a geographic coordinate system but this program only works
    on projected data. Please convert your raster to a projected coordinate
    system like UTM."*
  - *"Your units are not in metres. Are you using feet? Please convert your
    raster to a unit that is not only used by lunatics."*

---

## CLI executables

Current CLI programs are:

- **`lsdtt-raster-preprocessing.cpp`**  
  Raster preprocessing workflow: projection enforcement/reprojection, masking,
  nodata filling, breaching, pit filling, hillshade, and related outputs.

- **`lsdtt-synthetic-surfaces.cpp`**
  Synthetic raster generation and testing utilities.

- **`lsdtt-synthetic-networks.cpp`**  
  Generates synthetic Optimal Channel Networks (OCNs) with tunable branching
  statistics (tunable Tokunaga exponents) via simulated-annealing edge flips that
  minimise `E = Σ Aᵢ^gamma` (Metropolis acceptance `exp(-ΔE/T)`; the cooling
  schedule sweeps from random/Shreve-like at high T to strongly aggregated as
  T→0). Engine: `lsdtt_synthetic_network_cpu`. Elevations use the **chi-based
  fluvial steady state** via the *same* shared core
  (`compute_fluvial_steady_state_on_tree`) as lsdtt-lem's steady-state solver, so
  an OCN and a LEM steady state on the same graph give the same surface; relief
  is controlled exactly as in the LEM (`steady_state_mode` K/relief,
  `target_relief`, `uplift_rate`, `erodibility`, `m_exponent`, `n_exponent`,
  `baselevel_elevation`); fluvial relief only. The domain can be the full
  rectangle or an **irregular island** (`domain_type` = efd_island /
  noise_island / polygon): EFD = Elliptical Fourier random coastline,
  noise_island = thresholded radial red-noise, polygon = a user vector
  (`polygon_sizing` = scaled to the grid, or true_area = grid derived from the
  polygon's real extent + cell_size, preserving area and the polygon's CRS)
  (`lsdtt_synthetic_domain_cpu`). Outlets follow `outlet_mode`
  (`auto` → islands drain to their whole coastline (all_coast), the rectangle to
  a single outlet). Generated DEMs carry the default synthetic georeference (see
  below: centred on "Watopia"/Teanu Island, UTM zone 58S), with sea pixels as
  nodata. Outputs a synthetic DEM, the OCN drainage-area raster, and a summary
  CSV that includes the **D8 round-trip fidelity** (fraction of interior pixels whose D8 receiver equals the OCN
  receiver — ~1 for relaxed low-T networks, lower for random high-T networks, an
  inherent OCN self-consistency property).

- **`lsdtt-surface-metrics.cpp`**  
  Surface-analysis workflow that does not require flow routing: hillshade,
  polynomial surface fit, Otsu segmentation, Hough line extraction. Also writes
  a TopoCurve-style morphological surface-shape classification (SMAP) via
  `classify_surface_shape`: the polyfit principal curvatures are thresholded by
  `shape_curvature_threshold` (kt) and each pixel is labelled by the signs of
  Gaussian and mean curvature (peak / saddle / antiform / synform / basin /
  plane). Output: `_<method>_shape_class.tif`.

- **`lsdtt-drainage-analysis.cpp`**  
  Computes drainage-related rasters: optionally fills pits, then computes drainage
  area using D8 (`lsdtt_flowinfo`), D-infinity, or MFD (`lsdtt_flowgraph_cpu`),
  and derives the Topographic Wetness Index (TWI).

- **`lsdtt-internal-drainage.cpp`**  
  Detects closed depressions, round-trips sink candidates through CSV or vector
  review, and routes accepted basins either to non-destructive forced outlets or
  legacy-compatible nodata holes. Writes sink labels, drainage area, statistics,
  basin polygons, optional sink channel networks, and full-raster chi.

- **`lsdtt-channel-extraction.cpp`**  
  Builds or reuses a pit-filled DEM, constructs D8 flow routing, extracts source
  pixels from a contributing-pixel threshold, builds the junction/channel
  network, and writes channel-network output to CSV, GeoPackage, or FlatGeobuf.

- **`lsdtt-basin-extraction.cpp`**  
  Builds D8 flow routing, extracts a channel/junction network, and optionally
  selects and writes junction-aware basins (outlines, raster labels, statistics).

- **`lsdtt-basin-stats.cpp`**  
  Reads a DEM and a CSV of basin outlet points (latitude, longitude; optional gauge_id,
  gauge_name columns) and computes a suite of hydromorphic statistics per basin:
  elevation percentiles, slope (gradient magnitude) percentiles, flow-distance-from-outlet
  percentiles, centroid lat/lon, drainage area, nodata-edge flag, and channel length by
  Strahler stream order (orders 1–6). Basins touching nodata or the raster edge are
  optionally rejected (reject_nodata_edge_basins, default true).
  Optional toggles (each independent): `calculate_junction_area_ratios`
  (network-wide junction area-ratio CSVs) and `calculate_centre_of_mass_bearing`
  (adds com_bearing_deg / com_distance_m — bearing from the user outlet to the
  basin centre of mass; the outlet is used directly, the channel is not followed).
  Optionally (independent toggle `calculate_shape_metrics`, like `calculate_channel_stats`)
  computes per-basin planform shape metrics from the traced outline polygon
  (`compute_basin_shape` in `lsdtt_polygon_geometry.hpp`): elongation ratio, form factor,
  circularity ratio, Gravelius compactness, solidity, moment eccentricity, orientation, and
  hypsometric integral. Basin length L_b is selected by `basin_length_definition`
  (`flow_path` default, or `geometric` = max Feret diameter); the perimeter for
  circularity/Gravelius is measured on a Douglas–Peucker-simplified ring
  (`shape_perimeter_simplify_cells`) to mitigate the staircased divide-grid boundary.
  Theory and references: `docs/tex_technical/basin_shape_metrics.tex`.
  Also has an independent `calculate_junction_area_ratios` switch that writes the
  per-junction drainage-area ratio table + distribution summary (see below).
  `write_basin_polygons` writes one basin-outline polygon per basin to a vector
  layer (`write_vector_format` flatgeobuf/gpkg) with the per-basin hydromorphic,
  shape, and centre-of-mass metrics attached as typed attributes — i.e. the basin
  metrics joined to their basin polygons.

- **`lsdtt-feature-finder.cpp`**  
  Feature-detection workflows on DEMs (e.g. terrace, piedmont, hillslope morphology).

- **`lsdtt-mountain-front.cpp`**  
  Identifies mountain fronts using multi-scale relief, piedmont surface fitting, and
  residual thresholding.  Supports two input modes:
  - Synthetic DEM (use_synthetic_dem: true) — generates a test range from one of four
    geometries (linear, halfspace, circular, arc) and writes the ground-truth front_mask
    alongside the detected mountain_mask for algorithm validation.
  - Real DEM (use_synthetic_dem: false) — loads a projected DEM and runs the pipeline.
  Outputs: low_ground_mask, fitted_surface, residual, mountain_mask, and optional
  per-scale relief diagnostic rasters.

- **`lsdtt-flat-surface-finder.cpp`**  
  **Pure morphometry: which ground is flat?** Identifies flat (low-texture)
  surfaces using MADk2, residual IQR, and/or NVD metrics. Outputs per-metric
  texture rasters, binary masks, and a combined mask. Thresholds are set by
  percentile; masks can be combined with AND or OR logic. Also offers plane
  analysis and the VEC break-of-slope skeleton as diagnostics.

  No flow routing, no channel network, no height reference — everything that
  needed one moved to `lsdtt-fluvial-surfaces`. Its `_combined_flat_mask.tif`
  is that tool's `flat_mask_fname` input.

- **`lsdtt-fluvial-surfaces.cpp`**  
  **Which fluvial level is each flat surface on?** Takes a flat mask (or falls
  back to a built-in relief/slope admission), resolves a channel scope, builds a
  REFERENCE FIELD, bands it, and labels connected components WITHIN each band.

  Two reference fields, one ladder. `hand` is height above the in-scope channel
  network — the right field for a terrace flight, where treads sit at a roughly
  constant height above the channel. `fitted_surface` is the signed residual to
  a surface fitted to the flat ground — the right field for a fan, which climbs
  continuously in HAND and would be cut into meaningless arcs by a HAND ladder.
  The fitted surface may be a `cone` (the default, and the fan option), a
  RANSAC `plane`, or a `poly2`.

  **Band first, then label.** Two pixels can only share a label if they share a
  band, so two treads at different heights cannot merge — and the riser between
  them never has to be found or be continuous. This replaces four separate
  repair mechanisms that existed to detect that merge after the fact.

  **There is no valley mode and no terrace mode**, only a level count: the
  valley floor is the lowest band of a two-band cut, a terrace flight is an
  n-band ladder on the same population. Set `ladder_level_count` to 2, to n, or
  leave it at 0 and let the resolvability floor decide how many levels the DEM
  can carry — reporting UNRESOLVED where it cannot.

  Outputs `_band_index.tif`, `_banded_labels.tif`, `_valley_mask.tif` (band 0,
  the input to `lsdtt-valley-centreline`), `_surface_class.tif`,
  `_surface_patches.csv` (the join key for everything downstream) and
  `_ladder_windows.csv` (the resolvability report — read it first).
  `hand_fname` lets a precomputed HAND be supplied, which is what makes
  iterating on the level count affordable.

- **`lsdtt-graphflood.cpp`**
  Graph-based hydraulic routing on a fixed D-inf or MFD flow graph. Constructs
  `lsdtt_hydraulic_graph_cpu`, runs N timesteps, and writes water-depth,
  water-surface, and flux rasters (Gailleton et al. 2024).

- **`lsdtt-lem.cpp`**  
  Simple landscape evolution model (evolved from MuddPILE / DONKmodel). Builds a
  synthetic initial surface (flat / white-noise / red-noise / diamond-square) or
  loads a DEM, then runs one of two paths:
  - **Time-stepping LEM** — loops `run_lem_step`: implicit stream-power fluvial
    incision (D8 or D-inf), optional hillslope transport (cidre / roering / linear /
    adi / roering_implicit — the unconditionally-stable backward-Euler/TR-BDF2 solver,
    with channels as Dirichlet sinks and an optional `channel_mask_fname`), and uniform
    uplift, writing per-step elevation stats and optional intermediate rasters.
  - **Analytical steady state** (`run_steady_state: true`) — chi-based fluvial snap
    (optionally iterated over `steady_state_n_cycles` to let divides migrate), then
    optional hillslope steady-state solve (SOR linear/roering or pseudo-transient
    Cidre) with channels as Dirichlet BCs. Exits after writing.
  Supports sinuous domain masking, N-S parabola conditioning, breaching, and a
  `png` output format with worldfile/PAM sidecars. See the LEM stack section below.

- **`lsdtt-channel-network-analysis.cpp`**  
  Builds a channel/junction network and runs one or more network analyses:
  channel-link bearing fitting (orthogonal regression), junction-angle calculation
  (donor opening angle, donor–receiver angles, bending angle),
  Strahler/Tokunaga self-similarity analysis (with a count-weighted Poisson GLM
  fit alongside the OLS fit), per-junction drainage-area ratios
  (`write_junction_area_ratio_data` / `_summary` — a continuous, high-sample-count
  complement to order-binned Tokunaga, with branch-bootstrap CIs and an
  area-decade collapse table), and per-junction basin centre-of-mass bearings
  (`write_junction_centre_of_mass_bearings` — for each junction the basin is the
  upslope area of the pixel just before the next downstream junction; outputs the
  outlet→centre-of-mass bearing). Writes CSV, GeoPackage, or FlatGeobuf.

- **`lsdtt-chi-analysis.cpp`** *(requires Eigen3)*  
  End-to-end chi workflow: chi extraction on threshold-defined channel networks,
  piecewise-linear chi-network fitting, disorder metrics with jackknife/bootstrap
  uncertainty, and spatial k_sn/u* inversion with m/n model selection.
  Outlets can be supplied as a CSV of lat/lon points; otherwise the largest basin is used.

- **`lsdtt-spectral-analysis.cpp`** *(requires FFTW3)*  
  Spectral surface analysis workflows using the FFTW double API.

- **`lsdtt-stream-burning.cpp`**  
  Loads a polyline vector channel network and a DEM, finds degree-1 source tips in the
  vector network (`find_vector_network_source_tips`), optionally writes those tips as a QC
  point layer, then seeds a D8 junction/channel network from those tips. Writes the
  DEM-extracted network as lines and optionally as points and junctions. The two networks
  (vector-supplied vs DEM-extracted) are intended for subsequent comparison to decide
  where to burn the vector channel directly into the DEM vs follow DEM routing.
  With `carve_before_fill: true` it conditions the DEM by breaching first
  (`breach_algorithm_choice`: `lindsay2016` or `least_cost`; least-cost uses
  `breach_max_dist_pixels` (1000), `breach_max_cost`, `breach_minimize_dist`
  (false) and `breach_edge_policy` (outlet), validated before the DEM is read).

- **`lsdtt-cosmo-tool.cpp`**  
  Cosmogenic nuclide erosion-rate inversion. Loads measured nuclide concentrations at
  channel sample points (CSV with latitude/longitude), computes the Stone 2000 production
  raster, snaps points to the D8 channel network, and inverts each measured concentration
  to a basin-averaged effective erosion rate using Newton-Raphson + bisection.
  Optionally applies AMS standardisation and uses ERA40 atmospheric / palaeomagnetic
  corrections when a `crn_parameters.nc` file is supplied.

- **`lsdtt-cosmo-prediction.cpp`**  
  Predicts cosmogenic nuclide concentrations across a DEM.  Three modes, selectable
  independently and combinable in a single run:
  - `calculate_crn_concentration_raster` — steady-state concentration at every pixel
    given an erosion-rate field and shielding rasters
    (`calculate_crn_concentration_raster_whole`).
  - `calculate_step_change_crn_concentration_raster` — transient concentration after
    an instantaneous step change in erosion rate; per-pixel `time_since_change` raster,
    uniform pre- and post-change effective erosion rates
    (`calculate_crn_concentration_raster_step_change`).
  - `calculate_accumulated_crn_concentration` — erosion-weighted average concentration
    accumulated downstream via D8 flow routing (`calculate_accumulated_crn_raster`).
  All shielding inputs (topographic, snow, self) and the erosion-rate raster can be
  loaded from disk or specified as uniform scalar defaults.
  Outputs: `_<nuclide>_crn_conc.tif`, `_<nuclide>_step_crn_conc.tif`,
  `_<nuclide>_crn_acc.tif`.

- **`lsdtt-cosmo-column-tool.cpp`**  
  1-D cosmogenic nuclide column simulation (no DEM required).  Three modes:
  `steady_state`, `step_change`, and `full_transient` (user-supplied erosion-history CSV).
  Uses `lsdtt_crn_column` and `lsdtt_crn_parameters`.

- **`lsdtt-cosmo-landslides.cpp`**  
  Cosmogenic nuclide analysis for earthquake-triggered landslides. Models CRN 
  concentration in landslide samples accounting for burial, incomplete exhumation, 
  and depth-integrated production. Supports ensemble models for undated samples.

- **`lsdtt-basin-surface-metrics.cpp`**  
  Per-basin surface metrics extraction. Computes hillshade, relief, roughness, 
  curvature, and texture metrics for each basin in a network, writing per-basin 
  statistics CSV and optional raster maps.

- **`lsdtt-swath-profile.cpp`**  
  Swath profile extraction. Computes summary statistics (min, max, mean, percentiles) 
  across a swath corridor perpendicular to a polyline baseline. Outputs swath 
  CSV and optional visualization rasters.

- **`lsdtt-lithology-erodibility.cpp`**  
  **Lithology map to an erodibility (K) raster for `lsdtt-lem`.** Reads a polygon
  map with a lithology field (e.g. a Macrostrat download) or an integer code
  raster, maps values to codes (crosswalk, integer, or table name), and writes K
  from a property table in the table's own convention, on a DEM's grid, a
  synthetic grid (with `lsdtt-lem`'s placements) or the map's own grid. Also
  writes the code raster, a legend and the `lsdtt-lem` parameter lines that state
  the table's calibration. Uses `lsdtt_lithology_map_cpu`,
  `rasterise_polygon_field_values` and `lsdtt_lem_input_grid_cpu`.

- **`lsdtt-cross-sections.cpp`**  
  Cross-section extraction along a polyline. Samples elevation, slope, and metrics 
  perpendicular to a baseline, outputs profiles as CSV or rasters.

- **`lsdtt-valley-centreline.cpp`**  
  **Trace a centreline through a supplied valley.** Takes a valley mask
  (`valley_mask_fname`, REQUIRED — normally `_valley_mask.tif` from
  `lsdtt-fluvial-surfaces`), carves a trough, and traces the medial-axis /
  least-cost centreline through it. Also fits the pinned monotonic thalweg
  baseline, extracts cross sections, and derives cross-section hydraulics.

  It does **not** identify the valley. That moved to `lsdtt-fluvial-surfaces`,
  where a valley is simply band 0 of a two-band ladder — a formulation that also
  works on a fan, which a corridor HAND break does not. Running without a mask
  is a hard error naming that tool.

  `refine_valley_with_baseline` refines the supplied mask against the fitted
  baseline and writes the result for the user; it is never fed back to
  `lsdtt-fluvial-surfaces`. It requires an explicit `refine_max_height_metres`,
  because there is no longer a derived break to reuse.

- **`lsdtt-water-centreline.cpp`**  
  Water/floodplain centreline extraction. Combines multiple methods to identify 
  river and floodplain centrelines. Outputs polyline networks.

- **`lsdtt-water-line-repair.cpp`**  
  Repairs vector water-line networks. Fixes topological issues, gaps, and 
  inconsistencies in digitised stream networks.

- **`lsdtt-ridgeline-analysis.cpp`**  
  Ridgeline extraction and analysis. Identifies mountain ridges from high-curvature 
  topography. Outputs ridgeline polylines and statistics.

- **`lsdtt-divide-migration.cpp`**  
  Divide migration analysis. Tracks changes in drainage divides through time or 
  across parameter space. Uses multiple methods to quantify divide movement.

- **`lsdtt-sediment-routing.cpp`**  
  Sediment transport and routing through a flow network. Models detachment-limited 
  and transport-limited sediment flux.

- **`lsdtt-linear-stats.cpp`**  
  Linear feature statistics. Extracts statistics from linear features (channels, 
  ridges) in raster or vector form, including bearing analysis and sinuosity.

- **`lsdtt-raster-calculator.cpp`**  
  Map algebra calculator. Element-wise operations on rasters: addition, subtraction, 
  multiplication, division, powers, and other mathematical functions.

- **`lsdtt-fetch-raster.cpp`**  
  Raster data fetching from remote sources. Downloads and processes elevation data 
  from cloud providers (COP30, Planetary Computer, etc.) with optional tiling, 
  reprojection, and mosaicking.

- **`lsdtt-fetch-lithology.cpp`**  
  **Geological map units from Macrostrat for a bounding box.** Reads Macrostrat's
  `carto` vector tiles (GDAL MVT driver), keeps the units overlapping the box,
  looks up each map unit's linked lithologies and source maps in the API
  (optional; a failed lookup is a warning), and writes a polygon layer in the
  local UTM zone with `lith_primary` for `lsdtt-lithology-erodibility`, plus a
  lithology vocabulary, a crosswalk template and the CC BY 4.0 attribution.
  Shares the bounding-box, bounds-file and HTTP parameters with
  `lsdtt-fetch-raster` (`lsdtt_common_parameters`). Uses `lsdtt_macrostrat_fetch`.

- **`lsdtt-fetch-erosion-rates.cpp`**  
  **Basin-averaged cosmogenic erosion rates from OCTOPUS for a bounding box.**
  Counts the intersecting basins on the keyless WFS (`resultType=hits`), pages
  the basins and outlets, joins them on `OBSID1`, nulls the no-data values, and
  writes an erosion-rate sample CSV (outlet, rate, 1-sigma uncertainty, units,
  reported area), the basin polygons, a suggested DEM box and the CC BY 4.0
  attribution. Shares the bounding-box, bounds-file and HTTP parameters with the
  other fetch CLIs. Uses `lsdtt_octopus_fetch`.

- **`lsdtt-marsh-topography.cpp`**  
  Marsh and tidal platform analysis. Identifies marsh platforms and intertidal zones, 
  computes platform statistics, and extracts edges.

- **`lsdtt-slope-area.cpp`**  
  Slope-area statistics. Builds D-infinity (default) or D8 drainage area,
  computes direct-polyfit slope for each value in `polyfit_lambda_metres_list`,
  optionally clips to a basin or polygon mask, and writes one log-binned
  drainage-area CSV per lambda plus a summary of the drainage area at maximum
  median slope. Makes no plots.

- **`lsdtt-advection-diffusion.cpp`**  
  A deliberately small fixed-timestep landscape evolution model: red-noise
  initial surface, stream-power fluvial incision on D8 or D-infinity routing,
  linear ADI hillslope diffusion and uniform uplift. Optional hillshade and
  log10 drainage-area diagnostic frames (`write_hillshade_rasters`,
  `write_drainage_area_rasters`, `diagnostics_write_output_format`). Additional
  solvers, steady-state snapping, adaptive timesteps and forcing histories are
  reserved for `lsdtt-lem`.

- **`lsdtt-bankfull-width.cpp`**  
  Manning's-equation bankfull channel width from a high-resolution DEM. Extracts
  cross sections along a single-source channel network, then writes a channel
  profile, sampled cross-section elevations, Manning bank points and per-reach
  widths. Every output shares one cross-section extraction, so
  `cross_section_id` matches across files.

- **`lsdtt-distributed-hydrology.cpp`**  
  A raster TOPMODEL: a time-evolving, rainfall-driven simulation built on the
  same wetness index as `lsdtt-drainage-analysis`. The CLI contains no hydrology
  of its own. It assembles library modules (`lsdtt_wetness_index_cpu`,
  `lsdtt_topmodel_cpu`, `lsdtt_topmodel_dynamic_cpu`,
  `lsdtt_canopy_interception_cpu`, `lsdtt_infiltration_excess_cpu`,
  `lsdtt_hydro_routing_cpu`, ...). Routing is D-infinity or MFD through
  `lsdtt_flowgraph_cpu`.

- **`lsdtt-hillslope-channel-coupling.cpp`**  
  Couples hillslope form to channel form. For every ridgeline pixel it measures
  hilltop curvature, traces the flow path to the channel it drains to, and
  reports hillslope length, relief and the non-dimensional E* / R* pair
  (Hurst et al. 2012; Grieve et al. 2016).

- **`lsdtt-outlet-matching.cpp`**  
  Matches outlets on a reference DEM to cells on a target DEM by comparing D8
  contributing-basin areas, preferring the closest-area candidate that is not
  influenced by the DEM edge. Writes matched-outlet and diagnostic CSVs.

### Polyfit method rule

Wherever a CLI or algorithm lets the user choose a polyfit method, the options are:
- `direct` — rolling-moment separable fit (`direct_polyfit`); **default and preferred**; faster (O(n) per pixel) and better at resolving fine-scale features
- `wavelet` — FFT/Gaussian-smoothed fit (`wavelet_quadratic_fit`); available as a fallback

**Always default to `direct`.** Both methods return the same `poly_fit_coeffs` type and are interchangeable as input to `compute_polyfit_metrics`.

### CLI design rule

Every CLI uses `lsdtt_param_parser`. Use `cli/lsdtt-channel-extraction.cpp` and
`cli/lsdtt-surface-metrics.cpp` as canonical reference implementations.

CLI source files follow this layout, in order:

1. file-header comment block (copyright, brief description, outputs, usage)
2. `#include` block
3. `static const` status, citation, and tool-name strings (`k_status`,
   `k_citation`, `k_tool_name`)
4. `print_banner()` — fixed-width box to stdout. Required content: a
   brief summary of what the CLI does; then a versioning box, directly
   below the summary, containing `|| lsdtt3 version : ` followed by
   `lsdtt::lsdtt3_package_version()`, and `|| Status         : <k_status>`,
   with " -- not yet fully tested" appended when the status is not
   `stable`; then the citation and a line pointing to `-h` for
   parameter help. CLIs have no version of their own: lsdtt3 has one
   version number. `k_status` is `"experimental"` or `"stable"`; set it to
   `"stable"` only once a CLI is tested (a review decision; the guard in
   `unit_cli_parameter_conventions` checks only the value)
5. helper free functions (if unavoidable — see helper rule below)
6. `main()`:
   a. call `print_banner()` and `lsdtt_param_parser::parse_args_deferred()`
   b. declare all five maps and fill defaults + help (see below)
   c. if `cli_args.help_requested`, call
      `lsdtt_param_parser::print_terminal_help(...)` and return 0
      (help is rendered from the declarations, so no parameter file is needed)
   d. construct `lsdtt_param_parser` from `cli_args.path` / `cli_args.fname`
   e. call `parse_all_parameters`
   f. read results back from the parser
   g. cross-parameter validation / derived values
   h. load DEM and verify projection
   i. build flow routing and channel network (if needed)
   j. run analysis pipeline
   k. write outputs
   l. call `pp.print_help(...)` then `pp.print_parameters()`

All CLIs use this deferred flow (the legacy `parse_args()` was removed
2026-07-15 once the last CLI migrated). See `cli/lsdtt-raster-preprocessing.cpp`
and `cli/lsdtt-channel-extraction.cpp` for reference implementations.

**When adding or renaming a parameter:**

1. Name it according to [`cli_parameter_conventions.md`](cli_parameter_conventions.md).
2. If the CLI has been released, keep the old name working for one release
   through the deprecated/removed key maps (see
   [Deprecated, removed and unknown parameter keys](#deprecated-removed-and-unknown-parameter-keys)).
   A CLI not yet in any release tag can be renamed outright.
3. Rebuild the CLI and regenerate the CLI reference pages:

   ```bash
   python3 tools/generate_cli_reference.py --update --bin-dir build
   ```

   This rewrites the parameter tables between the
   `<!-- BEGIN/END GENERATED CLI PARAMETERS -->` markers in
   `docs/cli_reference/`. The CTest case `unit_cli_reference_is_current` runs
   the same script with `--check` and fails if any page is out of date. The
   hand-written tables and example parameter files outside the markers are not
   regenerated, so update the old names there by hand.

### Parameter map pattern (mandatory)

Every CLI declares exactly **five maps**:

```cpp
std::map<std::string, bool>        bool_default_map;
std::map<std::string, float>       float_default_map;
std::map<std::string, int>         int_default_map;
std::map<std::string, std::string> string_default_map;
lsdtt::lsdtt_help_registry         help_map;
```

(`lsdtt_help_registry` is a drop-in replacement for the older
`std::map<std::string, std::vector<std::string>>` that additionally records
declaration order so `-h` can group parameters in pipeline order.)

For every parameter, the default and the help entry are declared on **consecutive lines**:

```cpp
// --- Section label ---
float_default_map["my_param_metres"] = 10.0f;
help_map["my_param_metres"] = {"float","10.0",
    "Short one-line description including units [m].",
    "Longer explanation: valid range, physical meaning, cross-references to "
    "related parameters, caveats. Be specific about units and defaults.",
    "Section label", "param", "write_my_raster"};

bool_default_map["write_my_raster"] = false;
help_map["write_my_raster"] = {"bool","false",
    "Write the result raster.",
    "Output file: <write_fname>_my_result.tif.",
    "Section label", "switch", ""};
```

The `help_map` entry format is
`{type, default, short_desc, long_desc, group, role, parent}`:

- `type`: `"bool"`, `"int"`, `"float"`, or `"string"`.
- `group`: must match the `// --- Section label ---` comment above the block.
- `role`: `"switch"` (bool/enum that enables an analysis or output),
  `"param"` (tunes an analysis), or `"io"` (paths, formats, threading).
- `parent`: the switch key this parameter belongs to; `""` if none.

Entries with only the first four elements remain valid (legacy CLIs) but all
new parameters must declare all seven.

### Parsing and reading parameters back

```cpp
// Always call with this argument order: float, int, bool, string
pp.parse_all_parameters(
    float_default_map,
    int_default_map,
    bool_default_map,
    string_default_map
);

// Read results back into short-named maps for convenience
std::map<std::string, bool>        bf  = pp.get_bool_parameters();
std::map<std::string, float>       ff  = pp.get_float_parameters();
std::map<std::string, int>         ifs = pp.get_int_parameters();
std::map<std::string, std::string> sf  = pp.get_string_parameters();
```

Parameters are not declared `const` because cross-parameter logic may adjust them.

### Writing the help file and parameter echo

At the end of every successful run:

```cpp
pp.print_help(help_map, pp.get_write_fname() + "_README",
              k_status, k_citation);
pp.print_parameters();
```

Note: `print_help` prepends `write_path_` internally, so pass
`get_write_fname()` (not `get_full_write_path()`, which would double the
path).

`print_help` writes a CSV and an HTML reference file alongside the outputs so
users can look up every parameter, its type, default, and description without
reading source code.

### Flow routing in CLIs

The concrete class defined in `lsdtt_flowinfo_cpu.hpp` is named **`lsdtt_flowinfo`**
(the file name and the class name differ — the class has no `_cpu` suffix).

CLIs that need only the abstract flow interface (e.g. writing flow-direction rasters,
then dynamic-casting for threshold source extraction) should use the factory:

```cpp
// Include lsdtt_flow_base.hpp, not lsdtt_flowinfo_cpu.hpp
auto flow = lsdtt::make_flow_base(dem_data, dem_mask, meta, fi_opts);
// returns unique_ptr<lsdtt_flow_base>

// Cast when concrete API is needed
const lsdtt::lsdtt_flowinfo* flow_cpu =
    dynamic_cast<const lsdtt::lsdtt_flowinfo*>(flow.get());
```

CLIs that need the full `lsdtt_flowinfo` API throughout (e.g. sediment routing,
LEM) may construct the concrete object directly:

```cpp
// Include lsdtt_flowinfo_cpu.hpp — the class inside is lsdtt::lsdtt_flowinfo
lsdtt::lsdtt_flowinfo flow(
    dem.get_data_array(), dem.get_mask_array(), dem.get_meta(), fi_opts);
```

Never write `lsdtt::lsdtt_flowinfo_cpu` — that type does not exist.

### Junction network in CLIs

Use the concrete class directly when topological-order and s_vector access are needed:

```cpp
lsdtt::lsdtt_junction_network_options jn_opts;
jn_opts.build_topological_order = true;
jn_opts.build_s_vector          = true;

lsdtt::lsdtt_junction_network_cpu jn(flow, source_flats, jn_opts);
```

Source flats come from `flow.get_sources_index_threshold(threshold)`.

### Helper function rule

Avoid helper free functions at the top of CLI files. If a helper is non-trivial,
it almost certainly belongs in a `.hpp/.cpp` pair in `src/cpu/` or `src/utils/`
because the same logic will likely be needed by another CLI. A function declared
at file scope in a CLI is invisible to the rest of the project.

### Comma-separated string parameters

Parameters such as grain-size vectors that are supplied as comma-separated values
should be stored as `string_default_map` entries and parsed with:

```cpp
std::vector<double> vals = pp.parse_double_vector("my_param");
```

The `parse_double_vector` / `parse_float_vector` / `parse_int_vector` helpers are
on `lsdtt_param_parser` and expect the string to already have been registered via
`string_default_map`.

### Shared parameter declaration layer

Files: `src/core/lsdtt_common_parameters.hpp/.cpp`

**Architectural role:** A CLI declares no shared parameter name of its own. Every
parameter shared by two or more CLIs is declared exactly once in the shared layer,
and CLIs adopt it via a function call.

This eliminates the root cause of parameter drift: each shared concept had been
re-typed independently by each CLI, so one author wrote 0 and another wrote 1
with no intent behind the difference, 20+ CLIs re-declared the same help text with
silent corrections, and threshold values diverged without reason.

**The pattern:** Each independently selectable semantic group gets one function
that populates a caller's default and help maps. Flow parameters are deliberately
granular because a CLI may support routing, one of two distinct fill concepts,
or channel initiation without supporting all of them:

```cpp
lsdtt::declare_routing_core_parameters(
    std::map<std::string, int>&         int_map,
    std::map<std::string, std::string>& string_map,
    lsdtt::lsdtt_help_registry&         help_map,
    const std::string& prefix = ""      // optional, for multi-stage CLIs
);

lsdtt::declare_raster_is_filled_parameter(bool_map, help_map);
lsdtt::declare_fill_before_routing_parameter(bool_map, help_map);
lsdtt::declare_fill_min_slope_parameter(float_map, help_map);
lsdtt::declare_channel_initiation_parameters(int_map, help_map);
```

On return, the maps contain only the selected shared parameters and their help
entries, ready for a CLI to call `parse_all_parameters()`.

**Adoption in a CLI:** Instead of hand-declaring shared parameters, a CLI
composes only the declarations for concepts it actually supports:

```cpp
lsdtt::declare_routing_core_parameters(
    int_default_map, string_default_map, help_map);
lsdtt::declare_raster_is_filled_parameter(bool_default_map, help_map);
lsdtt::declare_fill_min_slope_parameter(float_default_map, help_map);
lsdtt::declare_channel_initiation_parameters(int_default_map, help_map);
```

After adoption, renaming that parameter is a one-line change in one file, not a
twenty-CLI edit. This is why parameter renames are low-risk:
the shared layer makes them self-enforcing.

- **`lsdtt-distributed-hydrology.cpp`**  
  A raster TOPMODEL: a time-evolving, rainfall-driven simulation built on the
  same wetness index `lsdtt-drainage-analysis` reports. The driver contains no
  hydrology of its own — it assembles `lsdtt_wetness_index_cpu`,
  `lsdtt_landcover_scenario_cpu`, `lsdtt_channel_extent_cpu`,
  `lsdtt_hydro_roughness_cpu`, `lsdtt_hydro_velocity_field_cpu`,
  `lsdtt_canopy_interception_cpu`, `lsdtt_infiltration_excess_cpu`,
  `lsdtt_topmodel_cpu`, `lsdtt_topmodel_dynamic_cpu` and
  `lsdtt_hydro_routing_cpu` in order. Routing is D-infinity or MFD through
  `lsdtt_flowgraph_cpu`; D8 is deliberately not offered, because sending all of
  a cell's area to one neighbour stripes the wetness index and therefore the
  deficit field built on it.

- **`lsdtt-hillslope-channel-coupling.cpp`**  
  Couples hillslope form to channel form. For every ridgeline pixel it measures
  hilltop curvature, traces the flow path down to the channel that pixel drains
  to, and reports hillslope length, relief and the non-dimensional pair
  E\* / R\* (Hurst et al. 2012; Grieve et al. 2016). Divides are traced on the
  pixel-corner divide grid, so a divide runs *between* DEM pixels: hilltop
  pixels are the pixels a ridgeline separates, one either side, and because the
  ridgeline is ordered each pixel carries its position along the divide and
  which side it sits on. That side partitions a divide's hilltops into the two
  opposing hillslopes, so E\*/R\* can be compared across a divide rather than
  only pooled.

- **`lsdtt-outlet-matching.cpp`**  
  Matches reference-DEM outlets to target-DEM cells by comparing reference and
  target D8 contributing-basin areas. Loads both DEMs and a reference outlet
  CSV, validates CRS and transformed-extent overlap, builds D8 flow plus
  target contributing-pixel and edge-influence rasters, then selects the
  closest-area non-edge candidate. Writes matched-outlet and diagnostic CSVs.
  Uses the outlet-resolution stack (`lsdtt_outlet_bounds_discovery`,
  `lsdtt_outlet_clusters`, `lsdtt_outlet_resolution_cpu`,
  `lsdtt_outlet_basin_verification_cpu`).

**Migration status:** Not all CLIs have adopted the shared layer yet. Legacy CLIs
that predate the shared layer still hand-declare their parameters, and are migrated
one at a time.

**Naming standard:** Every shared parameter follows the conventions in
`docs/developer_notes/cli_parameter_conventions.md`. That document is the durable record of what
names are canonical, what units they use, and when to use `min_` vs `minimum_`,
`threshold_` vs `min_`, and so on.

---

## lsdtt_crs_suggest

Files: `src/core/lsdtt_crs_suggest.hpp/.cpp`

Inspects a `lsdtt_raster_meta` and returns a `lsdtt_crs_suggestion` struct
describing the appropriate projected CRS for analysis. Does **not** reproject —
that is the job of `lsdtt_raster_reproject`.

**Decision logic:**
1. Already projected (non-geographic, non-Web Mercator) → return `already_projected`
2. All four corners are unprojected to WGS84 lon/lat via `OGRCoordinateTransformation`
3. `|centre_lat| > 75°` → NSIDC polar stereographic (EPSG:3413 north / EPSG:3976 south)
4. Otherwise → UTM zone of centre longitude + hemisphere
5. If raster spans ≥ 3 UTM zones: set `warn_multi_zone`, populate `laea_wkt`

**Key types:**
```cpp
enum class CrsRecommendationType {
    already_projected, polar_north, polar_south, utm, laea_custom
};

struct lsdtt_crs_suggestion {
    CrsRecommendationType type;
    std::string recommended_wkt;    // WKT for primary recommendation
    std::string recommended_name;   // human-readable label
    int         recommended_epsg;   // 0 for custom CRS
    double centre_lon, centre_lat;
    double min_lon, max_lon, min_lat, max_lat;
    int    utm_zone_min, utm_zone_max, utm_zone_count;
    bool   already_projected;
    bool   warn_multi_zone;
    std::string warning_message;
    std::string laea_wkt;           // populated when warn_multi_zone && offer_laea
    std::string laea_name;
};
```

**Public API:**
- `lsdtt_suggest_crs(meta, offer_laea=true)` → `lsdtt_crs_suggestion`
- `lsdtt_print_crs_suggestion(sug)` — pretty-print to stdout

Web Mercator (EPSG:3857) is treated as unsuitable alongside EPSG:4326 — both
trigger the projection recommendation path.

---

## lsdtt_raster_reproject

Files: `src/core/lsdtt_raster_reproject.hpp/.cpp`

Reprojects a raster file on disk to a suitable projected CRS using `GDALWarp()`.
Calls `lsdtt_suggest_crs()` internally. Returns a `lsdtt_reproject_result`
containing the output path, the CRS applied, and the grid spacing used.

**Key types:**
```cpp
struct lsdtt_reproject_options {
    bool        enforce_projection      = true;
    std::string resampling_method       = "cubic";
    double      grid_spacing            = -1.0;  // ≤ 0 = auto
    bool        allow_lambert_projection = false;
};

struct lsdtt_reproject_result {
    bool        reprojected;
    std::string out_path;        // empty if !reprojected
    lsdtt_crs_suggestion suggestion;
    double      grid_spacing_used;
    std::string target_wkt;
};
```

**Grid spacing snapping** (`snap_grid_spacing(pixel_m)`):

| Estimated pixel size | Snapping rule |
|---|---|
| ≥ 20 m | nearest 5 m |
| 5–20 m | nearest 2 m |
| < 5 m  | nearest 1 m |
| rounds to 0 | preserve raw pixel size (sub-metre guard) |
| ≤ 0 after guard | clamp to 0.01 m (malformed metadata safety net) |

**Output naming:** `write_path + write_stem + suffix + ".tif"`
where suffix is `_utm`, `_polar`, or `_lambert`.

**Resampling methods** (any GDAL-accepted string):
`near`, `bilinear`, `cubic` (default), `cubicspline`, `lanczos`,
`average`, `mode`, `max`, `min`, `med`, `q1`, `q3`, `sum`.

**LAMBERT / LAEA path:**
Only activated when `allow_lambert_projection == true` AND the raster spans
≥ 3 UTM zones. The custom LAEA is centred on (centre_lon, centre_lat) of the
raster. Output suffix: `_lambert`.

**Public API:**
- `lsdtt_reproject_if_needed(in_path, write_path, write_stem, meta, opts)` → `lsdtt_reproject_result`
- `snap_grid_spacing(pixel_m)` → `double`

---


## lsdtt_raster_warp

Files: `src/core/lsdtt_raster_warp.hpp/.cpp` (header GDAL-free).

- `clip_raster_to_reference_extent` and `warp_raster_to_reference_grid`: file to
  file, against a reference raster.
- `warp_raster_to_grid_in_memory(in_path, target_grid, opts)`: warps onto an
  explicit grid description, so a grid that exists only in memory (a synthetic
  model grid) can be the target. It writes a plain GeoTIFF to `/vsimem/`, reads it
  back with `lsdtt_data_loader` and removes it. Shares its request with
  `warp_raster_to_reference_grid`.
- `crs_equal_semantically(a_wkt, b_wkt)`: public; OGR `IsSame` with traditional
  GIS axis order.

The `-te` extent is passed to GDAL at 17 significant digits. Until 2026-10-16 it
went through a default `ostringstream` (6 digits), which turned a northing of
5000202 into `5.0002e+06` and shifted any warp or clip with 7-digit coordinates
(every UTM northing) by up to several metres.

## lsdtt_raster_base

File: `src/core/lsdtt_raster_base.hpp/.cpp`

Pure dispatch interface. Does NOT expose `get_data()` or `get_mask()`, which
return CPU-specific `std::vector` types that have no meaning on a GPU backend.
Metadata is always host-side.

CPU analysis free functions take `const lsdtt_raster_cpu&` directly — they
never go through the base class for data access.

**Pure virtual** (all of them; the base implements nothing):

- state: `empty()`, `full_grid_size()`
- metadata: `meta()`, `n_rows()`, `n_cols()`, `cell_size()`
- contiguous storage for core-layer raster writing: `data_ptr()`, `mask_ptr()`
- `clone()`
- loader-to-backend transfer: `set_from_vectors(data&&, mask&&, meta)`

The factory `make_raster_base(const lsdtt_data_loader&)` calls
`make_raster_engine(loader, backend_type::auto_detect)` in
`lsdtt_backend_factory.cpp`. There is no GPU raster engine: the GPU branch of
`make_raster_engine` always throws, and GPU work goes through per-kernel
selective dispatch instead (see "GPU compute kernels" above). Nothing in `cli/`
or the bindings calls either factory today. `auto_detect` picks the GPU branch
whenever `gpu_available()` is true. That function tests `LSDTT_GPU_ENABLED`,
which no build defines, so today it is always false and the factory always
returns an `lsdtt_raster_cpu`. If `LSDTT_GPU_ENABLED` is ever defined before a
GPU raster engine exists, `make_raster_base` will throw instead of falling back
to the CPU. There are no analysis methods on the base:
flow routing, FFT, filtering and smoothing are free functions or engine
objects in `src/cpu/`.

---

## lsdtt_raster_cpu

File: `src/cpu/lsdtt_raster_cpu.hpp/.cpp`

**Storage:** `lsdtt_array2d_f` and `lsdtt_array2d_u8` (each wraps `std::vector`).

**Accessors:**
- `get_data_array()` → `const lsdtt_array2d_f&` — preferred for analysis, zero-copy
- `get_mask_array()` → `const lsdtt_array2d_u8&` — zero-copy
- `get_data()` → `const std::vector<float>&` — zero-copy into array's vector
- `get_mask()` → `const std::vector<uint8_t>&` — zero-copy
- `get_meta()` → `const lsdtt_raster_meta&` (`meta()` is the base-class override)

**Constructors:**
- `lsdtt_raster_cpu(const lsdtt_data_loader&)` — read full raster from GDAL
- `lsdtt_raster_cpu(const lsdtt_raster_tile&)` — from loader tile
- `lsdtt_raster_cpu(lsdtt_array2d_f&&, lsdtt_array2d_u8&&, meta)` — preferred for results
- `lsdtt_raster_cpu(std::vector<float>&&, std::vector<uint8_t>&&, meta)` — backwards compat
- `lsdtt_raster_cpu(lsdtt_raster_result&&)` — adopt a core-layer engine result

---

## lsdtt_param_parser

File: `src/core/lsdtt_param_parser.hpp/.cpp`

Parameter file format: one `key: value` per line (snake_case keys,
case-sensitive, `#` for comments). A value ends at the first whitespace, so a
value containing spaces (typically a list typed as `a, b, c`) is a hard error
that suggests the space-free spelling, rather than silently keeping only `a,`.
Lists are comma-separated with no spaces. Booleans accept `true`/`True`/`TRUE`/
`t`/`T`/`1`; anything else is false.

**File I/O parameters:**
- `read_fname` — full filename **including extension** (e.g. `my_dem.tif`)
- `write_fname` — stem only, no extension (defaults to read stem)
- `read_path` — directory (default: parameter file directory)
- `write_path` — directory (default: read_path). It is **not created**.

`check_paths()` throws if `read_path` or `write_path` is set to a non-empty value
that is not an existing directory. Before this check, a missing `write_path`
made every output fail to open while the run still printed "done" and exited 0.
Unset or blank paths are not checked: they fall back to the parameter file's
directory (`read_path`) or to `read_path` (`write_path`). A blank value is the
empty string, which resolves relative to the current directory.
- `boundary_conditions` — string of characters: `o`=open, `n`=no-flux, `p`=periodic

`get_full_read_path()` returns `read_path + read_fname + read_ext`.

**Other input files** go through `add_read_path_if_needed(name)`: an absolute
path (`/...`, `\\...`, `C:...`) is returned unchanged, and every relative path
(`sites.csv` or `inputs/sites.csv`) is prefixed with `read_path`. Use it for every
`*_fname` input rather than concatenating `get_read_path()` by hand, which breaks
absolute paths. Until 2026-10-14 any name containing a separator was returned
unchanged, so `inputs/sites.csv` was read relative to the working directory.
`get_full_write_path()` returns `write_path + write_fname` (no extension).

**Static helpers:** `stem(filename)`, `extension(filename)` (returns dot + ext),
`parse_args_deferred(argc, argv, tool_name, status, citation)` — handles
`-v`/`-c` (print and exit) and records `-h` variants in the returned
`lsdtt_cli_args` struct so the CLI can render help from its parameter
declarations via `print_terminal_help()` (see the CLI design rule above).
Every CLI uses this instead of a per-file
`parse_args` free function. Each CLI declares
`static const std::string k_tool_name = "lsdtt-foo";`
alongside `k_status` and `k_citation`.

**Command-line contract.** A CLI accepts only `<parameter_file>`,
`<path> <parameter_file>`, or one of the flags `-h` (with an optional
`<parameter_name>`, `html` or `json`), `--help-json`, `-v` and `-c`. Every
analysis option goes in the parameter file. `parse_args_deferred` enforces
this:

| invocation | behaviour |
|---|---|
| no arguments | "No parameter file specified", then the parameter-file guide; exit 1 |
| an unrecognised `-...` argument (e.g. `--input dem.tif`) | "Unrecognised option", then the guide; exit 1 |
| more than two positional arguments | "Too many arguments", then the guide; exit 1 |
| a parameter file that does not exist | the parser's exception says what a parameter file is and points to `-h` |

Without the flag check, `--input dem.tif` would be read as `<path>` =
`--input` and parameter file = `dem.tif`.

**Self-describing binaries.** The release packages ship the binaries and a
README but not these docs, and the person reading the output may be an LLM with
nothing else to go on. So the binary must explain itself: that options live in
a parameter file, how that file is written, and how it is passed. The text
lives in one place, `write_parameter_file_guide()` in
`lsdtt_param_parser.cpp`, and is printed by:

- the no-argument and bad-argument paths above (to stderr, with no example key,
  because the CLI has not declared its parameters yet);
- plain `-h`, before the grouped parameter list (to stdout, with an example key
  taken from the CLI's first switch, e.g. `write_hillshade_raster: true`).

`-h <parameter_name>` adds a line showing how to set that key in the file, and
`-h html` repeats the format summary in its Usage section. `--help-json` is
unchanged, because the QGIS plugin and `tools/generate_cli_reference.py`
consume it.

When editing the guide:

- The end-user README written by `tools/package_{linux,macos}_backend.sh` and
  `tools/package_windows_backend.ps1` ("How the tools work: parameter files")
  says the same things at greater length. Keep the two consistent.
- Every CLI prints this text in its `-h` output, so it must pass every CLI's
  help tests. For example, `rigorous_fetch_raster_no_credential_parameters`
  fails if `lsdtt-fetch-raster -h` contains the word `token`.

**Version reporting (`-v` / `--version`):** prints the package version and
the CLI's status, e.g. `lsdtt-chi-analysis (lsdtt3 0.7.0, stable)`, then the
GDAL/PROJ lines. `--help-json` carries the same two values as `"version"` and
`"status"`, and `-h html` / the end-of-run `_README` files record both. The
package version is injected by CMake as the `LSDTT3_VERSION` compile
definition on `lsdtt3_core` (used only in `lsdtt_param_parser.cpp`, and
exposed as `lsdtt::lsdtt3_package_version()`). A build made exactly on the
release tag `v<PROJECT_VERSION>` reports the plain version; any other build
reports a development string derived from `git describe`, e.g.
`0.7.0-dev, git v0.7.0-14-g04460f4-dirty` — so distributed development
builds self-identify. Captured at configure time (rerun cmake to refresh).
The Python wheel compiles `lsdtt_param_parser.cpp` in its own CMake project,
which sets `LSDTT3_VERSION` to the package version scikit-build-core passes
(`X.Y.Z`, or `X.Y.Z.postN` for a rebuilt wheel); `_core.lsdtt3_version()`
returns it.

### Deprecated, removed and unknown parameter keys

The typed parse functions iterate over a CLI's **defaults maps** and pull
matching keys out of the file, so a key present in the file but absent from
every defaults map is never looked at by anything — silently a no-op. That is
the worst available behaviour when a parameter is renamed or retired: an old
file saying `find_floodplain: true` would run happily having done none of what
it asked for.

`parse_all_parameters` therefore scans the **raw keys present in the file**
(never the defaults maps, which would silently re-accept everything) against
three maps a CLI may populate:

| setter | behaviour | for |
|---|---|---|
| `set_removed_keys` | prints the exact migration message, **exits 1** | keys that TURNED AN ANALYSIS ON, where ignoring one leaves the user believing it ran |
| `set_deprecated_aliases` | copies the value to the new name, warns naming it | a parameter that still exists under a new name |
| `set_ignored_keys` | prints a specific message, run continues | a parameter that moved to another tool, or whose mechanism no longer exists |

Anything in none of the three and in no defaults map is reported as an unknown
key — a warning, never fatal, since a CLI may legitimately read a key by another
route. All three maps start empty, so a CLI that sets none behaves as before.

A removed key **exits rather than throws**: almost no CLI guards
`parse_all_parameters` with try/catch, so throwing would abort through
`std::terminate` and bury the migration message under a crash trace. Exiting
from the parser is what `parse_args_deferred` already does for `-v` and `-c`.

**These maps are a one-release compatibility window and must actually be
emptied**, or they become a permanent second spelling for every parameter.

---

## lsdtt_data_loader

File: `src/core/lsdtt_data_loader.hpp/.cpp`

`lsdtt_data_loader` handles both raster and vector I/O through GDAL/OGR.

`open_raster` refuses geotransforms the meta cannot represent (the meta keeps
only an origin and one cell size): rotated rasters and south-up rasters always,
and non-square pixels unless the CRS is geographic. Geographic rasters must still
open, because `lsdtt-raster-preprocessing` reads them to reproject them; analysis
tools refuse geographic rasters themselves. `band_count()` reports the bands of
the open raster.

Its header includes no GDAL header: `GDALDataset` is forward declared, the
dataset handles sit behind a deleter defined in the `.cpp`, and the OGR layer
argument of `write_record_to_layer` is an opaque `void*`. See **The GDAL
boundary** below.

### Raster output

**`raster_output_format` enum:**
- `cog` — Cloud Optimised GeoTIFF (default; falls back to gtiff for small rasters)
- `gtiff` — plain compressed GeoTIFF
- `same_as_input` — use the GDAL driver of the input raster

**Key raster methods:**
- `write_full(data, meta, path, format=cog)`
- `finalise_output(path)`
- `create_output_raster(path, meta, format=cog, n_bands=1)`

### Vector output

Various file formats are available for vector output but preference is for either csv
or FlatGeobuf. FlatGeobuf should be the default output since it is more storage efficient
than csv, includes georeferencing, and does not have the same file writing issues
as gpkg. GPKG struggles wien mapping across drives, which coould be a common problem if 
users are on a WSL2 implementation of linux writing to `/mnt/c`.  

Supported vector outputs are:
- GeoPackage (`.gpkg`)
- FlatGeobuf (`.fgb`)
- CSV

lsdtt3 holds vector layers in memory using an `lsdtt_vector_layer` structure, which contains:
- an `lsdtt_vector_schema` structure to hold attribute field names and types
- a number of features, each stored as an `lsdtt_vector_record`. Records contain an `lsdtt_geom` geometry (which can be point, line or polygon) and a vector of attribute table values for that feature.

The `loader` object in `lsdtt_data_loader.hpp/cpp` manages vector writes, using:
- `begin_vector_write(...)` to start a writing session, takes parameters for paths, layer names, the schema, the crs, the format (gpkg/flatgeobuf) and whether or not to write a spatial index. Including the spatial index makes loading data and spatial queries faster in a GIS, but to calculate the (optional) spatial index, the entire vector layer needs to sit in memory (no feature-by-feature writing) which results in large overheads when writing large vector datasets.
- `write_point/line/polygon(..)` write individual features, which are streamed if no spatial index required and buffered in memory if spatial index required.
- `finish_vector_write()` closes the writing session, and calculates the spatial index if required.

GeoPackage and FlatGeobuf are written through GDAL/OGR; no direct SQLite coding
is required in `lsdtt3`.

Fields in FlatGeobuf should be typed and not strings.


---


## Vector I/O

Vector I/O is handled through `lsdtt_data_loader` and shared helpers in
`src/core/lsdtt_vector_utils.hpp/.cpp`.

### Preferred output formats

See above, FaltGeobuf is the preferred format. 

Supported vector formats are:

- **GeoPackage (`.gpkg`)**  
  Default binary single-file GIS output.

- **FlatGeobuf (`.fgb`)**  
  Lightweight binary single-file GIS output.

- **CSV**  
  Mainly for backward-compatible tabular export and plotting workflows.

### Attribute typing

Vector fields should preserve their natural types:
- integer attributes → integer fields
- floating-point attributes → real fields
- text attributes → string fields

Do not coerce all vector attributes to strings unless there is a compelling
backward-compatibility reason.

This is important for fields such as:
- stream order
- junction id
- receiver id
- row / col / flat index
- elevation
- flow distance

### CSV spatial convention

If point data are written to CSV for interoperability with older scripts, the
CSV should contain:
- `latitude`
- `longitude`

in **EPSG:4326**, using GDAL/OGR coordinate transformation from the source CRS.

### Shared helper split

- **`src/core/lsdtt_types.hpp`**  
  all lightweight shared primitive types (spatial, raster geometry, flow, junction)

- **`src/utils/lsdtt_raster_geometry.hpp`**  
  header-only grid geometry helpers such as row/col to cell-centre x/y;
  primitive types (lsdtt_xy, lsdtt_extent, lsdtt_rowcol_window) come from lsdtt_types.hpp

- **`src/core/lsdtt_vector_utils.hpp/.cpp`**  
  GDAL/OGR-backed vector/projection helpers such as reusable coordinate
  transformations

---


## Flow routing — how to navigate the files

There are several flow-routing files. This table explains which one you need.

| You are... | Include this | Notes |
|---|---|---|
| Writing analysis code (chi, steepness, drainage area, etc.) | `lsdtt_flowinfo_cpu.hpp` | Only file you need. Construct `lsdtt_flowinfo`, call its methods directly. |
| Writing a CLI that builds a flow object | `lsdtt_flow_base.hpp` | Call `make_flow_base()`. Returns `unique_ptr<lsdtt_flow_base>`. GPU-compatible. |
| Writing D-inf / MFD distributed-flow analysis | `lsdtt_flowgraph_cpu.hpp` | Fully implemented (CSR layout, Kahn BFS). |
| Adding a GPU flow routing backend | `lsdtt_flow_base.hpp` + new `lsdtt_flowinfo_gpu.hpp` | Implement `lsdtt_flow_base`, add branch in `src/cpu/lsdtt_flow_base.cpp`. |
| Understanding cell-index types | `lsdtt_types.hpp` | Flow-routing cell indices use `std::int64_t` directly throughout; there is no `flat_index_type` typedef. |

**Analysis code** never needs to know about `lsdtt_flow_base`
or `lsdtt_flow_base.cpp`. Those are infrastructure.

---

## lsdtt_types

File: `src/core/lsdtt_types.hpp` (header-only, no `.cpp`)

**Single home for all lightweight primitive types shared across lsdtt3.**

Dependencies: only `<cstdint>`, `<cstddef>`. No GDAL, no OGR, no other lsdtt3 header.

**Rule:** Do not create a new `*_types.hpp` file. If you need a new lightweight
shared type, add it here. A separate types file is only justified when including
the type everywhere would create an unwanted compile-time dependency — that
situation should be extremely rare in this codebase.

**Include chain (acyclic):**
```
lsdtt_types.hpp                          (no lsdtt3 includes — just <cstdint> etc.)
       ↑
lsdtt_flow_base.hpp                      (includes lsdtt_types.hpp; owns lsdtt_flowinfo_options)
       ↑
lsdtt_flowinfo_cpu.hpp                   (includes lsdtt_flow_base.hpp)
       ↑
lsdtt_flow_base.cpp                      (in lsdtt3_cpu — includes lsdtt_flowinfo_cpu.hpp, no cycle)

lsdtt_types.hpp
       ↑
lsdtt_junction_network_base.hpp          (includes lsdtt_types.hpp; owns lsdtt_junction_network_options)
       ↑
lsdtt_junction_network_cpu.hpp           (includes lsdtt_junction_network_base.hpp)

lsdtt_types.hpp
       ↑
lsdtt_raster_geometry.hpp               (includes lsdtt_types.hpp + lsdtt_raster_meta.hpp)
```

**Contents:**

*Spatial primitives*
- `lsdtt_pixel_rc` — integer row/col location; negative = unset/invalid
- `lsdtt_xy` — double x/y projected coordinate pair

*Raster geometry*
- `lsdtt_extent` — axis-aligned bounding box (lower-left origin)
- `lsdtt_rowcol_window` — half-open `[row_begin, row_end)` / `[col_begin, col_end)` window

*Flow-routing types*
- `k_nodata_receiver = -1`, `k_outlet_receiver = -2` — `std::int64_t` sentinels
- `k_no_direction = 255` — `std::uint8_t` sentinel
- `lsdtt_flat_index_slice` — zero-copy view into s_vector; `value_type = std::int64_t`

*Junction-network types*
- `k_nodata_junction = -1`, `k_outlet_junction = -2` — `std::int32_t` sentinels
- `lsdtt_junction_index_slice` — zero-copy view into junction s_vector; `value_type = std::int32_t`

**`lsdtt_flowinfo_options`** lives in `src/core/lsdtt_flow_base.hpp` — it is
domain-specific to flow routing and would only clutter the shared types file.

**`lsdtt_flowinfo_options` fields:**
```
boundary_conditions       = "oooo"   # north east south west
interior_nodata_is_outlet = true
require_acyclic_graph     = true
build_topological_order   = true
build_s_vector            = true
allow_pits_as_outlets     = false    # true = treat transient pits as outlets (used by the LEM)
n_threads                 = 0
```
Boundary codes: `o`=open, `n`=no-flux, `p`=periodic (symmetric pairs),
`b`=forced outlet on that entire edge.

**`lsdtt_junction_network_options`** lives in `src/core/lsdtt_junction_network_base.hpp`.

---

## lsdtt_flow_base

File: `src/core/lsdtt_flow_base.hpp` (header-only — factory impl in src/cpu/lsdtt_flow_base.cpp)

Abstract base class for D8 flow routing backends.

**Factory function `make_flow_base()`**:
Declared in `lsdtt_flow_base.hpp`, implemented in `lsdtt_flow_base.cpp` (lsdtt3_cpu).
Currently returns `std::make_unique<lsdtt_flowinfo_cpu>`.

**Adding GPU flow routing in future:**
1. Write `src/gpu/lsdtt_flowinfo_gpu.hpp/.cpp` implementing `lsdtt_flow_base`
2. Add `gpu_available()` branch in `lsdtt_flow_base.cpp`
3. Add `lsdtt_flowinfo_gpu.cpp` to `lsdtt3_gpu` in CMakeLists
4. All existing CLIs work unchanged

---

## lsdtt_flowinfo_cpu

File: `src/cpu/lsdtt_flowinfo_cpu.hpp/.cpp`

CPU-only D8 flow routing. Inherits `lsdtt_flow_base`. All storage is
host-side (`lsdtt_array2d`, `std::vector`).
Requires a pit-free, flat-resolved DEM — use `lsdtt_raster_flow_preprocessing`
first.

**Direction code convention (receiver_dir_):**
```
7  0  1
6  x  2
5  4  3
```

**Build pipeline (in order):**
1. `validate_inputs()`
2. `parse_boundary_conditions()`
3. `initialise_storage()`
4. `build_receiver_graph_d8()` — parallelisable with OMP collapse(2)
5. `build_donor_count_d8()` — serial
6. `build_outlets()` — serial
7. `build_topological_order()` — Kahn's BFS, serial
8. `build_s_vector()` — iterative DFS, serial

**s_vector contract:**
- `s_vector_[s_vector_start_[i]] == i` — node is always first in its own slice
- `upslope_nodes(flat)` — returns zero-copy `lsdtt_flat_index_slice`
- `is_upslope_of(candidate, target)` — O(1) interval membership test

---

## lsdtt_flowgraph_cpu

Files: `src/cpu/lsdtt_flowgraph_cpu.hpp/.cpp`

DAG-based CPU flow routing for D-infinity and MFD. No s_vector (upslope sets
are fractional and overlapping). Use for distributed-flow analyses: drainage
area (D-inf/MFD variant), TWI, sediment flux, specific catchment area on
hillslopes. Also provides the routing snapshot consumed by `lsdtt_hydraulic_graph_cpu`.

**Routing methods:**
- `d_infinity` — Tarboton (1997). At most two receivers per cell (one cardinal,
  one diagonal). Proportions from within-facet angle `r`. Uses `atan2(s2, s1)`
  for numerical stability; no artificial slope perturbation.
- `mfd` — Quinn et al. (1991) / Freeman (1991). Up to eight receivers,
  weighted by `slope^exponent`. `mfd_exponent` defaults to 1.0 (Quinn);
  Freeman recommends 1.1 for natural terrain.

**Data model (CSR):**
- `receiver_offsets_`, `receiver_indices_`, `receiver_weights_`
- Cell `i`'s receivers: `[receiver_offsets_[i], receiver_offsets_[i+1])`
- Weights sum to 1.0 per valid cell that has at least one receiver.

**Topological order:** Kahn's BFS with Anand et al. (2020) ready-queue
extension for multi-receiver DAGs. Sources first, outlets last.

**Key methods:**
- `accumulate(cell_weights = {})` → `std::vector<double>` — forward push sweep
- `accumulate_raster(cell_weights = {})` → `lsdtt_raster_result`
- `dinf_slopes()` → `lsdtt_raster_result` — steepest-descent slope magnitude
- `dinf_flow_direction()` → `lsdtt_raster_result` — compass bearing (0=N, CW)
- `n_outlets()` — count of sink cells in topological order

**Options (`lsdtt_flowgraph_options`):**
```
boundary_conditions       = "oooo"   # north east south west
interior_nodata_is_outlet = true
method                    = mfd      # d_infinity or mfd
mfd_exponent              = 1.0f
n_threads                 = 0        # 0 = OMP default
```

**Boundary condition codes:**
- `o` / `n` — open / no-flux: neighbours outside bounds are skipped; edge
  cells become sinks if no in-bounds downslope neighbours exist.
- `p` — periodic: wraps to opposite edge (north↔south, east↔west must match).
- `b` — forced outlet: entire edge is forced to be sinks.

---


## Junction-network architecture

The junction-network layer is the compact channel-graph analogue of the pixel
flow-routing graph.

### Files

- `src/core/lsdtt_types.hpp` — shared primitive types (sentinels, index slices)
- `src/core/lsdtt_junction_network_base.hpp` — abstract interface + `lsdtt_junction_network_options`
- `src/cpu/lsdtt_junction_network_cpu.hpp/.cpp`

### Type policy

- **pixel / raster flat indices** use `std::int64_t`
- **junction indices** use `std::int32_t`

This split is deliberate:
- full rasters may exceed the safe range of 32-bit indexing
- junction counts are much smaller than pixel counts

### Construction model

`lsdtt_junction_network_cpu` is built on top of an existing flow-routing object.

The flow object is not owned by the junction network and must outlive it.

### Persistent storage

The persistent raster-scale storage is intentionally small:
- channel mask
- junction-id raster

The compact graph stores:
- junction flat index
- receiver junction
- donor count
- stream order
- outlets
- topological order
- junction `s_vector`

Temporary per-pixel products, such as pixel stream order used during
construction or export, should remain local and not be stored permanently unless
there is a clear reason.

---


## CPU analysis free functions (src/cpu/)

All take `const lsdtt_raster_cpu&` directly — not `lsdtt_raster_base`.
All return `lsdtt_raster_cpu`. Input is never modified.

**`lsdtt_raster_algebra` (`src/cpu/lsdtt_raster_algebra.hpp/.cpp`):**

Lives in `namespace lsdtt::algebra`. All operations are element-wise; output inherits
metadata from the left-hand operand. A cell is nodata if either input cell is nodata.
Compatibility (same n_rows, n_cols, cell_size, x_min, y_min, no_data_value) is
enforced by `assert_compatible()` — throws `std::invalid_argument` on mismatch.
Inner loops are OpenMP-parallelised.

- `assert_compatible(a, b)` — throws if rasters are not element-wise compatible
- `raster_add(a, b)` — result = a + b
- `raster_subtract(a, b)` — result = a − b
- `raster_multiply(a, b)` — result = a × b
- `raster_divide(a, b)` — result = a / b; additionally nodata where b == 0
- `raster_log(a)` — natural log; additionally nodata where a ≤ 0
- `raster_log10(a)` — base-10 log; additionally nodata where a ≤ 0
- `raster_negate(a)` — unary minus: result = −a
- `raster_scalar_subtract(c, a)` — scalar minus raster: result = c − a  (e.g. DEM inversion: `raster_scalar_subtract(max_val, dem)`)

**`lsdtt_raster_surface_metrics_cpu`:**
- `compute_hillshade(raster, azimuth=315, altitude=45, z_factor=1.0)`
- `moving_window_mean(raster, radius_metres, shape=circular)`
- `moving_window_local_relief(raster, radius_metres, shape=circular)` — max minus min in window
- `moving_window_stddev(raster, radius_metres, shape=circular)` — population stddev
- `moving_window_median(raster, radius_metres, shape=circular)`
- `moving_window_mad(raster, radius_metres, shape=circular)` — median absolute deviation
- `moving_window_tpi(raster, radius_metres, shape=circular)` — Topographic Position Index (Weiss 2001)
- `calculate_madk2(...)` — **moved to `lsdtt_surface_texture_cpu`**
- `direct_polyfit(raster, lambda_metres=100.0, n_threads=0)` → `poly_fit_coeffs` — **default and preferred polyfit method**; separable two-pass rolling-moment, O(n) per pixel, faster than wavelet and better at resolving fine-scale features
- `wavelet_quadratic_fit(raster, lambda=100.0, n_threads=0)` → `poly_fit_coeffs` — FFT-based alternative; available but not the default
- `compute_polyfit_metrics(polycoeff, request)` → `polyfit_metric_result` — batch computation of curvature and derivative rasters from polyfit coefficients (works with output of either method above)
- `window_shape` enum: `square`, `circular`
- `polyfit_metric` enum: slope, aspect, curvatures (mean, Gaussian, max, min, profile, plan, tangential, etc.), rotor

### Curvature sign convention

lsdtt3 uses the **strict mathematical sign convention** for all curvature metrics:

- **Concave-up terrain (channels, hollows, bowl-shaped) → positive curvature**
- **Convex-up terrain (ridges, noses, dome-shaped) → negative curvature**

This is the convention that follows directly from the second derivatives of the fitted surface:
a valley bottom has positive ∂²z/∂x² (concave up), so curvature there is positive.

**This is NOT the Zevenbergen & Thorne (1987) convention**, which defines plan and profile
curvature with an extra leading negative sign, making channels negative and ridges positive.
The Z&T convention is common in older GIS software (e.g., ArcGIS) and many geomorphology
papers; lsdtt3 deliberately does not follow it.

The strict mathematical convention means:
- `plan_curvature > 0` at a channel (convergent flow, concave in plan)
- `plan_curvature < 0` at a ridge (divergent flow, convex in plan)
- `profile_curvature > 0` at a concave hillslope profile (decelerating flow)
- `profile_curvature < 0` at a convex hillslope profile (accelerating flow)
- `mean_curvature > 0` at net concave terrain
- `laplacian > 0` at net concave terrain

`compute_polyfit_metrics` prints a one-line reminder of this convention to stdout each
time it is called so that the sign is never ambiguous in logged output.

**`lsdtt_raster_surface_fitting_cpu`:**

Plane and polynomial surface fitting to masked or labelled raster regions.
Use for piedmont / fan fitting, terrace treads, floodplain surfaces.
For spectral (long-wavelength) surface fitting see `lsdtt_raster_spectral_cpu`.

- `fit_plane_to_mask(dem, mask)` → `lsdtt_plane_fit_result` — OLS plane, z = a*x + b*y + c
- `fit_plane_to_mask_ransac(dem, mask, inlier_tolerance, max_iter=1000, seed=1)` → `lsdtt_plane_fit_result` — robust RANSAC fit; tolerates incised channels and outliers
- `evaluate_fitted_plane(meta, coeffs, n_threads=0)` → `lsdtt_raster_cpu`
- `residual_from_fitted_plane(dem, coeffs, n_threads=0)` → `lsdtt_raster_cpu`
- `fit_planes_to_labels(dem, label_raster)` → `std::map<int, lsdtt_plane_fit_result>` — one plane per label
- `fit_poly2_to_mask(dem, mask)` → `lsdtt_poly2_fit_result` — OLS degree-2 polynomial, z = a*x² + b*y² + c*x*y + d*x + e*y + f
- `evaluate_fitted_poly2(meta, coeffs, n_threads=0)` → `lsdtt_raster_cpu`
- `residual_from_fitted_poly2(dem, coeffs, n_threads=0)` → `lsdtt_raster_cpu`
- `poly2_max_principal_curvature(coeffs)` → `double` — max Hessian eigenvalue for curvature screening

**`lsdtt_raster_flow_preprocessing_cpu`:**
- `mask_to_nodata_using_threshold(raster, threshold, is_below_nodata)`
- `fill_pits(raster, min_slope=0.0001f, [boundary_conditions,]
  [flat_resolution,] n_threads=0, queue_backend=binary_heap)` — dispatches on
  `min_slope`: > 0 → Priority-Flood with the gradient imposed inline (Barnes
  et al. 2014); <= 0 → grayscale reconstruction (Vincent 1993), same fill
  levels, faster. The `flat_resolution_choice` overloads run the Barnes 2014b
  resolver after a flat fill.
- `fill_pits_with_flat_resolution(...)` → `flat_resolution_result` — the same,
  returning the resolver's report (rounds, cells raised, final-pass counts).
- `fill_pits_priority_flood(raster, min_slope, n_threads, queue_backend)` and
  `fill_pits_reconstruction(raster, ...)` — the two backends, called directly.
- `priority_flood_queue_backend` enum: `binary_heap` (default), `radix_hybrid`.
  Same fill levels; at `min_slope > 0` the imposed gradients on flats can
  differ, because ties are broken in a different order.
- `fill_pits_with_boundary(...)` — serial only
- `breach_depressions(raster, n_threads=0)` — Lindsay 2016 complete breaching

**`lsdtt_raster_breach_least_cost_cpu`** (Lindsay 2020 least-cost breaching):
- `breach_depressions_least_cost(raster, const breach_least_cost_options&,
  breach_least_cost_diagnostics* = nullptr)` — the form every caller uses.
- `breach_least_cost_options`: `max_dist` (path length in cells; default
  1000), `max_cost` (inf), `minimize_dist` (false), `flat_increment` (<= 0 →
  auto, 2 × ULP_f32(max(|z|max, 1))), `fill_remainder` (true),
  `fill_min_slope` (1e-4), `fill_flat_resolution` (`inline_gradient`),
  `edge_policy` (`outlet`).
- `breach_lc_edge_policy` enum: `outlet` (an off-grid neighbour ends the search
  like a nodata neighbour; the user default) and `whitebox` (off-grid
  neighbours skipped; the bitwise oracle mode). `default_breach_lc_max_dist()`
  gives 1000 for `outlet` and 100 for `whitebox`; callers resolve an unset
  max_dist through it, so an explicit value always wins.
- `parse_breach_lc_edge_policy` / `breach_lc_edge_policy_name`,
  `validate_breach_least_cost_options` (throws `std::invalid_argument` naming
  the parameter; the CLIs call it before reading any raster).
- `breach_least_cost_diagnostics`: pits detected / searched / skipped / solved
  / unsolved, cells without descent before the remainder fill, the increment
  used.
- A legacy positional overload remains (`max_dist=1000, max_cost, min_dist=false,
  flat_increment=0, fill_remainder=true, n_threads`), forwarding to the options
  form; no CLI or binding uses it.
- Numerics: a float64 working copy, rounded once to float32; the search heap
  is a hand port of Rust's `BinaryHeap`, so ties pop in Whitebox's order; the
  TU is compiled with `-ffp-contract=off` (see `CMakeLists.txt`). With
  `edge_policy = whitebox`, an explicit increment and `fill_remainder = false`
  the output is float32(Whitebox 2.0.6's output), bit for bit. Serial.
- Exposed in Python as `breach_least_cost` and `least_cost_carve_depth`
  (`max_distance=None` → the edge-policy default, `prefer_shortest_path=False`,
  `edge_policy="outlet"`, `fill_min_slope=1e-4`).

**`lsdtt_raster_flat_resolution_cpu`** (Barnes et al. 2014b):
- `flat_resolution_choice` enum: `inline_gradient` (default; no separate
  pass), `barnes2014b_gradient`, `barnes2014b_float_step`;
  `parse_flat_resolution_choice` / `flat_resolution_choice_name`.
- `compute_barnes2014b_flat_mask(filled, ...)` — the towards-low /
  away-from-high mask.
- `resolve_flats_barnes2014b(filled, flat_resolution_options)` →
  `flat_resolution_result` — up to four refill-and-resolve rounds and a final
  strict-rise Priority-Flood pass; never lowers a cell; output independent of
  thread count.
- Algorithm and guarantees: `docs/design_and_algorithms/flow_routing.md`.

**`lsdtt_raster_nodata_filling_cpu`:**
- `identify_interior_nodata(raster)` → `lsdtt_array2d<uint8_t>`
- `fill_interior_nodata(raster, method=laplacian, max_iter=10000, tolerance=1e-6, n_threads=0)`
- `clip_to_maximum_rectangle(raster)`
- `nodata_fill_method` enum: `iterative_weighted_average`, `laplacian`

**`lsdtt_raster_spectral_cpu`:**
- `wavelet_quadratic_fit(raster, lambda=100.0, n_threads=0)`
  — FFTW double API, Gaussian kernel; **not the default** — prefer `direct_polyfit` (in surface_metrics)
  — Returns `poly_fit_coeffs{a, b, c, d, e, f}` — six `lsdtt_raster_cpu` members
- `fftw2D_masked_surface_recovery(dem, mask, lambda_cutoff_metres, window=hann, n_threads=0)`
  — Recovers smooth long-wavelength surface from pixels where mask > 0.5 only
  — Detrends using all valid pixels, zeros excluded pixels, lowpass filters, restores trend
  — Output valid everywhere the DEM is valid (interpolates under the masked-out gap)
  — Intended use: piedmont surface recovery for mountain front detection

**`lsdtt_streamlines_cpu`:**

Kinematic streamline tracing (steepest-descent flow paths) on raw topography.
Implements the Lea (1992) algorithm with Tarboton (1997) improvements.

Key types:
- `streamline_outcome` enum: `reached_sink`, `hit_edge`, `loop_detected`, `no_aspect`
- `lsdtt_streamline_node{x, y, z, distance_m}` — one sub-pixel path point
- `lsdtt_streamline{source_row, source_col, sink_row, sink_col, length_m, relief_m, outcome, nodes}`
- `lsdtt_streamline_options{thinning, max_steps, save_nodes, n_threads}`
- `lsdtt_streamline_result{streamlines, n_reached_sink, n_hit_edge, n_loop_detected, n_no_aspect}`

`trace_streamlines(dem, source_mask, sink_mask, opts)`:
- Precomputes aspect (steepest-descent bearing) from DEM central differences
- Traces from every source pixel (source_mask > 0.5) until hitting a sink pixel, raster edge, or loop
- OMP-parallelised over source pixels; loop detection is per-thread (unordered_map)
- `thinning` controls trace density (every Nth source pixel in row-major order)
- `save_nodes` enables sub-pixel path storage for visualisation

Intended as the foundation for `lsdtt_hilltop_flow_routing_cpu` (ridge statistics, E*, R*).

**`lsdtt_mountain_front_cpu`:**

Three-stage pipeline for identifying mountain fronts from a DEM.

```
Stage 1 — low-ground candidate mask (multi-scale relief)
Stage 2 — piedmont surface fitting (plane OLS/RANSAC, poly2 OLS, or spectral lowpass)
Stage 3 — residual thresholding → mountain mask
```

Key types:
- `relief_scale{radius_metres, max_relief_metres}` — one moving-window scale
- `mountain_front_fit_method` enum: `plane_ols`, `plane_ransac`, `poly2_ols`, `spectral_lowpass`
- `mountain_front_request` — all parameters; see header for full field documentation
- `mountain_front_result{low_ground_mask, fitted_surface, residual, mountain_mask, relief_rasters}`

`relief_rasters` is populated only when `request.return_diagnostics == true`; contains one raster per entry in `relief_scales`, sorted by radius (smallest first).

`residual_percentile_threshold` (0.0–1.0): pixels above this percentile of the valid residual distribution are labelled mountain.  E.g. 0.80 → top 20 % → mountain.

Stage 1 mask logic:
- `require_all_scales_pass = false` (default): pixel passes if largest-scale relief < threshold
- `require_all_scales_pass = true`: pixel must pass every scale's threshold

Ground-truth validation: use `lsdtt_raster_synthetic_surfaces_cpu` test surfaces
(`generate_linear_range_surface`, `generate_halfspace_range_surface`,
`generate_circular_range_surface`) whose `front_mask` raster is the oracle.

**`lsdtt_raster_synthetic_surfaces_cpu`:**
- Takes `lsdtt_raster_meta` or raw dimensions (n_rows, n_cols, cell_size)
- Basic generators: uniform random, red noise, blue noise, diamond-square, periodic wave,
  global quadratic, random quadratic
- Mountain range generators (return `lsdtt_mountain_range_surface{dem, front_mask}`):
  - `generate_linear_range_surface` — bilateral piedmont, ridge in DEM interior
  - `generate_halfspace_range_surface` — ridge at DEM edge (Himalaya half-space analogy)
  - `generate_circular_range_surface` — elliptical dome in DEM interior
  - `generate_arc_range_surface` — parabolic curved front at one edge
- `lsdtt_range_edge` enum: `north`, `south`, `east`, `west`

---


## Utility-layer design (`src/utils/`)

The `src/utils/` layer is for generic mathematical, numerical, and lightweight
helper code that is reusable across the project and is not tied to GDAL/OGR,
file I/O, or a specific backend implementation.

Typical contents include:
- lightweight containers
- metadata-only raster geometry helpers
- descriptive statistics
- histogram / PDF / CDF tools
- binning utilities
- nonparametric statistical tests
- circular-data comparison tools

### Valid-data contract

Utility functions in `src/utils/` assume that the values passed to them are
already valid.

This means `utils` code should **not** perform:
- nodata filtering
- raster flattening
- moving-window extraction
- basin extraction
- backend-specific mask handling

Those tasks belong in:
- `src/core/` when they are shared infrastructure or GDAL/OGR-aware helpers
- `src/cpu/` for concrete CPU-side extraction and traversal
- `src/gpu/` for concrete GPU-side extraction and traversal

In practice, raster, moving-window, and basin code should extract valid values
first and then pass the resulting flat sample into `src/utils/` routines.




## lsdtt-raster-preprocessing CLI

File: `cli/lsdtt-raster-preprocessing.cpp`  |  `k_status`: stable

**Preprocessing pipeline (in order):**
1. **Projection** — detect and reproject geographic / Web Mercator rasters
2. **Threshold masking** — remove pixels outside elevation range
3. **Interior nodata filling** — fill enclosed voids
4. **Clip to rectangle** — largest valid rectangle
5. **Depression breaching** — Lindsay 2016 complete breaching (default) or
   Lindsay 2020 least-cost (`breach_algorithm_choice: least_cost`, with its own
   remainder fill)
6. **Pit filling** — Barnes 2014 Priority-Flood (default), Vincent 1993
   reconstruction, or Fill-Spill-Merge; flats by the inline gradient or Barnes
   2014b (`fill_flat_resolution_choice`)
7. **Hillshade** — Zevenbergen & Thorne 1987

Steps 2–7 always operate on the (possibly reprojected) raster. If reprojection
occurs the new file is reopened before step 2.

**Threshold masking semantics:**
`min_elevation` masks pixels **at or below** the threshold (i.e. `value <= min_elev`).
Default is `0.0`, which catches ocean pixels that raster providers commonly set to exactly 0.
`maximum_elevation` masks pixels above the threshold (`value > max_elev`); default `30000.0`.

**Projection parameters:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enforce_projection` | bool | true | Reproject geographic/WebMercator input automatically. If false, exit with error. |
| `resampling_method` | string | cubic | Any GDAL warp resampling string. |
| `grid_spacing` | float | -1 | Target grid spacing (m). ≤ 0 = auto-snap from input pixel size. |
| `allow_lambert_projection` | bool | false | Use LAEA instead of UTM when raster spans ≥ 3 UTM zones. |

**All other parameters:**
```
read_fname, write_fname, read_path, write_path
write_output_format: cog / gtiff / same_as_input
write_hillshade_raster, write_fill_raster, write_fill_diff_raster
write_breached_raster, write_carve_diff_raster
write_filled_nodata_raster, write_clipped_rectangle_raster
breach_before_fill, fill_remainder_after_breach
fill_algorithm_choice: barnes2014 / vincent1993 / barnesetal2021
priority_flood_queue: binary / radix
fill_flat_resolution_choice: inline_gradient / barnes2014b_gradient / barnes2014b_float_step
breach_algorithm_choice: lindsay2016 / least_cost
nodata_fill_method_choice: laplacian / iterative_weighted_average
hillshade_azimuth, hillshade_altitude, hillshade_z_factor
fill_min_slope, min_elevation, maximum_elevation
breach_lc_max_dist (default 1000 outlet / 100 whitebox), breach_lc_max_cost,
breach_lc_minimize_dist (false), breach_lc_edge_policy (outlet / whitebox),
breach_lc_flat_increment (0 = auto)
n_threads, nodata_fill_max_iter
```

**Diff raster outputs** (use `lsdtt::algebra::raster_subtract` from `lsdtt_raster_algebra.hpp`):

| Parameter | Output suffix | Formula | Interpretation |
|---|---|---|---|
| `write_fill_diff_raster` | `_fill_diff.tif` | filled − original | Positive = surface raised by filling |
| `write_carve_diff_raster` | `_carve_diff.tif` | breached − original | Negative = surface lowered by breaching |

**Help / parameter printing:**
- `-h` flag prints usage to stdout and exits immediately (before any processing).
- `print_parameters()` is called at the end of every successful run to echo resolved parameters.
- `print_help()` on `lsdtt_param_parser` writes CSV and HTML help files, but this CLI does not call it during `-h` handling or normal runs.

---

## lsdtt-drainage-analysis CLI

File: `cli/lsdtt-drainage-analysis.cpp`

Computes drainage-related rasters. Uses `make_flow_base()` for D8 or
`lsdtt_flowgraph_cpu` for D-infinity / MFD.

On a raw DEM (`raster_is_filled: false`) it runs `fill_pits` first and prints a
`Compute time:` line for the fill, the D8 build and each accumulation; this is
the default route from a raw DEM to a routed DEM described in
`docs/design_and_algorithms/flow_routing.md`.

**Key parameters:**
```
read_fname, write_fname, read_path, write_path
mfd_exponent: 1.0
boundary_conditions: oooo
raster_is_filled (false: fill_pits on the raw DEM first)
fill_min_slope: 0.0001, fill_flat_resolution_choice: inline_gradient
build_flow_s_vector: true
write_d8_drainage_area_raster, write_dinf_drainage_area_raster,
write_mfd_drainage_area_raster, write_twi_raster, write_sti_raster,
write_connectivity_raster, write_height_above_channel_raster
write_output_format: cog
n_threads: 0
```

**Flow direction codes (D8):**
`0`=N `1`=NE `2`=E `3`=SE `4`=S `5`=SW `6`=W `7`=NW; outlets = `-2`

---

## lsdtt_stat_csv

File: `src/core/lsdtt_stat_csv.hpp/.cpp`

Shared CSV writers for statistical outputs. Used by multiple CLIs to avoid
repeating boilerplate. All functions take a `lsdtt_data_loader&` for I/O.

- `write_qq_csv(loader, path, qq)` — QQ-plot table (standard normal variate, data quantile, normal reference line)
- `write_smoothed_cdf_csv(loader, path, scdf, fits)` — empirical CDF + best-fitting parametric CDF + divergence column
- `write_dist_fits_csv(loader, path, fits)` — one row per fitted distribution family with KS statistic and validity flag

---

## lsdtt_hydraulic_base / lsdtt_hydraulic_graph_cpu

### lsdtt_hydraulic_base

File: `src/core/lsdtt_hydraulic_base.hpp` (header-only)

Narrow abstract interface for graph-based hydraulic routing engines, analogous
to `lsdtt_flow_base`. Conceptually inspired by Gailleton et al. (2024).

**`lsdtt_hydraulic_options` fields:**
```
dt                              = 1.0e-3   # timestep (s)
min_slope                       = 1.0e-4
receiver_width                  = -1.0     # <= 0 = auto (cell_size)
friction_depth_exponent         = 5/3      # q / width = coefficient * h^a * S^b
friction_slope_exponent         = 0.5      # (Manning: coefficient = 1/n)
adaptive_dt_enabled             = false
adaptive_dt_min / adaptive_dt_max                 = 1.0e-4 / 1.0
adaptive_target_median_abs_dhw  = 1.0e-3
adaptive_dt_growth_factor / adaptive_dt_shrink_factor = 1.5 / 0.5
print_startup_message           = true
```

Rainfall and roughness are per-cell fields, not options. The resistance law is
written in a generic form so that Chézy or other laws can be added without
changing the hydraulic balance loop.

**Key virtual methods:**
- `run_one_step()`, `run_n_steps(n)`, `current_time()`, `step_count()`
- forcing: `set_rainfall_rate_field(m/s per cell)`,
  `set_friction_coefficient_field(...)`, `add_point_inflow(...)`,
  `clear_point_inflows()`
- state: `set_bed_elevation(...)`, `set_water_depth(...)`, `reset_water_depth(...)`,
  `total_water_volume()`
- outputs: `get_bed_raster()`, `get_water_depth_raster()`,
  `get_water_surface_raster()`, `get_qin_raster()`, `get_qout_raster()`,
  `get_dhw_raster()`
- `refresh_routing()` and `routing_is_dynamic()`

### lsdtt_hydraulic_graph_cpu

File: `src/cpu/lsdtt_hydraulic_graph_cpu.hpp/.cpp`

CPU implementation of `lsdtt_hydraulic_base`. Prints a startup message crediting
Gailleton et al. (2024).

**Run modes** (`set_run_mode()`):
```cpp
enum class lsdtt_hydraulic_run_mode {
    transient_fixed_graph,      // default
    stationary_fixed_graph,
    transient_dynamic_graph,
    stationary_dynamic_graph
};
```

All four are implemented. In the fixed-graph modes the routing graph is a
snapshot copied from an `lsdtt_flowgraph_cpu`. In the dynamic-graph modes the
routing is rebuilt from the current water surface. `refresh_routing()` does
whichever applies: it rebuilds from the water surface when a hydraulic domain is
set, and otherwise re-copies the source flowgraph. Constructing the engine on a
hydraulic domain selects `stationary_dynamic_graph`.

---

## lsdtt_surface_texture_cpu

File: `src/cpu/lsdtt_surface_texture_cpu.hpp/.cpp`

Surface texture and roughness metrics. MADk2 was moved here from
`lsdtt_raster_surface_metrics_cpu`; `moving_window_iqr` and `compute_nvd` are new.

**`madk2_result`** — Trevisani, Teza & Guth (2023):
- `omnidirectional_roughness` — second-order directional MAD; slope-invariant
- `direction_maximum_continuity` — compass bearing of maximum continuity (°, [0, 180))
- `roughness_anisotropy` — directional anisotropy index [0, 1]

`calculate_madk2(raster, window_radius_metres, lag_pixels=1, n_threads=0)` → `madk2_result`

**`moving_window_iqr(raster, radius_metres, shape=circular, n_threads=0)`** → `lsdtt_raster_cpu`  
IQR (Q3 − Q1) of valid pixel values within a moving window. Pass a residual
raster from `lsdtt_raster_surface_fitting_cpu` for detrended roughness.

**`nvd_result`** — Normal Vector Dispersion (Woodcock 1977):
- `s1`, `s2`, `s3` — normalised eigenvalues of the 3×3 orientation tensor (sum to 1)
- Low `s3` → co-planar normals (flat or sloping smooth surface)
- High `s3` → scattered normals (rough terrain)

`compute_nvd(d_coeffs, e_coeffs, radius_metres, shape=circular, n_threads=0)` → `nvd_result`  
Input `d` and `e` are ∂z/∂x and ∂z/∂y rasters from any polyfit workflow.

---

## Chi analysis stack (`src/cpu/lsdtt_chi_*`)

Four files implement the full chi workflow in a pipeline:

```
lsdtt_chi_extraction_cpu      — extract chi values
       ↓
lsdtt_chi_network_analysis_cpu — piecewise-linear fit
       ↓
lsdtt_chi_disorder_cpu         — disorder metrics & uncertainty
lsdtt_chi_data_maps_cpu        — output: per-node points, polylines, basin polygons
lsdtt_chi_ksn_inversion_cpu    — spatial k_sn / u* inversion (requires Eigen3)
```

### lsdtt_chi_extraction_cpu

**`lsdtt_chi_options`** — `m_over_n = 0.45`, `A_0 = 1.0`

Three extraction modes:
1. **Full raster** — whole-raster chi from all baselevel outlets at chi = 0
2. **Outlet-seeded, zero base** — upslope domain of supplied outlets; each starts at chi = 0
3. **Outlet-seeded, inherited base** — outlets inherit chi values from a precomputed chi raster (for nested analyses)

Each mode has three domain variants: full upslope raster / threshold-defined channel / supplied junction network.

**Key result types:**
- `lsdtt_chi_basin_result` — outlet flat, per-node `flats`, `chi`, `receiver_flats`
- `lsdtt_chi_collection_result` — vector of `lsdtt_chi_basin_result`

Custom-denominator overloads accept a `denominator_raster` in place of drainage area (e.g. discharge).

### lsdtt_chi_network_analysis_cpu

Recursive piecewise-linear fitting of chi-elevation profiles basin by basin.

**`lsdtt_chi_network_fit_options`:**
- `critical_divergence_m` — maximum residual (m) before a segment is split (default 5.0 m)
- `mainstem_choice` — `max_flow_length` (default) or `max_chi`

**Key result types:**
- `lsdtt_chi_network_fit_result` — per-node: `chi`, `observed_elevation`, `fitted_elevation`, `residual_elevation`, `branch_id_of_node`, `is_main_stem`, `is_split_point`, plus segment list
- `lsdtt_chi_network_fit_segment` — segment endpoints in chi/elevation space, slope (k_sn = dz/dchi), `max_abs_residual`

**Key functions:**
- `fit_chi_network(chi_network, dem, options)` → `lsdtt_chi_network_fit_result`
- `export_chi_network_fit_nodes(...)` → `lsdtt_chi_network_fit_node_table`
- `export_chi_network_fit_segments(...)` → `lsdtt_chi_network_fit_segment_table`
- `write_chi_network_fit_segments_vector(...)` — writes polylines with k_sn attributes

### lsdtt_chi_disorder_cpu

Disorder metric (non-collinearity of chi-elevation branches).

**Analysis tools:**
- `get_chi_disorder(fit_result, options)` → `float` — full-network disorder
- `get_chi_branch_subset_sensitivity(fit_result, options)` — main stem + N side-branch combinations
- `get_chi_disorder_jackknife(fit_result, options)` — leave-one-branch-out jackknife
- `get_chi_disorder_bootstrap(fit_result, options)` — bootstrap over non-mainstem branches

**`lsdtt_chi_disorder_sweep_row`** — one row per basin × m/n for a sweep table written to CSV.

All result types have collection overloads (`_collection_result`) that operate basin-by-basin.

### lsdtt_chi_data_maps_cpu

High-level output for chi analysis.

- `write_chi_data_maps(chi_networks, flow, dem, output_prefix, fmt)` — writes per-node points, basin outlets, channel heads, and basin polygons. Incomplete basins (upslope area touches nodata) are automatically excluded.
- `write_chi_data_map_lines(...)` — writes channels as polylines keyed by `source_key` / `basin_key`

**Source-key convention:** the longest channel in each basin retains its headwater key to the outlet; shorter tributaries keep their own keys.

### lsdtt_chi_ksn_inversion_cpu *(requires Eigen3)*

Spatial k_sn / u* inversion with m/n model selection via regularised sparse least squares.

For each trial m/n:
1. recompute chi on the fixed channel network
2. assemble sparse design matrix from downstream chi increments
3. solve regularised least squares for a coarse spatial u* field
4. compute predicted elevations and MAE/RMSE misfit

**`lsdtt_chi_ksn_inversion_options`:** `A_0=1.0`, `coarse_cell_size=2000 m`, `roughness_weight=1e4`, `diagonal_jitter=1e-8`, `use_mean_absolute_error_for_model_selection=true`

**Key functions:**
- `make_chi_ksn_inversion_network(chi_basin, dem)` — build fixed network input
- `invert_chi_ksn(flow, network, tested_m_over_n, options)` → `lsdtt_chi_ksn_inversion_result`

Result includes per-trial `MAE`/`RMSE`, best-fit m/n, node-wise `chi`, `predicted_elevation`, `node_u_star`, and raster outputs `u_star_field` + `support_count_field`.

---

## lsdtt_channel_network_analysis_cpu

File: `src/cpu/lsdtt_channel_network_analysis_cpu.hpp/.cpp`

Channel link bearing analysis and junction angle analysis.

### Link sampling

`lsdtt_link_sampling_options` — controls how much of a link is used for bearing fitting:
- `full_link` — all nodes
- `fixed_n_pixels` — up to N pixels from the anchor end
- `fixed_length` — up to L metres from the anchor end

`trace_channel_link_between_junctions(flow, network, upstream_junction)` → `lsdtt_traced_channel_link`

### Channel link bearings

`calculate_channel_link_bearings(flow, network, junctions=all, sampling_options, anchor=upstream_end)` → `std::vector<lsdtt_channel_link_bearing>`

Each record includes: upstream/downstream junctions and flats, stream order, drainage area, bearing (rad and deg), `fit_linearity`, and a `representative_flat` for visual output.

Writers: `write_channel_link_bearings(...)`, `write_channel_link_bearing_vectors(...)` (visual QA lines)

### Junction angles

`calculate_junction_angles(flow, network, junctions=all, sampling_options)` → `std::vector<lsdtt_junction_angle_record>`

Each record includes (using descriptive names that map to Strong & Mudd 2023):
- `donor_opening_angle` = α (angle between two donor channels)
- `donor_1_receiver_angle` = β₁, `donor_2_receiver_angle` = β₂
- `bending_angle_signed/abs` — direction of the donor-opening midpoint from the inward receiver, `wrap(180 + (beta_1 - beta_2)/2)` (symmetric junction = ±180); NOT Strong & Mudd's beta; new diagnostic, excluded from the statistics

Writer: `write_junction_angle_data(...)`, `write_junction_angle_bearing_vectors(...)` (visual QA)

---

## lsdtt_strahler_tokunaga_analysis_cpu

File: `src/cpu/lsdtt_strahler_tokunaga_analysis_cpu.hpp/.cpp`

**Terminology:** A *Strahler channel* is a maximal downstream chain of compact
links sharing the same Strahler order. It differs from a compact link.

**Key types:**
- `lsdtt_strahler_channel` — geometry (flats, xy, cumulative distance), metrics (length, elevation drop, mean slope), and Tokunaga receiver info
- `lsdtt_strahler_channel_network` — container with `channel_id_of_upstream_junction` lookup

**Pipeline:**
1. `build_strahler_channel_network(flow, network, elevation)` — aggregates compact links, computes longitudinal metrics
2. `assign_tokunaga_receiver_channels(strahler_network)` — populates `receiver_channel_id`, `receiver_order`, `tokunaga_offset`
3. `analyse_tokunaga_self_similarity(strahler_network)` → `lsdtt_tokunaga_self_similarity_test`
   - pair counts → offset statistics → fit of T_k = a · c^(k−1)

**Writers:**
- `write_strahler_channels(...)` — Strahler-focused attribute schema (order, length, slope)
- `write_tokunaga_channels(...)` — Tokunaga-focused attribute schema (donor order, receiver order, offset); same geometry
- `write_tokunaga_summary_data(...)` — pair counts, offset statistics, and fit results to CSV

---

## lsdtt_cosmo_raster_cpu

File: `src/cpu/lsdtt_cosmo_raster_cpu.hpp/.cpp`

Raster-level cosmogenic nuclide operations.  All free functions in `namespace lsdtt`.

### Production raster

`calculate_production_raster(dem, nc={}, age_yr=0)` → `lsdtt_raster_cpu`

Computes the Stone (2000) spallation scaling factor at every valid pixel.
`nc` is an optional `CrnNcData` struct loaded via `lsdtt_crn_parameters::load_nc_data()`;
when `nc.loaded == false` the standard-atmosphere formula is used with no palaeomagnetic
correction.  `age_yr > 0` activates the Muscheler VDM time-integrated correction.

### Steady-state concentration — single basin

`calculate_crn_concentration_raster(nuclide, muon_scheme, erosion, prod, topo, snow, self, flow, outlet_node)` → `lsdtt_raster_cpu`

Computes steady-state concentration at every pixel **upslope** of `outlet_node`.
Used by the erosion-rate inversion workflow in `lsdtt-cosmo-tool`.

### Steady-state concentration — whole raster

`calculate_crn_concentration_raster_whole(nuclide, muon_scheme, erosion, prod, topo, snow, self)` → `lsdtt_raster_cpu`

Same physics as the single-basin version but operates on **every valid pixel** with no
flow-routing dependency.  OpenMP-parallelised over all pixels.  Used by
`lsdtt-cosmo-prediction`.

### Step-change concentration — whole raster

`calculate_crn_concentration_raster_step_change(nuclide, muon_scheme, time_since_change, prod, topo, eff_erate_old, eff_erate_new)` → `lsdtt_raster_cpu`

Transient surface concentration after an instantaneous step change in erosion rate.
Per-pixel `time_since_change` raster (yr); uniform scalar `eff_erate_old` and
`eff_erate_new` (g cm⁻² yr⁻¹).

At each pixel where `time_since_change > 0`:
1. Starting depth = `eff_erate_new × time_since_change`
2. Initialise SS at `eff_erate_old` at that depth
3. Evolve for `time_since_change` yr at `eff_erate_new` → surface concentration

At pixels where `time_since_change ≤ 0`: returns SS at `eff_erate_new`.

### Accumulated concentration — single outlet (scalar)

`calculate_accumulated_crn_concentration(conc, erosion, flow, outlet_node)` → `double`

Erosion-weighted mean concentration across all upslope pixels of `outlet_node`.

### Accumulated concentration — whole raster

`calculate_accumulated_crn_raster(conc, erosion, flow)` → `lsdtt_raster_cpu`

Erosion-weighted mean concentration at **every pixel** using D8 flow accumulation:
`C_acc[i] = Σ_upslope(C × ε) / Σ_upslope(ε)`.  Calls `flow.accumulate()` twice
(O(N) each).  Used by `lsdtt-cosmo-prediction`.

### Erosion-rate inversion

`calculate_eff_erate_from_conc(measured_conc, nuclide, muon_scheme, prod, topo, snow, self, flow, outlet_node, predicted_conc&)` → `double`

Inverts a measured nuclide concentration to an effective erosion rate using
Newton-Raphson (analytic derivative) followed by bisection if N-R does not converge.
Precomputes a basin summary in O(N) so each N-R and bisection step costs O(4).

### Supported nuclides and muon schemes

Nuclides: `"Be10"`, `"Al26"`, `"C14"`

Muon schemes: `"Braucher"`, `"BraucherBorchers"` (default), `"Granger"`, `"Schaller"`, `"newCRONUS"`

---

## lsdtt_crn_parameters / lsdtt_crn_particle / lsdtt_crn_column

**`lsdtt_crn_parameters`** (`src/cpu/lsdtt_crn_parameters.hpp/.cpp`)

Physics constants for CRN production and decay.  Workflow per pixel:
1. `reset_scaling(muon_scheme)` — restore standard F values
2. `scale_F_values(total_shielding, nuclide)` — fold shielding into F arrays
3. Pass `const&` to `lsdtt_crn_particle`

Also provides: `stone2000sp()`, `stone2000sp_paleomag()`,
`atmospheric_pressure_standard()`, `atmospheric_pressure_era40()`,
`load_nc_data()`, `apply_production_scaling()`, standardisation maps.

**`lsdtt_crn_particle`** (`src/cpu/lsdtt_crn_particle.hpp`, header-only)

Lightweight particle tracking nuclide concentrations.  Methods:
- `update_10Be_ss / update_26Al_ss / update_14C_ss` — steady-state at point or depth-integrated
- `init_ss_all` — SS for all three nuclides at current depth
- `update_conc(dt, eff_erate, params)` — analytic transient update for all three nuclides
- `remove_mass(depth_m, density)` — advance particle to shallower depth

**`lsdtt_crn_column`** (`src/cpu/lsdtt_crn_column.hpp/.cpp`)

1-D transient column simulation.  Reads an erosion-history CSV (`time,rate,removal`).
Methods: `run_transient(depths, scaling, muon_scheme, initial_*)`, `steady_state(...)`,
`step_change(d_eff, eff_erate_old, eff_erate_new, t_step, scaling, muon_scheme)`.

---

## lsdtt_burn_streams_cpu

File: `src/cpu/lsdtt_burn_streams_cpu.hpp/.cpp`

Burns a user-supplied polyline vector into a DEM by extracting the DEM elevation
profile beneath each polyline and enforcing a monotone-decreasing (downhill)
elevation constraint along the profile.

**`stream_burn_mode`:** `burn_in_place` | `snap_to_network` | `snap_to_accumulation`

- `burn_in_place` lowers the DEM under each polyline so that its profile is
  non-increasing downstream. Used by both `lsdtt-channel-extraction` and
  `lsdtt-stream-burning`.
- `snap_to_network` is implemented in `lsdtt-stream-burning`, where it is the
  default. It compares the supplied vector network's junctions with junctions
  extracted from the DEM flow network, writes diagnostics of where they agree or
  diverge, and repairs gaps (`build_vector_networks`,
  `find_first_layer_junction_pairs`, `repair_junction_gaps`,
  `find_candidate_junctions_phase`, `compute_reach_divergence`,
  `repair_candidate_junctions`, `repair_diverged_reaches`).
  `lsdtt-channel-extraction` does not support it and exits with an error if it
  is requested.
- `snap_to_accumulation` is for accurate vector layers (e.g. OS VectorMap
  Local). It is implemented in `lsdtt-stream-burning`: it traces D8 along the
  vector, keeps reaches that follow the DEM thalweg, and bridges obstruction
  gaps (roads, embankments) with least-cost-path carving.

`burn_streams_in_place(dem, vector_path, n_threads=0)` → `std::vector<burned_stream_profile>`

Each `burned_stream_profile` contains `burned_stream_node` records ordered outlet-first by `flow_distance`, with fields: `flow_distance`, `dem_elevation`, `burned_elevation`, `x`, `y`.

**Writers:** `write_burned_stream_profiles_csv(...)`, `write_burned_stream_profiles_points(...)` (gpkg / flatgeobuf)

---

## lsdtt_flat_surface_finder

File: `src/cpu/lsdtt_flat_surface_finder.hpp/.cpp`

Identifies flat (low surface texture / roughness) surfaces from a DEM. "Flat"
means low texture — not necessarily level. Alluvial fans, terraces, and piedmont
surfaces qualify.

**Three complementary texture metrics:**
- **MADk2** — second-order directional differences; slope-invariant; recommended for all cases
- **Residual IQR** — IQR of residuals to a locally fitted reference surface; three detrending methods: `fftw_highpass`, `wavelet_quadratic`, `moving_window_mean` (FFTW-free fallback)
- **NVD s3** — smallest eigenvalue of the surface-normal orientation tensor; requires FFTW for the wavelet quadratic fit; disabled by default

**Thresholding:** each metric is thresholded at a user-specified percentile of its own distribution. The QQ-plot CSV output helps identify the break-point between flat and rough terrain in the lower tail. Metrics are combined with AND or OR logic.

**Scope:** this module answers "which ground is flat?" and nothing more. Sorting flat ground into fluvial levels needs a channel and a height reference, and lives in `lsdtt_reference_field_ladder_cpu` behind `lsdtt-fluvial-surfaces`.

---

## lsdtt_reference_field_ladder_cpu

File: `src/cpu/lsdtt_reference_field_ladder_cpu.hpp/.cpp`

Stratifies a REFERENCE FIELD into height bands and labels connected components
**within each band**. The engine behind `lsdtt-fluvial-surfaces`.

**One operation, two tasks.** Finding a valley floor and segmenting a terrace
flight are the same computation at different band counts: the valley floor is
the lowest band of a two-band cut, a terrace flight an n-band ladder on the same
population. There is deliberately no "valley mode" and no "terrace mode" — only
`level_count`.

**Why band before labelling.** The incumbent ordering labelled once over the
whole flat mask and consulted height afterwards, which forces the riser network
to be continuous: a road or gully through a riser lets two treads at different
heights merge into one component, and nothing downstream can separate them
again. Banding first makes that merge impossible rather than detectable.

**Key API:**
- `build_reference_field_ladder(field, admitted, request, source_id?, coord?)`
  → `ladder_result` — the band index raster plus one `ladder_window` per
  along-channel window, each carrying its status, breaks and band centres.
- `label_components_within_bands(band_index, admitted, field, min_pixels, conn)`
  → `banded_label_result` — labels are a GLOBAL RUNNING COUNTER 1..N assigned
  band by band. The band is **not** packed into the label; it travels in the
  band-index raster and in `patch_record`. Throws if the total patch count would
  reach 2^24, above which a label cannot round-trip through the raster's float32
  backing store and every join keyed on it would silently corrupt.
- `build_fitted_surface_reference_field(dem, fit_mask, request)` — fits a
  `cone`, RANSAC `plane` or `poly2` and returns the signed residual.

**`ladder_status` is a reported output, not an internal detail:** `resolved`,
`unresolved`, `forced`, `degenerate`. A ladder that could not resolve says so
rather than returning a number, because what may be claimed about it depends on
which branch was taken.

**Two numerics decisions that must not be "tidied":** band centres are MEDIANS
of each band's values rather than midpoints of the break edges, and every level
count in range is scanned rather than searched downward (admissibility is not
monotone under median centres). Both are argued in the header. This file WRAPS
`find_population_breaks` from `lsdtt_flat_surface_channel_classify_cpu` and must
not modify it — that function is shared, so a change there moves two tools at
once.

---

## lsdtt_raster_binary_mask_cpu / lsdtt_raster_skeleton_graph_cpu

**`lsdtt_raster_binary_mask_cpu`** (`src/cpu/lsdtt_raster_binary_mask_cpu.hpp/.cpp`)  
Generic binary mask operations (foreground = value > 0.5):
- `binary_morphological_dilation`, `erosion`, `opening`, `closing` (structuring element radius in metres)
- Reconstruction-based morphology, directional line-based closing
- Hole filling, size-based component filtering and labelling
- `filter_binary_components_by_seed_reachability` — a morphological opening
  CONSTRAINED BY DISTANCE TO A SEED SET: erode, drop the patches a seed cannot
  reach (either not adjacent to one, or only reachable by a detour longer than
  `max_detour_ratio` times the straight-line distance), then regrow with each
  patch capped at its own farthest distance from a seed. The seed is whatever
  the caller supplies — a channel network, a coastline, a road centreline — and
  nothing in the operation knows what the mask represents. Note it is **not** a
  filter: the regrow fills enclosed holes, so the result can contain more set
  pixels than the input. An in-place raw-vector overload exists for callers that
  already hold a byte mask and a distance field in an inner loop.

  This arrived from `lsdtt_valley_centreline_cpu`, where it had been cleaning
  valley floors against the channel that carved them. **General raster
  operations belong here and should be reusable, not duplicated in
  task-specific module pairs** — that move is why the valley engine went from
  1484 lines to 626.

Does NOT perform DEM feature scoring, polygon conversion, or skeleton analysis — those are in the three separate files below.

**`lsdtt_raster_binary_mask_vectorise_cpu`** — binary mask to vector polygon conversion  
**`lsdtt_raster_skeleton_graph_cpu`** — 1-pixel skeleton graph: degree classification, spur pruning, segment splitting, endpoint tangent estimation, endpoint pairing, bridge rasterisation

---

## Utility-layer additions (`src/utils/`)

**`lsdtt_distributions.hpp/.cpp`** — parametric distribution fitting on valid flat samples (no nodata filtering). Supported families: `normal`, `exponential`, `weibull`, `pareto`. All use MLE. KS goodness-of-fit test computed for each. `max_divergence_x` marks where the empirical CDF departs most from the fitted CDF — useful as a threshold candidate for roughness distributions.

**`lsdtt_linear_algebra.hpp`** — header-only lightweight linear algebra (no GDAL, no lsdtt3 headers):
- `sym3x3_eigenvalues(M)` — analytical Cardano eigenvalues of a real symmetric 3×3 matrix
- `sym2x2_eigen(a, b, c)` — eigenvalues + major-axis orientation of a symmetric 2×2 matrix (inertia-tensor / shape-orientation analysis)
- `invert3x3(M)` — Gauss–Jordan inversion
- `angle_between_vectors_degrees(a, b)`
- `point_to_segment_distance(p, a, b)`

**`lsdtt_polygon_geometry.hpp`** — header-only planar polygon geometry on rings of `lsdtt_xy` (depends only on `lsdtt_types.hpp` + `lsdtt_linear_algebra.hpp`):
- `polygon_signed_area` / `polygon_area`, `polygon_perimeter`, `polygon_centroid`
- `polygon_central_moments` — Green's-theorem mu20/mu11/mu02 (feeds `sym2x2_eigen`)
- `simplify_douglas_peucker` — strips the pixel-staircase before perimeter-based metrics
- `convex_hull` (Andrew monotone chain), `min_area_bounding_rect` (rotating calipers), `feret_diameters`, `equivalent_circle_diameter` (Heywood)
- `compute_basin_shape(...)` → `lsdtt_basin_shape` — headline drainage-basin shape indices (elongation, form factor, circularity, Gravelius compactness, solidity, moment eccentricity/orientation, hypsometric integral). Theory + references: `docs/tex_technical/basin_shape_metrics.tex`.

**`lsdtt_thread_utils.hpp`** — header-only: `resolve_thread_count(n_threads)` returns `n_threads` if > 0, otherwise the OpenMP or hardware maximum.

**`lsdtt_binary_array2d_tracing.hpp/.cpp`** — binary array tracing utilities (used by skeleton graph and related morphology code).

---

## lsdtt_landscape_drape_cpu

File: `src/cpu/lsdtt_landscape_drape_cpu.hpp/.cpp`

Computes an upper landscape envelope surface that drapes over the topography,
resting on ridges and summits, and derives a pixel-wise incision depth raster
(`envelope − DEM`).

**Algorithm: Verlet cloth simulation**

A flat virtual cloth is initialised at `max_elevation` and allowed to fall under
gravity onto the terrain. 8-connected spring forces keep the cloth smooth between
contact points. A particle is pinned when it descends within `class_threshold`
metres of the terrain surface. At equilibrium the cloth satisfies ∇²f = 0
between contact points — identical to harmonic interpolation through the ridge
points.  Verlet time-stepping naturally enforces the above-terrain constraint that
a direct solver cannot.

This is the upper-envelope equivalent of the Zhang et al. (2016) Cloth Simulation
Filter (CSF) for LiDAR ground classification, and replicates the geomorphological
intent of the Carretier & Valenti EnvelopSurface GRASS script without requiring
explicit crest extraction via inverted-DEM flow routing.

**Key types:**

```cpp
struct lsdtt_cloth_drape_options {
    float gravity           = 0.3f;    // downward force per step (m/step²)
    float rigidity          = 1.0f;    // spring stiffness — higher = smoother envelope
    float time_step         = 0.65f;   // Verlet dt
    float class_threshold   = 0.5f;    // contact tolerance (m); raise for noisy DEMs
    int   max_iterations    = 500;
    float convergence_tol   = 0.005f;  // max displacement (m) to declare settled
    int   n_threads         = 0;
};

struct lsdtt_cloth_drape_result {
    lsdtt_raster_cpu envelope;              // cloth equilibrium = upper landscape envelope
    lsdtt_raster_cpu incision;             // envelope − DEM (negative = DEM above envelope)
    lsdtt_raster_cpu crest_mask;           // 1 = pinned pixel (ridge / summit contact)
    int              iterations_taken;
    float            final_max_displacement;
};
```

**Public API:**
- `compute_cloth_drape(dem, opts)` → `lsdtt_cloth_drape_result`

**Parameter guidance:**
- `rigidity`: equivalent to the `Tension` parameter in the Carretier/Valenti script
  and to spring stiffness in CSF. Start at 1.0; increase (5–20) for broader,
  smoother envelopes over highly dissected terrain.
- `class_threshold`: increase (1–5 m) for DEMs with residual artefact pits or
  coarse resolution.
- `convergence_tol`: 0.005 m is appropriate for metre-scale DEMs; tighten to
  0.001 m if incision values are noisy at convergence.

**References:**
- Zhang W. et al. (2016) Remote Sensing 8(6):501. DOI: 10.3390/rs8060501
- Carretier S. & Valenti C. — EnvelopSurface.bash (IRD/GET, France)

---

## Landscape evolution model (LEM) stack

The LEM is a layered set of CPU engines driven by the `lsdtt-lem` CLI. The CLI
owns surface generation, conditioning, I/O, and the run loop; each physical
process lives in its own free-function module so the same engines can be reused
outside the LEM (e.g. the steady-state solvers are useful standalone).

```
cli/lsdtt-lem.cpp
   │  resolves the model grid, builds/loads DEM, conditions it, owns the time loop
   │  or steady-state branch
   ├── lsdtt_lem_input_grid_cpu        model grid (parameters / DEM / K raster) and
   │                                   placement of input rasters onto it
   ├── lsdtt_lem_parameter_rasters_cpu K values: calibration, conversion, validation,
   │                                   per-step scaling by the forcing history
   ▼
lsdtt_lem_cpu            run_lem_step(): one timestep = fluvial → hillslope → uplift → restore boundaries
   ├── lsdtt_stream_power_cpu        implicit stream-power incision (transient)
   ├── lsdtt_hillslope_transport_cpu transient hillslope transport
   └── (uplift + boundary restore handled inline in lsdtt_lem_cpu)

cli/lsdtt-lem.cpp  (run_steady_state branch)
   ├── lsdtt_steady_state_cpu            analytical chi-based fluvial steady state
   └── lsdtt_hillslope_steady_state_cpu  hillslope steady state on top of the fluvial result
```

All five engines are CPU free functions in `namespace lsdtt`, take a
`const lsdtt_raster_cpu&`, and return a new `lsdtt_raster_cpu` (input never
modified). They build their own flow object internally (D8 `lsdtt_flowinfo` or
D-inf `lsdtt_flowgraph_cpu`) from the elevation passed in.

### lsdtt_lem_cpu (orchestrator)

File: `src/cpu/lsdtt_lem_cpu.hpp/.cpp`

`run_lem_step(dem, opts, boundary_flats)` → `std::pair<lsdtt_raster_cpu, lsdtt_lem_step_stats>`.

Per-step sequence:
1. **Fluvial** — D8 (`compute_d8_stream_power_incision`) or D-inf
   (`compute_dinf_stream_power_incision`); D-inf takes priority and the CLI forces
   D8 off if D-inf is requested.
2. **Hillslope** (optional) — `adi` runs without a flow graph; `cellular` uses conservative
   soft-critical face transfers with optional nonlocal routing and sparse particles;
   `linear` / `roering` / `roering_implicit` need a D-inf graph only for channel identification, and skip
   building one entirely when an explicit channel mask is supplied (it is otherwise
   rebuilt every step, which dominates a long fixed-network run); `cidre` always
   builds one, because it routes sediment on the D-inf topological forward sweep.
3. **Uplift** (optional) — uniform `uplift_rate * timestep` added to valid pixels.
4. **Statistics** — computed full-domain *before* boundary restoration. When
   `compute_erosion_rate_field` is set, the per-pixel erosion rate (m/yr, uplift
   backed out, deposition negative) is returned in `lsdtt_lem_step_stats::erosion_rate`.
5. **Boundary restore** — pixels in `boundary_flats` are reset to their pre-step
   elevation (fixed-elevation baselevels), *less* `baselevel_fall_rate * timestep`.

`lsdtt_lem_options` carries all per-step knobs (fluvial K/m/n, threshold area,
hillslope model + diffusivity + Sc, uplift rate, base-level fall rate, process
switches, boundary conditions). `erodibility_field` is an optional per-pixel K
(flat raster order, absolute K in the run's convention; empty = the scalar).
`run_lem_step` checks it on every call (size, and finite K > 0 on every valid
pixel, serially before any parallel region) and caches nothing derived from it,
so a caller may rewrite it between steps. `lsdtt_lem_step_stats` reports `delta_max_elevation`,
`delta_mean_elevation`, `max_local_elevation_change`, `mean_erosion_rate`,
`max_erosion_rate`, and optionally the `erosion_rate` field.

The fluvial flow object is built with `allow_pits_as_outlets = true` because the
transient LEM surface may contain pits that must not stall routing.

**Channel mask.** `compute_lem_channel_mask(dem, opts)` returns the Dirichlet
channel set the hillslope step would use for a surface — the explicit mask when
`opts.explicit_channel_mask` is supplied, else the pixels whose D-inf drainage area
reaches `incision_threshold_specific_area`. It shares its implementation with the hillslope
step, so the mask a run *writes* is exactly the one it *used*. The CLI exposes it as
`write_channel_mask_raster` (with an independent `channel_mask_threshold_specific_area`, so a
run that erodes fluvially everywhere can still emit a network), which lets a fluvial
run freeze its network and feed it back as `channel_mask_fname` to a hillslopes-only
run. `hold_channel_mask_fixed` additionally adds the masked pixels to `boundary_flats`,
so uplift cannot lift the frozen network away from the hillslopes relaxing toward it.

**Snap helper.** `snap_lem_surface_to_fluvial_steady_state(dem, ss_opts, use_dinf,
fi_opts, fg_opts, fill_min_slope, boundary_conditions, n_threads, before_solve)`
fills, routes with the caller's flow options and solves the analytical steady
state; it returns the result and the drainage area. Both the CLI's initial snap
and each cycle of its legacy steady-state branch call it. `before_solve` receives
the filled surface and a mutable copy of the options just before the solve: the
hook for a per-step K source (lithology).

### Time-varying forcing (pulses of incision)

`lsdtt_lem_forcing_history` (in `lsdtt_lem_cpu`) makes `uplift_rate` and
`baselevel_fall_rate` functions of model time, read from a CSV
(`forcing_history_fname`) with columns `time_yr, uplift_rate, baselevel_fall_rate`.
It is applied as a **piecewise-constant step function** — the row in force at time
`t` is the last one with `time_yr <= t` — so a pulse has sharp on/off edges rather
than being smeared across the surrounding rows by interpolation. The CLI resolves
the forcing at the start of each step, and records what it applied in the
`uplift_rate` / `baselevel_fall_rate` columns of the stats CSV.

`baselevel_fall_rate` lowers the fixed-elevation boundary pixels instead of
restoring them, which drops base level away from the domain and sends a wave of
incision up the channel network. Together with `write_erosion_rate_raster` and
`write_channel_mask_raster` (both written per frame at every
`frequency_of_output_writing`), this is what `lsdviztools3`'s
`render/lem_pulse.py` animates: the pulse is essentially invisible in elevation
but obvious in the erosion-rate field.

An `erodibility` column changes the scalar K in time. With a K field active the
CLI keeps the loaded field as a base and, at the start of each step, sets
`opts.erodibility_field = scale_lem_erodibility_field(base, K_param, K_forcing(t))`,
i.e. K(x, t) = K_base(x) K_forcing(t) / K_param. The CLI marks reserved hooks
H1-H4 (time step, initial snap, start of the steady-state branch, each legacy
cycle) where a per-step K source rewrites the base field. An exception inside a
step is rethrown to one handler that prints `FATAL: time step failed ...` and
exits 1.

### lsdtt_lem_parameter_rasters_cpu (per-pixel K values)

File: `src/cpu/lsdtt_lem_parameter_rasters_cpu.hpp/.cpp`. What K values mean and
whether they are acceptable; never where pixels sit (that is the input-grid
module). No GDAL headers; returns notes and warnings and never prints.

- `check_lem_erodibility_calibration(stated, run, label)`: the single comparison
  of a source's stated m, n, area convention and calibration cell size against
  the run's. Returns the conversion factor (exactly 1.0 when nothing converts),
  notes and warnings; a different m or n is fatal.
- `validate_lem_erodibility_values` / `validate_lem_erodibility_field`: values on
  the DEM grid → calibration → conversion → checks (nodata on a DEM-valid pixel
  and K ≤ 0 are refused, naming the pixel). The `_field` form first applies
  `check_lem_field_grid_alignment`.
- `load_lem_erodibility_field(path, dem, placement, opts)`: prepare onto the grid,
  then validate.
- `check_lem_erodibility_field`: the per-step check `run_lem_step` uses.
- `scale_lem_erodibility_field(base, reference, applied, out)`:
  `float((double(b) · applied) / reference)`, which leaves `base` bitwise unchanged
  when applied == reference; throws `range_error` naming the pixel on overflow or
  underflow.

The units and conversion factor themselves (`describe_erodibility_units`,
`erodibility_area_conversion_factor`, `parse_erodibility_area_convention`) live in
`lsdtt_stream_power_cpu`, as properties of the law.

### lsdtt_lem_input_grid_cpu (model grid and input-raster placement)

File: `src/cpu/lsdtt_lem_input_grid_cpu.hpp/.cpp`. The LEM step, the steady state
and every validator only see rasters already on the model grid; getting an input
raster there happens here, before the LEM, never inside a step. Uses only
GDAL-free core headers.

- Model grids: `make_lem_synthetic_grid` (Watopia when no WKT is given; the CLI's
  synthetic branch calls it), `make_lem_model_grid_from_raster` (mode 2: the
  source's grid; refuses no CRS, geographic, non-metre and distorting CRSs),
  `make_lem_model_grid_from_raster_dimensions` (mode 1 (a)),
  `apply_lem_source_footprint` (mode 2: source nodata becomes DEM nodata).
- Grid functions take an `lsdtt_lem_source_grid` (meta + band count), never a
  path, so a multi-band source can supply a grid. `describe_lem_source_grid`
  builds one from an open loader.
- `check_lem_field_grid_alignment`: the strict check (`same_grid` at 1e-6 m plus
  `crs_equal_semantically`).
- `check_lem_projection_distortion`: refuses EPSG:3857 by name and any CRS whose
  point scale factor at the grid's centre or corners is more than 1 % from 1
  (`lsdtt_raster_meta::max_linear_scale_error`).
- Clip windows (`compute_lem_clip_window`, `compute_lem_clip_extent`,
  `check_lem_clip_extent_within`): the only implementation of the centre /
  lower-left rules. Equal cell size gives a pixel-aligned window counted from the
  source's south-west corner (a point on a grid line belongs to the pixel to its
  north / east); a different cell size gives an extent the source is resampled
  onto. A window off the source is refused.
- `prepare_lem_input_raster(path, model_grid, opts, label)`: the single entry
  point for the five placements; the result's meta is the model grid.

### lsdtt_lithology_map_cpu (lithology maps to erodibility)

File: `src/cpu/lsdtt_lithology_map_cpu.hpp/.cpp`. No GDAL; never prints.

- `load_lithology_property_table(path)`: `#` header lines (`erodibility_convention`,
  `m_exponent`, `n_exponent` required; `calibration_cell_size_metres`, `provenance`
  optional) then `code,name,erodibility[,reference]`. Codes 0..65534 and names
  unique; K finite and > 0; errors name the line (and column). The stated
  calibration is an `lsdtt_erodibility_calibration`, checked later by the LEM's
  `check_lem_erodibility_calibration` when the raster is used.
- `load_lithology_crosswalk(path)`: `value,code`.
- `lithology_codes_for_values(values, counts, crosswalk*, table)`: crosswalk, else
  integer code, else table name (case-insensitive); lists every unmapped value with
  its pixel count.
- `make_lithology_erodibility(grid, value_index, values, codes, table, opts)`: the
  code and K rasters on `grid` (whose mask is the domain) and a legend. Unmapped
  domain pixels are refused, given `unmapped_code`, or (footprint mode) left
  nodata.

Lithologic units with no property table, for analyses categorised by unit
(a DEM and a geologic map are all they need):

- `read_lithology_units_on_grid(grid, source)`: the unit code of every valid pixel
  (-1 = unmapped) and a label per code (name, map values, pixel count). Vectors go
  through `rasterise_polygon_field_values` (the field must exist); a code raster is
  read as it is when already on the grid, else warped onto it with nearest
  neighbour. Only map values reaching a valid pixel need a code.
- `lithology_unit_codes_for_values(values, counts, crosswalk, names)`: with a
  crosswalk, its codes, named by its optional `name` column (else its values joined
  with "; "); without, integer values are their own codes and otherwise the sorted
  values are numbered 1..N. Codes are local to a run; the unit name is the key
  between runs.
- `lithology_source_from_path`, `load_lithology_crosswalk_names`,
  `lithology_unit_code_raster`.

The parameters are declared once in `lsdtt_common_parameters`
(`declare_lithology_map_parameters`, `declare_write_lithology_code_raster_parameter`);
`lsdtt-lithology-erodibility` adopts them and overwrites the guidance text with its
table-name matching, so its `-h json` is unchanged.

The categorical polygon rasteriser it relies on is
`rasterise_polygon_field_values` in `lsdtt_vector_rasterise_cpu`: cell-centre
assignment (overlaps to the lower sorted value), then a coverage pass on an
oversampled subgrid, in 512² chunks, that fills sliver gaps (fragments of one
value add up). The older `rasterise_polygon_string_field` (STI) and
`rasterise_macrostrat_geology` are separate implementations still.

### lsdtt_lithology_channel_stats_cpu (k_sn by lithologic unit)

File: `src/cpu/lsdtt_lithology_channel_stats_cpu.hpp/.cpp`. GDAL-free; never prints.
Used by `lsdtt-chi-analysis` when `lithology_fname` is set (needs
`fit_chi_networks: true`).

- `assign_chi_fit_nodes_to_segments(fit)`: segments share end nodes; a node belongs
  to the segment it lies on except that segment's downstream end, which belongs to
  the next segment down (a confluence to the receiving branch); the outlet to the
  last segment.
- `compute_lithology_contact_distance(grid, units)`: distance to the nearest
  contact pixel (a valid pixel with a valid 4-neighbour of another unit; unmapped
  counts as a unit), through `euclidean_distance_transform`.
- `build_chi_lithology_node_table`: every fitted node with its segment k_sn, unit,
  flow length to its receiver and distance to contact. Observation only: the fits
  are not changed and the existing fit tables stay byte-identical.
- `summarise_chi_ksn_by_lithology`: per unit, unmapped and all: node, segment and
  basin counts, channel length, median, quartiles, mean, range and length-weighted
  mean of the node k_sn; `low_sample` below 5 segments.
- `make_chi_ksn_histogram_by_lithology`: common bins (linear or log10) over the
  pooled 1st-99th percentile, with below and above rows, so units overlay.

k_sn per node is its fitted segment's; segments are not split at contacts.

### lsdtt_macrostrat_fetch and lsdtt_http_fetch (Macrostrat map units)

Files: `src/core/lsdtt_macrostrat_fetch.hpp/.cpp`, `src/core/lsdtt_http_fetch.hpp/.cpp`.
GDAL-free headers; never print (notes and warnings are returned).

Macrostrat's map API answers per point, so an area comes from its vector tiles
(`https://tiles.macrostrat.org/carto/{z}/{x}/{y}.mvt`, layer `units`, EPSG:3857).

- `http_request_body(url, post_body, opts, body, error)`: one GET (or POST) through
  GDAL's curl with timeouts and retries. `lsdtt_raster_fetch` uses it too.
- `slippy_tiles_for_bbox`, `slippy_tile_bounds_3857`, `utm_epsg_for`: tile maths.
- `parse_macrostrat_lithologies` / `_map_unit_liths` / `_sources`: the API's
  `success/data` JSON (malformed JSON or an `error` response throws).
- `macrostrat_primary_lithology(lith_text, lith_ids, vocabulary)`: the first
  lithology of a unit: the first `Major:{...}` term, else the first linked
  vocabulary name, else the whole text, else empty; lower case, whitespace
  collapsed. `lith_primary_source` records which.
- `fetch_macrostrat_units(opts)`: downloads the tiles (or reads
  `<tile_dir>/{z}/{x}/{y}.mvt|.pbf`), each through a `/vsimem` path carrying z/x/y
  because that is how the MVT driver georeferences a single tile; keeps features
  whose envelope overlaps the box; batches the API lookups 100 map ids at a time;
  reprojects to the UTM zone of the box centre (or `target_epsg`) and writes layer
  `units`. Too many tiles, no readable tile or no units are errors; a failed API
  lookup is a warning and `lith_primary` falls back to the tiles' text.
- `macrostrat_attribution_text(result)`: the CC BY 4.0 citation (Peters, Husson &
  Czaplewski 2018) and each source map.

### lsdtt_octopus_fetch (OCTOPUS basin erosion rates)

Files: `src/core/lsdtt_octopus_fetch.hpp/.cpp`. GDAL-free header; never prints.

- The endpoint (`http://geoserver.octopusdata.org/geoserver/wfs`, HTTP only) and
  the layers (`be10-denude:crn_int_basins`, `be10-denude:crn_int_outlets`) are
  constants.
- Pure recipe: `octopus_wfs_bbox_parameter` (latitude first,
  `S,W,N,E,urn:ogc:def:crs:EPSG::4326`; other forms return no features),
  `octopus_wfs_hits_url`, `octopus_wfs_page_url` (`count` / `startIndex`, sorted
  by `OBSID1`; never `maxFeatures`, which the server ignores),
  `octopus_wfs_outlets_by_id_url`, `octopus_wfs_page_starts`,
  `parse_octopus_wfs_number_matched`, `octopus_wfs_exception_text`,
  `is_octopus_no_data_sentinel` (negative all-nines values with at least three
  nines) and `octopus_rate_fields(scheme, nuclide)`.
- `fetch_octopus_erosion_rates(opts)`: hits, then pages, each body staged in
  `/vsimem` and read with the GeoJSON driver; or the GeoJSON bodies of
  `response_dir` offline. Returns the basins with their typed attributes, the
  per-field counts of nulled sentinels and the calibration rows.
- Writers: `write_octopus_erosion_rate_csv`, `write_octopus_basin_polygons`,
  `octopus_dem_bounds` / `write_octopus_dem_bounds`, `octopus_attribution_text`.

`bounds_from_vector_file` (`lsdtt_raster_fetch`) also takes a raster: when the
file opens as no vector layer, the WGS84 envelope of its densified footprint is
the box (`geometry_kind` `raster`).

### lsdtt_erodibility_calibration_cpu (stream-power K from basin erosion rates)

Files: `src/cpu/lsdtt_erodibility_calibration_cpu.hpp/.cpp`. Never prints; takes
`const lsdtt_flow_base&`. Used by `lsdtt-chi-analysis` (`erosion_rate_fname`).

- Input: `canonical_erosion_rate_units` / `erosion_rate_units_to_m_per_yr` (`mm_per_kyr`
  = `m_per_myr` 1e-6, `cm_per_kyr` 1e-5, `mm_per_yr` 1e-3, `m_per_yr` 1; no default,
  mass fluxes refused) and `read_erosion_rate_samples` (any point input through
  `load_point_layer_with_attributes`; a per-row `erosion_rate_units` column must agree).
- `delineate_erosion_rate_basins` (m/n independent): snapping
  (`snap_pixel_by_method`, `reported_area` by default), basin completeness, DEM vs
  reported area, nested and duplicate samples, nesting groups (union-find), the
  nested-basin policy, lithologic composition.
- `compute_erosion_rate_basin_ksn` (per m/n): the reach median of each sample basin's
  own fitted chi network on a single-unit grid (reaches split at confluences only), its
  reach bootstrap, and the diagnostic statistics.
- `calibrate_erodibility_from_erosion_rates` (per m/n): K per sample (Monte Carlo over
  E and the k_sn replicates), K at fixed n (weighted log mean), a fitted (n, K)
  (`fit_york_line`), K per unit (`mono`, `unit_in_basin`); intervals by a bootstrap
  over nesting groups from `erodibility_calibration_bootstrap_engine` (stream tag 2).
- Writers: the samples, basin K, fit, by-lithology and basin-lithology CSVs, the basin
  polygons, and the draft property table (`write_lithology_erodibility_property_table`).
- Parameter groups: `declare_erosion_rate_input_parameters`,
  `declare_erodibility_calibration_parameters`, `declare_stream_power_n_values_parameter`
  (`lsdtt_common_parameters`). Shared: `stream_power_erodibility_from_erosion_rate`
  (`lsdtt_lithology_relative_erodibility_cpu`).

### lsdtt_stream_power_cpu (transient incision)

File: `src/cpu/lsdtt_stream_power_cpu.hpp/.cpp`

Implicit stream-power incision, `U = K·A^m·S^n` (Braun & Willett 2013 generalised
to D-inf per Anand et al. 2020). Both functions sweep the topological order in
**reverse** (outlets first) so a receiver elevation is already updated when its
donor is solved.

- `compute_d8_stream_power_incision(dem, flowinfo, K, m, n, dt, threshold_area, divide_area_by_dx, slope_clamp, ..., erodibility_field = {})`
  — single-receiver D8; serial reverse-topological sweep. A non-empty
  `erodibility_field` (n_rows·n_cols values) gives each pixel its own K in the
  one coefficient K·A^m every update path reads, so a uniform field is bitwise
  the scalar run.
- `compute_dinf_stream_power_incision(..., flowgraph, ..., n_threads, ..., erodibility_field = {})`
  — D-inf with up to two receivers; level-parallel (`build_reverse_levels`) sweep,
  OMP within each level.

`n = 1` uses the closed-form implicit update; `n ≠ 1` uses Newton-Raphson with an
overexcavation clamp to `receiver_elevation + slope_clamp·dx_eff`. Diagonal
directions use `dx_eff = dx·√2`.

### lsdtt_steady_state_cpu (analytical fluvial steady state)

File: `src/cpu/lsdtt_steady_state_cpu.hpp/.cpp`

Chi-based analytical solution of `U = K·A^m·S^n` (Mudd et al. 2014). `z[i] =
baselevel + f·chi[i]`, where chi is integrated up-network and is independent of
`f` (so it is computed once from geometry).

- `compute_d8_fluvial_steady_state(flowinfo, opts)` — reverse topological sweep.
- `compute_dinf_fluvial_steady_state(flowgraph, opts)` — level-parallel sweep;
  two-receiver chi uses `chi[i] = chi[e1] + dx·√(A^(-2m/n) − dchi²)` with a
  cardinal-only fallback when the discriminant is negative.

`lsdtt_fluvial_steady_state_options` has two modes: **K-mode** (user supplies K,
`f = (U/K)^(1/n)`) and **relief-mode** (user supplies `target_relief`, K is
back-calculated so `f·MaxChi = target_relief`). Result reports `K_used`,
`max_unit_height` (= MaxChi), and `computed_relief`. The CLI can iterate this over
`steady_state_n_cycles`, pit-filling between cycles, so divides migrate one pixel
per cycle and the single-pass elevation breaks across divides are smoothed.

`opts.erodibility_field` (optional, flat on the flow grid) turns chi into a
K-weighted chi: each increment is multiplied by `w_i = (K/K_i)^(1/n)` (K being
`opts.K`, the reference), precomputed once and validated serially (finite
K_i > 0 on every valid pixel) before the parallel sweep. With no field the weight
is exactly `1.0f`, so the scalar arithmetic is unchanged. In relief mode the
field is scaled by `erodibility_field_scale = K_used / opts.K`, reported in the
result.

### lsdtt_hillslope_transport_cpu (transient hillslope)

File: `src/cpu/lsdtt_hillslope_transport_cpu.hpp/.cpp`

Three transient hillslope laws, selected by `hillslope_transport_model`:

- `compute_cidre_hillslope_transport(dem, flowgraph, kappa, dt, nonlinear, Sc, threshold_area, divide_area_by_dx, deposit_on_hillslope_baselevels, n_threads)`
  — Carretier et al. (2016) finite-volume forward topological sweep over the D-inf
  graph; flux accumulated via atomic adds. Nonlinear branch uses
  `L = 1/(1−S²/Sc²)` and switches to slope-restoration when `S ≥ Sc`
  (unconditionally stable near Sc).
- `compute_roering_hillslope_transport(dem, is_channel, kappa, dt_total, nonlinear, Sc, n_threads)`
  — explicit FD of `D·∇·[∇z/(1−(|∇z|/Sc)²)]` on the 4-connected grid; conservative
  divergence form; auto-subcycles to satisfy the CFL bound `dt ≤ h²/(4·D·α_max)`.
  `nonlinear == false` reduces to the linear Laplacian.
- `compute_adi_hillslope_transport(dem, kappa, dt, boundary_conditions)`
  — Crank-Nicolson linear diffusion via uncoupled ADI (Thomas algorithm) sweeps;
  no flow graph required.

Channel pixels (drainage area at or above `incision_threshold_specific_area`) act as hillslope sinks
/ Dirichlet nodes; the caller builds the `is_channel` mask for the FD laws.

### lsdtt_hillslope_steady_state_cpu (steady-state hillslope)

File: `src/cpu/lsdtt_hillslope_steady_state_cpu.hpp/.cpp`

Solves the hillslope steady state on top of a fluvial result, with channel pixels
fixed as Dirichlet boundaries. Two interchangeable solvers (same signature):

- `compute_hillslope_steady_state(dem, is_channel, opts)` — Red-Black SOR solving
  `D·∇²z = −U` (linear) or the Roering nonlinear law with Picard-lagged
  coefficients. Auto-tuned ω (capped 1.95); O(N) convergence.
- `compute_cidre_hillslope_steady_state(dem, is_channel, opts)` — pseudo-transient
  Cidre iterated to a fixed point, rebuilding the D-inf graph each iteration
  (Picard update of routing). At steady state the nonlinear Cidre and Roering
  solvers produce the same field (Carretier et al. 2016, §3.3).

`lsdtt_hillslope_steady_state_options`: `D`, `U`, `Sc`, `nonlinear`,
`boundary_conditions`, `tolerance`, `max_iterations`, `sor_omega` (0 = auto),
`n_threads`. Result reports `n_iterations`, `final_residual` (max|Δz|), and
`converged`. In the CLI the solver is chosen by `hillslope_ss_model`
(`linear`/`roering` → SOR, `cidre` → pseudo-transient); for `cidre` the nonlinear
flag comes from `enable_nonlinear_hillslopes`.

---

## Graph-based hydraulic routing (no vendored code)

There is **no `graphflood/` directory and no vendored third-party code**. An
earlier version of lsdtt3 carried a clone of the graphflood headers, ported in
from upstream. That is gone: Boris Gailleton, graphflood's author, is a
developer of lsdtt3 and wrote a native implementation for it.

The routing now lives in `src/cpu/lsdtt_hydraulic_graph_cpu.hpp/.cpp` as
first-party AGPL-3.0-or-later code implementing `lsdtt_hydraulic_base`, like any
other CPU module. It includes nothing outside `src/`. The method is still
Gailleton et al. (2024) and is cited as such in the CLI and the module header —
an attribution question, not a licensing one.

This matters for packaging: the whole tree is single-licence, which is why
`debian/copyright` needs no third-party stanza.

---




## Key design principles

1. No `using namespace std` anywhere
2. All code in `namespace lsdtt`
3. `check_projected_metres()` called at entry of metric analysis functions
4. Mask array (`uint8_t`, 0=nodata, 1=valid) separate from data array
5. `meta.is_nodata(v)` handles both sentinel values and NaN
6. No warnings — `-Wall -Wextra -Wpedantic` enforced
7. FFTW double (`fftw_*`), not float; not threaded FFTW
8. Analysis free functions take `const lsdtt_raster_cpu&` — not the base class
9. The concrete D8 flow class is `lsdtt_flowinfo` (defined in `lsdtt_flowinfo_cpu.hpp`
   — note: class name has no `_cpu` suffix; file name does).
   CLIs that only need the abstract interface use `make_flow_base()` and hold
   `unique_ptr<lsdtt_flow_base>`. CLIs that need the full concrete API throughout
   (LEM, sediment routing) construct `lsdtt_flowinfo` directly.
   Never write `lsdtt_flowinfo_cpu` — that type does not exist.
10. `src/core/` = abstract interfaces, factories, shared types, I/O only
11. `src/cpu/` = all CPU implementations (raster, flow, and algebra)
12. `src/gpu/` = all GPU implementations
13. All lightweight shared primitive types live in `src/core/lsdtt_types.hpp`.
    Do not create a new `*_types.hpp` file — add to lsdtt_types.hpp instead.
    A separate types file is only justified when inclusion would create an
    unwanted compile-time dependency; that situation should be extremely rare.
14. `src/core/lsdtt_crs_suggest` inspects CRS; `src/core/lsdtt_raster_reproject`
    acts on it — these are deliberately separate so suggestion logic can be used
    without triggering I/O
15. Reprojected rasters always receive a projection suffix (_utm / _polar / _lambert)
    so the provenance of every output file is unambiguous from its filename
16. `lsdtt_raster_algebra` lives in `src/cpu/lsdtt_raster_algebra.hpp/.cpp`; include from `src/cpu/`.
    Use `lsdtt::algebra::raster_subtract(a, b)` etc. for diff rasters in CLIs.
17. Every CLI uses `lsdtt_param_parser`. See the "CLI design rule" section above
    for the mandatory five-map pattern, `parse_all_parameters` call signature,
    and `print_help` / `print_parameters` conventions.
18. CLI parameters are grouped by workflow theme; each default and its `help_map`
    entry are declared on consecutive lines (default first, then help).
19. `parse_all_parameters` argument order is always: float, int, bool, string.
    Read results back via `pp.get_bool_parameters()` etc. into `bf`, `ff`, `ifs`, `sf`.
20. CLI parameter variables are intentionally not `const` because cross-parameter
    logic may adjust them after parsing.
21. Pure grid-geometry helpers belong in `src/utils/`; GDAL/OGR-backed vector
    and projection helpers belong in `src/core/`
22. Vector outputs should prefer typed attributes and single-file binary formats
    such as FlatGeobuf
23. `lsdtt_chi_ksn_inversion_cpu` requires Eigen3 and is compiled only when
    `LSDTT_HAS_EIGEN=1`. Do not include it in code that does not need inversion.
24. Hydraulic routing is first-party code in `lsdtt_hydraulic_graph_cpu`; there
    is no vendored graphflood directory to include from anywhere.
25. `src/utils/lsdtt_linear_algebra.hpp` and `lsdtt_thread_utils.hpp` are header-only
    with no lsdtt3 or GDAL dependencies — safe to include anywhere in the project.
26. GPU compute kernels use selective dispatch, not the raster backend. The
    `LSDTT_GPU_KERNELS` axis (kernels for heavy free functions, e.g. `direct_polyfit`)
    is separate from the future `LSDTT_GPU_ENABLED` raster-engine backend. The
    global toggle `set/get_compute_backend()` lives in core as pure host state;
    the dispatch decision and the `lsdtt3_cpu → lsdtt3_gpu` link live in `src/cpu`,
    so `src/core` stays CUDA-free. GPU kernels take POD pointers only (no
    `lsdtt_raster_cpu`/GDAL), accumulate in double for CPU parity, and always have
    a CPU fallback. See `docs/developer_notes/gpu_compute.md`.
