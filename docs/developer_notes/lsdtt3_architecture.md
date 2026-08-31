# lsdtt3 Architecture Document

A guide to the architecture of lsdtt3. 
Use for context when interacting with AI  during code development. 
Can also be used as a developer guide.

---

## Project identity

- **Package name:** lsdtt3
- **Repository:** https://github.com/LSDtopotools/lsdtt3
- **Copyright:** Simon M. Mudd
- **Licence:** GPL-3.0
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
├── cli/                                   — command-line executables (lsdtt-* dashed names)
│   ├── lsdtt-basin-extraction.cpp
│   ├── lsdtt-basin-stats.cpp
│   ├── lsdtt-basin-surface-metrics.cpp
│   ├── lsdtt-bearing-comparison.cpp
│   ├── lsdtt-channel-extraction.cpp
│   ├── lsdtt-channel-network-analysis.cpp
│   ├── lsdtt-chi-analysis.cpp             *(requires Eigen3)*
│   ├── lsdtt-cosmo-column-tool.cpp
│   ├── lsdtt-cosmo-landslides.cpp
│   ├── lsdtt-cosmo-prediction.cpp
│   ├── lsdtt-cosmo-tool.cpp
│   ├── lsdtt-cross-sections.cpp
│   ├── lsdtt-divide-migration.cpp
│   ├── lsdtt-drainage-analysis.cpp
│   ├── lsdtt-feature-finder.cpp
│   ├── lsdtt-fetch-raster.cpp
│   ├── lsdtt-flat-surface-finder.cpp
│   ├── lsdtt-graphflood.cpp
│   ├── lsdtt-internal-drainage.cpp
│   ├── lsdtt-lem.cpp
│   ├── lsdtt-linear-stats.cpp
│   ├── lsdtt-marsh-topography.cpp
│   ├── lsdtt-mountain-front.cpp
│   ├── lsdtt-raster-calculator.cpp
│   ├── lsdtt-raster-preprocessing.cpp
│   ├── lsdtt-raster-synthetic-tests.cpp
│   ├── lsdtt-ridgeline-analysis.cpp
│   ├── lsdtt-sediment-routing.cpp
│   ├── lsdtt-spectral-analysis.cpp        *(requires FFTW3)*
│   ├── lsdtt-stream-burning.cpp
│   ├── lsdtt-surface-metrics.cpp
│   ├── lsdtt-swath-profile.cpp
│   ├── lsdtt-synthetic-networks.cpp
│   ├── lsdtt-valley-centreline.cpp
│   ├── lsdtt-water-centreline.cpp
│   └── lsdtt-water-line-repair.cpp
├── docs/                                   — organised by audience
│   ├── index.md, overview.md               — user-facing, and the mkdocs site root
│   ├── install_lsdtt3.md                   — user-facing install guide
│   ├── testing_and_benchmarking.md
│   ├── cli_reference/                      — one parameter reference per CLI
│   ├── developer_notes/                    — architecture, building, packaging
│   │   ├── lsdtt3_architecture.md          — this file
│   │   ├── developer_onboarding.md
│   │   ├── building_from_source.md
│   │   ├── making_release_packages.md
│   │   ├── debian_packaging.md
│   │   ├── lsdtt_utils_guide.md
│   │   └── using_lsdtt_flowinfo.md
│   └── tex_technical/                      — LaTeX method write-ups (untracked)
├── graphflood/                             — third-party graph-based flood routing (Gailleton et al.)
│   ├── graphflood.hpp
│   ├── graphflood_v2.hpp
│   ├── fastflood.hpp
│   ├── fastflood_recorder.hpp
│   ├── graphflood_enums.hpp
│   └── graphflood_parts.hpp
├── python_helpers/
│   ├── dem_downloader_instructions.md
│   └── download_dem_tiles.py
├── src/
│   ├── core/
│   │   ├── lsdtt_types.hpp                  — all lightweight shared primitive types (see below)
│   │   ├── lsdtt_raster_meta.hpp            — metadata struct + CRS validation
│   │   ├── lsdtt_raster_base.hpp/.cpp       — abstract raster dispatch interface
│   │   ├── lsdtt_raster_result.hpp          — plain data-transfer struct for rasters returned by core-layer engines
│   │   ├── lsdtt_data_loader.hpp/.cpp       — GDAL/OGR I/O, RAII, COG output
│   │   ├── lsdtt_backend_factory.hpp/.cpp   — make_raster_engine(), gpu_available()
│   │   ├── lsdtt_param_parser.hpp/.cpp      — parameter file parser
│   │   ├── lsdtt_basin_base.hpp             — abstract basin interface
│   │   ├── lsdtt_junction_network_base.hpp  — abstract channel network interface + lsdtt_junction_network_options
│   │   ├── lsdtt_flow_base.hpp              — abstract D8 flow interface + lsdtt_flowinfo_options
│   │   ├── lsdtt_hydraulic_base.hpp         — abstract graph-based hydraulic routing interface + lsdtt_hydraulic_options
│   │   ├── lsdtt_crs_suggest.hpp/.cpp       — CRS inspection and projection recommendation
│   │   ├── lsdtt_raster_reproject.hpp/.cpp  — CRS reprojection
│   │   ├── lsdtt_raster_warp.hpp/.cpp       — GDAL warp helpers
│   │   ├── lsdtt_point_input.hpp/.cpp       — point/coordinate input handling
│   │   ├── lsdtt_stat_csv.hpp/.cpp          — CSV writers for QQ, CDF, and distribution-fit tables
│   │   └── lsdtt_vector_utils.hpp/.cpp      — GDAL/OGR-backed vector helpers
│   ├── cpu/                                   — 99 CPU-specific analysis modules (flow routing, surface metrics, channel analysis, LEM, cosmogenic nuclides, etc.)
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
│   │   ├── lsdtt_raster_algebra.hpp/.cpp                     — element-wise raster algebra
│   │   ├── lsdtt_raster_surface_metrics_cpu.hpp/.cpp         — hillshade, moving window, polyfit metrics, surface-shape (SMAP) classification
│   │   ├── lsdtt_surface_texture_cpu.hpp/.cpp                — MADk2, moving_window_iqr, NVD roughness metrics
│   │   ├── lsdtt_raster_surface_fitting_cpu.hpp/.cpp         — plane and polynomial surface fitting to masks
│   │   ├── lsdtt_raster_flow_preprocessing_cpu.hpp/.cpp      — fill, breach, nodata
│   │   ├── lsdtt_raster_nodata_filling_cpu.hpp/.cpp          — interior void filling
│   │   ├── lsdtt_raster_spectral_cpu.hpp/.cpp                — wavelet quadratic fit; spectral filters; fftw2D_masked_surface_recovery
│   │   ├── lsdtt_raster_opencv_cpu.hpp/.cpp                  — OpenCV-based raster processing
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
│   │   ├── lsdtt_stream_power_cpu.hpp/.cpp                   — implicit stream-power incision (D8 + D-inf)
│   │   ├── lsdtt_steady_state_cpu.hpp/.cpp                   — analytical chi-based fluvial steady state (shared core compute_fluvial_steady_state_on_tree, reused by the OCN generator)
│   │   ├── lsdtt_hillslope_transport_cpu.hpp/.cpp             — transient hillslope transport (cidre, roering, linear, adi)
│   │   ├── lsdtt_hillslope_implicit_cpu.hpp/.cpp              — implicit transient Roering solver (Eigen; backward-Euler/TR-BDF2 + damped Picard, channels as Dirichlet sinks)
│   │   ├── lsdtt_hillslope_steady_state_cpu.hpp/.cpp          — steady-state hillslope solver (SOR + pseudo-transient Cidre)
│   │   ├── lsdtt_hillslope_newton_cpu.hpp/.cpp                — steady/transient Newton Roering solver
│   │   ├── lsdtt_hillslope_cellular_cpu.hpp/.cpp              — soft-critical face transport, nonlocal routing, sparse particles
│   │   ├── lsdtt_hillslope_particle_io.hpp/.cpp               — particle CSV restart/output
│   │   ├── lsdtt_landscape_drape_cpu.hpp/.cpp                — Verlet cloth drape: upper envelope surface and incision depth
│   │   ├── lsdtt_landscape_equilibration_cpu.hpp/.cpp        — divide equilibration and landscape relaxation
│   │   ├── lsdtt_flat_surface_finder.hpp/.cpp                — flat surface identification (MADk2, residual IQR, NVD)
│   │   ├── lsdtt_flat_surface_plane_analysis.hpp/.cpp        — plane analysis on detected flat surfaces
│   │   ├── lsdtt_flat_surface_channel_classify_cpu.hpp/.cpp  — classify channels in flat surface areas
│   │   ├── lsdtt_hydraulic_graph_cpu.hpp/.cpp                — CPU hydraulic routing (lsdtt_hydraulic_base impl, Gailleton et al. 2024)
│   │   ├── lsdtt_chi_extraction_cpu.hpp/.cpp                 — chi coordinate extraction (full raster, outlet-seeded, inherited-base)
│   │   ├── lsdtt_chi_network_analysis_cpu.hpp/.cpp           — piecewise-linear fitting of chi channel networks
│   │   ├── lsdtt_chi_disorder_cpu.hpp/.cpp                   — chi disorder metric, jackknife, bootstrap, m/n sweep
│   │   ├── lsdtt_chi_data_maps_cpu.hpp/.cpp                  — chi data-map output (per-node points, basin outlines, polylines)
│   │   ├── lsdtt_chi_ksn_inversion_cpu.hpp/.cpp              — spatial k_sn/u* inversion with m/n model selection *(requires Eigen3)*
│   │   ├── lsdtt_channel_network_analysis_cpu.hpp/.cpp       — channel link bearings, junction angles
│   │   ├── lsdtt_channel_sources_wiener_cpu.hpp/.cpp         — channel source detection via Wiener filtering
│   │   ├── lsdtt_channel_scope_cpu.hpp/.cpp                  — channel scope analysis and classification
│   │   ├── lsdtt_strahler_tokunaga_analysis_cpu.hpp/.cpp     — Strahler channel aggregation, Tokunaga self-similarity, junction area-ratio distribution
│   │   ├── lsdtt_burn_streams_cpu.hpp/.cpp                   — burn polyline vector into DEM with monotone profile
│   │   ├── lsdtt_ridgeline_cpu.hpp/.cpp                      — ridgeline extraction and analysis
│   │   ├── lsdtt_valley_centreline_cpu.hpp/.cpp              — valley centreline extraction and tracing
│   │   ├── lsdtt_valley_cross_sections_cpu.hpp/.cpp          — valley cross-section analysis
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
│   └── utils/
│       ├── lsdtt_array2d.hpp                — header-only 2D array
│       ├── lsdtt_raster_geometry.hpp        — header-only raster geometry helpers
│       ├── lsdtt_polygon_geometry.hpp       — header-only: polygon area/perimeter/moments, convex hull, min-area rect, Feret, basin shape metrics
│       ├── lsdtt_linear_algebra.hpp         — header-only: sym3x3/sym2x2 eigenvalues, 3×3 inversion, 2D geometry helpers
│       ├── lsdtt_thread_utils.hpp           — header-only: resolve_thread_count() (OpenMP / hardware fallback)
│       ├── lsdtt_binary_array2d_tracing.hpp/.cpp — binary array tracing utilities
│       ├── lsdtt_distributions.hpp/.cpp     — parametric distribution fitting: normal, exponential, Weibull, Pareto (MLE + KS)
│       ├── lsdtt_statistics.hpp/.cpp        — descriptive statistics (median, MAD, stddev, etc.)
│       ├── lsdtt_histogram.hpp/.cpp         — histogram, PDF, CDF
│       ├── lsdtt_statistical_tests.hpp/.cpp — nonparametric tests, KS, circular statistics
│       ├── lsdtt_splines.hpp/.cpp           — 1D spline interpolation
│       └── lsdtt_regression.hpp/.cpp        — linear and polynomial regression helpers
├── CMakeLists.txt
├── diagnose_fftw.sh
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



**`lsdtt3_core` sources** (`src/core/`):
- `lsdtt_raster_base.cpp`
- `lsdtt_data_loader.cpp`
- `lsdtt_backend_factory.cpp`
- `lsdtt_param_parser.cpp`
- `lsdtt_crs_suggest.cpp`
- `lsdtt_raster_reproject.cpp`
- `lsdtt_raster_warp.cpp`
- `lsdtt_vector_utils.cpp`
- `lsdtt_point_input.cpp`
- `lsdtt_stat_csv.cpp`

**`lsdtt3_utils` sources** (`src/utils/`):
- `lsdtt_statistics.cpp`
- `lsdtt_histogram.cpp`
- `lsdtt_statistical_tests.cpp`
- `lsdtt_distributions.cpp`
- `lsdtt_splines.cpp`
- `lsdtt_regression.cpp`
- `lsdtt_binary_array2d_tracing.cpp`

**`lsdtt3_cpu` sources** (`src/cpu/`):
- `lsdtt_raster_cpu.cpp`
- `lsdtt_raster_algebra.cpp`
- `lsdtt_raster_surface_metrics_cpu.cpp`
- `lsdtt_surface_texture_cpu.cpp`
- `lsdtt_raster_surface_fitting_cpu.cpp`
- `lsdtt_raster_flow_preprocessing_cpu.cpp`
- `lsdtt_raster_nodata_filling_cpu.cpp`
- `lsdtt_raster_spectral_cpu.cpp`
- `lsdtt_raster_synthetic_surfaces_cpu.cpp`
- `lsdtt_synthetic_network_cpu.cpp`
- `lsdtt_synthetic_domain_cpu.cpp`
- `lsdtt_raster_feature_finder_cpu.cpp`
- `lsdtt_flowinfo_cpu.cpp`
- `lsdtt_flow_base.cpp`
- `lsdtt_junction_network_cpu.cpp`
- `lsdtt_junction_network_base.cpp`
- `lsdtt_depression_hierarchy_cpu.cpp`
- `lsdtt_fill_spill_merge_cpu.cpp`
- `lsdtt_flow_snapping_cpu.cpp`
- `lsdtt_basin_cpu.cpp`
- `lsdtt_junction_basin_cpu.cpp`
- `lsdtt_junction_basin_manager_cpu.cpp`
- `lsdtt_basin_outline_cpu.cpp`
- `lsdtt_flat_surface_finder.cpp`
- `lsdtt_flat_surface_plane_analysis.cpp`
- `lsdtt_vector_rasterise_cpu.cpp`
- `lsdtt_proximity_cpu.cpp`
- `lsdtt_swath_cpu.cpp`
- `lsdtt_streamlines_cpu.cpp`
- `lsdtt_lem_cpu.cpp`
- `lsdtt_stream_power_cpu.cpp`
- `lsdtt_steady_state_cpu.cpp`
- `lsdtt_hillslope_transport_cpu.cpp`
- `lsdtt_hillslope_steady_state_cpu.cpp`
- `lsdtt_flowgraph_cpu.cpp`
- `lsdtt_raster_vectorise_cpu.cpp`
- `lsdtt_raster_mosaic_cpu.cpp`
- `lsdtt_connectivity_cpu.cpp`
- `lsdtt_hydraulic_graph_cpu.cpp`
- `lsdtt_chi_extraction_cpu.cpp`
- `lsdtt_chi_network_analysis_cpu.cpp`
- `lsdtt_chi_disorder_cpu.cpp`
- `lsdtt_chi_data_maps_cpu.cpp`
- `lsdtt_channel_network_analysis_cpu.cpp`
- `lsdtt_strahler_tokunaga_analysis_cpu.cpp`
- `lsdtt_raster_skeleton_graph_cpu.cpp`
- `lsdtt_raster_binary_mask_vectorise_cpu.cpp`
- `lsdtt_raster_binary_mask_cpu.cpp`
- `lsdtt_burn_streams_cpu.cpp`
- `lsdtt_mountain_front_cpu.cpp`
- `lsdtt_landscape_drape_cpu.cpp`
- `lsdtt_cosmo_raster_cpu.cpp`
- `lsdtt_crn_parameters.cpp`
- `lsdtt_crn_column.cpp`
- `lsdtt_raster_opencv_cpu.cpp` *(when OpenCV is available)*
- `lsdtt_chi_ksn_inversion_cpu.cpp` *(when Eigen3 is available)*

`lsdtt3_cpu` depends on `lsdtt3_core` and `lsdtt3_utils`. All CLI executables
link against `lsdtt3_cpu` and `lsdtt3_core`. OpenMP, FFTW, OpenCV, and Eigen3
link to `lsdtt3_cpu`.

**Dependencies:**
- GDAL — required
- OpenMP — optional
- FFTW3 — optional (enables `lsdtt-spectral-analysis`)
- OpenCV — optional
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


---

## Naming conventions

- All snake_case, `lsdtt_` prefix on all classes and files
- No camelCase anywhere
- Opening braces always on their own line
- All parameters use snake_case with no spaces
- Output raster parameters: `write_*_raster`
- Algorithm selection parameters: `*_algorithm_choice`
- Output format parameter: `write_output_format` (`cog` / `gtiff` / `same_as_input`)

**CLI executables use dashes, not underscores:**
- Object/library: `lsdtt_flowinfo_cpu`, `lsdtt_raster_cpu`, `lsdtt_param_parser`
- CLI executable: `lsdtt-test-flowinfo`, `lsdtt-raster-preprocessing`

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

**`src/gpu/`** — GPU implementations. Currently `lsdtt_raster_gpu` placeholder.
Future: `lsdtt_flowinfo_gpu`.

**`src/utils/`** — header-only utilities with no backend dependency.

**Abstract interface pattern:**
```
src/core/lsdtt_raster_base    ← abstract interface for raster backends
src/cpu/lsdtt_raster_cpu      ← CPU implementation
src/gpu/lsdtt_raster_gpu      ← GPU implementation (placeholder)

src/core/lsdtt_flow_base      ← abstract interface for flow routing
src/cpu/lsdtt_flowinfo_cpu    ← CPU implementation
src/gpu/lsdtt_flowinfo_gpu    ← GPU implementation (future)
```

**Factory pattern:**
```
lsdtt_backend_factory  →  lsdtt_raster_cpu or lsdtt_raster_gpu
make_flow_base()       →  lsdtt_flowinfo_cpu or lsdtt_flowinfo_gpu (future)
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

## lsdtt_raster_meta

Plain struct in `src/core/lsdtt_raster_meta.hpp`.
Requires: `<cmath>`, `<stdexcept>`, `<string>`, `<gdal.h>`, `<ogr_spatialref.h>`

Key fields: `n_rows`, `n_cols`, `x_min`, `y_min`, `cell_size`, `no_data_value`,
`crs_wkt`, `epsg_code`, `data_type`.

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

- **`lsdtt-raster-synthetic-tests.cpp`**  
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
  layer (`basin_polygon_format` flatgeobuf/gpkg) with the per-basin hydromorphic,
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
  Identifies flat (low-texture) surfaces using MADk2, residual IQR, and/or NVD
  metrics. Outputs per-metric texture rasters, binary masks, and a combined mask.
  Thresholds are set by percentile; masks can be combined with AND or OR logic.

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

- **`lsdtt-cross-sections.cpp`**  
  Cross-section extraction along a polyline. Samples elevation, slope, and metrics 
  perpendicular to a baseline, outputs profiles as CSV or rasters.

- **`lsdtt-valley-centreline.cpp`**  
  Valley centreline extraction and tracing. Identifies valley networks from DEMs 
  using SMAP (surface morphological classification) and traces centrelines through 
  the identified valleys. Outputs centreline polylines.

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

- **`lsdtt-marsh-topography.cpp`**  
  Marsh and tidal platform analysis. Identifies marsh platforms and intertidal zones, 
  computes platform statistics, and extracts edges.

- **`lsdtt-bearing-comparison.cpp`**  
  Bearing and orientation analysis. Computes channel, ridge, and fracture bearings; 
  outputs rose diagrams and bearing statistics.

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
3. `static const` version, citation, and tool-name strings
4. `print_banner()` — fixed-width box to stdout. Required content: a
   brief summary of what the CLI does; `|| CLI version: <k_version>`
   (labelled as the CLI's own version, not the package version);
   the citation; a line pointing to `-v` / `--version` for the full
   lsdtt3 release version; and a line pointing to `-h` for parameter help
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
float_default_map["my_param_m"] = 10.0f;
help_map["my_param_m"] = {"float","10.0",
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
              k_version, k_citation);
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

## lsdtt_raster_base

File: `src/core/lsdtt_raster_base.hpp/.cpp`

Pure dispatch interface. Does NOT expose `get_data()` or `get_mask()` —
those are CPU-specific types that have no meaning on a GPU backend.
Only `get_meta()` is on the base (metadata is always host-side).

CPU analysis free functions take `const lsdtt_raster_cpu&` directly — they
never go through the base class for data access.

**Pure virtual:** `get_meta()`, `write_to_loader()`, `read_from_loader()`,
`smooth_gaussian()`, `smooth_bilateral()`

**CPU fallback virtuals** (base implements, GPU overrides when ready):
`compute_flow_direction()`, `compute_flow_accumulation()`, `compute_fft()`,
`apply_wavelet_filter()`

---

## lsdtt_raster_cpu

File: `src/cpu/lsdtt_raster_cpu.hpp/.cpp`

**Storage:** `lsdtt_array2d_f` and `lsdtt_array2d_u8` (each wraps `std::vector`).

**Accessors:**
- `get_data_array()` → `const lsdtt_array2d_f&` — preferred for analysis, zero-copy
- `get_mask_array()` → `const lsdtt_array2d_u8&` — zero-copy
- `get_data()` → `const std::vector<float>&` — zero-copy into array's vector
- `get_mask()` → `const std::vector<uint8_t>&` — zero-copy
- `get_meta()` → `const lsdtt_raster_meta&` — implements base pure virtual

**Constructors:**
- `lsdtt_raster_cpu(const lsdtt_data_loader&)` — read full raster from GDAL
- `lsdtt_raster_cpu(const lsdtt_raster_tile&)` — from loader tile
- `lsdtt_raster_cpu(lsdtt_array2d_f&&, lsdtt_array2d_u8&&, meta)` — preferred for results
- `lsdtt_raster_cpu(std::vector<float>&&, std::vector<uint8_t>&&, meta)` — backwards compat

---

## lsdtt_param_parser

File: `src/core/lsdtt_param_parser.hpp/.cpp`

Parameter file format: `key: value` (snake_case keys, `#` for comments).

**File I/O parameters:**
- `read_fname` — full filename **including extension** (e.g. `my_dem.tif`)
- `write_fname` — stem only, no extension (defaults to read stem)
- `read_path` — directory (default: parameter file directory)
- `write_path` — directory (default: read_path)
- `boundary_conditions` — string of characters: `o`=open, `n`=no-flux, `p`=periodic

`get_full_read_path()` returns `read_path + read_fname + read_ext`.
`get_full_write_path()` returns `write_path + write_fname` (no extension).

**Static helpers:** `stem(filename)`, `extension(filename)` (returns dot + ext),
`parse_args_deferred(argc, argv, tool_name, version, citation)` — handles
`-v`/`-c` (print and exit) and records `-h` variants in the returned
`lsdtt_cli_args` struct so the CLI can render help from its parameter
declarations via `print_terminal_help()` (see the CLI design rule above).
Every CLI uses this instead of a per-file
`parse_args` free function. Each CLI declares
`static const std::string k_tool_name = "lsdtt-foo";`
alongside `k_version` and `k_citation`.

**Version reporting (`-v` / `--version`):** prints both the CLI's own
`k_version` and the package-level version, e.g.
`lsdtt-chi-analysis version 0.2.0 (lsdtt3 0.2.0)`. The package version is
injected by CMake as the `LSDTT3_VERSION` compile definition on
`lsdtt3_core` (used only in `lsdtt_param_parser.cpp`). A build made exactly
on the release tag `v<PROJECT_VERSION>` reports the plain version; any
other build reports a development string derived from `git describe`, e.g.
`0.3.0-dev, git v0.2.0-14-g04460f4-dirty` — so distributed development
builds self-identify. Captured at configure time (rerun cmake to refresh).
Per-CLI `k_version` numbers are versioned independently of the package.

---

## lsdtt_data_loader

File: `src/core/lsdtt_data_loader.hpp/.cpp`

`lsdtt_data_loader` handles both raster and vector I/O through GDAL/OGR.

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

**Analysis code** never needs to know about `lsdtt_flow_base`, `lsdtt_flow_types`,
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
- `fill_pits(raster, min_slope=0.0001f, n_threads=0)` — Barnes et al. 2014
- `fill_pits_with_boundary(...)` — serial only
- `breach_depressions(raster, n_threads=0)` — Lindsay 2016
- `breach_depressions_least_cost(raster, max_dist=1000, max_cost=max_float,
  min_dist=true, flat_increment=0, fill_remainder=true, n_threads=0)` — Lindsay 2020

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

File: `cli/lsdtt-raster-preprocessing.cpp`  |  Version: 0.2.0

**Preprocessing pipeline (in order):**
1. **Projection** — detect and reproject geographic / Web Mercator rasters
2. **Threshold masking** — remove pixels outside elevation range
3. **Interior nodata filling** — fill enclosed voids
4. **Clip to rectangle** — largest valid rectangle
5. **Depression breaching** — Lindsay 2016 or Lindsay 2020 least-cost
6. **Pit filling** — Barnes 2014 Priority-Flood
7. **Hillshade** — Zevenbergen & Thorne 1987

Steps 2–7 always operate on the (possibly reprojected) raster. If reprojection
occurs the new file is reopened before step 2.

**Threshold masking semantics:**
`minimum_elevation` masks pixels **at or below** the threshold (i.e. `value <= min_elev`).
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
fill_algorithm_choice: barnes2014 / barnesetal2021
breach_algorithm_choice: lindsay2016 / least_cost
nodata_fill_method_choice: laplacian / iterative_weighted_average
hillshade_azimuth, hillshade_altitude, hillshade_z_factor
fill_min_slope, minimum_elevation, maximum_elevation
breach_lc_max_dist, breach_lc_max_cost, breach_lc_flat_increment
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

**Key parameters:**
```
read_fname, write_fname, read_path, write_path
flow_algorithm_choice: d8 / d_infinity / mfd   (default: d8)
mfd_exponent: 1.0
boundary_conditions: oooo
fill_before_routing, fill_algorithm_choice, breach_algorithm_choice
fill_min_slope, minimum_elevation, maximum_elevation
write_n_contributing_pixels_raster, write_flow_direction_raster
write_twi_raster
write_fill_raster, write_breached_raster
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
dt                          = 1.0e-3   # fixed timestep (s)
mannings_n                  = 0.033
min_slope                   = 1.0e-8
uniform_rainfall_rate       = 0.0      # m/s
receiver_width              = -1.0     # ≤ 0 = auto (cell_size)
adaptive_dt_enabled         = false
print_startup_message       = true
```

**Key virtual methods:**
- `run_one_step()`, `run_n_steps(n)`
- `set_uniform_rainfall_rate()`, `add_point_inflow(flat, m3/s)`
- `get_water_depth_raster()`, `get_water_surface_raster()`, `get_qin_raster()`, `get_qout_raster()`, `get_dhw_raster()`
- `refresh_routing()` — no-op for fixed-routing impls; intended extension point for dynamic routing

### lsdtt_hydraulic_graph_cpu

File: `src/cpu/lsdtt_hydraulic_graph_cpu.hpp/.cpp`

CPU implementation of `lsdtt_hydraulic_base`. Routing graph is copied from an
`lsdtt_flowgraph_cpu` at construction; the hydraulic solver evolves water depth
while routing remains fixed. Prints a startup message crediting Gailleton et al. (2024).

**Run modes:**
```cpp
enum class lsdtt_hydraulic_run_mode {
    transient_fixed_graph,
    stationary_fixed_graph,
    transient_dynamic_graph,    // future
    stationary_dynamic_graph    // future
};
```

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
- `bending_angle_signed/abs` — offset between inward receiver and bisector of donor opening (new diagnostic)

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

**`stream_burn_mode`:** `burn_in_place` (implemented) | `snap_to_network` (not yet implemented)

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

---

## lsdtt_raster_binary_mask_cpu / lsdtt_raster_skeleton_graph_cpu

**`lsdtt_raster_binary_mask_cpu`** (`src/cpu/lsdtt_raster_binary_mask_cpu.hpp/.cpp`)  
Generic binary mask operations (foreground = value > 0.5):
- `binary_morphological_dilation`, `erosion`, `opening`, `closing` (structuring element radius in metres)
- Reconstruction-based morphology, directional line-based closing
- Hole filling, size-based component filtering and labelling

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
   │  builds/loads DEM, conditions it, owns the time loop or steady-state branch
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
switches, boundary conditions). `lsdtt_lem_step_stats` reports `delta_max_elevation`,
`delta_mean_elevation`, `max_local_elevation_change`, `mean_erosion_rate`,
`max_erosion_rate`, and optionally the `erosion_rate` field.

The fluvial flow object is built with `allow_pits_as_outlets = true` because the
transient LEM surface may contain pits that must not stall routing.

**Channel mask.** `compute_lem_channel_mask(dem, opts)` returns the Dirichlet
channel set the hillslope step would use for a surface — the explicit mask when
`opts.explicit_channel_mask` is supplied, else the pixels whose D-inf drainage area
reaches `incision_threshold_area`. It shares its implementation with the hillslope
step, so the mask a run *writes* is exactly the one it *used*. The CLI exposes it as
`write_channel_mask_raster` (with an independent `channel_mask_threshold_area`, so a
run that erodes fluvially everywhere can still emit a network), which lets a fluvial
run freeze its network and feed it back as `channel_mask_fname` to a hillslopes-only
run. `hold_channel_mask_fixed` additionally adds the masked pixels to `boundary_flats`,
so uplift cannot lift the frozen network away from the hillslopes relaxing toward it.

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

### lsdtt_stream_power_cpu (transient incision)

File: `src/cpu/lsdtt_stream_power_cpu.hpp/.cpp`

Implicit stream-power incision, `U = K·A^m·S^n` (Braun & Willett 2013 generalised
to D-inf per Anand et al. 2020). Both functions sweep the topological order in
**reverse** (outlets first) so a receiver elevation is already updated when its
donor is solved.

- `compute_d8_stream_power_incision(dem, flowinfo, K, m, n, dt, threshold_area, divide_area_by_dx, slope_clamp)`
  — single-receiver D8; serial reverse-topological sweep.
- `compute_dinf_stream_power_incision(..., flowgraph, ..., n_threads)`
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

Channel pixels (drainage area ≥ `incision_threshold_area`) act as hillslope sinks
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

## graphflood/ (third-party)

The `graphflood/` directory contains third-party header-only graph-based flood routing code by Gailleton et al. It is not a CMake library target — it is included directly by `lsdtt_hydraulic_graph_cpu.cpp`. Do not add it as a CMake target or include it from other lsdtt3 modules.

---

## Known compile errors fixed (for reference)

| File | Issue | Fix |
|------|-------|-----|
| `lsdtt-raster-preprocessing.cpp` | `const lsdtt_raster_base&` to cpu fn | Changed to `const lsdtt_raster_cpu&` |

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
24. `graphflood/` headers are consumed only by `lsdtt_hydraulic_graph_cpu.cpp`.
    Do not include them from other lsdtt3 sources.
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
