# lsdtt3 tests

What each test does and how to run it. For **methodology, benchmark results, and
the documented differences between lsdtt3 and RichDEM**, see
`docs/developer_notes/testing_and_benchmarking.md` — read that before concluding a mismatch
against another tool is a bug.

CTest-driven test suite. Every test carries labels on four independent axes;
select with `-L` / `-LE`, which take **regular expressions**, so anchor them
(`-L '^hermetic$'`): `-L network` also matches `synthetic_networks`.

**1. Runtime tier — exactly one per test: what it needs to start.**

| Label | Registered | Needs |
|---|---|---|
| `tier_native` | always (with `BUILD_TESTING`) | the C++ build only |
| `tier_python_stdlib` | `-DLSDTT_ENABLE_PYTHON_TESTS=ON` | Python ≥ 3.10, standard library only; the `needs_gdal_tools` ones also need GDAL's command-line programs |
| `tier_python_geo` | also `-DLSDTT_ENABLE_PYTHON_GEO_TESTS=ON` | Python plus the GDAL bindings, NumPy, Rasterio, GeoPandas, Shapely |
| `tier_pixi` | `-DLSDTT_ENABLE_BENCHMARKS=ON` | `pixi run` into the `benchmark` environment |

**2. Environment flags — zero or more.**

| Label | Meaning |
|---|---|
| `hermetic` | native, no skip path, no network, no host data: the reproducible floor |
| `needs_network` | downloads data (registered only with `-DLSDTT_ENABLE_NETWORK_TESTS=ON`; skips when offline) |
| `needs_gdal_tools` | runs `gdalinfo` / `ogrinfo` / `gdal_translate` (Ubuntu: `gdal-bin`, which `libgdal-dev` does not install) |

**3. Component.** `lem_component`: every `lsdtt-lem` CLI test and the unit
tests of the engines it drives.

**4. Cost.** `heavy`: too slow for the inner loop.

**Kind and domain** labels say what is tested and where:

| Label | What it covers |
|---|---|
| `smoke` | each CLI runs end to end on a tiny fixture |
| `rigorous` | exact output checks — checksums, sizes, feature counts, hashes |
| `failure` | bad parameters are rejected with the expected message |
| `unit` | C++ unit tests (parameter parser, cosmogenic physics, flow properties) |
| `divide_diagnostics` | divide-diagnostics compute-layer unit tests |
| `landscape_equilibration` | landscape equilibration cycle unit test |
| `cosmo` | cosmogenic-nuclide physics vs frozen CRONUS + numerics sentinel |
| `flow` | flow routing and depression-handling properties |
| `cosmo_realworld` | inverted erosion rates vs frozen published/OCTOPUS rates on real DEMs (network) |
| `benchmark` | comparison against RichDEM, TopoToolbox, Whitebox on fetched real DEMs (opt-in) |
| `performance` | in-process timings (opt-in, reported not asserted) |
| `configure_contract` | configure the whole project in a scratch directory and check what it reports for a disabled or missing tier |

Configure prints what it registered, per tier, and which options left tests
out and how many.

### What each tier needs installed

- **Native**: the C++ build dependencies only.
- **Python**: Python ≥ 3.10; GDAL's command-line programs for the
  `needs_gdal_tools` tests (`sudo apt install gdal-bin`).
- **Python geospatial**: easiest from the Pixi `test` environment — configure
  *outside* Pixi and pass its interpreter (see `pixi.toml`):

  ```bash
  pixi install -e test
  cmake -S . -B build -DLSDTT_ENABLE_PYTHON_TESTS=ON -DLSDTT_ENABLE_PYTHON_GEO_TESTS=ON \
        -DPython3_EXECUTABLE="$(pixi run -e test python -c 'import sys; print(sys.executable)')"
  pixi run -e test ctest --test-dir build --output-on-failure
  ```

- **Network**: `-DLSDTT_ENABLE_NETWORK_TESTS=ON` and a connection.
- **Benchmarks**: `-DLSDTT_ENABLE_BENCHMARKS=ON` and the Pixi `benchmark`
  environment.

A missing interpreter, program or package is a configure error that names it.

### Execution contract

- **No test reads the terminal.** CTest passes the terminal it was started
  from to every test as stdin, so anything that prompts would hang an
  interactive run while passing in CI. Every test runs with
  `LSDTT_NONINTERACTIVE=1` (the credential lookup then never prompts), and the
  Python scripts start every program through `tests/scripts/lsdtt_test_env.py`,
  which gives it `/dev/null` as stdin.
- **Every test has a timeout.** Unless a test sets its own: `failure_*` 60 s,
  `heavy` 3600 s, the rest 300 s; benchmarks and performance tests are left to
  their options.
- **OpenMP tests do not fight under `-j`.** Every test except the serial,
  benchmark and performance ones runs with `OMP_WAIT_POLICY=passive`.
- **Pixi isolation.** Run inside a Pixi environment, the launcher strips that
  environment's GDAL/PROJ data and library variables from lsdtt3 programs built
  against a different GDAL (`unit_test_env_isolation` checks it).
- **The audit.** `unit_test_inventory_audit` checks every test's labels,
  environment and timeout against the command it actually runs, and that test
  scripts launch programs only through `lsdtt_test_env`.

## Running the suite

```bash
ctest --test-dir build -L '^hermetic$' --output-on-failure   # native, offline, reproducible
ctest --test-dir build -LE '^heavy$'                          # everything registered except the slow tests
ctest --test-dir build -L '^cosmo$' -LE '^heavy$'             # one domain
ctest --test-dir build -LE '^lem_component$'                  # everything except the LEM
ctest --test-dir build -L '^lem_component$' -LE '^heavy$'     # the LEM, cheap part
ctest --test-dir build -L '^cosmo_realworld$'                 # real-landscape cosmo validation
                                                              # (needs -DLSDTT_ENABLE_NETWORK_TESTS=ON)
```

For long runs, the logging wrapper keeps a complete log in `build/test-logs/`
even when a run is killed (CTest's own `LastTest.log.tmp` is overwritten by the
next `ctest` call of any kind):

```bash
cmake -P build/run_tests.cmake -LE '^heavy$' -j4
```

**When tests fail**, build the report and send it with the bug report:

```bash
cmake --build build --target test_report
```

It writes one plain-text file in `build/test-logs/`: build options and
versions, the machine, the GDAL/PROJ actually loaded, the Python environment,
the run's arguments and counts, and each failing test's output and logs. The
home directory appears as `~` and only a short fixed list of environment
variables is included; read it before sending it.

The offline run is dominated by `unit_hillslope_implicit`, whose wall time swings
widely under WSL2 host load (68 s to 219 s observed on the same machine). If that
becomes a nuisance, gate it by label rather than loosening its tolerances.

`cosmo_realworld` is registered only with `-DLSDTT_ENABLE_NETWORK_TESTS=ON`
(since 2026-09-18; before that a bare `ctest` pulled it in). It downloads real
DEMs and dominated the wall time when it ran by default (~745 s of ~865 s, one
test alone ~540 s). Those tests are validation against published erosion rates, not smoke tests: each
one's DEM extent is tied to the study area of its samples, so they cannot be made
faster by shrinking the DEM without destroying what they measure. Run them when
touching cosmogenic code, not on every iteration.

## CLI smoke and rigorous tests

- `smoke_raster_preprocessing`: runs `lsdtt-raster-preprocessing` on a tiny DEM and checks `tiny_fill.tif` exists and is non-empty.
- `smoke_raster_preprocessing_tif`: runs `lsdtt-raster-preprocessing` on `tiny30_dem.tif` and checks `tiny30_fill.tif` exists and is non-empty.
- `smoke_channel_extraction`: runs `lsdtt-channel-extraction` on the same DEM and checks `tiny_channel_network_lines.fgb` exists and is non-empty.
- `smoke_basin_extraction`: runs `lsdtt-basin-extraction` and checks `tiny_junction_basin_outlines.fgb` exists and is non-empty.
- `smoke_nodata_fill`: runs `lsdtt-raster-preprocessing` on a tiny DEM with an interior nodata pixel and checks `tiny_nodata_nodata_filled.tif` exists and reports `STATISTICS_VALID_PERCENT=100` via `gdalinfo -stats`.
- `smoke_nodata_edge_unfilled`: runs `lsdtt-raster-preprocessing` on a tiny DEM with an edge nodata pixel and checks `STATISTICS_VALID_PERCENT=96` (edge nodata remains unfilled).
- `smoke_drainage_analysis`: runs `lsdtt-drainage-analysis` and checks `tiny_d8_drainage_area.tif` exists.
- `smoke_chi_analysis`: runs `lsdtt-chi-analysis` and checks `tiny_chi_fit_nodes_mn_0p500.csv` exists.
- `smoke_raster_calculator_difference`: runs `lsdtt-raster-calculator` and checks `tiny_difference.tif` exists.

More rigorous checks currently included:

- `rigorous_channel_lines_feature_count`: requires channel lines output to contain at least 1 feature.
- `rigorous_basin_outlines_feature_count`: requires basin outlines output to contain at least 1 feature.
- `rigorous_tif_fill_raster_size`: requires the `tiny30` fill raster to be exactly `32x39` pixels.
- `rigorous_raster_calculator_difference_checksum`: checks the GDAL checksum of `tiny_difference.tif`.
- `rigorous_drainage_d8_checksum`: checks the GDAL checksum of `tiny_d8_drainage_area.tif`.
- `rigorous_chi_nodes_sha256`: checks SHA256 of `tiny_chi_fit_nodes_mn_0p500.csv`.
- `failure_channel_extraction_invalid_threshold`: expected-failure test for invalid `threshold_contributing_pixels`.
- `failure_raster_calculator_missing_second`: expected-failure test for missing `second_raster_fname`.

## Divide migration

`unit_divide_migration` (labels `unit;divide_migration`) runs
`tests/cpp/test_divide_migration.cpp` — 30 checks over
`src/cpu/lsdtt_divide_migration_cpu.cpp`. The headline property is **N6**: a head
is never matched to a divide between two *other* units, even when that divide is
80× closer in space. That is what the per-unit spatial index exists to guarantee.
Also covered: mid-segment projection (every fixture puts the answer strictly
between two vertices, so a nearest-vertex bug cannot hide), a NaN — not infinite
— gradient at zero distance, negative relief passed through unchanged, nodata
sampling, and per-fragment chi re-zeroing under a common base level.

The same binary writes the **CLI fixture** into its output directory (a CTest
`FIXTURES_SETUP`): a synthetic three-valley DEM plus outlet layers in FGB and
CSV. Generated in C++, never committed.

- `smoke_divide_migration_no_clusters` — one unit per outlet: 3 basins, **2**
  divide segments (1|2 and 2|3; the outer valleys are not adjacent so 1|3 never
  exists), 69 channel heads, all matched.
- `smoke_divide_migration_clusters` — the same DEM with outlets 1 and 2 grouped:
  the 1-2 divide becomes internal and the topology collapses to **1** segment.
- `failure_divide_migration_single_unit` — all outlets in one cluster is rejected
  with `outlets=3, units=1` rather than silently producing no divides.

## Distributed hydrology: the S-HYPE column

`column_model: hype` of `lsdtt-distributed-hydrology` is tested against
**HYPE 5.36.0's own output**. The files are committed under `tests/expected/hype/`
and were written by `tests/scripts/generate_hype_oracle.py`. That script is run
by hand: CTest only reads its files and never downloads or compiles HYPE.

```bash
ctest -R "hype|subbasin_partition|distributed_hydrology_forcing" -LE heavy
```

**Library suites** (C++ drivers):

- `unit_hype_parameters`, `unit_hype_column`, `unit_hype_routing`,
  `unit_hype_model`: compare against the oracle within a tolerance about five
  times the measured single-vs-double precision floor of HYPE itself (1.08e-4 mm).
- `unit_subbasin_partition`, `unit_hype_raster_setup`: the set-up from rasters.
- `unit_hype_thread_invariance`: the model must be bitwise identical at 1 and 8
  threads.

**CLI tests:**

- `fixture_hype_cli` generates the input set (a 60×60 Y-network DEM, land use,
  soil, parameter tables and 120 days of forcing with snow). Nothing is committed.
- The `smoke_distributed_hydrology_hype_*` tests cover:
  - the three outlet modes;
  - `expanding_square` and `downstream` snapping of off-channel gauges;
  - nodata holes in the soil raster;
  - a run shorter than the forcing;
  - the daily default for `dt_seconds`.
- The `failure_distributed_hydrology_hype_*` tests assert the FATAL message for
  each rejected input.

The timing harness (`lsdtt-test-hype-performance`) is always built, but
it is registered as the test `performance_hype_model` only under
`LSDTT_ENABLE_PERF_TESTS`. What each suite asserts, and the measured numbers,
are in §10q of `docs/developer_notes/testing_and_benchmarking.md`.

## lsdtt-lem

`smoke_lem_strip_fixture` generates the hillslope strip through
`lsdtt-synthetic-surfaces`; `smoke_lem_roering_implicit_strip` then runs
`lsdtt-lem` with `hillslope_model: roering_implicit` on it and asserts the
**converged steady state** — peak 4.9429 m and mean 3.2719 m, both two-sided.

This is the integration test §8.4 of `docs/developer_notes/testing_and_benchmarking.md` argues
for. The semi-implicit fast path is only taken when the field is well below `Sc`,
which is a regime the standalone driver's steady-state tests never *finish* in —
so the omega bug there biased coupled LEM steady states ~40 % high while every
unit check passed. Reinstating that bug moves this test's peak to 6.862 m
(+38.8 %), far outside the band. The strip starts **flat**, so "the solver did
nothing" leaves it at 0 m and cannot be mistaken for the answer. ~1.5 s for both.

### Spatially varying erodibility and the model grid

`fixture_lem_erodibility_cli` writes every fixture these tests read, with
`lsdtt-test-lem-erodibility --write-fixtures`, into
`build/tests/lem_erodibility_cli_fixture/`: a 40 × 20 DEM, K rasters on and off
its grid, placement sources whose pixels hold their own indices
(`1 + 1000 r + c`), the expected results of each placement (computed from the
placement rules, not by calling the code under test), channel masks and forcing
CSVs. The `smoke_`/`failure_lem_erodibility_*`, `..._lem_grid_*` and
`..._lem_channel_mask_*` runs read them through `--read-path`, and the
`regression_lem_*` tests compare written rasters with `--compare-rasters` (dims,
cell size, origin, CRS, mask, data bitwise) or `--compare-raster-values` (dims,
mask, data). The library tests are `unit_lem_erodibility` (fast) and
`rigorous_lem_erodibility_analytic` (the analytic steady states); see section 8c
of `docs/developer_notes/testing_and_benchmarking.md`.

### Lithology to erodibility

`fixture_lem_lithology_erodibility_cli` writes the fixtures for
`lsdtt-lithology-erodibility` (`lsdtt-test-lithology-erodibility --write-fixtures`)
into `build/tests/lem_lithology_erodibility_cli_fixture/`: GeoPackage maps (the
three-unit map, one with a hole, one in EPSG:4326, one with an unknown
lithology), a code raster, a 20 m DEM, property tables (valid, K = 0, no m), a
crosswalk, and the expected K rasters. `smoke_lem_lithology_erodibility_chain`
then runs `lsdtt-lem` on the K raster the own-grid test wrote.

### k_sn by lithologic unit" in tests/README.md with this. -->
### k_sn by lithologic unit

`fixture_chi_lithology_cli` (`lsdtt-test-lithology-channel-stats --write-fixtures`)
writes `build/tests/chi_lithology_cli_fixture/`: a steady-state DEM with two units
(k_sn about 20 west, 5 east), the unit polygons, a crosswalk naming the units
(`lithks_crosswalk.csv`) and one missing a value (`lithks_crosswalk_gap.csv`), the
unit polygons and a code raster 200 km east of the DEM in the same CRS
(`lithks_units_far.gpkg`, `lithks_codes_far.tif`, for the no-overlap refusals), and a
two-unit code raster beside a copy of `tiny_dem.asc`. `rigorous_chi_lithology_observation_only`
reuses the hash of `rigorous_chi_nodes_sha256`: if that hash is ever re-pinned, re-pin both.

`rigorous_chi_lithology_ksn_ratio` is not a CLI run: it runs the driver's
`--check-ksn-by-lithology` mode on the summary that `smoke_chi_lithology_ksn` wrote
(FIXTURES_REQUIRED on that test), and checks the sandstone/shale ratio of the reach
median k_sn against 4. The same check mode is exercised inside
`unit_chi_lithology_channel_stats` on a summary written in process.
`rigorous_chi_lithology_bootstrap_per_mn` compares the 0.45 tables of
`smoke_chi_lithology_ksn` and `smoke_chi_lithology_multi_mn` (m/n 0.4 and 0.45) byte
for byte, with `tests/cmake/compare_output_sets.cmake`.

### Macrostrat lithology fetch

`unit_fetch_lithology` and `smoke_fetch_lithology_offline` read one committed,
synthetic Mapbox Vector Tile, `data/macrostrat_tiles/10/540/364.pbf`, through
`macrostrat_tile_dir`, so they need no network. `data/macrostrat_tiles/README.md`
lists its units and the `ogr2ogr` command that regenerates it from
`units_source.geojson`. `network_fetch_lithology_cli` downloads from Macrostrat
and is registered only with `-DLSDTT_ENABLE_NETWORK_TESTS=ON`.

### Channel-link bearing and junction-angle statistics

`lsdtt-test-channel-bearing-statistics` (`tests/cpp/test_channel_bearing_statistics.cpp`)
runs one mode per test, all labelled `unit;network` and native (no Python):
`unit_channel_bearing_statistics_circular` (the analytic circular Wasserstein and
Kuiper cases, including the four the retired `lsdtt-bearing-comparison` smoke tests
pinned), `_summaries`, `_comparisons`, `_junction_angles`, `_trellis` (a synthetic
trellis network built in memory: N-S trunks, E-W tributaries, known bearings and
junction angles), `_grouping` (region purity, order pooling, order-k and nested
outlet basins) and `_isotropic_null` (a seeded diamond-square surface; the seed and
the measured value are recorded in the source).

`fixture_channel_bearing_statistics` (`--write-fixtures`) writes
`build/tests/channel_bearing_statistics_fixture/`: the trellis DEM and its
transpose, a north/south region code raster on the trellis grid, an off-grid region
raster with no CRS, a code raster and a one-polygon GeoPackage (field `unit`) 200 km
east of the trellis in the same CRS (`cbs_region_far.tif`, `cbs_region_far.gpkg`, for
the no-overlap refusals), and a copy of `data/channel_bearing_statistics/ns_trunk_bearings.csv`.
The `*_channel_network_analysis_bearing_*` CLI tests (Python harness) read it.

### Channel-network statistics by lithologic unit

`unit_lithology_channel_network_stats` (`tests/cpp/test_lithology_network_stats.cpp`,
`unit;network;lithology`, native) builds its own copy of the trellis in memory and
gives it lithologic units on the grid in process: a contact parallel to the trunks
(no channel crosses it), one across them, and one with an unmapped strip. The
`*_channel_network_analysis_lithology*` CLI tests read
`build/tests/channel_bearing_statistics_fixture/` and use its north/south code
raster as the lithology map. `rigorous_channel_network_analysis_lithology_observation_only`
compares the statistics outputs of a run with an unused lithology map against
`smoke_channel_network_analysis_bearing_statistics_order` byte for byte, with
`tests/cmake/compare_output_sets.cmake`.

## Cosmogenic nuclide tests

`unit_cosmo_physics` (labels `unit;cosmo`) runs `tests/cpp/test_cosmo_physics.cpp`,
the validation driver for the CRN physics core (methodology:
`docs/developer_notes/testing_and_benchmarking.md`). It runs offline against committed reference
data and prints one PASS/FAIL line per check.

> **The CRN numeric core is locked.** The analytical solutions exercised here
> encode years of manual verification. Any change to `lsdtt_crn_parameters.*`,
> `lsdtt_crn_column.*`, `lsdtt_crn_particle.*` or `lsdtt_cosmo_raster_cpu.*`, or
> to the reference values in `tests/expected/cosmo/`, requires Simon Mudd's
> explicit review **before** merging. See docs/developer_notes/testing_and_benchmarking.md

What it checks:

- **Numerics sentinel** — pinned values in
  `tests/expected/cosmo/numerics_sentinel.csv` that fail on any drift in the
  numeric core. This is the guard that makes the lock enforceable.
- **Equation fidelity (E1–E9)** — the library is checked against independent
  re-implementations written from the papers rather than by calling the library,
  so agreement is meaningful. Covers the CAIRN steady-state and depth-integrated
  solutions (Mudd et al. 2016 Eq. 12/14), the step-change point solution and its
  analytic depth integral, limiting behaviour (saturation, fast erosion), and
  flux-weighted basin accumulation.
- **Exact parameter pins** — Be-10 (19) and Al-26 (7) standardisation factors,
  and Braucher/Borchers production fractions summing to 1.
- **CRONUS comparison (T-C1…T-C3)** — against a frozen CRONUS v3 reference in
  `tests/expected/cosmo/cronus_v3_reference.csv`. Pressure is submitted directly
  so atmosphere-model differences are excluded from the comparison by
  construction.
- **Geomagnetic / time-varying production (G1–G4, §6c(A))** — the Muscheler/
  Knudsen VDM time-integrated Stone/Lal correction (`stone2000sp_paleomag`) and
  the ERA40 atmosphere. G1 static/degenerate limits (a constant-VDM record must
  reproduce the time-independent factor exactly; `!loaded`/`age≤0` fallbacks),
  G2 physics invariants (high-latitude insensitivity, equatorial sign and
  monotonicity, ERA40 grid interpolation), G4 golden trajectories against the
  committed `data/crn_parameters.nc` (Knudsen VDM), and the equatorial LSDTT-1
  geomag ON/OFF delta. These run against the bundled `crn_parameters.nc`,
  **auto-discovered** relative to the executable (`src/core/lsdtt_exe_path.hpp`).
- **LSDn scaling (L0–L6, L10, L11)** — the Lifton–Sato–Dunai (2014) production
  scaling that runs additively alongside Stone. L0 pins the packed grid axes and
  the Pavón-Carrasco/Lifton record constants by exact equality; **L1 is a
  SEPARATE numerics sentinel**, `tests/expected/cosmo/numerics_sentinel_lsdn.csv`
  (360 values, own gated flag `--regenerate-lsdn-sentinel-i-have-simons-approval`
  — the Stone `numerics_sentinel.csv` must stay byte-identical); L2 reproduces
  the frozen v3 grid reference `lsdn_grid_reference.csv` (generated offline by
  `python_helpers/extract_lsdn_reference.py --emit-reference` from the untracked
  v3 `.mat` files) including the axis-orientation value 4.880884; L3 checks
  LSD2014 Eq. 2 against an independent re-implementation, its negative leak above
  80° and the genuine longitude dependence of the moving paleopole; L4 the clamps
  and the flat sub-1 GV columns; L5 the degenerate cases — a missing LSDn grid is
  a **hard failure, never a silent Stone fallback** (checked against a pre-LSDn
  `crn_parameters.nc` in `tests/data`); L6 the physics invariants and the `S_tot`
  interaction with `scale_F_values`; L10 that the precomputed time-averaged
  surface equals the direct per-pixel loop (both paths invoked explicitly via
  `force_build`, including an age past the end of the long record so the
  final-step stretch is exercised); L11 a measured horizontal-resolution sweep
  emitting `lsdn_surface_resolution_sweep.csv`. The Stone half of L11 — and
  therefore the shipped `h_step_deg` — is still open.
- **Plottable CSVs** — `production_vs_depth.csv`, `geomag_scaling.csv`,
  `geomag_on_off_LSDTT1.csv`, `lsdn_vs_stone.csv`,
  `lsdn_surface_vs_direct.csv`, `lsdn_surface_resolution_sweep.csv` and the
  earlier ones, written to the test output directory for inspection with
  lsdviztools3.

The three cosmo CLIs have offline integration tests (labels `smoke;cli;cosmo`
and `failure;cli;cosmo`): `smoke_cosmo_{tool,prediction,column_tool}_help` check
that every parameter declaration renders, and `failure_cosmo_tool_geographic_crs`
checks that a non-projected (EPSG:4326) DEM is rejected with a clean FATAL rather
than a crash or a silent wrong-units run.

### Real-landscape validation tier (download → run → delete)

`unit_cosmo_validation_harness` (labels `unit;cosmo`) runs the **offline**
self-test of `tests/scripts/run_cosmo_validation.py` — the "download → run →
delete" harness for published-basin validation
(`tests/expected/cosmo/published/README.md`). The self-test proves the cleanup
trap removes the temporary working directory on a *failed* run (with a mock
workdir, no network), that `--keep` is honoured, that the manifest/sample-table
schema guards trip, and that summary-append leaves the frozen skeleton
byte-unchanged.

`cosmo_realworld_dethier2014` (label `cosmo_realworld`, registered only with
`-DLSDTT_ENABLE_NETWORK_TESTS=ON`) is the real thing: it fetches a DEM **keyless** via `lsdtt-fetch-raster`
(Copernicus GLO-30 / Planetary Computer, no API key), clips it, runs
`lsdtt-cosmo-tool`, compares to the tracked published rates, appends the frozen
summary CSV, then deletes the temporary workdir including on failure. It needs
network and is run on request (`ctest -L '^cosmo_realworld$'` in a build configured
with the network option); it exits 77 (skip)
when the network/DEM source is unavailable or the sample table is still a
template, so it never turns CI red. **No DEM is ever committed** — only the small
manifest, sample table, and frozen summary under `tests/expected/cosmo/published/`.

`cosmo_realworld_lsdn_pixel_vs_lumped` (label `cosmo_realworld`, registered only
with `-DLSDTT_ENABLE_NETWORK_TESTS=ON`) is the **generation** half of the per-pixel-vs-lumped LSDn result.
It fetches all three published datasets' DEMs keyless, runs `lsdtt-cosmo-tool`
under `scaling_scheme: lsdn` at three sample ages plus `stone`, and records each
basin's mean production scaling over the exact upslope node set the inversion
uses. Those numbers are **frozen** into
`tests/expected/cosmo/lsdn_pixel_vs_lumped_perpixel.csv`; the cheap comparison
against the lumped LSDn kernel then runs in the **default** `unit;cosmo` suite
inside `unit_cosmo_physics`, which emits `lsdn_pixel_vs_lumped.csv` and asserts
the frozen bands. Expensive step manual, cheap check always. Like its siblings it
is **check-only** by default and rewrites the tracked artefact only under an
explicit `--freeze`.

Refreshing the CRONUS reference is a **manual, reviewable** step and is never
done by CI:

```bash
python tests/scripts/fetch_cronus_reference.py \
    --samples build/tests/cosmo_test_output/cronus_fixture_samples.csv
```

It sends a single batched request, caches the raw XML alongside the CSV, and the
regenerated CSV should be committed with a note explaining why it moved.

## Native DEM fetch tests

`unit_raster_fetch` (labels `unit;fetch`) runs `tests/cpp/test_raster_fetch.cpp`
in offline mode: 31 checks over `src/core/lsdtt_raster_fetch.*`, needing no
network. See `docs/developer_notes/testing_and_benchmarking.md`.

- UTM zone derivation, including zone boundaries and the ±180 antimeridian
  (an off-by-one silently reprojects into the neighbouring zone rather than
  failing outright)
- Copernicus tile naming, including bboxes crossing 1° boundaries in either or
  both axes, and southern/western hemisphere naming
- Option validation — a malformed bbox is an **error**, never a "skip"
- **Graceful degradation** — an unobtainable tile yields `skipped=true`, so
  callers can skip rather than fail
- **No credential leakage** — asserts nothing credential-shaped appears in any
  loggable result field. It passes trivially today (Phase 1 has no keyed
  sources); its job is to fail the day a key is routed through these structs.

`network_raster_fetch` actually downloads and is **off by default**:

```bash
cmake -S . -B build -DLSDTT_ENABLE_NETWORK_TESTS=ON
ctest --test-dir build -L '^fetch$' --output-on-failure
```

It fetches a 128² window, checks the crop, EPSG, nodata and that the elevations
are plausible alpine values, and exits 77 (skip) when the source is unreachable.

## Credential handling tests

`unit_secret` (labels `unit;fetch`) runs `tests/cpp/test_secret.cpp` — 24
offline checks over `src/core/lsdtt_secret.*`, the only place in lsdtt3 that
handles secrets. See `docs/developer_notes/testing_and_benchmarking.md`.

These are security properties, not conveniences:

- **Masking removes every occurrence** of a key from a string, leaves the rest
  of a URL intact for diagnosis, and treats an empty secret as a no-op (an empty
  needle would otherwise match at every position and shred the text)
- **Resolution order** is exactly `OT_API_KEY` → `OT_API_KEY_FILE` → TTY prompt,
  and a named-but-unreadable key file is reported rather than silently skipped
- **No prompt when stdin is not a terminal** — prompting under CTest would hang
  the run forever instead of failing cleanly
- **A world-readable key file warns** (and the warning does not contain the key)
- **The secret never appears** in `source`, `warning` or `message` — the fields
  intended for logging

The `chmod 600` check is skipped where `chmod` is a no-op (NTFS under WSL); the
test detects that rather than reporting a filesystem limitation as a defect. On
a real POSIX filesystem all 24 checks assert.

Why this matters concretely: the OpenTopography API **echoes the submitted key
back in its own error message**. Without masking, every mistyped or expired key
would be printed to the terminal and into any redirected log.

### Bounds from a vector file

`unit_raster_fetch` also covers deriving the fetch area from an `.fgb`/`.gpkg`
containing a polygon or two corner points (fixtures in `tests/data/bounds/`,
committed because they are a few kB each).

The headline property: **both diagonals give the identical box.** Two points at
lower-left/upper-right and at lower-right/upper-left produce the same
`min`/`max`, so there is no orientation to detect and no way to hold it wrong.
A test asserts the two arrangements agree to 1e-12.

Also covered: reprojection from UTM (the geometry is transformed *before* the
envelope is taken, since transforming envelope corners under-covers the true
extent), rejection of three points, of mixed geometry types, and of missing
files.

The `lsdtt-fetch-raster` CLI adds four more offline tests plus one networked
one: `smoke_fetch_raster_help`, `failure_fetch_raster_missing_bbox` (a missing
bbox edge must be a hard error — a silently zeroed box fetches the Gulf of
Guinea), `failure_fetch_raster_bad_source`, and
`rigorous_fetch_raster_no_credential_parameters`.

That last one scans the CLI's help output and fails if any parameter is named
like a credential (`api_key`, `password`, `token`, `secret`). Parameter files
get committed and `print_parameters()` echoes every parameter at the end of a
run, so **a credential must never be a parameter**. The test exists to stop a
future change quietly introducing one.

## Flow-routing performance

`performance_flow_routines` runs `tests/cpp/test_flow_performance.cpp`, behind
`-DLSDTT_ENABLE_PERF_TESTS=ON`. It is **reported, never asserted** — wall-clock
thresholds in CI are flaky on shared runners — and every row of
`flow_performance.csv` carries a git SHA so runs are comparable across commits.

```bash
cmake -S . -B build -DLSDTT_ENABLE_PERF_TESTS=ON
ctest --test-dir build -L '^performance$' --output-on-failure

# or directly, including a real raster for a like-for-like RichDEM comparison
build/tests/lsdtt-test-flow-performance out_dir --full \
    --raster tests/data/dem_cache/alpine_small.tif
```

Timings are **in-process**: the raster is already in memory and nothing is
written to disk. This matters — the `lsdtt_seconds` column of
`benchmark_results.csv` is whole-process CLI wall time (startup, parameter
parsing, GDAL open and write) and is **not** comparable with RichDEM's
in-routine numbers. Timed properly, lsdtt3 is 1.8–6.6× *faster* than RichDEM
across every routine.

## Flow-routing property tests

`unit_flow_properties` (labels `unit;flow`) runs
`tests/cpp/test_flow_properties.cpp` — Phase 1 of
`docs/developer_notes/testing_and_benchmarking.md`. It asserts invariants of depression
handling and flow routing on deterministic synthetic DEMs generated in-process,
so it needs no network, no downloaded DEM, and no RichDEM:

- **P1** filled DEM has no interior local minimum (each paired with a
  non-triviality check that the fixture *had* depressions to begin with)
- **P2** fill never lowers the surface
- **P3** breaching raises only pit/flat cells, and never above the lowest
  original neighbour (`breach_depressions`
  pit-flattens before carving, so a plain "never raises" assertion is wrong)
- **P4** total D8 accumulation reaching outlets equals the valid-cell area
- **P5** accumulation increases downstream
- **P6** filling an already-filled DEM is a no-op
- **P7** D-inf / MFD receiver weights sum to 1; D-inf has ≤ 2 receivers per cell
- **P8** D8 accumulation is **exact on a planar ramp**, and every cell routes due
  south rather than diagonally — an analytic oracle needing no reference
  implementation. This is the check that distinguishes true steepest descent
  (drop ÷ distance) from "lowest neighbour by elevation"; RichDEM fails it, which
  is why accumulation is validated here rather than by comparison.

It also writes `fill_algorithm_characterisation.csv` (timings and divergence
between fill algorithms) into its output directory. That table is reported, not
asserted.

```bash
ctest --test-dir build -L '^flow$' --output-on-failure
```

## Flow-routing benchmark vs RichDEM (opt-in)

Phases 0 and 2 of `docs/developer_notes/testing_and_benchmarking.md`. Scores lsdtt3's
depression handling against RichDEM on DEM tiles fetched at test time — **no DEM
is committed to the repo**.

Off by default (needs network + optional RichDEM), and skips cleanly via CTest's
exit code 77 when either is unavailable, so it can never turn CI red.

Python tooling uses **pixi** (`pixi.toml` at the repo root), matching
lsdviztools3. Install RichDEM from conda-forge via pixi, not pip: `pip install
richdem` gives 0.3.4, which pins `numpy<2` and downgrades numpy in whatever
environment it lands in. conda-forge ships 2.4.3 with numpy-2 builds.

```bash
# Fixtures are fetched automatically by lsdtt-fetch-raster on first run.
pixi run -e benchmark python tests/scripts/run_flow_benchmark.py \
    --exe build/lsdtt-raster-preprocessing

# via CTest
cmake -S . -B build -DLSDTT_ENABLE_BENCHMARKS=ON
ctest --test-dir build -L '^benchmark$' --output-on-failure
```

Comparisons run against **two** external references: RichDEM for depression
handling, and **TopoToolbox for D8 accumulation** — RichDEM's D8 is not steepest
descent and is recorded but not asserted. lsdtt3 and TopoToolbox agree *exactly*
on high-relief D8 accumulation (rank correlation 1.000000).

Two terrain regimes are fetched, because they diverge very differently —
`alpine_small` (σ = 510 m) and `lakeland_small` (σ = 2.9 m). Epsilon-fill
disagreement with RichDEM is 0.5 % of cells on the former and 91 % on the
latter, so the low-relief fixture is not optional.

Tolerances are **pinned per fixture** in `run_flow_benchmark.py` and enforced
under `--check` (which the CTest wiring passes). They are set just above
observed divergence, so they assert something real — Barnes 2014 fill is pinned
to **exactly 0.0** against RichDEM, since it is bit-exact and both routines are
deterministic.

Two guards make that safe:

- Tolerances are per fixture, because breach divergence is fixture-dependent
  (alpine breach measured 4.54 m on one fixture build and 0.00 m on another).
- `EXPECTED_FIXTURE_SHA` records the fixture each tolerance set was measured
  against. A different GDAL build produces a different raster, and the run
  **skips** rather than applying numbers that no longer mean anything.

The benchmark also asserts the **depression-free property** on the real
fixtures — no interior minima after filling — with or without `--check`. That
check needs no reference, so it still has teeth if both implementations were
wrong in the same way.

## Flat resolution: Barnes et al. (2014b)

The opt-in `fill_flat_resolution_choice` options, `barnes2014b_gradient` and
`barnes2014b_float_step` (the default `inline_gradient` is covered by the
existing fill tests). The method is in
`docs/design_and_algorithms/flow_routing.md`, and the test record is in
`docs/developer_notes/testing_and_benchmarking.md` (section 10r).

```bash
ctest --test-dir build -R 'flat_resolution' -LE '^heavy$' --output-on-failure
```

**C++ drivers** (`tests/cpp/test_flat_resolution.cpp`, hermetic, in-process):

- `unit_flat_resolution` (`unit;flow`): the mask, the step, both resolvers, the
  rounds and the final pass. Groups F1-F27 and F29 (there is no F28); F1 is the
  bitwise `step_up_float32` check, F5 pins the valley profile, F6 the exact
  gradient, F21-F23 the stranding cases (interior rim, two flats a few ULPs
  apart, cascades) and F29 the single-pass output frozen as hashes.
- `unit_flat_resolution_real_terrain` (`unit;flow`): the same executable in
  `--real-terrain <crop>` mode, on `tests/data/flat_resolution_lakeland_crop_dem.tif`.
  The single pass strands cells there and the router refuses them. The defaults
  must route.

**Fill CLI, `lsdtt-raster-preprocessing`:**

- `smoke_fill_flat_resolution_float_step` and `smoke_fill_flat_resolution_gradient`:
  both options run and write (plumbing only; `tiny_dem.asc` has no flats).
- `failure_fill_flat_resolution_unknown_choice`,
  `failure_fill_flat_resolution_gradient_zero_slope` and
  `failure_fill_flat_resolution_fill_spill_merge`: the exact refusal messages.
- `rigorous_fill_flat_resolution_valley_float_step` and
  `rigorous_fill_flat_resolution_valley_gradient`: exact fill-difference ranges
  and report lines on the valley fixture (`flat_valley_dem.asc`).
- `rigorous_fill_flat_resolution_knob_rounds`: the one-ULP knob
  (`flat_valley_knob_dem.asc`) needs two rounds, and the result routes.

**Routing CLIs with a barnes option:**

- `smoke_chi_analysis_flat_resolution_gradient`,
  `smoke_channel_extraction_flat_resolution_gradient` and
  `smoke_drainage_analysis_flat_resolution_float_step`: one run each, asserting
  a line only the barnes path prints.
- `smoke_chi_analysis_flat_resolution_filled_float_step`: `raster_is_filled: true`
  takes the resolve-only path.
- `failure_chi_analysis_flat_resolution_filled_inline`: the control; the default
  on the same filled input must fail with "could not assign a receiver".
- `failure_chi_analysis_flat_resolution_unfilled_input`: a valley with a pit
  (`flat_valley_pit_dem.asc`) and `raster_is_filled: true` is refused with the
  no-outlet message.
- `smoke_basin_extraction_flat_resolution_gradient`,
  `smoke_basin_extraction_flat_resolution_filled_float_step` and
  `failure_basin_extraction_flat_resolution_filled_inline`.
- `smoke_basin_stats_flat_resolution_float_step` and
  `failure_basin_stats_flat_resolution_no_fill_note` (the
  `fill_before_routing: false` NOTE, on `flat_valley_dem.asc`).
- `smoke_basin_surface_metrics_flat_resolution_gradient` (outlets in
  `flat_valley_outlets.csv`).
- `smoke_channel_network_analysis_flat_resolution_float_step`.
- `smoke_cross_sections_flat_resolution_gradient`.
- `smoke_bankfull_width_flat_resolution_inline_warning`,
  `smoke_bankfull_width_flat_resolution_gradient_warning`,
  `smoke_bankfull_width_flat_resolution_float_step_warning` and
  `smoke_bankfull_width_flat_resolution_filled_float_step_warning`: the bankfull
  WARNING, one per branch (source points in `flat_valley_source.csv`).
- `smoke_distributed_hydrology_flat_resolution_hype_gradient` and
  `smoke_distributed_hydrology_flat_resolution_topmodel_float_step`.
- `smoke_fluvial_surfaces_flat_resolution_float_step`,
  `failure_fluvial_surfaces_flat_resolution_no_fill_note` and
  `failure_fluvial_surfaces_flat_resolution_gradient_zero_slope` (the latter names
  `fill_min_slope_for_routing` in its message).
- `smoke_ridgeline_flat_resolution_gradient`.
- `smoke_valley_centreline_flat_resolution_float_step` (mask in
  `flat_valley_mask.asc`).

**Real terrain, through `lsdtt-drainage-analysis`** (the committed crop, the
command that failed before the rounds were added):

- `rigorous_drainage_analysis_flat_resolution_real_terrain_float_step` and
  `rigorous_drainage_analysis_flat_resolution_real_terrain_gradient`: exit 0, a
  D8 raster, `Rim cells without descent: 0`, the pinned rounds and final pass,
  and no "could not assign a receiver".
- `smoke_drainage_analysis_flat_resolution_real_terrain_inline_gradient`: the
  control on the same crop.
- `smoke_flat_resolution_real_terrain_routing`: a network-tier script,
  `tests/scripts/check_flat_resolution_real_terrain_routing.py`. It runs the three
  options on the flow-benchmark alpine and lakeland small and medium fixtures (fetched
  when not cached) and on the 500 and 1000 m benchmark rasters when cached. Every
  run must exit 0 with no cell without descent. Registered only with
  `-DLSDTT_ENABLE_NETWORK_TESTS=ON`; it skips (77) offline.

**RichDEM oracles** (tier `pixi`, `-DLSDTT_ENABLE_BENCHMARKS=ON`):

- `benchmark_flat_resolution_richdem` (O1,
  `tests/scripts/check_flat_resolution_vs_richdem.py`): the first round of
  `barnes2014b_float_step` against RichDEM `ResolveFlats` on float32, bit for bit
  except on near-zero-guarded flats and on the cells the later rounds changed. It
  runs on synthetic fixtures and the committed crop.
- `benchmark_flats_float_step_alpine_small` and
  `benchmark_flats_float_step_lakeland_small` (O2, `run_flow_benchmark.py`
  comparison `flats_float_step`): the same comparison on the fetched small DEMs.

## Run locally

```bash
cmake -S . -B build -DLSDTT_ENABLE_PYTHON_TESTS=ON -DLSDTT_ENABLE_PYTHON_GEO_TESTS=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Without the two Python options only the native tests are registered; most
`smoke`, `rigorous` and `failure` tests are Python-driven. Run by label:

```bash
ctest --test-dir build -L '^smoke$' --output-on-failure
ctest --test-dir build -L '^rigorous$' --output-on-failure
ctest --test-dir build -L '^failure$' --output-on-failure
```

Label conventions:

- `smoke_*` tests are labeled `smoke;cli`
- `rigorous_*` tests are labeled `rigorous;cli`
- `failure_*` tests are labeled `failure;cli`

## Add a new smoke test

1. Add a param template in `tests/params/`.
2. Add expected output filename(s).
3. Register it in `tests/CMakeLists.txt` with `add_lsdtt_cli_smoke_test(...)`,
   **above** the DEPENDENCY TIERS block near the end of the file (a test
   registered after it gets no tier label, which the audit reports). The
   helper registers it only when `LSDTT_ENABLE_PYTHON_TESTS` is ON and marks
   it as a Python test; set its properties with `lsdtt_set_tests_properties`,
   which skips a test that is switched off.
4. A test registered with `add_test` and Python directly must call
   `lsdtt_mark_test_runtime(python <name>)` right after it, inside
   `if(LSDTT_ENABLE_PYTHON_TESTS)`; one whose script imports third-party
   packages also goes in `LSDTT_PYTHON_GEO_TESTS`.
5. Scripts start programs through `tests/scripts/lsdtt_test_env.py`
   (`lsdtt_test_env.run(...)`), never `subprocess` directly.

The runner script is `tests/scripts/run_cli_smoke.py`.

## Fixtures and test data

Committed, tiny, offline:

- `tests/data/tiny_dem.asc` and variants (`_nodata`, `_edge_nodata`, `_plus1`) —
  the DEMs behind the smoke and failure tests.
- `tests/data/tiny30_dem.tif` — 30x30 GeoTIFF clipped from
  `../data/scotland/scotland.tif`.
- `tests/data/flat_valley_*` — the 11 × 24 valley at 10 m used by the
  flat-resolution tests. `flat_valley_dem.asc` is the valley after a
  `min_slope: 0` fill (a flat floor at 100 m, walls at 110 m, outlet at row 5,
  column 0); `flat_valley_knob_dem.asc` raises one floor cell by one float32
  ULP; `flat_valley_pit_dem.asc` is unfilled, with one floor cell at 99.5 m.
  `flat_valley_outlets.csv`, `flat_valley_source.csv` and `flat_valley_mask.asc`
  give the outlet, source and mask inputs of the routing CLIs that take them.
- `tests/data/flat_resolution_lakeland_crop_dem.tif` — 128x128 float32 GeoTIFF
  (DEFLATE, predictor 3), 30 m, WGS 84 / UTM 15N, no nodata cells; Canadian
  Shield lakeland, northwestern Ontario (origin 564930 E, 5604600 N). It is
  pixel window xoff 256, yoff 192 of the benchmark suite's `lakeland_500` rung
  (`run_benchmark_suite.py`: the centred 500 px crop of the `lakeland` 1000 px
  source fetched from Copernicus GLO-30 on AWS by `lsdtt-fetch-raster`), i.e.
  window xoff 506, yoff 442 of that source; cut with no resampling.
  sha256 `95864d0f08bcb2fce2ba72ee9b948aff4bc358ed43ee58a501cf367bf14cc4fc`.
  Lake surfaces a few float32 ULPs apart, on which a single Barnes et al.
  (2014b) flat-resolution pass leaves cells without a lower neighbour; used by
  `unit_flat_resolution_real_terrain` and the
  `*_drainage_analysis_flat_resolution_real_terrain_*` CLI tests.
  Attribution: produced using Copernicus WorldDEM-30 © DLR e.V. 2010-2014 and
  © Airbus Defence and Space GmbH 2014-2018 provided under COPERNICUS by the
  European Union and ESA; all rights reserved.
- `tests/expected/cosmo/` — frozen CRONUS v3 reference, raw XML, and the
  numerics sentinel.

Fetched on demand, never committed:

- `tests/data/dem_cache/` (gitignored) — real DEM tiles cut from Copernicus
  GLO-30 on AWS Open Data. Public, keyless, and read as Cloud-Optimised GeoTIFFs
  so only the requested window is transferred.

Fixtures are fetched by **`lsdtt-fetch-raster`**, the lsdtt3 CLI itself — there
is no Python fetcher any more. `run_flow_benchmark.py` calls it automatically
the first time a fixture is needed, so there is no separate step to remember.
The fixture registry (centre point, EPSG, pixel size) lives in
`run_flow_benchmark.py`.

Fixtures are reprojected to a metric CRS **at fetch time** rather than by the
CLI under test — when two implementations are compared, both must read a
byte-identical raster or the differences mean nothing.

To fetch one by hand:

```bash
build/lsdtt-fetch-raster <dir> <param file>     # see docs/cli_reference/cli_reference_lsdtt_fetch_raster.md
```

Named fixtures can be SHA256-pinned (`--print-hashes` prints a paste-ready
block). A hash mismatch is a hard failure, not a silent re-download: it means
the upstream tile or the fetch parameters changed, which voids every downstream
expectation.

The pins are only meaningful **together with the committed `pixi.lock`**.
Fetching is deterministic for a given GDAL build — repeated `--force` runs
reproduce the hash exactly — but warp output is GDAL-version sensitive: a pixi
re-solve that moved GDAL changed the alpine fixture by up to 2.6 m per cell
(while leaving lakeland byte-identical). Regenerating `pixi.lock` therefore
means re-pinning, and the GDAL version should be recorded with the hashes.
