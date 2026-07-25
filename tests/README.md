# lsdtt3 tests

What each test does and how to run it. For **methodology, benchmark results, and
the documented differences between lsdtt3 and RichDEM**, see
`docs/testing_and_benchmarking.md` — read that before concluding a mismatch
against another tool is a bug.

CTest-driven test suite, selected by label:

| Label | Needs | What it covers |
|---|---|---|
| `smoke` | nothing | each CLI runs end to end on a tiny fixture |
| `rigorous` | nothing | exact output checks — checksums, sizes, feature counts, hashes |
| `failure` | nothing | bad parameters are rejected with the expected message |
| `unit` | nothing | C++ unit tests (parameter parser, cosmogenic physics, flow properties) |
| `cosmo` | nothing | cosmogenic-nuclide physics vs frozen CRONUS + numerics sentinel |
| `flow` | nothing | flow routing and depression-handling properties |
| `benchmark` | network + RichDEM | comparison against RichDEM on fetched real DEMs (opt-in) |
| `performance` | nothing | in-process timings of the flow routines (opt-in, reported not asserted) |

Everything except `benchmark` runs offline with committed fixtures. `benchmark`
is off unless configured with `-DLSDTT_ENABLE_BENCHMARKS=ON`, and skips rather
than fails when its dependencies are absent.

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
- `failure_channel_extraction_invalid_threshold`: expected-failure test for invalid `source_pixels_threshold`.
- `failure_raster_calculator_missing_second`: expected-failure test for missing `second_raster_fname`.

## Cosmogenic nuclide tests

`unit_cosmo_physics` (labels `unit;cosmo`) runs `tests/cpp/test_cosmo_physics.cpp`,
the validation driver for the CRN physics core (methodology:
`docs/testing_and_benchmarking.md`). It runs offline against committed reference
data and prints one PASS/FAIL line per check.

> **The CRN numeric core is locked.** The analytical solutions exercised here
> encode years of manual verification. Any change to `lsdtt_crn_parameters.*`,
> `lsdtt_crn_column.*`, `lsdtt_crn_particle.*` or `lsdtt_cosmo_raster_cpu.*`, or
> to the reference values in `tests/expected/cosmo/`, requires Simon Mudd's
> explicit review **before** merging. See docs/testing_and_benchmarking.md

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
delete" harness for published-basin validation (plan §5b/§6c(B),
`tests/expected/cosmo/published/README.md`). The self-test proves the cleanup
trap removes the temporary working directory on a *failed* run (with a mock
workdir, no network), that `--keep` is honoured, that the manifest/sample-table
schema guards trip, and that summary-append leaves the frozen skeleton
byte-unchanged.

`cosmo_realworld_dethier2014` (label `cosmo_realworld`, **not** in the default
suite) is the real thing: it fetches a DEM **keyless** via `lsdtt-fetch-raster`
(Copernicus GLO-30 / Planetary Computer, no API key), clips it, runs
`lsdtt-cosmo-tool`, compares to the tracked published rates, appends the frozen
summary CSV, then deletes the temporary workdir including on failure. It needs
network and is invoked manually (`ctest -L cosmo_realworld`); it exits 77 (skip)
when the network/DEM source is unavailable or the sample table is still a
template, so it never turns CI red. **No DEM is ever committed** — only the small
manifest, sample table, and frozen summary under `tests/expected/cosmo/published/`.

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
network. Plan: `docs/testing_and_benchmarking.md` Phase 1.

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
ctest --test-dir build -L fetch --output-on-failure
```

It fetches a 128² window, checks the crop, EPSG, nodata and that the elevations
are plausible alpine values, and exits 77 (skip) when the source is unreachable.

## Credential handling tests

`unit_secret` (labels `unit;fetch`) runs `tests/cpp/test_secret.cpp` — 24
offline checks over `src/core/lsdtt_secret.*`, the only place in lsdtt3 that
handles secrets. Plan: `docs/testing_and_benchmarking.md` §1.

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
ctest --test-dir build -L performance --output-on-failure

# or directly, including a real raster for a like-for-like RichDEM comparison
build/tests/lsdtt-test-flow-performance out_dir --full \
    --raster tests/data/dem_cache/alpine_small.tif
```

Timings are **in-process**: the raster is already in memory and nothing is
written to disk. This matters — the `lsdtt_seconds` column of
`benchmark_results.csv` is whole-process CLI wall time (startup, parameter
parsing, GDAL open and write) and is **not** comparable with RichDEM's
in-routine numbers. Timed properly, lsdtt3 is 1.8–6.6× *faster* than RichDEM
across every routine (plan §7e).

## Flow-routing property tests

`unit_flow_properties` (labels `unit;flow`) runs
`tests/cpp/test_flow_properties.cpp` — Phase 1 of
`docs/testing_and_benchmarking.md`. It asserts invariants of depression
handling and flow routing on deterministic synthetic DEMs generated in-process,
so it needs no network, no downloaded DEM, and no RichDEM:

- **P1** filled DEM has no interior local minimum (each paired with a
  non-triviality check that the fixture *had* depressions to begin with)
- **P2** fill never lowers the surface
- **P3** breaching raises only pit/flat cells, and never above the lowest
  original neighbour (see plan §7a Finding 1 — `breach_depressions`
  pit-flattens before carving, so a plain "never raises" assertion is wrong)
- **P4** total D8 accumulation reaching outlets equals the valid-cell area
- **P5** accumulation increases downstream
- **P6** filling an already-filled DEM is a no-op
- **P7** D-inf / MFD receiver weights sum to 1; D-inf has ≤ 2 receivers per cell
- **P8** D8 accumulation is **exact on a planar ramp**, and every cell routes due
  south rather than diagonally — an analytic oracle needing no reference
  implementation. This is the check that distinguishes true steepest descent
  (drop ÷ distance) from "lowest neighbour by elevation"; RichDEM fails it, which
  is why accumulation is validated here rather than by comparison (plan §7d).

It also writes `fill_algorithm_characterisation.csv` (timings and divergence
between fill algorithms) into its output directory. That table is reported, not
asserted.

```bash
ctest --test-dir build -L flow --output-on-failure
```

## Flow-routing benchmark vs RichDEM (opt-in)

Phases 0 and 2 of `docs/testing_and_benchmarking.md`. Scores lsdtt3's
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
ctest --test-dir build -L benchmark --output-on-failure
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

## Run locally

```bash
cmake -S . -B build
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Run by label:

```bash
ctest --test-dir build -L smoke --output-on-failure
ctest --test-dir build -L rigorous --output-on-failure
ctest --test-dir build -L failure --output-on-failure
```

Label conventions:

- `smoke_*` tests are labeled `smoke;cli`
- `rigorous_*` tests are labeled `rigorous;cli`
- `failure_*` tests are labeled `failure;cli`

## Add a new smoke test

1. Add a param template in `tests/params/`.
2. Add expected output filename(s).
3. Register it in `tests/CMakeLists.txt` with `add_lsdtt_cli_smoke_test(...)`.

The runner script is `tests/scripts/run_cli_smoke.py`.

## Fixtures and test data

Committed, tiny, offline:

- `tests/data/tiny_dem.asc` and variants (`_nodata`, `_edge_nodata`, `_plus1`) —
  the DEMs behind the smoke and failure tests.
- `tests/data/tiny30_dem.tif` — 30x30 GeoTIFF clipped from
  `../data/scotland/scotland.tif`.
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
