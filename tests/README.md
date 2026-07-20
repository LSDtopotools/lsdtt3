# lsdtt3 tests

CTest-driven test suite. Four tiers, selected by label:

| Label | Needs | What it covers |
|---|---|---|
| `smoke` | nothing | each CLI runs end to end on a tiny fixture |
| `rigorous` | nothing | exact output checks — checksums, sizes, feature counts, hashes |
| `failure` | nothing | bad parameters are rejected with the expected message |
| `unit` | nothing | C++ unit tests (parameter parser, cosmogenic physics, flow properties) |
| `cosmo` | nothing | cosmogenic-nuclide physics vs frozen CRONUS + numerics sentinel |
| `flow` | nothing | flow routing and depression-handling properties |
| `benchmark` | network + RichDEM | comparison against RichDEM on fetched real DEMs (opt-in) |

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
the validation driver for the CRN physics core. Plan:
`docs/cosmo_cli_testing_plan.md`. It runs offline against committed reference
data and prints one PASS/FAIL line per check.

> **The CRN numeric core is locked.** The analytical solutions exercised here
> encode years of manual verification. Any change to `lsdtt_crn_parameters.*`,
> `lsdtt_crn_column.*`, `lsdtt_crn_particle.*` or `lsdtt_cosmo_raster_cpu.*`, or
> to the reference values in `tests/expected/cosmo/`, requires Simon Mudd's
> explicit review **before** merging. See plan Sec3b.

What it checks:

- **Numerics sentinel** (plan Sec3b) — pinned values in
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
- **Plottable CSVs** — written to the test output directory for inspection with
  lsdviztools3.

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
network. Plan: `docs/fetch_raster_cli_plan.md` Phase 1.

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
handles secrets. Plan: `docs/fetch_raster_cli_plan.md` §1.

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

## Flow-routing property tests

`unit_flow_properties` (labels `unit;flow`) runs
`tests/cpp/test_flow_properties.cpp` — Phase 1 of
`docs/flow_routing_benchmark_plan.md`. It asserts invariants of depression
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

It also writes `fill_algorithm_characterisation.csv` (timings and divergence
between fill algorithms) into its output directory. That table is reported, not
asserted.

```bash
ctest --test-dir build -L flow --output-on-failure
```

## Flow-routing benchmark vs RichDEM (opt-in)

Phases 0 and 2 of `docs/flow_routing_benchmark_plan.md`. Scores lsdtt3's
depression handling against RichDEM on DEM tiles fetched at test time — **no DEM
is committed to the repo**.

Off by default (needs network + optional RichDEM), and skips cleanly via CTest's
exit code 77 when either is unavailable, so it can never turn CI red.

Python tooling uses **pixi** (`pixi.toml` at the repo root), matching
lsdviztools3. Install RichDEM from conda-forge via pixi, not pip: `pip install
richdem` gives 0.3.4, which pins `numpy<2` and downgrades numpy in whatever
environment it lands in. conda-forge ships 2.4.3 with numpy-2 builds.

```bash
pixi run fetch-test-dem --suite flow      # cache fixtures into tests/data/dem_cache/
pixi run -e benchmark python tests/scripts/run_flow_benchmark.py \
    --exe build/lsdtt-raster-preprocessing

# via CTest
cmake -S . -B build -DLSDTT_ENABLE_BENCHMARKS=ON
ctest --test-dir build -L benchmark --output-on-failure
```

Two terrain regimes are fetched, because they diverge very differently —
`alpine_small` (σ = 510 m) and `lakeland_small` (σ = 2.9 m). Epsilon-fill
disagreement with RichDEM is 0.5 % of cells on the former and 91 % on the
latter, so the low-relief fixture is not optional.

Tolerances in `run_flow_benchmark.py` are **not pinned yet** (Phase 3); the
driver reports metrics and only enforces limits under `--check`.

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

`tests/scripts/fetch_test_dem.py` is shared by any test suite needing real
topography, not just one. Fixtures are reprojected to a metric CRS **at fetch
time** rather than by the CLI under test — when two implementations are
compared, both must read a byte-identical raster or the differences mean
nothing.

```bash
pixi run fetch-test-dem --list                    # registry
pixi run fetch-test-dem --suite flow              # everything one suite needs
pixi run fetch-test-dem --fixture alpine --size small

# Ad-hoc bbox for the download-run-delete validation tier: transient DEM,
# only a small summary CSV is kept (docs/cosmo_cli_testing_plan.md Sec5b).
pixi run fetch-test-dem --bbox 10.3 46.4 10.5 46.6 \
    --dest /tmp/work/dem.tif --no-cache
```

Bounding boxes spanning several source tiles are mosaicked automatically, and
the UTM zone is derived from the bbox centre unless `--epsg` is given. To add a
fixture, add an entry to `REGISTRY` in that script rather than writing a new
fetcher.

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
