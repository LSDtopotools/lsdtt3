# Benchmarking and testing in lsdtt3

How lsdtt3 is tested, what it is benchmarked against, and — importantly — the
places where lsdtt3 and the external reference implementations **deliberately
differ**. If you are comparing lsdtt3 output against another tool and the
numbers do not match, read §5 before concluding anything is wrong.

- Practical "what each test does and how to run it": `tests/README.md`
- This document: methodology, external references, measured findings, and how
  to extend the suite.

---

## 1. Principles

A few rules that the suite is built on. They are worth stating because they
explain why some obvious-looking tests are absent and some unusual ones exist.

**Property tests beat comparisons.** A property test asserts something that must
be true of any correct implementation — a filled DEM has no interior minimum,
accumulation increases downstream, weights sum to 1. It keeps its value even
when the reference implementation is wrong, which (see §5.1) is not a
hypothetical concern.

**An analytic oracle beats both.** Where a case has a known closed-form answer —
a planar ramp, a cone — assert the exact answer. No reference required, and the
test says something precise about correctness rather than about agreement.

**A missing optional dependency must skip, never fail.** Network access,
RichDEM, and downloaded fixtures are all optional. When absent the driver exits
with CTest's skip code (77) and the test reports as skipped. CI must never go
red because a machine is offline.

**A malformed request is an error, not a skip.** The converse of the above: only
an unavailable *dependency* may skip. A bad bounding box or a missing parameter
must fail loudly, or a typo silently registers as "skipped" and the suite passes
while testing nothing.

**Reported is not asserted.** Performance numbers and comparisons against a
reference we do not fully trust are recorded to CSV but never assert. A test
that cannot fail is noise; a wall-clock threshold on a shared runner is flaky.
Both are worse than an honest record.

**Tolerances are pinned just above observed divergence.** A tolerance set to a
comfortable round number asserts nothing. Where a result is exact, it is pinned
to exactly zero.

---

## 2. Test tiers

Everything except `benchmark` and `performance` runs offline against committed
fixtures.

| Label | Needs | Covers |
|---|---|---|
| `smoke` | — | each CLI runs end to end on a tiny fixture |
| `rigorous` | — | exact output checks: checksums, sizes, feature counts, hashes |
| `failure` | — | bad parameters are rejected with the expected message |
| `unit` | — | C++ unit tests |
| `cosmo` | — | cosmogenic physics vs frozen CRONUS + numerics sentinel |
| `cosmo_realworld` | network | published-basin runs against frozen summaries |
| `surface_metrics` | — | surface-metrics validation (slope, curvature, SMAP, MADk2) |
| `flow` | — | flow routing and depression-handling properties |
| `fetch` | — | native DEM fetching; offline logic by default |
| `network` | network | tests that actually download |
| `benchmark` | network + RichDEM/TopoToolbox | comparison against external references on real DEMs |
| `performance` | — | in-process timings, reported and never asserted |

```sh
cmake -S . -B build && cmake --build build -j
ctest --test-dir build --output-on-failure          # default: offline suite

ctest --test-dir build -L unit --output-on-failure
ctest --test-dir build -L flow --output-on-failure
```

### Compartments: run what your change can affect

Running the whole suite after every edit is wasteful. Most tests cannot be
affected by a given change, and a few are heavy enough to dominate the clock.
So every test carries **two orthogonal kinds of label**: a *kind*
(`unit` / `smoke` / `rigorous` / `failure`) and a **domain** saying which part of
the codebase it exercises. The expensive ones additionally carry `heavy`.

| Domain | Tests | Covers |
|---|---|---|
| `cosmo` | 29 | cosmogenic physics, the three cosmo CLIs, CRONUS comparisons |
| `fetch` | 24 | DEM fetching, offline logic by default |
| `network` | 14 | chi, channels, drainage, basins, flow properties |
| `divide_migration` | 11 | divide migration |
| `raster` | 11 | preprocessing, filling, nodata, the raster calculator |
| `landslide` | 8 | the triggered-landslide module |
| `surface_metrics` | 5 | slope, curvature, SMAP, MADk2 |
| `hillslope` | 3 | the implicit hillslope solver, wherever it is driven from |
| `hydraulic` | 3 | hydraulic routing |
| `lem` | 2 | landscape evolution |
| `step_change` | 2 | the analytical step-change raster generator |
| `core` | 2 | parameter parsing, credential handling |

The intended inner loop:

```sh
ctest --test-dir build -L cosmo -LE heavy     # I changed cosmogenic code
ctest --test-dir build -L network -LE heavy   # I changed routing / chi / basins
ctest --test-dir build -LE heavy              # everything cheap
ctest --test-dir build                        # everything, before committing
```

Measured 2026-08-04 on an 8-core machine, `-j4`:

| Run | Tests | Wall |
|---|---|---|
| `-L step_change -LE heavy` | 2 | 1.5 s |
| `-L cosmo -LE heavy` | 21 | 34 s |
| `-LE heavy` | 93 | 71 s |
| full suite | 103 | ~30 min |

So the cheap tier is **~25x faster than the full suite** and still covers 93 of
the 103 tests. Nothing is lost by using it in the inner loop as long as the full
suite runs before a commit.

#### What is `heavy`, and why `RUN_SERIAL` matters

`heavy` is `unit_hillslope_implicit`, `smoke_lem_roering_implicit_strip`, and
the eight `cosmo_realworld_*` tiers. The first two **saturate every core** — the
implicit hillslope test takes 9 min 24 s wall for 73 min of CPU, a parallelism
of 7.8 on 8 cores.

Both are marked `RUN_SERIAL`. This is not a nicety. Under `ctest -j4` a test
that wants 8 threads and gets 2 does not merely run slower, it can be starved
past its timeout: `unit_hillslope_implicit` "failed" at CTest's default 1500 s
ceiling while passing in 9 minutes when given the machine to itself. A timeout
under `-j` is a scheduling result, not a defect — check for `RUN_SERIAL` before
investigating one as a bug.

#### The domain label is derived, and an unclassified test is a configure error

Domains are derived from the **test name** by a rule table at the foot of
`tests/CMakeLists.txt`, rather than declared at each of the hundred-odd call
sites, so there is one table to read and one place to change. A test whose name
matches no rule **fails the CMake configure step** with a message naming it.

That guard is the point: without it, a newly added test silently belongs to no
compartment and stops being run by anyone working on that subsystem — which is
worse than no compartmentalisation at all, because the compartment run looks
green. When you add a subsystem, add a rule.

#### Do not build while a suite is running

The build replaces binaries the running tests are about to exec, and on Linux
that surfaces as `[Errno 26] Text file busy` inside whichever test was unlucky.
It looks exactly like a real failure. `rigorous_cosmo_cli_bridge` reported
21 of 35 checks failing that way once, entirely because of a concurrent rebuild.

Opt-in tiers, each behind its own CMake option so enabling one does not drag in
the others. **All three default to `OFF`**:

```sh
cmake -S . -B build \
  -DLSDTT_ENABLE_NETWORK_TESTS=ON \   # tests that download
  -DLSDTT_ENABLE_BENCHMARKS=ON \      # RichDEM comparisons (needs pixi)
  -DLSDTT_ENABLE_PERF_TESTS=ON        # performance timings
```

### Keep the timing tiers in a separate build directory

The opt-in tiers are slow by design — they exist to measure, not to check. Turn
them on in the build directory you develop in and every routine `ctest` run pays
for them: on a typical workstation the offline suite is **~14 s**, the same suite
with `-DLSDTT_ENABLE_BENCHMARKS=ON -DLSDTT_ENABLE_PERF_TESTS=ON` is **~17
minutes**, almost all of it the six `benchmark_surface_*` timings.

So use two build directories, and never enable the options in the first:

```sh
cmake -S . -B build && cmake --build build -j            # develop here: fast
ctest --test-dir build --output-on-failure               # ~14 s

cmake -S . -B build-perf \                               # measure here
  -DCMAKE_BUILD_TYPE=Release \
  -DLSDTT_ENABLE_PERF_TESTS=ON -DLSDTT_ENABLE_BENCHMARKS=ON
cmake --build build-perf -j
```

This is better than filtering at run time: the benchmark targets are not even
compiled in `build/`.

If you have inherited a build directory with the options already on — check with
`grep '^LSDTT_ENABLE' build/CMakeCache.txt` — you can either reconfigure it with
`-DLSDTT_ENABLE_BENCHMARKS=OFF -DLSDTT_ENABLE_PERF_TESTS=OFF`, or exclude the
slow tiers **by label** for one run:

```sh
ctest --test-dir build -LE "benchmark|performance|cosmo_realworld|network"
```

Exclude by label (`-LE`), not by test-name regex (`-E`). Labels are assigned
deliberately; names drift — `perf_cosmo_lsdn_production` is labelled
`performance` but is not named `performance_*`, so a name-based filter silently
runs it.

---

## 3. Test data

**No DEM is committed to this repository.** Committed fixtures are tiny
(`tests/data/tiny*.asc`, a 30×30 GeoTIFF, a few kB of vector bounds files) plus
frozen reference CSVs under `tests/expected/`.

Anything needing real topography is fetched on demand by **`lsdtt-fetch-raster`**
— the lsdtt3 CLI itself, not a Python helper — into `tests/data/dem_cache/`,
which is gitignored. The benchmark driver does this automatically the first time
a fixture is needed.

Two terrain regimes, because they behave very differently:

| Fixture | Relief (σ) | Why |
|---|---|---|
| `alpine_small` | 510 m | glacial overdeepenings; deep depressions |
| `lakeland_small` | 2.9 m | broad shallow basins and large flats |

Using only high-relief terrain would have hidden most of §5's findings. The
Alpine tile alone reports near-perfect agreement with RichDEM on epsilon fill;
the low-relief tile reports a factor of 40.

### Determinism caveat

Fetching is deterministic for a **given GDAL build** — repeated fetches
reproduce byte-identical rasters. It is *not* stable across GDAL versions: a
version change altered a fixture by one raster row (a 30 m registration shift,
worth up to 91 m of elevation in Alpine terrain), because the exact-centre crop
computes `(rows − size) / 2` and the warp extent differed by one row.

Consequences, both already implemented:

- Pinned benchmark tolerances record the fixture hash they were measured
  against (`EXPECTED_FIXTURE_SHA`). On mismatch the run **skips** rather than
  applying numbers that no longer mean anything.
- Cached reference outputs are keyed on the fixture's **content hash**, not its
  name. Keying on the name silently compared new output against a stale
  reference and once reported "99.9 % of cells differ" where the truth was 0 %.

> Any cache keyed on a name rather than on content will eventually answer a
> question it was not asked.

---

## 4. External references

RichDEM is the external reference for depression handling. It is installed from
**conda-forge via pixi**, never pip:

```sh
pixi run -e benchmark python tests/scripts/run_flow_benchmark.py \
    --exe build/lsdtt-raster-preprocessing
```

`pip install richdem` yields version 0.3.4, which pins `numpy<2` and will
downgrade numpy in whatever environment it lands in. conda-forge ships 2.4.3
with numpy-2 builds.

### DAGGER / GraphFlood - hydraulic parity reference

Hydraulic parity runs also go through the same pixi benchmark environment
(`numpy` + `gdal` compatibility, and activated PROJ/GDAL data paths):

```sh
pixi run -e benchmark python tests/scripts/run_hydraulic_dagger_benchmark.py \
    --exe build/lsdtt-hydraulic-routing \
    --dagger-root ../DAGGER \
    --dagger-ref tests/scripts/dagger_graphflood_reference.py \
    --scenario stationary_fixed \
    --work-dir build/hydraulic_benchmark_stationary \
    --write-json \
    --check
```

The benchmark writes `hydraulic_benchmark_results.csv` and exits 77 (skip)
when optional DAGGER dependencies are unavailable. In `--check` mode it enforces
both nodata-mask compatibility and pinned parity tolerances.

Real-terrain parity (fetched fixture in `tests/data/dem_cache/`) uses the same
harness. Pass `--fetch-exe` so the driver can fetch on demand with
`lsdtt-fetch-raster`:

```sh
pixi run -e benchmark python tests/scripts/run_hydraulic_dagger_benchmark.py \
    --exe build/lsdtt-hydraulic-routing \
    --fetch-exe build/lsdtt-fetch-raster \
    --dagger-root ../DAGGER \
    --dagger-ref tests/scripts/dagger_graphflood_reference.py \
    --fixture alpine_small \
    --scenario stationary_fixed \
    --work-dir build/hydraulic_benchmark/alpine_small_stationary \
    --write-json \
    --check
```

The pinned tolerance set currently covers `tiny_dem.asc` and `alpine_small`.

CTest wiring (`-DLSDTT_ENABLE_BENCHMARKS=ON`) registers
`benchmark_hydraulic_dagger_*` with the same skip-77 contract.

Hydraulic module test outcomes (2026-07-26):

| Test set | Command | Outcome |
|---|---|---|
| unit properties | `ctest --test-dir build -R unit_hydraulic_routing --output-on-failure` | pass |
| CLI smoke/failure | `ctest --test-dir build -R "smoke_hydraulic_routing|failure_hydraulic_routing_invalid_mode" --output-on-failure` | pass |
| DAGGER parity (checked) | `ctest --test-dir build -R benchmark_hydraulic_dagger --output-on-failure` | pass (4/4: tiny + alpine, stationary + transient) |

Current checked parity fixtures/scenarios:

| fixture | scenario | result row status |
|---|---|---|
| `tiny_dem.asc` | `stationary_fixed`, `transient_dynamic` | `phase4_checked` |
| `alpine_small` | `stationary_fixed`, `transient_dynamic` | `phase4_checked` |

The latest checked benchmark rows are written under:

- `build/tests/hydraulic_benchmark/stationary_fixed_tiny_dem.asc/hydraulic_benchmark_results.csv`
- `build/tests/hydraulic_benchmark/transient_dynamic_tiny_dem.asc/hydraulic_benchmark_results.csv`
- `build/tests/hydraulic_benchmark/stationary_fixed_alpine_small/hydraulic_benchmark_results.csv`
- `build/tests/hydraulic_benchmark/transient_dynamic_alpine_small/hydraulic_benchmark_results.csv`

### TopoToolbox — the D8 reference

Because RichDEM cannot be used to check D8 (§5.1), TopoToolbox's Python package
was evaluated as a second opinion. **It passes the test RichDEM fails.**

`libtopotoolbox/src/flow_routing.c` computes true steepest descent:

```c
float chamfer[8] = {1.0, sqrtf(2.0), 1.0, sqrtf(2.0), ...};
float g = (z - neighbor_z) / chamfer[neighbor];
if (g > max_gradient) { max_gradient = g; direction = query_direction; }
```

Drop divided by distance, diagonals weighted √2 — exactly the rule RichDEM
omits. Verified empirically: on the planar ramp it reproduces the analytic
answer cell for cell, and on a surface engineered so the two rules must disagree
(`z = 100 − row − 0.3·col`, where S is the steepest slope but SE is the lower
cell) all 2450 interior cells route due south, never to the lower diagonal.

```python
import topotoolbox as ttb
dem = ttb.read_tif("dem.tif")                       # or assign .z / .cellsize
fd  = ttb.FlowObject(dem, method="d8", sink_resolution="carve")
acc = fd.flow_accumulation()
```

It also resolves flats by imposed auxiliary topography rather than by perturbing
elevations — closer in spirit to lsdtt3 than to RichDEM's ULP fill.

Practical caveats before wiring it in:

- **PyPI only** (`topotoolbox`; not on conda-forge), so it must enter through
  pixi `[pypi-dependencies]` in the optional feature.
- Wheels exist for **cp311–cp313 only**. The repo's `python = ">=3.10"` can
  resolve to 3.14 and silently build from sdist — which works on a dev box with
  a compiler and fails in a bare CI container. Pin `>=3.11,<3.14` in the
  benchmark feature.
- Version **0.0.12** — pre-1.0, so pin the exact version; the underlying
  algorithms are the mature MATLAB TopoToolbox lineage.
- GPL-3.0, same as RichDEM, and invoked as an optional subprocess dependency —
  no change to lsdtt3's licence posture.

**Now wired in.** `accum_d8` is **asserted against TopoToolbox**; the RichDEM row
for the same comparison is recorded but not asserted. Both references score the
*same* lsdtt3 output, so any difference between the two rows is a difference
between the references, not between lsdtt3 runs.

Measured, on the epsilon-filled fixtures:

| Fixture | reference | rank correlation | top 1 % channels | peak accumulation |
|---|---|---|---|---|
| alpine | **TopoToolbox** | **1.000000** | **0.0000 %** | 41 437 = lsdtt3's 41 437 |
| alpine | RichDEM | 0.508883 | 16.56 % | 40 446 |
| lakeland | **TopoToolbox** | 0.940185 | 15.11 % | — |
| lakeland | RichDEM | 0.792751 | 6.58 % | — |

On high-relief terrain lsdtt3 and TopoToolbox agree **exactly** — same
steepest-descent rule, and steep ground leaves almost no ties to break. That is
about as strong a corroboration of lsdtt3's D8 as an independent implementation
can give.

Low relief does *not* reproduce exactly, and that is expected rather than
alarming: after epsilon filling, flat ground still leaves many near-ties, and D8
must pick one neighbour. Tie-breaking moves *which* cell carries a value without
changing the network — precisely what rank correlation tolerates and a cell-wise
diff would not. The binding pin there is the rank correlation (0.93), not the
channel agreement.

TopoToolbox's `fillsinks()` is also bit-exact against lsdtt3's `min_slope = 0`
fill on both fixtures, which makes the fill result independently confirmed
rather than confirmed only by RichDEM — a sibling port of the same paper.

Speed, on the same 256² raster (all in-routine):

| | fill | D8 accumulation |
|---|---|---|
| lsdtt3 (in-process) | 0.0073 s | 0.0127 s |
| TopoToolbox | 0.0034 s | 0.0106 s (build 0.0102 + accumulate 0.00034) |
| RichDEM | 0.0134 s | 0.0250 s |

TopoToolbox is the faster of the three at filling, and comparable to lsdtt3 at
D8. Note the phase split is not identical: TopoToolbox's `FlowObject` build
includes sink resolution (a no-op on an already-filled DEM), whereas RichDEM's
single call does directions and accumulation without resolving sinks. The
isolated `flow_accumulation()` figure is the one comparable with lsdtt3's
in-process accumulation timing.

Two environment traps worth knowing:

- CTest must launch the benchmark through `pixi run`, not the environment's
  python directly — pixi exports `PROJ_DATA`/`GDAL_DATA` on activation, and
  without it GDAL cannot open the PROJ database.
- Conversely, the lsdtt3 **C++ binaries** must be launched with pixi's library
  paths *stripped*. Left alone they load pixi's `libgdal.so.38`, which needs a
  newer libstdc++ than the system-built binaries have (`CXXABI_1.3.15 not
  found`). The driver's `system_env()` handles this.

---

## 5. Where lsdtt3 and RichDEM deliberately differ

**This section is the reason this document exists.** These are measured,
source-verified differences. If you compare lsdtt3 against RichDEM (or against
another tool with similar conventions) and see disagreement, check here first.

### 5.1 RichDEM's D8 is not steepest descent — do not use it as a D8 reference

`richdem/flowmet/d8_flowdirs.hpp` selects the neighbour with the **lowest
elevation**, with a tie-break that explicitly prefers diagonals. True D8
(O'Callaghan & Mark 1984, which both tools cite) divides the drop by the
**distance**, so a diagonal neighbour must be √2 lower to win.

Demonstrated on a 50×50 plane tilted purely north–south, where every cell must
drain straight down its column (row *r* accumulates *r*+1 cells):

| Row | Expected | lsdtt3 | RichDEM |
|---|---|---|---|
| 25 | 26 | **26** | 1, 2, 3, 4 … |
| bottom | 50 | **50** | 1, 2, 3, 4 … |

On a plane the S, SE and SW neighbours sit at the same elevation, so RichDEM
ties and breaks toward a diagonal, producing diagonal stripes instead of
columns. **lsdtt3 is exactly right here and RichDEM is not.** Planar hillslopes
are ubiquitous in real terrain, so the disagreement is pervasive rather than
confined to flats — measured rank correlation on a real DEM is 0.51 for D8,
0.74 for D-inf, 0.99 for MFD (MFD spreads flow across neighbours, so a
single-direction tie-break barely registers).

Consequently **RichDEM's accumulation comparisons are reported, not asserted.**
Correctness of lsdtt3's accumulation is established instead by property tests,
the analytic ramp oracle, and comparison against **TopoToolbox** (§4), which
implements the correct rule and agrees with lsdtt3 exactly on high-relief
terrain.

A corroborating data point: peak D8 accumulation on `alpine_small` is 41 437
cells for lsdtt3 and 41 437 for TopoToolbox, against 40 446 for RichDEM on the
same raster.

### 5.2 The two "epsilon fills" are different operations

Both tools offer what is loosely called an epsilon fill. They are not the same.

- **lsdtt3** imposes a true **slope**: each cell in a flat is raised by
  `fill_min_slope × cell_size` (×√2 diagonally). On a 30 m grid the default
  gives 3 mm per cell. Uniform across the raster, scales with resolution,
  independent of elevation and data type.
- **RichDEM** raises each cell to the next representable float (`std::nextafter`
  — one ULP). That step depends on elevation and dtype: ~1.2e-4 m at 2000 m in
  float32, ~7.6e-6 m at 100 m, and a denormal ~1e-45 m near sea level. RichDEM
  refuses integer rasters for this reason.

`fill_min_slope: 0` **is** equivalent to `FillDepressions(epsilon=False)` and is
bit-exact against it on both terrain regimes.

Measured share of total fill volume contributed by the flat-gradient term:

| Terrain | Depression fill | Flat-gradient term | Share |
|---|---|---|---|
| Alpine | 6091 m | 0.56 m | 0.01 % |
| Lakeland | 188 m | 78.8 m | 29.5 % |

So the two agree closely in steep terrain and diverge by ~40× in material moved
in flat terrain — not because the step sizes coincide anywhere, but because the
term is negligible where depressions are metres deep and dominant where they are
centimetres deep.

**Practical advice:** if you work in low-relief terrain (coastal plain,
glaciated shield, marsh, floodplain), `fill_min_slope` is doing substantial work
and is worth reducing. The default is tuned for steeper landscapes.

The ULP approach was evaluated for lsdtt3 and **rejected**: it collapses to a
denormal near sea level (exactly where `lsdtt-marsh-topography` operates), varies
by orders of magnitude within a single raster, is resolution-blind, and is
destroyed by any integer-typed output or by warping a filled DEM.

### 5.3 Breaching pit-flattens before carving — in both tools

`breach_depressions` is not a pure breach. Before carving, every pit/flat cell
is raised to its lowest neighbour, so "breaching never raises the surface" is
false by design — measured up to 80 m on a high-relief DEM.

RichDEM does exactly the same thing (verified: 1.0328 m rise on a shared test
fixture against lsdtt3's 1.032768 m, both raising the same 2 cells). This is a
faithful port, not drift. The invariant that *does* hold, and which the property
tests assert, is the bounded form: only pit/flat cells may rise, and never above
the elevation of their lowest original neighbour.

### 5.4 lsdtt3's D8 router requires a gradient on flats

Feeding lsdtt3 a `fill_min_slope: 0` DEM produces:

```
lsdtt_flowinfo: could not assign a receiver at row 12, col 129.
The DEM likely contains an unresolved sink or flat.
```

This is by design — the imposed gradient is what makes flats routable — but it
means the bit-exact `min_slope = 0` fill is not a usable input to lsdtt3's own
flow routing. `fill_min_slope` is load-bearing, not cosmetic.

---

## 6. Performance

Timed **in-process**: raster already in memory, nothing written to disk, best of
3, every row carrying a git SHA so results are comparable across commits.
Reported, never asserted.

### Running the performance suite on your own machine

Use a **separate, Release** build directory (§2) so the timing tiers never slow
down routine testing:

```sh
cmake -S . -B build-perf -DCMAKE_BUILD_TYPE=Release -DLSDTT_ENABLE_PERF_TESTS=ON
cmake --build build-perf -j
ctest --test-dir build-perf -L performance --output-on-failure
```

`-DCMAKE_BUILD_TYPE=Release` is not optional for timing work: a Debug or
unspecified build measures the compiler's inlining decisions, not the algorithms.

That gives two tests, both of which write CSV into `build-perf/tests/perf_output/`:

| Test | Emits | Covers |
|---|---|---|
| `performance_flow_routines` | `flow_performance.csv` | fill, breach, D8/D-inf/MFD build + accumulate |
| `perf_cosmo_lsdn_production` | `lsdn_performance.csv` | Stone vs LSDn per-pixel production; LSDn at 0 / 20 ka / 1 Ma; time-averaged surface vs direct loop, with the measured crossover pixel count; time-varying Stone (`stone2000sp_paleomag`) direct vs scalar-API vs surface |

A third, `performance_surface_metrics`, additionally needs
`-DLSDTT_ENABLE_BENCHMARKS=ON`, pixi and a fetched fixture; it skips with code
77 when any of those is missing.

Every row carries a git SHA, so results are comparable across commits and you
can keep a server's CSVs as a time series.

**Run the drivers directly** for more control than CTest gives — this is the
useful mode on a server, because `--full` and a real raster are where the
interesting numbers are:

```sh
# flow routines: --quick is what CTest runs; --full sweeps more sizes
build-perf/tests/lsdtt-test-flow-performance out_dir --full \
    --raster tests/data/dem_cache/alpine_small.tif

# per-pixel cosmogenic production
build-perf/tests/lsdtt-test-cosmo-physics out_dir --performance
```

Both take the output directory as their first positional argument and create it
if needed.

#### On a shared server

- **Pin the threads.** `OMP_NUM_THREADS=8 ctest --test-dir build-perf -L performance`.
  Left unset, OpenMP takes every core it can see, so a busy node produces numbers
  that are not comparable with anything. Note §6's finding that thread scaling is
  limited and accumulation is *slower* with more threads — more threads is not
  automatically faster here.
- **Run CTest serially** (no `-j`) for the `performance` label. Parallel test
  execution and in-process timing measure each other.
- **Nothing is asserted.** These tests are reported, never asserted — they can
  only fail by not running, so a loaded runner cannot turn the suite red. Correctness
  claims live in the `unit` tier: that the LSDn time-averaged surface and the direct
  per-pixel loop agree, for instance, is asserted by L10 in `unit_cosmo_physics`,
  not by the timing test.
- Timings are **not** committed. Compare against the numbers recorded in this
  document, and treat a large divergence as a question about the machine before
  it is a question about the code.

### Time-varying scaling on a precomputed surface (2026-07-27)

Both time-varying production paths — LSDn, and Stone with the palaeomagnetic
correction — walk a geomagnetic record (108 steps for the Muscheler/Knudsen
record, ~1800 for the 1 Ma LSDn case) at *every pixel*. That per-pixel record
walk is the entire cost, and it is the same walk at every pixel bar the
latitude, so both now precompute a time-averaged surface once per raster and
interpolate it per pixel.

For the Stone palaeomagnetic path the surface is built on **latitude nodes only,
with an exact pressure dependence**: the underlying interpolation table has
pressure-independent breakpoints, so the factor separates exactly into
`Σ_j β_j(Rc)·S_j(P)` and only the latitude axis is ever approximated. The
shipped latitude step of 0.01° puts the horizontal error at ~8e-8 relative —
below a float32 half-ulp, so it is invisible in the output raster. Correctness
is asserted in the `unit` tier (surface versus the direct definition, at and
between nodes; scalar API versus raster surface), not here.

Measured with `lsdtt-test-cosmo-physics --performance`, medians of 5 repeats on
an 8-core WSL2 laptop. Point level, 4000 sites, single-threaded, 20 ka:

| path | median | vs direct |
|---|---|---|
| direct per-site record walk | 220.3 ms | 1× |
| scalar API (surface snapped per call) | 123.6 ms | ~1.8× |
| surface built once, then evaluated | **1.75 ms** | **≈ 126×** |

At raster level, 10⁶ px of time-varying Stone at 20 ka costs **0.937 s**
single-threaded and **0.435 s** at 8 threads — only ~1.40× / ~1.30× a *static*
Stone raster of the same size (0.671 s / 0.335 s). The surface itself is
negligible: 41 latitude nodes and 12 kB for the 4000-site extent, built in
2e-05 s. For LSDn the same pairing gives a measured crossover of ~130 px at
20 ka and ~170 px at 1 Ma, which is why the cost model builds a surface for
anything larger than a very small window.

One caveat, stated because it limits the claim: the pre-surface *raster* path no
longer exists, so a like-for-like raster before/after cannot be measured. The
126× above is a point-level measurement; scaling it to 10⁶ px would suggest
~59×, but that is an extrapolation, not a measurement.

### Do not read performance from the CLI timings

`benchmark_results.csv` records lsdtt3 as **whole-process CLI wall time** —
startup, parameter parsing, GDAL open, the routine, GDAL write — against
RichDEM's in-routine call. Those numbers are not comparable and make lsdtt3 look
~300× slower than it is.

Timed like for like on the same raster:

| Routine | lsdtt3 | RichDEM | |
|---|---|---|---|
| fill (`min_slope=0`) | 0.0073 s | 0.0134 s | **1.8× faster** |
| fill (`min_slope=1e-4`) | 0.0066 s | 0.0209 s | **3.2× faster** |
| breach | 0.0080 s | 0.0250 s | **3.1× faster** |
| D8 build + accumulate | 0.0127 s | 0.0250 s | **2.0× faster** |
| D-inf build + accumulate | 0.0057 s | 0.0379 s | **6.6× faster** |
| MFD build + accumulate | 0.0089 s | 0.0239 s | **2.7× faster** |

### The reconstruction backend (2026-07-20)

Investigating why TopoToolbox filled ~2× faster than lsdtt3 found that it is not
an implementation-quality gap but a **different algorithm class**:

- lsdtt3's `fill_pits` used Barnes et al. (2014) Improved Priority-Flood —
  **O(N log N)**, because on continuous terrain nearly every cell passes through
  the `std::priority_queue` rather than the flat-tie FIFO.
- TopoToolbox's default `fillsinks(hybrid=True)` is Vincent (1993) grayscale
  morphological reconstruction — **O(N)**: two raster scans plus a bounded FIFO
  drain, almost entirely sequential memory access.

The decisive control: TopoToolbox's *non*-hybrid path tracks lsdtt3 almost
exactly at every size, so their implementation quality was never the
differentiator — the hybrid reconstruction specifically was.

lsdtt3 now has both. `fill_pits_reconstruction()` implements Vincent 1993, and
`fill_pits()` dispatches on `min_slope`:

| `min_slope` | backend | why |
|---|---|---|
| `<= 0` | reconstruction | 2–3× faster, bit-identical |
| `> 0` | priority-flood | reconstruction cannot impose a gradient |

Measured on synthetic surfaces (in-process, best of 3):

| size | priority-flood | reconstruction | speed-up |
|---|---|---|---|
| 256² | 0.0051 s | 0.0024 s | 2.1× |
| 512² | 0.0275 s | 0.0091 s | 3.0× |
| 1024² | 0.1142 s | 0.0563 s | 2.0× |
| 2048² | 0.6791 s | 0.2251 s | 3.0× |

On the real `alpine_small` fixture the gain is larger — 0.0073 s → 0.0014 s,
about **5×** — because alpine terrain has few depressions (0.9 % of cells), so
the FIFO drain stays tiny while the heap still has to process everything.

At 2048² lsdtt3's reconstruction (0.225 s) is now slightly **faster than
TopoToolbox's** (0.250 s), which was the original benchmark.

#### A negative result: flat resolution as a second pass does not pay

The obvious way to extend the fast path to `min_slope > 0` is to fill by
reconstruction and then impose the gradient in a separate step. That was
implemented (`resolve_flats_gradient`) and **it is slower**, by 8–88 % against
the priority-flood it was meant to replace, across every size tested. Two
independent implementations were tried — fast-sweeping raster passes, then a
work queue — and both lost.

The reason is structural rather than a tuning problem. Imposing the gradient is
a shortest-path problem over the surface: each cell must end one step above
whichever neighbour it drains to. That is exactly the ordered propagation
priority-flood already performs inline while flooding. Doing it as a second pass
means paying for that propagation separately, on top of the fill, and the heap
it removes costs less than the second pass adds. Solving the second pass
"properly" with a priority queue would just be priority-flood again.

Two bugs found along the way, both caught by asserting **routability** rather
than only depression-freedom — worth recording because a depression-free check
alone passes both:

1. Adding `min_slope × cell_size × distance` post hoc can lift a flat cell
   *above* a neighbouring cell that was previously lower, inverting that pair
   and leaving the neighbour nowhere to drain. A new pit, and the D8 router
   refuses the surface. The fix is to solve the monotone fixed point
   `S(p) = max(F(p), min over q of [S(q) + step])` instead, where every cell
   necessarily ends above whichever neighbour supplied its minimum.
2. In that fixed point, pinning every already-draining cell to its filled
   elevation reintroduces the same failure: if the neighbour it drains to is
   later raised, the pinned cell must rise too, and pinning prevents it. Only
   the drainage *exits* may be pinned; everything else must be computed.

`resolve_flats_gradient` is retained — it is correct, tested, and defines the
gradient by propagated distance rather than FIFO first-arrival, which may be
independently useful — but it is not wired into any CLI and is not a default.

#### The `min_slope > 0` output depends on heap tie-breaking

Measured, not inferred. Changing **only** the tie-break in `pf_node::operator>`
— adding `(row, col)` as a secondary key, algorithm otherwise untouched — and
comparing filled rasters:

| Case | Output |
|---|---|
| `min_slope = 0`, alpine and lakeland | **identical** |
| `min_slope = 1e-4`, alpine (0.9 % flats) | identical |
| `min_slope = 1e-4`, lakeland (flat-dominated) | **different** |

On the flat-dominated fixture: max difference **0.053 m** (≈18 steps of
`min_slope × cell_size`), **22.8 %** of cells differing at all, 4.2 % by more
than 1 cm, mean 0.0016 m.

The mechanism: a cell's filled value is assigned **once, on first arrival**
(`idx != k_open` blocks reprocessing). At `min_slope = 0` the assigned value is
the parent's elevation, and since priority-flood pops in increasing elevation,
every tied parent carries the same value — the result is the spill level, a
uniquely determined quantity. That is why it is bit-exact against RichDEM,
TopoToolbox and reconstruction. At `min_slope > 0` the assigned value is
`parent + step`, where the parent's value already includes its own accumulated
increments, so two tied parents can carry different values and whichever arrives
first wins. `std::priority_queue` leaves the order among equal elements
**unspecified**.

Two consequences worth knowing:

- **A heap swap is not output-neutral for `min_slope > 0`.** An indexed or radix
  heap would break ties differently and shift gradient values on flats by up to
  a few centimetres in low-relief terrain. Not *wrong* — the current values are
  one arbitrary choice among many valid ones — but it would change published
  numbers for anyone working on flat ground, for an estimated 10–20 % on a path
  that is not the bottleneck. Judged poor value; not done.
- **The current output is already unspecified-order dependent.** Nothing
  guarantees `std::priority_queue` behaves identically across libstdc++
  versions, compilers or optimisation levels. Adding a defined tie-break would
  make it deterministic at essentially no cost. Considered and deliberately
  **not** adopted (Simon, 2026-07-20) — recorded here so the exposure is known
  rather than forgotten. It matters only for gradient values on flat terrain;
  depression fill levels, flow routing and everything at `min_slope = 0` are
  unaffected.

#### Why it is not the default everywhere

Reconstruction fills to exactly the spill level and leaves **true flats**. It has
no mechanism for the gradient that `fill_min_slope` imposes, and lsdtt3's own D8
router *throws* on a truly flat surface. Since the default `fill_min_slope` is
`1e-4`, the default path still uses priority-flood — there is no faster option
for the gradient case. The dispatch is therefore an optimisation of the
`min_slope == 0` case only, and changes no output.

`fill_algorithm_choice: vincent1993` selects it explicitly in
`lsdtt-raster-preprocessing`, and **refuses to run** with a positive
`fill_min_slope` rather than silently ignoring it and returning undrainable
flats.

#### How the equivalence is verified

**P9** in `tests/cpp/test_flow_properties.cpp` asserts the two backends agree
**exactly** — `max |difference| == 0.0`, not a tolerance, because both are
deterministic algorithms over the same float32 input and anything less than
exact would mean one of them is wrong. Covered:

- all four standard fixtures, including the interior-nodata one
- six degenerate cases: single row, single column, 3×3, **perfectly flat**
  (every cell ties with every neighbour — the worst case for any tie-breaking
  difference), a monotonic ramp with nothing to fill, and a mostly-nodata raster
  with a valid island containing its own pit

P9 calls `fill_pits_priority_flood()` explicitly rather than `fill_pits()`,
since the latter now dispatches and the comparison would otherwise be circular.

Independently, the RichDEM and TopoToolbox benchmark comparisons run at
`min_slope = 0` and so exercise the new default path — both still report
bit-exact agreement.

### Scaling

Constant throughput from 65k to 1M cells confirms the O(n) behaviour
Priority-Flood advertises: `fill_pits` ~9.6 Mcell/s, `breach_depressions`
~7.0 Mcell/s. `breach_least_cost` is the exception, falling 3.4 → 2.3 Mcell/s
with size, consistent with a per-pit least-cost search rather than a sweep.

### Thread scaling is limited, and accumulation is slower with more threads

At 1024², going from 1 to 8 threads:

| Routine | speed-up |
|---|---|
| `flowinfo_build` (D8) | 1.37× |
| `flowgraph_build` (D-inf) | 1.73× |
| `flowgraph_build` (MFD) | 1.12× |
| `flowinfo_accumulate` (D8) | **0.72×** |
| `flowgraph_accumulate` (D-inf) | 0.83× |
| `flowgraph_accumulate` (MFD) | 0.91× |

Graph construction barely scales because only `build_receiver_graph_d8` is
OpenMP-parallel; donor counts, outlet detection, Kahn's BFS and the s_vector
build are all serial (Amdahl). Accumulation is a serial topological sweep, so
`n_threads > 0` buys nothing and costs thread setup — passing `n_threads: 8` to
a drainage-area run makes accumulation ~35 % slower than passing 1.

**Open:** accumulation should probably ignore `n_threads`, or the parameter
should be documented as applying to graph construction only.

---

## 7. Cosmogenic nuclide testing

> **The CRN numeric core is locked.** The analytical solutions exercised by
> `unit_cosmo_physics` encode years of manual verification. Any change to
> `lsdtt_crn_parameters.*`, `lsdtt_crn_column.*`, `lsdtt_crn_particle.*` or
> `lsdtt_cosmo_raster_cpu.*`, or to the reference values in
> `tests/expected/cosmo/`, requires Simon Mudd's explicit review before merging.

The cosmogenic tests live in `tests/cpp/test_cosmo_physics.cpp` — a standalone
C++ driver built from `tests/CMakeLists.txt`, linked against `lsdtt3_cpu` (the
cosmo sources live there, not in `lsdtt3_core`), registered with CTest under the
labels `unit;cosmo`. Test data lives inside `tests/` and is version-controlled:
parameter files in `tests/params/`, tiny fixtures in `tests/data/cosmo/`, frozen
values in `tests/expected/cosmo/`. Version control is load-bearing here — the
safeguards in §7.1 only work if sentinel changes show up in diffs. Python is
used for orchestration and the one-off CRONUS fetch only; never for domain
generation.

Correctness is pinned at four levels: equation fidelity, internal consistency,
community comparison, and the geomagnetic correction.

### 7.1 The numerics lock and what enforces it

The analytical solutions underpinning CAIRN and the muon-based step change were
checked and re-checked over several years. The suite is built to make it
impossible to alter their numerical behaviour unnoticed.

| Safeguard | What it does |
|---|---|
| **Sentinel pinning suite** | A matrix of exact concentration/ε values — every (nuclide × scheme) at several (ε, d, d_s, d_t, t, ε_old/ε_new, S_tot) points, covering CAIRN Eq. (14), knickpoint Eqs. (7)–(8), the Λ_v solve and the inversion — frozen to full double precision (17 significant digits) in `tests/expected/cosmo/numerics_sentinel.csv` and compared at 1e-12 relative tolerance. Deliberately stricter than physical accuracy requires: its job is to trip on *any* change — refactor, algebraic rearrangement, compiler flag, parameter tweak. **LSDn has a second, separate sentinel** (`numerics_sentinel_lsdn.csv`, §7.9) covering its own data path; keeping the two apart means a change to one cannot be masked by the other. Both are additionally pinned by sha256. |
| **Exact parameter pins** | F, Γ, P0 and λ for every muon scheme asserted with exact `==` equality on the doubles, not a tolerance. A changed constant is a deliberate act and must fail loudly. |
| **Loud, attributed failure** | A sentinel or parameter-pin failure prints an unmissable block naming what the values encode, instructing that the sentinel must not be regenerated, and that any change requires review and sign-off by Simon Mudd. The same text lives in the sentinel CSV header and at the top of the driver. |
| **Gated regeneration** | The sentinel records its generating commit and date. There is no regenerate-on-failure convenience path; regeneration sits behind a separate, deliberately awkward flag (`--regenerate-sentinel-i-have-simons-approval`) so it cannot be run absent-mindedly and is unambiguous in a diff. |
| **Tolerances are pinned too** | Loosening a tolerance to make a test pass is the same class of event as changing the sentinel and gets the same gate. Every frozen bound carries a comment recording the value and date it was measured from. |

Currently frozen bounds, each pinned just above a measurement:

| Bound | Frozen at | Measured |
|---|---|---|
| numerics sentinel | 1e-12 rel | — (drift detector, not an accuracy claim) |
| E7 analytic vs particle integrator | 5e-7 | 9.25e-8 |
| E6 forward→inverse round-trip | 3e-3 | 9.93e-4 (the inversion solver's convergence tolerance) |
| T-C2 vs frozen CRONUS v3 St | 10 % | 6.65 % worst |

The protected surface is `lsdtt_crn_parameters.*`, `lsdtt_crn_column.*`,
`lsdtt_crn_particle.*` and the concentration/inversion functions of
`lsdtt_cosmo_raster_cpu.*`. Reviewers — and any AI-assisted edit — must treat a
sentinel-affecting diff the way golden-value changes are treated in the LEM
benchmark: a reviewable event, never an automatic pass, with the added rule that
here the reviewer must be Simon.

### 7.2 Equation fidelity — reference equations, not the library

The tests do not check the code against itself. The driver carries its own
double-precision re-implementations of the published equations, written from the
papers (Mudd et al. 2016, *Earth Surf. Dynam.* 4, 655–674, the CAIRN method; Mudd
& Roda-Boluda, the ¹⁴C–¹⁰Be knickpoint paper), so agreement means something. The
checks are properties of the *paper equations the code claims to implement*, so a
failure is a defect rather than a modelling disagreement.

| Check | What is asserted |
|---|---|
| E1 | Depth-integration limits: Eq. (14) → Eq. (12) as `d_t → 0`; Eq. (8) → Eq. (7) as the slab shrinks |
| E2 | Step-change end members: `t = 0` → steady state at ε_old; `t → ∞` → steady state at ε_new; `ε_old = ε_new` → steady state for all t; the trajectory stays inside `[C_ss(new), C_ss(old)]` and approaches monotonically |
| E3 | Exponential identity in Eq. (7) evaluated two algebraically equivalent ways — guards against mis-transcribed exponents |
| E4 | Vermeesch Λ_v Newton solve: `Σ_j S_j F_j = S_tot`; `S_tot = 1` ⟹ `Λ_v = 0` and all `S_j = 1`; `S_tot < 1` ⟹ every `S_j < 1` with spallation attenuated least; `S_tot > 1` ⟹ sum identity holds with the ordering reversed (spallation amplified most) |
| E5 | Linearity in P0 and S_tot; `ε → 0` saturation `C → (P0/λ)Σ S_j F_j`; `ε ≫ Γ₀λ` gives `C ∝ 1/ε` |
| E6 | Forward Eq. (14) at imposed ε → invert → recover ε, all nuclides and schemes, with and without snow/self shielding and AMS standardisation |
| E7 | The column/particle integrator agrees with the analytic formula across an ε sweep, depths, shielding and nuclides |
| E8 | Averaging identities: unweighted mean on the inversion path; flux-weighted mean `Σ(Cε)/Σ(ε)` on the prediction path, on a two-rate basin where the weighting actually matters, plus the uniform-ε reduction |
| E9 | AMS standardisation factors match the CRONUS table (Be-10 KNSTD → 07KNSTD 0.9042; Al-26 ZAL94 → KNSTD 0.9134), and applying a factor scales the inverted ε consistently |

Parameter invariants also asserted: the knickpoint Table 1 F-sums are exactly 1
for ¹⁰Be/²⁶Al/¹⁴C under BraucherBorchers; the CAIRN Table 2 ²⁶Al F sums to
0.97525 and must **not** be normalised; ¹⁴C uses P0 = 15.533, F = {0.788, 0.212,
0, 0} (Lupker 2015) and λ = 121e-6 yr⁻¹; λ(¹⁰Be) = 500e-9 and λ(²⁶Al) = 980e-9
yr⁻¹. Schemes without a paper anchor (Granger, Schaller) are frozen as exact
in-code constants with a comment naming their literature lineage (COSMOCALC,
Schaller et al. 2009/2016).

> **Internal cross-checks are not independent oracles.** E7 established that the
> analytic basin formula and the CRN-column particle path evaluate the *same*
> closed form — measured agreement 9.25e-8, which is pure float-raster storage
> error, not a depth-discretisation gap. Independence comes from the driver's own
> reference equations and from the external comparisons in §7.4. The S_tot ≥ 1
> bug (§7.6) proved the point: E6 and E7 were blind to it because both sides
> shared the clamped scaling; the external CRONUS comparison caught it.

### 7.3 Fixture sites

| Fixture | What it is | Why |
|---|---|---|
| Hypothetical point site | 70° N, 1007 hPa, Lal/Stone, no shielding, ε sweep 0–0.20 g cm⁻² yr⁻¹ | The CAIRN paper's own comparison fixture (Figs. 8–9), so the forward curve is directly comparable to the published CRONUS-2.2 (+20→41 %) and CRONUScalc (−3 to −5 %) difference curves |
| Synthetic basin | 40×31 V-valley DEM, Watopia/Teanu UTM 58S via `set_default_synthetic_georeference`, single forced outlet at (row 0, centre column), every pixel provably draining to it; variants with uniform ε and with a two-rate upstream/downstream split | Exercises basin averaging and the flux weighting; generated in C++, so no user ever goes into Python to build an input |
| Real-coordinate sample set | 15 synthetic samples: 5 latitudes (0–70°) × 3 elevations (250/1500/3000 m), ε cycling 3×10⁻⁴–3×10⁻² g cm⁻² yr⁻¹, both nuclides | Maps the difference surface against CRONUS. CAIRN found the differences are a function of S_tot rather than ε — confirmed, and it was the elevation axis that exposed the S_tot ≥ 1 bug |

### 7.4 Community comparison against CRONUS

The answers are **expected to differ**. Each comparison records the percent
difference and asserts only that it stays inside a frozen band.

| Comparison | Difference | Cause |
|---|---|---|
| vs CRONUS-2.2 (St) | up to ~40 % in C, a function of S_tot | higher P0 (4.49) plus Heisinger muon overproduction |
| vs CRONUScalc | 2.5–5 % in ε | muon scheme plus spatial averaging |
| vs CRONUScalc running newCRONUS parameters | < 2 % | residual spatial averaging and fixed Λ₀ = 160 |
| **vs CRONUS v3 St, default `BraucherBorchers`** | **measured 2026-07-16: +1.8 % to +6.65 %, mean +3.6 %; band frozen at 10 %** | Balco (2017) model-1A muons; grows with ε as the muon share rises |
| **vs CRONUS v3 St, `cronusEmulatorSt`** | **measured 2026-08-01: mean −0.23 %, max 1.03 % (Be-10)** | purpose-built to reproduce v3; see §7.4.1 |
| vs CRONUS v3 Lm | mean +14.6 %, max +45 % against our time-independent rates | different VDM record |
| vs CRONUS v3 LSDn | mean +11.8 % against our time-independent rates | LSD framework |

**Getting the numbers, and freezing them.** CRONUS v3's MATLAB source is not
public, but the calculator is scriptable: an HTTP POST of the plain-text input
block to `…/cgi-bin/matweb` with form fields `mlmfile=erosion_input_v3`,
`reportType=XML`, `resultType=long` (endpoints: hess.ess.washington.edu, with a
stoneage.hzdr.de mirror) returns XML with St/Lm/LSDn erosion rates and internal
and external uncertainties. This is the same route the open-source `riversand`
package uses.

The workflow is **fetch once, freeze, compare offline**.
`tests/scripts/fetch_cronus_reference.py` (orchestration only, per house rules)
submits the fixture samples, caches the XML and writes
`tests/expected/cosmo/cronus_v3_reference.csv` — currently 90 rows, St/Lm/LSDn ×
Be-10/Al-26 — recording the generating date and inputs. **CI never touches the
network**: the comparison test reads the frozen CSV. Refreshing the reference is
a manual, reviewable event, and diffs to it are reviewed like golden-value
changes. Be polite to the server: batch, cache and throttle.

Secondary offline oracles, both scriptable and network-free: the downloadable
**CRONUS v2.3 MATLAB** source (runs under Octave; per-function GPL downloads) and
**expage** (Heyman; open MATLAB/Octave, LSDn-based). These are useful for
isolating *causes* rather than just measuring totals — v2.3's `stone2000.m` lets
the Stone scaling factor be compared directly, separating scaling differences
from muon-scheme differences.

Difference attribution is run in three configurations: (i) our default scheme vs
CRONUS St, giving the total; (ii) our newCRONUS scheme vs CRONUScalc/v3, with the
muon-scheme difference removed; (iii) our Stone factor vs `stone2000.m` at the
same latitude and pressure, pure scaling. Atmosphere differences (ERA40 vs NCEP)
are isolated by running in standard-atmosphere mode against `pre` (direct
pressure) CRONUS inputs.

The comparison tests are **T-C1** (forward-concentration difference curve at the
70° N fixture, reproducing the CAIRN Fig. 9b shape and bounds), **T-C2**
(inversion difference for the real-coordinate sample set against the frozen v3
St columns; band asserted), **T-C3** (Stone scaling factor against independent
Stone (2000) pins — tight, near-exact) and **T-C4** (newCRONUS mode vs
CRONUScalc/v3, expected < 2–5 %; folded into the T-C2 machinery and enabled by
adding a newCRONUS concentration column to the fixture).

#### 7.4.1 The CRONUS v3 emulator schemes

`cronusEmulatorSt` and `cronusEmulatorLSDn` exist for one purpose: to reproduce
the CRONUS v3 online calculator's erosion rates. **Use them when agreement with
that calculator is the goal.** For everything else the default
`BraucherBorchers` remains appropriate — it is a physical parameter set, not an
emulation target, and its few-percent offset from v3 is a real difference
between production models rather than an error.

**What they emulate, and why the obvious approach fails.** CRONUS v3 ships two
muon models. One is a single-exponential "fast approximation" with a
pressure-dependent effective attenuation length; the other is a numerically
integrated model-1A muon inventory, tabulated on a (pressure, log erosion rate)
grid. **v3 uses the exponential only for exposure ages, saturation checks and
plotting — its erosion-rate root-find interpolates the table.** The exponential's
attenuation length is described in its own source as a *near-surface* value, and
an eroding column integrates muon production to a depth of order ε/λ, which for
Be-10 at 0.05 g cm⁻² yr⁻¹ is ~10⁵ g cm⁻² — kilometres of rock, dominated by the
deep fast-muon tail. Measured against v3's own table the exponential delivers
only **17.6–78.5 %** of the true muon inventory across the erosion rates that
matter for erosion-rate work (Be-10, ε = 0.01–0.55 g cm⁻² yr⁻¹), worsening with
altitude; over the full tabulated grid the ratio spans 0.175–1.037, approaching
parity only in the near-saturation corner where the attenuation length barely
enters. Frozen at `tests/expected/cosmo/nmu_exponential_error_report.csv`.

lsdtt3 therefore fits the **tabulated** inventory using the three muon pathways
the concentration kernel already provides:

    N_mu(ε, p) = Σ_k P_k(p) Λ_k(p) / (ε + Λ_k(p) λ)     k = 1..3

with `log P_k` and `log Λ_k` quadratic in pressure, `Λ_1 < Λ_2 < Λ_3` enforced
structurally so the ordering cannot swap with pressure. Fit quality against the
table over 450–1020 hPa and 2.7e-6–2.7 g cm⁻² yr⁻¹, frozen in
`tests/expected/cosmo/cronus_muon_emulator_fit_report.csv`:

| exponentials | Be-10 | Al-26 | C-14 |
|---|---|---|---|
| 1 (a single-exponential scheme) | 22.4 % | 19.8 % | 10.6 % |
| 2 | 2.8 % | 2.9 % | 1.9 % |
| **3 (shipped)** | **0.77 %** | **0.87 %** | **1.15 %** |

Coefficients are frozen in `tests/expected/cosmo/cronus_muon_emulator_coeffs.csv`
and regenerated by `python_helpers/fit_cronus_muon_emulator.py`, which reads the
frozen table extraction rather than the GPL reference tree. The surface
production is **not** anchored to Balco's published `P_mu,0(p)` formula: the
table itself departs from that formula by up to 3.0 % (Be-10), and anchoring
measurably pulls the fit off ground truth (Be-10 0.77 % → 1.57 %).

**Structure.** v3 applies the spallation scaling and topographic shielding to
the spallation pathway only — its muon lookup takes erosion rate and pressure
and nothing else. The emulator reproduces that split. Snow and self shielding
are genuine depth attenuation and apply to all pathways. The muon terms are
identical between the St and LSDn variants, because v3 calls the same muon
lookup for both; only the spallation reference production differs.

**Three things to know before using them.**

1. They are **pressure- and nuclide-dependent**, so they are not reachable
   through `reset_scaling()`, which has neither argument. It throws a message
   naming the right entry point. The raster paths dispatch to them per pixel.
2. They are **steady-state only**. The step-change and transient concentration
   paths reject them explicitly.
3. `lsdtt-cosmo-column-tool` does not offer them: it works on a single column
   with no pressure raster.

**Validation against the live calculator (2026-08-01).** 210 points, latitude
{−60, −30, −1, 30, 60} × elevation {0, 500, 1000, 1500, 2500, 3500} m × erosion
rate {0.001, 0.01, 0.05, 0.1, 0.3, 0.55, 1.0} g cm⁻² yr⁻¹, `sample_age_yr = 0`,
sign `100 (lsdtt3 − v3)/v3`:

| run | mean | max abs | latitude trend |
|---|---|---|---|
| `cronusEmulatorSt`, Be-10 | **−0.233 %** | 1.028 % | 0.114 pp |
| `cronusEmulatorSt`, Al-26 | −0.255 % | 1.083 % | 0.108 pp |
| `cronusEmulatorSt`, C-14 | −0.252 % | 4.234 % | 0.255 pp |
| `BraucherBorchers`, Be-10 (baseline) | −4.638 % | 10.534 % | 2.767 pp |
| `cronusEmulatorLSDn`, Be-10 | −3.366 % | 26.573 % | 9.154 pp |

C-14's 4.2 % worst case is confined to ε = 0.001, where τ = Λ_sp/ε is ~28 C-14
half-lives: the nuclide is saturated and the inversion is ill-conditioned, so a
sub-percent concentration difference maps to a several-percent erosion rate. At
every other rate C-14 stays inside 1.25 %.

**`cronusEmulatorLSDn` does not meet the same standard, and the cause is not the
muons.** v3's LSDn root-find integrates a time-varying spallation production
history that lsdtt3's steady-state inversion has no analogue for, and no single
`sample_age_yr` makes the two agree by construction. The residual's shape
confirms it — −8 % at both poles, +1 % at the equator, asymmetric between −30°
and +30° — which is spallation- and field-side structure, not muon-side. **The
muon coefficients must not be retuned to absorb it**: they are shared with the
St variant by construction, which passes at −0.23 %.

**A comparison trap worth knowing.** The scan driver defaults to
`--sample-age-yr 20000`. Under `scaling_scheme: stone` that switches lsdtt3 to
the Muscheler VDM paleomagnetic correction, while **CRONUS v3's St is by
definition non-time-dependent**. Comparing at age 20000 is apples to oranges and
carries a latitude shape, because paleomagnetic corrections vanish as cutoff
rigidity goes to zero: the same emulator grid reads mean +1.00 %, max 4.87 %,
trend 2.99 pp at age 20000 against −0.23 %, 1.03 %, 0.11 pp at age 0. **Any St
comparison against CRONUS must use `--sample-age-yr 0`.** Earlier scans carry
this term; scheme-versus-scheme comparisons at equal age remain valid, but
absolute residuals quoted from them are inflated.

**Test coverage.** `P2-T1` asserts the structural invariants (Γ₀ = 160, ordered
Γ, ΣF = 1, spallation production equal to v3's reference, pressure clamped to
the fitted range, guards, and that St and LSDn share bit-identical muon terms).
`P2-T2` reproduces the frozen tabulated inventory at all 560 grid points per
nuclide inside the fit gates — that is what catches the C++ drifting from the
Python fit. Eight cases in `rigorous_cosmo_cli_bridge` close the forward↔inverse
loop through both CLIs across erosion rate, nuclide, shielding and both
variants, to ≤ 0.27 %.

#### Muon-scheme isolation decomposition (M7)

`phase8_muon_isolation_comparison()` is the end-to-end decomposition that led to
the emulator described in §7.4.1: the same 15-sample fixture, age-matched
inversion (`t = Γ_sp/ε`), compared against the frozen CRONUS v3 **St** and
**LSDn** columns and emitted as `muon_isolation_comparison.csv`. Reported sign is
`100*(CRONUS − lsdtt)/lsdtt`. Nuclide coverage is Be-10, Al-26 and C-14.

The CSV carries explicit solve diagnostics per row — `solve_status`,
`hit_lower_bound`, `hit_upper_bound`, `eps0_gcm2yr` — so bracket-limited rows are
machine-identifiable rather than inferred from huge percentages.
`tests/scripts/summarize_muon_isolation.py` prints stable per-variant counts and
min/max/mean summaries for review diffs.

**As of 2026-08-01 this harness carries one variant, `bb_vermeesch`**
(BraucherBorchers constants plus Vermeesch pathway redistribution), and serves as
an end-to-end regression guard on the default scheme. Three comparison variants
were removed along with the setters they exercised:

- two variants built on a setter reproducing v3's single-exponential muon
  approximation — the wrong muon model for erosion rates (§7.4.1);
- one built on a fit to the tabulated inventory made at a **single pressure**
  (1013.25 hPa), which froze its attenuation lengths at their sea-level values
  and degraded from ~6–10 % at sea level to 37–47 % at 700 hPa. Its fit report is
  retained at `tests/expected/cosmo/balco_1A_exp_fit_report.csv` as the record of
  what a pressure-independent fit costs.

Both failure modes are superseded by the pressure-dependent emulator, which is
validated far more strongly — against the tabulated inventory at every grid point
(`P2-T2`) and against the live calculator (§7.4.1) — so the decomposition's role
as a *diagnostic* is complete. What it retains is the `bb_vermeesch` regression
guard, `muon_run_e_bb_worst_pct = 15.0`.

**What the decomposition established**, and why it is worth recording even though
the variants are gone: the residual against CRONUS v3 was dominated not by the
choice of muon constants but by two structural differences — which muon model is
being reproduced at all, and whether the site scaling is redistributed across
every pathway (Vermeesch) or applied to spallation alone (as v3 does). Applying
Vermeesch redistribution with Balco's constants was severe enough to collapse
some high-elevation inversions onto the lower bisection bracket entirely. Both
lessons are built into the emulator.

### 7.5 Geomagnetic correction — no published benchmark

The Muscheler-VDM + ERA40 machinery (`crn_nc_file` + `sample_age_yr`) is
post-CAIRN and appears in neither paper, so it is tested through invariants
rather than against an external answer.

- **G1 — static-field and degenerate limits.** `sample_age_yr = 0`, and a
  constant-VDM record equal to the modern value at any age, must reproduce the
  time-independent Lal/Stone result to tight tolerance. This validates the
  time-integration loop itself. The degenerate fallbacks (`nc.loaded == false`,
  `age_yr <= 0`) are covered here too.
- **G2 — physics invariants.** (a) At high geomagnetic latitude production is
  insensitive to VDM — time-integrated scaling at 70° N changes ≪ 1 % across the
  record; (b) at low latitude higher VDM ⟹ lower production ⟹ higher inverted ε
  for the same concentration (sign and monotonicity checked at the equator); (c)
  for an age spanning a Laschamp-type low-field interval, integrated production
  exceeds modern-field production at low latitude; (d) ERA40 pressure at the
  fixture sites versus standard atmosphere: differences small, bounded, frozen.
- **G3 — comparison against CRONUS Lm/LSDn.** The T-C2 sample set rerun with
  `sample_age_yr > 0` against the frozen Lm and LSDn columns. Differences are
  systematic because the VDM records differ. The more meaningful check is
  relative: the *spread between our static and geomagnetic runs* should be
  comparable in magnitude and sign to the *spread between CRONUS St and Lm*.
- **G4 — golden trajectories.** Time-integrated scaling factors for a few (lat,
  pressure, age) triplets frozen with their generating commit.

`CrnNcData` is a plain struct, so a substantial subset of these needs **no
external file**: the driver builds synthetic VDM records and uniform
pressure/temperature grids in code, consistent with the no-external-data rule.
That covers G1 in full, all of G2's sign and monotonicity invariants, and ERA40
interpolation against exactly predictable grids. Only G4's goldens against the
real Muscheler record and validation against the real ERA40 climatology need
`crn_parameters.nc` — a file that was **not found anywhere under
/mnt/c/LSDTopoTools on 2026-07-16**; its location and redistributability remain
open.

### 7.6 Findings — two numerics defects the suite caught

Both were found on the suite's first day, both fixed the same day, and both fixes
were reviewed and approved by Simon Mudd.

**`advance_conc()` transient update (found by E2).** The implementation computed
`exp_decay * (exp(alpha*dt) − 1) / (alpha + lambda)`, algebraically
`(exp((alpha−lambda)·dt) − exp(−lambda·dt))/(alpha+lambda)`, where CAIRN Eq. (4) /
knickpoint Eq. (7) — and the function's own header comment — require
`(exp(alpha·dt) − exp(−lambda·dt))/(alpha+lambda)`. The growth exponent was
missing `+lambda·dt`. Verified numerically: with `eps_old == eps_new` at
t = 23 456 yr the concentration was low by exactly `1 − exp(−λt)` = **1.166 %**.
Consequences: every transient result biased low by ~λt on the accumulation term —
negligible for small λt, **94 % error at t ≈ 23 kyr for ¹⁴C**; as `t → ∞` the
result decayed to 0 instead of to steady state at ε_new; and the
`exp(−d0/Γ)·exp(alpha·dt)` grouping overflowed to `0 × inf = NaN` for
`t·ε/Γ ≳ 700`. Six call sites were affected
(`calculate_crn_concentration_raster_step_change` ×3, `lsdtt_crn_column.cpp` ×3);
steady state and inversion were untouched. The fix adopts the numerically stable
regrouping `exp(−(d0 − ε·dt)/Γ) − exp(−d0/Γ)·exp(−λt)`, which also removes the NaN
overflow. Sentinel regenerated under the approval flag; only the six `t > 0`
`conc_step` rows changed, exactly as predicted.

**Production scaling silently dropped when scaling × shielding ≥ 1 (found by
T-C2).** `scale_one_nuclide` returned immediately for `total_shielding >= 1.0`
("no shielding — nothing to do"). Correct for pure shielding factors in [0,1],
but every raster-pipeline call site passes **Stone production scaling ×
topographic shielding** through it, so wherever that product was ≥ 1 — i.e. most
sites above a few hundred metres — production scaling was dropped and
concentrations computed as if S_tot = 1. The CRONUS comparison exposed it because
inverted-rate differences grew with elevation exactly in the ratio of the Stone
factor: **250 m: 5–32 %; 1500 m: 90–370 %; 3000 m: 330–1020 %** (e.g. sample
LSDTT-12 at 60°/3000 m: ratio 11.2 against a Stone factor of 11.01). The
Vermeesch redistribution is well defined for S_tot > 1 (negative effective depth;
Newton converges — CAIRN Eqs. 15–16 place no upper bound on S_tot). The guard
became `if (total_shielding == 1.0) return;`, E4 was extended with S_tot > 1
cases, and the CRONUS reference was re-fetched because the first submission's ten
Stone > 1 samples carried clamped concentrations. **The sentinel was
bit-unchanged after the fix** — all pinned cases use S_tot < 1 — which is itself
the evidence that the fix was surgical.

The lesson generalises: *internal cross-checks cannot see a defect that both
sides share. That is exactly what external references are for.*

### 7.7 Real-landscape validation tier (download → run → delete)

Published basin-averaged erosion rates computed with other calculators — the
CAIRN paper's seven datasets (Bierman 2005, Dethier 2014, Kirchner 2001, Munack
2014, Palumbo 2010/2011, Safran 2006, Scherler 2014) — are the ultimate
comparison target. They need real DEMs: **hundreds of MB to GB that must never be
stored in the repository or tracked in version control.** The rule that follows:

- **Tracked (small):** per-dataset sample tables (lat/lon, concentration,
  uncertainty, standard, published rate and original-calculator metadata) as CSVs
  under `tests/expected/cosmo/published/`; a manifest per dataset naming the DEM
  source (provider, dataset ID, bounding box, resolution); and the frozen
  *summary* of each validation run (per-sample lsdtt rate vs published rate,
  percent difference). A few kB each, reviewable in a diff.
- **Not tracked (large):** DEMs and every intermediate raster. The harness
  downloads each dataset's DEM into a temporary working directory, runs the
  comparison, appends the small summary CSV, then **deletes the working
  directory — including on failure, via a cleanup trap**, so no bloat survives a
  crashed run. A `--keep` flag exists for debugging one dataset; the default is
  always delete.
- **Not in CI.** This tier is manually invoked; it needs network, time and
  possibly an OpenTopography API key (from the environment — never in the repo;
  see §9). CI-tier cosmo tests use only the tiny tracked fixtures and frozen CSVs.
  Comparing a *new run* against the *frozen summary* of the last accepted run is
  itself a cheap tracked regression check, so the expensive download tier only
  needs rerunning when the code or a manifest changes.
- **Match the published product.** CAIRN used SRTM 90 m throughout, so validation
  runs default to the same product for comparability rather than to the best
  available data. Per-sample regression targets for the CAIRN datasets live in the
  LSDTopoTools GitHub data repositories, not in the paper tables; the manifests
  point there.

#### 7.7.1 riversand cross-check (manual, not in CI)

The Phase 5 riversand harness lives at
`tests/scripts/run_cosmo_riversand_comparison.py` and follows the same
download -> run -> delete discipline. In default mode it runs riversand directly
from the harness (with `lsdtt-basin-extraction`), then compares against the
tracked frozen summary `<dataset>_riversand_summary.csv`.

```bash
# check against frozen summary (default mode)
python3 tests/scripts/run_cosmo_riversand_comparison.py dethier2014 \
    --exe /mnt/c/LSDTopoTools/lsdtt3/build/lsdtt-cosmo-tool \
    --fetch-exe /mnt/c/LSDTopoTools/lsdtt3/build/lsdtt-fetch-raster \
    --samples S169WTS035

# freeze/update the tracked summary
python3 tests/scripts/run_cosmo_riversand_comparison.py dethier2014 \
    --exe /mnt/c/LSDTopoTools/lsdtt3/build/lsdtt-cosmo-tool \
    --fetch-exe /mnt/c/LSDTopoTools/lsdtt3/build/lsdtt-fetch-raster \
    --samples S169WTS035 --freeze
```

You can still bypass built-in execution with `--riversand-csv` (out-of-band
results) or `--riversand-cmd` (external invocation template).

### 7.8 Plottable outputs

Every cosmo test family writes tidy long-format CSVs into the test output
directory, with consistent column names (`nuclide, scheme, lat, pressure_hPa,
erate_gcm2yr, conc_atoms_g, source, pct_diff, t_yr` as applicable) so they plot
directly with lsdviztools3. The driver prints `PASS`/`FAIL` per test and a final
`RESULT:` line regardless; the CSVs are a side product for inspection and for
documentation figures.

| File | Contents |
|---|---|
| `conc_vs_erate.csv` | C(ε) at the synthetic-basin outlet |
| `conc_vs_erate_70N.csv` | C(ε) per scheme at the 70° N / 1007 hPa CAIRN fixture — the Fig. 9a analogue |
| `step_change_relaxation.csv` | C(t)/C_ss for ¹⁰Be and ¹⁴C after a step change — the knickpoint Fig. 1 analogue, showing ¹⁴C re-equilibrating orders of magnitude faster |
| `analytic_vs_integrator.csv` | the full E7 comparison grid with relative differences |
| `cronus_fixture_samples.csv` | the 15-sample comparison fixture (inputs) |
| `cronus_comparison.csv` | per-sample ε from lsdtt vs CRONUS St/Lm/LSDn with percent differences |
| `cronus_comparison_lsdn.csv` | Phase 7/5 LSDn-column comparison vs frozen CRONUS reference at `t=0` and `t=Gamma/eps` |
| `cronus_comparison_geomag_full.csv` | Phase 5 G3 full-inversion comparison vs CRONUS `Lm`/`LSDn` with explicit `solve_status` |
| `geomag_on_off_LSDTT1.csv` | equatorial (`LSDTT-1`) geomag ON/OFF delta sweep for §6c(A) |
| `geomag_era40_fixture_validation.csv` | real-atmosphere ERA40 pressure/scaling deltas vs standard atmosphere over the 15-point fixture |
| `geomag_scaling.csv` | six frozen Knudsen-record G4 trajectory triplets (plus modern ratio context) |
| `lsdn_vs_stone.csv` | LSDn against Stone across the fixture grid |
| `lsdn_scaling_surface.csv` | the LSDn scaling factor over the (lat, p) fixture grid |
| `lsdn_surface_vs_direct.csv` | time-averaged surface against the direct per-pixel loop |
| `lsdn_surface_resolution_sweep.csv` | surface-vs-direct error against latitude step |
| `stone_paleomag_surface_vs_direct.csv`, `stone_paleomag_resolution_sweep.csv` | the same pair for the Stone palaeomagnetic path |
| `g4_paleomag_refreeze.csv` | per-triplet old/new/shift for the six geomagnetic goldens — written every run, as a report; the driver never edits the goldens |
| `lsdn_pixel_vs_lumped.csv` | per-basin per-pixel LSDn against LSDn evaluated once at the basin's published centre |

Tracked/frozen realworld artefacts (under `tests/expected/cosmo/published/`) also
include riversand summaries, e.g.
`dethier2014_riversand_summary.csv` from
`tests/scripts/run_cosmo_riversand_comparison.py`.

### 7.9 LSDn production scaling, and what "validated" means for it

LSDn (Lifton–Sato–Dunai 2014) is the nuclide-specific, spectrum-based
alternative to Lal/Stone. lsdtt3 computes it **per pixel by default**: each pixel
gets its own latitude, longitude and pressure, and when a sample age is supplied
the geomagnetic record is time-averaged over the exposure window. That default is
the point of the implementation, and it is also what makes validation awkward —
almost every published LSDn rate was produced by evaluating the scaling **once**
at a single lumped position for the whole basin, so a naive comparison conflates
"is our kernel right?" with "does per-pixel differ from lumped?".

The suite therefore keeps those two questions apart, and a third one too. For a
given basin there are three quantities:

- **(a) the per-pixel LSDn basin mean — this is the product**, what the tool
  actually returns;
- **(b) LSDn evaluated once** at the basin's published centre latitude and
  pressure — how OCTOPUS's **LSDn-equivalent comparison rate** is produced, and
  how RIVERSAND works generally;
- **(c) the published LSDn rate** itself.

> Scope, because this is easy to over-read: (b) describes OCTOPUS's **LSDn**
> column only. OCTOPUS's headline erosion rates are computed **per pixel using
> CAIRN**, in the same spatially distributed way lsdtt3 computes them, and the
> per-pixel-versus-lumped term below does not apply to them.

Comparing **(b) against (c) validates the scaling kernel** and needs no DEM at
all. Comparing **(a) against (b) is a scientific result, not an error term** — it
is the size of the per-pixel correction. Comparing **(a) against (c) is a
monitoring band and nothing more**, because it mixes both of the above with DEM
and delineation differences.

#### Sign conventions — read the emitting test before combining any two numbers

The cosmo comparisons publish percent differences in **two opposite
conventions**, and figures from one must never be added to figures from the
other:

| comparison | expression | so a positive number means |
|---|---|---|
| kernel oracle, (b) vs (c) | `100·(reference − lsdtt3)/lsdtt3` | the **reference** is higher |
| per-pixel vs lumped, (a) vs (b) | `100·(per_pixel − lumped)/lumped` | per-pixel is higher |
| real-landscape tiers, (a) vs (c) | `100·(lsdtt3 − published)/published` | **lsdtt3** is higher |

This is not hypothetical bookkeeping. During development a residual band was
predicted by subtracting a per-pixel figure from a kernel figure, which silently
mixed the first two rows and produced a number with the wrong sign. Before
combining or quoting any two of these, read the actual expression in the test
that emits it — do not infer it from the sign, from a column name, or from a
summary table. When quoting one, say which direction it runs.

#### The sentinel

LSDn has its own full-double-precision sentinel,
`tests/expected/cosmo/numerics_sentinel_lsdn.csv`, alongside the pre-existing
`numerics_sentinel.csv`. Both are pinned by sha256 in the driver and neither may
be regenerated without explicit sign-off (§7.1). They are deliberately separate:
the LSDn grids are a different data path, and keeping them apart means a change
to one lights up without being masked by the other.

#### Measured

**Kernel, (b) vs (c), DEM-free, 50 OCTOPUS basins.** Mean **+4.34 %** (min −1.18,
max +16.08) using Balco's `P0`; **+6.54 %** on lsdtt3's own `P0`. Independently,
a direct comparison against CRONUS v3 gives **+4.13 %** — two routes to the same
few-percent kernel offset, which is the cross-check that matters. The spread
tracks erosion rate (r = +0.54 against log ε) rather than position (latitude
r = −0.18), which is the signature of a muon-scheme difference; that is offered
as an attribution, not asserted as a test.

**Per-pixel versus lumped, (a) vs (b) — the scientific result.** Over 47
delineation-reliable basins at 5 kyr: **median −0.525 %**, IQR
[−1.590, +0.141], mean −0.827, min −9.211, max +6.735. Inverted to erosion rate,
median −0.440 %. This sits squarely on the independent published estimate of the
same quantity — Codilean & Munack (2025, *Geochronology* 7:113, p. 117) report an
LSDn median of −1.0 % with IQR −1.8 to −0.1 % over 290 points: same sign, same
order of magnitude, overlapping IQRs. Age sensitivity is negligible (median
−0.540 / −0.525 / −0.532 % at 0 / 5 k / 20 kyr), because the age is matched on
both sides by construction and what remains is purely the spatial term.

A second result falls out of the same machinery: inverting our per-pixel **Stone**
basin mean for the pressure that would reproduce it at the published centre
recovers OCTOPUS's own published pressure to a median of **+0.73 hPa** (IQR
−0.18 to +2.19). That is ~0.1 % of a ~900 hPa pressure. It confirms our reading
of their reduction, and bounds the caveat that their pressure was fitted under
Stone but is being consumed by an LSDn kernel — the effect is roughly a fifth of
the spatial term, not a competitor to it.

**End-to-end, (a) vs (c) — the monitoring band.** Three OCTOPUS-derived datasets,
sign `100·(lsdtt3 − EBEL)/EBEL`, delineation-unreliable basins excluded. The tier
is **age-matched** — it runs with `self_consistent_age`, so both sides integrate
over the same erosional timescale:

| dataset | n | median | mean | min | max |
|---|---|---|---|---|---|
| Clapp 2002 (Yuma Wash) | 21 | −4.59 | −4.78 | −8.27 | −3.57 |
| Dethier 2014 (Front Range) | 12 | −5.85 | −5.85 | −8.08 | −3.71 |
| DiBiase 2010 (San Gabriel) | 14 | −8.23 | −8.76 | −21.22 | −2.65 |

**These are bands to watch, not accuracy figures, and nothing may be tuned to
shrink them.** Each still mixes the kernel residual, the per-pixel-versus-lumped
term, and the DEM and delineation difference.

The age matching is worth recording, because it is a prediction that was tested.
The tier originally ran at the CLI default `sample_age_yr = 0` (modern field)
while the published `EBEL` is time-averaged over Γ/ε, and it measured **+3.82 /
+4.28 / +10.14 %** — positive. Expressed in the real-landscape sign and on
lsdtt3's own `P0`, the age-matched decomposition was predicted, in advance, as
kernel −6.14, per-pixel-versus-lumped −0.53, ERA40-rather-than-NCEP atmosphere
−0.64, topographic shielding +0.70 ⇒ ≈ **−6.6 %**. Re-running age-matched flipped
the sign and landed on it, for two of the three datasets; DiBiase sits wider, as
it does in every other comparison, because those basins erode up to 10× faster
and the spread tracks erosion rate — the muon-fraction signature. Reading the
original +4 to +10 % as "accuracy" would have been wrong in both magnitude and
sign, and the corrected figure is still not an accuracy number.

Two basins remain excluded as mis-delineated rather than wrong: `S169WTS036`
(0.07 km²) and `S055WTS044` (0.08 km²), both below the snap-threshold footprint.
`S055WTS044` is not fully explained even by that rule and is worth raising with
the OCTOPUS maintainers.

#### Closing the age loop

The age mismatch above is not a defect but it is avoidable, and `lsdtt-cosmo-tool`
can now avoid it. `self_consistent_age: true` solves the age and the erosion rate
together — invert, form `t = Γ_sp/ε`, rebuild production, re-invert, repeat — which
is the CRONUScalc / Balco v3 convention that published LSDn rates are computed
under. It is off by default, because it changes what the reported rate means and
because it costs a production-raster rebuild per iteration **per sample**.

Three CLI tests cover it, and each was verified to fail on a deliberate
perturbation before being kept:

- the iteration really iterates — asserted by requiring a *second* pass in the
  output, so a loop that degenerated to the fixed-age path cannot pass; forcing a
  single pass fails it;
- non-convergence is **reported, not absorbed** — a capped run must emit the
  warning and write `self_consistent_converged = 0`; making the flag optimistic
  fails it;
- the diagnostic columns are present in the CSV, so the converged age is a
  recorded per-sample result rather than console narration; dropping them fails it.

A non-positive tolerance or a zero iteration cap is fatal up front, because either
would let the loop exit on its first comparison and report a "converged" age that
is nothing of the kind.

#### The CLI bridge: forward against inverse

Everything above drives the library. `rigorous_cosmo_cli_bridge` drives the two
CLIs against **each other**: `lsdtt-cosmo-prediction` forward-models a
concentration for a known erosion rate and writes the points CSV it advertises as
ingestible; `lsdtt-cosmo-tool` inverts that CSV; the imposed rate must come back.

This is the only test that can catch *plumbing* faults — a column renamed on one
side, a knob honoured asymmetrically, a shielding term applied twice. Every one
of those produces perfectly plausible numbers and passes all single-CLI tests.
Seventeen round-trip cases apply the same knob to both sides (rate sweep, three
nuclides, five muon schemes, both scaling schemes, aged and unaged, four
shielding combinations), plus checks on the density unit-conversion identity, the
standardisation factor, and ±1σ bracketing (more nuclide means slower erosion, so
the +σ rate must be the *lower* one — a sign error no single-value test can see).

Worst recovery error is **6.2e-4 relative**. The floor is not the inversion
solver but the **float32** accumulated-concentration raster the forward model
hands over, at ~7 significant digits.

Two design points keep it honest. A knob-effect self-check asserts every knob
actually *moves* the forward concentration (measured 2.2 % to 628 %), so no case
can pass merely because its knob is ignored on both sides. And `atmosphere` is
**deliberately not covered**: `lsdtt-cosmo-tool` declares it but
`lsdtt-cosmo-prediction` does not, so no matching atmosphere can be imposed on
both sides. Setting it on one side only closed to 0.38 % on this near-sea-level
fixture — which looks like a pass and established nothing. `atmosphere` was
added to `lsdtt-cosmo-prediction` on 2026-07-27 and all three settings are now
covered on both sides. An unknown atmosphere name is now **fatal** in both CLIs
rather than silently falling back to ERA40: since the atmosphere sets the
pressure at every pixel, a typo used to bias every rate in a run while looking
entirely normal. The documented ERA40/NCEP → standard fallback for *absent
grids* is unchanged — that is a data condition, not a typo.

#### Column and raster must be the same physics

The goal for `lsdtt-cosmo-prediction` is to model several isotopes responding
together to a change in erosion rate, **per pixel**. That only works if the 1-D
column model and the raster step-change kernel are the same physics driven two
ways — otherwise nothing prototyped or inverted in a column transfers to the
grid. They are separate implementations, so it has to be measured, not assumed.

Tests X1 and X3 drive both paths over three nuclides, four times, two shielding
values and rates spanning two orders of magnitude. Worst relative difference:
**5.96e-08** for steady state, **5.78e-08** for the step change — the float32
floor of the raster's own storage, i.e. exact agreement. The step is compared
away from `t = 0`, because the raster documents `t <= 0` as "no change" while the
column evaluates the step there; equating them would erase a real documented
difference that another test asserts.

**They also share one production-scaling convention, as of 2026-07-28.** The
raster folds production scaling and shielding *together* through the Vermeesch
redistribution — including the `S_tot > 1` branch that reproduces CRONUS v3 at
altitude — and the column now does the same. It previously applied production
scaling as a flat `P0` multiply, which put it −14.9 % out at a scaling of 0.5 and
+5.4 % at 1.4, and a guard meant factors above 1 were silently ignored, so the
column could not express altitude enhancement at all. Both are gone: one factor,
one mechanism, verified equal to the raster's own values. The test asserts the
sharp form — production scaling and shielding are **interchangeable**, because
only their product reaches the physics.

#### Arbitrary erosion histories

`run_transient` is the path a pixel prediction would drive for anything beyond a
single step: a column of particles carried through a user-supplied erosion
history, with instantaneous removal events and optional inherited inventories.

It is **not** given a frozen reference file of its own. A two-row history
reproduces `step_change` to **1.09e-15** — the evolution loop degenerates to the
single call `step_change` makes, from the same starting depth — and `step_change`
is already pinned at full double precision by the numerics sentinel. So
`run_transient` is pinned transitively, against the kernel the raster tools use,
rather than merely against its own past. A second frozen file would have been the
weaker guarantee.

The rest are identities no frozen file could express. Subdividing a constant
phase into 2, 5 or 20 sub-steps must change nothing (measured 2.01e-16) — that is
what catches an off-by-one in the step loop or a depth decrement applied at the
wrong point. A constant history must reproduce steady state at every sampled
depth. The particle must return to exactly its sampling depth, so the record's
erosion and the mass-removal bookkeeping cancel. Removals are covered twice,
because they take two different code paths depending on whether they sit at the
most recent row or inside the record.

**Unordered records are rejected, as of 2026-07-28.** The history validator
originally checked the background-row sentinel, rate/removal exclusivity and
removal timing, but not that the record was *ordered in time*. Index 0 must be
the most recent step; a reversed record fed negative timesteps into the solver
and returned a confidently wrong number — 72 357 against 178 467, 147 % apart —
with no error, no warning and no NaN. It now throws, with a message naming the
rule and how to fix the file. This deliberately rejects records that previously
ran, because what they returned was meaningless.

The check earns its keep immediately: it caught a mis-ordered fixture inside this
very test suite, which had been passing only because its erosion rate was
constant, so the negative timesteps cancelled.

#### The multi-isotope step-change raster

`lsdtt-cosmo-prediction`'s `calculate_multi_isotope_step_change` carries Be-10,
Al-26 and C-14 through the *same* step change at every pixel and writes the two
isotope-ratio rasters alongside the three concentrations.

It is not a wrapper around three runs. All three come from **one particle
evolution per pixel**, with the three F arrays scaled independently — which
matters for correctness, not just cost: a ratio is only meaningful if both
members saw an identical scenario, and computing them together makes that true by
construction rather than by careful parameter matching. The cost is close to a
single-nuclide run.

The load-bearing test is **bitwise equality**: each nuclide from the multi kernel
must equal the single-nuclide kernel exactly — not to a tolerance — across four
shielding values, five times and three nuclides. Anything less would mean a
second implementation of the same physics, which is precisely the situation that
produced the column-versus-raster convention bugs. The remaining checks confirm
the three really are one scenario (the ratios evolve, and C-14's moves further
over the same window), that all three masks agree so a ratio can never divide a
valid pixel by a nodata one, and that the documented `t <= 0` branch still
returns steady state at the new rate.

Ratios are formed at the output layer rather than in the library, so the nodata
policy for a zero denominator lives with the writer; pixels with zero Be-10 are
reported and written as nodata rather than silently producing infinities.

#### Arbitrary erosion histories, per pixel

`calculate_multi_isotope_transient` is the general case of which the step-change
mode is the two-row special case: the three nuclides carried through a
user-supplied erosion history at every pixel, with the history uniform in time
and the production scaling and shielding varying per pixel.

The history walk is **not** reimplemented in the raster kernel — each pixel calls
`lsdtt_crn_column::run_transient`, the same code the column tool drives. A raster
kernel with its own copy of the walk would be a second implementation of the same
physics, which is the failure mode this codebase has already been bitten by
twice.

The load-bearing test exploits the special-case relationship: a two-row history
must reproduce the multi-isotope step-change raster, which is itself
bitwise-pinned against the single-nuclide kernel, which the sentinel pins.
Measured agreement is **bitwise zero** across three nuclides, two shielding
values and three times — so this path inherits the existing pinning rather than
needing a frozen file of its own. A separate check confirms a three-phase history
gives a materially different answer from any single step, without which the
general kernel would be doing nothing the step kernel could not.

The history is validated **once, before** the parallel region, for two reasons: an
exception escaping an OpenMP loop is undefined behaviour, and a malformed history
should fail immediately rather than after most of a long raster has been
computed. A time-reversed history is rejected at the CLI boundary with the
message naming the ordering rule.

Cost is `O(pixels × history rows)` and therefore much higher than the
step-change kernel. That is inherent — an arbitrary history has no closed form —
and the CLI says so when the mode runs.

#### Two-nuclide diagnostics: burial, and inverting a ratio pair

One nuclide gives one number, so it resolves one unknown — an apparent erosion
rate, assuming steady erosion and no burial. A *pair* with different half-lives
gives two numbers and can resolve two unknowns. `lsdtt_crn_diagnostics` provides
the two cases that matter: burial/inheritance from a fast/slow pair, and
(erosion rate, time since change) from a step-change pair.

These diagnostics are now exposed on `lsdtt-cosmo-column-tool` in both single
and batch forms (`burial_pair`, `step_pair`, `burial_pair_batch`,
`step_pair_batch`). Pair choice is generic: users can choose fast/slow members
from `Be10`, `Al26`, `C14`. Solver status is surfaced explicitly as
`status`/`status_code`; non-identifiable outcomes (`not unique`, `no solution`)
are first-class outputs, never silent drops.

These **invert** the forward models the rest of the suite pins, so the decisive
test is a **round trip** — forward-model at known inputs, invert, require the
inputs back. That is much stronger than freezing outputs: a frozen value only
proves the answer has not changed, whereas a round trip proves the inverse
actually inverts *this* forward model. Measured worst errors: **5.2e-16** on the
pre-burial erosion rate, **9.3e-15** on the burial age, **6.4e-15** on the new
erosion rate and **1.0e-14** on the time since change.

Both solvers avoid a two-dimensional Newton iteration, which would be fragile
here. Burial reduces to **one** dimension because the burial time can be
eliminated analytically from the ratio equation. The step-change pair is
genuinely two-dimensional but is solved as *nested* one-dimensional problems: for
a trial time, the slow nuclide is inverted for the rate that reproduces it
(monotone, bracketed), leaving the fast nuclide's residual a function of time
alone.

**Non-identifiability is reported, not fitted.** The time axis is scanned to find
sign changes *before* bisecting, so more than one root is detected rather than
hidden. This matters: an Al-26/Be-10 pair cannot date a 6 kyr transient — Al-26's
0.7 Myr half-life makes it nearly insensitive on that timescale, and its residual
crosses zero twice — so the solver returns "not unique". Disabling that safeguard
makes it converge on *the true answer*, 6000 yr, which is exactly why the test
asserts the status rather than the value: being right by luck is not the same as
being identifiable, and a solver must not present a coin flip as a measurement.

The burial diagnostic also reports what a single-nuclide study would have
concluded. Ignoring 1 Myr of burial on a sample eroding at 0.005 g cm⁻² yr⁻¹
makes a Be-10-only analysis report a rate **69.5 % too fast**, because burial
removes nuclides that a steady-erosion model can only explain by faster erosion.
The direction of that bias is pinned, so the diagnostic cannot be wired up
backwards.

Pair modes also report uncertainty summaries from concentration-corner re-solves
(`fast ± σ_fast`, `slow ± σ_slow`): lower/upper bounds (`*_lo_*`, `*_hi_*`),
an internal uncertainty (`*_unc_int_*`, half-range), and an external
uncertainty (`*_unc_ext_*`, internal plus optional production/decay terms in
quadrature). The number of converged corners used is reported (`*_unc_n_valid`)
so a partially non-convergent corner set is visible rather than hidden. The
`unc_int`/`unc_ext` naming matches the CRONUS internal/external convention used
throughout this test regime.

#### Multi-isotope transients, and two expectations the data corrected

Isotope *ratios* are where a transient signal lives — a nuclide silently sharing
another's decay constant would leave every individual concentration looking
plausible. Test X4 asserts monotone relaxation and endpoint bounding for all
three nuclides, both steady-state ratio limits, and the ratio's transient
behaviour.

Writing it corrected two plausible-sounding expectations, both worth recording
because they are easy to repeat:

- **Relaxation is not ordered by half-life.** Measured at 8 kyr after a tenfold
  erosion increase: C-14 0.850 of the way, Al-26 0.551, Be-10 0.562 — Al-26
  *trails* Be-10 despite the shorter half-life. At ε = 0.02 g cm⁻² yr⁻¹ the
  erosional timescale Γ_sp/ε is 8000 yr, orders of magnitude shorter than either
  0.7 or 1.39 Myr half-life, so erosion sets the relaxation for both and decay is
  only a small correction; what separates them is muon partitioning. C-14 leads
  because its 5.7 kyr half-life is the only one comparable to the erosional
  timescale.
- **The C-14/Be-10 ratio does not overshoot its final value.** Erosion increased,
  so every concentration falls; C-14 reaches its new low value first while Be-10
  is still draining from the old high one, so the ratio is *suppressed* and
  approaches its endpoint from below (0.332 → 0.452 at 8 kyr → 1.831). The
  datable signal is precisely that lag — C-14 is 85 % relaxed while the ratio has
  covered about 8 % of its range — not an excursion past the endpoint.

#### The column tool

`lsdtt-cosmo-column-tool` is the 1-D case: one column, no DEM, no flow routing,
no basin averaging. That makes it the natural place to compare lsdtt3 against
point calculators such as CRONUS v3 and CRONUScalc, because nothing in the
pipeline can mask a scaling or muon discrepancy.

`rigorous_cosmo_column_tool` covers it with limit identities that need no
external reference — cases where the physics forces the answer, so a wrong one
is unambiguous. A step change evaluated at t = 0 must equal steady state at the
*old* rate, and at t → ∞ the *new* rate. A constant-rate transient run must
relax to the steady state the first mode computes: two independent code paths
meeting at the same number, measured to 1e-8. Concentration must fall
monotonically with depth and with erosion rate. C-14 must re-equilibrate faster
than Be-10 after a step change — the reason C-14 is interesting at all.

Two checks pin *which mechanism each knob uses*, not merely that it does
something: `production_scaling` must be exactly linear in `P0`, while
`total_shielding` must **not** be, because it goes through the Vermeesch
virtual-attenuation redistribution that attenuates the pathways unequally. A
silent swap of the two therefore fails, even though both would still move the
answer.

**A defect this found, and how it had hidden.** The two point modes originally
applied `total_shielding` by different mechanisms — `steady_state` through the
Vermeesch redistribution, `step_change` as a flat `P0` scale — so they disagreed
by 11.4 % at a shielding of 0.5, in a case where the limit identity says they
must agree exactly. It survived because the identity was only ever checked at a
shielding of 1, where the two mechanisms coincide. A test that exercises an
invariant at a single point often cannot see the fault it was written for.

Fixed by giving `lsdtt_crn_column` **additive** overloads that take the site
production scaling and the shielding separately and apply each by its own
mechanism; the existing single-factor signatures are untouched, so no other
caller can move, and both numerics sentinels stayed bit-identical. The CLI now
calls the library in both modes rather than one of them rebuilding the parameter
object inline — which is how the two drifted apart to begin with. The identity is
now asserted across four `(shielding, scaling)` combinations including a scaling
above 1, so the altitude-enhancement branch is covered too.

#### Where the expensive parts run

The per-pixel side of (a) needs real DEMs, so it follows the §7.7 pattern: a
`cosmo_realworld` script (`run_cosmo_validation.py` / `run_cosmo_lsdn_pixel_vs_lumped.py`)
fetches keyless COP30, runs 50 basins, and freezes the
per-basin means into a small tracked CSV; the default `unit` suite then does the
cheap comparison against that frozen artefact. The frozen artefacts are
**check-only by default** and regenerate only under an explicit `--freeze`. A
missing artefact **fails** — it never reports a skip, because the files are
tracked and their absence means a broken checkout, not an optional step.

---

## 8. Implicit hillslope solver validation

The backward-Euler + damped-Picard solver for the Roering (1999) nonlinear
hillslope flux law (`compute_roering_hillslope_transport_implicit`) is validated
by the unit-test driver `tests/cpp/test_hillslope_implicit.cpp` (CTest target
`unit_hillslope_implicit`, label `unit;hillslope`), plus in-LEM
hillslope-strip runs driven entirely from `lsdtt-lem` parameter files. The driver
takes its output directory as `argv[1]` and emits the plottable `*.csv` there —
under CTest that is `build/tests/hillslope_implicit_output/`. The CSVs are
diagnostics, not assertions: an unwritable directory warns and the run continues.

The same binary carries the Phase 4 linear-solve **timing sweep** behind
`--performance`; it is registered separately as `performance_hillslope_implicit`
(label `performance;hillslope`, gated on `LSDTT_ENABLE_PERF_TESTS`) and is
reported, never asserted, so it stays out of the default suite.

**Test domains are generated in C++, never in Python.** The hillslope-strip
generator in `lsdtt_raster_synthetic_surfaces_cpu` (exposed through
`lsdtt-raster-synthetic-tests`) returns `{dem, channel_mask}` with flat, analytic
parabola and raised-plateau initial profiles; the network tests use the
`diamond_square` initial surface already built into `lsdtt-lem` and let the
model's own stream-power incision produce the branching network. Python supplies
only analytic reference solutions and diagnostic plots — never input DEMs.

### 8.1 The analytic oracle

The primary reference is the **Roering (1999) nonlinear steady state**, not the
linear parabola. Mass balance gives flux `q_s(x) = U·x` from the divide, and
`q_s = D·S/(1 − (S/Sc)²)`, so

```
S(x) = (Sc²/(2Ux))·[√(D² + 4U²x²/Sc²) − D]
```

integrated to `z(x)`: a rounded crest with planar mid-slopes and `max|S| < Sc`.
The linear parabola `(U/2D)(L² − x²)` is **only** the `Sc → ∞` limit and is used
solely for the linear-mode check. Choose `D`, `Sc` and domain length so the
channel-end slope `U·L/D` is a good fraction of `Sc`, or the nonlinear term is
never exercised.

Domains are a few columns wide with **periodic column edges**, so the columns
must stay identical to machine precision — that is the cross-slope/lateral bug
check, and it is the reason for using more than one column. Two boundary
scenarios: *Full* (`"opop"`, both end rows fixed, symmetric hill, central divide)
and *Half* (`"opnp"`, one fixed end and one no-flux end, divide at the no-flux
end) — the latter being the no-flux-BC check. Fixed ends are imposed through the
`'o'` boundary code, not through `is_channel`.

Note that `'o'` means **fixed elevation** for the implicit solver, whereas the
existing SOR and explicit hillslope solvers special-case only `'p'` and treat
both `'o'` and `'n'` as Neumann. The SOR comparison therefore restricts itself to
`'n'`/`'p'` edges with `is_channel` Dirichlet sinks rather than extending the SOR
edge handling — the lower-risk of the two options, and it works.

### 8.2 What is validated, and the measured numbers

| Test | Check | Result |
|---|---|---|
| 1 | Roering steady state, Full and Half scenarios, under grid refinement | passes at **2nd-order** convergence in `h` |
| 2 | Step change in uplift (U 5×10⁻⁴ → 1×10⁻³) — the field must land on the Roering steady state for each rate | strongly supercritical driving (`U·L/D = 2.5·Sc`, 40 m domain) lands on Roering(U₂) at Linf/peak ≈ 4.1e-4, `max_slope_ratio ≈ 0.816`, 51 rejected sub-steps |
| 3 | Coupled lowering-channel dynamic equilibrium — the Dirichlet hand-off check | new-time relief matches the Roering `(U+E)` profile to ~2e-5 and is `dt`-independent |
| 4 | Match the nonlinear SOR steady state on the same DEM and `is_channel` | agrees to the SOR float32 floor: linear rel **3.9e-6**, nonlinear rel **2.4e-6** |
| 5 | Linear-mode contrast (`nonlinear = false`) against the exact parabola | **NOT IMPLEMENTED** — no such check exists in the driver. The row describes intent, not a result. The only linear-mode check present is test 4's linear half, which compares the implicit solver against SOR, not against the analytic parabola |
| 6 | Timestep independence across `dt_init = tau/{10, 40, 160}` | rel ≈ **2.4e-6** |
| 7 | Conservation (U = 0, `'nnnn'`, no sinks) — the discrete volume budget must close | max relative volume drift **2.3e-8** |
| 8 | Performance benchmark, one implicit solve, flat, U = 0 | see §8.3; run via `--performance` as `performance_hillslope_implicit`, not part of the default suite |
| 9 | `Sc` enforcement stress (U = 2×10⁻²) | caps at `max_slope_ratio = 0.99` with `limiter_fired` set and bounded rejections |
| 12 | TR-BDF2 temporal accuracy | reaches Roering(U) at Linf/peak 5.0e-4 with **observed order 2.00** against 0.96 for backward Euler — ~35× more accurate at the coarsest `dt` |

A **backend-agreement check** compares `SimplicialLDLT` against
`ConjugateGradient` + `IncompleteCholesky`: the steady fields agree to **0.0** —
identical operator, different solver.

On test 9, note that a large per-step channel drop instantly over-steepens the
channel-margin face, so the `r_max` cap is *expected* to fire on Picard iteration
0 there and relax as the solve converges. That is benign, not a limiter failure.

Two network-scale tests remain outstanding: a dendritic-network domain
(hillslopes only, on a frozen fluvial-steady-state channel mask) and a fully
coupled fluvial + hillslope diamond-square run. Both run entirely from
`lsdtt-lem` parameter files.

### 8.3 Solver-backend benchmark

One implicit solve, flat field, U = 0, run as back-to-back interleaved A/B pairs
on the same machine — WSL2 host-load noise makes isolated timings unreliable.
Run it with `ctest -L performance -R hillslope` after configuring with
`-DLSDTT_ENABLE_PERF_TESTS=ON`, or directly as
`lsdtt-test-hillslope-implicit <out_dir> --performance`.

| Grid | LDLT | CG + IncompleteCholesky |
|---|---|---|
| 128² | 43 ms → 42 ms | 20 ms → 13 ms |
| 256² | 332 ms → 228 ms | 62 ms → 43 ms |
| 512² | 2554 ms → 2043 ms | 265 ms → 214 ms |

LDLT scales roughly as O(N^1.5) — the factorisation cost — while CG+IC is near
O(N), which is the empirical confirmation of the LDLT→CG crossover. The
arrow-pairs are before and after the assembly/scan rework, a consistent ~1.3× for
LDLT and ~1.2–2× (noisier) for CG+IC, at **identical results**: the whole suite
passed unchanged, with sentinel values matching (test 2 strong 4.103e-4, ratio
0.8157, 51 rejects; conservation drift 2.27e-8; SOR agreement at the float32
floor; TR-BDF2 observed order 2.00; column spread 0), and the in-LEM strip run
reproducing its recorded steady peak (4.943 vs 4.94).

### 8.4 Cautionary findings

**The semi-implicit omega bug — the case for an integration test.** The
semi-implicit fast path applied the *damped* Picard `omega` to its single IMEX
solve, under-relaxing and therefore under-diffusing the step. This was **invisible
to the unit-test driver**, whose steady-state tests all finish on full Picard, so
every unit-level check passed — while coupled LEM steady states were biased about
**40 % high**. It was only exposed by the in-LEM integration run. A lone
semi-implicit solve now uses `omega = 1`.

> A fast path that is only taken in regimes the unit tests do not finish in is a
> fast path the unit tests cannot validate. Where a code path exists specifically
> for the cheap case, an end-to-end test that actually spends time in it is not
> optional.

That end-to-end test is now registered and runs in the default suite:
`smoke_lem_strip_fixture` generates the strip through
`lsdtt-raster-synthetic-tests`, and `smoke_lem_roering_implicit_strip` runs
`lsdtt-lem` on it (`tests/params/lem_roering_implicit_strip.param.in`) and pins
the converged steady state two-sided — peak **4.9429 m** and mean **3.2719 m**,
the recorded in-LEM strip values. `U·L/D = 0.2` is a quarter of `Sc`, so the run
spends essentially every step inside the semi-implicit fast path. Reinstating the
omega bug moves the peak to **6.862 m (+38.8 %)** and the mean to 4.516 m, both
far outside the bands; and because the strip starts **flat**, a solver that does
nothing leaves it at 0 m rather than somewhere plausible. Both tests together
take ~1.5 s.

**The explicit solver's CFL estimator was wrong — and silently so.**
`compute_roering_hillslope_transport` sub-cycled on `dt_stable = h²/(4·D·α_max)`
with `α = 1/(1 − r²)`, `r = |∇z|/Sc` — the *secant* diffusivity from the flux law.
The coefficient that actually governs explicit stability is the **tangent**
modulus `dq_s/dS = D·(1 + r²)/(1 − r²)²`, which diverges faster by a factor
`(1 + r²)/(1 − r²)` — **~13× at r = 0.93**, and unbounded as `r → 1`. Sub-steps
were therefore far above the true limit on near-critical slopes: the step
overshoots, drives faces supercritical, and runs away through the singularity.

It surfaced only when the explicit model was used as the *reference* for the
dendritic-network test: at `D = 1e-2` with max `S/Sc` ≈ 0.93 it blew up at
`dt = 500 yr` (single-step elevation changes of 85–170 m, `delta_max_z` →
2×10⁴ m) *despite* sub-cycling, while `roering_implicit` sat at
`delta_max_z = 0.02 m` with `max_erosion_rate` exactly U. With the tangent modulus
in place the explicit model is stable again, at a much smaller — i.e. honest —
sub-step. Perron (2011) makes the same point about the linearised coefficient.
The implicit solver was never affected; its `r_max` cap, damped Picard and
post-step rejection are built around the singularity.

> The reference path deserves the same scrutiny as the path under test. A
> reference that fails loudly is a nuisance; one that is quietly too permissive
> corrupts every comparison drawn against it.

**Float versus double.** Assembly and solve are in double for conditioning;
raster I/O stays float. Test 4's agreement with the SOR solver bottoms out at the
SOR float32 floor, which is why the pinned tolerance sits there and not lower.

---

## 9. Credentials in tests

`lsdtt-fetch-raster` is the only lsdtt3 tool that touches secrets, and the test
suite encodes the rules that keep them out of the repository.

**An API key must never be a parameter.** Every lsdtt3 CLI is parameter-file
driven, users routinely commit their `.param` files alongside analyses, and
`print_parameters()` echoes every parameter at the end of a successful run. A key
in a parameter file would be committed *and* printed. The tool therefore refuses
to read a key from the parameter file at all — this is a structural rule, not a
documentation convention. Nor may a key arrive in argv, which is visible to any
user on the box via `ps` and lands in shell history.

Two tests enforce this:

| Test | What it asserts |
|---|---|
| `rigorous_fetch_raster_no_credential_parameters` (`rigorous;cli;fetch`, offline) | Scans the CLI's own help output and fails if any parameter is named like a credential (`api_key`, `password`, `token`, `secret`). This is what stops a future change quietly introducing a key parameter. |
| `test_no_credential_leakage` (in `unit_raster_fetch`) | Asserts nothing credential-shaped (`API_Key`, `Bearer`, `password`, `token`) appears in any loggable result field. It passes trivially today; its job is to fail the day someone routes a key through `lsdtt_fetch_options`. |

`lsdtt_fetch_options::api_key` exists but is filled **at the call site** from
`lsdtt_secret`, never from a parameter — the header says so at the field.
Resolution order is `OT_API_KEY` (environment variable) → `OT_API_KEY_FILE`
(environment variable naming a file) → TTY prompt with echo disabled, skipped
entirely when stdin is not a terminal so CI fails cleanly instead of hanging.
`src/core/lsdtt_secret.hpp/.cpp` is the only place in lsdtt3 that handles
secrets, which makes it the only place to audit; it has 24 offline checks in
`tests/cpp/test_secret.cpp`.

**A missing key is an error, not a skip.** This is the §1 rule applied
specifically: a skip would hide a misconfiguration behind a "test skipped"
result. Only an unavailable *dependency* — no network — may skip.

**Masking is not theoretical.** Run against the live OpenTopography API with a
deliberately fake key, the provider **echoes the submitted key back inside its own
error message**:

```
COULD NOT FETCH: ... the OpenTopography API returned an error:
<error>Error: Not a valid format API Key: ***. Please register ...</error>
```

Without `mask_secret()` every mistyped or expired key would be printed to the
terminal and into any redirected log. GDAL error messages embed the failing URL
for the same reason, so `quiet_gdal_errors` is the masking hook point. Planetary
Computer's SAS tokens are credentials too — time-limited, in the URL query string,
granting blob read access even though obtaining one needs no account — and
because the signature is generated per request there is no fixed string to mask;
`mask_url_query()` replaces the whole query with `?***`. Verified that no SAS
material (`sig=`, `st=`, `se=`) appears anywhere in the CLI's output.

Two practical notes. The key-file permission warning (group/world-readable) is a
warning, not a failure, and is skipped on filesystems where `chmod` is a no-op —
NTFS under WSL, where this build tree lives. The test detects that and reports
SKIP; on a real POSIX filesystem all 24 checks assert and pass. And the
OpenTopography **happy path has not been verified** — no key was available in the
environment. The error paths, resolution order, masking and transport are all
exercised; `cop30_aws` remains the only source proven to return real data.

---

## 10. Fetch cross-validation

Three routes to the same 128²/256² alpine window — explicit `west/south/east/north`,
a polygon in UTM 32N reprojected on the way in, and two points on the *reversed*
diagonal — produce a **byte-identical raster**, sha256 `899cc313e27cfcab…`. So do
the C++ library path, the CLI, and the now-deleted `fetch_test_dem.py`, all on
GDAL 3.4.1, and so does the `planetary_computer` source against `cop30_aws` (the
same underlying Copernicus data by different routes, which makes it a useful
independent check on both). This byte-identity is what justified retiring the
Python fetcher without re-pinning anything — *provided the GDAL version is held
constant* (§3).

Bounding-box handling is deliberately unforgiving: a missing bounding-box edge is
a hard error, because a silently zeroed box fetches the Gulf of Guinea rather
than failing (`failure_fetch_raster_missing_bbox`). Giving both `bounds_fname`
and any of the four edges is likewise a hard error rather than a precedence rule,
so the user is never left unsure which area was fetched. Two corner points need
no orientation autodetection: both diagonals give the identical box, since it is
just min/max of the two points either way, and a test asserts the two
arrangements produce bit-identical bounds.

Live timings, for reference: a 128² alpine fetch completes in **1.95 s**, a 256²
Ortler window through the CLI in **2.3 s**, and a 256² Planetary Computer fetch
in **2.8 s**.

### Bugs the end-to-end runs caught that unit tests could not

- **`GDALAllRegister()` in the wrong place.** It was called just before the fetch,
  but reading bounds from a vector file happens during parameter resolution — so
  no OGR drivers were registered and every vector file failed to open. The unit
  tests missed it because the test driver registers drivers itself; only the
  end-to-end CLI run exposed it.
- **The packaged bundle.** `package_linux_backend.sh --cli-profile core` stages
  libraries by SONAME. Running the binary out of the *extracted bundle* and
  checking it produces the same byte-identical raster was worth doing explicitly:
  a missing `libcurl.so.4` would have left `/vsicurl/` silently broken for every
  user of a release bundle while working perfectly on the build machine.
- **Packaging gotcha.** `resolve_exe_path()` hard-exits 1 when a CLI listed in
  `tools/cli_profiles.txt` is absent from the build directory, and that file feeds
  all three packaging scripts plus `make_release.sh`. A profile entry must never
  be added before the CLI exists — doing so breaks every release build
  immediately, on all platforms, with no partial-success mode.
- **`CPL_VSIL_CURL_ALLOWED_EXTENSIONS=.tif`** would have blocked *every* Planetary
  Computer asset, because their hrefs carry no file extension
  (`.../Copernicus_DSM_COG_10_N46_00_E010_00_DEM`). The option was removed;
  `GDAL_DISABLE_READDIR_ON_OPEN` already provides most of the latency win.
  Relatedly, unsigned access to Planetary Computer assets returns HTTP 409 — the
  signing endpoint needs no account but *is* required.

---

## 10b. Surface metrics (`lsdtt-surface-metrics`)

### How correctness is established

**Not** by comparison with another tool. The oracle is analytic:
`generate_global_quadratic_test_surface()` and
`generate_periodic_wave_test_surface()` return exact local Taylor coefficients
at every pixel, and on a globally quadratic surface the fitted model matches the
true surface class exactly. Since all 27 derived metrics are algebraic functions
of the six coefficients, this validates the whole chain. Sphere-cap, cylinder,
saddle and tilted-plane surfaces with textbook closed-form curvatures give a
second, independent oracle.

`tests/cpp/test_surface_metrics.cpp` runs 124 checks in ~1 s, offline, no
committed data (`unit;surface`). Note that **23 of the 27 metrics have no
TopoToolbox counterpart at all**, so an external reference could not have
covered them.

### Measured accuracy, and why 1e-6 is not the right target everywhere

| Quantity | Measured worst | Pinned tolerance |
|---|---|---|
| Quadratic, 2nd-derivative coeffs (a/b/c), W=1 | 1.05e-3 | 5e-3 |
| Quadratic, same at W=5 | 9.6e-6 | 5e-5 |
| Quadratic, 1st-derivative coeffs (d/e) | 2.5e-6 | 1e-5 |
| Sphere / cylinder curvatures | 3.5e-5 – 7.0e-5 | 2e-4 |
| Plane, all curvatures (exact zero) | 5.3e-15 1/m | 1e-6 |
| Algebraic identities (mean, Gaussian) | 9.4e-8 / 1.45e-7 | 1e-5 / 1e-4 |

The second-derivative coefficients cannot reach 1e-6 from a float32 DEM at W=1,
and this is arithmetic, not a defect: `z_xx` is recovered from elevation
differences over `(W*L)^2`, so a float32 quantum `eps*Z` propagates to
`eps*Z/(W*L)^2` — with Z~1e3 m, W=1, L=10 m that is ~1e-3 relative. At W=5 it
falls to ~1e-5, as predicted. First-derivative and constant coefficients do meet
1e-6. Tolerances are pinned to measurement with ~5x headroom.

### Angular conventions (project-wide, enforced 2026-07-22)

**Every bearing, aspect and orientation lsdtt3 emits uses the COMPASS
convention: 0 degrees = NORTH, increasing CLOCKWISE** (90 = east, 180 = south,
270 = west). This matches GDAL, ArcGIS and the geomorphology literature.

It is **not** the mathematical convention (0 = +x axis, anticlockwise), which
in map view puts 0 at **east**. To convert: `compass = fmod(90 - maths + 360,
360)`. Note `atan2(y, x)` is the maths form; a compass angle is
`atan2(east, north)`, arguments in that order - getting them the wrong way round
is the classic defect and accounted for three of the five found here.

**Aspect** additionally means the **down-slope** azimuth - the direction the
slope FACES - i.e. `atan2(-p, -q)` with p = dz/dx east-positive, q = dz/dy
north-positive. `atan2(p, q)` is the up-slope azimuth and is 180 degrees wrong.

**Axial** quantities (line orientations with no head or tail: MADk2 DMC, basin
principal axis, skeleton dominant orientation, feature azimuths) use the same
zero and sense but fold to [0, 180).

This note is repeated in the headers of every module that emits such a value.

#### The five defects fixed on 2026-07-22

An audit of every angular expression under `src/`, `cli/` and `tests/` found
five wrong sites. All are fixed; each is guarded by a test or verified against
a known-answer surface.

| # | Site | Was | Now |
|---|---|---|---|
| 1 | `aspect_radians` / `aspect_degrees` | `atan2(p,q)` - up-slope, so 0 meant "faces south" | `atan2(-p,-q)`, 0 = faces north |
| 2 | `compute_hillshade` sun azimuth | `A = az - pi/2`, lighting as though the sun were at az+90 | `A = az - pi` |
| 3 | MADk2 `DMC` | `+pi/2` applied, emitting the compass bearing of the maximum-**roughness** axis | halved resultant directly; now the maximum-**continuity** axis in compass axial degrees |
| 4 | `lsdtt-flat-surface-finder` `aspect_degrees` | `atan2(b,a)` - maths frame **and** up-slope | `atan2(-a,-b)`, [0,360) |
| 5 | basin `orientation_deg`; skeleton `dominant_orientation_degrees` | maths frame; and `atan2(dc,dr)` mirrored (NE-SW read 135 where compass says 45) | compass axial [0,180) |

Sites 1 and 2 shared a root cause: the textbook hillshade formula
`cos(az - pi/2 - aspect)` is only valid when `aspect` is the down-slope azimuth,
and lsdtt3's was up-slope. Fixing either alone would have left them
inconsistent.

Verified correct and unchanged by the audit: streamline steepest-descent
bearings (the reference implementation in the repo), channel-link and
centre-of-mass bearings, junction angles, D-infinity flow direction,
cross-section and tidal tangent/normal bearings, linear-feature and spectral
azimuths, and basin-extraction fault geometry.

**Compatibility.** Aspect rasters from earlier versions need rotating by 180
degrees; hillshades were lit from 90 degrees off the requested azimuth; DMC
rasters need rotating by 90 degrees; flat-surface aspect and basin
`orientation_deg` change frame entirely.

### Conventions that differ from other tools

### Conventions that differ from other tools

Three traps for anyone comparing lsdtt3 against ArcGIS, GDAL or TopoToolbox:

1. ~~**Aspect is the azimuth of the UP-slope direction**~~ **FIXED 2026-07-22.**
   `aspect_radians` / `aspect_degrees` are now `atan2(-p, -q)` - the **down-slope**
   azimuth, clockwise from north, the direction the slope FACES. 0 = faces north.
   This is the GDAL / ArcGIS convention and the universal meaning of "aspect" in
   the geomorphology literature (Smith & Bookhagen 2020 calls `gdaldem Aspect`
   directly; Pelletier et al. 2018, Poulos et al. 2012 and Yetemen et al. 2015
   all use pole-/equator-facing framing that presumes it).

   Before the fix, 0 meant "faces south", so anyone binning lsdtt3 aspect as
   north-facing got exactly the opposite hemisphere - silently, with plausible
   output. Verified after the fix on planes descending toward 0/45/90/180/270
   degrees: aspect returns the descent azimuth exactly in every case.

   The only internal consumer, `lsdtt_raster_feature_finder_cpu.cpp`, added 180
   degrees to recover a downslope bearing; that rotation was removed in the same
   change. `lsdtt_streamlines_cpu.cpp` computes its own `atan2(-p,-q)` descent
   bearing and was always correct.

2. ~~**Hillshade azimuth is offset by 90 degrees from TopoToolbox**~~
   **THIS WAS A BUG IN lsdtt3, NOT A CONVENTION DIFFERENCE. FIXED 2026-07-22.**
   An earlier version of this section recorded the measured relation
   `ttb_azimuth = lsdtt3_azimuth + 90 (mod 360)` as a difference of convention.
   It was not: `compute_hillshade` used `A = az - pi/2` where the physics gives
   `A = az - pi`.

   Illumination is the dot product of the unit normal `n ~ (-p, -q, 1)` with the
   unit sun vector `s = (sin(az)cos(alt), cos(az)cos(alt), sin(alt))`, giving
   `hs ~ sin(alt) - cos(alt)*(p*sin(az) + q*cos(az))`, which requires
   `sin_A = -sin(az)`, `cos_A = -cos(az)`. The old code lit the terrain as
   though the sun were at **az + 90 degrees** - which is precisely the measured
   offset, so the benchmark had detected the defect and it was mis-filed as a
   convention. The root cause was the base formula `cos(az - pi/2 - aspect)`,
   which is only valid when `aspect` is the down-slope azimuth; the aspect in
   use was the up-slope one. Both errors had the same origin.

3. **Curvature sign**: convex-up (ridge, peak) gives NEGATIVE mean curvature.
   This is not the Zevenbergen & Thorne convention.

Also: on low-relief terrain lsdtt3 masks substantially more ground as flat.
On a Lakeland fixture the direction-dependent metrics gave lsdtt3 246,450 valid
cells against TopoToolbox's 518,400 — 52% fewer. Comparisons there are not
like-for-like and the harness warns above a 10% asymmetry.

### Speed vs TopoToolbox

**TopoToolbox is used here as a performance reference only, not a correctness
one.** `evansslope` and `curvature` use the same Evans 3x3 kernels lsdtt3's
polyfit uses at W=1 — verified from source — but they are **single-threaded
`scipy.ndimage.convolve` calls**, so the 1-thread column is the honest
comparison. `aspect` is plain `np.gradient` central differences, which does less
work than a quadratic fit, and is report-only.

alpine_medium (720x720 @ 30 m), W=1 forced via `polyfit_lambda_metres: 60`,
best-of-5 after a warm-up, both sides float32, native ext4.

**These ratios are NOT publication-grade on this hardware.** Re-running the
identical sweep moved TopoToolbox's *own* reference times by up to 4x between
runs (slope 10.1 -> 40.4 ms; curvature_all4 98.7 -> 159.8 ms) on a 15 W laptop
part that throttles under sustained load. The lsdtt3 absolute times below are
stable and were cross-checked against controlled in-process microbenchmarks;
**the ratio column should be treated as indicative only** until the sweep is
re-run on a machine that is not power- and thermally-limited.

Measured after the 2026-07-22 optimisations (`direct_polyfit` column blocking,
2.07x; `compute_hillshade` transcendental removal, 5.2x):

| comparison | lsdtt3 1t | topotoolbox | ratio (indicative) |
|---|---|---|---|
| slope | 35.4 ms | 40.4 ms | 0.5 - 1.0x |
| profile curvature | 31.5 ms | 27.0 ms | ~0.9x |
| mean curvature | 27.4 ms | 28.6 ms | ~1.0x |
| plan curvature | 28.1 ms | 27.8 ms | ~1.0x |
| **all four curvatures** | **30.2 ms** | **159.8 ms** | **~5x** |
| hillshade (fused, compiled) | 11.2 ms | 3.0 ms | ~0.3x |

For reference, the same table before those optimisations read: slope 41.2 ms,
all-four 37.6 ms, hillshade 18.0 ms.

**Read this honestly.** For a *single* metric lsdtt3 is roughly at parity or
somewhat slower, because the direct polyfit always solves all six quadratic
coefficients whatever you ask for; it wins decisively once you want more than
one product of the same fit (~5x on four curvatures). The hillshade row is the
only compiled-vs-compiled comparison (libtopotoolbox) and lsdtt3 still loses it
by ~3x even after the 5.2x rewrite - the remaining cost there is the per-pixel
3x3 nodata check, not the arithmetic.

Threading buys only ~1.3-1.7x on 8 cores, not 8x — worth investigating
separately.

**Memory is a consistent win**: lsdtt3 83-86 MB peak against TopoToolbox's
207-233 MB on the same fixture; 64 MB vs 208 MB for hillshade.

### Two structural properties worth knowing

**Adding a metric is nearly free.** The fit is unchanged between a 4-metric and
a 27-metric request (47.9 ms vs 41.7 ms — the difference is noise, and the
27-metric run was faster). What you pay is the raster write, plus memory: peak
RSS goes 86.4 MB to 143.7 MB, about 2.5 MB per extra metric, which is one
float32 raster of that size. TopoToolbox has no equivalent — each extra
curvature type is another five convolutions before any write.

**The direct polyfit really is O(n) independent of window size.** Fit time
across W = 1, 5, 15, 40 (9 to 6561 window cells, a 729x increase):
39.8 / 43.3 / 42.5 / 47.5 ms — a 1.19x spread. The separable rolling-moment
claim in the CLI reference is correct.


### Why hillshade is ~4x slower, and why threading buys so little

Investigated 2026-07-22 after the benchmark raised both. Neither is what it
looks like.

**Hillshade: five transcendentals per pixel that cancel algebraically.**
`compute_hillshade` evaluates, per pixel, in double precision:
`atan(sqrt(p*p+q*q))`, `atan2(p,q)`, `cos(slope)`, `sin(slope)` and
`cos(az - pi/2 - aspect)`. That is where essentially all the time goes, and it
is avoidable. With `r^2 = p^2+q^2` and `A = az - pi/2`, using
`cos(slope) = 1/sqrt(1+r^2)`, `cos(aspect) = q/r`, `sin(aspect) = p/r`, the
expression collapses to a dot product of the surface normal with the sun
vector:

```
hs = 1 + 254 * (sin_alt + cos_alt*(sinA*p + cosA*q)) / sqrt(1 + p^2 + q^2)
```

One square root per pixel, no transcendentals; the flat-pixel case `r = 0` is
handled identically by both forms. Measured on 720x720, single thread:

| | time | agreement |
|---|---|---|
| current implementation | 31.5 ms | - |
| closed form | **2.74 ms** | **max abs difference 0.000e+00 over 515,524 px** |

**APPLIED 2026-07-22.** Verified bit-identical against the previous
implementation across **252 cases** - 2 terrain types (random quadratic,
diamond-square) x 7 azimuths x 3 altitudes x 2 z-factors x 3 seeds - with zero
differing pixels and a maximum absolute difference of 0.000e+00.

In place, `compute_hillshade` went from **31.5 ms to 7.58 ms single-threaded on
720x720: 4.2x faster.** It does not reach the 2.74 ms of the bare kernel above
because the real function also runs a 3x3 nodata check per pixel (nine mask
lookups) and writes an output mask; that check, not the arithmetic, is now the
dominant remaining cost and is the next lever if hillshade ever matters more.

Note on the end-to-end benchmark: it is too noisy on this thermally-constrained
laptop to confirm the ratio cleanly - TopoToolbox's *own* reference time moved
from 4.07 ms to 2.61 ms between two runs of the identical code, and an 8-thread
figure came out worse than the 1-thread one. The controlled in-process
microbenchmark (same process, same thermal state, best-of-7) is the number to
trust. Re-measure the end-to-end comparison on a machine that is not power- and
thermally-limited before quoting a lsdtt3-vs-TopoToolbox hillshade ratio.

**Threading: this machine's ceiling is 3.3x, not 8x.** Measured on the
benchmark box (Intel i5-1145G7: **4 physical cores, 8 hyperthreads, 15 W laptop
part** that cannot hold turbo across all cores):

| kernel | 1t | 2t | 4t | 8t |
|---|---|---|---|---|
| pure compute, negligible memory traffic | 1.00x | 1.76x | 3.00x | **3.30x** |
| pure streaming (triad, ~12.4 GB/s) | 1.00x | 1.31x | 1.30x | 1.13x |
| `compute_hillshade` | 1.00x | 1.77x | 2.45x | **3.10x** |
| `direct_polyfit` | 1.00x | 1.43x | 1.57x | 0.99x |

Two ceilings, and each kernel sits on the right one. 3.30x is the compute
ceiling - hyperthreads are not cores. Memory bandwidth saturates at **2
threads**, so anything bandwidth-bound is capped near 1.3x however many threads
it gets.

`compute_hillshade` reaches 3.10x against a 3.30x ceiling: **its threading is
working correctly.** `direct_polyfit` tops out at 1.43-1.57x and degrades past
4 threads, matching the streaming profile - it is **not
parallelism-limited**, and more threads cannot help it.

**It is cache-line-utilisation bound, not bandwidth bound.** (An earlier draft
of this section said `direct_polyfit` allocated two complete sets of six
coefficient rasters. That was wrong: the `g_*` set is inside
`#ifdef LSDTT_GPU_KERNELS` behind a live-device check and is not compiled in a
CPU build. Peak is 58 B/pixel, measured - `VmHWM` delta 58.1 B/pixel at n =
1440, 2880 and 4096 against an analytic 58; two live sets would give 83.)

The cause is traversal order. Pass 2 parallelises over columns and walks *down*
a column of row-major arrays, so every 8-byte read pulls a 64-byte cache line
and every 4-byte write pulls one for ownership, across seven strided streams.
That is how ~58 bytes of data per pixel generates ~1.6 kB of traffic, and why
per-pixel cost **grows with raster size** - 104.7 / 196.5 / 321.6 ns/pixel at
n = 1440 / 2880 / 4096 - which a genuinely streaming kernel does not do. Pass 2
is 55-65% of the function; allocation is ~10% and pass 1 ~3%.

Processing columns in blocks of 8, with the prefix arrays interleaved so reads
and writes run along rows, measures **2.8x on pass 2 and ~1.7x on the whole
function at 2880x2880, bit-identical** - the arithmetic and its order are
untouched, only the iteration order over independent columns changes. Blocked
pass 2 at one thread beats the current pass 2 at four. There is no gain at
1440x1440, where the working set still fits in cache: the expected signature.

Second item: the zero-fill of `H0`/`H1`/`H2`/`Hn` (28 B/pixel, 232 MB at
2880x2880) is dead, because pass 1 writes every element; it costs a measured
5-10% of runtime. The `out_*` nodata fill is *not* dead - border and
skipped-window pixels are never written.

Re-check the 3.30x ceiling before reading any thread-scaling number measured on
a different machine.

### Defects this work found

All fixed 2026-07-22. Recorded because each was a case of "the run succeeds
while doing less than you asked", which is the failure mode tests exist for.

1. **MADk2 treated axial directions as vectors.** The four sampling directions
   are axes (period 180 degrees), so the anisotropy resultant must be formed at
   the doubled angle and halved. Using the raw angle gave RA = 0.653 on a
   perfectly isotropic surface instead of 0 - so RA never reached its own lower
   end and "low RA means isotropic" could not be read off the output - and put
   DMC 35.26 degrees out for one corrugation orientation while leaving it exact
   for the perpendicular. Guarded by the C5 checks.
2. **Parameter values were truncated at the first space**, silently:
   `polyfit_metrics: slope, aspect, mean_curvature` computed only `slope`,
   exit 0, no warning - and the CLI reference's own example used that spaced
   form. Affected every space-containing value in every CLI. The parser now
   rejects it with a copy-pasteable correction.
3. **Three silent no-ops**: a sub-cell-size `polyfit_lambda_metres` warned on
   stdout and exited 0 with no outputs; requesting no outputs at all exited 0;
   an over-large MADk2 window wrote three 100%-nodata rasters. All now fail
   with an explanation.

### Running it

```bash
ctest -L surface_metrics               # unit tests, ~1 s, no network
ctest -R benchmark_surface             # vs TopoToolbox, needs pixi + network
ctest -R performance_surface_metrics   # full sweep, ~4 min
```

The benchmark and performance tests are opt-in
(`LSDTT_ENABLE_BENCHMARKS`, `LSDTT_ENABLE_PERF_TESTS`) and skip (77) when pixi,
topotoolbox or the network is missing. Timings are reported, never asserted;
the only assertion is a correlation guard that fails when the two tools turn
out not to be computing the same kind of quantity.

---

## 11. Adding a test

1. **Can it be a property or an analytic oracle?** Prefer that over a
   comparison. `tests/cpp/test_flow_properties.cpp` is the model.
2. **Does it need data?** Add a tiny committed fixture, or add an entry to the
   fixture registry in `tests/scripts/run_flow_benchmark.py` and let
   `lsdtt-fetch-raster` fetch it.
3. **Does it need the network or an optional package?** Put it behind a CMake
   option and make the driver skip with exit code 77.
4. **Registering it**: CLI tests via `add_lsdtt_cli_smoke_test(...)` in
   `tests/CMakeLists.txt`; C++ drivers as their own executable linked against
   `lsdtt3_core` or `lsdtt3_cpu`.
5. **If you pin a tolerance**, measure first and pin just above what you
   measured — and record what you measured it against.

---

## 12. Status of this document

**This document is the durable record.** The working documents that produced
the findings above are uncommitted local scratch files: they may not be present
in any given checkout, and nothing in the tracked tree cites them. Everything
needed to understand, reproduce, or challenge a result above is in this file.
