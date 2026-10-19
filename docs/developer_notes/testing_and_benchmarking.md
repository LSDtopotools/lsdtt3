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

### Dependency tiers come first

What a test needs **to start** is labelled on its own axis, separate from what
it tests. Every test has exactly one runtime tier:

| Tier | Registered by | Needs |
|---|---|---|
| `tier_native` | `BUILD_TESTING=ON` (always) | the C++ build only |
| `tier_python_stdlib` | `LSDTT_ENABLE_PYTHON_TESTS=ON` | Python ≥ 3.10, standard library; some also GDAL's command-line programs |
| `tier_python_geo` | `LSDTT_ENABLE_PYTHON_GEO_TESTS=ON` | Python + GDAL bindings, NumPy, Rasterio, GeoPandas, Shapely |
| `tier_pixi` | `LSDTT_ENABLE_BENCHMARKS=ON` | the Pixi `benchmark` environment |

and zero or more flags: `hermetic` (native, no skip path, no network, no host
data — the reproducible floor), `needs_network` (downloads; registered only
with `LSDTT_ENABLE_NETWORK_TESTS=ON`), `needs_gdal_tools` (runs `gdalinfo`,
`ogrinfo` or `gdal_translate`). `lem_component` marks the landscape evolution
model's tests; `heavy` marks the slow ones.

A default configure registers **only the native tests**. Most CLI coverage is
in the Python tiers, so a native-only green run says the libraries are sound,
not that every CLI is. Configure prints what it registered and what each OFF
option left out, and a missing interpreter, program or Python package is a
configure error naming it.

```sh
# native only: the C++ build dependencies and nothing else
cmake -S . -B build && cmake --build build -j
ctest --test-dir build -L '^hermetic$' --output-on-failure

# add the Python tiers
cmake -S . -B build -DLSDTT_ENABLE_PYTHON_TESTS=ON -DLSDTT_ENABLE_PYTHON_GEO_TESTS=ON
ctest --test-dir build -LE '^heavy$' --output-on-failure
```

`-L` and `-LE` take regular expressions: anchor them. `-L network` matches the
`network` domain *and* `synthetic_networks`.

**Pixi and the system GDAL.** The C++ build and the Python tests can use
different GDALs. Configure the build *outside* Pixi, give CMake the Pixi
`test` environment's interpreter, and run the tests inside Pixi:

```sh
pixi install -e test
cmake -S . -B build -DLSDTT_ENABLE_PYTHON_TESTS=ON -DLSDTT_ENABLE_PYTHON_GEO_TESTS=ON \
      -DPython3_EXECUTABLE="$(pixi run -e test python -c 'import sys; print(sys.executable)')"
pixi run -e test ctest --test-dir build --output-on-failure
```

Inside Pixi, `GDAL_DATA`, `PROJ_DATA` and `GDAL_DRIVER_PATH` belong to Pixi's
GDAL. A CLI built against the system GDAL would read them without complaint —
another installation's PROJ database — so the test launcher
(`tests/scripts/lsdtt_test_env.py`) removes them, with `LD_LIBRARY_PATH` and
`CONDA_PREFIX`, from lsdtt3 programs whose GDAL is not the environment's.
`unit_test_env_isolation` checks this, with a control that shows the leak
without it. Do not activate Pixi around a system-built suite any other way.

**Execution contract.** Every test runs with `LSDTT_NONINTERACTIVE=1` and,
from Python, with `/dev/null` as stdin — CTest passes the terminal it was
started from to every test, so a prompt would hang an interactive run while
passing in CI. Every test has a timeout (`failure_*` 60 s, `heavy` 3600 s,
otherwise 300 s, unless it sets its own). OpenMP tests run with
`OMP_WAIT_POLICY=passive` so that they do not spin against each other under
`-j`; the serial, benchmark and performance tests are exempt.
`unit_test_inventory_audit` checks all of this, and every label, against the
commands the tests actually run.

**Logs and reports.** CTest's per-test log lives in
`Testing/Temporary/LastTest.log.tmp` until a run completes and is overwritten
by the next `ctest` call of any kind; it also records a timeout as "Test
Failed." and a skip as "Test Passed.". For a durable record, run through the
wrapper, and after a failing run build the report a user can send:

```sh
cmake -P build/run_tests.cmake -LE '^heavy$' -j4     # tees into build/test-logs/
cmake --build build --target test_report            # one plain-text file to send
```

### Kind labels

| Label | Covers |
|---|---|
| `smoke` | each CLI runs end to end on a tiny fixture (mostly Python-driven) |
| `rigorous` | exact output checks: checksums, sizes, feature counts, hashes |
| `failure` | bad parameters are rejected with the expected message |
| `unit` | unit tests, mostly C++ |
| `cosmo` | cosmogenic physics vs frozen CRONUS + numerics sentinel |
| `cosmo_realworld` | published-basin runs against frozen summaries (`needs_network`) |
| `surface_metrics` | surface-metrics validation (slope, curvature, SMAP, MADk2) |
| `flow` | flow routing and depression-handling properties |
| `fetch` | native DEM fetching; offline logic by default |
| `benchmark` | comparison against external references on real DEMs |
| `performance` | in-process timings, reported and never asserted |
| `configure_contract` | configures the project in a scratch directory and checks what it reports for a disabled or missing tier |

Downloading tests carry `needs_network`; `network` is the flow-routing
**domain** below, not network access.

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
ctest --test-dir build -L '^cosmo$' -LE '^heavy$'     # I changed cosmogenic code
ctest --test-dir build -L '^network$' -LE '^heavy$'   # I changed routing / chi / basins
ctest --test-dir build -LE '^heavy$'              # everything cheap
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

`heavy` is `unit_hillslope_implicit` and the eight `cosmo_realworld_*` tiers.
The implicit hillslope test **saturates every core**: 9 min 24 s wall for
73 min of CPU, a parallelism of 7.8 on 8 cores (2026-08-04).
`smoke_lem_roering_implicit_strip` was also `heavy` until 2026-09-18, when it
measured 2 s and was taken off the list.

It is marked `RUN_SERIAL`. This is not a nicety. Under `ctest -j4` a test
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
    --exe build/lsdtt-graphflood \
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
    --exe build/lsdtt-graphflood \
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

`richdem/flowmet/d8_flowdirs.hpp`, `d8_FlowDir()`, selects the neighbour with the
**lowest elevation**. Distance never enters the comparison. True D8
(O'Callaghan & Mark 1984, which both tools cite) divides the drop by the
**distance**, so a diagonal neighbour must be √2 lower to win. lsdtt3 and
TopoToolbox both implement drop-over-distance; RichDEM does not.

The divergence appears whenever the lowest neighbour is diagonal but a cardinal
neighbour is steeper. On a 1 m grid, with cardinal south 1.0 m lower and diagonal
south-east 1.2 m lower:

| | drop | distance | gradient | chosen by |
|---|---|---|---|---|
| cardinal S | 1.0 | 1.0 | **1.000** | lsdtt3, TopoToolbox |
| diagonal SE | 1.2 | 1.414 | 0.849 | **RichDEM** |

Verified by calling `d8_FlowDir` directly against RichDEM at commit `af61c65`:
it returns 6 (diagonal SE) for this configuration. No tie is involved, which is
why the disagreement is pervasive on real terrain rather than confined to flats —
measured rank correlation on a real DEM is 0.51 for D8, 0.74 for D-inf, 0.99 for
MFD (MFD spreads flow across neighbours, so the single-direction choice barely
registers).

#### Correction, 2026-08-21

An earlier version of this section made two claims about the *mechanism* that are
wrong. Both were checked against the RichDEM source and by running `d8_FlowDir`:

1. It said the tie-break "explicitly prefers diagonals". **It prefers cardinals.**
   In RichDEM's numbering `n8_diag = {false,false,true,false,true,false,true,
   false,true}`, so **odd `n` is cardinal**; the clause
   `flowdir>0 && flowdir%2==0 && n%2==1` switches *from* a diagonal *to* a
   cardinal at equal elevation. RichDEM's own source carries
   `//TODO: What is this modulus stuff for?` beside it.
2. It presented a table showing RichDEM producing `1, 2, 3, 4 …` accumulation on
   a 50×50 north–south plane. **RichDEM passes that test.** Running `d8_FlowDir`
   over the interior of exactly that ramp returns direction 7 — cardinal south —
   for **all 2304 cells**, because the tie-break above pulls it off the diagonal.
   Those RichDEM figures were never measured; the ramp test in
   `test_flow_properties.cpp` is an analytic oracle and does not invoke RichDEM.

The conclusion is unchanged — RichDEM's D8 is not steepest descent and is not a
valid D8 reference — but it follows from the *selection rule*, not from tie
behaviour, and not from the ramp.

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

The ULP approach was evaluated for lsdtt3 and is **not the default**: it collapses
to a denormal near sea level (exactly where `lsdtt-marsh-topography` operates),
varies by orders of magnitude within a single raster, is resolution-blind, and is
destroyed by any integer-typed output or by warping a filled DEM. It is now
available as an opt-in, `fill_flat_resolution_choice: barnes2014b_float_step`,
which is RichDEM's `ResolveFlats` on float32 in its first pass. Its slopes
are still ULP-based, and §10r describes how the difference from RichDEM is
confined and tested.

### 5.3 Breaching pit-flattens before carving — in both tools

`breach_depressions` is not a pure breach. Before carving, every pit/flat cell
is raised to its lowest neighbour, so "breaching never raises the surface" is
false by design — measured up to 80 m on a high-relief DEM.

RichDEM does exactly the same thing (verified: 1.0328 m rise on a shared test
fixture against lsdtt3's 1.032768 m, both raising the same 2 cells). This is a
faithful port, not drift. The invariant that *does* hold, and which the property
tests assert, is the bounded form: only pit/flat cells may rise, and never above
the elevation of their lowest original neighbour.

### 5.3b Pit filling is periodic-aware only through the boundary overload (2026-08-14)

`fill_pits(raster, min_slope, n_threads)` treats **all four edges as spill
points** and its neighbour walk does not wrap. That is right for open boundaries
and wrong for periodic ones, and the failure is not subtle: fill with open
edges, then route with `'p'`, and cells next to the seam have no downhill
neighbour, so routing aborts with

```
lsdtt_flowinfo: could not assign a receiver at row 3, col 0.
```

There is now an overload that takes the NESW boundary string:

```cpp
lsdtt_raster_cpu fill_pits(const lsdtt_raster_cpu&, float min_slope,
                           const std::string& boundary_conditions,
                           int n_threads = 0);
```

Two things change on a periodic axis: its edges are **not seeded** as spill
points, and the neighbour walk **wraps** across the seam so a depression that
spills round the edge fills to the level it actually spills at. Measured on a
40x60 noisy fixture, the periodic and open fills differ on **2182 cells, up to
37.7 m**.

Guard rails, all pinned by `unit_flow_properties` case PP:

- periodic on **both** axes is refused -- the grid would have no drainage exit
  at all and every cell would fill to the global maximum;
- mismatched pairs (`'opoo'`) are refused, matching the router's own validation;
- with no periodic flag the overload defers to the plain one, so it is
  bit-identical and keeps the reconstruction fast path;
- with a periodic flag it always takes the priority-flood backend, because the
  reconstruction backend has its own neighbour walk and has **not** been given
  the wrap.

PP3 is the check worth knowing about: it asserts that the OPEN fill still fails
under periodic routing. Without it, PP2 ("the periodic fill routes") could pass
on a fixture where nothing was ever broken.

#### The CLIs had the same latent bug

21 CLIs call `fill_pits`, and every one that also exposes `boundary_conditions`
was filling open and then routing with whatever the user asked for. Ten were
wired to pass the boundary string through:

`lsdtt-basin-extraction`, `lsdtt-basin-surface-metrics`,
`lsdtt-channel-extraction`, `lsdtt-channel-network-analysis`,
`lsdtt-chi-analysis`, `lsdtt-cross-sections`, `lsdtt-divide-migration`,
`lsdtt-drainage-analysis`, `lsdtt-lem`, `lsdtt-valley-centreline`.

`lsdtt-ridgeline-analysis` was wired later (2026-08-14), when a periodic run of
the paired-flow-path metric on the synthetic noise surfaces hit exactly the
failure above — `could not assign a receiver at row 23, col 399`. Non-periodic
behaviour is unchanged: the same DEM under `oooo` gives an identical median
ΔL and 9721 pairs against 9755 for `opop`, the difference being the seam.
Note there is **no CTest coverage of this CLI at all**, so that check was
manual; the wiring is not regression-guarded.

The change is a no-op for every non-periodic run, which is why the whole suite
is unchanged at 134/134. **Still unwired**, because their call sites did not
match the mechanical pattern and each needs a look: `lsdtt-cosmo-landslides`,
`lsdtt-cosmo-prediction`, `lsdtt-cosmo-tool`, `lsdtt-flat-surface-finder`,
`lsdtt-graphflood`, `lsdtt-sediment-routing`,
`lsdtt-stream-burning`. Those will still fail on a periodic run.
`lsdtt-basin-stats`, `lsdtt-feature-finder` and `lsdtt-mountain-front` do not
expose `boundary_conditions` at all, so they cannot route periodically anyway.

`breach_depressions` was **not** made periodic-aware; its seeding helper is a
separate path.

---

### 5.4 lsdtt3's D8 router requires a gradient on flats

Feeding lsdtt3 a `fill_min_slope: 0` DEM produces:

```
lsdtt_flowinfo: could not assign a receiver at row 12, col 129.
The DEM likely contains an unresolved sink or flat.
```

This is by design — the imposed gradient is what makes flats routable — but it
means the bit-exact `min_slope = 0` fill is not a usable input to lsdtt3's own
flow routing. `fill_min_slope` is load-bearing, not cosmetic.

#### That message used to kill the process (fixed 2026-08-14)

Until 2026-08-14 the `throw` above sat inside the `#pragma omp parallel for` that
builds the receiver graph. **An exception that leaves an OpenMP structured block
calls `std::terminate`**, so the error could not be caught by anyone: a correctly
written `try`/`catch` around the routing call never ran, and the process died
with a bare libstdc++ message. Under several threads it was worse — multiple
threads threw at once and the output interleaved into garbage:

```
terminate called recursively
terminate called after throwing an instance of '  what():  lsdtt_flowinfo: could not assign a receiver at row 0, col 2. ...
std::runtime_error'
```

This affected **every CLI**, not one tool. Any DEM with an unresolved sink or
flat — a diverged hillslope solve, a closed outlet boundary, a bad
`fill_min_slope` — aborted rather than reporting.

The fix records the failing cell inside the region and throws once the region has
closed. Two properties worth knowing:

- **It is catchable now.** The same failure that used to abort now surfaces as a
  normal exception, so callers can add context. `lsdtt-synthetic-networks`
  reports `FATAL: landscape equilibration failed: lsdtt_flowinfo: could not
  assign a receiver at row 0, col 2 ...` and exits 1.
- **It is deterministic.** The cell reported is the lexicographically first
  (row, then col) failing cell, not whichever thread arrived first. Verified
  identical at 1, 2, 4 and 8 threads; before the fix the same run reported rows
  47, 25, 37 and 16 on different threads of one execution.

A tree-wide scan found this was the **only** `throw` reachable from inside an
OpenMP region. `lsdtt_cosmo_raster_cpu.cpp` has parallel regions that call
`reset_scaling()`, which can throw, but every call site is guarded by
`if (!is_cronus_emulator(...))` so it never does for the schemes that reach it —
a latent risk if that guard is ever weakened, not a live defect.

> Writing the scan is easy to get wrong. The first version reported the tree
> clean, including the file with the known defect still in it, because it closed
> each region at the first line that did not increase brace depth — and the
> `#pragma` sits inside `#ifdef _OPENMP`, so the `#endif` closed it immediately.
> **Validate a checker against a known positive before trusting a clean result.**
> Matching the bare word `throw` also hits comments; every remaining hit in the
> corrected scan was a comment containing "throws".

---

## 5b. The standalone benchmark suite (cross-implementation)

Everything else in this document is testing: it runs under CTest, it asserts,
and it can fail a build. **The benchmark suite is not.** It asserts nothing,
CTest never invokes it, and its only product is a measurement of how lsdtt3
compares against other implementations on the machine it is run on. It exists
so that a user, or a reviewer, can reproduce a published comparison without
knowing anything about our test wiring.

```sh
pixi run bench-laptop            # reduced fixtures, minutes
pixi run bench-server            # the peer-review run
pixi run bench-laptop --dry-run  # print the plan, measure nothing
```

Driver: `tests/scripts/run_benchmark_suite.py`. Read its module docstring
before quoting any number from it -- it states what each measurement means and
where the comparison is still not apples to apples.

### The two profiles

| Profile | Fixtures | Threads | Repeats | Purpose |
|---|---|---|---|---|
| `laptop` | `*_small` | 1, 2, 4 | 3 | a real baseline on a laptop, **and** the smoke test before dispatching a server campaign |
| `server` | small + medium + large | 1, 2, 4, 8 | 5 | the profile whose output is quotable |

Run `bench-laptop` before `bench-server`. It exercises every code path in a few
minutes, so a broken executable path or a missing fixture surfaces immediately
rather than after a queue wait.

### Slurm is deliberately outside the benchmark

`tests/scripts/slurm_benchmark.sbatch` contains **no measurement logic**. It
arranges for a scheduler node to run `pixi run bench-server`, the identical
command a laptop runs, and nothing in the suite knows Slurm exists. Sites with
PBS, LSF, SGE or no scheduler at all run the command directly and get an
identical campaign directory. Keep it that way: put logic in the sbatch file
and Slurm sites stop being comparable with everyone else.

The one thing it does enforce is `--exclusive`. A benchmark sharing a node is
measuring the co-tenant's memory bandwidth as much as its own, and it refuses
to start if `--cpus-per-task` is below the largest thread count in the sweep --
an "8 thread" row on a 4-core allocation is not a slow result, it is a wrong
one, and it is indistinguishable from a genuine scaling ceiling.

### Where things are written, and how to run from a data directory

Three locations, with three different lifetimes:

| what | default | override | notes |
|---|---|---|---|
| campaign output | `./benchmark_campaigns/` | `--out-dir` | small; the thing you keep and send |
| fixture cache | `<repo>/tests/data/dem_cache/` | `--fixture-dir` | ~512 MB for `server`; reused by every later campaign |
| staged fixtures + intermediates | `$TMPDIR` | `--scratch-dir` | node-local, disposable, must be fast — measurements happen here |

A relative `--out-dir` resolves against the **current directory**, so running
from a data directory puts the campaign in that data directory.

**The one trap:** `pixi run` sets the working directory to the manifest root,
so `pixi run bench-server --out-dir results` lands in the repo regardless of
where you invoked it. To collect results elsewhere, pass an **absolute** path,
or call the script directly:

```sh
cd /path/to/my/data
pixi run --manifest-path /path/to/lsdtt3/pixi.toml bench-server \
    --out-dir "$PWD" --fixture-dir /scratch/$USER/lsdtt3_fixtures
```

**On a cluster, move the fixture cache off the repo.** It defaults inside the
repository, which usually sits on a quota'd home filesystem, and the `server`
profile stores ~512 MB of region sources there. Point `--fixture-dir` at
project or scratch storage — somewhere persistent and shared, because it is
reused by every later campaign and re-downloading it is the slowest part of a
cold run.

### Running under Slurm without editing anything

`slurm_benchmark.sbatch` needs **no edits**. Everything site-specific is an
sbatch flag or an environment variable, so the file can stay exactly as
committed.

**Slurm's own settings** go to `sbatch`, which overrides the `#SBATCH` lines,
or through Slurm's native `SBATCH_*` variables. There is deliberately **no
`#SBATCH --account` line** — a placeholder would silently become the default
the one time you forgot to override it.

```sh
sbatch --account=PROJECT --partition=compute --cpus-per-task=32 \
       tests/scripts/slurm_benchmark.sbatch
export SBATCH_ACCOUNT=PROJECT SBATCH_PARTITION=compute   # same effect
```

**The script's settings** are environment variables:

| variable | purpose |
|---|---|
| `BENCH_SETUP_CMD` | shell eval'd first — where `module load` / `conda activate` go |
| `BENCH_SETUP_SCRIPT` | a file to source instead (both run if both are set) |
| `BENCH_RUNNER` | `pixi` (default) or `python` |
| `BENCH_PYTHON` | interpreter for `BENCH_RUNNER=python` |
| `BENCH_BUILD_DIR` | where the `lsdtt-*` binaries are; default `<repo>/build` |
| `BENCH_REPO` | the checkout; default derived from the script's own path |
| `BENCH_OUT_DIR` | campaign output; default `<submit dir>/benchmark_campaigns` |
| `BENCH_FIXTURE_DIR` | fixture cache; default `<submit dir>/benchmark_fixtures` |
| `BENCH_PROFILE` | `server` (default) or `laptop` for a shakedown |

Arguments after the script name reach the suite unchanged. A full conda
example, with results in the directory you submit from and fixtures on scratch:

```sh
cd /path/to/my/data
export BENCH_SETUP_CMD="module load anaconda; conda activate lsdtt3"
export BENCH_RUNNER=python
export BENCH_BUILD_DIR=/path/to/lsdtt3/build
export BENCH_FIXTURE_DIR=/scratch/$USER/lsdtt3_fixtures
sbatch --account=PROJECT --cpus-per-task=32 \
       /path/to/lsdtt3/tests/scripts/slurm_benchmark.sbatch --threads 1,2,4,8,16,32
```

**`conda activate`, not a bare interpreter path.** Pointing `BENCH_PYTHON` at
an environment's `python` binary without activating gives you GDAL with no
`PROJ_DATA`, and GDAL then falls back to CRS parameters from the GeoTIFF keys
rather than the EPSG registry. The script checks and says so before starting.

**A missing reference is silent.** The suite skips any reference it cannot
import (exit 77) rather than failing, so an environment without RichDEM,
TopoToolbox or Whitebox produces a campaign that completes "successfully" with
lsdtt3-only rows. This is not hypothetical: an early version of the sbatch ran
`pixi run` without `-e benchmark`, landed in the default environment, and
skipped all three. Both runners now print what is importable before the
campaign starts — **read that list before trusting the tables.**

### What each campaign writes

One self-contained, tar-and-send directory per run, under
`benchmark_campaigns/<profile>_<host>_<timestamp>/`:

| File | Contents |
|---|---|
| `benchmark_raw.csv` | one row per measurement, every column below |
| `SUMMARY.md` | rendered tables, medians, scaling and memory |
| `hardware_manifest.json` | CPU, cores, RAM, governors, compiler, git SHA, reference versions |
| `run_manifest.json` | profile settings, fixture hashes, measured overheads |

### Why the numbers are defensible

- **Threads are pinned for every engine.** Whitebox's `max_procs` defaults to
  `-1` (all cores) and lsdtt3 takes every OpenMP thread it can; an unpinned
  comparison measures core count. `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS` and
  `MKL_NUM_THREADS` are all set, the last two because a transitively-loaded
  BLAS will otherwise spin up its own pool and the "1 thread" row is not one
  thread.
- **Reference subprocesses are pinned too**, via `reference_env()`, which sets
  the same variables without stripping the pixi paths the Python references
  need to import at all. `system_env()` (for the C++ CLIs) and
  `reference_env()` (for the Python references) need opposite treatment of
  `LD_LIBRARY_PATH` and must stay separate.
- **RichDEM and TopoToolbox are serial** for these routines. They are measured
  once, at the **baseline** (lowest) thread count in the sweep, and flagged
  `serial_reference`. A flat scaling curve for them is not a finding, and
  re-running them at 8 threads would record the same number under a label
  implying it had been tested there. (Baseline, not a hard-coded 1: with
  `--threads 2,4` an earlier version produced no reference rows at all and they
  vanished from the campaign silently.)
- **Every measurement has a deadline**, `--timeout-seconds` (default 1800). A
  killed measurement becomes a `timeout` row -- a distinct status from `error`,
  because the engine was working, just too slow -- and the campaign carries on.
  This is not hypothetical: Whitebox's breach on low-relief terrain measured
  ~300x slower than lsdtt3's, and one such cell would otherwise stall an
  unattended server campaign for hours.
- **Oversubscription is refused, not measured.** The runner warns when the
  thread sweep exceeds the machine's logical cores, and the Slurm wrapper
  refuses to start when `--cpus-per-task` is below the sweep. An "8 thread" row
  on 4 cores is a mislabelled result indistinguishable from a real scaling
  ceiling.
- **The exact invocation is recorded** in `run_manifest.json` (`argv`, cwd,
  interpreter). Every profile setting is overridable from the command line, so
  the profile name alone does not say what was actually run.
- **Repeats, not single shots**: a discarded warm-up then N timed repeats, with
  best, median, stdev and the full list all in the CSV. `SUMMARY.md` quotes the
  **median** -- best-of-N answers "how fast can this go" and a user wants "what
  will I get".
- **I/O is separated from compute.** The engines have very different I/O layers
  and a bundled number compares file formats.
- **Memory** is child-process peak RSS via `wait4(2)`, so it sees the C++ and
  Rust allocation a Python-level probe cannot.
- **Fixture identity is hashed** into every row. Two campaigns on differently
  built fixtures are comparing different rasters.
- **A non-Release build is called out** in `SUMMARY.md`, as is a CPU governor
  that is not uniformly `performance`.

### Timing scope: the thing to check before quoting a ratio

Rows are only ever compared with others measured the same way, recorded in the
`timing_scope` column:

| Scope | Meaning |
|---|---|
| `cli_reported_compute` | the compute time the lsdtt3 CLI printed itself; excludes startup and I/O |
| `in_process_routine` | a reference script timing the routine in-process |
| `whole_process` | wall clock for an entire lsdtt3 CLI run |

`lsdtt-surface-metrics`, `lsdtt-raster-preprocessing` and
`lsdtt-drainage-analysis` all print `Compute time:` (and `Write time:`) lines,
which the driver parses and sums, so the flow, surface and depression suites are
like-for-like. The drainage CLI prints the routing build and the accumulation
separately, matching how TopoToolbox splits the same work; the driver sums them.

That printing was added *for* this suite, and it matters more than it sounds:
lsdtt3's fill on `alpine_small` is **7.6 ms of compute inside a 279 ms
process**. Whole-process wall clock was 97% startup and I/O. Any CLI that
prints no timing falls back to `whole_process`, and those rows are tabled
separately in `SUMMARY.md` and never enter a ratio.

### Fixed per-call cost is measured for every engine

Do **not** assume an engine has no fixed overhead because it is "timed
in-process". That was assumed here once and it was wrong. Each campaign
measures every engine's per-call cost on an 8×8 crop (64 cells), where the
algorithm itself is microseconds. Median of 7, `fill`, on crops of
`alpine_small`:

| grid | richdem | topotoolbox | whitebox (1t) |
|---|---|---|---|
| 8² | **9.83 ms** | 0.036 ms | 0.038 ms |
| 16² | **11.56 ms** | 0.029 ms | 0.032 ms |
| 64² | **10.19 ms** | 0.111 ms | 0.126 ms |
| 256² | **13.73 ms** | 1.966 ms | 1.440 ms |

TopoToolbox and Whitebox extrapolate to ~0.03 ms at zero cells. **RichDEM
extrapolates to ~10 ms.** Its own internal timer prints
`Zhou2016 wall-time = 1.2e-05 s` for the 8×8 fill, so ~99.8% of what a caller
measures there is RichDEM's Python binding: a per-call `rdarray` construction
plus the citation block it formats and prints on *every* invocation. On
`alpine_small` that is roughly half of RichDEM's reported fill time.

`overhead_dominated` therefore compares each row against **its own engine's**
floor. Holding lsdtt3 alone to a startup floor would have quietly let RichDEM
off, and RichDEM is the reference lsdtt3 is most often compared against.

### Compute-time printing is a parsed interface, not a log message

`lsdtt-surface-metrics`, `lsdtt-raster-preprocessing` and
`lsdtt-drainage-analysis` print their own compute times, and **two scripts parse
those lines**: `run_surface_benchmark.py` (`FIT_TIMING_RE`) and
`run_benchmark_suite.py` (`parse_reported_times`). Two phrasings exist:

```
  Compute time: 0.0026 s                      the common form
  Write time: 0.0090 s
  Computed hillshade... done (0.0035 s)       surface-metrics only
  Running direct fit... done (0.0041 s)
```

The `... done (X s)` phrasing looks like something to tidy up. **It is not.**
Normalising hillshade to `Compute time:` while building this suite broke
`run_surface_benchmark.py`, which raises rather than silently falling back to
wall time when its pattern stops matching -- a good guard that caught it
immediately. The CLI wording was restored and the new parser taught to accept
both forms. When a CLI's output is parsed, the parser accommodates the CLI, not
the other way round.

### Fixtures: one download per region, then crop the ladder

The suite fetches **one large real DEM per region** with `lsdtt-fetch-raster`
(COP30) and crops every rung of the size ladder out of it.

| profile | regions | ladder (px) | source fetched at | download |
|---|---|---|---|---|
| `laptop` | alpine, lakeland | 500, 1000 | 1000 px (30 km) | ~4 MB/region |
| `server` | alpine, lakeland | 1000, 2000, 4000 | 4000 px (120 km) | ~64 MB/region |

The `server` ladder went to 8000 px until 2026-09-13. On the forth02 campaign
of 2026-09-09 that rung took 4.2 h of 5.5 h. For the surface suite, D8 and
D-infinity, every engine's 8000² throughput was 0.95–1.04× its 4000² figure.
The fill family still fell further (0.59–0.82×), but that fall is confounded
with terrain, because each rung is a wider crop. A cached larger source is
reused rather than re-fetched, so rung hashes do not change.

| region | centre | UTM | terrain |
|---|---|---|---|
| `alpine` | 86.9 E, 27.9 N | 32645 | High Himalaya, Khumbu/Makalu |
| `lakeland` | 92.0864 W, 50.5734 N | 32615 | Canadian Shield, NW Ontario |

**Why crop rather than fetch each rung.** A ladder built from independent
downloads confounds size with terrain: at 30 m, 2400 px is a 72 km box and
8000 px is 240 km, so scaling the box up means "alpine" stops being alpine.
Cropping every rung from one source keeps terrain as constant as it can be and
needs one download per region instead of one per rung.

**Why not synthetic.** `lsdtt-synthetic-surfaces` generates a
diamond-square surface of any size in ~1.3 s with no network, which is
tempting, and it is not usable here. Depression density (fraction of cells a
no-epsilon RichDEM fill raises, 256²):

| surface | relief | depression cells |
|---|---|---|
| REAL alpine (Italian Alps) | 1976 m | **0.94 %** |
| REAL lakeland (Finland) | 51 m | **1.29 %** |
| SYNTH diamond-square, no noise | 4.3 m | 11.54 % |
| SYNTH roughness 0.70 + noise | 6.8 m | 34.24 % |
| SYNTH roughness 0.80 + noise 10 | 15.4 m | 37.44 % |

10–30× more depressions than real terrain, at 1 m cell size with a few metres
of relief. Since depression structure is the largest single driver of the fill
and breach timings, a fractal surface would inflate every depression-handling
measurement for every engine. Synthetic surfaces belong to the analytic oracles
in the CTest suite, not here.

**The confound that remains, and how it is exposed.** Centred crops of
increasing size necessarily reach into different landscape — a 240 km box in
the Himalaya includes both the Gangetic plain and the Tibetan plateau. This
cannot be engineered away with real terrain, so it is *measured*:
`describe_fixture()` records relief and true depression density (RichDEM
priority-flood) for every rung into `run_manifest.json`. A scaling curve across
rungs whose depression density doubles is not a clean scaling curve, and the
reader can see that rather than having to assume otherwise.

Do not substitute a cheap proxy for that number. An earlier version used a 5×5
`grey_closing` and reported 29–55% where the true figures were 6.4% and 18.5% —
it flags any local concavity, not a closed basin. A statistic wrong by 3× is
worse than no statistic, because it looks authoritative in the manifest.

**These regions are much harder than the ones they replace**, which is the
point:

| region | relief | depression cells |
|---|---|---|
| Himalaya (new alpine) | 1662 m | **6.40 %** |
| Ontario (new lakeland) | 94 m | **18.51 %** |
| Italian Alps (old) | 1976 m | 0.94 % |
| Finland (old) | 51 m | 1.29 % |

`--regions` and `--ladder` override the profile. `--fixtures` bypasses the
system entirely and takes raster names from the cache verbatim — the escape
hatch for legacy or hand-placed fixtures.

### Fixtures are staged to local disk, and this is not optional

**The single largest measurement error found while building this suite.**
Whitebox touches its working directory inside a tool call, so the filesystem
the fixture sits on lands inside the timed region. Same 256² raster, same
routine, `fill`:

| filesystem | whitebox | richdem | topotoolbox |
|---|---|---|---|
| `9p` (WSL2 `/mnt/c`) | **16.4 ms** | 11.0 ms | 1.96 ms |
| local `ext4` | **2.6 ms** | 12.5 ms | 2.15 ms |

A **6x penalty on Whitebox alone**, entirely attributable to the filesystem.
RichDEM and TopoToolbox are unaffected. A campaign run straight off `/mnt/c`
reported that lsdtt3 beat Whitebox on fill and breach; after staging, Whitebox
is several times *faster* on both. The conclusion had inverted.

The suite now copies every fixture, and the shared filled DEM, to local storage
before measuring anything. `--scratch-dir` chooses where (default: a
subdirectory of the system temp directory, which is node-local on a cluster and
therefore the right place anyway). The campaign records the filesystem it
measured on, and `SUMMARY.md` refuses to endorse the run if that filesystem is
one of `SLOW_FILESYSTEMS` (9p, drvfs, cifs, nfs, fuseblk, overlayfs).
`--no-stage` exists only to demonstrate what staging is worth.

**If you take one thing from this section:** check the "measured on" line at
the top of any `SUMMARY.md` before quoting a number from it.

### Two measurement artefacts found while building this

Both are recorded because they are the kind of thing that silently produces a
wrong published table.

**Whitebox tool-handle resolution cost ~190 ms per attribute access.**
`wbe.hydrology.sink` is not an attribute read; it is a registry resolution
through the Rust layer. 1000 accesses took 191.4 s, while 1000 accesses of the
namespace `wbe.hydrology` took 0.8 ms. Resolving the handle inside the timing
lambda charged every Whitebox measurement a fixed ~190 ms:

| fill on `alpine_small` (256²) | seconds |
|---|---|
| handle resolved inside the timing loop | 0.207 |
| handle bound once, outside it | 0.019 |

That is a **10.9x artefact, and it flattered lsdtt3.**
`whitebox_reference.py` now binds every handle before timing. Each campaign
re-measures the cost as `whitebox_dispatch_seconds` and `SUMMARY.md` prints a
warning refusing to endorse the Whitebox rows if it exceeds 100 ms, so the
regression cannot come back quietly.

**Small fixtures are dominated by fixed cost.** Whitebox's fill is 0.29 µs/cell
on 256² and 0.09 µs/cell on 2400²; the difference is overhead amortising, not
an algorithm. Every row is checked against its own engine's measured floor and
flagged `overhead_dominated` if it is within 3x of it. The `laptop` profile
proves the pipeline works; the `server` profile, with the large fixtures,
produces the numbers.

### Two like-for-like failures found on the first server campaign (2026-08-31)

The first `laptop`-profile campaign on forth03 reported lsdtt3 **slower than
TopoToolbox on every single-core comparison** — 20 of 20, with fill 2.9× behind.
Investigating against the TopoToolbox source showed **both headline losses were
the benchmark comparing unlike things**, not lsdtt3 being slow. Recorded here
because the corrected cases are only meaningful alongside what they replaced.

#### 1. `fill` compared two different algorithms

lsdtt3 implements both algorithm families; TopoToolbox implements only
reconstruction, and the two default to a *different* family:

| | hybrid reconstruction (O(N), FIFO queue) | sequential reconstruction (repeated scans) | priority-flood (O(N log N), heap) |
|---|---|---|---|
| lsdtt3 | `vincent1993` | — | `barnes2014` ← **CLI default** |
| TopoToolbox | `hybrid=True` ← **default** | `hybrid=False` | — |

**TopoToolbox has no priority-flood fill.** `fillsinks(hybrid=False)` calls
libtopotoolbox `reconstruct`, Vincent's (1993) *sequential* algorithm (forward
and backward raster scans repeated until nothing changes); `hybrid=True` calls
`reconstruct_hybrid`. Neither calls the `pq_*` heap routines (call graph of the
installed `topotoolbox/_grid*.so`). The sequential pass count is data-dependent:
33.5 s on lakeland_4000 against 3.2 s on alpine_4000.

The suite used both defaults, so it ran lsdtt3's heap against TopoToolbox's
reconstruction. On alpine_1000 at `fill_min_slope: 0`, output **bit-identical**
(0 differing cells of 1,000,000):

| | compute | throughput |
|---|---|---|
| lsdtt3 `barnes2014` | 0.1282 s | 7.80 Mcells/s |
| lsdtt3 `vincent1993` | 0.0446 s | **22.42 Mcells/s** |
| TopoToolbox `hybrid=True` | — | 22.84 Mcells/s |

**Matched, lsdtt3 is within 2% of TopoToolbox.** The reported 2.9× gap was
entirely the algorithm choice.

**Fix:** the `fill` case now pins `fill_algorithm_choice: vincent1993` against
`hybrid=True`, and a `fill_priority_flood` case pins `barnes2014` against
RichDEM (Zhou 2016) and Whitebox. Both families are measured; neither is
compared across.

Until 2026-10-08 `fill_priority_flood` also carried TopoToolbox
`hybrid=False`, on the mistaken belief that it was a heap. It was sequential
reconstruction against lsdtt3's priority-flood, so those rows are not
like-for-like. Their timings remain valid measurements of sequential
reconstruction. The routine is removed: lsdtt3 has no sequential-reconstruction
backend to pair it with.

`fill_epsilon` has no reconstruction row, and cannot: reconstruction fills to
exactly the spill level and has no mechanism for a gradient across flats.

#### 2. `accum_d8` charges lsdtt3 for the S-vector

lsdtt3's D8 build produces an S-vector the references do not build. Measured on
a 1000² filled tile, single-threaded:

| configuration | build |
|---|---|
| all on (what the CLIs request) | 186.4 ms |
| without the acyclicity check | 190.0 ms (free) |
| without the topological order | 178.3 ms |
| **without the S-vector** | **80.1 ms** |
| receivers only | 52.4 ms |

The S-vector is **~106 ms of the 186 ms**. `get_drainage_area_raster()` then
just reads its precomputed subtree sizes, which is why the accumulation phase
measures 1.4 ms — the work was already done. Accumulating over the topological
order gives the identical raster (0 differing cells of 1,000,000) for
0.1341 s against 0.1951 s, **1.45× faster**.

**This is deliberately NOT changed**, and the CLI keeps building the S-vector.
It is what makes `is_upslope_of` an O(1) interval test, and chi analysis, basin
summaries, source finding and the cosmogenic routines all depend on it. See
`docs/design_and_algorithms/flow_routing.md`.

**Fix:** the `accum_d8` case note now states that lsdtt3 performs work the
references do not, and quantifies it. The asymmetry is recorded rather than
engineered away — but do not quote that row as a like-for-like accumulation
comparison.

### Surface metrics: what the shared fit costs and when it pays

lsdtt3 fits one 6-coefficient quadratic and derives every metric from it;
TopoToolbox and Whitebox recompute their derivatives per metric. The suite runs
`polyfit_method: direct` throughout. Measured on alpine_1000:

| phase | time |
|---|---|
| shared `direct` fit | **82.5 ms** |
| + slope | 28.0 ms |
| + aspect_degrees | 44.3 ms |
| + profile_curvature | 22.3 ms |
| + plan_curvature | 20.4 ms |

The fit is ~40–55% of the total at four metrics and is paid **once**, so
lsdtt3's position improves with every metric added. End to end on one machine:

| metrics | lsdtt3 | if each were computed independently |
|---|---|---|
| 1 (`slope_only`) | 93.7 ms | 93.7 ms |
| 4 (`polyfit_all`) | ~199 ms | ~375 ms |
| 8 (`polyfit_many`) | 275.5 ms | ~749 ms |

lsdtt3 exposes ~25 metrics from the same fit.

#### Marginal cost, and where lsdtt3 becomes cheaper

Measured as a series rather than extrapolated from two points, **curvatures
only** on alpine_1000, single-threaded, best of 3–4. Matching the metric family
matters — mixing slope and aspect in distorts both sides (see the correction
below).

| | fixed | per metric |
|---|---|---|
| **lsdtt3** | 50.6 ms (the quadratic fit) | **24.3 ms** |
| **TopoToolbox** | 8.0 ms | **50.2 ms** (5 convolutions each) |

TopoToolbox's marginal cost is **2.1×** lsdtt3's, and it carries almost no fixed
cost because it has no shared fit to amortise.

| N curvatures | lsdtt3 | TopoToolbox | winner |
|---|---|---|---|
| 1 | 79.5 ms | 62.4 ms | TopoToolbox |
| **2** | **98.6 ms** | **102.7 ms** | **lsdtt3** |
| 3 | 117.2 ms | 157.0 ms | lsdtt3 |
| 4 | 151.6 ms | 211.6 ms | lsdtt3 |
| 5 | 162.2 ms | [258.8 ms] | lsdtt3 |
| 6 | 204.5 ms | [309.0 ms] | lsdtt3 |

`[]` extrapolated; TopoToolbox exposes only four curvature types. Crossover is
at **N = 1.65**, so **lsdtt3 is cheaper from two metrics onward** — and N = 2, 3
and 4 are directly measured rather than fitted. By four metrics lsdtt3 is 1.4×
ahead, and the gap widens linearly.

!!! warning "Correction: an earlier estimate of ~6 metrics was wrong"
    This document previously put the crossover at "roughly six metrics onward",
    from a two-point extrapolation over a MIXED metric set (slope + aspect + 2
    curvatures). Mixing families distorts both sides: TopoToolbox's
    `evansslope` is 2 convolutions rather than 5 — and much of its 68 ms is a
    `distance_transform_edt` NaN infill that lsdtt3 does not perform — while
    `aspect_degrees` is lsdtt3's *dearest* derivation at 44 ms, roughly double
    its other metrics. Compare curvature against curvature.

**The crossover depends on which metrics you want.** Two is the figure for
curvatures. An aspect-heavy set moves it right, because aspect is lsdtt3's most
expensive per-metric operation. Quote the metric family alongside the number.

**Three reasons the surface rows are `report_only`:**

1. **`slope_only` is lsdtt3's worst case by construction.** It pays the full
   quadratic fit to produce one number, where TopoToolbox runs 2 convolutions.
   Reported because it is the honest counterweight to `polyfit_all`; quoting
   only the favourable case would be cherry-picking.
2. **The metric sets differ.** `polyfit_all` yields slope + aspect + 2
   curvatures; TopoToolbox's `curvature_all4` yields 4 curvatures. Not the same
   products.
3. **`polyfit_many` has no counterpart at all** — no reference exposes an
   eight-metric call. It was an lsdtt3-only row, included to show the fit
   amortising. It was removed from the campaign on 2026-09-13: it cost ~44 min
   of a server campaign, mostly writing eight rasters. The table above remains
   the record of the amortisation.

!!! warning "Correction: window size is free, contrary to an earlier claim"
    This document previously stated that lsdtt3's default
    `polyfit_lambda_metres: 70` is a wider window than the references' fixed
    3×3, so lsdtt3 does more work per cell and its surface numbers are
    "conservative". **That is wrong.** The `direct` method is a separable
    rolling-moment fit whose cost is independent of window size — measured flat
    from lambda 35 m to 280 m (0.0191 s to 0.0193 s on 1000², a 729× range of
    window cells). Lambda does not tilt the comparison in either direction.

### Historical: why the original parallel scaling was poor (2026-08-31)

!!! warning "Historical implementation"
    This section describes the original two-pass `direct_polyfit`. On
    2026-09-01 the prefix storage was bounded and the passes were fused over row
    strips. The full-raster `H0/H1/H2/Hn` arrays and full-height prefix buffers
    discussed below no longer exist. Keep these measurements as provenance, but
    do not use them to diagnose the current fused loop. The current diagnosis
    and next measurements follow this historical section.

The first server campaign showed **every lsdtt3 scaling row at or below
1.00×** — adding threads made lsdtt3 slower on every operation, on a 192-core
machine. Worst cases were `slope_only` and `polyfit_all` at 0.64–0.67× on two
threads. This section records what that turned out to be, because the shape is
easy to misread as a threading bug and it is not one.

#### First establish what the machine can do

Do not interpret a scaling curve without a ceiling to compare it against. Two
OpenMP kernels at the same working-set size (1M doubles), on the development
machine — an i5-1145G7, **4 physical cores, 8 threads**:

| kernel | 1t | 2t | 4t | 8t |
|---|---|---|---|---|
| compute-bound (60 `sqrt` per element) | 1.00× | **1.83×** | **3.09×** | 2.96× |
| memory-bound (streaming triad, ~24 MB) | 1.00× | **1.03×** | **1.14×** | **0.27×** |

The machine parallelises compute fine. A bandwidth-bound kernel gets nothing —
and at 8 threads it is **4× WORSE than serial**, because hyperthreads sharing a
core's load/store ports simply contend.

`direct_polyfit` measures 1.17× at 2t and degrades after. That places it
squarely in the memory-bound band. **The OpenMP is correct**: the pragmas,
`num_threads` clauses and `resolve_thread_count` all do what they should. The
kernel is not compute-limited, so there is nothing for extra threads to win.

#### Which pass, and why

Instrumenting the two passes of `direct_polyfit`, alpine_1000, best of 3:

| pass | 1t | 2t | 4t | speedup |
|---|---|---|---|---|
| 1 — horizontal rolling sums (parallel over rows) | 0.0154 s | 0.0112 s | 0.0122 s | 1.37× → 1.26× |
| **2 — vertical prefix + solve (parallel over column blocks)** | **0.0264 s** | **0.0269 s** | **0.0286 s** | **0.98× → 0.92×** |

**Pass 2 is 63% of the fit and has negative speedup.** Pass 1 at least gains at
two threads.

The cause is pass 2's per-thread working set. Each 8-column block allocates
seven full-height prefix arrays:

```cpp
const std::size_t np = (nrows + 1) * BLK;
std::vector<double> ph0(np), ph0r(np), ph0r2(np), ph1(np), ph1r(np), ph2(np), phn(np);
```

7 × (nrows+1) × BLK doubles **per thread**, growing with the raster:

| nrows | per-thread buffers |
|---|---|
| 1000 | 438 KB |
| 2000 | 875 KB |
| 4000 | 1.7 MB |
| **8000** | **3.4 MB** |

against 1280 KB of L2 per core and 8 MB of shared L3. The 438 KB buffer at
1000² fits L2; it exceeds 1280 KB only above about 2900 rows. At the server
profile's 8000² rung, N threads need 3.4N MB and
thrash L3 collectively. **More threads make it worse, which is the observed
shape.** It also means the scaling gets worse as the DEM gets bigger — the
opposite of the usual expectation, and a trap when reading the ladder.

The arrays are only ever consumed as a difference of two rows:

```cpp
const double raw_S_z = ph0[pb + j] - ph0[pa + j];
```

So a full-height prefix array was materialised to compute a sliding vertical
window. A naive rolling sum would use **O(BLK) scalars**, but it changes
floating-point addition order and is not bit-identical. The landed replacement
instead retains the prefix recurrence in bounded strip storage and carries its
exact accumulator between strips.

#### Historical NUMA hypothesis for the server result (0.64× vs 1.17×)

Untested — forth03 was not available for instrumented runs — but the likely
additional factor is **NUMA first-touch**. `H0`, `H1`, `H2` and `Hn` (~28 MB at
1000²) are allocated and zero-filled by the master thread *before* the parallel
region:

```cpp
std::vector<double> H0(npix, 0.0);
```

Linux binds those pages to the master thread's NUMA node. On a 192-core box —
certainly multi-socket — threads on other sockets then read across the
interconnect on every access. For a bandwidth-bound kernel that turns added
threads into a net loss, which fits 0.64× better than a single-socket 1.17×.

Testable on the server with no code change:

```sh
numactl --hardware                                   # confirm multi-socket
numactl --cpunodebind=0 --membind=0 <the 2-thread run>
```

If pinning to one node recovered the scaling, it would identify NUMA as a
factor. The full-raster `H` arrays have since been removed. The six output
coefficient rasters remain possible first-touch candidates, but only after an
exclusive server run confirms that NUMA still matters.

#### A hypothesis that was tested and REJECTED

Pass 2 constructs 875 `std::vector<double>` per call (7 per block × 125
blocks), ~56 MB of allocation and zeroing, all inside the parallel loop body.
Allocator contention across threads looked like an obvious culprit.

**It is not.** Hoisting the buffers to thread-private — allocated once per
thread instead of once per block — produced **identical timings**:

| | 1t | 2t | 4t | 8t |
|---|---|---|---|---|
| original (per-block allocation) | 0.0588 s | 0.0504 s | 0.0530 s | 0.0600 s |
| hoisted (per-thread allocation) | 0.0587 s | 0.0523 s | 0.0534 s | 0.0594 s |

Recorded so the experiment is not repeated. glibc's per-thread arenas already
make these allocations cheap; the working-set size is what matters, not the
allocation count.

#### Current fused-loop result and leading candidates (2026-09-01)

!!! warning "Superseded on 2026-09-21"
    The server measurements in this section were taken on a session restricted
    to two CPUs by a login-node affinity mask, which was not known at the time.
    The fused loop **does** scale: 3.22× at 4000² on four threads, measured in
    a real allocation. See
    [Server campaign 2026-09-21](#server-campaign-2026-09-21-forth01-release-052-thread-scaling-measured).
    The laptop measurements below are unaffected — that machine was never
    restricted — and the investigation order at the end of this section is
    still the right order for the residual gap between 3.2× and 4×.

The bounded-prefix and fused-strip changes were bit-identical on all 324
reference rasters. They made `direct_polyfit` 20-73% faster and reduced peak
memory by 13-25%, but did not produce useful thread scaling: one and four
threads remain within a few percent. The old pass profile above cannot explain
the current result because those passes and working sets changed.

The leading current candidate was the shared, full-width `H` strip. Its size is
`(RSTRIP + 2W + 1) * ncols * 28` bytes. With the original fused-loop
`RSTRIP=256` and `W=1` it occupied
7.25 MB at 1000 columns and 29.0 MB at 4000 columns, against 8 MB shared L3 on
the development machine. Pass 1 writes the complete strip before pass 2 reads
it, so wide rasters can still evict `H` and recreate substantial DRAM traffic.
Pass 1 distributes rows while pass 2 distributes column blocks, so the thread
that produced a cache line generally is not the thread that consumes it;
private L2 residency does not carry cleanly between the phases.

Initial follow-up measurements confirmed that this working set mattered.
Changing `RSTRIP` from 256 to 64 improved best complete profiled time by 18-36%
at 1000² and 16-26% at 4000² across one and four threads. Generated outputs were
byte-identical at 1/2/4 threads and the five surface-metric tests passed.
`RSTRIP=64` is now provisional pending the full 324-raster exact-comparison
matrix.

The following serial-memory step removed the redundant `H` zero-fill, initialized
the independent output arrays concurrently, and copied the five required masks
concurrently. On seven unprofiled generated-raster runs, the resulting whole-fit
curve reached 1.79x at 1000² and 2.20x at 4000² on four threads. Outputs were
byte-identical to the `RSTRIP=256` baseline at 1/2/4 threads and all five
surface-metric tests passed. This remains provisional until the full 324-raster
matrix and uninitialized-memory checks are available.

Explicit OpenMP SIMD across the block-column loops was subsequently tested and
rejected. GCC 11 vectorized carry seeding, prefix construction and carry storage,
but not the nodata-guarded solve. After discarding a batch taken while the laptop
was in power-saving mode, nine scalar/SIMD pairs were run interleaved. Results
were inconsistent, and the important 4000² four-thread case regressed by 14%
on best time and 4% on median. Outputs remained byte-identical, but the SIMD
directives were removed under the no-regression rule. Do not repeat the same
directive-only experiment.

Before the serial-memory step, the whole-function fit timer included six serial
coefficient fills, one mask fill and five serial mask copies outside OpenMP:
approximately 30 bytes per pixel. The implementation described above now runs
those independent fills and copies concurrently when multiple threads are
requested; the allocations themselves remain serial so failures propagate
normally.

GCC 11 with the project's Release flags also leaves the inner prefix-build and
solve loops scalar. The prefix recurrence is serial down each column, but
columns and blocks are independent; at width 1000 there are about 125 blocks,
so the recurrence does not inherently serialize the OpenMP workers. Shared
cache pressure, scalar per-core execution and the timed serial fraction are
stronger candidates than a shortage of parallel blocks.

Current investigation order:

1. Repeat representative cases with `OMP_PLACES=cores` and
   `OMP_PROC_BIND=close|spread` to exclude WSL2 migration and sibling placement.
2. Complete `BLK` tuning around the provisional `RSTRIP=64`.
3. Check load balance on asymmetric nodata rasters.
4. Test padded carry storage last. Static scheduling limits cross-thread carry
   sharing to partition boundaries, and carries are written only once per strip.

These are outstanding work rather than findings, so they are tracked
elsewhere rather than here. This document records
what was measured; it is not a to-do list.

#### Sequential cases are no longer swept (changed 2026-09-01)

Depression filling (both algorithm families), least-cost breaching, and
D8/D-inf/MFD graph construction are **ordered traversals**: each step depends on
the result of the one before it. There is no parallelism in them to measure, in
lsdtt3 or in any of the references.

Sweeping threads over them measured nothing. It multiplied the campaign length
and — worse — filled the scaling table with flat rows that invite a reader to
conclude something is broken. These cases are marked `sequential` in the case
registry (list as of 2026-10-11) and are **measured once, at the baseline thread
count**:

`fill`, `fill_priority_flood`, `fill_priority_flood_radix`, `fill_epsilon_radix`, `fill_epsilon`, `fill_flats_float_step`, `fill_flats_gradient`, `breach`, `breach_least_cost`, `breach_least_cost_default`, `breach_then_fill_lindsay2016`, `breach_then_fill_least_cost`, `accum_d8`, `accum_d8_no_svector`, `route_d8_default`, `route_d8_epsilon`, `accum_dinf`, `accum_mfd`

They carry `sequential_case = True` in the CSV, are excluded from the scaling
table, and `SUMMARY.md` states why in that section. This cut the `laptop`
profile from 308 measurements to 200 with no loss of information.

The surface-metric cases stay in the sweep. They are genuinely parallel code, so
their poor scaling is a defect worth tracking rather than an expected result.

!!! warning "Do not add a case to `SEQUENTIAL_CASES` because its scaling is poor"
    The flag means "this algorithm has no parallelism to measure", not "this
    scales badly". Marking a genuinely parallel kernel sequential would hide
    exactly the finding the benchmark exists to surface.

Whether these algorithms *could* be parallelised bit-identically is a separate
question, and an open one.

#### Consequences for reading and running the benchmark

- **`fill` and `breach` are no longer swept at all**, for the reason above.
  Before 2026-09-01 they appeared in the scaling table at ~1.00×, which was the
  correct result rather than a finding.
- **Do not sweep hyperthreads.** On a bandwidth-bound kernel they cost 4× on
  the development machine. Cap `--threads` at the physical core count;
  `--cpus-per-task` should match, and the Slurm wrapper's `--exclusive` matters
  more here than for a compute-bound benchmark, because a co-tenant competes
  for the same memory bandwidth.
- **Scaling curves need the ceiling stated alongside them.** An lsdtt3 curve of
  1.17× means something quite different on a machine that can do 3.09× than on
  one that can do 1.14×. Re-run the two-kernel ceiling test on any machine
  whose scaling numbers are going to be quoted.

### What is compared, and what must never be

The suite (Simon, 2026-10-10; revised 2026-10-11) shows (a) which lsdtt3
methods are fastest, lsdtt3 against itself, and (b) lsdtt3 against other
packages only where the algorithm and the output are like-for-like.

1. **Surface fitting** -- `slope_only`, `polyfit_all`, `hillshade`
   (report-only, kept).
2. **Single-step fills** -- `fill` (vs TopoToolbox) and `fill_priority_flood`
   (vs RichDEM and Whitebox Wang & Liu). Real comparisons since 2026-10-11:
   same algorithm, outputs asserted identical (§4-5), even though the output
   keeps its flats. lsdtt3's `fill_priority_flood_radix` is lsdtt3-only.
3. **Breaching** -- `breach` (Lindsay 2016, vs RichDEM) and
   `breach_least_cost` (Lindsay 2020, vs Whitebox 2.0.6 under matched
   parameters, outputs asserted bitwise; re-instated 2026-10-11). lsdtt3-only:
   `breach_least_cost_default`, `breach_then_fill_lindsay2016`,
   `breach_then_fill_least_cost`.
4. **Connected flow networks** -- `fill_epsilon`, the D8 accumulation cases
   and `route_d8_epsilon` (raw DEM to D8 by the same pipeline in lsdtt3 and
   Whitebox: epsilon priority flood + steepest descent). lsdtt3-only:
   `fill_epsilon_radix`, `fill_flats_float_step`, `fill_flats_gradient`.
   `route_d8_default` (each package's default pipeline) is report-only.
   `SUMMARY.md` also ranks lsdtt3's own fill options against each other
   (below).
5. **Multi-direction flow** -- `accum_dinf`, `accum_mfd` (report-only).

| Suite | Cases | Status |
|---|---|---|
| `flow` | fill_epsilon, breach, accum_d8, accum_d8_no_svector | timing; values asserted separately in §4-5 |
| `flow` | fill, fill_priority_flood | timing; values asserted bit-exact in §4-5. Report-only 2026-10-10 .. 2026-10-11. Run at `fill_min_slope: 0`, so the output keeps its flats and is **not a connected flow network**: quote them as fill-algorithm comparisons, never as a drained-surface comparison |
| `flow` | breach_least_cost | timing; lsdtt3 least_cost (edge policy whitebox, max_dist 100, minimize_dist false, increment 0.001, no remainder fill) vs Whitebox 2.0.6 with the same parameters (flat_increment = float32(0.001)); outputs asserted bitwise (`benchmark_breach_least_cost_*`, `benchmark_breach_least_cost_fuzz`, `unit_breach_least_cost_real_terrain`) |
| `flow` | breach_least_cost_default, breach_then_fill_lindsay2016, breach_then_fill_least_cost | **timing only, report-only, lsdtt3-only** -- least_cost with the CLI defaults; and breach + barnes2014 fill (min_slope 1e-4) per breaching method, the time to a routable surface |
| `flow` | route_d8_epsilon | timing -- raw DEM to routed D8 accumulation: lsdtt3 drainage-analysis (Barnes 2014 priority flood, fill_min_slope 1e-4, steepest-descent D8) vs Whitebox Wang & Liu (fix_flats, flat_increment 1e-4 × cell size) + d8_flow_accum. Differs only in the diagonal increment (√2 in lsdtt3) |
| `flow` | fill_priority_flood_radix, fill_epsilon_radix | **timing only, report-only, lsdtt3-only** -- the `fill_priority_flood` and `fill_epsilon` lsdtt3 commands with `priority_flood_queue: radix` |
| `flow` | route_d8_default | **timing only, report-only** -- raw DEM to routed D8 accumulation using each package's default pipeline; algorithms and outputs differ by design |
| `flow` | fill_flats_float_step | **timing only, report-only, lsdtt3-only since 2026-10-11** -- fill + Barnes (2014b) flat resolution in float32 steps (`barnes2014b_float_step`). RichDEM FillDepressions + ResolveFlats was paired until then; it does one pass (bitwise equal to lsdtt3's round 1, asserted in §5), lsdtt3 loops to a routable surface |
| `flow` | fill_flats_gradient | **timing only, report-only, lsdtt3-only** -- fill + Barnes (2014b) flat resolution with a 1e-4 gradient (`barnes2014b_gradient`); no reference implements it |
| `mfd` | accum_dinf, accum_mfd | **timing only** -- Tarboton vs Quinn vs Freeman proportioning genuinely differ |
| `surface` | slope_only, polyfit_all, hillshade | **timing only** -- curvature definitions and aspect handedness differ |
| ~~`flow`~~ | ~~fill_whitebox_default~~ | **removed 2026-10-11** -- Whitebox's default `fill_depressions` (pit-centric, different surface) is not like-for-like with any lsdtt3 fill |

Two of these deserve their reasoning stated wherever they are quoted:

`polyfit_all` vs `slope_only` are both reported on purpose. lsdtt3 fits one
quadratic and derives every metric from it; TopoToolbox and Whitebox recompute
their derivatives per metric. `polyfit_all` is therefore where lsdtt3's design
pays, and `slope_only` is where it does not -- lsdtt3 pays for a full fit to
produce a single number. Quoting only the first would be cherry-picking. Note that the
window size does NOT tilt this: the `direct` fit's cost is independent of
`polyfit_lambda_metres` (measured flat over a 729× range of window cells).

**Pairings corrected 2026-09-13.** Three reference rows ran a different
algorithm from the lsdtt3 row they were compared against:

- `fill` is the grayscale-RECONSTRUCTION case, but Whitebox's `fill_depressions`
  is a priority-flood tool (its own tool description). Whitebox now appears in
  `fill_priority_flood`; `fill@whitebox` is `n/a` because whitebox_workflows
  2.0.6 has no reconstruction fill.
- `fill_epsilon@whitebox` ran `fill_depressions` (the tool that disagrees with
  three other implementations by up to 9.78 m, see §10o) at Whitebox's default
  `flat_increment` of 1e-4 m per cell, against lsdtt3's 1e-4 GRADIENT (3e-3 m
  per cell at 30 m) — 30x less flat work. It now runs
  `fill_depressions_wang_and_liu(fix_flats=True)` with
  `flat_increment = fill_min_slope * cell size`, read from each fixture.
- `breach@whitebox` ran `breach_depressions_least_cost` (Lindsay 2020) against
  lsdtt3's and RichDEM's Lindsay (2016). `breach` is now lsdtt3 vs RichDEM.
  A `breach_least_cost` case (lsdtt3 least_cost vs Whitebox) was added the same
  day. It was not like-for-like: lsdtt3 weights the cut by step distance while
  Whitebox ran `minimize_dist=False`, and the automatic flat increments differ
  by about seven orders of magnitude. Campaigns from 2026-09-13 to 2026-10-10
  carry those `breach_least_cost` rows; do not quote them. The case was
  **removed on 2026-10-10**: lsdtt3's least-cost breach outputs looked
  incorrect on inspection, and it returns only once that implementation has
  passed validation. Breaching is benchmarked on Lindsay (2016) only.

Campaigns before 2026-09-13 carry the old pairings in those three rows.


**Pairings corrected 2026-10-08.** `fill_priority_flood@topotoolbox` ran
`fillsinks(hybrid=False)`, which is Vincent's sequential reconstruction, not a
priority flood. TopoToolbox has no priority-flood fill, so it is now `n/a` in
that case. Campaigns before 2026-10-08 carry this row; it must not be read as a
priority-flood comparison.

`route_d8_default` (added 2026-10-08) answers the user's question rather than
the algorithm question: how long each package takes, with its defaults, to get
from a raw DEM to a routed D8 accumulation surface. The algorithms are not
matched and the rasters differ; no value is compared and no ratio is built.
Every pipeline must return a **fully drained** surface, because that is what
lsdtt3 provides. Checked 2026-10-08 on alpine_1000 and lakeland_1000: interior
cells with no strictly lower neighbour (TopoToolbox: no receiver in the
`FlowObject`) were 0 for all four engines, against 5,532 and 128,931 in the
raw DEMs.

| engine | fill | flats | physical slope on flats in the DEM? | D8 rule |
|---|---|---|---|---|
| lsdtt3 | `lsdtt-drainage-analysis`, `raster_is_filled: false`: `fill_pits` with `fill_min_slope` 1e-4 and `fill_flat_resolution_choice: inline_gradient`, i.e. Barnes et al. (2014) priority flood | gradient imposed inside the flood: each cell is raised to at least 1e-4 x step length (cell size, x sqrt2 on diagonals) above the cell it was reached from, so flats drain in parallel lines along the flood paths | yes, >= 1e-4 m/m | steepest descent (drop / distance) |
| TopoToolbox | `FlowObject(dem)` defaults (`method="d8"`, `sink_resolution="carve"`): `fillsinks_hybrid` (Vincent 1993 reconstruction, no gradient) | `identifyflats`, then gwdt costs from fill depth (`gwdt_computecosts`) and a geodesic-time auxiliary surface (`gwdt`); `flow_routing_d8_carve` routes flats and filled depressions by least-cost carving (Schwanghart et al. 2013, ESPL 38:1576-1586) | no: the DEM is not modified, flats stay at zero slope | steepest descent on the filled DEM; on flats, from the auxiliary surface |
| RichDEM | `FillDepressions(epsilon=True)`: Barnes (2014) Priority-Flood+Epsilon. `epsilon=True` is NOT the default; the default (`epsilon=False`) leaves flats that its D8 does not route | each flat cell raised to `nextafter` of the cell it was reached from, one float step | no: steps are ULP(z) of the raster dtype, slopes of order ULP(z) / cell size | the LOWEST neighbour, not steepest descent |
| Whitebox | `fill_depressions()` defaults: pit-centric (finds interior pits, sorts them, works pit by pit), `fix_flats=True`, `flat_increment` 1e-4 | `flat_resolution` `garbrecht_martz`: per flat, rank = dist_low x (max_high + 1) + (max_high - dist_high), raise = 1e-4 m x rank -- towards lower terrain first, away from higher terrain second | yes, in metres (1e-4 m per rank step), so a small real gradient | steepest descent (`d8_flow_accum`) |

Time is each engine's in-process call; for lsdtt3 the sum of its fill, D8
build and accumulation `Compute time:` lines. Whitebox's
`breach_depressions_least_cost` is not used: its default `fill_deps=False` can
leave depressions unresolved.

Where each row was checked (2026-10-10):

- **lsdtt3**: `cli/lsdtt-drainage-analysis.cpp` (default `fill_min_slope`
  0.0001; the fill block calls `lsdtt::fill_pits(dem, min_slope, bc, threads)`
  under `inline_gradient`); `src/cpu/lsdtt_raster_flow_preprocessing_cpu.cpp`
  (`fill_pits` sends `min_slope > 0` to `fill_pits_priority_flood`; the flood
  raises to `centre.zeta + min_slope * distance`, `distance` = L or L x sqrt2);
  `src/cpu/lsdtt_flowinfo_cpu.cpp` (receiver = largest `dz / path_length`).
- **TopoToolbox** 0.0.12, installed sources in the pixi `benchmark` env:
  `flow_object.py` (`FlowObject.__init__` defaults, `_d8_carve` ->
  `grid.auxiliary_topography` -> `_grid.flow_routing_d8_carve(..., filled_dem,
  aux, flats, ...)`) and `grid_object.py` (`auxiliary_topography`:
  `fillsinks_hybrid`, `identifyflats`, `gwdt_computecosts(costs, ..., dem,
  filled_dem)`, `gwdt`); no DEM is returned. Steepest descent off flats is
  asserted against lsdtt3 in §4.
- **RichDEM** 2.3.0: `richdem/__init__.py` (`FillDepressions(dem,
  epsilon=False, ...)`; `epsilon=True` -> `rdPFepsilonD8`), the installed
  headers `include/richdem/depressions/Barnes2014.hpp`
  (`PriorityFloodEpsilon_Barnes2014`, `std::nextafter`) and
  `flowmet/OCallaghan1984.hpp` (`FM_OCallaghan`: keeps the neighbour with the
  lowest elevation). `tests/scripts/run_flow_benchmark.py` records the same
  lowest-neighbour finding.
- **Whitebox** 2.0.6: the stub `whitebox_workflows.pyi` in the pixi
  `benchmark` env (`fill_depressions(dem, fix_flats=True,
  flat_increment=0.0001, max_depth=None, flat_resolution=None, ...)`), and
  `../whitebox_next_gen/crates/wbtools_oss/src/tools/hydrology/mod.rs` at
  commit cd24675: `parse_flat_resolution_mode` (`None` ->
  `"garbrecht_martz"`), the tool defaults (`fix_flats` true, `flat_increment`
  0.0001, `flat_resolution` `garbrecht_martz`), `fill_depressions_core` (pits
  collected in the interior, sorted, popped one at a time) and
  `resolve_flats_garbrecht_martz` (the rank formula above);
  `tools/flow_algorithms/mod.rs` `d8_dir_from_dem`, used by
  `D8FlowAccumTool` (largest `(z0 - z) / length`). Steepest descent was also
  measured, §10o (`tests/scripts/check_whitebox_d8_rule.py`).

`fill_spill_merge` (lsdtt3's nested depression hierarchy, Barnes 2020, plus
Fill-Spill-Merge, Barnes 2021, lsdtt3 only) was **removed on 2026-10-10**: the
suite now compares like-for-like algorithms only, and no binding benchmarked
here offers Fill-Spill-Merge (RichDEM's Python package 2.3.0 exposes only
`FillDepressions` and `BreachDepressions`). Until 2026-09-13 it had been paired
with Whitebox's `sink`, which only labels depression cells.

**lsdtt3 fill ranking (added 2026-10-10).** `SUMMARY.md` has a section "lsdtt3
fill methods, fastest first": per fixture, at the baseline thread count, the
lsdtt3 rows of `fill`, `fill_priority_flood`, `fill_priority_flood_radix`,
`fill_epsilon`, `fill_epsilon_radix`, `fill_flats_float_step`,
`fill_flats_gradient`, sorted by median,
with a one-phrase algorithm, whether the output is a connected flow network,
the median, Mcells/s and the ratio to the fastest. Only `cli_reported_compute`
rows with status ok enter it. It compares lsdtt3 with itself on one machine, so
the report-only flags of the cross-package table do not apply. The radix variants are
honoured by the CLI for `barnes2014` under the default `inline_gradient` flat
resolution: `fill_pits_priority_flood` receives the queue backend directly at
`fill_min_slope` 0 and > 0, inside the single `Compute time:` bracket.

Two cases reuse rows rather than re-measuring (`reuse_from` in the case
registry, `reused_from` in the CSV). `accum_d8_no_svector` copies the
reference rows of `accum_d8`, because the references run the identical call in
both. `fill_whitebox_default` copies lsdtt3's `fill_priority_flood` row, because
its lsdtt3 side is the same barnes2014 command.

**Pairings corrected 2026-10-10.** `fill_priority_flood@whitebox` ran
`fill_depressions`, which is pit-centric (finds interior pits, sorts them, works
pit by pit) rather than an edge-seeded priority flood. It returns a different
surface (6.45 m high on alpine_small, §10o) and its cost depends strongly on the
number of pits: 0.63 s on alpine_4000 against 54.7 s on lakeland_4000
(forth01, 2026-09-21, 1 thread). The row now runs Whitebox's
`fill_depressions_wang_and_liu(fix_flats=False)`, an edge-seeded priority flood
that is bit-exact against lsdtt3 (3.1 s alpine, 4.2 s lakeland on the same
campaign). The former `fill_wang_liu` case made exactly that pairing and was
removed; `fill_depressions` is kept, report-only, as `fill_whitebox_default`.
Campaigns before 2026-10-10 carry the old `fill_priority_flood@whitebox` row; do
not quote it as a priority-flood comparison. The CSV also gained
`seconds_build_only` / `seconds_accumulate_only`: TopoToolbox's `FlowObject`
refills its input and runs its gwdt flat routing even on the pre-filled
`accum_d8` DEM, so its accumulate-only time is the like-for-like part of that
row.

**Memory.** Every command is launched through `tests/scripts/rusage_launcher.py`.
A Linux child inherits its parent's lifetime peak RSS and `wait4` reports it as
the child's own. Measured on WSL2: after a Python parent touched and freed
1.5 GiB, `/bin/true` launched from it reported 1509.5 MiB; from a fresh process,
9.0 MiB. The 2026-09-09 forth02 campaign launched directly from a driver that
had peaked at 3.5 GB in `describe_fixture()`, and 441 of its 504 rows read
exactly 3556.98 MiB. Memory columns from campaigns before 2026-09-13 are invalid
except where a value exceeds that floor. Each campaign now records `rss_floor`
(`/bin/true` through the launcher: 11.0 MiB on WSL2, against 93.0 MiB launched
directly) and fails a tripwire above 64 MiB. **Confirmed on the server
2026-09-13 (forth01):** floor 10.5 MiB through the launcher against 942.5 MiB
launched directly; 233 distinct peak values across 318 rows, none near the
floor, scaling with grid size.

### Accumulation runs on a shared, pre-filled DEM

Flow accumulation cannot run on a raw DEM -- `lsdtt-drainage-analysis` aborts
with *"could not assign a receiver ... unresolved sink or flat"*, and each other
engine would silently apply its **own** depression handling first, turning an
accumulation comparison into a comparison of fills. The suite therefore fills
once and every engine accumulates over the identical raster.

lsdtt3's fill is the shared input, which needs justifying since lsdtt3 is one of
the tools under test: it is bit-exact against both RichDEM's Barnes 2014 fill
and TopoToolbox's `fillsinks` (§4), so it favours nobody. The fill uses
`fill_min_slope: 0.0001`, not `0` -- a fill with true flats is unroutable,
because a flat cell has no downslope neighbour for D8 to pick.

### A Whitebox fill result we could not explain

On `alpine_small` (256², no nodata) the raw minimum 1177.1334 sits at row 254,
an **interior** cell, and the minimum along the border is 1178.3820. A
depression fill should therefore bring the global minimum to the border
minimum. RichDEM and Whitebox's `fill_depressions_wang_and_liu` both return
exactly 1178.3820. Whitebox's `fill_depressions` -- its Barnes tool -- returns
**1184.8361**, 6.45 m higher, and the two Whitebox fills differ by up to 9.78 m.

lsdtt3 is bit-exact against RichDEM and TopoToolbox on this operation, so the
odd one out is Whitebox's Barnes tool, most likely a boundary convention
(edge cells treated as closed). Until that is understood, `fill@whitebox` is a
**timing** row and `depression_fill_wl` carries the Whitebox value comparison.

### Whitebox licensing

The suite constructs `wb.WbEnvironment()`, the **open tier**, which needs no
licence key. 733 of Whitebox's 775 tools are open tier and every routine used
here is one of them (verified against 2.0.6). Never add a routine requiring a
Pro entitlement, and never call `from_floating_license_id` or
`from_signed_entitlement_*`: a benchmark a reviewer cannot run without buying a
licence is not a benchmark.

Whitebox is a ~22 MB wheel, ~60 MB installed, with **zero Python
dependencies**, so it cannot perturb the numpy/GDAL solve. It has **no osx-64
wheel**, which is why the `whitebox` feature declares only
`linux-64, osx-arm64, win-64` and the `benchmark-lite` environment exists as the
Intel-Mac route:

```sh
pixi run -e benchmark-lite python tests/scripts/run_benchmark_suite.py --profile laptop
```

Drivers skip any reference they cannot import, so `benchmark-lite` simply
produces fewer rows.

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
ctest --test-dir build-perf -L '^performance$' --output-on-failure
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

- **Pin the threads.** `OMP_NUM_THREADS=8 ctest --test-dir build-perf -L '^performance$'`.
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

TopoToolbox's *non*-hybrid path was measured tracking lsdtt3's priority-flood
almost exactly at every size, and was read at the time as a control showing
that implementation quality was not the differentiator. That reading is
withdrawn: `hybrid=False` is Vincent's *sequential* reconstruction, not a
priority flood (TopoToolbox has none), so the two curves come from different
algorithms and their agreement is not a control. The measurement stands. The
conclusion that the hybrid reconstruction is the faster algorithm class rests on
the matched `vincent1993` vs `hybrid=True` comparison, not on this.

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
The opt-in Barnes et al. (2014b) flat resolution (§10r) is a separate code path
in `lsdtt_raster_flat_resolution_cpu`: it does not use this fixed point.

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

### The D8 flow object build: refactored for speed, byte-identical (2026-09-13)

`lsdtt_flowinfo::create()` builds every D8 flow object in lsdtt3 (every CLI,
and every LEM timestep that re-routes). On the 2026-09-09 forth02 campaign,
lsdtt3's D8 drainage area as shipped (s-vector on) took 4.0-5.0x as long as
Whitebox's `d8_flow_accum`, and 2.2x with the s-vector off. Profiling found
that most of the difference was bookkeeping rather than the algorithm. It was
refactored under one rule: **not one output byte may change.**

**The pin, written first.** `unit_flow_d8_identity`
(`tests/cpp/test_flow_d8_identity.cpp`) was committed (`c9aecb1`) BEFORE any
change. It pins hashes of every array the build produces: receivers,
direction codes, donor counts, outlets, topological order, the full s-vector
with its start/size arrays, drainage area, and a weighted accumulation (the
path the cosmogenic code uses). The fixtures are 17 surfaces covering every
boundary code (`o`, `n`, `p`, `b`, mixed), forced outlets, pits as outlets,
nodata holes with nodata-as-outlet on and off, a nodata border, a non-integer
cell size, the periodic 2-row and 1-column degenerate cases, and 1x1 / 2xN
grids. Each is built in 4 option configurations at 1 and 3 threads. It also
asserts that drainage area is identical with the s-vector on and off, and that
each larger fixture actually contains tied steepest slopes (76-383 cells), so
the tie-break is exercised.

The fixtures are integer terrain (`z = 1000 * steps-to-outlet + {0,250,500,750}`),
routable by construction and exactly reproducible on any machine:
`tests/data/dem_cache` is gitignored, and libm-built surfaces can differ in the
last bit across platforms. Mutation-tested: flipping the receiver tie-break
`>` to `>=` fails 277 checks; reordering s-vector donors fails exactly the
s-vector pins.

**What changed (`eb65499`), each change run against the pin on its own:**

| change | why it is identical by construction |
|---|---|
| interior cells find neighbours by precomputed flat offset, skipping the per-direction boundary-condition resolution | an interior cell's neighbours are inside the grid under every boundary code; same comparisons in the same direction order, so the strict `>` first-wins tie-break is unchanged |
| donor counting and outlet collection merged into one pass | same increments; outlets still row-major |
| the Kahn queue is used as the topological order | the queue's contents in push order ARE the order; it used to be copied cell by cell into a second n-cell vector |
| `accumulate()` seeds weights, and `accumulate_raster()` copies, in row-major order | the topological order holds exactly the valid cells, so the same values land in the same cells; no addition is reordered |
| `class lsdtt_flowinfo final` | removes indirect calls to the virtual index helpers from the build loops (objdump: 0 left) |
| s-vector donor search on interior cells visits neighbours in ascending flat order | yields the donors already sorted and duplicate-free, which the general path achieves by dedup and insertion sort |
| `s_vector_start_`/`s_vector_size_` and `donor_count_` no longer filled twice | the second fill wrote the same values |

The mutation tests were repeated against the new interior paths (tie-break:
252 checks fail; donor order: 22 s-vector checks fail). `ctest -L flow -LE heavy`
18/18, `ctest -L network -LE heavy` 92/92.

**`validate()` moved from every build to CTest.** `create()` ended with
`validate()`, a full audit of the object it had just built. None of its 36
checks can fail on a correct build, whatever the DEM: pits, flats, cycles and
bad inputs are caught earlier by `build_receiver_graph_d8`,
`build_topological_order`, `build_s_vector` and `validate_inputs`, which are
unchanged. So it only ever detected bugs in `lsdtt_flowinfo_cpu.cpp`. It cost
time on every D8 build (see the laptop note below) and ~10 B/cell of transient
peak memory (an
8 B/cell `position` array plus two 1 B/cell planes). It is now called
explicitly by `unit_flow_d8_identity` (every fixture, configuration and thread
count) and `unit_flow_properties` (the fractal, nodata and pit surfaces, with
the s-vector on and off), and no longer by `create()`. Checked by injecting a
donor-count bug: `create()` no longer throws on it, and the CTest `validate()`
checks fail on it. `validate()` remains public, documented as test-time only.

**Rejected, and why:**

- *Stack (depth-first) accumulation, as Whitebox does.* It was the largest
  single gain profiled, but it adds a receiver's donors in a different order,
  so the s-vector-off path would compute drainage area in a different order
  from the s-vector-on path. lsdtt3 must not behave differently depending on
  the s-vector switch, and `accumulate(weights)` feeds the cosmogenic
  calculations, where addition order can change the last bits.
- *32-bit flat indices.* 64-bit indices are a deliberate design choice so very
  large rasters work: 32-bit caps a raster at 2,147,483,647 cells (46,340²),
  and a 50 x 50 km 1 m lidar mosaic is 2.5e9 cells. It was also a small gain.
- *Single-precision slopes.* Changes output bytes.

**Speed (server).** Measured in the 2026-09-13 forth01 campaign against the
2026-09-09 forth02 campaign; see "Server campaign 2026-09-13" below for the
method and caveats. Single thread, old time / new time: D8 drainage area as
shipped (s-vector on) **2.07-2.38x** faster, s-vector off **1.33-1.71x**,
growing with grid size. forth01 ran unchanged code 5-24% slower than forth02,
so these understate the gain slightly.

!!! warning "Laptop figures only -- not results"
    The only measurements so far are from a thermally limited laptop (4-core,
    15 W, WSL2), where speeds are unreliable and drift with temperature. They
    show the direction of the change and nothing more. With that warning:
    against an in-process copy of the old code, the new build was 1.3x faster
    at 1000² and 1.5-1.6x at 4000² with the s-vector off, while still running
    `validate()`. `validate()` alone took roughly 10-16% of a build. Against
    the old library in earlier runs, roughly 1.5-1.9x
    with the s-vector on or off. The receiver build remains the largest phase,
    and is slower on flat terrain (data-dependent branching in the slope
    comparison).

**The junction network followed (2026-09-13).**
`lsdtt_junction_network_cpu::create()` also ended with `validate()`. It was a
much smaller share of that build than the flow object's audit (laptop only,
thermally limited, direction only: roughly 1-11% of the junction network build
and ~1-2% of flow plus network), but the reasoning is the same, so it was moved
the same way. `unit_flow_d8_identity` now builds junction networks on every
flow fixture, at dense, moderate and sparse channel thresholds and with each
option combination, and calls `validate()` on each. Checked by injecting a
junction donor-count bug: 6 fixtures fail.

**Still open:** the receiver-build branch cost on flat terrain (any change
must keep the strict first-wins tie-break). The junction network BUILD itself,
not its audit: with a dense channel network it cost more than the D8 build
(laptop, same warning as above: 4.0 s against 1.6 s on a 4000² alpine grid at a
50-pixel threshold).

### Server campaign 2026-09-13 (forth01): speedups confirmed, thread scaling not measured

`server` profile (1000/2000/4000², threads 1/2/4/8, 5 repeats) at `39b5293`
(dirty tree), on **forth01**. The previous campaign (2026-09-09) ran on
**forth02**. Same CPU model (AMD EPYC 9654P, one socket, 192 logical cores),
compiler (GCC 13.3.0, Release), reference package versions, and byte-identical
fixtures (same SHA256 prefixes). Neither run was an exclusive allocation. 318
rows, all `ok`. Raw data: `../data/benchmarking/server_forth01.geos.ed.ac.uk_20260913_162715`.

**The machines were not equivalent during the run.** The reference engines did
not change between campaigns, so their movement measures the machine:

| drift control (old time / new time; < 1 = forth01 slower) | median | range |
|---|---|---|
| RichDEM (serial) | 0.89 | 0.76-1.00 |
| TopoToolbox (serial) | 0.92 | 0.82-1.00 |
| Whitebox, serial tools (D-inf, MFD, Wang & Liu) | 0.75-0.93 | |
| Whitebox, internally multi-threaded tools (slope, hillshade, curvature, D8) | 0.25-0.63 | |
| lsdtt3, code unchanged since 09-09 (fill, breach, fill_epsilon, fill-spill-merge, D-inf, MFD; 1 thread) | 0.85 | 0.76-0.95 |

The slowdown was worst on the largest grids, which points to memory-bandwidth
contention.

**Thread scaling was NOT measured.** Every multi-threaded measurement in the
campaign took the same time at 2, 4 and 8 threads as at 1. Process CPU time
divided by wall time was 1.3-1.7 at every thread count, for lsdtt3 and for
Whitebox alike. For example, `polyfit_all` on alpine_4000 took 2.29 s at 1 thread
and 2.18 s at 8 (forth02: 2.13 s and 0.62 s). Something limited the whole job to
about 1.5 cores. Candidates, not yet resolved: a per-user CPU limit or load on
forth01, or (not ruled out) an interaction with the `rusage_launcher.py` launch
path, which is the only harness change since 09-09. No scaling figure from this
campaign may be quoted, and no Whitebox timing from its multi-threaded tools.

**lsdtt3 speedups (single thread, old time / new time, all six fixtures):**

| case | speedup | code change |
|---|---|---|
| `accum_d8` (s-vector on, as shipped) | **2.07-2.38x** | D8 flow object refactor, `validate()` removed from `create()` |
| `accum_d8_no_svector` | **1.33-1.71x** | same |
| `hillshade` | **1.56-1.65x** | no input copies, invalid-cell plane |
| `slope_only`, `polyfit_all` (1 thread) | 0.90-0.94 | unchanged; matches the machine |

**Where lsdtt3 now stands (same run, 1 thread; reference time / lsdtt3 time):**

| case | RichDEM | TopoToolbox | Whitebox (serial tools only) |
|---|---|---|---|
| D8 as shipped | 1.04-1.17 | 1.56-1.78 | (multi-threaded, not comparable) |
| D8, s-vector off | 1.41-1.71 | 2.09-2.43 | (multi-threaded, not comparable) |
| D-inf | 0.93-1.16 | n/a | 0.70-0.89 |
| MFD (Quinn vs Freeman) | 0.64-0.66 | n/a | 0.53-0.63 |
| fill, reconstruction | - | 0.80-1.09 | n/a |
| fill, priority-flood | 0.50-0.55 alpine, 0.64-0.79 lakeland | 0.71-0.95 alpine, 3.3-9.5 lakeland | 0.13-0.19 alpine, 2.4-18 lakeland |
| fill_epsilon | 0.94-1.04 | n/a | 0.88-1.13 alpine, 1.18-1.71 lakeland |
| breach (Lindsay 2016) | 0.96-1.06 | - | n/a |
| breach, least-cost | n/a | - | 0.53-1.27 alpine, 0.20-0.21 lakeland |
| hillshade | n/a | 0.51-0.56 | (multi-threaded, not comparable) |
| slope_only (lsdtt3 at 1 thread) | 0.43-0.49 | 0.95-1.01 | (multi-threaded, not comparable) |
| polyfit_all (lsdtt3 at 1 thread; metric sets differ) | 0.20-0.26 | 1.05-1.21 | (multi-threaded, not comparable) |

On 2026-09-09 the D8 build as shipped was 0.46-0.48x RichDEM and 0.60-0.73x
TopoToolbox (4000²); it now beats both. Hillshade went from 0.31x TopoToolbox
to 0.55x. The surface-metric ratios against RichDEM are single-thread lsdtt3
against serial RichDEM; lsdtt3's multi-threaded figures would change them, but
this campaign could not measure those.

**Memory (peak RSS of the whole process, 4000²; reference MiB / lsdtt3 MiB, > 1 =
lsdtt3 uses less).** Python reference processes carry ~100-200 MiB of
interpreter and NumPy before any work, so the small-grid ratios flatter lsdtt3;
4000² is the rung to read.

| case | RichDEM | TopoToolbox | Whitebox |
|---|---|---|---|
| D8 as shipped (lsdtt3 1015 MiB) | 0.99 | 1.33 | 1.33-1.34 |
| D8, s-vector off (769 MiB) | 1.31 | 1.75 | 1.76 |
| D-inf (1320-1364 MiB) | 0.74-0.76 | n/a | 0.99-1.02 |
| MFD (2734-2736 MiB) | 0.37 | n/a | 0.49 |
| fill, priority-flood (570-585 MiB) | 0.61-0.62 | 0.81-0.84 | 3.03-3.11 |
| fill_epsilon (570-585 MiB) | 0.59 | n/a | 3.01-3.08 |
| breach (539-555 MiB) | 0.76-0.78 | - | n/a |
| breach, least-cost (434-533 MiB) | n/a | - | 3.33-4.07 |
| hillshade (285 MiB) | n/a | 1.57 | 4.67 |
| slope_only / polyfit_all (742-743 MiB) | 0.45 / 0.59 | 1.01 / 1.69 | 1.79 / 2.29 |

The s-vector costs 246 MiB at 4000² (1015 against 769 MiB), 15.4 B/cell of
peak memory. This is the first campaign in which that is visible.

**To do before quoting scaling:** establish why the job had ~1.5 cores (check
`cpu.max` / `Cpus_allowed_list` / load on the node, and one direct
`/usr/bin/time -v` run of a surface CLI at `n_threads: 8` outside the harness),
then re-run `--suites surface` where 8 threads really get 8 cores.

> **Resolved 2026-09-21.** The cause was a login-node CPU affinity mask, not a
> code or harness fault. See the next section.

### Server campaign 2026-09-21 (forth01, release 0.5.2): thread scaling measured

`server` profile (1000/2000/4000², threads 1/2/4/8, 5 repeats) at `b8b6cdf`,
**clean tree**, on **forth01**. `b8b6cdf` is the commit tagged `v0.5.2`, so
this is the first campaign that measures a tagged release rather than a
development tree. 318 rows, all `ok`. Raw data:
`../data/benchmarking/server_forth01.geos.ed.ac.uk_20260921_153611`.

Note that the suite records the git SHA and not the release version, so the
"0.5.2" in this heading is asserted here by hand; `git rev-parse
"v0.5.2^{commit}"` resolves to `b8b6cdf` and is the check.

#### What had been limiting the earlier campaigns

An interactive session on forth01 is restricted to **two CPUs** —
`Cpus_allowed_list: 1,97`, the two SMT siblings of a single physical core —
while `lscpu` and `os.cpu_count()` both report 192. forth01 is a Slurm head
node and this is ordinary login-node policy. Every campaign run from an
interactive shell therefore swept 1/2/4/8 threads on two CPUs, which is why the
09-13 job appeared to have "about 1.5 cores".

It was identified with `tests/scripts/omp_ceiling.cpp`, which reported
`logical cores reported by OpenMP: 2` and an identical 1.62× at **both** 2 and
4 threads for a pure-arithmetic kernel containing no lsdtt3 code. In an
allocation the same kernel gives 2.00× and 3.99×.

Two consequences, both acted on:

- The suite's oversubscription guard used `os.cpu_count()`, which ignores
  affinity. It now uses `len(os.sched_getaffinity(0))`, records `usable_cores`
  and `affinity_restricted` alongside `logical_cores`, and warns in
  `SUMMARY.md`. **This campaign predates that change**, so its manifest carries
  no `usable_cores` field; the evidence that it had real cores is in the data
  below rather than in the manifest.
- `docs/what_is_lsdtt3/running_benchmarks.md` documents the check and the
  `srun --cpus-per-task=8` route.

#### Thread scaling (lsdtt3, alpine; lakeland within 0.03× throughout)

| case | rung | 1t | 2t | 4t | 8t |
|---|---|---|---|---|---|
| `slope_only` | 1000² | 1.00× | 1.85× | 1.96× | 0.50× |
| | 2000² | 1.00× | 1.87× | 1.96× | 2.13× |
| | 4000² | 1.00× | 1.87× | **3.18×** | **3.62×** |
| `polyfit_all` | 1000² | 1.00× | 1.77× | 1.99× | 1.01× |
| | 2000² | 1.00× | 1.91× | 2.02× | 2.74× |
| | 4000² | 1.00× | 1.85× | **3.22×** | **3.48×** |
| `hillshade` | 1000² | 1.00× | 1.50× | 1.49× | 0.76× |
| | 2000² | 1.00× | 1.52× | 1.61× | 1.04× |
| | 4000² | 1.00× | 1.52× | **2.05×** | **2.09×** |

The previous campaign measured 1.07× at 4000², four threads, for
`polyfit_all`. The polyfit kernel is byte-identical between the two runs — the
only source difference between `a8529f6` and `b8b6cdf` in
`lsdtt_raster_surface_metrics_cpu.cpp` is two comment paths.

Scaling **improves with grid size**, which reverses the behaviour recorded for
the historical two-pass implementation. The small rungs still lose at 8
threads (0.50× at 1000²); there is not enough work per thread to amortise the
parallel region, and this is the expected shape rather than a finding.

#### The 1-thread rows are the control

lsdtt3's single-thread times are unchanged from the affinity-limited campaign
six hours earlier — `polyfit_all` 4000² 2.042 s → 2.064 s, `slope_only` 0.805 s
→ 0.805 s. Across all 1-thread rows over 50 ms, every lsdtt3, RichDEM and
TopoToolbox case moved by ≤1.5% (median). Same machine, same binaries, same
work: only the threaded rows changed.

Peak RSS is identical to the previous campaign — **zero rows differ by more
than 10%**, on any engine — so the memory table from 09-21 10:43 stands
unchanged, and memory figures are insensitive to the affinity problem.

#### Whitebox's thread pinning does not take effect

The 09-13 section already refused to quote Whitebox's multi-threaded tools.
This campaign gives the mechanism. `whitebox_reference.py:504` sets
`wbe.max_procs = int(args.threads)` and the suite does pass `--threads 1`, but
five cases ran **40-72% faster** at "1 thread" once the allocation had more
CPUs, while every other engine reproduced to within 1.5%:

| case | 1-thread change, 2 CPUs → allocation |
|---|---|
| `slope_only` | −72.3% |
| `hillshade` | −59.8% |
| `polyfit_all` | −55.3% |
| `accum_d8`, `accum_d8_no_svector` | −41.9% |

Whole-process CPU seconds are unchanged across both runs and both thread counts
(22.9 s for `slope_only` 4000² at 1 and 2 threads, in both campaigns). The work
is constant; only the wall time moved, because more cores were available to
absorb it. `max_procs` is not constraining these tools in
whitebox_workflows 2.0.6.

**`accum_d8` and `accum_d8_no_svector` are in the like-for-like ratio table**,
where Whitebox now reports 15.95-24.68 Mcells/s against lsdtt3's 9.25-11.75.
That comparison is one-thread lsdtt3 against many-thread Whitebox and **must
not be quoted**. The other three affected cases are `report_only` and carry no
ratio. Until the pinning is understood, treat every Whitebox row for slope,
hillshade, curvature and D8 accumulation as an unpinned measurement.

**Fixed 2026-10-11 by pinning the process, not the setting.** The forth01
campaign of that date showed the same thing from the other side: on a 2-CPU
mask, Whitebox used 1.30 (`accum_d8`), 1.48 (hillshade), 1.54 (`polyfit_all`)
and 1.65 (slope) CPUs per process at `max_procs=1`. `run_benchmark_suite.py`
now launches every Whitebox measurement under `taskset -c <cpus>`, with as many
CPUs as the row's thread count, one per physical core first
(`whitebox_cpu_pin()`). The mask is inherited by every thread Whitebox starts.
A laptop check at 1000² after the change (laptop, thermally limited: the
direction is meaningful, the values are not): Whitebox slope 130 ms at 1 thread and
82 ms at 2, hillshade 91 and 53 ms, so the tools do scale once the 1-thread row
really is one core. `REFERENCE_MAX_THREADS` was emptied at the same time,
because the cap at 2 threads rested on flat curves measured on a restricted
mask. `SUMMARY.md` gained a CPU-use table that lists any engine and case
above 1.2 CPUs at the baseline thread count, so an unpinned row cannot pass
unnoticed.

### Server campaign 2026-10-11 (forth01): Whitebox pinned, least-cost breach back

`server` profile (1000/2000/4000², 5 repeats after a warm-up, medians) at
`4d8d844`, clean tree, `--threads 1,2,4`. Raw data:
`../data/benchmarking/server_forth01.geos.ed.ac.uk_20261011_140209`. 372 rows,
all `ok`.

- **Machine:** forth01 (AMD EPYC 9654P), GCC 13.3 `-O3`. The session had 8
  logical CPUs (`Cpus_allowed_list: 75,90-92,171,186-188`), which are **4
  physical cores** and their SMT siblings. **Not exclusive**; governor
  `schedutil`.
- **Whitebox pinned:** every Whitebox measurement ran under `taskset`, one
  physical core per thread. The CPU-use table flags nothing: every 1-thread row
  of every engine is at about 1.0 CPU/wall.
- **Earlier run the same day** (`…_20261011_120654`, `a3d90e6`) had a 2-CPU
  mask and unpinned Whitebox. Its lsdtt3, RichDEM and TopoToolbox 1-thread rows
  agree with this campaign within 0.89-1.08 (median per case). Its Whitebox
  rows do not: pinning made Whitebox's 1-thread times longer by 1.42×
  (`accum_d8`, hillshade), 1.32× (`route_d8_epsilon`), 1.25× (fills), 1.19×
  (slope, `route_d8_default`), 1.10-1.12× (polyfit, `breach_least_cost`) and
  1.04-1.08× (D-inf, MFD). **Whitebox rows from any campaign before this one
  are unpinned and must not be compared with these.**

#### 1-thread comparison

Reference time / lsdtt3 time, range over the six fixtures (above 1 means lsdtt3
took less time). Report-only cases compare different work.

| case | RichDEM | TopoToolbox | Whitebox |
|---|---|---|---|
| `accum_d8` (S-vector on) | 1.05-1.23 | 1.21-1.72 | 0.95-1.17 |
| `accum_d8_no_svector` | 1.50-1.74 | 1.77-2.26 | 1.32-1.87 |
| `fill_epsilon` | 0.92-0.99 | n/a | 1.03-2.75 |
| `fill_priority_flood` | 0.52-0.80 | n/a | 1.05-2.10 |
| `fill` (reconstruction) | n/a | 0.84-1.03 | n/a |
| `breach` (Lindsay 2016) | 0.94-1.49 | n/a | n/a |
| `breach_least_cost` (bitwise-equal output) | n/a | n/a | 0.69-0.83; 1.45 on alpine 4000² |
| `route_d8_epsilon` | n/a | n/a | 0.99-2.00 |
| `accum_dinf` (report-only) | 0.93-1.21 | n/a | 0.69-0.89 |
| `accum_mfd` (report-only) | 0.57-0.63 | n/a | 0.51-0.63 |
| `slope_only` (report-only) | 0.43-0.51 | 0.85-1.02 | 2.11-2.47 |
| `hillshade` (report-only) | n/a | 0.45-0.59 | 9.7-12.2 |
| `polyfit_all` (report-only) | 0.20-0.26 | 1.04-1.20 | 0.97-1.14 |

Raw DEM to D8 accumulation with each package's defaults (`route_d8_default`,
report-only), seconds:

| package | alpine 1000² | lakeland 1000² | alpine 2000² | lakeland 2000² | alpine 4000² | lakeland 4000² |
|---|---|---|---|---|---|---|
| lsdtt3 | 0.241 | 0.240 | 1.103 | 1.100 | 4.909 | 5.310 |
| RichDEM | 0.247 | 0.245 | 1.150 | 1.097 | 5.095 | 5.272 |
| TopoToolbox | 0.147 | 0.285 | 0.578 | 1.335 | 2.790 | 6.729 |
| Whitebox | 0.127 | 0.494 | 0.446 | 6.648 | 2.361 | 81.273 |

The timings of every lsdtt3 conditioning option on the same campaign, and the
decomposition of the TopoToolbox / lsdtt3 difference on alpine terrain (fill
family about 1.7 s, S-vector about 0.5 s at 4000²), are in
`docs/design_and_algorithms/flow_routing.md`, "From a raw DEM to a fully routed
DEM".

#### Least-cost breaching

`breach_least_cost` runs lsdtt3 and Whitebox with matched parameters (edge
policy `whitebox`, max_dist 100, increment float32(0.001), no remainder fill);
the outputs are asserted bitwise in CTest. Seconds:

| | alpine 1000² | lakeland 1000² | alpine 2000² | lakeland 2000² | alpine 4000² | lakeland 4000² |
|---|---|---|---|---|---|---|
| lsdtt3 | 0.040 | 0.263 | 0.135 | 1.051 | 0.498 | 4.695 |
| Whitebox | 0.028 | 0.199 | 0.093 | 0.779 | 0.723 | 3.920 |
| lsdtt3 / Whitebox | 1.41 | 1.32 | 1.45 | 1.35 | 0.69 | 1.20 |
| lsdtt3, CLI defaults (`breach_least_cost_default`) | 0.035 | 0.196 | 0.133 | 0.931 | 0.491 | 10.055 |

On identical output lsdtt3 took 1.20-1.45× Whitebox's time on five of six
fixtures. A laptop profile (direction only) put 87-90% of lsdtt3's time in the
search loop and found no algorithmic difference; heap reuse, int16 path
lengths, no validity lookup, row/col carried in the heap node and
`-march=native` each changed nothing beyond noise. The compiler (GCC here,
LLVM for Whitebox's Rust) is the remaining untested candidate. Peak memory at
4000²: lsdtt3 644-657 MiB, Whitebox 1759-1774 MiB.

With the CLI defaults (edge policy `outlet`, max_dist 1000, automatic
increment, remainder fill on) lakeland 4000² took 10.1 s against 4.7 s for the
matched run: the longer searches cross large lakes.

#### Thread scaling on 4 physical cores

Speedup over each engine's own 1-thread time, 4000² (alpine and lakeland
within 0.07× of each other):

| case | lsdtt3 2t | lsdtt3 4t | Whitebox 2t | Whitebox 4t |
|---|---|---|---|---|
| `slope_only` | 1.80-1.82× | 2.14-2.21× | 1.00× | 2.40-2.42× |
| `polyfit_all` | 1.83-1.85× | 2.12-2.13× | 1.02-1.03× | 1.80-1.83× |
| `hillshade` | 1.48× | 1.50-1.53× | 1.01× | 2.05× |

- Whitebox gains nothing from a second core and then 1.8-2.4× from four. Its
  2-thread rows here used 0.99-1.08 CPU/wall. The 0.6.0 rerun the next day
  (below) measured 1.51-1.81× at 2 threads, with 1.39-1.57 CPU/wall, on the
  same suite code, so the flat 2-thread curve belongs to this session, not
  to Whitebox.
- lsdtt3's 4-thread figures are below the 3.18-3.22× of the 2026-09-21
  allocation. lsdtt3's OpenMP threads were not pinned here, so two can share a
  physical core's SMT siblings, and the node was shared. Not separated.
- lsdtt3's 4-thread rows at 1000² are bimodal (hillshade alpine repeats of 8 ms
  and 20 ms); quote the 4000² rows.

### Server campaign 2026-10-12 (forth01, release 0.6.0): rerun at the tag

The same `server` profile and suite as 2026-10-11, rerun at `e09e809`, clean
tree. `git rev-parse "v0.6.0^{commit}"` resolves to `e09e809`, so this
campaign measures release 0.6.0. Raw data:
`../data/benchmarking/server_forth01.geos.ed.ac.uk_20261012_115629`. 372 rows,
all `ok`, the same row set as 2026-10-11.

- **Machine:** forth01, 8 logical CPUs in the session's mask (which 8 is not
  recorded in the manifest), `schedutil`, not exclusive. The CPU-use table
  flags nothing at 1 thread.
- **Code:** between `4d8d844` and `e09e809` the only changes are the version
  number and the default of `min_dist` in the positional overload of
  `breach_depressions_least_cost` (now `false`, matching the options struct
  and Whitebox). The benchmark scripts are unchanged.
- **Every engine took less time.** New / old median at 1 thread, median over
  rows: lsdtt3 0.87, RichDEM 0.86, TopoToolbox 0.89, Whitebox 0.86. RichDEM,
  TopoToolbox and Whitebox are the same builds as before, so the shift is in
  the session, not the code. Absolute seconds from the two campaigns should not
  be mixed. Whitebox's fills and `route_d8_epsilon` moved more than the rest
  (0.56-0.79 at 2000² and 4000²), so their ratios with lsdtt3 changed (table
  below).

#### 1-thread comparison

Reference time / lsdtt3 time, range over the six fixtures (above 1 means lsdtt3
took less time). The 2026-10-11 range is in brackets where it moved by more
than 0.1 at either end.

| case | RichDEM | TopoToolbox | Whitebox |
|---|---|---|---|
| `accum_d8` (S-vector on) | 0.98-1.12 | 1.33-1.72 (1.21-1.72) | 1.01-1.16 |
| `accum_d8_no_svector` | 1.39-1.64 | 1.92-2.43 (1.77-2.26) | 1.42-1.83 |
| `fill_epsilon` | 0.96-1.07 | n/a | 0.92-1.93 (1.03-2.75) |
| `fill_priority_flood` | 0.52-0.80 | n/a | 0.90-1.74 (1.05-2.10) |
| `fill` (reconstruction) | n/a | 0.85-1.00 | n/a |
| `breach` (Lindsay 2016) | 0.93-1.07 (0.94-1.49) | n/a | n/a |
| `breach_least_cost` (bitwise-equal output) | n/a | n/a | 0.60-0.79; 1.36 on alpine 4000² |
| `route_d8_epsilon` | n/a | n/a | 0.89-1.34 (0.99-2.00) |
| `accum_dinf` (report-only) | 0.95-1.16 | n/a | 0.74-0.93 |
| `accum_mfd` (report-only) | 0.60-0.64 | n/a | 0.55-0.65 |
| `slope_only` (report-only) | 0.46-0.51 | 0.96-1.05 | 2.35-2.63 |
| `hillshade` (report-only) | n/a | 0.47-0.56 | 9.8-12.0 |
| `polyfit_all` (report-only) | 0.21-0.27 | 1.03-1.15 | 0.96-1.13 |

The `breach` upper end of 1.49 on 2026-10-11 was RichDEM alpine 4000² at
5.17 s; here it took 3.26 s.

Raw DEM to D8 accumulation with each package's defaults (`route_d8_default`,
report-only), seconds:

| package | alpine 1000² | lakeland 1000² | alpine 2000² | lakeland 2000² | alpine 4000² | lakeland 4000² |
|---|---|---|---|---|---|---|
| lsdtt3 | 0.213 | 0.216 | 0.949 | 0.945 | 4.227 | 4.464 |
| RichDEM | 0.213 | 0.208 | 0.972 | 0.924 | 4.278 | 4.329 |
| TopoToolbox | 0.133 | 0.257 | 0.513 | 1.149 | 2.369 | 5.640 |
| Whitebox | 0.103 | 0.411 | 0.379 | 4.963 | 1.980 | 60.053 |

#### Least-cost breaching

Matched parameters as on 2026-10-11. Seconds:

| | alpine 1000² | lakeland 1000² | alpine 2000² | lakeland 2000² | alpine 4000² | lakeland 4000² |
|---|---|---|---|---|---|---|
| lsdtt3 | 0.035 | 0.234 | 0.118 | 0.930 | 0.434 | 4.069 |
| Whitebox | 0.021 | 0.163 | 0.071 | 0.640 | 0.590 | 3.210 |
| lsdtt3 / Whitebox | 1.68 | 1.43 | 1.67 | 1.45 | 0.74 | 1.27 |
| lsdtt3, CLI defaults (`breach_least_cost_default`) | 0.031 | 0.176 | 0.115 | 0.813 | 0.425 | 8.189 |

Against 2026-10-11, Whitebox's times fell to 0.74-0.82 and lsdtt3's to
0.87-0.89, so lsdtt3 / Whitebox rose from 1.20-1.45 to 1.27-1.68 on the five
fixtures where Whitebox is faster. Peak memory at 4000² is unchanged: lsdtt3
644-657 MiB, Whitebox 1757-1774 MiB.

#### Thread scaling on 4 physical cores

Speedup over each engine's own 1-thread time, 4000², alpine / lakeland:

| case | lsdtt3 2t | lsdtt3 4t | Whitebox 2t | Whitebox 4t |
|---|---|---|---|---|
| `slope_only` | 1.50 / 1.51× | 2.59 / 2.61× | 1.81 / 1.81× | 3.02 / 2.82× |
| `polyfit_all` | 1.71 / 1.71× | 2.96 / 2.96× | 1.51 / 1.53× | 1.98 / 1.95× |
| `hillshade` | 1.56 / 1.54× | 2.04 / 1.56× | 1.66 / 1.62× | 2.45 / 2.42× |

- Whitebox now gains from a second core: 1.51-1.81× at 2 threads, with
  1.39-1.57 CPU/wall, against 1.00-1.03× and 0.99-1.08 CPU/wall on
  2026-10-11. The suite code is the same, and `whitebox_cpu_pin()` picks the
  CPUs from the session's affinity mask, which is the one input that differs
  between sessions and was not recorded.
- lsdtt3's 4-thread speedups rose (slope 2.14-2.21 → 2.59-2.61×,
  `polyfit_all` 2.12-2.13 → 2.96×) and its 2-thread slope speedup fell
  (1.80-1.82 → 1.50-1.51×). The 2-thread slope rows themselves are about the
  same in seconds (0.51 s then, 0.53 s now) while the 1-thread rows
  fell (0.93 → 0.79 s), so the change is in the 1-thread baseline.
- 4-thread lsdtt3 rows still have outlier repeats: `polyfit_all` 4000² has
  one repeat of 1.02 s among four of 0.69-0.76 s; hillshade alpine 4000² has
  one of 0.100 s among four of 0.075-0.078 s, and lakeland did not speed up
  past 2 threads. Hillshade at 4 threads is not quotable from one campaign.

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
for ¹⁰Be/²⁶Al/¹⁴C under BraucherBorchers; the CAIRN Table 2 ²⁶Al F values sum
to exactly 1 in code, because lsdtt3 uses 0.0275 where the published table
prints 0.00275 (which would sum to 0.97525) — see the do-not-"fix" note in
`tests/cpp/test_cosmo_physics.cpp`; ¹⁴C uses P0 = 15.533, F = {0.788, 0.212,
0, 0} (Lupker 2015) and λ = 121.58e-6 yr⁻¹ under BraucherBorchers (5700 yr
half-life, IAEA 2024); λ(¹⁰Be) = 500e-9 and λ(²⁶Al) = 980e-9
yr⁻¹. Granger is frozen as exact in-code constants with a comment
naming its literature lineage (COSMOCALC). `Schaller` is NOT an unanchored
scheme: its ¹⁰Be production combines Braucher et al. (2011) and Borchers et al.
(2016), transcribed by Mirjam Schaller, so the scheme name identifies the
transcriber rather than a source. An earlier comment citing "Schaller et al.
2009/2016" for an updated Gamma[1] was wrong on both counts — no such
publication was found, and Gamma[1] = 1500 g cm⁻² is shared with Braucher and
BraucherBorchers rather than being specific to this scheme.

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

The Muscheler-VDM + ERA40 machinery (`crn_nc_fname` + `sample_age_yr`) is
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

### 7.10 Landslides and the Niemi (2005) ensemble

`unit_niemi_ensemble` (`tests/cpp/test_niemi_ensemble.cpp`, 24 checks, ~6 s)
reproduces Niemi et al. (2005) on a synthetic 110×90 fixture at 30 m, over 200
landslide seeds and 21 basins spanning 1.5 decades of area. It emits
`niemi_basins.csv`, `niemi_basins_volumetric.csv`, the matching area-bin files,
`niemi_ensemble.csv` and `niemi_timeseries.csv`. The acceptance script
`tests/scripts/niemi_ensemble_acceptance.py` reads those CSVs and evaluates the
paper's claims one at a time, so a failure names the claim rather than the run.

#### The denominator is the whole problem

**A landsliding landscape has three defensible "true" erosion rates and the
answer to "is the cosmogenic estimate biased?" depends entirely on which one you
divide by.** Getting this wrong cost this module a fictitious "unexplained
factor of 1.7", recorded as an open question for a fortnight.

| Rate | What it is | Spread across basins |
|---|---|---|
| imposed | the long-term rate the scenario was asked to deliver | none |
| volumetric, `E_v` | the mass leaving the basin **at the sampling instant** | large |
| apparent, `E_c` | the flux-weighted basin concentration, inverted as an analyst would | large |

Niemi's comparison is `E_c` against `E_v`, and **both are snapshots**. His
§2.4.2 defines the volumetric rate as "the total depth of material removed from
the watershed divided by the watershed area multiplied by the model **time
step**" — the same per-timestep depth array that weights his Eq. 6 mixture. His
Fig. 5 shows `E_v` with a *larger* spread than the cosmogenic rates at every
catchment scale, which a long-term average could not have.

Measured, 200 seeds:

| quantity | median |
|---|---|
| `E_c / E_imposed` | 0.46 |
| **`E_v / E_imposed`** — volumetric alone, no cosmogenics in it at all | **0.44** |
| **`E_c / E_v`** — Niemi's comparison | **0.96** (IQR 0.56–1.08) |

So essentially the whole apparent shortfall is in the denominator, and it is
Niemi's own §3.1 result: "the median 50% of the CRN-derived **and volumetric**
erosion rates generally fall below the total imposed erosion rate." The mass
arrives in rare large events, so the median snapshot carries less than its
share. The area trend is his too — `E_c/E_v` falls monotonically, 1.02 → 0.93 →
0.87 → 0.81, reproducing "at small catchment scales, CRN-derived rates are
typically higher than volumetric rates … but converge at larger catchment
scales."

The `E_c/E_v` median converges slowly (0.85 at 20 seeds, 0.94 at 100, 0.96 at
200, 0.97 at 400) while its quartiles barely move. That is the shape of the
distribution, not noise: a basin either has debris in transit at the sampling
instant or it does not, so the ratio is close to bimodal and the median sits in
the sparse trough between the modes. **The quartiles are the stronger
statistic**, and the assertion is that the IQR brackets unity.

#### The defect this caught, and the check that caught it

`trigger_landslide_event` accumulated the excavation raster (`exc = prev +
depth`) while `time_of_last_failure` recorded only the most recent event. Every
consumer reads the two as **one event** — the scar history strips the excavation
depth instantaneously at the failure time, and the flush weight delivers that
slab over the window that failure opened. A twice-failed pixel therefore
re-delivered its whole history's mass in the last flush window and was stripped
of the total-ever depth in one instant at the wrong time.

It was invisible to every existing test because the mass-balance check that did
exist is a **single-event** test, where cumulative and per-event are the same
number. What found it was a check nobody had written: the **domain mean of the
delivery field against what the series itself reports having excavated**, over
its own active area.

| | domain-mean delivery | expected |
|---|---|---|
| accumulating | 0.00852 g/cm²/yr | 0.00542 — **57 % excess** |
| per-event (fixed) | 0.00528 | 0.00542 — 2.7 % short |

The residual shortfall is expected and is asserted as a two-sided 10 % band:
single-failure memory legitimately loses the mass of repeat failures, and a
flushing pixel delivers its slab *instead of* its background erosion rather than
in addition to it. The excavation raster now replaces; nothing is lost, because
the cumulative record is the DEM.

**Generalisable lesson: a conservation check needs a regime where the quantity
can actually differ.** A mass balance asserted only on one event cannot see an
accumulation bug, and reads as coverage while providing none.

#### Two checks removed for being unfalsifiable

Both had passed since the ensemble was written, and neither could have failed
for the right reason:

- `|median(r) − 1| < median(|r − 1|)`. When every ratio falls below 1 — the
  normal case, and Niemi's own result — these are **algebraically identical**,
  since `|r−1| = 1−r` gives `median(|r−1|) = 1 − median(r)`. It passed by 0.001,
  on the handful of samples above 1.
- spread-of-medians `< 0.5 ×` `median(|r − 1|)`. That compares a **spread**
  against a **bias**: the right-hand side is dominated by the systematic offset
  from 1, so it passed by a factor of six while saying nothing about averaging.

The replacement compares spread against spread, and the honest answer is weaker
than the old one implied: the IQR of per-realisation medians is **0.77×** that
of single samples, not the ~6× the old check suggested. **The sampled basins are
nested** — a large one contains the small ones and all share one landslide field
— so they are nowhere near independent and the median cannot average down like
√n. Niemi's 46 sub-basins are nested the same way, so the caveat applies to his
practical advice as much as to ours: "sample several catchments" buys much less
than it looks unless the catchments are genuinely independent.

#### What is validated, and what is not

Validated against the paper: the `E_c/E_v` comparison and its area trend; the
volumetric rate's own under-read of the imposition; the depth-averaging factor
of his Eq. 5, pathway by pathway, to 7.5×10⁻¹⁵; and the landslide module's rate
closure at **1.86 %** against the ~3 % he reports. The **control** — no
landslides, apparent recovers imposed at every basin, worst |ratio − 1| = 0.0016
— is what makes any of the rest meaningful.

Not validated, and not to be quoted as if it were:

- **Nothing has been run on a real DEM.** Every fixture is synthetic.
- **Niemi Eq. 7 cannot be measured on this fixture.** `A_avg = 100 / E` predicts
  ~4700 km² at the fixture's rate, against a sampled range of 0.04–1.2 km², so
  the criterion never binds. The acceptance script reports it and says so rather
  than manufacturing a pass. Testing the *scaling* needs a DEM spanning more
  decades of area.
- **No realised landslide inventory has been compared with a mapped one.** The
  size distribution is imposed, so reproducing it proves the sampler works, not
  that the distribution suits any landscape.

#### Running it

```bash
ctest --test-dir build -R unit_niemi_ensemble --output-on-failure

cd build/tests/niemi_test_output
python3 ../../../tests/scripts/niemi_ensemble_acceptance.py \
    --basins niemi_basins.csv \
    --volumetric niemi_basins_volumetric.csv --strict
```

Like `landslide_distribution_acceptance.py`, the script is a diagnostic tool
rather than a CTest target: it is the thing to reach for when moving to a real
DEM, where the claims are no longer guaranteed to hold. `--strict` makes an
unmet claim exit non-zero; without it every claim is reported and the exit
status reflects only whether the files could be read.

---

### 7.10d The isotope-ratio output policy (`unit_cosmo_ratio`, 2026-08-05)

Ratios are deliberately **not** returned by the concentration kernels — that is
a settled convention, on the grounds that a ratio is a pure per-pixel division
and a nodata policy for zero denominators belongs at the output layer. What was
missing was *the* output layer: the policy lived as a lambda inside
`lsdtt-cosmo-prediction`, and the landslide coupling and basin-averaged paths
would each have had to re-decide it.

`make_isotope_ratio_raster` is now that layer, with the policy stated:

| case | result |
|---|---|
| both inputs valid, denominator ≠ 0 | the quotient |
| denominator exactly 0 | nodata, **counted** |
| one input masked, the other not | nodata, **counted** |
| both masked | outside the domain; not a failure |
| non-finite quotient | nodata, counted **separately** |

Two of those are corrections rather than transcriptions. The original lambda
tested only the **denominator's** mask, so a nodata numerator would have been
divided and the result emitted — harmless for `MultiNuclideRasters`, where all
three nuclides share a mask by construction, and wrong for any other caller.
And a non-finite result had no handling at all; it is unreachable given the
other rules, which is exactly why it is counted separately rather than lumped
in. If that count ever fires, something upstream produced an inf or NaN
concentration, and the alternative is a NaN propagating into a basin average.

**A zero numerator is a legitimate ratio of zero, not nodata.** Easy to get
wrong when writing the guard, and pinned.

`unit_cosmo_ratio` R8 requires the helper to reproduce the original lambda
**bitwise** on shared-mask input — not "close". A promotion that shifts a
shipped raster by one ulp is still a shipped raster that moved. The refactor
oracle (§7.11) was unchanged at 34/34, so the promotion is verified at both
levels.

#### A silent-overwrite bug found in the same code

`calculate_multi_isotope_step_change` and `calculate_multi_isotope_transient`
write **identical filenames** — both use `_step_crn_conc.tif`,
`_R2610_step.tif` and `_R1410_step.tif`. Nothing prevented both being enabled,
the transient ran first, and the step change then overwrote all five of its
outputs. The run reported success having produced half of what it was asked
for, under names that said "step" either way.

Now fatal, pinned by `failure_cosmo_prediction_both_multi_modes`. **Guarded
rather than renamed**: the `_step` suffix on the transient outputs is a naming
defect too, but changing shipped output filenames is a separate decision with a
migration cost, and the guard removes the data loss without imposing one.

---

### 7.10e Named unit conversions and the redundancy guard (2026-08-05)

`x * rho * 0.1` — a real thickness in metres to a mass depth in g/cm² — was
open-coded at **nine** sites across the library and two CLIs, each independently
choosing which density to use.

That is not a style complaint. The density split of 2026-08-04 exists precisely
because one of those choices was wrong: a soil scar carried at the bedrock
density (2650 rather than ~1500) is **1.77× too deep in mass**, which changes
spallation production at the scar base by **2.05×**. The arithmetic was never
wrong; the density was.

`src/cpu/lsdtt_crn_units.hpp` now owns it, and **the density argument has no
default, deliberately**. A defaulted density would restore exactly the failure
the file prevents: a call site that compiles, looks complete, and silently
assumes bedrock.

#### Expression order is preserved, not tidied

`0.1 * rho * m` and `m * rho * 0.1` are not the same double, and the refactor
oracle compares file hashes. Each function reproduces the multiplication order
of the sites it replaced, and `unit_cosmo_units` U7 asserts that **bitwise**
rather than "close". The oracle was unchanged at 34/34 across the migration.

`unit_cosmo_units` also recomputes the 1.77× and 2.05× figures from the
conversions rather than asserting them from memory, so the documented numbers
and the code cannot drift apart.

#### The guard, and what it caught

The pattern accumulated one call site at a time, so a one-off cleanup would
accumulate again. A **CMake-time configure error** now rejects any new
hand-written conversion in the cosmogenic and landslide sources — matching the
"unclassified test is a configure error" precedent, on the principle that a lint
nobody runs is a lint that does not exist. Verified to fire: reintroducing the
pattern in `lsdtt_landslide_size_cpu.cpp` fails the configure and names the
line.

Writing the guard immediately found a **ninth site the manual inventory had
missed** — `ice_effective_depth_g_cm2` in `lsdtt_quartz_fertility_cpu.cpp` —
and then a glob gap: the first draft scanned only `.cpp`, while two live
conversions sit in `.hpp` files. Both lessons are the same one: an inventory
assembled by reading is not a guarantee, and only the automated check is.

**Exempt, deliberately**: `lsdtt_crn_column.{hpp,cpp}` and
`lsdtt_crn_particle.hpp`. `to_eff()`, `remove_mass()` and the four `depth_m`
recoveries are inside the locked numeric core and must not be touched.

#### One uniform-raster helper, with the difference made explicit

Four local variants of "a raster of one value shaped like this one" had drifted
into **two different behaviours**: `make_constant_raster` (cosmo-prediction) and
`like` (cosmo-landslides) filled the *whole array* including pixels outside the
mask, while `make_uniform` (cosmo-tool) wrote nodata outside the mask.
`get_data()` returns the raw array and the mask is not applied on write, so for
any raster that reaches disk — `_LS_EXCAVATION.tif` is built this way — those
are different files.

`make_uniform_raster(reference, value, MaskedFill)` replaces all three, and the
choice is a **parameter rather than a convention**. Standardising on either
would have been a silent behaviour change to a shipped output; the point of the
helper is to make a difference that was accidental into one that is stated.

The `uniform()` fixtures in the test files are deliberately **not** migrated:
they build from metadata alone with an all-valid mask, which is unambiguous, and
the guard covers `src/` and `cli/` where the semantics actually differed.

---

### 7.10c The `lsdtt-cosmo-landslides` rename, and the compartment trap in it (2026-08-05)

`lsdtt-landslides` became **`lsdtt-cosmo-landslides`**, with no compatibility
alias. Parameter files are unaffected; scripts invoking the old executable name
are not.

Recorded here for one non-obvious reason. The domain compartments (§2) derive a
test's domain **from its name, first match wins**, and the rule table has
`"cosmo"` *before* `"landslide"`. So renaming `smoke_landslides_*` to
`smoke_cosmo_landslides_*` would silently have moved all six CLI tests out of
the `landslide` compartment — `ctest -L landslide` would have stopped seeing the
CLI for the module it is named after, and nothing would have failed to say so.

The fix follows the existing `roering|implicit` precedent, which exists for the
same reason: a test can legitimately belong to two compartments.

```cmake
if(_test MATCHES "landslide" AND NOT _domains STREQUAL "landslide")
    list(APPEND _domains "landslide")
endif()
```

Verified: all six renamed tests now carry both `cosmo` and `landslide`, and
`ctest -L landslide -LE heavy` lists 14 tests as before.

**The generalisable point:** a name-derived classification scheme makes renames
load-bearing. Any rename that adds a word matching an earlier rule silently
reclassifies the test, and the "unclassified test is a configure error" guard
cannot catch it, because the test *is* classified — just not where it was.

**The refactor oracle (§7.11) was re-frozen** for this phase, and the diff is
the evidence it was behaviour-preserving: every `sha256`, `kind` and
`exit_code` identical, with only the informational `exe` column changing on the
four `L1` rows. A re-freeze whose diff touches a hash is a different kind of
change and needs a different justification.

---

### 7.10f The shared cosmogenic site component (2026-08-05)

`src/cpu/lsdtt_cosmo_site_cpu.{hpp,cpp}` now owns everything between "I have a
DEM and some parameters" and "I have a production raster and three shielding
rasters": loading `crn_parameters.nc`, resolving the atmosphere and scaling
scheme, building the pressure raster and the CRONUS emulator inputs, computing
or loading production, and resolving shielding from a raster or a scalar.

That stretch was written out separately in `lsdtt-cosmo-prediction` and
`lsdtt-cosmo-tool`, ~150 lines each. Both now call the component; neither
contains a single reference to `calculate_production_raster`,
`calculate_pressure_raster`, `load_nc_data` or `AtmSource` any more.

#### What "migrating lsdtt-cosmo-tool" means — and what it does not

`lsdtt-cosmo-tool` **inverts measured concentrations into erosion rates**, so it
is worth being exact about what was touched.

Its inversion call is one line:

```cpp
double erate = lsdtt::calculate_eff_erate_from_conc(
    conc, nuclide, muon_scheme,
    *prod_for_inversion,                                        // production
    topo_shield_raster, snow_shield_raster, self_shield_raster, // shielding
    flow, outlet_node, predicted_conc, emu_inputs);             // emulator inputs
```

Five of those arguments are **fields the tool must build before it can invert
anything**. Building them is what moved. Nothing else did:

| moved into the component | untouched |
|---|---|
| load `crn_parameters.nc` | reading the sample CSV |
| resolve atmosphere / scaling scheme | snapping samples to the channel network |
| check the LSDn grids are present | basin accumulation |
| build the pressure raster | `calculate_eff_erate_from_conc` itself |
| compute the production raster | the self-consistent-age iteration |
| build the three shielding rasters | uncertainty propagation, output CSV |

**Why an inversion tool needs the same setup as a forward tool** — and therefore
why the duplication existed at all: an inversion is *"find the erosion rate
whose forward prediction matches the measurement"*. `calculate_eff_erate_from_conc`
is Newton–Raphson plus bisection over the forward model, so it needs exactly the
production and shielding fields a forward run needs. The two CLIs had converged
on near-identical setup because the physics requires it.

**How that claim is held to account.** The refactor oracle (§7.11) hashes
`_CN_ERATE.csv` — *the inverted erosion rates themselves* — for four
`lsdtt-cosmo-tool` cases: `stone`, `lsdn`, non-unit shielding, and
`self_consistent_age`. It came back **34/34 bitwise**. That the oracle genuinely
reaches the inverted rates is not assumed: the 1-in-10⁶ production perturbation
run when the oracle was frozen moved `t1_CN_ERATE.csv` and `t2_CN_ERATE.csv`.
Both CLIs' `-h` parameter sets and **ordering** are also unchanged.

#### Capabilities, not a union — and every predicted difference was real

The two CLIs are **not** the same shape, and a naive merge would have been a
silent feature change. Each difference below was found by reading or by
migrating, and each is now an explicit capability flag:

| capability | prediction | tool | what it guards |
|---|---|---|---|
| `Production` | ✓ | ✓ | — |
| `ProductionFromDisk` | ✓ | ✗ | the tool has no `load_production_raster` |
| `ShieldingScalars` | ✓ | ✓ | — |
| `ShieldingRasters` | ✓ | ✗ | **the tool has no shielding-raster parameters at all** |
| `PressureAlways` | ✗ | ✓ | the tool reports basin-mean pressure in its CSV, so pressure is a *result* not an intermediate |

The last one also carries a policy difference: a pressure failure is **fatal**
where something depends on it, and a **warning → nan** where it is only
reported. Both behaviours were already in the tree.

Granting `ShieldingRasters` to `lsdtt-cosmo-tool` as a refactor side effect
would have added four parameters and an untested code path to an inversion tool
**while passing every numeric check**, because the prefixes default to `NULL`.
That is the specific accident the capability design exists to prevent.

#### Two design points worth keeping

- **`emulator_inputs()` is built on demand, not stored.**
  `CronusEmulatorInputs` holds a raw pointer to the pressure raster; a stored
  copy would dangle the moment the site were moved.
- **Production is re-computable, not just readable.**
  `self_consistent_age` rebuilds production over *one basin's pixels* at each
  iterate, and the original comment records that whole-DEM evaluation "was the
  dominant cost of this mode". `CosmoSite::recompute_production(dem, age,
  eval_only)` exposes the configured compute so that survives; a component
  owning only a finished raster would have forced the loop to bypass it and put
  the duplication straight back.

---

### 7.10h One landslide CRN kernel (2026-08-05)

`compute_landslide_surface_concentration` is now a **thin wrapper** over
`compute_landslide_step_change_concentration`. There was one erosion history in
the codebase; there are no longer two.

A steady landscape is the step-change kernel with equal rates, so the wrapper
expresses its scalar inputs as the uniform fields that kernel takes:

| scalar input | becomes |
|---|---|
| `production_scaling` | a uniform production raster |
| `total_shielding` | a uniform shielding raster — **kept separate**, since the column forms `S_tot = scaling × shielding` and solves the Vermeesch redistribution on the product |
| background rate | *both* `eff_erate_old` and `eff_erate_new` |
| — | time-since-change 10⁶ yr: steady for far longer than any scar age |
| — | flush window 0, so every scar takes the scar treatment |

`set_erosion_time_series` now has **one** call site in the file, down from two,
and the file is 55 lines shorter.

#### Collapsing two implementations destroys the test that compared them

T17 compared the two kernels across 384 combinations. Once one calls the other
that comparison is tautological — it cannot fail. If nothing replaced it, the
collapse could have moved every one of those 384 values undetected.

So T17 gained a **golden checksum**: an order-weighted sum over all 384
concentrations, asserted to 10⁻¹².

```
before the collapse   113969168222.3623
after                 113969169785.83105
```

The 1.4×10⁻⁸ shift is float32 raster storage — the scalar path used to compute
in double and store once, and now round-trips through the kernel's raster.
Nothing physical moved, and T17b independently checks the one value that would
have if it had.

**The refactor oracle was unchanged at 34/34**, which is worth understanding
rather than just noting: its `L1` case fires at `event_time_yr: 0` with
`crn_report_time_yr: 0`, so every scar sits at elapsed time zero and takes the
`steady_state` branch — identical through either route. The float32 shift lives
at *positive* scar ages, which `L1` has none of. The golden checksum is what
covers that range, and the oracle alone would not have.

**Generalisable point.** When two implementations are merged, the test that
proved they agreed becomes worthless at exactly the moment it is most needed.
Freeze the absolute values *before* merging, or the merge is unverifiable.

### 7.10f2 "The Niemi closure" — what that number actually is

Quoted throughout this document as **1.86437 %**, so it is worth saying plainly
what it measures. It is a **mass-balance check on the landslide generator**, and
it contains **no cosmogenics at all** — the test file has zero CRN references.

The loop it closes:

1. **Prescribe** a landslide erosion rate — say 50 mm/kyr — over the ACTIVE area
   (the hillslopes that can actually host a landslide, not the whole raster).
2. The module turns that into a volume budget per event, samples landslide areas
   from the size distribution, solves for how many landslides that implies,
   places them, and excavates them.
3. **Measure** the volume actually removed; divide by active area × duration;
   convert back to mm/kyr.
4. Does it come back as 50?

| requested (mm/kyr) | realised | ratio |
|---|---|---|
| 20 | 19.63 | 0.981 |
| 50 | 49.80 | 0.996 |
| 100 | 100.02 | 1.000 |
| 200 | 196.97 | 0.985 |

Worst departure **1.86 %**. Niemi et al. (2005) report **~3 %** for their own
model doing the same thing, which is the external comparison the name refers to.

**Why it is not zero.** The size distribution is heavy-tailed and the landslide
count is an integer, so one draw from the tail can move an event's whole budget.
Exact closure is impossible in principle; ~2 % is the achievable level.

**Why it is quoted to five decimal places.** It is a **sentinel, not a result
being re-established**. It exercises size sampling, placement and the volume
budget, and nothing cosmogenic — so it must not move when CRN code changes.
Every phase of the 2026-08-05 work was in CRN territory, and a shift in this
number would have meant the landslide generator had been perturbed by accident.
"About 1.9 %" would still read as fine if it had drifted to 1.94; `1.86437 %` is
a tripwire.

**Do not confuse it with the Niemi ENSEMBLE result** (§7.10, `E_c/E_v`), which
is a different check that unfortunately also carries his name. The closure is
volumetric bookkeeping in the landslide module; the ensemble is about what
landsliding does to a basin-averaged cosmogenic erosion rate.

**What it does not establish.** That the imposed size distribution suits any
real landscape, or anything about the cosmogenic coupling.

### 7.10g A landslide that just happened has happened (2026-08-05)

**A convention change, decided deliberately.** It is documented at length
because it changes `lsdtt-cosmo-prediction` output and because the convention it
replaces was written down in a header, so a future reader will find the old rule
and need to know it was superseded on purpose.

#### The question

*When a landslide fires at exactly the moment you are sampling, has it
happened?* Until 2026-08-05 the two tools answered differently:

| tool | answer | reported at a 2 m scar |
|---|---|---|
| `lsdtt-cosmo-landslides` | **yes** | 237,058 atoms/g — the depleted rock the scar exposed |
| `lsdtt-cosmo-prediction` | **no** | 1,295,640 atoms/g — background, as if nothing happened |

A factor of 5.5 on the same scar at the same instant.

The old rule was explicit, in `time_since_failure_from_time_of_failure`:
*"A failure recorded at or after `t_report_yr` has not happened yet at the
reporting time and is also returned as never."* It was arbitrary rather than
reasoned, and it has been replaced.

#### The rule now

**Elapsed time zero means the landslide just happened.** The surface is the rock
that was at the scar depth, carrying the concentration it had there — inheritance,
with no re-accumulation. Only a failure in the **future** of the reporting time
is "never".

Three sites changed, all in `lsdtt_landslide_crn_cpu.cpp`: the conversion helper
(`elapsed > 0` → `>= 0`), the raster kernel's failed test, and the flush
selector. The negative sentinel `k_landslide_never_failed` means `>= 0` still
excludes "never" correctly.

#### Why it needed its own branch, and a number that proves it

At `t_f == 0` the four-row erosion history would have a **zero-length
most-recent segment**, and that does not reduce to the inherited value: it gives
**362,047 atoms/g against the correct 237,058**, a factor of 1.53. So the
instant-of-failure case evaluates the pre-failure column at the scar depth
directly, via the same `steady_state()` call `compute_landslide_surface_concentration`
has always used.

Scaling and shielding are handed to that call as `(site, 1.0)`. That is exact,
not a shortcut: the column forms `S_tot = scaling × shielding` and solves the
Vermeesch redistribution on the product, so passing the already-multiplied value
is the same thing.

**One approximation, stated:** under a genuine step change this uses the new
rate's steady state rather than the mid-transient profile, so it is exact only
where the incision wave passed long ago. That is the same assumption the scalar
path has always made.

#### Why this was worth stopping for

`crn_report_time_yr` defaults to 0 and a `single_event` run defaults to
`event_time_yr: 0` — so **the default configuration of
`lsdtt-cosmo-landslides` sits exactly on this boundary.** Unifying the two
kernels without settling it would have silently turned *"the landslide exposed
depleted rock"* into *"the landslide did nothing"*, while every aggregate check
still passed, because the affected pixels would simply have read background.

The tool now says so at run time when the reporting time is at or before the
last event, rather than leaving the user to infer it.

#### The two kernels, and the test that pins them

There are still two landslide → concentration implementations, sharing no code
260 lines apart: `compute_landslide_surface_concentration` (scalar production
scaling and shielding) and `compute_landslide_step_change_concentration`
(per-pixel rasters). Nothing asserted they agreed — the density split of
2026-08-04 had to be applied to both by hand, and had it reached only one,
no test would have noticed.

`unit_landslide_crn` T17 now compares them across **384 combinations**: depth ×
age × production scaling × shielding × scar density × nuclide. They agree to
**1.2×10⁻⁷** — float32 raster storage, nothing more.

`total_shielding != 1` is in that grid deliberately: scaling and shielding reach
the column as a product, so a change that merged or reordered them would be
invisible at shielding = 1. This codebase has already paid 11.4 % for a
scaling/shielding mistake once.

**Generalisable point, and the reason this section exists.** The first version
of T17 used ages 50–50,000 yr and reported clean agreement. The disagreement
lived at `age = 0` — the boundary the grid excluded, and the shipped default. An
agreement test written only over the "interesting" range can certify an
equivalence that does not hold.

#### A gap this exposed in the refactor oracle

The oracle was **unchanged at 34/34** through this change, which is not the
reassurance it looks like: no oracle case exercises the *raster* landslide
coupling with a scar at elapsed time zero. The behaviour that moved is real and
simply is not covered. Worth a case when the landslide coupling next changes.

---

### 7.10b Validity is not capability (`lsdtt-cosmo-landslides`, 2026-08-05)

A shipped defect, found while planning the cosmogenic CLI rationalisation and
fixed as part of that rationalisation. It is recorded because the *shape* of it recurs.

**What it was.** `lsdtt-cosmo-landslides` caught any failure inside its cosmogenic
block, printed a stderr warning, and fell through to **exit 0** having written
no concentration raster. Measured on the shipped binary:

```
$ lsdtt-cosmo-landslides t.param          # muon_scheme: ThisIsNotARealScheme
EXIT=0 — "Warning: CRN concentration failed: unknown muon scheme ..."  (no raster)

$ lsdtt-cosmo-landslides t2.param         # muon_scheme: cronusEmulatorSt
EXIT=0 — "Warning: CRN concentration failed: ... cannot be set through
          reset_scaling() ..."                                        (no raster)
```

To a scripted workflow both are **successes that produced nothing** — a worse
failure mode than a crash, because nothing notices.

**The second case is the interesting one.** `cronusEmulatorSt` is a perfectly
valid lsdtt3 muon scheme; `is_valid_muon_scheme()` accepts it. But it is
pressure- *and* nuclide-dependent, and this tool reaches the column through
`reset_scaling()`, which deliberately throws for exactly those two schemes
because it is handed neither. **Validity is a property of the name; capability
is a property of the tool.** A single shared validator that checked only the
name would accept something that fails later — which is precisely what
`lsdtt-cosmo-column-tool` already documents and guards against with its own
`column_muon_scheme_is_supported()`.

That has a direct consequence for the planned shared cosmogenic site component:
**its validator must take a capability set, not just a scheme name.**

**The fix.** Validation up front and fatal, before the DEM is read or a single
landslide placed, plus the catch-all downgraded from a warning to a fatal. It is
**conditional on the cosmogenic output actually being requested**, because
`muon_scheme` and `nuclide` are inert for a pure-geomorphology run and
rejecting a stale line that changes nothing would be officious.

Three tests pin all three behaviours, and the third is the one that stops the
natural over-correction:

| test | pins |
|---|---|
| `failure_landslides_bad_muon_scheme` | unknown scheme ⇒ nonzero exit |
| `failure_landslides_emulator_scheme` | valid-but-unserviceable scheme ⇒ nonzero exit, message names the tool that does support it |
| `smoke_landslides_geomorph_only_ignores_muon_scheme` | a nonsense scheme with the cosmogenic output OFF still runs |

Both failure tests would have failed before the fix with "CLI was expected to
fail but returned 0", so they pin the change rather than merely accompanying it.

**The refactor oracle (§7.11) was unchanged by this**, at 34/34 — the proof that
a fix to invalid-input handling disturbed no valid run. That is the intended
division of labour: the oracle covers valid configurations, the `failure_` tier
covers invalid ones.

---

### 7.10i Goal 3 unlocked: multi-isotope step change + landslide coupling (2026-08-06)

Until 2026-08-06, `lsdtt-cosmo-prediction` rejected this combination with

```
FATAL: landslide_coupling is not yet available for
       calculate_multi_isotope_step_change.
```

That guard is now removed. The combined path is implemented by routing through
the existing multi-nuclide landslide kernel
`compute_landslide_step_change_concentration_multi`, so Be-10/Al-26/C-14 are
evolved once under one composed scenario and written consistently with the
existing ratio policy.

#### What changed in behaviour

- `calculate_multi_isotope_step_change` can now be used with
  `landslide_coupling: scar|flush|both`.
- One run writes the three concentration rasters and both ratio rasters for the
  composed step-change + landslide case.
- The sampled accumulated outputs continue to carry `R2610_ACC`/`R1410_ACC` as
  **ratios of accumulated means** (not means of per-pixel ratios).

#### Acceptance evidence

`regression_cosmo_multi_landslide_coupling` is the end-to-end guard for the
formerly fatal path. It:

1. generates landslide rasters with `lsdtt-cosmo-landslides`,
2. runs `lsdtt-cosmo-prediction` with
   `calculate_multi_isotope_step_change + landslide_coupling: both`,
3. asserts the three nuclide concentration rasters and two ratio rasters exist,
4. asserts sampled accumulated `R2610_ACC` and `R1410_ACC` equal
   `Al26_CONC/Be10_CONC` and `C14_CONC/Be10_CONC` within output precision.

Run used for sign-off:

```bash
ctest --test-dir build -R "regression_cosmo_multi_landslide_coupling|rigorous_cosmo_multi_accum_ratio" --output-on-failure
```

Both tests passed on 2026-08-06.

---

### 7.10j SC divide-step gate (A'6) reproduced (2026-08-13)

This is the acceptance measurement that gates the divide-equilibration work: the
known divide-cliff signature on the generated SC DEM must be reproducible with
the corrected divide diagnostics before anything built on top of them is
trusted. It is recorded here because it is the baseline every later divide
measurement is compared against.

Dataset used:

- DEM: `/mnt/c/LSDTopoTools/data/cosmo_transient_test/lsdtt3_version/mega_sc5x_conc_SC_DEM.tif`
- tags: `/mnt/c/LSDTopoTools/data/cosmo_transient_test/lsdtt3_version/mega_sc5x_conc_SC_TAG.tif`

Method used (same counting convention as the original transient report):

1. Divide-step baseline on directed 4-neighbour tag-4 vs tag-(1,2)
   adjacencies, measured as `z_tag4 - z_tag12`.
2. Cross-divide drainage fraction by D8 routing rebuilt on the SC DEM,
   counting tag-4 pixels whose D8 receiver is tag-1 or tag-2.

Measured values:

| metric | measured | reference (2026-08-12 transient report) |
|---|---:|---:|
| tag-4 vs tag-(1,2) adjacency count | 12,589 | 12,589 |
| median divide step (m) | +55.9469 | +55.95 |
| mean divide step (m) | +64.7419 | +64.74 |
| p10 / p90 (m) | +6.4191 / +140.8479 | +6.42 / +140.85 |
| fraction with affected side higher | 0.92398 | 0.924 |
| tag-4 -> tag-(1,2) crossings | 8,919 / 165,523 | 8,919 / 165,523 |
| cross-divide fraction | 0.0538838 | 0.0539 |

Verdict: **gate reproduced**. The corrected diagnostic path recovers the SC DEM
signature at the expected values, so A'6 is numerically satisfied and the
Phase C work can proceed on a validated baseline.

---

### 7.10k Transient (rarefaction) pixels and hillslope `S_c` (2026-08-14)

#### Naming: transient, not fan

The step-change tags for partially adjusted pixels were called `fluvial_fan` and
`hillslope_fan`. They are now `fluvial_transient` and `hillslope_transient`.

The zone is a **rarefaction fan** in the kinematic-wave sense — the spreading
solution, as opposed to the shock. It is emphatically **not an alluvial fan**,
and nothing about it is depositional, so in a geomorphology toolbox the old
names invited exactly the wrong reading of a tag raster. The mathematical
vocabulary survives only where it names the mathematical object itself:
`chi_fan_lower` / `chi_fan_upper` are the fan's edges in χ, and the `summary()`
line still reports "fan in chi".

#### When transient pixels exist at all

Only when `n != 1`, and then only for one direction of the uplift step.
Celerity goes as `ksn^(n-1)`, so:

| `n` | behaviour | stretches when |
|---|---|---|
| `n == 1` | celerity independent of `ksn`; the wave translates rigidly | never — no transient pixel is ever emitted |
| `n > 1` | steeper travels faster | the step **falls** |
| `n < 1` | steeper travels slower | the step **rises** |

The code decides on `c_new < c_old` rather than on `n` and the sign of `dU`
separately, so both directions fall out without a special case. `n == 1` is
excluded explicitly rather than left to that comparison: it is what guarantees
no transient tag at `n == 1`, and it keeps the fan algebra — which divides by
`(n - 1)` — off a singular branch.

#### Hillslope transience is inherited, and used to be flat

A hillslope is tagged transient because the channel it drains to is
mid-transition, not because anything was solved on the hillslope. The channel
side is resolved properly: a transient channel pixel takes its erosion rate from
its position within the fan, strictly between `U_old` and `U_new`.

That position used to be discarded at the hillslope boundary. Every transient
hillslope received the arithmetic midpoint `0.5*(S_c_old + S_c_new)` regardless
of how far through the transition its channel actually was, so a hillslope above
a 5%-adjusted channel and one above a 95%-adjusted channel got identical `S_c`.

Since 2026-08-14 `S_c` is interpolated on the pixel's **own erosion rate**:

```
f   = (e - U_old) / (U_new - U_old)
S_c = S_c_old + f * (S_c_new - S_c_old)
```

Erosion rate, not position in χ, is the interpolant: `S_c` is paired with the
rate the landscape carries (`S_c_old` with `U_old`, `S_c_new` with `U_new`).
Interpolating on χ would be a geometric proxy that does *not* reduce to the
endpoints correctly, because `e` is not linear in χ inside the fan
(`S ~ χ^(1/(n-1))`, hence `e ~ χ^(n/(n-1))`).

No new plumbing was required: pass 2a already inherits the receiving channel's
erosion rate onto the hillslope, which is exactly the channel that hillslope
hangs off. `f` is clamped to `[0,1]` because `e` is stored as float and the
endpoints must not overshoot; `dU` cannot be zero on this branch, since that is
the `no_step` case and it emits no transient pixels.

Measured on the `n = 0.7` fixture (47 transient hillslope pixels, `S_c_old`
0.25, `S_c_new` 0.40):

| quantity | before | after |
|---|---:|---:|
| realised `S_c` range | 0.325 only | 0.2504 – 0.3909 |
| largest departure from the midpoint rule | — | 3.17 m |

#### How it is tested, and what the test deliberately does not do

Restating the interpolation formula in the test would be near-tautological, so
the per-pixel assertion checks only the **geometry** — that the imposed
elevation really is `rec_elev + dx*S_c` for that pixel's own `f`, i.e. that the
per-cell `gradient_scale` reaches pass 2b intact and lands on the right pixel.

The physics is asserted separately and formula-independently, from `S_c`
recovered out of the elevations:

- **bounded** — every realised `S_c` lies between the two endpoints;
  interpolation, never extrapolation;
- **spans a range** — guards against collapse to *any* single constant, which
  every other check here would tolerate;
- **monotone in erosion rate**, checked pairwise across the whole population
  rather than at the extremes; the required sign follows the direction of the
  step, so a falling `S_c` is handled too (0 inversions measured);
- **departs from the midpoint rule by > 1 m**, so the interpolation cannot
  silently revert.

The previous guard — "every transient pixel is more than 1 m from `S_c_old`" —
was removed rather than repaired. Under interpolation, pixels near the leading
edge legitimately sit close to `S_c_old`; that is the point of the change.

#### The oracle is silent on this branch

`regression_cosmo_oracle` passes across this change, and that is **not**
evidence the change is correct. Its step-change parameter file does not set `n`,
so it takes the default `n = 1.0`, `is_shock` is forced, and zero transient
pixels are produced — the rasters are bit-identical because the altered branch
never executes. The oracle confirms nothing else was disturbed; it says nothing
about the transient treatment.

Closing that gap means an `n != 1` oracle tier. Until it exists,
`unit_step_change` is the only guard on transient hillslope `S_c`.

---

### 7.11 The refactor oracle (`regression_cosmo_oracle`)

`tests/data/cosmo_oracle.csv` is a **bitwise hash of every raster and CSV the
four cosmogenic CLIs produce**, across a matrix of scaling scheme, muon scheme,
shielding path and prediction mode. 34 artefacts over 15 cases; 8 kB tracked.

It exists because the production, pressure and shielding setup is CLI-local and
duplicated between `lsdtt-cosmo-prediction` and `lsdtt-cosmo-tool`, and is due
to be extracted into a shared component. **The oracle for a refactor is the
pre-refactor binary, and that binary stops existing the moment the refactor
lands** — so it is captured first, deliberately, as a tracked artefact.

Every subsequent change is one of exactly two kinds and must declare which:

| kind | obligation |
|---|---|
| behaviour-**preserving** | `regression_cosmo_oracle` stays green. That is the proof the refactor was safe. |
| behaviour-**changing** | re-freeze with `--freeze`, and state in the commit message which artefacts moved and why. |

**Do not re-freeze to make a red test pass.**

#### How the hash is exact

Rasters go through `gdal_translate -of ENVI` and the raw bytes are sha256'd —
bitwise on the pixel values, verified by round-tripping through numpy against
the source band. CSVs are hashed with absolute paths stripped so the run
location is irrelevant.

Two alternatives were tried and rejected:

- **`gdalinfo -checksum`** is a 16-bit value, and an all-zero band (a landslide
  excavation raster) checksums to 0 — it cannot tell "all zeros" from "all
  zeros in a different shape".
- **`gdalinfo -stats`** truncates `STATISTICS_MEAN` and `_MAXIMUM` to two
  decimal places, nowhere near enough for a bitwise claim, and it writes a
  `.aux.xml` sidecar that then caches stale statistics.

stdout is deliberately **not** hashed: it carries absolute paths and free text
and is not a stable contract. Everything it reports that matters is expressed in
the rasters. Exit codes **are** recorded, because a CLI that starts failing — or
starts succeeding — is a real change.

#### Verified to fail

A 1-part-in-10⁶ perturbation of the production raster moves **16 of the 22
cheap-tier artefacts**, including `lsdtt-cosmo-tool`'s inverted erosion-rate
CSVs. So the oracle reaches through the whole forward/inverse chain, not just
the raster it perturbs. The 6 that survive are the landslide and column-tool
cases, which take the scalar production path — itself a useful signal about
which cases exercise which path.

#### A branch the matrix does not reach: `n != 1`

The step-change case does not set `n`, so it runs at the default `n = 1.0`. That
forces the shock branch, and **no transient (rarefaction) pixel is ever
produced** in any frozen artefact. A change confined to transient pixels — the
hillslope `S_c` interpolation of 2026-08-14 is exactly that — leaves every hash
untouched, so the oracle goes green without having exercised the change at all.

This is the general trap: a green oracle proves only that the *covered* branches
did not move. Before reading it as a safety result, check that the matrix
actually reaches the branch you edited. See §7.10k.

#### A guard on the oracle itself

The parameter parser **ignores unknown names silently**, so a case with a typo
freezes the CLI's *default* rather than the behaviour it names — an oracle that
looks fine and guards nothing. This was not hypothetical: the first draft of
case C1 wrote `effective_erosion_rate` where `lsdtt-cosmo-column-tool` declares
`effective_erosion_rate_g_cm2_yr`, and froze a run at the default 0.1 g/cm²/yr
instead of the intended 0.002.

The script now validates every case's parameter names against the CLI's own
`-h` output — which is rendered from the declarations, so it is authoritative
rather than a second list to keep in step — and **fails the run** on an
undeclared name. An unparseable help block is also a hard error rather than an
empty "nothing to check" set; that caught a missing `re.MULTILINE` immediately.

#### Cross-checks that fall out of it

Freezing the matrix surfaced three consistency properties worth knowing:

- `lsdtt-cosmo-prediction` and `lsdtt-cosmo-tool` produce **byte-identical
  production rasters** under `stone`. That is direct evidence the shared-site
  extraction is legitimate rather than merely plausible.
- The multi-isotope Be-10 raster is **byte-identical** to the single-nuclide
  step-change raster, confirming at CLI level the equality the unit tier asserts
  internally.
- Changing shielding alone leaves the production raster identical and moves the
  concentration — so the two really are applied separately, which is the
  property an 11.4 % disagreement once turned on.

#### Running it

```bash
ctest --test-dir build -R regression_cosmo_oracle            # cheap, ~4 s
ctest --test-dir build -R regression_cosmo_oracle_heavy      # + LSDn/emulator

# re-freeze — a deliberate act
python3 tests/scripts/freeze_cosmo_oracle.py --build-dir build \
    --read-path tests/data --oracle tests/data/cosmo_oracle.csv \
    --work-dir /tmp/oracle --tier all --freeze
```

The heavy tier carries the `lsdn` and CRONUS-emulator rows: they need the LSDn
scaling grids in `crn_parameters.nc` and carry the most setup state, which is
exactly why they are worth keeping rather than dropping.

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
`lsdtt-synthetic-surfaces`) returns `{dem, channel_mask}` with flat, analytic
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
Run it with `ctest -L '^performance$' -R hillslope` after configuring with
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
`lsdtt-synthetic-surfaces`, and `smoke_lem_roering_implicit_strip` runs
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

## 8b. Synthetic networks, routing fidelity, and divide equilibration

`lsdtt-synthetic-networks` generates an Optimal Channel Network (OCN) by
simulated annealing and then paints elevations on it as a chi-based fluvial
steady state. Everything downstream — the LEM, the cosmogenic scenarios, any
basin analysis — throws the graph away and re-routes the DEM. So the question
that decides whether the generator is usable at all is: **does the DEM it writes
still carry the network it was asked to generate?**

The measurements below answer that, and then record what the in-cycle hillslope
solve does about the answer. All of them were made on 2026-08-13.

### 8b.1 D8 routing fidelity is bounded well below 1 — and that is structural

Fidelity here is the fraction of interior, non-outlet land pixels whose D8
receiver on the synthesised DEM equals their OCN receiver. It is computed by the
CLI itself (`verify_routing`, reported as `routing_fidelity` in
`<write_fname>_ocn_summary.csv`) and asserted in-process by
`unit_synthetic_networks` (case `E`), which regenerates the OCN, synthesises the
DEM, routes it with D8 and compares receiver for receiver.

Measured on a 48×48 grid, cell 50 m, γ = 0.5, seed 42, `single_south` outlet,
2116 interior pixels:

| `n_iter` | Tscale = 0 | 0.01 | 0.1 | 0.3 | 1 |
|---|---:|---:|---:|---:|---:|
| 5 000 | 0.4660 | 0.3374 | 0.3077 | 0.3119 | 0.3053 |
| 20 000 | 0.7032 | 0.5723 | 0.4693 | 0.4031 | 0.3620 |
| **92 160 (production default, 40·N)** | **0.916352** | 0.879962 | 0.856805 | 0.873819 | **0.845936** |
| 400 000 | 0.9353 | 0.9069 | 0.9008 | 0.8993 | 0.8974 |

On the production default schedule the endpoints are **0.916352 (1939/2116) at
`temperature_scale = 0`** against **0.845936 (1790/2116) at
`temperature_scale = 1`**, same seed and same grid.

Three things follow, and the first is the load-bearing one.

**1. Fidelity plateaus near 0.935 and never reaches 1, however long the anneal
runs.** Even at 400 000 iterations — more than four times the production
schedule — roughly 6–8 % of interior pixels re-route under independent D8. This
is not a defect in the annealer: a local minimum of `sum(A^γ)` under single-flip
moves is simply not the same object as the steepest-descent graph of the chi
surface built from it. The consequence for the whole synthetic-landscape chain
is that **the OCN chi DEM does not exactly carry the graph that generated it**,
so re-routing always perturbs the network — most visibly at divides, where the
chi surface has one-pixel discontinuities by construction (fluvial relief only,
no hillslope component). That makes divide equilibration a **structural
requirement of the generator, not cosmetic polish**: there is no anneal length
at which it becomes unnecessary. Anything that reads a generated DEM and cares
about which basin a pixel belongs to has to deal with this.

**2. Anneal length dominates temperature.** Truncating `n_iter` from 92 160 to
20 000 costs more fidelity (0.916 → 0.703 at `Tscale = 0`) than the entire
temperature range costs at fixed `n_iter` (0.916 → 0.846). A branching study
that varies `temperature_scale` must therefore hold `n_iter` at the default or
scale it with N, or it confounds "more side branching" with "under-relaxed". At
`n_iter = 5000` with `Tscale ≥ 0.1` the anneal ends *above* its starting energy
(`E_final/E_0` = 1.06–1.11) — it never cools at all. Put `n_iter` and
`E_final/E_0` in any screening table so that failure is visible.

**3. Fidelity is not monotone in the interior of the temperature range** (0.8568
at `Tscale = 0.1` against 0.8738 at 0.3). Only the endpoint inequality is safe
to assert, and that is all the test asserts.

The gate is therefore a **floor pinned at 0.90**, about 1.8 % under the measured
0.916352, plus the endpoint inequality `fidelity(Tscale = 1) < fidelity(Tscale = 0)`.
Negative controls were verified: truncating `n_iter` to 20 000 and flattening the
hot run's temperature each fail the relevant assertion and nothing else.

A weaker property is asserted unconditionally and should not be mistaken for
this one: because `z[i] = baselevel + f·chi[i]` and chi strictly increases
upstream, **every land pixel is strictly higher than its OCN receiver at any
temperature**, so the OCN DEM has no interior local minimum. That is necessary
but nowhere near sufficient for fidelity — a pixel can be higher than its OCN
receiver while a *different* neighbour is lower still, which is exactly the
failure the D8 round trip catches and monotonicity does not.

### 8b.1b The energy is length-weighted by default (2026-08-14)

`lsdtt-synthetic-networks` minimises

```
E = sum_i A_i^gamma * L_i
```

where `L_i` is the length of the link from node i to its receiver: `cell_size`
for a cardinal step, `cell_size*sqrt(2)` for a diagonal one. **This is the
default.** `energy_length_weighted: false` drops `L_i` and reproduces OCNet
(Carraro et al. 2020) bit-for-bit.

#### Why the length term belongs there

The energy expenditure of a channel link is `E_i ~ L_i Q_i S_i` — Strong & Mudd
(2022, *PNAS*, "Explaining the climate sensitivity of junction geometry in
global river networks", eq. 4), who note the problem is "length minimization
where the junction link lengths are weighted according to `A_i^gamma`". Goren &
Shelef (2024, *Earth Surf. Dynam.* 12, 1347–1369,
doi:10.5194/esurf-12-1347-2024, eq. 5) use the same form with `l_i` explicit,
and observe that with uniform `l_i` the energy reduces to a function of the
number of nodes.

On a **4-connected** lattice every `L_i` is equal, so it factors out and the
literature abbreviates it to `sum(A^gamma)`. On the **8-connected** lattice used
here it does not factor out: a diagonal link spans 1.414 cell widths for the
same energy, so lateral movement is free.

#### What that cost, measured

400x200 strip at 30 m draining north and south, periodic sides, gamma 0.5,
temperature_scale 0, fixed `n_iter`. Flow directions from
`lsdtt-channel-extraction`, **area-weighted over channels > 2 km²** — i.e.
where the water actually goes:

| variant | diagonal | pure N/S |
|---|---:|---:|
| unweighted (OCNet parity) | **0.937** | 0.062 |
| length-weighted, BFS seed | 0.827 | 0.151 |
| length-weighted + Dijkstra seed (default) | **0.056** | **0.944** |

The main channels ran diagonally in 94% of the flow, on a strip whose only base
levels are north and south. Weight by area or the effect hides: the *count* of
diagonal links only moves 0.45 -> 0.22, because headwater links point in every
direction in any network. `unit_synthetic_networks` case L pins the
area-weighted statistic (0.937 -> 0.029 on its own smaller fixture).

#### The seed tree mattered more than the energy term

Length weighting alone moved the area-weighted diagonal fraction 0.937 -> 0.827
only. The rest came from the **initial spanning forest**. An unweighted
breadth-first search over eight neighbours is the **Chebyshev** metric, in which
a diagonal and a cardinal step are the same distance — the wrong metric for a
length-weighted energy, and it hands the annealer a network whose shortest paths
are already diagonal. Single-flip annealing at low temperature is greedy enough
not to undo it. Seeding instead with **Dijkstra on Euclidean link length** took
the fraction to 0.056, and also reached a *lower* final energy
(8.27e6 vs 8.62e6) and higher round-trip fidelity (0.843 vs 0.805). Both fixes
are selected together by `energy_length_weighted`.

#### Consequences to know about

- **Round-trip fidelity is LOWER under the default**, 0.837 vs 0.916 on the
  48x48 reference fixture, and the pinned floor moved from 0.90 to 0.82. This is
  not a regression: the default drives the graph toward straight cardinal trunks
  while the chi surface is still re-routed by D8 steepest descent, which favours
  diagonals, so the two disagree more. The round-trip gap is a property of the
  discretisation either way (see 8b.1).
- **The snap cycle converges sooner** — the 32x32 equilibration fixture reaches
  a period-1 fixed point at cycle 5 instead of 10, because the length-weighted
  network starts closer to self-consistent.
- **The A3 annealing check was weakened to its two-point form.** Under the
  length-weighted energy `E(0) < E(0.1) < E(1)` no longer holds
  (357686 / 361794 / 361234); the endpoints are still ordered and that is what
  is pinned. The test's own note anticipated this weakening.
- **Networks become quite regular.** The corrected strip is close to a comb of
  parallel N-S trunks. Whether that is under-branched for a real mountain belt
  is open, and is the next thing to look at rather than a settled result.

#### Both energies are pinned

`rigorous_synthetic_networks_oracle` now sets `energy_length_weighted: false`
and its frozen values are **unchanged** across this work — same drainage
checksum 3275, same CSV sha256. That unchanged hash is the evidence that parity
mode is bit-for-bit OCNet rather than merely close.
`rigorous_synthetic_networks_oracle_lengthweighted` pins the new default
(checksum 3905). `unit_synthetic_networks` case E3 additionally asserts that
parity mode reports the unweighted sum and that the two forms are
distinguishable, so the flag cannot silently become a no-op.

---

### 8b.1c Drainage-area exceedance exponent beta, and an unresolved gap (2026-08-14)

`compute_area_scaling` in `lsdtt_strahler_tokunaga_analysis_cpu` fits
`P(A > a) ~ a^-beta`, the standard characterisation of how drainage area is
partitioned across a network (Rodriguez-Iturbe et al. 1992 WRR 28:1095; Rigon et
al. 1996 WRR 32:3367). Exposed as `write_area_scaling_summary` /
`write_area_scaling_curve`; pinned by `unit_synthetic_networks` case B.

#### Three things that had to be got right, each of which was wrong first

**The fit must be two-sided.** The first implementation trimmed only from below
and fitted the tail, which is the standard heavy-tail recipe and completely wrong
here: at the large-a end sits the finite-size rollover, since no basin can exceed
the domain. Worse, for a complete single basin the endpoints PIN the full-range
slope near 1 regardless of the network, because N(>= smallest) = N and
N(>= A_max) = 1. Measured on a 400x200 OCN: full-range slope **-0.986**, while
half-decade local slopes in the scaling regime give beta = **0.31-0.53**. The
tail fit returned a confident **0.76**. `upper_exclude_decades` (default 1.0)
now discards the rollover before any cutoff is selected.

**KS selection alone favours short ranges.** A short segment of any smooth curve
looks like a power law, so the search happily returned r^2 = 0.999 over half a
decade. Without a floor the chosen range collapsed to 0.53-1.10 decades and beta
swung 0.29-0.52 with it. `min_fit_decades` (default 1.5) is enforced inside the
candidate loop, so a short-and-tidy segment is never a candidate; if the
constraint cannot be met the result is NaN, which is the correct answer when
there is not enough range.

**Grid size matters more than expected.** At 96x96 the fit is ill-conditioned
(r^2 = 0.83, beta = 0.57). The test fixture is 200x200; below that beta should
not be read.

The exponent is the **MLE** (Clauset, Shalizi & Newman 2009, SIAM Review 51:661),
not least squares on the log-log curve, which those authors show to be biased.
Both are reported: `beta` (MLE) and `beta_ols` with its `r_squared`, kept as
shape diagnostics. On the test fixture they agree to 0.02, and case B asserts
they agree within 0.1 -- a wide gap means the fitted range is not a power law and
the number should not be read as one.

#### The unresolved part, stated plainly

**lsdtt3 measures beta ~0.54-0.58 where the literature reports 0.43-0.45.**

The obvious suspect was the length-weighted energy that became the default on
2026-08-14. It is not: measured on the same grid, seed, gamma and schedule,

| network | beta |
|---|---:|
| OCNet parity (`energy_length_weighted: false`) | 0.540 |
| length-weighted (default) | 0.581 |

The parity network is bit-identical to OCNet and still sits ~0.1 above the
published value, so the gap is **not** caused by our energy change. It is either
in how beta is measured here versus in the literature, or in grid size, and
neither has been established. Candidates not yet tested: the literature may fit
over a different range, may use a different exceedance normalisation, or may
require substantially larger grids than 200x200.

Case B therefore asserts a **broad** band (0.3-0.8) and prints the value every
run, rather than a band tuned to the literature. It is a drift guard, not a
validation of agreement. **Do not cite beta from lsdtt3 as agreeing with
Rodriguez-Iturbe et al. until this gap is explained**, and do not narrow the band
to make a point -- the gap is the finding.

---

### 8b.2 Regression artefacts for the OCN — oracle, band, pin, and no sentinel

These four words are not interchangeable, and the OCN is a good example of why.

| Artefact | Kind | Why |
|---|---|---|
| GDAL checksum of `*_ocn_drainage_cells.tif` | **oracle** (refactor tripwire) | The area field is integer-valued and exactly representable in float32, so the checksum fingerprints the whole tree and is robust across toolchains. Primary. |
| sha256 of `*_ocn_summary.csv` | **oracle** | Energies print at six significant figures — robust to last-bit drift, sensitive to any real change. Secondary. |
| `--expect-raster-stat` min/max/mean on `*_ocn_dem.tif` | **band**, not an oracle | Float chi accumulation is `pow`-implementation sensitive; a hash here would be flaky across libm versions. |
| pinned `final_energy()` in the C++ driver, 1e-9 relative | **characterisation pin** | It records what the code does today so a refactor cannot move it silently. It is not a validated numeric standard. |

**There is no numerics sentinel for the OCN, and there must not be one.** A
sentinel asserts a *validated* numeric standard. No such standard exists here:
the OCNet R package uses a different RNG and a different initial tree, so it can
validate *statistics* (Horton and Tokunaga ratios, energy scaling) but never
values. Calling any of the above a sentinel would claim a guarantee that has
never been established.

### 8b.3 The elliptic hillslope solve in the snap cycle: contraction, not one-shot stabilisation

The design claim under test was this. Critical-slope hillslope construction is
path integration (`z_i = z_receiver + Sc·L_i`): the elevation depends on the
receiver and the receiver depends on the elevations, a fixed-point iteration on
a piecewise-constant, non-variational map with no Lyapunov functional forcing it
downhill — so period-2 orbits are its natural failure mode. A steady-state
elliptic solve is a different kind of object: with the channel network fixed as
a Dirichlet set, `D∇²z = −U` (or the Roering form, which stays elliptic and
monotone for `|∇z| < Sc`) has a unique solution, and the divide is wherever
`∇z = 0` in that solution — an *output* rather than an imposed integration path.
The prediction was therefore that re-routing on that surface reproduces the same
divide, and the second cycle is a no-op on hillslope pixels: **the solve
"stabilises on the first try"**.

The literal criteria were: **S1** — cycle-2 `frac_receiver_changed` on hillslope
pixels is exactly 0.0; **S2** — cycle-2 `divide.step_p90_m` within 1e-6 of cycle
1; **S3** — the limit-cycle detector reports `limit_cycle_period == 1` by cycle
3.

**Measured result: S1 and S2 FAIL. S3 passes.** The mechanism is right; the
literal claim is wrong. Measured on a 48×49 noisy cone with a nodata-ring
outlet, `'nnnn'` edges, D = 0.3, Sc = 0.8, tolerance 1e-4, `sor_omega = 1.3`:

| cycle | `frac_receiver_changed_hillslope` (linear) | `divide.step_p90_m` |
|---|---:|---:|
| 2 | 0.773429 (1502/1942) | 2.01537 |
| 3 | 0.0216272 (42/1942) | 2.01052 |
| 4 | 0.0015448 (3/1942) | 2.01069 |
| 5 | 0.0 (0/1942) | 2.01069 |

Roering behaves the same way (0.0293512 at cycle 2→3). So the elliptic route
does **not** stabilise in one cycle. What it does instead is contract by better
than 10× per cycle and reach a genuine period-1 fixed point at cycle 4 —
qualitatively unlike a limit cycle, which is the distinction that mattered.
**The operational expectation is 3–4 cycles, not 2.** That is a revision forced
by the measurement, not the original prediction. `hillslope_every_cycle = true`
remains the right default.

The committed tests assert what is *true* — >10× contraction per cycle, exactly
0.0 by cycle 5, period 1 — and **print** the literal S1/S2 values marked
`NOT MET` rather than quietly rewriting the criteria to match the outcome.

**What the residual instability is not.** Three candidate explanations were
measured and two were eliminated:

- **Not solver error.** Sweeping the hillslope tolerance across four decades
  (1e-2 → 1e-5) moves the cycle-2→3 number by less than 5 %: 0.0227 / 0.0216 /
  0.0216 / 0.0216.
- **Not a flickering Dirichlet set.** The channel count was constant at 38
  pixels across cycles, so the channel network is not moving underneath the
  solve.
- **It scales with diffusivity**: D = 0.3 → 0.0216, D = 1.0 → 0.0509,
  D = 3.0 → 0.157. The smoother the hillslope, the more pixels flip.

**Near-routing-ties are a contributing but incomplete cause.** For hillslope
pixels, the ratio of the second-steepest to the steepest D8 descent candidate
(1.0 being an exact tie):

| | median | p10 | p90 | within 1 % of a tie |
|---|---:|---:|---:|---:|
| flipped pixels | 0.917 | 0.736 | 0.9915 | 13.0 % |
| stable pixels | 0.717 | 0.425 | 0.921 | 2.1 % |

Flipped pixels sit far closer to a routing tie than stable ones — an 8 % margin
against 28 % — and the top decile of flips are effectively exact ties. But the
median flipped pixel still has a real 8 % margin, so the surface itself is
genuinely still moving a little between cycles. The consequence for anyone
tempted by the obvious fix: deterministic D8 tie-breaking in the router would
reduce the flip count and is worth doing on its own merits, but **it would not
by itself drive the hillslope receiver-change fraction to zero**.

### 8b.3b The snap cycle needs tens to hundreds of cycles, and does not approach smoothly (2026-08-14)

`steady_state_n_cycles` defaults to 1, the equilibration smoke fixture pins 5,
and the synthetic-network sweeps were first run at 8 and then 20. **Measured
cycles to convergence on a 400x200 strip, gamma 0.5, length-weighted OCN:**

| hillslope model | T = 0 | T = 0.5 | T = 1.0 |
|---|---:|---:|---:|
| linear | 68 | 28 | 30 |
| Roering | 73 | 36 | 31 |
| critical slope | never (500) | never (500) | never (500) |

Three things follow.

**Both elliptic models converge; critical slope does not, at any budget tried
(8, 20, 500).** Its median divide step sits at exactly 14.06 m throughout, and
29-50% of pixels change receiver on every cycle. Linear leads Roering slightly
because max hillslope gradients here are ~0.05 against Sc = 0.8, so the
nonlinear term is nearly inactive and Roering pays for Picard lagging it does
not need.

**The approach is NOT monotone.** On T = 0 linear, max|dz| per cycle runs 17.6
(cycle 10), 18.4 (20), 29.1 (30), **45.9 (40)**, 14.6 (50), 6.7 (60), 13.9 (65),
then **0.10 at cycle 66** and converged at 68. The excursion at cycle 40 is
larger than the value at cycle 10, and the collapse happens inside one cycle.
**A truncated run therefore lands somewhere on the wandering path, not near the
answer, and its residual gives no indication of how far from the end it is.**
Simon reports the same behaviour in lsdtopotools2, sometimes needing ~1000
cycles to reach oscillating divides.

Do not read "reached max_cycles" as "the landscape will not settle" without
having given it hundreds of cycles. Two earlier reports in this project made
exactly that error at 8 and at 20 cycles.

**A more aligned network is harder to equilibrate.** T = 0 needs 68-73 cycles
against 28-36 at T = 0.5-1.0, and it is the greedy case that produces the
strongly axially aligned comb (R2 = 0.879 against ~0.32; see 8b.1b). Trunks that
are near-parallel and nearly equal in elevation leave the divide between them
weakly constrained, so the set of flipping pixels empties slowly. The same
comparison at fixed settings across the ENERGY change: the unweighted OCNet
network converges at cycle 17, the length-weighted one at 68 -- so the
better-oriented default costs roughly 4x in equilibration cycles.

#### Why it wanders: the cycle is a closed feedback

`lsdtt_landscape_equilibration_cpu.cpp:820-1226`. Each cycle fills last cycle's
**delivered** surface (hillslopes included), re-routes it, **recomputes drainage
area**, **re-snaps channel elevations** by chi integration on those areas,
rebuilds the Dirichlet channel set with hysteresis, solves hillslopes against
the new channel elevations, and then sets `cycle_dem = delivered` (`:1226`).

That last assignment is what closes the loop -- were it `snap.elevation`,
hillslopes could never move a divide. So one pixel flip perturbs downstream
area, which shifts channel elevation, which moves the hillslope solution and can
flip further pixels. The wandering is that feedback, not solver noise.

---

### 8b.4 What the hillslope solve buys: the divide step

This is the practical payoff and the real deliverable. Same fixture, hillslope
solve on against off:

| metric | solve off | solve on | change |
|---|---:|---:|---|
| divide `step_p90_m` | 8.32 m | 3.27 m | 2.5× smaller |
| divide gradient, median | 0.0961 | 0.0235 | 4× smaller |

Read this against §8b.1: the generated DEM never exactly carries its own
network, and the divide is where that shows up as a cliff. The elliptic solve
does not remove the effect, but it reduces the divide step by a factor of 2.5
and the divide gradient by a factor of 4 — which is why the solve is on by
default in the equilibration cycle even though it costs an extra solve per
cycle.

Note also that convergence and acceptability are separate questions. A landscape
can reach a period-1 fixed point that still has a cliff across the divide — the
OCN chi surface with no hillslope treatment does exactly that.
`divide_gradient_target` reports the acceptance criterion ("is the landscape
smooth across divides") as distinct from "has the iteration stopped moving", and
it is reported rather than enforced by default.

### 8b.5 Two solver traps, recorded because they cost real time

Both were found while wiring the hillslope solve into the equilibration cycle,
and both are warned about rather than fixed — the hillslope steady-state solver
is a protected numeric core.

**(i) The nonlinear (Roering) SOR solve diverges under auto-omega.** With
`sor_omega = 0` — the default, meaning "auto-estimate from grid size" — and
`nonlinear = true`, the solve ran to the 50 000-iteration cap at a residual of
~480 m and produced ~145 m of divide-step nonsense, **at every diffusivity and
tolerance tried** on the 48×49 cone fixture. `sor_omega = 1.3` converges in
~1000 iterations. The plausible mechanism is that the auto-omega formula is the
optimal-SOR result for the *linear* Poisson problem with Dirichlet boundaries on
a rectangle, and has no derivation behind it for the Picard-lagged nonlinear
operator, whose effective coefficient `1/(1 − (|∇z|/Sc)²)` changes between
sweeps. Nothing currently caps omega when `nonlinear = true`. The equilibration
cycle warns on stderr and writes `DID NOT CONVERGE` into its `stop_reason`.

> If a solve is used as the *reference* for another solve, assert that the
> reference converged before comparing fields. A reference that silently hit its
> iteration cap makes the comparison pass or fail for the wrong reason.

**(ii) An absolute tolerance that cannot be reached.** The residual floors near
**2e-4 m on a ~200 m elevation range** — consistent with float32 round-off. The
solver's `tolerance` is an *absolute* `max|Δz|` criterion, so
`hillslope_ss_tolerance: 1e-5` cannot succeed on a landscape of that relief and
simply burns `max_iterations`. The symptom looks like "the solver is slow"; the
cause is "the criterion is unreachable".

---


## 8c. Spatially varying erodibility in the LEM (2026-10-16)

Per-pixel erodibility K(x) in the stream-power solvers and the analytical steady
state, the `lsdtt-lem` model-grid and input-raster placements, and the Python
`LandscapeModel` field. Tests:

| Test | Labels | What it covers |
|---|---|---|
| `unit_lem_erodibility` | `unit;lem` | kernels (K1-K6), units and conversion (U1-U3), calibration (V11-V12), steady state (S1-S6), values and scaling (V1-V14, F1-F3), `run_lem_step` (L1-L4), model grids and placements (G1-G16). Driver `tests/cpp/test_lem_erodibility.cpp`, `--suite fast` |
| `rigorous_lem_erodibility_analytic` | `rigorous;lem` | analytic steady states A1, A1b, A2, A3 (`--suite analytic`) |
| `unit_raster_loader_geotransform_guard` | `unit;gdal` | the loader refuses rotated, south-up and projected non-square rasters (G1-G8) |
| `unit_lem_forcing_history` (F5) | `unit;lem` | an explicit forcing K <= 0 is refused, naming the row |
| `smoke_/failure_lem_erodibility_*`, `smoke_/failure_lem_grid_*`, `..._channel_mask_*`, `regression_lem_*` | `smoke`/`failure`/`regression`, `cli;lem` | 39 CLI runs on fixtures from `fixture_lem_erodibility_cli`, and native comparisons of their outputs |
| `smoke_lem_synthetic_default_georeference`, `smoke_lem_synthetic_supplied_crs` + `regression_lem_synthetic_*_crs` | `smoke`/`regression`, `lem` | the default synthetic georeference (Watopia, EPSG:32758, grid centre within one cell) and a supplied CRS |
| `test_landscape_model.py::LandscapeModelErodibilityTest` | pytest | the Python field: bitwise equality with the scalar, `set_erodibility_field`, refusals, the strip identity |

### How correctness is established

* **Bitwise, within one build.** A uniform field equal to the scalar K must give
  the scalar result bit for bit: one `run_lem_step`, five steps (D8 and D-inf,
  n = 1 and 2), the steady state (D8 and D-inf, K and relief mode), the CLI snap
  and the CLI run (compared with `lsdtt-test-lem-erodibility --compare-rasters`).
  No test pins an absolute checksum of evolved topography.
* **The pixel's own K, not its neighbour's.** K3 recomputes the n = 1 closed form
  with each pixel's K to 1 ulp. K5 and K6 evaluate the Newton residual with the
  pixel's K (within `slope_clamp · dx + 2 ulp`) and with the receiver's K, which
  must exceed 10× the tolerance on at least 90 % of crossing pixels that incise
  at least 0.01 m (100 % measured). These are the tests that catch K applied on
  the wrong side of the implicit update; the steady-state identities do not.
* **Analytic steady states.** A1: on an 80 × 3 strip (m = 0) the steady slopes
  either side of a vertical contact are U/K_i: S1 = 0.1000 and S2 = 0.39999
  (expected 0.1 and 0.4) at n = 1. A1b and A2: K_i a^m S^n = U on every pixel to
  1e-3/n (measured worst 9.5e-4 at n = 1, 4.9e-4 at n = 2), for D8 and D-inf.
  S2: the K-weighted chi on a 200-node chain matches a double reference within
  200 ulp and U = K_i a^m S^n link by link to ~5e-6.
* **Placements, value by value.** Source rasters hold `1 + 1000 r + c` (× a
  constant), so every clip and warp is checked pixel by pixel against indices
  computed in the test, including points on grid lines.

### Findings

* **Steady identities must be measured on the surface the step solved.** Boundary
  pixels that are not flow outlets erode inside the implicit sweep and are reset
  only after it, so their donors are solved against the eroded value. Measured on
  the uplifted, restored surface, the identity is off by up to a factor of 4 next
  to the boundary. The analytic tests evaluate it on one kernel call's output.
* **p999 is not a per-pixel steady-state criterion.** Stopping at p999 |dz| <
  1e-3 U dt leaves 0.1 % of pixels moving, enough for a per-pixel identity to fail
  at 2e-3. The analytic tests stop on the maximum.
* **The D-infinity two-receiver Newton stalls when a pixel ends below its cardinal
  receiver** (f' changes sign); 33-40 of ~1100 two-receiver pixels end there and 6
  had residuals of 0.1-1.4 m, with scalar K as much as with a field. K6 excludes
  and counts them.
* **The analytical snap is not a fixed point of the transient model** (see the
  landscape evolution design notes): A3 compares the first step after a field snap
  with the first after a scalar snap (99.9th-percentile change 1.368 m and 1.365 m)
  rather than demanding a steady state the scalar snap does not reach either.
* **The GDAL warp extent lost precision.** G6 (a 4 m raster clipped onto a 10 m
  grid) failed by exactly one row: the warp passed `-te` with 6 significant digits,
  so a northing of 5000202 became 5000200. Fixed for every warp and clip.
* **The geotransform guard first broke reprojection.** The only test raster with
  non-square pixels is a geographic DEM (1.5" × 1"), which
  `lsdtt-raster-preprocessing` must open to reproject. The guard therefore exempts
  geographic rasters from the square-pixel rule.


### 8c.1 Lithology maps to erodibility (`lsdtt-lithology-erodibility`)

Every test whose name contains `lithology` also carries the label `lithology`, beside its
domain (`lem`, `fetch`, later `chi` and `network`), so `ctest -L lithology -LE heavy` runs them all.

| Test | Labels | What it covers |
|---|---|---|
| `unit_lem_lithology_erodibility` | `unit;lem;lithology` | property table (valid, duplicates, bad numbers, K <= 0, missing m / n / convention, code range, extra columns); crosswalk; values to codes; the polygon rasteriser; the K raster; lithologic units with no table (`read_lithology_units_on_grid`: text, integer and mixed values, crosswalk with and without names, a missing value, one code with two names, an unknown field, unmapped and grid-nodata pixels, a code raster on the grid and warped to 20 m, a non-integer code). Driver `tests/cpp/test_lithology_erodibility.cpp` |
| `smoke_/failure_lem_lithology_erodibility_*`, `regression_lem_lithology_erodibility_*` | `smoke`/`failure`/`regression`, `cli;lem` | the CLI on fixtures from `fixture_lem_lithology_erodibility_cli`: own grid, DEM grid, centre clip, code raster, crosswalk, a named unmapped code; seven refusals; the K rasters compared bitwise against expected rasters |
| `smoke_lem_lithology_erodibility_chain` | `smoke;cli;lem` | `lsdtt-lem` on the produced K raster with the printed parameter lines (converts the map-area K by dx^m) |

The fixture map's unit edges lie on pixel boundaries and the clip point is
chosen so no pixel centre lies on a unit edge, so the expected K rasters follow
from the layout alone. The rasteriser tests use a 20 × 20 grid with three
overlapping pixel centres (resolved to the lower value), a 0.6-cell sliver gap
(10 pixels filled by coverage), a one-pixel hole (left unmapped), two fragments
of one value, reference nodata, and the same layout in EPSG:4326 (interior pixels
exact on the UTM grid).

### 8c.3 k_sn by lithologic unit (`lsdtt-chi-analysis` with `lithology_fname`)

| Test | Labels | What it covers |
|---|---|---|
| `unit_chi_lithology_channel_stats` | `unit;chi;lithology` | `compute_weighted_quantile` (lower weighted quantile: fixed vectors at q = 0, 0.25, 0.5, 1; zero weights; empty = NaN; zero total, negative weight, size mismatch and q outside [0, 1] throw; equal weights = `compute_quantile_floor_in_place` at q = 0, 0.5, 1 for n = 1..12); contact distance on a two-unit grid (exact, incl. nodata and unmapped); node -> segment on a hand-built fit (confluence, outlet); the lithology of a set of nodes (majority, purity with dx and dx sqrt2 weights, tie to the lower code, all unmapped, zero weights, empty); lithologic reaches on a hand-built network (3 reaches, 12 / 18 / 9 nodes, end types contact/outlet, source/contact, source/junction; with a 20 m buffer 4 nodes excluded and the ends become buffer); reach OLS on z = 3 + 7 chi (k_sn 7, R^2 1, SE 0) and `lithology_min_reach_pixels`; segment purity (180/290); summary bookkeeping (rows, sums, `low_sample` below 5 qualifying reaches with a NaN interval); the bootstrap (same seed bitwise, other seed differs, constant values, basin clusters, 0 replicates); the fit result unchanged; the known-answer steady states for n = 1 and n = 2 in two layouts; the `--check-ksn-by-lithology` mode; relative erodibility (R1-R4, below). Driver `tests/cpp/test_lithology_channel_stats.cpp` (156 checks) |
| `fixture_chi_lithology_cli` | `fixture;chi;lithology` | writes the steady-state DEM, its unit polygons, two crosswalks, and a two-unit code raster beside a copy of `tiny_dem.asc` |
| `smoke_chi_lithology_ksn` | `smoke;cli;chi;lithology` | the CLI on the steady state: all seven outputs, the code raster holds 1 and 2, the run-local codes 1 and 2 in `_lithology_units.csv`, the unmapped row is empty, the reach columns, and each unit has no node in the other's histogram peak bin |
| `rigorous_chi_lithology_ksn_ratio` | `rigorous;chi;lithology` | the reach-median ratio sandstone/shale from that run against 4, within 2.42 % |
| `smoke_chi_lithology_crosswalk` | `smoke;cli;chi;lithology` | codes 10 and 30 and names from a crosswalk's `name` column, a 60 m contact buffer (52 nodes set aside), the bootstrap by basin |
| `rigorous_chi_lithology_observation_only` | `rigorous;cli;chi;lithology` | the `chi_analysis_smoke` run plus a lithology map keeps the fit-nodes hash that `rigorous_chi_nodes_sha256` pins |
| `failure_chi_lithology_{without_fit,missing_file,bad_field,crosswalk_gap,bad_cluster}` | `failure;cli;chi;lithology` | `lithology_fname` without `fit_chi_networks: true`; a map that does not exist; a field the map lacks (the fields are listed); a value missing from the crosswalk (with its pixel count); `lithology_bootstrap_cluster` other than reach or basin |
| `failure_chi_lithology_no_overlap_{raster,vector}` | `failure;cli;chi;lithology` | a code raster and a polygon map 200 km east of the DEM (same CRS) label no valid DEM pixel: refused, the message giving both extents and CRSs |
| `smoke_chi_lithology_multi_mn`, `rigorous_chi_lithology_bootstrap_per_mn` | `smoke`/`rigorous`, `cli;chi;lithology` | the bootstrap stream is per m/n (`std::seed_seq{seed, round(m/n * 10000)}`): with m/n 0.4 and 0.45 the 0.45 summary, unit x basin and reach tables equal those of the single-m/n run byte for byte |
| `smoke_chi_lithology_relative_erodibility` | `smoke;cli;chi;lithology` | `write_relative_erodibility` on the steady state (n = 1 and 2, E = 2e-4 m/yr): the relative-erodibility CSV header and rows, K units strings, no unmapped row, both draft property tables with their header lines and provenance, the default reference and the conditions printed |
| `rigorous_chi_lithology_relative_erodibility_observation_only` | `rigorous;cli;chi;lithology` | the seven k_sn-by-lithology outputs of that run equal those of `smoke_chi_lithology_ksn` byte for byte (the relative-erodibility stream is separate) |
| `rigorous_chi_lithology_relative_erodibility_{ratio,reference}` | `rigorous;chi;lithology` | `--check-relative-erodibility` on that CSV: shale K/K_ref against 4 within 2.42 % and K against 4e-5 within 6.4 %; sandstone (the reference) 1 and K against 1e-5 within 4 % |
| `smoke_chi_lithology_relative_erodibility_k_raster` | `smoke;cli;chi;lithology` | `lsdtt-lithology-erodibility` reads the n = 1 draft table unchanged and, with the fixture DEM and map, writes a 160 x 120 K raster holding the table's two K values |
| `smoke_chi_lithology_relative_erodibility_multi_mn`, `rigorous_chi_lithology_relative_erodibility_per_mn` | `smoke`/`rigorous`, `cli;chi;lithology` | with m/n 0.4 and 0.45 the 0.45 relative-erodibility CSV and both property tables equal those of the single-m/n run byte for byte |
| `failure_chi_lithology_relative_erodibility_{no_map,bad_reference}` | `failure;cli;chi;lithology` | `write_relative_erodibility` without `lithology_fname`; a reference name that is not a unit (the units are listed) |

**Known answer.** A 120 x 160 grid at 30 m (a tilted red-noise surface, seed 11, D8,
open edges, draining north) with K1 = 1e-5 and K2 = 4e-5, U = 2e-4, map area, steady
state from `compute_d8_fluvial_steady_state` with an erodibility field; chi (m/n =
theta, A0 = 1) and the fits on the same flow tree. Two layouts: A, the contact a
vertical line (column 80) parallel to the drainage; B, a horizontal line (row 60)
across it, so every trunk crosses it. The expected medians are (U/K)^(1/n): 20 and 5
for n = 1 (m = 0.45), sqrt(20) and sqrt(5) for n = 2 (m = 0.9).

Measured 2026-10-17 (the relative errors are identical for n = 1 and n = 2):

| Statistic | Layout | Unit K1 | Unit K2 | Ratio | Pins |
|---|---|---|---|---|---|
| node median (first cut) | A | +1.14 % | +1.25 % | 0.108 % | medians 2 % (cap), ratio 0.25 % |
| reach median (OLS per reach, 10 nodes, buffer 0) | A | +1.759 % | +1.310 % | 0.443 % | medians 2 % (cap), ratio 0.9 % |
| reach median | B | +1.192 % | +1.581 % | 0.383 % | medians 2 % (cap), ratio 0.8 % |

The medians sit high because the steady state integrates chi at each cell's own
area and `lsdtt-chi-analysis` trapezoidally. Twice the measured median error would
exceed the 2 % cap, so the medians are pinned at the cap; the ratios are pinned at
about twice the measured error, under the 1 % cap. Reach counts are pinned exactly:
layout A 46 (K1) and 42 (K2) reaches (16 and 14 qualifying), layout B 57 and 38 (18
and 20 qualifying); B has 20 contact-ended reaches against 8 for A, and 1 (n = 1) or
7 (n = 2) of its fitted segments span the contact.

Through the CLI, which re-routes flow on the steady DEM: node medians 20.78 and 5.22;
reach medians 20.403 and 5.163, ratio 3.9515, 1.212 % from 4, pinned at 2.42 %.

**West Virginia (Allegheny Front, live data, not a CTest).** `wvfront_DEM.tif` (COP30,
EPSG:32617, 30 m, -79.45/-79.15/38.85/39.05), Macrostrat SGMC units, 105 outlets from
`lsdtt-basin-extraction`, threshold 200 pixels, `min_channel_pixels` 30, m/n 0.45,
A0 = 1, 1000 replicates by reach, seed 5489. `lith_primary` with a crosswalk (codes
10/20/30/40/90): 15847 channel nodes, 0 unmapped; 967 reaches, 531 qualifying.
Reach median k_sn [95 % interval], buffer 0 / buffer 120 m:

| Unit | Reaches (qualifying), buffer 0 | Median [interval], buffer 0 | Reaches (qualifying), 120 m | Median [interval], 120 m |
|---|---|---|---|---|
| sandstone | 358 (190) | 65.34 [60.01, 69.59] | 265 (126) | 64.56 [58.64, 69.43] |
| siltstone | 68 (44) | 54.09 [45.86, 60.09] | 60 (39) | 55.99 [40.84, 66.77] |
| shale | 334 (177) | 41.14 [36.66, 48.97] | 277 (134) | 36.23 [31.30, 42.34] |
| limestone | 196 (113) | 57.84 [51.62, 64.87] | 167 (68) | 44.31 [38.45, 54.22] |
| alluvium | 11 (7) | 8.63 [4.65, 14.69] | 9 (5) | 7.02 [3.00, 15.55] |
| all | 967 (531) | 54.58 [51.05, 58.22] | 778 (372) | 48.43 [43.72, 53.23] |

(Intervals from the per-m/n bootstrap stream, seed 5489; the medians do not depend
on the seed.)

The 120 m buffer sets aside 4677 of 15847 nodes. With `name` (no crosswalk), 28 map
units are on the grid, 27 with channel nodes, 21 with at least 5 qualifying reaches.
111 of the 967 reaches (2 qualifying) have |dz| < 0.5 m.

**Relative erodibility from k_sn (`write_relative_erodibility`).** K_i / K_ref =
(k_ref / k_i)^n from the units' reach-median k_sn, and with an erosion rate E,
K_i = E / (A0^m k_i^n) (m = (m/n) n, map area, units m^(1 - 2m) yr^-1). The interval
resamples the ratio inside each bootstrap replicate (reaches drawn per unit, or basins
drawn once for every unit), from `std::seed_seq{seed, round(m/n * 10000), 1}`, a
stream separate from the k_sn summary's. Driver checks (`unit_chi_lithology_channel_stats`):

* R1 (exact): the default reference is the unit with the most qualifying reaches, then
  the larger channel length, then the lower code; unmapped and all are never chosen and
  are refused by name; an unknown name throws listing the units; a named unit with no
  qualifying reach, and no unit with one, throw; the units string (m^(0.1), m^(-0.8),
  m^(0) yr^-1) and the n file suffix.
* R2 (exact, hand-built reaches): the point ratio (k_ref/k_i)^n for n = 1 and 2; same
  seed, bitwise-identical intervals; the interval brackets the point ratio; the n = 2
  interval is the n = 1 replicates squared (ends at or above the squared n = 1 ends,
  within 0.1 %); the reference row 1 [1, 1]; K = E/k^n with an interval bracketing it;
  the stream differs from the k_sn summary's and from m/n 0.4, and the 0.45 interval is
  unchanged by running 0.4 first; constant reaches give [4, 4] and K [E/5, E/5]; basin
  clusters with the second unit in one of four basins use 672 of 1000 replicates for
  it (expected about 684) with the interval [4, 4]; fewer than 5 reaches gives
  low_sample with NaN intervals; a low_sample reference gives NaN ratio intervals and
  K intervals; n <= 0 throws.
* R3: the draft table loads with `load_lithology_property_table` (map_area, m 0.9,
  n 2, K equal to the CSV as float, the provenance text); a unit with no reach has no
  row and is listed in `# not_estimated`, and the table still loads; two names equal
  ignoring case and E = 0 are refused.
* R4: K2/K1 and K with E = U on the T10/T11 steady states, both layouts, default
  reference (west in A, east in B, by qualifying reaches 16 / 14 and 18 / 20).

Measured 2026-10-19 (relative errors; K/K_ref is 4 against west, 0.25 against east):

| n | Layout | K/K_ref | K west (K1) | K east (K2) | Pins: ratio, K |
|---|---|---|---|---|---|
| 1 | A | 4.0177 (0.443 %) | 9.827e-06 (1.728 %) | 3.948e-05 (1.293 %) | 0.9 %, 2 % |
| 1 | B | 0.25096 (0.384 %) | 9.882e-06 (1.178 %) | 3.938e-05 (1.556 %) | 0.8 %, 2 % |
| 2 | A | 4.0355 (0.887 %) | 9.657e-06 (3.426 %) | 3.897e-05 (2.570 %) | 1.8 %, 4 % |
| 2 | B | 0.25192 (0.770 %) | 9.766e-06 (2.342 %) | 3.876e-05 (3.088 %) | 1.6 %, 4 % |

K and K/K_ref are powers -n of the reach medians, so the median errors (above) are
multiplied by n. The pins are twice the measured error, capped at n times the median
caps (1 % on a ratio, 2 % absolute): K for n = 1 sits at the 2 % cap, K for n = 2 at
4 %, and the n = 2 ratios at 1.8 % and 1.6 %, above the 1 % median-ratio cap. The
K intervals (1000 replicates by reach) do not contain the generating K: for n = 1,
layout A, K west [9.760e-06, 9.877e-06] against 1e-5; the reach medians are biased
high by the rectangle-versus-trapezoid chi difference above, and the bootstrap
does not see a bias.

Through the CLI (n = 1, E = U, reference sandstone by 15 against 11 qualifying
reaches): shale K/K_ref 3.951501 (1.212 %), K shale 3.8734e-05 (3.165 %), K sandstone
9.8024e-06 (1.976 %); pinned at twice: 2.42 %, 6.4 %, 4 %.

West Virginia, the same runs as above with relative erodibility on (`wvrk_*` in the
data directory): the k_sn-by-lithology outputs, fit nodes and fit segments are byte-identical
to the runs without it. `lith_primary`, reference sandstone (190 qualifying reaches,
the most), 1000 replicates by reach; with an illustrative uniform E = 3e-5 m/yr
(30 m/Myr; not a measurement):

| Unit | Qualifying reaches | K/K_ref, n = 1 | K/K_ref, n = 2 | K (n = 1), m^(0.1) yr^-1 | K (n = 2), m^(-0.8) yr^-1 |
|---|---|---|---|---|---|
| sandstone | 190 | 1 [1, 1] | 1 [1, 1] | 4.591e-07 [4.309e-07, 4.966e-07] | 7.026e-09 [6.188e-09, 8.219e-09] |
| siltstone | 44 | 1.208 [1.032, 1.488] | 1.459 [1.065, 2.214] | 5.546e-07 [4.965e-07, 6.536e-07] | 1.025e-08 [8.217e-09, 1.424e-08] |
| shale | 177 | 1.588 [1.303, 1.864] | 2.523 [1.698, 3.476] | 7.292e-07 [6.126e-07, 8.370e-07] | 1.772e-08 [1.251e-08, 2.335e-08] |
| limestone | 113 | 1.130 [0.9755, 1.281] | 1.276 [0.9516, 1.642] | 5.186e-07 [4.646e-07, 5.850e-07] | 8.966e-09 [7.195e-09, 1.141e-08] |
| alluvium | 7 | 7.573 [4.196, 14.65] | 57.35 [17.61, 214.5] | 3.477e-06 [2.043e-06, 6.445e-06] | 4.030e-07 [1.391e-07, 1.384e-06] |

The ratio columns of the E run equal those of the run without E byte for byte.
`lsdtt-lithology-erodibility` read the n = 1 table with the same map and crosswalk and
wrote the K raster (legend: the five K values above). With `name` (n = 1), the default
reference is Greenbrier Group (51 qualifying reaches); 27 of 28 units have a
qualifying reach, 21 have at least 5.

### 8c.4 Channel-network statistics by lithologic unit (`lsdtt-channel-network-analysis` with `lithology_fname`)

Lithology is a grouping key of the bearing and junction-angle statistics
(`channel_bearing_group_by` / `_stratify_by`, `junction_angle_group_by` /
`_stratify_by`: `lithology`), so its rows are in the same summary, rose and
comparison files as every other key. A link takes the length-weighted majority unit
of the nodes its bearing was fitted on (each node weighted by its D8 step, ties to the
lower code). No record is dropped for impurity: a link below
`lithology_min_majority_fraction` (default 0.8) goes into the group
"mixed (<80% one unit)" (the percentage follows the parameter; code -3), and every
record keeps its majority unit and fraction in the per-record columns. A junction takes
the unit of its pixel (`group_key` lithology) and, as a second selection
(`lithology_interior`), the same unit when its three sampled arms are all in it with
purity at or above the threshold, else the group "arms disagree" (code -4). Mixed and
arms disagree are **descriptive only** (`channel_record_categories::descriptive_only`):
they have their own summary and rose rows (note "descriptive only: not compared") but
are never compared, never part of "rest", and never compared with a reference, because
they have no single unit and pooling them into "rest" would put a contact effect into a
lithology contrast. One stdout note per run gives how many records that leaves out of
the comparisons. Unmapped records are a group like any unit.
Lithology-keyed junction rows carry a 95 % percentile bootstrap interval of the median
(the junction is the resampling unit, none below 5 junctions, one random stream per
run). In `_junction_angle_comparisons.csv`, which has no group-key column, lithology
rows carry `group_key=lithology` or `group_key=lithology_interior` in `note`. With
`lithology_fname`, `_channel_link_bearings` and `_junction_angles` gain the unit
columns after all others; without it every record output is unchanged.

A lithology or region map that labels no valid DEM pixel is refused, with the map's
and the DEM's extents and CRSs (for a raster, the warp's "do not overlap" is caught and
reported the same way). A partial cover gives a note with the unlabelled fraction and a
WARNING above 10 %.

| Test | Labels | What it covers |
|---|---|---|
| `unit_lithology_channel_network_stats` | `unit;network;lithology` | on the in-memory trellis (2400 junctions, 2396 links): with the contact parallel to the trunks every link takes the unit of its side with purity 1 (1198 / 1198); with it across them the 12 crossing links' unit and purity equal a length-weighted computation in the test, the 4 links below 0.8 are the group "mixed (<80% one unit)" (none excluded), and a 4-node sample at the upstream end makes the 8 crossing links longer than 4 nodes pure (the tag follows the sampled nodes); mixed has its own summary row, the group n sum to every link, and mixed is in no comparison row and not in "rest" (north vs rest compares with the south links only); unmapped links are the group unmapped; 1192 junctions in lithology and in lithology_interior (1180 in their unit, 12 arms disagree); the axial rose by unit equals counts made in the test (596 / 602); the junction median bootstrap; the writers (appended columns, comma quoting); the analysed-area mask; the block-order rule on synthetic counts (auto picks the largest k with >= min, override with a note, nothing reaching min = no blocks with the note, nan p); block permutations (an extreme separation over 40 blocks gives p = 1/(1 + 199); same seed, same p; within-block nan when no block holds both groups; on the trellis by order with order-1 blocks, order 2 has no blocked record = nan with a note; with order-2 blocks the statistic equals Kuiper V on the blocked records, 4 unblocked links left out). Driver `tests/cpp/test_lithology_network_stats.cpp` (30 checks) |
| `smoke_channel_network_analysis_lithology` | `smoke;cli;network;lithology` | the trellis fixture with its north/south code raster as the map: appended record columns, the interval and block columns, the mixed summary row, the arms-disagree junction rows, no comparison row naming mixed, the stdout note, the unit legend (24000 / 23636 pixels in analysed basins), the block-order counts file |
| `smoke_channel_network_analysis_lithology_observation`, `rigorous_channel_network_analysis_lithology_observation_only` | `smoke`/`rigorous`, `cli;network;lithology` | a lithology map that no key uses leaves the six statistics outputs byte-identical to `smoke_channel_network_analysis_bearing_statistics_order` (`tests/cmake/compare_output_sets.cmake`) |
| `failure_channel_network_analysis_lithology_{bad_source,key_without_map}` | `failure;cli;network;lithology` | `lithology_source: shapefile`; `channel_bearing_group_by: lithology` without `lithology_fname` |
| `failure_channel_network_analysis_{lithology,region}_no_overlap_{raster,vector}` | `failure;cli;network` | lithology and region maps 200 km east of the trellis (same CRS), as a code raster and as a polygon layer: refused with both extents and CRSs |

Without the statistics switches, 32 outputs of two full runs (the lake and Allegheny
Front DEMs, every record, Tokunaga, area-ratio, area-scaling and vector output) are
byte-identical to those before the lithology and permutation work. With them, the
summaries gain `block_order, block_order_mode, min_blocks, n_unblocked, note`, the
comparisons fill `permutation_p` and gain `permutation_statistic`, and n_blocks follows
the automatic block order (section 10s.6).

**West Virginia (live data, not a CTest).** `wvbig_DEM.tif` (-79.80/-79.05/38.85/39.15),
Macrostrat SGMC `lith_primary` with the five-class crosswalk (codes 10/20/30/40/90),
threshold 300 pixels, all complete-basin junctions, full-link sampling with at least
5 nodes, automatic block order (order 3), purity 0.8, 1000 replicates, seed 5489,
`one_vs_rest`. Of 2974 links, 732 are "mixed" (by majority unit: sandstone 266,
siltstone 40, shale 250, limestone 159, alluvium 17); no channel node is unmapped. Of
953 junctions, 489 are "arms disagree" in `lithology_interior`.

| Unit (all orders) | Links | n_blocks | Axial mean / R (count) | Axial mean / R (length) | D8 fraction | One-vs-rest axial W1 (deg) / Kuiper V / p_iid / permutation p |
|---|---|---|---|---|---|---|
| sandstone | 585 | 66 | 110.6 / 0.047 | 166.3 / 0.018 | 0.056 | 2.7 / 0.069 / 0.215 / 0.562 |
| siltstone | 318 | 27 | 114.8 / 0.138 | 121.4 / 0.207 | 0.053 | 3.0 / 0.083 / 0.291 / 0.959 |
| shale | 945 | 82 | 113.4 / 0.124 | 116.4 / 0.149 | 0.048 | 2.5 / 0.065 / 0.155 / 0.678 |
| limestone | 343 | 46 | 63.7 / 0.023 | 17.2 / 0.044 | 0.076 | 4.5 / 0.118 / 0.008 / 0.863 |
| alluvium | 51 | 5 | 110.9 / 0.322 | 114.3 / 0.415 | 0.039 | 9.7 / 0.236 / 0.062 / 0.200 |
| mixed (<80% one unit) | 732 | 85 | 110.7 / 0.262 | – | – | not compared |

| Unit | Junction alpha median [interval], lithology (n) | lithology_interior (n) |
|---|---|---|
| sandstone | 62.3 [60.3, 68.2] (284) | 64.8 [60.1, 72.5] (126) |
| siltstone | 66.2 [58.9, 72.4] (113) | 70.7 [60.3, 76.9] (76) |
| shale | 60.3 [57.0, 65.5] (385) | 58.4 [54.5, 60.9] (206) |
| limestone | 63.9 [59.3, 69.6] (155) | 66.9 [61.2, 80.0] (53) |
| alluvium | 43.6 [28.0, 83.1] (16) | 19.8, no interval (3) |
| arms disagree | – | 62.0 [60.0, 66.2] (489) |

### 8c.2 Macrostrat lithology fetch (`lsdtt-fetch-lithology`)

| Test | Labels | What it covers |
|---|---|---|
| `unit_fetch_lithology` | `unit;fetch` | tile maths against hand-computed tiles (zoom 0, Berlin 550/335 at zoom 10, polar clamp, bounds); UTM zones; the three API parsers incl. malformed JSON, an `error` response and a missing `success/data`; the `lith_primary` rule; a fetch from the committed tile: the unit outside the box dropped, EPSG:32632, areas within 1.5 % of the source rectangles projected independently, the vocabulary, the attribution; an explicit EPSG to GeoPackage; too many tiles, no readable tile and a box with no units refused; an unreachable API (127.0.0.1:9) is a warning with the output still written. Driver `tests/cpp/test_macrostrat_fetch.cpp`; skips (77) without the MVT driver |
| `smoke_fetch_lithology_offline` | `smoke;cli;fetch` | the CLI on the committed tile with `lookup_map_units: false`: all four outputs, vocabulary rows and sources, no unit from outside the box |
| `failure_fetch_lithology_{no_bbox,too_many_tiles,web_mercator}` | `failure;cli;fetch` | a missing edge, a box needing more than `max_tiles`, `output_epsg: 3857` |
| `smoke_fetch_lithology_help` | `smoke;cli;fetch` | `-h` lists the Macrostrat parameters |
| `network_fetch_lithology_cli` | `network;cli;fetch` | a live fetch near Edinburgh (two tiles, BGS 1:625k), only with `LSDTT_ENABLE_NETWORK_TESTS` |

The fixture is one synthetic vector tile, `tests/data/macrostrat_tiles/10/540/364.pbf`,
generated with `ogr2ogr -f MVT` from `units_source.geojson` beside it (the
command is in that directory's README). Its vertices are quantised to 1/4096 of
the tile (9.6 m), which sets the area tolerance. MVT and OpenFileGDB are in the
vector drivers `unit_format_support` requires.

### 8c.5 OCTOPUS erosion rates and the first pieces of the erodibility calibration (2026-10-19)

Calibrating the stream-power K from basin-averaged erosion rates is being built in
steps. Built and tested so far: the OCTOPUS fetch (`lsdtt-fetch-erosion-rates`,
`lsdtt_octopus_fetch`), area-matched snapping, the York regression, the erosion-rate
sample reader with explicit units, the parameter groups (declared, not yet used by a
CLI), and the K formula and property-table writer shared with relative K. Every test
below also carries the label `erodibility_calibration`.

| Test | Labels | What it covers |
|---|---|---|
| `unit_fetch_erosion_rates` | `unit;fetch` | 42 checks, no HTTP. The endpoint and layer constants; the latitude-first bbox string; the hits, page and outlets-by-id URLs as exact strings, none with `maxFeatures`, all on `http://geoserver.octopusdata.org`; the paging plan (28 at page 10: startIndex 0, 10, 20); `numberMatched` from the canned hits XML (28) and from GeoJSON; the ExceptionText; the sentinel predicate (true for -9999, -9999.99, -99.99, -99.999, -9.999, -9.99999, -999; false for -9998.9, 0, -9.9, -99, -9, data, NaN); the rate-field map for all four schemes × two nuclides. Offline reads of the large box: 28 basins, 23 outlets and rows, 94 fields (93 attributes + `id`), 896 values nulled in 32 fields (28 each, including every Al-26 field, `GLA_PCNT`, `SMP_*`; none in `EBE_MMKYR`); the S186WTS046 row (38.88961897, -79.40200545, 12.88 ± 2.62 mm/kyr, 20.14 km²) and its internal uncertainty `ERRBE_AMS × EBE_MMKYR / EBE_GCMYR`; S186WTS001 left out (outlet outside); S046WTS009-012 four rows of 9.45 km²; the CSV header, no sentinel in the CSV; the polygons (28 in EPSG:32617, 97 fields, Al-26 rates null); the DEM box (21 basins: -79.719167, 38.270833, -78.51, 39.729167 with a 0.02° pad); the attribution (every DOI, both study lines, DBVERNO v2.3(4)); lsdn (EBEL 12.8 ± 0.83); al26 (all 23 dropped, no rate); the small box (5 basins, 1 row); `max_features` 5 against the hits fixture and an ExceptionReport body both fail. Driver `tests/cpp/test_octopus_fetch.cpp` |
| `smoke_fetch_erosion_rates_offline` | `smoke;cli;fetch` | the CLI on `wv_large_box`: the four outputs, 28 polygons, the exact CSV header and S186WTS046 row, no S186WTS001 and no `-9999` in the CSV, the gchron DOI and the S186 study line, `Read 28 OCTOPUS basins (23 with an outlet in the box)`, `samples written: 23`, EPSG:32617 |
| `smoke_fetch_erosion_rates_offline_lsdn` | `smoke;cli;fetch` | `octopus_rate_scheme: lsdn`: the S186WTS046 row with EBEL 12.8, EBEL_EXT 0.83, EBEL_INT 0.29 |
| `smoke_fetch_erosion_rates_offline_empty` | `smoke;cli;fetch` | an empty FeatureCollection: header-only CSV, the attribution, the warning, exit 0 |
| `smoke_fetch_erosion_rates_{raster,vector}_bounds` | `smoke;cli;fetch` | `bounds_fname` as a GeoTIFF (`raster, EPSG:32632`) and as a polygon FlatGeobuf (`polygon, EPSG:4326`) |
| `failure_fetch_erosion_rates_{no_bbox,too_many,bad_scheme,web_mercator,exception}` | `failure;cli;fetch` | no box; `max_features: 5` against the canned hits (28); `octopus_rate_scheme: cronus`; `output_epsg: 3857`; an ExceptionReport body |
| `smoke_fetch_erosion_rates_help` | `smoke;cli;fetch` | `-h` lists `octopus_rate_scheme` |
| `network_fetch_erosion_rates_cli` | `network;cli;fetch` | a live fetch of the wvbig box: S186WTS046 in the CSV, at least 5 polygons, the attribution; only with `LSDTT_ENABLE_NETWORK_TESTS` |
| `unit_chi_erodibility_calibration` | `unit;chi` | 32 checks. York (below); the sample reader (unit factors exact, `m_per_myr` = `mm_per_kyr`, unknown and empty units listed, mass fluxes refused; a units column that agrees reads, one that disagrees throws naming row 2; missing rate / uncertainty columns throw listing the columns; no id column gives `row_<n>`; `bad_rate` and `missing_uncertainty` rows kept; km² to m²; a UTM GeoPackage point layer reads the same pixels and statuses as the CSV); `stream_power_erodibility_from_erosion_rate` (value, NaN cases, A0 = 1 vs 1e6 equal to 1e-10 relative for n = 1, 1.5, 2); the K units strings for m = 0.45 and 0.9; the shared property-table writer byte for byte and its case-insensitive name check; the two parameter groups (29 keys, defaults, groups, no n-list key); basin completeness at a grid edge (below). Driver `tests/cpp/test_erodibility_calibration.cpp` |
| `unit_flow_snapping` | `unit;flow;snapping` | 13 new checks of `snap_pixel_to_reported_drainage_area_expanding_square` on a synthetic junction (below) |
| `unit_raster_fetch` | `unit;fetch` | 2 new checks: `bounds_fname` pointing at `tests/data/bounds/bounds_raster_utm.tif` returns its WGS84 footprint |

**York regression.** `fit_york_line` (York et al. 2004, r = 0) on Pearson's data with
York's weights, transcribed from Cantrell (2008, ACP 8, 5477, Table 1): slope
-0.480533407446, intercept 5.479910224033, analytic errors 0.057985009 /
0.294970735, MSWD 1.483294149, corr(slope, intercept) -0.963088, 10 iterations. The
published values (Cantrell 2008, Table 2, as printed to 5 significant figures) are
slope -0.48053 ± 0.0706 and intercept 5.4799 ± 0.359; the printed errors equal the
analytic ones times sqrt(MSWD) (0.0706203, 0.3592465), which `fit_york_line` returns
as `se_*_scaled`. The test asserts the published values to their printed precision
and the analytic results against scipy.odr (ODRPACK, scipy 1.18.0: slope
-0.480533392649, intercept 5.479910153327, errors 0.0579850066 / 0.2949707257,
reduced chi-square 1.4832941493) within 5e-8 (slope; ODRPACK's own convergence leaves
1.5e-8), 2e-7 (intercept) and 1e-8 (errors, MSWD). York et al. (2004) itself
(doi:10.1119/1.1632486) could not be read (HTTP 403), so its tabulated values are not
asserted. Also: an exact line y = 2 + 1.5x gives the slope and intercept to 1e-12
and MSWD < 1e-20; sigma_x = 0 equals weighted least squares to 1e-10; fewer than 3
points or a zero x range give `converged = false` without an exception.

**Area-matched snapping** (`unit_flow_snapping`). A 30 × 21 grid at 10 m, "oooo":
a trunk running south down column 10 and a tributary running east along row 12.
Measured contributing areas: tributary (12, 9) 80 px, trunk (13, 10) below the
confluence 310 px, trunk (11, 10) above it 207 px, (12, 10) 218 px. From (13, 9) with a
20 px channel threshold and R = 3: expanding-square snapping takes the trunk (13, 10);
`reported_area` with the tributary's area takes (12, 9). A reported area of 50 × the
trunk's takes the cell closest in area, (16, 10), still with success and a ratio
outside the tolerance. No channel cell within R gives
`threshold_not_found_in_search_window`; off the grid, `input_out_of_bounds`; equally
near in-tolerance cells go to the smaller |ln ratio| (from (11, 9): 150 px picks the
trunk (11, 10), 110 px the tributary (12, 9)); a tolerance <= 1, a non-positive or NaN
area and a zero threshold throw; the dispatcher overload passes `reported_area`
through, and the old overload refuses it.

**Basin completeness at a grid edge.** A 12 × 9 plane falling south, "oooo": the basin
of (8, 4) is 9 pixels including the north-edge cell (0, 4), and
`check_basin_completeness()` is false; the basin of (8, 6) below a 2-pixel nodata block
(5 pixels) is false too.

**The seams shared with relative K.** `stream_power_erodibility_from_erosion_rate`,
`write_lithology_erodibility_property_table` (the relative-K table is a wrapper over
it) and the bootstrap stream registry (k_sn summary: `seed_seq{seed, round(m/n·1e4)}`;
relative K: tag 1; calibration: tag 2) moved into
`lsdtt_lithology_relative_erodibility_cpu`. Every output of the 6 `smoke_chi_lithology_*`
runs (76 files) and of the three West Virginia relative-K runs (`wvrk_lithprimary`,
`wvrk_lithprimary_E`, `wvrk_name`: 44 files, two property tables) is byte-identical
before and after, and so is their stdout.

**`bounds_fname` as a raster.** `bounds_from_vector_file` falls back to a raster when
the file opens as no vector layer: the WGS84 envelope of its footprint, every edge
densified (64 steps). The fixture GeoTIFF's corners transformed by `gdaltransform`
(PROJ, GDAL 3.4.1) give west 10.3352468907, south 46.4553477705, east 10.4649546110,
north 46.5446520523; the test asserts them to 1e-6. The helper now registers the GDAL
drivers itself: before, `bounds_fname` failed with "could not open ... as a vector
dataset" in `lsdtt-fetch-lithology` (which resolves its bounds before any other GDAL
call) even for a vector file.

**Live fetches (2026-10-19, laptop, WSL2; wall times are dominated by the server and the network, and laptop timings are thermally limited, so they are not benchmarks).** wvbig box (-79.80/-79.05/38.85/39.15): 5
basins, 1 with its outlet in the box (S186WTS046, 12.88 ± 2.62 mm/kyr, 20.14 km²), 4
requests, 825,209 bytes, 3.6 s wall. Large box (-80.5/-78.5/38.3/39.7): 28 basins, 23
rows (`EBE_MMKYR` 3.74 to 24.08 mm/kyr, `AREA` 1.06 to 3310.36 km²), 4 requests,
1,028,013 bytes, 3.3 s; its CSV holds the same rows as the offline read of the
committed fixture. The `network_fetch_erosion_rates_cli` command run by hand against
the build: pass, 3.3 s.

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

**Both levers pulled, 2026-09-13, byte-identical.** Two further costs sat in
front of the arithmetic:

1. `data_view()` / `mask_view()` in `lsdtt_raster_surface_metrics_cpu.cpp` were
   documented as non-owning views but built `lsdtt_array2d` objects with the
   `from_ptr` constructor, which COPIES. Every kernel using them duplicated
   the DEM and mask (5 B/cell) serially before doing any work. They now return
   a `raster_view2d<T>` (pointer plus dimensions). This covers all 11 users in
   the file.
2. The per-pixel 9-cell nodata test is replaced, in `compute_hillshade`,
   `compute_gradient_magnitude_3x3_square` and
   `compute_plan_curvature_3x3_square`, by a 1-byte invalid-cell plane built
   once in parallel and a branch-free 9-byte OR. When the raster has no invalid
   cell, the plane is released and the test is skipped entirely.

Single precision was also tried and rejected: it changes output bytes.

Verified byte-identical before/after on 19 calls x 15 inputs (every view user;
nodata holes, nodata border, value==nodata with mask 1, mask 0 with a finite
value, finite and NaN nodata values, 1x1/2x5/3x3 grids, a real DEM) at 1 and 4
threads. Pinned permanently by section C6 of `unit_surface_metrics`, which
compares the three 3x3 kernels with verbatim reference copies of the old
algorithms (memcmp, data and mask). C6 was mutation-tested: dropping one cell
from the window OR fails it. Server result (2026-09-13 forth01 against
2026-09-09 forth02, single thread): **1.56-1.65x** faster on every fixture; see
"Server campaign 2026-09-13" in section 6. Multi-threaded hillshade was not
measurable in that run.

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
ctest -L '^surface_metrics$'               # unit tests, ~1 s, no network
ctest -R benchmark_surface             # vs TopoToolbox, needs pixi + network
ctest -R performance_surface_metrics   # full sweep, ~4 min
```

The benchmark and performance tests are opt-in
(`LSDTT_ENABLE_BENCHMARKS`, `LSDTT_ENABLE_PERF_TESTS`) and skip (77) when pixi,
topotoolbox or the network is missing. Timings are reported, never asserted;
the only assertion is a correlation guard that fails when the two tools turn
out not to be computing the same kind of quantity.

---

## 10c. Multipart geometry, and a silent-drop defect (2026-08-26)

### The defect

`lsdtt_data_loader::read_vector_layer` dispatched on `wkbPoint`, `wkbLineString`
and `wkbPolygon` and had **no branch for any multipart type**. An unrecognised
geometry returned an empty variant with no throw and no warning, and every
consumer then did

```cpp
const auto* line = std::get_if<lsdtt_line_geom>(&record.geometry);
if (!line || line->points.size() < 2) continue;    // silently skipped
```

so the feature vanished and the tool exited **zero**. The same hole existed for
`wkbMultiPoint`, `wkbMultiPolygon` and `wkbGeometryCollection`, and in a SECOND
reader — `load_vector_point_layer_in_raster_crs`, which the point-sampling path
uses instead of `read_vector_layer`.

Measured on shipped binaries, with geometrically identical input in two
encodings:

| CLI | probe | before | after |
|---|---|---:|---:|
| `lsdtt-linear-stats` | output features | **0** (empty layer) | 2 |
| `lsdtt-raster-preprocessing` (points) | sampled points | **0** (0-byte file) | 3 |
| `lsdtt-stream-burning` | burned profile rows | **0** | 68 |
| `lsdtt-channel-extraction` (`stream_vector_path`) | burned profile rows | **0** | 68 |
| `lsdtt-swath-profile` | non-empty outputs | **0** | 2 |

Every "before" run exited 0. Stream burning is the worst of these: burning that
silently does nothing leaves a DEM that looks correct and routes wrongly.

**Scope.** All-single-part files were never affected, which is why ordinary
point layers from QGIS always worked. A file is affected only if it contains a
multipart feature — or is in a format that promotes to multipart (below).

### Format behaviour, measured

| format | all single-part | already multipart | mixed content |
|---|---|---|---|
| Shapefile | `LineString` | `MultiLineString` | **keeps both** |
| GeoPackage | `LineString` | `MultiLineString` | promotes all to `Multi` |
| FlatGeobuf | `LineString` | `MultiLineString` | promotes all to `Multi` |
| GeoJSON | `LineString` | `MultiLineString` | keeps both |

Shapefile has no `LineString`/`MultiLineString` distinction in the format, so a
multi-part PolyLine reads back as `MultiLineString` — and since most public
hydrography and road data is distributed as shapefile, that is the commonest
route by which multipart data arrives. GeoPackage and FlatGeobuf declare one
geometry type per layer, so any multipart content makes the whole layer read
back as `Multi`.

### The fix

Both readers explode multipart geometry into one record per part, recursing
through nested `GeometryCollection`s, with attributes copied onto each part and
the grouping recoverable from `lsdtt_vector_layer::source_feature_index`.

`source_feature_index` is a **parallel vector, not a schema field**. It was
written as a field first, and `unit_internal_drainage_detect` caught it: adding
a synthetic column makes the schema stop describing the file, so every
read-then-write round trip emits something the source did not have.

Unsupported types are now reported by name and count instead of returning
silently.

### Related: `link_length` under periodic boundaries

Independent of the above, and wrong in shipped output. `link_length` was
accumulated as `sqrt(dx² + dy²)` between successive cell centres. Under periodic
boundary conditions a receiver can be on the opposite edge of the grid — the
flow step is one cell, the two centres are a domain apart — so every wrap added
a spurious domain width. On a 250 × 250 periodic landscape, 3 of 97 links
wrapped:

| reported | true | error |
|---:|---:|---:|
| 1730.8 m | 485.8 m | 3.6× |
| 1317.4 m | 72.4 m | 18× |
| 1325.7 m | 80.7 m | 16× |

Any drainage density, total stream length or Hack's-law result on a periodic
domain was wrong by those factors. The length now comes from the row/column
step: a D8 receiver is always adjacent, so the step is one cell width or
`dx√2`, and a wrap is detectable exactly (`|Δrow| > 1 || |Δcol| > 1`) with no
threshold heuristic.

Channel links are now written as uniform `MultiLineString`, split at wraps, so
no line is drawn across the domain. One vertex is lost per link whose *first*
step wraps — a lone vertex cannot be a line — measured at 1 of 1819 on that
landscape; `link_length` is unaffected because it no longer derives from
coordinates.

`*_channel_network_points` gained `flow_distance` and `elevation`, so a
longitudinal profile can be built from attributes alone. Before this it could
not be built from lsdtt3's channel output at all, on any domain: the points
carried neither.

### Tests

| test | covers |
|---|---|
| `unit_multipart_geometry` | 22 checks: all four types, attributes, `source_feature_index`, single-part unchanged, GeoPackage, and a genuinely mixed shapefile |
| `unit_link_length_periodic` | 11 checks: the length invariant, non-periodic unchanged, uniform MultiLineString, no jump inside any part |
| `unit_channel_profile_attributes` | 10 checks: profile reconstructed from attributes on a surface whose gradient is known, recovered as 0.050000 against 0.05 |
| `smoke_multipart_cli_sweep` | 5 shipped CLIs, single vs multipart vs mixed |

**All four were verified to fail with the fix reverted** — the sweep at 0/5.
That mattered: three of its five cases were initially probing something the
defect did not touch (stream burning via `_vector_edges.fgb`, which a
MultiLineString-aware code path already produced; channel extraction via a
point count driven by the area threshold rather than the vector; and a "mixed"
shapefile fixture whose multipart half had only one part and so round-tripped as
a plain `LineString`). A passing check proves nothing until it has been shown to
fail.

---

## 10d. The steady-state SOR tolerance sat below the float noise floor (2026-08-26)

A 250 x 250 steady-state hillslope solve at 5 m was taking ~13 minutes and reporting
`converged: NO`. Optimal SOR needs O(n) iterations, so n = 250 should reach a 1e-4 m
tolerance in roughly 364 sweeps. The gap between ~364 and 250,000 was the whole finding.

**Measurement.** Running the LINEAR solver (which removes the Roering nonlinearity from
suspicion) at increasing iteration caps on the `ocn` landscape:

| iterations | max\|dz\| (m) |
|---:|---:|
| 200 | 2.5e-3 |
| 500 | 1.98e-4 |
| 1000 | 1.83e-4 |
| 2000 | 1.91e-4 |
| 5000 | 1.98e-4 |

The residual falls normally to ~2e-4 by iteration 500 and then stops -- and *wanders*,
non-monotonically. That is an arithmetic floor, not a convergence rate. Confirmed directly
by tolerance bisection against the real binary:

* tolerance 3e-4 -> converged in **243** iterations
* tolerance 1e-4 -> ran **20,000** iterations, finishing at 2.29e-4, i.e. *worse* than it
  had been at iteration 243.

**Cause.** The working array `z` was `std::vector<float>`. The convergence test is
`max|z_new - z_old|`, an increment tending to zero measured on a field of O(100 m), formed
by differencing a four-neighbour sum of O(400). Float spacing at that magnitude is ~6e-5 m,
so the increment could not be resolved below ~2e-4 m. The default tolerance of 1e-4 m was
therefore unreachable in principle: `converged` could never be set, and every solve ran to
`hillslope_ss_max_iterations` no matter how well converged the field was.

**Fix.** `src/cpu/lsdtt_hillslope_steady_state_cpu.cpp` carries the SOR iteration, its
constants and its residual accumulators in `double`; the raster is converted back to `float`
on the way out. Same landscape, same parameters, after the fix:

| case | before | after |
|---|---|---|
| linear, tol 1e-4 | 20,000 iters, not converged | **261 iters, 2.7 s, converged** |
| linear, tol 1e-6 | unreachable at any iteration count | **350 iters, 1.8 s, converged** |

261 and 350 both sit where SOR theory says they should for n = 250. See
[Landscape evolution modelling](../design_and_algorithms/landscape_evolution_modelling.md) for the numerics.

**Separate defect, since RESOLVED: the nonlinear solve did not converge.** Double precision did
not change the Roering case, correctly -- its residual plateaued at ~1.3 m, nine orders above any
precision floor, and the plateau scaled with the relaxation factor (omega 0.3 -> 1.33, 0.8 -> 4.64,
1.3 -> 18.8, auto -> 133). The field was converged (RMS change 0.023 m between 20k and 60k sweeps)
but **106 cells of 62,500 (0.17%)**, at median `S/Sc` 0.735 against 0.365 grid-wide, oscillated at
~1.3 m and held the max-norm up. No cell reached the `s2` clamp.

**Cause: the face coefficient was not a property of the face.** `alpha` for the N face of `(r,c)`
used the E-W cross-slope evaluated at `(r,c)`, while `alpha` for the S face of `(r-1,c)` -- the
same physical face -- used it at `(r-1,c)`. The two differed, making the discrete operator
non-symmetric and the scheme non-conservative: flux leaving one cell through a face did not equal
flux entering its neighbour through it. Non-symmetric operators admit complex eigenvalue pairs,
i.e. modes that oscillate rather than decay.

**Fix.** The tangential slope is averaged *across* the face from the two cells sharing it, which
required splitting the `nbr_z` accessor into an index resolver plus a value read so the
neighbour's own centred slope gets identical periodic-wrap / Neumann / masked-cell handling.

| Roering nonlinear, omega 0.3 | before | after |
|---|---|---|
| result | 30,000 sweeps, plateau 1.33 m, not converged | **5,120 sweeps, 23 s, converged (9.99e-5)** |
| max `S/Sc` in the solution | 0.970 | **0.882** |
| cost per sweep | 2.8 ms | 4.6 ms (1.63x) |

The over-steepening (0.970 -> 0.882) is the conservation defect appearing directly in the answer.
**Physical validation against known forcing:** for a 1-D steady hillslope the divide curvature is
exactly `U/D` as `S -> 0`. The converged field gives a median Laplacian over the flattest 1% of
interior pixels of **-0.007985** against a predicted **0.008** (convex-up is negative in lsdtt3) --
99.8%, tightening toward the prediction as `S -> 0`: 0.007892 at the 5th percentile, 0.007963 at
the 2nd, 0.007985 at the 1st.

**Still open: the nonlinear solver remains unusable at competitive relaxation.** After the fix,
30,000 sweeps still gives 5.02 at omega 0.8 and 23.1 at omega 1.3, and auto-omega (~1.95) diverges.
The auto-omega formula is the optimal-SOR result for the *linear* Poisson problem and has no
derivation behind it for a Picard-lagged nonlinear operator. Being pinned at omega 0.3 costs
roughly 500x against optimal SOR (`rho_GS ~ 0.99984` at n=250 gives ~59,000 Gauss-Seidel sweeps
and ~364 optimal-SOR sweeps; omega 0.3 gives ~196,000). Picard cannot fix this -- it has no
damping on the coefficients at all, so under-relaxing `z` leaves the alpha feedback untouched.
A Newton scheme is the intended resolution.


## 10e. `'o'` was a Neumann boundary in two solvers, and the SOR solver had no test (2026-08-26)

Found while resolving a discrepancy flagged during the hillslope-numerics work: the codebase held
**two contradictory meanings of the `'o'` boundary code**.

* `lsdtt_hillslope_implicit_cpu.cpp:192` had it right -- `is_o_edge()` classifies `'o'` edge cells
  as Dirichlet (fixed elevation at base level), alongside channel cells.
* `compute_hillslope_steady_state` (SOR) and `compute_adi_hillslope_transport` special-cased only
  `'p'`, and treated `'o'` **and** `'n'` alike as Neumann reflection.

So the implicit solver and the SOR solver were solving *different boundary-value problems from the
same input* -- which is why the implicit solver's "match the SOR steady state" test had to dodge
`'o'` edges entirely to pass.

**Decision (Simon): `'o'` is Dirichlet and must never be Neumann.** Both solvers now mirror the
implicit solver's predicate -- SOR skips `'o'` edge cells in the red-black sweep exactly as it
skips channel cells; ADI gives them the identity row in both tridiagonal sweeps.

**This changed answers, and the old ones were wrong.** The hillslope test landscapes declare
`boundary_conditions: opop` -- periodic E-W, *fixed elevation* N-S, which is what was intended --
but were silently getting no-flux N-S. After the fix the N/S edge rows pin at 0.0 m with a proper
gradient rising into them (row 1 spans 0.23-2.39 m). **Any landscape generated with an `'o'` edge
before 2026-08-26 was solved with the wrong boundary condition and should be regenerated.**

The default is also `boundary_conditions = "oooo"`, which previously meant Neumann on all four
edges. With no channels that is an ill-posed problem -- no Dirichlet data anywhere -- and the field
grows without bound. The test below measures it drifting **698 m** in 500 sweeps. The fix makes the
default well posed.

Re-validated on the OCN landscape (Roering nonlinear, omega 0.3): 4,650 sweeps (down from 5,120 --
the extra Dirichlet constraint helps), `max S/Sc` 0.883, divide curvature **-0.007989 against a
predicted 0.008**. Full suite 192/192.

### `unit_hillslope_steady_state` -- the solver's first dedicated test

`compute_hillslope_steady_state` runs every steady hillslope landscape lsdtt3 generates and, until
2026-08-26, **had no test of its own** -- it was exercised only indirectly, as the reference the
implicit solver is compared against. Three defects were found in it in one sitting (10d and this
section). `tests/cpp/test_hillslope_steady_state.cpp` (17 checks, label `unit;hillslope`) pins all
three:

| | check | pins |
|---|---|---|
| H1 | converged field satisfies a **face-symmetric** operator written independently in the test | 10d face asymmetry |
| H2 | linear solve gives `lap(z) = -U/D` at every interior cell | operator correctness |
| H3 | a 1e-6 tolerance is **reachable** on O(100 m) relief | 10d precision floor |
| H4 | `'o'` cells are held; `'n'` cells are not | this section |
| H5 | periodic invariance: `solve(roll(z)) == roll(solve(z))` | seam handling |
| H6 | nonlinear divide curvature recovers `-U/D` | known forcing |

H1 is the important one and the reason it is written as it is. A test that recomputed the
coefficient the way the solver does would only confirm the solver agrees with itself. Instead the
test implements the face coefficient **face-first** -- taking the two cells that share a face and
averaging the tangential slope across it, so `alpha(N face of P) == alpha(S face of N)` by
construction -- then asserts the solver's converged field satisfies *those* equations. The fixture
carries a deliberate cross-slope that varies along the face direction, because the defect is
invisible on a field with a uniform tangent (both cells then report the same value and the
asymmetry cancels).

**Negative controls -- each fix reverted in turn, confirming the tests catch the real defects:**

| reverted | result |
|---|---|
| `'o'` Dirichlet -> Neumann | H4 fails; edge drifts **267,657 m** |
| face-symmetric tangent -> tangent at the updating cell | H1 fails; flux imbalance **575% of U** |
| working array `double` -> `float` | H3 fails; 200,000 sweeps, residual floors at 1.14e-5, never converges |

The 575%-of-U figure is worth keeping in mind: the non-conservative scheme was leaving a flux
imbalance nearly six times the forcing it was supposed to be balancing.

**Gate added to `unit_hillslope_implicit`.** Its test 4 compares the implicit solver against a SOR
reference requested at `tolerance = 1e-6` -- which was **unreachable in float**, so that reference
ran to its 200,000-iteration cap every time and was compared against unconverged. Two checks
(`test4_linear_sor_reference_converged`, `test4_nonlinear_sor_reference_converged`) now assert
`converged == true && final_residual <= tolerance` before any field comparison is made.


## 10f. Newton for the nonlinear hillslope, and one shared operator (2026-08-27)

Follows 10d/10e. Two additions, both driven by the same finding: Picard could not converge the
nonlinear steady hillslope at any usable relaxation factor.

### One shared operator

`src/cpu/lsdtt_hillslope_operator_cpu.hpp` now holds the discrete operator every hillslope solver
must share -- `lsdtt_hillslope_domain` (grid, mask, boundary codes, `neighbour()`,
`is_dirichlet_edge()`) and `lsdtt_hillslope_stencil<T>` (`at`, `slope_x/y`, `face_alpha`,
`face_alphas`). Three private copies of the neighbour accessor had grown up across the
steady-state, transport and implicit solvers, and they had drifted in **three** ways, each a real
defect: face symmetry (10d), the meaning of `'o'` (10e), and what a short boundary string pads
with -- `'o'` in the implicit solver, `'n'` in the other two, so identical input meant different
things. The shared contract pads with `'o'`.

Rewiring the SOR solver onto it is **not bit-identical**, and the reason is understood: the old
code formed `s_n^2/Sc^2 + s_t^2/Sc^2`, the shared operator forms `(s_n^2 + s_t^2)/Sc^2` --
algebraically identical, different rounding association. On the OCN landscape that moves the solve
from 4,653 to 4,650 sweeps and the field by at most **6.98e-4 m on 99 m relief (7e-6 relative)**,
far inside the 1e-4 m convergence tolerance. Batching the four face coefficients so the centre
cell's tangential slopes are evaluated once rather than once per face made the solve *faster*:
**41.5 s -> 20.0 s**.

### The Newton solver

`src/cpu/lsdtt_hillslope_newton_cpu.{hpp,cpp}`, `compute_hillslope_newton`. Parameterised by
`inv_dt`, so the steady problem (`inv_dt = 0`) and one backward-Euler step are **one
implementation** and cannot drift apart. The Jacobian is 9-point and non-symmetric -- the
cross-slope terms reach the diagonals, and the `d(alpha)` terms multiply `(z_nb - z_P)`, which is
antisymmetric -- so `SimplicialLDLT` does not apply; it uses `SparseLU`, with BiCGSTAB+ILUT for
large grids.

`compute_hillslope_newton_transient` adds adaptive sub-stepping. There is deliberately **no CFL
condition**: backward Euler is unconditionally stable, so sub-step size is chosen for accuracy and
nonlinear solvability only.

**Measured, 41x41 cross-slope fixture:**

| | result |
|---|---|
| SOR at omega = 1.5 | fails after 4,000 sweeps, residual 26.6 m |
| **Newton** | **converges in 42 iterations, `norm(F)_inf` = 4e-17 m/yr** |
| outer iterations vs grid | 10 / 42 / 6 for n = 21 / 41 / 81 -- does not grow with n, against SOR's O(n) |
| linear mode | exact in **one** step, matching SOR to 3.8e-6 m |
| periodic invariance | **exact** (0 m) |
| transient -> steady | 15 sub-steps, 49 Newton iterations, Linf 1.9e-6 m |

Newton is judged on the **residual** `norm(F)_inf` -- an erosion-rate imbalance in m/yr -- not on
`max|dz|`. That is a direct consequence of 10d: a criterion built on a difference of large
nearly-equal numbers has a high noise floor, while the residual is independent of elevation
magnitude and suffers no such cancellation. 1e-12 m/yr is comfortably reachable.

**A defect found and fixed during the work: a super-critical initial field deadlocked the line
search.** A strict "every trial step must be sub-critical" filter cannot be satisfied when the
*starting* field is already above Sc -- backtracking drives the trial toward that same infeasible
field, so every step is rejected and Newton gives up at iteration 0 without moving. Users supply
arbitrary DEMs and an initial guess above Sc is ordinary. The filter is now conditional: from a
feasible iterate it refuses to leave the feasible region; from an infeasible one it runs a
restoration step, accepting whatever reduces the slope ratio until the field is back under Sc.
This is why the n=41 case in `unit_hillslope_newton` W3 exists -- it was the case that failed.

`unit_hillslope_newton` (28 checks, label `unit;hillslope`) covers agreement with SOR, known
forcing, quadratic convergence, grid independence, periodic invariance, the transient limit, the
adaptive controller, and the new-time-level Dirichlet contract. That last one (W9) pins doc A's
warning quantitatively: with a channel lowering 50 m over one step, the new-time and old-time
conventions differ by exactly that 50 m and the old-time convention biases the hillslope **high**
by `O(E*dt)` -- a small, entirely plausible-looking error that would be very hard to spot in
output.


## 10g. `'n'` was not a zero-flux boundary, and the explicit solver had no BCs (2026-08-27)

Two more boundary defects, found by writing `unit_hillslope_transport_bc` (8 checks, label
`unit;hillslope`).

### The explicit Roering solver had no boundary-condition concept at all

`compute_roering_hillslope_transport` took no `boundary_conditions` parameter: it clamped indices
unconditionally, so every edge was Neumann, a periodic domain could not be expressed, and an `'o'`
edge was silently no-flux. It is reached from the LEM's explicit path and -- more dangerously --
from the implicit solver's subcycling fallback, so an `'o'` or `'p'` run that dropped to
subcycling **changed boundary condition mid-run with nothing in the output to say so**.

It now takes the string, with **no default** -- a default would let a caller keep the old wrong
behaviour without noticing -- and reads through the shared operator. Both call sites pass the
run's BCs through.

`compute_max_slope_ratio` had the same clamping problem. At a periodic seam it took a one-sided
(zero) difference instead of wrapping and **understated** the slope. Measured on a linear ramp:
periodic reads **4.875**, clamped reads **0.25** -- a 19x understatement, in the direction that
reports a domain as safer than it is, exactly where near-`Sc` behaviour matters. It now takes the
BC string too (and is declared in the header rather than living in an anonymous namespace).

### `'n'` did not implement zero flux

The more interesting one, and the shared operator's fourth drift. A no-flux face requires the
ghost value to equal the cell's own, so the face-normal slope is zero. The four solvers disagreed:

| solver | no-flux rule | correct? |
|---|---|---|
| implicit (`get_face`) | face skipped entirely | yes |
| ADI | face conductance zeroed | yes |
| explicit transport (old) | clamped index -> returns the cell itself | yes |
| **SOR** | `max(0, r - dr)` -> the **opposite** neighbour | **no** |

The SOR rule mirrors about the cell centre. That yields a zero *centred* slope at the boundary
cell, which looks right, but a **non-zero face flux**: row 0 receives flux from row 1 twice, once
through its real S face and once through its fictitious N face. The same accessor was serving both
the centred-slope estimate and the face flux, and it cannot be correct for both.

This was pre-existing (it predates the 2026-08-26 work), and rewiring the explicit solver onto the
shared operator in N1 propagated it *into* a solver that had been correct. The conservation test
caught that immediately, which is the point of writing it.

The shared operator now resolves a no-flux neighbour -- and a masked one -- to the cell itself.

**Measured, on a closed no-flux domain with no channels and no uplift, where mass must be
conserved exactly:**

| | mean elevation drift |
|---|---|
| mirror rule (before) | **0.278 m** |
| cell-itself rule (after) | **7.2e-8 m** |

The OCN test landscape is unaffected (4,649 sweeps, divide curvature -0.007989 against a predicted
0.008) because `opop` contains no `'n'` edge. Runs that used `'n'` edges do change, and the new
answers are the correct ones. Full suite 195/195.

**Test coverage** (`unit_hillslope_transport_bc`): periodic invariance under a column roll (exact,
0 m); `'o'` holds the edge while `'n'` evolves; `'p'` and `'n'` differ at the seam, so `'p'` is
genuinely honoured; mass conservation on a closed domain (the face-symmetry and zero-flux guard);
and `max_slope_ratio` finding a seam-crossing gradient.


## 10h. Newton becomes the default solver; SOR vs Newton on the LINEAR problem (2026-08-27)

`lsdtt-lem` steady-state runs now use Newton for **both** flux laws, and transient
`roering_implicit` steps use the Newton backward-Euler stepper. SOR + Picard is the automatic
fallback (used when Eigen is absent or Newton fails to converge) and can be forced with
`hillslope_force_sor: true`. The explicit `linear` and `roering` models are no longer offered as
user options -- they are redirected with a notice and retained for testing and Eigen-free builds.

**The linear case is the surprise, and it is why the default covers linear too.** A linear steady
solve has no critical-slope difficulty, so the natural expectation is that SOR wins. It does not:
`F` is already linear in `z`, so the Jacobian *is* the operator and Newton reduces to one sparse
direct solve. Measured at `tolerance = 1e-6`, single-threaded, `'nnnn'` edges with channel rows:

| grid | unknowns | SOR | sweeps | Newton | speed-up | `Linf(Newton - SOR)` |
|---:|---:|---:|---:|---:|---:|---:|
| 41 | 1,681 | 8.0 ms | 234 | 6.2 ms | 1.3x | 2.7e-5 m |
| 81 | 6,561 | 92.7 ms | 786 | 17.2 ms | 5.4x | 6.1e-5 m |
| 161 | 25,921 | 1218 ms | 2,188 | 146 ms | 8.3x | 1.8e-4 m |
| 251 | 63,001 | 8374 ms | 5,281 | 393 ms | 21x | 3.7e-4 m |
| 401 | 160,801 | 49,290 ms | 13,322 | 2124 ms | 23x | 9.8e-4 m |

Newton takes exactly **one** iteration at every size.

**Two findings behind those numbers.**

1. **The `omega` cap costs SOR its O(n) property.** The auto value is
   `min(2/(1+sqrt(1-rho_J^2)), 1.95)`, and the optimum passes 1.95 at about n = 160 (1.962 at
   n = 161, 1.975 at n = 251, 1.985 at n = 401). Above that the cap binds and the measured sweep
   count grows with an exponent near 1.97 in n -- essentially Gauss-Seidel, not optimal SOR. The
   cap is deliberate (it guards against the formula overestimating `omega` when channels are dense)
   but on a sparse-channel domain it is expensive.
2. **SOR's convergence criterion understates its own error by ~1000x.** Newton's direct solve is
   exact to round-off, so `Linf(Newton - SOR)` is SOR's error: 9.8e-4 m at n = 401 against a
   requested `tolerance` of 1e-6. This is the standard `error ~ step/(1-rho)` amplification, and
   `rho -> 1` is exactly the regime. It is the same lesson as 10d from the other direction -- a
   step-size criterion is not an error bound.

**Where a direct solve stops paying.** Newton's cost is sparse-LU fill-in dominated. The time
exponent in `N` stays near 1.5 up to a few hundred thousand unknowns, then climbs:

| unknowns | 251,001 | 491,401 | 1,002,001 |
|---|---|---|---|
| Newton | 4.3 s | 9.5 s | 38.2 s |

so `solver_choice` now defaults to AUTO: direct `SparseLU` below
`direct_solve_max_unknowns` (250,000), BiCGSTAB + IncompleteLUT above it.

**On the OCN landscape** (250x250, Roering nonlinear), Newton converges in **22 iterations to
9.6e-16 m/yr in 23 s with no `hillslope_sor_omega` set at all** -- against SOR's 4,649 sweeps which
required `omega = 0.3` by hand. Divide curvature -0.007965 against a predicted 0.008; max `S/Sc`
0.882; no cell at the limiter cap.

**A reporting defect fixed in passing.** `limiter_active` was set if the slope cap bound at any
point *during* the iteration, so a solve that started super-critical, clamped during the
restoration phase, and then converged to a field nowhere near the cap reported "limiter active".
That is a false alarm on a correctness flag -- the flag exists to say that near-`Sc` behaviour in
**the answer** is set by the cap rather than the flux law. It now reflects the converged field, and
the during-iteration signal is reported in the message instead.

### What could not be removed

The explicit solvers and Picard are hidden from the CLI, not deleted, and the reason is structural
rather than sentimental: **neither Newton nor the implicit Picard solver can run without Eigen.**
SOR, ADI and the explicit solvers have no Eigen dependency, `LSDTT_EIGEN` is a supported `OFF`
build option, and `compute_roering_hillslope_transport_implicit` already falls back to
`compute_roering_hillslope_transport` when Eigen is missing. So the explicit Roering solver is the
only nonlinear hillslope transport available in an Eigen-free build, and SOR is the only steady
solver. `cidre` is also retained regardless -- it is the only model that works with an emergent
channel network, since `roering_implicit` structurally requires imposed channels.


## 10i. Emergent coupling unblocked, and a factorisation cache (2026-08-27)

### Two guards blocked the mode they were meant to protect

`lsdtt-lem` can run fluvial incision and hillslope transport on **every cell**, with no channel mask
and no Dirichlet channel set, so that channels emerge where stream power outruns diffusion. Most of
the wiring already existed -- `incision_threshold_area: 0` applies stream power everywhere, and
`hillslope_channel_coupling: coupled` returns an empty channel set on purpose. Two guards blocked it:

* `lsdtt_lem_cpu.cpp` skipped the transient hillslope step whenever the channel set was empty;
* `cli/lsdtt-lem.cpp` made `coupled` + `roering_implicit` a **FATAL** error.

Both rested on "without Dirichlet channel sinks the solve degenerates to sink-free diffusion toward
the edges". That is true of a **steady** solve with no Dirichlet data anywhere -- an elliptic problem
with pure Neumann data is singular, and `unit_hillslope_steady_state` H4 measures it directly (698 m
of drift in 500 sweeps on `'nnnn'` with no channels). It is **not** true of a transient step, which
is what the guards blocked: a diffusion step over `dt` is well posed under any boundary condition,
the fluvial term is the sink in the coupled regime, and an `'o'` edge supplies base level regardless.
The hillslope operator needs boundary conditions, not channels.

Both are replaced by the narrow correct guard -- steady + coupled + no `'o'` edge is fatal.
`failure_lem_coupling_ri_steady_no_dirichlet` pins the singular case;
`smoke_lem_coupling_ri_transient` pins that the combination now runs. The old
`failure_lem_coupling_roering_implicit`, which asserted the wrong restriction, is deleted.

### The factorisation cache: 97% of runtime, removed

The first working coupled run took 37 minutes, and profiling put **97% of it in the hillslope
solve** -- 0.43 s for 20 steps with hillslopes off against 14.83 s with them on, i.e. 0.72 s/step on
40,000 cells for a *linear* solve.

The cause is structural. For a linear run with fixed `dt` and fixed boundary conditions, `alpha` is
1, the degree-of-freedom set never changes and the `I/dt` term never changes, so **the matrix is
identical at every timestep** and only the right-hand side moves. The code was rebuilding ~360k
triplets, sorting them, re-running `analyzePattern` and re-factorising on every one of 3,000 steps.

`lsdtt_hillslope_newton_workspace` (`src/cpu/lsdtt_hillslope_newton_cpu.hpp`, Eigen behind a pimpl
so the header stays Eigen-free) holds the cached matrix and factorisation **for the caller** -- a
cache inside the solve function would be useless, because the LEM makes a fresh call each timestep.
`run_lem_step` accepts one optionally and the CLI owns it across the loop. Reuse depends on the
problem: linear reuses the matrix and both factorisations; nonlinear reuses the sparsity **pattern**
only, since `alpha` depends on `z`. A validity key (grid, boundary codes, DOF hash, `D`, `Sc`, `dt`,
linear/nonlinear) is recomputed every call and the cache dropped if anything moved, so passing a
workspace cannot change an answer.

| 20 steps, 200x200, linear | before | after |
|---|---:|---:|
| wall | 14.83 s | **0.99 s** |

Verified directly rather than inferred: a 400-step run reports `1 pattern analyses,
1 factorisations, 400 solves`.

**Test coverage** (`unit_hillslope_newton` W10/W11). Both halves matter. W10 asserts the cached run
is **bit-identical** (`Linf = 0` over 12 steps) **and** that reuse actually happened (exactly 1
factorisation, 1 symbolic analysis, 12 solves) -- without the second assertion the test would pass
trivially on a cache that silently did nothing. W11 covers the dangerous direction: change `dt`,
require a refactorisation, and confirm the post-invalidation answer matches a fresh uncached solve.

**A bug found while building it.** A workspace-held `Eigen::SparseMatrix` is default-constructed
**0x0**, and `setFromTriplets` does not resize -- it writes outside the matrix, segfaulting. The
non-workspace path never hit this because `J_local(n_dof, n_dof)` is sized at construction. The
matrix is now sized explicitly on first use, with the cache invalidated whenever dimensions change.

### Solver-selection policy moved out of the CLI

`run_steady_hillslope()`, added to `cli/lsdtt-lem.cpp` when Newton became the default, was
solver-selection policy rather than presentation. It is now
`compute_hillslope_steady_state_auto()` in `src/cpu`, performing no I/O: it fills a
`lsdtt_hillslope_solve_report` and the CLI renders it. 65 lines removed from the CLI.

### A measurement trap worth recording

Per-step cost is ~5x higher during the transient than at equilibrium, and the **fluvial-only** cost
is ~3x higher too, so it is not specific to the hillslope solve:

| | fluvial only | with hillslopes | hillslope cost |
|---|---:|---:|---:|
| evolved landscape, 20 steps | 0.41 s | 1.39 s | 0.049 s/step |
| flat start, 400 steps | 24.3 s | 129.4 s | 0.263 s/step |

Plausibly flow routing over an unorganised surface full of pits and flats, but **not confirmed** --
recorded rather than explained. It does not affect the cache result, which is measured directly.

An earlier version of that comparison was **invalid**: the "evolved landscape" run had silently
failed to load its DEM (missing file extension) and done no work, making it appear 30x faster than
it is. `lsdtt-lem` prints a full, plausible configuration block before it discovers the file is
missing. When timing a LEM run, assert on the completion line rather than the absence of an error.


## 10j. Valley spacing, hillslope length, and coupled-timestep control (2026-08-27)

Three capabilities, each added because an emergent-drainage experiment failed without it.

### Hillslope length must be measured by distance transform, not row transects

A row-transect measure (mean run of convex cells along rows) is **not robust to branching**, and the
failure sharpens with the channelization index: oblique tributaries inflate every crossing by
`1/cos(theta)`; a single row mixes short first-order hillslopes with long inter-basin divides so the
mean describes neither; and both grow with `CI`.

`run_valley_spacing` (`src/cpu/lsdtt_valley_spacing_cpu`) measures hillslope length by a
two-pass chamfer **distance transform** from every convex cell to the nearest concave cell. It is
orientation-free and indifferent to branching order. On a five-point diffusivity sweep it beat the
row measure on both criteria:

| measure | exponent | R^2 | implied prefactor c |
|---|---:|---:|---|
| **distance transform** | **0.668** | 0.998 | **0.86** (0.78-0.93) |
| convex run along rows | 0.612 | 1.000 | 3.84 (3.18-4.62) |
| predicted 1/p | 0.690 | | |

`c` ~ 0.86 means `L_h ~ (D/K)^(1/p)` needs essentially no fudge factor; the row measure was inflated
~2.4x because it spans ridge-to-ridge including the valley. The `lsdtt-lem` hillslope-resolution
warning was recalibrated from c = 4.0 to 0.86 and now predicts within 7% of measured.

### Valley spacing is computed in C++ and only DRAWN in Python

`lsdtt-spectral-analysis` gains `run_valley_spacing: true`, reporting three measures side by
side because they fail differently and their disagreement is the diagnostic:

1. **Perron et al. (2008)** power-weighted mean of the spectral peaks -- weighted, not the argmax,
   because spacing occupies a band and an argmax hops between domain harmonics;
2. the **distance-transform hillslope length** above;
3. the small-scale spectral exponent alpha (**Bonetti et al. 2020**, **Hooshyar et al. 2021**),
   predicted to be `2m - 3`.

Alongside the **Bonetti channelization index** `CI = K l^(m+n) / (D^n U^(1-n))`, which for n = 1 is
independent of uplift; `CI_c ~ 32-37` marks first channelization. FFTW is optional in lsdtt3, so the
row-averaged spectrum uses a direct real DFT -- exact on a periodic domain, where zero-padding to a
power of two would not be. Outputs are `_valley_spacing.csv` and `_valley_spacing_spectrum.csv`.

`lsdviztools3` reads those and computes nothing: `plot_valley_spectrum`, `plot_spacing_vs_ci` and
`plot_landscape_panel` (hillshades a whole sweep with each panel's scale bar set to that landscape's
measured length, for a visual check before any number is believed).

### A coupled timestep can be badly wrong while both solvers are unconditionally stable

**This is the finding most likely to bite someone else.** Backward Euler for the hillslope and an
implicit Braun-Willett sweep for stream power at n = 1 are both unconditionally stable, so neither
complains at a large timestep. Stability describes a single operator; the ACCURACY of an
operator-split pair is a different question that nothing was asking.

Measured on a 200x200 grid at 5 m, U = 0.4 mm/yr, from an equilibrated landscape over 500 kyr:

| dt (yr) | grid-scale roughness | hillslope length | relief |
|---:|---:|---:|---:|
| 1,000 | 4.9 % | 34.9 m | 287.7 m |
| 5,000 | 4.4 % | 38.3 m | 280.7 m |
| 10,000 | 4.2 % | 37.8 m | 268.5 m |
| **20,000** | **12.2 %** | **19.7 m** | 211.9 m |

The 20,000 yr run **converged, reported ~99% of the uplift rate, and wrote a plausible landscape**
whose hillslope length was wrong by 5.7x. "Grid-scale roughness" is the fraction of adjacent cells
where curvature flips sign: a few per cent is normal, tens of per cent is an artefact.

**The criterion**: compare uplift per step with the relief across one cell,
`ratio = U*dt / (dx * S)`. Once one step of uplift approaches one cell's own relief, the split stops
resolving incision against diffusion. Keep the ratio at or below ~0.7 (the break above lies between
1.3 and 2.7).

`diagnose_lem_timestep` (`src/cpu/lsdtt_lem_timestep_cpu`) reports the ratio before every LEM run
and warns with a suggested timestep. It is silent on a flat initial surface -- no gradient means no
reference, and a spurious warning on every diamond-square run would train users to ignore the real
one.

`adaptive_timestep: true` enables a controller keyed on the **achieved**
`max_local_elevation_change` rather than `U*dt`, since baselevel fall, a migrating knickpoint or a
post-parameter-change transient all move the surface too. It shrinks faster than it grows (an
oversized step corrupts the landscape; an undersized one only costs time) and never rejects and
re-runs a step (one oversized step does not create the artefact -- it took ~25 -- and rejecting
would double the cost of every adjustment).

**Measured, same landscape and 500 kyr, dt STARTING at the known-bad 20,000 yr:**

| | roughness | hillslope length | relief |
|---|---:|---:|---:|
| start (clean) | 4.7 % | 36.4 m | 287.2 m |
| dt = 20,000 fixed | 12.2 % | 19.7 m | 211.9 m |
| **adaptive on** | **4.3 %** | **37.1 m** | 282.5 m |

It recovers the clean answer from a timestep 4x too large. The controller settled at 2,880 yr; the
static check independently suggested <= 5,232 yr. `unit_lem_timestep`, 11 checks.

### The emergent-spacing result these produced

600 x 600 at 5 m (3 km), OCN-seeded, swept over `CI` = 20-1000 at dt = 5,000. Grid-scale roughness
2.6-7.6% across all seven.

    CI >= 120 (channelized):  L_h ~ D^0.628,  R^2 = 0.996,  implied c = 0.74 (0.64-0.79)
    200 x 200, independent:   L_h ~ D^0.668,               implied c = 0.86 (0.78-0.93)
    predicted:                1/p = 0.690

Runs at `CI` = 20 and 35 sit at or below `CI_c` and legitimately do not follow the law -- an
unchannelized landscape has no spacing set by the D/K balance. `CI` = 60 is above `CI_c` and still
deviates; that is unexplained. **Do not fit across the threshold**: all seven points give `D^0.437`
with `R^2` = 0.88, a number that averages two regimes and means nothing. Every value is one draw
from a band; no repeats were run.


## 10k. A boundary-bookkeeping bug poisoned every LEM convergence statistic (2026-08-28)

**The most consequential defect found in this work, and the hardest to see.**

`run_lem_step` computed its step statistics BEFORE restoring the pinned base-level
boundaries. Base-level cells receive uplift like every other cell and are then reset to their
pinned value, so they never actually move -- but with the restore happening after the statistics,
every one of them contributed a spurious `|dz|` of exactly `U*dt`.

On a 200 x 200 grid with two base-level rows that is 400 cells, **1% of the domain**. Consequences:

* `max_local_elevation_change` sat at **exactly `U*dt`** for entire runs -- 0.4 m for 6 Myr;
* switching the convergence test to the 99.9th percentile changed **nothing**, because 1% is larger
  than 0.1% and the percentile was still inside the poisoned set;
* meanwhile a direct raster comparison of the same run showed the true p99.9 was **0.028 m** with
  **zero** cells anywhere moving by `U*dt`.

So no steady-state test built on either statistic could ever succeed, and the adaptive timestep
controller -- which keys on the same maximum -- was reacting to boundary bookkeeping rather than to
the landscape (223 timestep adjustments in one 600 x 600 run).

**Fix:** restore the boundaries before computing statistics (`lsdtt_lem_cpu.cpp`).

**Everything measured before this fix is suspect**, including equilibration assessments across the
emergent-spacing sweeps and the `mean_erosion` overshoots recorded in section 10j.

**The tell, in hindsight:** a convergence statistic that lands on `U*dt` to full floating-point
precision is describing bookkeeping, not physics. That should prompt a search for an exact
mechanism before any physical story is entertained -- I reached first for an explanation involving a
channel reach near base level, which was wrong.

### Steady state must be expressed as a RATE

With the bug fixed, the same relaxed landscape gives:

| dt (yr) | p99.9 per step | as a rate | % of U |
|---:|---:|---:|---:|
| 250 | 0.0063 m | 2.50e-5 m/yr | 6.3 % |
| 500 | 0.0103 m | 2.07e-5 m/yr | 5.2 % |
| 1000 | 0.0281 m | 2.81e-5 m/yr | 7.0 % |
| 2000 | 0.0290 m | 1.45e-5 m/yr | 3.6 % |

The per-step value scales with dt; the rate does not. **A per-step tolerance can therefore be
satisfied by shrinking the timestep** without the landscape being any steadier, which makes it
unusable as a convergence test whenever dt adapts. `steady_state_rate_tolerance` (m/yr) is
dt-independent and is the criterion to use; the per-step form is retained because it is what users
naturally specify, and both must be met when both are set.

**There is a floor.** The residual sits at **3.6-7.0% of the uplift rate**, independent of timestep
across a factor of eight, set by continual reorganisation of the drainage network. No tolerance
below that can be met however long a run goes.

### Run control added alongside

* `run_until_time` -- a fixed ITERATION count does not cover a fixed model time once dt adapts.
  Measured: a seven-run sweep intended to cover 6 Myr each reached 1.26, 1.12 and 0.70 Myr, because
  one run's dt fell from 5,000 to 326 yr. Runs in a sweep ended at different model times, which is
  an uncontrolled variable.
* `p999_local_elevation_change` in `lsdtt_lem_step_stats` -- a max-norm is hostage to a handful of
  cells (the same lesson as the SOR max-norm in 10d). The steady-state test uses the percentile; the
  TIMESTEP controller keeps the maximum, because an oversized step is a safety question and the
  worst cell is the right thing to be safe about.
* dt-collapse guard -- convergence is refused while dt is below a fraction of its initial value, so
  a shrinking timestep cannot manufacture a small per-step change.
* stall and oscillation detection (`lsdtt_lem_convergence_monitor`) -- distinguishes "converged"
  from "stopped improving", and flags sign-alternation as a period-2 cycle, which points at the
  numerics rather than at needing more model time. A stall WARNS rather than stops: detection is
  heuristic, a false positive discards a healthy run, a false negative only costs time. An earlier
  200-step window with a 10% decay threshold did stop a run that was converging perfectly well.

### Adaptive timestepping destroys a dt-keyed factorisation cache

The Newton workspace (10i) keys its cached factorisation on `inv_dt`, so every timestep adjustment
invalidated it: **28 factorisations over 60 steps against 1 with a fixed step**. Fixed by quantising
dt to a geometric ladder with ASYMMETRIC hysteresis -- descend a rung promptly when the achieved
ratio breaches target, ascend only after a full rung of travel.

Two failed attempts are recorded in the code because both are easy to repeat: applying the growth
factor to the QUANTISED value meant dt could never climb a 2x rung (54 consecutive grow decisions
with dt pinned and the ratio at 0.197 against a target of 0.7); and symmetric full-rung hysteresis
parked dt one rung too high, running at ratio 0.93 -- saving a refactorisation is not worth running
outside the criterion the controller exists to enforce.


## 10l. Compiler warnings, a warning that was load-bearing, and where threading actually is (2026-08-28)

### The build produced thousands of warnings; all were fixable at source

Suppressing them wholesale was rejected: `-Wmaybe-uninitialized` catches real defects elsewhere, and
a warning class turned off globally is one nobody ever turns back on. **The tree now builds with 0
warnings and 197/197 tests pass** (`ctest -LE heavy`).

Four causes, each fixed rather than silenced:

| Warning | Cause | Fix |
|---|---|---|
| `-Wmaybe-uninitialized` (the bulk) | `push_back(v)` on a `std::vector<lsdtt_vector_field_value>` builds a temporary variant and move-constructs it; GCC 13 cannot prove the discriminant is set through that move | `emplace_back` from the *alternative* type (`double`, `int64_t`, `monostate`) -- no variant temporary exists |
| `-Wstringop-overflow` | `checked_mul` only rules out `size_t` OVERFLOW, so its result is provably no tighter than `[0, SIZE_MAX]` | bound `n_pixels` at the allocation |
| `-Wunused-result` | `OGRLayer::CreateFeature` return ignored in a test fixture | checked -- a fixture that silently fails to write a feature makes every assertion below it vacuous |
| `-Wunused-variable` | dead local | removed |

Two details worth keeping, because both cost time:

**`emplace_back` is not automatically the fix.** It helps only when the argument is an *alternative*
of the variant. Where a helper returned `lsdtt_vector_field_value` by value, `emplace_back` would
still move-construct and change nothing; that case needed the helper rewritten to append in place.
Check what is being constructed, not which method is called.

**A guard the optimiser cannot see is not a guard.** The first attempt at the `-Wstringop-overflow`
fix asserted `nrows > 0 && ncols > 0` at the top of the function. It did not work, and the reported
bound (`0xE000000000000000`-`0xFFFFFFFFFFFFFFFF`) said why: the complaint was never about the
dimensions but about `n_pixels`, several hundred lines and an entire FFTW section away, where no
range information survived. Restate the invariant where it is needed.

### A warning was hiding a missing assertion

`-Wunused-but-set-variable` fired on `dt` in `tests/cpp/test_lem_timestep.cpp`. The comment said
"drive it in until it is inside target, then confirm it settles", but the only assertion was
`n_shrinks() >= 1`, which was already true before the loop ran. **The computed `dt` was never
checked** -- the test asserted nothing about the behaviour it claimed to cover.

Adding the assertion made it fail, and the test was wrong, not the controller: it fed a CONSTANT
movement regardless of step size. No real run does that -- movement scales with `dt` -- so the
achieved ratio could never rise to meet the target, and the EMA decaying from the initial overshoot
walked `dt` to the floor with nothing to pull it back. Driving it with movement proportional to `dt`
gives the loop a fixed point; it now recovers to exactly its starting step and holds.

**A `set but not used` on a value a test computed is worth reading as a missing assertion**, not as
lint. Deleting the variable would have silenced the warning and preserved the hole.

### `wavelet_quadratic_fit` has no ctest coverage

The 197/197 pass says nothing about it. Verified instead through the CLI's own oracles: on the
global quadratic -- the EXACT oracle for this fit -- all six coefficients recover to RMSE ~1e-10,
i.e. float round-off, with the periodic-wave approximate oracle giving the expected small errors.

That the check had to be run by hand is the finding. `lsdtt-synthetic-surfaces` already carries
`test_wavelet_on_global_quadratic` with an exact analytic answer; it is not wired into CTest, so
nothing guards this function against regression.

### Where threading is, and where it is not

Established by inspection, not timing, and it bounds what any tuning can achieve:

* `lsdtt_hillslope_newton_cpu.cpp` -- **0** OpenMP pragmas;
* `lsdtt_lem_cpu.cpp` -- **0**;
* nothing anywhere calls `Eigen::setNbThreads` or `initParallel`;
* the Newton path factorises with `Eigen::SparseLU`, which is **sequential**;
* `opts.n_threads` reaches that file only as `omp_set_num_threads` -- setting a count for pragmas
  that are not there.

The coupled transient step routes through `compute_hillslope_newton_transient` for BOTH the linear
and nonlinear cases (the threaded SOR solver is the steady-state one, which a transient sweep does
not use). **So the dominant cost of a coupled run is single-threaded by construction.** Raising the
thread count helps flow routing (`lsdtt_flowinfo_cpu.cpp`, 10 pragmas) and little else. Where runs
are independent -- a parameter sweep -- process-level parallelism is the lever, not threads.

### Two ways to measure nothing

Both were hit while trying to answer "does this scale?", and both produce a plausible number.

**Comparing two machines at the same thread count is not a scaling measurement.** A laptop-to-server
comparison of the same run gave 1.54x, which was read as poor scaling. Both ran at 8 threads:
the ratio is per-core speed and memory bandwidth, and thread count was the one variable that did not
change. It is not evidence about scaling in either direction.

**`nproc` is not the machine's core count.** It reports the CPUs available to the calling process --
affinity mask plus cgroup quota. On a cluster login node it read **2 while `nproc --all` read 192**,
a factor of 96. Every timing taken in that shell, including the server half of the 1.54x above, was
measured on 2 CPUs. Check `nproc` against `nproc --all` before believing any wall-clock number, and
treat a large discrepancy as an instruction to use the scheduler.

### Reproducibility of the sweep harness

Three defects were found in the harness itself, all of the same kind -- a failure that looks like
success:

* **Hand-maintained parameter files silently lost a case.** One `CI` value's files went missing and
  the sweep died at its first step; the only surviving record was an *ingested-parameter dump*,
  which is a CLI output rather than an input. The files are now GENERATED from a single table, with
  a `--check` mode that verifies on-disk files against it. Only `K` is a free input per `CI`: `D`
  follows exactly from `D = K*l^(m+n)/CI`, which reproduces all published values to six significant
  figures.
* **A `.t_*` cleanup glob raced under a job array.** Six tasks share one directory; a task finishing
  early deleted another's in-flight timing file, making it report a spurious failure. Cleanup is now
  scoped to the values that invocation handled.
* **Checking a binary is executable does not establish that it runs.** A binary can be present and
  executable and still fail to load its shared libraries. The preflight now EXECUTES each one
  (`-v`).

On the last point: Eigen is **header-only**, so it is a compile-time dependency that never appears in
`ldd` and cannot cause a runtime failure -- the library most obviously "in the environment" is the
one that leaves no runtime trace, which makes "the environment must be active to run" a natural but
wrong inference. The runtime dependencies are GDAL, PROJ, FFTW and the C++/OpenMP runtimes.


## 10m. Cellular soft-critical hillslopes and sparse tracing

Validated 2026-08-29 against Newton-Roering. On the `7 x 257` near-critical strip, where both
methods use the identical Roering branch, cellular differs from a Newton `dt=0.25 yr` reference by
`4.58e-5 m` Linf and `2.28e-5 m` RMS and is 12.6x faster. Closed volume mismatch is
`7.1e-11 m3`; ridgetop curvature recovers `-U/D` within 2%.

Generated near-critical one-step timings (cellular/Newton) are 0.0626/3.285 s at 250x250,
0.2649/36.133 s at 500x500, and 1.3466/25.825 s at 1000x1000. Four-core cellular speedup is
2.50x, 2.33x and 2.03x. These are thermally limited WSL2 laptop timings, not a native-server ceiling.
A coupled ten-step 250-grid run takes 1.91 s cellular versus 32.08 s Newton.

With 50% nonlocal excess flux, four-hop and eight-hop coupled runs take 2.06 and 2.24 s. Sparse
tracking leaves topography bit-identical. One particle per 100 active cells adds no measurable cost;
one per 10 adds 13% and remains 10.8x faster than bare Newton. Export probability is pinned against
topographic exported thickness, trajectories are exact across one/four threads, and particle CSV
restart round trips preserve all state exactly.

## 10n. The fluvial-surfaces split, and band-then-label (2026-09-04)

`lsdtt-flat-surface-finder` had grown from a morphometric detector into a
109-parameter pipeline that also did flow routing, channel scope resolution,
height above nearest drainage (HAND), floodplain/terrace classification and
vector network output. `lsdtt-valley-centreline` independently re-derived most
of the same channel-relative machinery to reach the part unique to it. The two
overlapped over roughly a third of their parameter surface and disagreed about
which one owned the valley.

The work split them at the point where a channel network becomes necessary:

| | before | after |
|---|---|---|
| `lsdtt-flat-surface-finder` | 109 params, 2696 lines | 52 params, 1422 lines |
| `lsdtt-fluvial-surfaces` | did not exist | ~75 params |
| `lsdtt-valley-centreline` | 68 params | 43 params |
| `lsdtt_valley_centreline_cpu` | 1484 lines | 1051 lines |

### The load-bearing change: band first, then label

The incumbent ordering labelled connected components once over the whole flat
mask and consulted height afterwards. That forces the riser network to be
continuous: a road or gully cutting through a riser lets two treads at different
heights merge into one component, and nothing downstream can separate them
again. Four independent mechanisms had accumulated to attack that one leak.

Banding the height reference FIRST and labelling within each band makes the
merge impossible rather than detectable — two pixels can only share a label if
they share a band, and two treads at different heights are by construction in
different bands. The riser never has to be found and never has to be continuous.

**Measured on a synthetic dissected staircase** (4 treads, 8 m risers, 3
road-like gaps 30 m wide), recorded here because the old path was deleted and
the comparison cannot be re-run:

| | labels | worst label spans | tread purity |
|---|---|---|---|
| old ordering (label, then band) | 2 | all 4 true treads | 0 |
| band-then-label | 8 in 4 bands | 1 tread | 1.0 |

### A fan needs a cone, not a plane or a poly2

HAND is the right reference field for a terrace flight, where treads sit at a
roughly constant height above the channel. It is the wrong field for a fan: a
fan climbs continuously in HAND from toe to apex, so a HAND ladder cuts one
surface into arcs banded by distance from the apex rather than by age. What is
constant on a fan is the depth of an inset lobe below the PROJECTION of the
older surface, which is what a residual to a fitted surface measures.

The obvious fitted surfaces do not work. A fan is radial about its apex, and a
cone is not a low-order polynomial in x and y; over a fan spanning a wide range
of azimuth the misfit is systematic, shaped like the fan, and larger than the
depositional offsets the residual exists to isolate.

**Measured on a synthetic fan** (2 km radius, 0.03 radial gradient, two lobes
inset 4 m below the older surface). The figure is the spread of the OLDER
SURFACE's own residual, which should be flat — so its spread is pure misfit:

| fitted surface | older-surface residual p5..p95 | spread | vs the 4 m inset |
|---|---|---|---|
| plane (RANSAC, 1 m tolerance) | -31.6 .. +2.7 m | ~34 m | swamped |
| poly2 (OLS) | -4.8 .. +5.2 m | ~10 m | swamped |
| **cone (OLS, auto apex)** | **-0.15 .. +0.15 m** | **~0.3 m** | **resolved** |

Consequence, with both arms running the SAME ladder and the same automatic level
selection so that only the field differs: under HAND the ladder chose 12 levels
and cut the older surface into 44 components; under the cone it chose 2, the
older surface stayed 2 components, and each inset lobe came out as exactly one
component at purity 1.0.

`fitted_surface_algorithm_choice` therefore offers `cone | plane_ransac | poly2`
and defaults to `cone`. Plane and poly2 remain right for a piedmont, a terrace
tread, or a broadly curved surface with no single centre.

### The antimode break finder cannot band a signed field

The antimode finder builds its histogram on `log10` of the field with a small
positive floor. That suits HAND, which is non-negative by construction, and
destroys a signed one: every negative value is floored into the lowest bin, so
on a fitted-surface residual the entire below-surface population — the inset
lobes, the thing being looked for — collapses into a single spike and the breaks
returned are meaningless.

`build_reference_field_ladder` throws rather than banding a field containing
negative values with the antimode finder, and `lsdtt-fluvial-surfaces` selects
`jenks` automatically when `fitted_surface` is chosen and no break algorithm was
named — while rejecting an explicit `antimode` with that explanation. Choosing
silently only where nothing was asked for is the line: quietly running something
other than what was requested is worse than refusing.

### How many levels: a constraint, not a trade-off

A level count is ADMISSIBLE when the tightest gap between consecutive band
centres clears `kappa * E_v`, with `E_v` the vertical accuracy of the elevation
product. The answer is the largest admissible count. Where none qualifies the
engine reports UNRESOLVED rather than returning a number, because what may be
claimed about the ladder depends on which branch was taken.

Two decisions inside that are easy to get wrong and are pinned by tests:

- **Band centres are MEDIANS of each band's values, not midpoints of the break
  edges.** The midpoint is geometric and depends only on where the breaks fell;
  the median is where the surface actually sits. Terrace correlation cares about
  the elevation of the tread. This is a deliberate divergence from the published
  source of the `3 * E_v` rule and is flagged in the header so it is not "fixed"
  back.
- **Every count in range is scanned, never searched downward.** Under median
  centres admissibility is not monotone: merging two well-separated bands can
  place the merged median next to a neighbour's, so a coarser partition can fail
  where a finer one passed.

There is deliberately **no additive score** of the form
`w1 * GVF + w2 * (gap / gap_min)`. Goodness-of-variance-fit is bounded on [0, 1]
while the band gap scales as `range / k` and is unbounded, so over a large
vertical range the gap term dominates by orders of magnitude and the criterion
collapses to the smallest count regardless of the data. `unit_fluvial_surface_
ladder` includes an arithmetic arm demonstrating that the additive score changes
its argmax on data where the shipped rule does not.

### The parameter key policy

Splitting the CLIs renamed or retired more than 50 parameters. The parser's
typed parses iterate over a CLI's DEFAULTS MAPS and pull matching keys out of
the file, so a key present in the file but absent from every defaults map was
never looked at by anything — silently a no-op. An old parameter file saying
`find_floodplain: true` would have run happily having done no classification at
all, with no error and no warning.

`lsdtt_param_parser` now scans the RAW KEYS PRESENT IN THE FILE (never the
defaults maps, which would silently re-accept everything) against three maps:

| group | behaviour | for |
|---|---|---|
| removed | exact migration message, exit 1 | keys that TURNED AN ANALYSIS ON, where ignoring one leaves the user believing it ran |
| deprecated alias | value copied to the new name, warning names it | a parameter that still exists under a new name |
| ignored | specific message, run continues | a parameter that moved to another tool, or whose mechanism no longer exists |

Anything in none of the three and in no defaults map is reported as an unknown
key — a warning, never fatal, since a tool may legitimately read a key by
another route.

**That scan found three stale committed parameter files on its first run**, each
carrying a key no CLI declares and that had therefore been doing nothing:
`hydraulic_run_mode` and `raster_is_filled` in the two hydraulic templates, and
`n_timesteps` in the `lem_coupling_*` templates. It also explains two long-
standing red tests: `failure_hydraulic_routing_invalid_mode` asserts a message
`lsdtt-graphflood` cannot produce, because that CLI does not declare
`hydraulic_run_mode` at all.

### Testing regime

Label `fluvial_surfaces`, 23 tests, run as
`ctest -L '^fluvial_surfaces$' -LE '^heavy$'`. What counts as a genuine regression:

- tread purity on the dissected staircase falling below 1.0;
- any ladder window reporting `resolved` when its minimum band-centre gap is
  below `kappa * E_v`;
- the label -> band map becoming multi-valued;
- the `hand_fname` path diverging from the computed path;
- `find_population_breaks` output moving (T5) — it is shared, so that moves two
  tools at once.

Not regressions: absolute component counts on real DEMs, and timings.

Three checks are worth singling out because no single-file assertion can express
them, and each needed its own runner:

- `rigorous_fluvial_surfaces_label_band_identity` checks the label JOIN across
  outputs — that every label in the raster has exactly one patch row, that the
  label -> band map is single-valued and agrees with the patch table, and that
  every per-pixel profile row joins back and agrees with both rasters. A schema
  assertion cannot see a join.
- `smoke_fluvial_surfaces_hand_input` runs the tool twice and requires the
  supplied-HAND path to produce bit-identical bands and labels to the computed
  path. Both arms use a global ladder, and that is load-bearing: a raster
  carries no HAND donor nodes, so a supplied HAND cannot be windowed, and
  comparing a windowed run against a supplied one would fail for a reason that
  has nothing to do with composability.
- `smoke_fluvial_surfaces_pipeline` runs both CLIs and asserts the HANDOVER —
  that stage 2 reports the valley was supplied, the one line separating "read
  the mask" from "quietly derived its own".

### A process note that cost real time

During the split, an edit to `tests/CMakeLists.txt` that deleted a range between
two named anchors silently removed six test registrations that happened to sit
between them. Every run afterwards reported `100% tests passed`, because a
deleted test cannot fail. The only visible symptom was the compartment count
falling from 21 to 17, and that was not noticed for two phases.

**A green run is not evidence when the run is a different set of tests.** Compare
the test COUNT, and preferably the test SET, across any change that was not
meant to remove tests. `tests/CMakeLists.txt` now carries a guard that refuses to
configure when a test matches no domain rule and would be unreachable by a
compartment run; that guard caught a badly-named test immediately in a later
phase and is exactly the protection this failure needed.

---

## 10o. Whitebox as a second D8 reference, and a raster-writer defect (2026-09-08)

Wiring `accum_d8@whitebox` into the flow benchmark's asserted comparisons. The
wiring itself was three lines; everything below is what it exposed.

### Whitebox uses steepest descent - measured, not assumed

D8 accumulation was asserted against ONE independent implementation
(TopoToolbox). RichDEM is excluded because it selects the LOWEST-ELEVATION
neighbour rather than the steepest DESCENT. Whitebox had never been
characterised, so it carried timing rows only.

`tests/scripts/check_whitebox_d8_rule.py` settles it against two analytic
fixtures.

**A planar ramp cannot answer this question**, and assuming otherwise is the
trap. On a ramp the S, SE and SW neighbours are TIED in elevation, so a
lowest-elevation rule that breaks ties toward cardinals passes it - RichDEM does
exactly that (corrected at `test_flow_properties.cpp:996`). The ramp is still
run, because it calibrates the counting convention: Whitebox's
`out_type="cells"` was measured to INCLUDE the cell itself.

The second fixture removes the tie. One cell is given a diagonal neighbour that
is lower and a cardinal neighbour that is steeper:

| receiver | drop | distance | gradient | |
|---|---|---|---|---|
| cardinal S | 8 m | 30.00 m | 0.2667 | steeper |
| diagonal SE | 10 m | 42.43 m | 0.2357 | lower |

The script computes the whole expected accumulation raster twice - once per rule
- with its own reference D8, rather than hand-deriving what the surrounding
cells do. The two rules differ on 11 of 441 cells. Whitebox matched steepest
descent on every one, and a `d8_pointer` read (direction codes calibrated on the
ramp, not assumed) confirmed it independently. The script reports INCONCLUSIVE
if accumulation and pointer ever disagree.

**Note the headline cells did not decide it.** At the two candidate receivers
the steepest-descent expectation is S 24 / SE 24 - equal, because SE inherits
its own feeder column. The whole-raster comparison is what discriminates. A
check reading only the two neighbours would have concluded nothing while
appearing to.

### The raster writer set a TIFF predictor that does not apply to float data

`lsdtt_data_loader.cpp` (2 sites) and `lsdtt_raster_warp.cpp` (1) hardcoded
`PREDICTOR=2` on every tiled DEFLATE output. Predictor 2 is horizontal
differencing, which TIFF 6.0 defines for INTEGER data; float data takes
predictor 3. lsdtt3 writes Float32 almost everywhere. GDAL accepts the pairing
and reads it back without complaint, so nothing inside this project noticed.

Measured on alpine_1000 (1000x1000 Float32, 3.81 MB raw), tiled DEFLATE:

| predictor | size | |
|---|---|---|
| 2 | 2.68 MB | what shipped; spec-invalid on float |
| 3 | 2.19 MB | correct, and 18 % smaller |
| none | 3.14 MB | |

Not a trade-off: the conforming choice is also the smaller one. Fixed by
`tiff_predictor_for()` in `lsdtt_data_loader.hpp`, which selects from the band
type and is used at all three sites. The other DEFLATE writers
(`lsdtt_raster_fetch.cpp`, `lsdtt_raster_reproject.cpp`) set no predictor and
were never affected - which is why fetched FIXTURES always read correctly while
lsdtt3-written outputs did not.

**This defect is in every raster lsdtt3 has written to date.** Any consumer
using a strict TIFF reader has been getting garbage from lsdtt3 outputs
silently.

### What the predictor fix did NOT fix

The obvious inference - invalid predictor causes the bad read - is wrong, and
was tested rather than assumed. A creation-option matrix on identical data:

| creation options | Whitebox read |
|---|---|
| `COMPRESS=DEFLATE` | correct |
| `COMPRESS=DEFLATE, PREDICTOR=2` | all NaN |
| `COMPRESS=DEFLATE, PREDICTOR=3` | all NaN |
| `COMPRESS=DEFLATE, TILED=YES` (no predictor) | correct |
| uncompressed | correct |

**Whitebox Workflows cannot read a DEFLATE-predictor GeoTIFF at all**, valid
predictor or not; tiling is irrelevant. That is a limitation of its TIFF reader.
Going predictor-free globally would cost every user 43 % on every float raster
to accommodate one reader, so `whitebox_reference.py` transcodes a
predictor-free temporary instead. Pixel values are bit-identical; only the
compression of the bytes changes. Staging precedes the timing and read time is
reported separately, so it costs the measurement nothing.

Two implementation notes worth keeping: Whitebox's `Raster` holds a file handle
and the TOOLS re-read from disk, so the staged copy must outlive every tool
call, not just `read_raster()`; and `assert_readable()` now refuses to time an
all-NaN raster, because that is the failure mode here.

### How it presented, which is the reason this section exists

The first run reported `accum_d8@whitebox` with a maximum accumulation of **123
cells against lsdtt3's 41,437**. It read as a catastrophic routing
disagreement. Whitebox was routing a grid of NaN, faithfully, and reporting a
number. No error, no warning.

The diagnostic that broke it open was not the accumulation raster but the
POINTER raster: 53,155 of 65,536 cells with no receiver, on a filled DEM
measured to contain zero flat cells. That is not a routing difference, and it
is what said "this tool is not seeing the surface".

**A silent all-NaN read is indistinguishable from a result once it reaches a
summary table.** Every reference wrapper should assert its input loaded.

### Measured agreement, and the pinned tolerances

With the DEM actually loading:

| fixture | reference | rank corr | top 1 % channels | total |
|---|---|---:|---:|---:|
| alpine_small | topotoolbox | 1.0000000000 | 0.0 | 0.0 |
| alpine_small | **whitebox** | **1.0000000000** | **0.0** | **0.0** |
| alpine_small | richdem *(report-only)* | 0.5088829645 | 16.56 % | 4.592 % |
| lakeland_small | topotoolbox | 0.9401849288 | 15.11 % | 0.1466 % |
| lakeland_small | **whitebox** | **0.9996777158** | **0.60 %** | **0.0051 %** |
| lakeland_small | richdem *(report-only)* | 0.7927508258 | 6.58 % | 0.0991 % |

Three implementations agree cell-for-cell on high relief. Pinned just above
observed: alpine `0.9999 / 1e-6 / 1e-6`, lakeland `0.999 / 0.0065 / 6e-5`.

**Whitebox tracks lsdtt3 more closely than TopoToolbox does on LOW relief**
(0.60 % vs 15.11 % on top channels), which is the opposite of what the alpine
row suggests. All three use the same selection rule, so the spread is
tie-breaking among the near-ties that survive an epsilon fill on flat ground.
That is a coincidence of implementation, not a correctness result, and must not
be quoted as one.

### Whitebox as a value reference for fill (same session)

With the raster now loading, Whitebox was added to two more value comparisons.

`fill@whitebox`, routed through `fill_depressions_wang_and_liu` (`fix_flats=
False`), is **bit-exact on both fixtures** — 0.0 m, 0 differing cells, volume
ratio 1.0. lsdtt3's fill now reproduces cell-for-cell against three independent
implementations.

`fill_epsilon` needed a units correction first. **`fill_min_slope` is a
GRADIENT** — `lsdtt_raster_flow_preprocessing_cpu.cpp` imposes
`min_slope * distance`, so on a 30 m grid 1e-4 raises a cardinal neighbour by
3e-3 m and a diagonal by sqrt(2) more — whereas Whitebox's `flat_increment` is
metres per cell. Comparing them directly compares increments 30x apart. With
`--flat-increment` set to `fill_min_slope * cell_size`:

| fixture | reference | max abs | cells >1cm | volume ratio |
|---|---|---:|---:|---:|
| alpine | richdem | 0.1052 m | 0.523 % | 1.0023 |
| alpine | **whitebox** | **0.0246 m** | **0.105 %** | **1.00037** |
| lakeland | richdem | 0.4995 m | 91.29 % | **40.54** |
| lakeland | **whitebox** | **0.1459 m** | 81.36 % | **1.39** |

The lakeland volume ratio is the point: RichDEM's one-ULP scheme moves 40.5x the
material lsdtt3 does on flat-dominated terrain; a matched constant increment
moves 1.39x. Same fixture, same lsdtt3 output — **the 40x was the epsilon
scheme, not lsdtt3.** 1.39 is not forced to 1.0 and should not be: lsdtt3 scales
the gradient by distance and Whitebox does not, which is a real operational
difference.

One discrepancy remains open. Whitebox's Barnes `fill_depressions` disagrees
with lsdtt3 by up to 9.78 m on alpine — the same offset by which it disagrees
with Whitebox's *own* `wang_and_liu`, which is bit-exact. So the disagreement is
that tool, not the increment scheme. Until 2026-09-13 `fill@whitebox` and
`fill_epsilon@whitebox` in the timing campaign still called it; the pairings
are now corrected (see "What is compared"). **Why it disagrees with three other
implementations is unresolved.**

### Drainage area without the s-vector (same session)

`get_drainage_area_raster()` now falls back to `accumulate_raster()` when the
s-vector was not built, rather than throwing, and `lsdtt-drainage-analysis`
exposes `build_flow_s_vector` (default true, ignored with a warning when
connectivity or height-above-channel is requested). **P11** in
`test_flow_properties.cpp` asserts the two paths give a bit-identical raster —
exact, not a tolerance — on all four fixtures.

This exists so the benchmark can report a like-for-like D8 row against tools
that build no equivalent structure. The timing that justifies it is NOT recorded
here: it must be measured on a server, not on a thermally limited laptop.

### Chi agrees with TopoToolbox to sub-percent (2026-09-08)

The first value comparison in lsdtt3's ANALYSIS layer rather than its
preprocessing layer. Both implement the Perron & Royden (2013) integral
`chi = int (A0/A)^(m/n) ds`; lsdtt3 follows Mudd et al. (2014) with the same
definition.

**Method.** Same DEM, same channel-defining threshold, `A0 = 1 m^2` and
`m/n = 0.45` on both sides. lsdtt3's chi fit node table is joined to
TopoToolbox's chi grid on the FLAT INDEX, then regressed per basin.

Per basin rather than as a difference, because **chi is baselevel-relative**:
each tool seeds from its own outlet, so a raw difference is dominated by a
constant offset that says nothing about the integrand. The regression separates
the two: the SLOPE is the ratio of the integrands and must be 1, the INTERCEPT
is the baselevel offset and is free.

Measured on a synthetic island (400^2, 57 % land), 9 complete basins,
804 matched nodes:

| basin | n | slope | intercept | R^2 | max abs residual |
|---:|---:|---:|---:|---:|---:|
| 0 | 50 | 1.006378 | -2.57999 | 0.9999863 | 0.004928 |
| 1 | 167 | 1.003793 | -0.41878 | 0.9999826 | 0.011923 |
| 2 | 58 | 1.004604 | -3.30436 | 0.9999706 | 0.008226 |
| 3 | 25 | 1.002134 | -1.35847 | 0.9999025 | 0.008481 |
| 4 | 42 | 1.003935 | -3.09087 | 0.9999748 | 0.007555 |
| 5 | 101 | 1.004117 | -1.53699 | 0.9999673 | 0.009011 |
| 6 | 69 | 1.001952 | -1.94635 | 0.9999708 | 0.009110 |
| 7 | 231 | 1.003246 | -0.41685 | 0.9999959 | 0.011733 |
| 8 | 61 | 1.003783 | -1.11684 | 0.9999903 | 0.006219 |

Slope within **0.2-0.6 % of 1** on every basin, `R^2 >= 0.9999`, residuals below
0.012 chi units. Two independent implementations of the chi integral agree to
sub-percent on the same terrain.

**The residual slope is SYSTEMATIC, not noise** -- all nine basins between
1.0020 and 1.0064, all positive. Run down 2026-09-08, and it is entirely the
QUADRATURE RULE:

| | rule | increment over one segment of length L |
|---|---|---|
| lsdtt3 | left-endpoint rectangle | `L * f(A_upstream)` |
| TopoToolbox | trapezoidal (`streamquad_trapz`) | `L * (f(A_up) + f(A_down)) / 2` |

with `f(A) = (A0/A)^(m/n)`. Drainage area increases downstream, so `f`
decreases, so `f(A_up)` is always greater than the trapezoidal average.
**lsdtt3's increment is therefore larger on EVERY segment** -- which is why the
bias is systematic and one-signed rather than scatter.

Proved rather than argued. `chitransform` was reimplemented in numpy from
TopoToolbox's own `source`/`target`/`node_to_node_distance` arrays; the
trapezoidal version reproduces `chitransform` to **max difference 0.0**, which
establishes the edge handling is right. Swapping only the quadrature rule to
left-endpoint and re-running the same per-basin regression:

| regression | slope range | mean |
|---|---|---|
| lsdtt3 vs TopoToolbox (trapezoidal) | 1.00195 - 1.00638 | 1.00377 |
| lsdtt3 vs left-endpoint variant | **1.00000 - 1.00000** | **1.00000** |

Exact on all nine basins. So the two implementations agree EXACTLY once the
quadrature rule is matched: same integrand, same drainage area, same path
lengths (including the sqrt(2) on diagonals), same `A0`, same `m/n`. The 0.4 %
was never a disagreement about chi.

**RESOLVED: lsdtt3 adopted the trapezoidal rule (2026-09-08).** Trapezoidal is
second-order accurate and left-endpoint first-order, so the old rule carried an
O(h) bias that grows with cell size -- ~0.4 % at 30 m, worse on a coarse DEM.
The package is a development version with no released results, so continuity
with past output was not a reason to keep a known first-order error.

`chi_increment()` in `lsdtt_chi_extraction_cpu.cpp` now holds the rule in ONE
place, used by all three extraction paths (full-raster, outlet-seeded
zero-base, outlet-seeded inherited-base) -- previously the expression was
written out three times.

Verified after the change, same island fixture, same 9 basins, 804 nodes:

| | slope range | mean | R^2 | max abs residual |
|---|---|---|---|---|
| before (rectangle) | 1.00195 - 1.00638 | 1.00377 | >= 0.9999 | 0.012 |
| **after (trapezoid)** | **0.999999 - 1.000001** | **1.000000** | **1.00000000** | **5e-6** |

The residual is float32 rounding in the CSV, not a difference in the value.
**lsdtt3 and TopoToolbox now agree on chi to the precision the output format
carries.**

Direction of the change, on the committed `tiny` fixture: every chi value fell
(max 4.0297 -> 3.6297), and chi = 0 is preserved at the outlet. That is the
predicted sign -- the old rule evaluated the integrand at the upstream end of
each segment, where it is largest.

`rigorous_chi_nodes_sha256` caught this, which is exactly what a hash pin is
for. It was re-pinned only after the direction and the magnitude were checked
against TopoToolbox; the CMakeLists entry records why the hash moved.

**Every chi and k_sn produced before this date is ~0.4 % high at 30 m
resolution**, more at coarser cell sizes. Anything carried forward from an
earlier run should be regenerated rather than compared against new output.

**What this rules out.** A mismatched `A0` would show as a slope of
`A0_ratio^(m/n)` -- at `m/n = 0.45`, TopoToolbox's own default of `a0 = 1e6`
against lsdtt3's 1.0 would appear as ~631. Slope 1.004 confirms both sides ran
at `A0 = 1 m^2`, which is the condition for a chi gradient to be k_sn.

**The fixture requirement is itself a finding.** lsdtt3 drops any basin whose
upslope area is `is_influenced_by_invalid` -- a cell with ANY of its 8
neighbours nodata, or on an open boundary. On a cropped DEM every basin's
headwaters reach an edge, so all chi outputs are empty. TopoToolbox does not
decline, and returns chi for edge-truncated basins whose drainage area is wrong
by construction. The comparison therefore needs a DEM containing a whole
catchment; a synthetic island supplies one, but only a fat, smooth one
(`island_fill_fraction: 0.95`, 3 harmonics gave 57 % land; 0.55 with the
defaults gave 16 % in lobes, and no complete basins at all).

### Chi value tolerances, pinned (2026-09-08)

`tests/scripts/run_svector_benchmark.py --cases chi_value --check` now enforces
the comparison. Value metrics only: nothing in that driver asserts a timing,
which would be meaningless off a dedicated machine.

**The fixture is generated, not committed.** `lsdtt-synthetic-networks` is
deterministic -- the same `random_seed` gives a byte-identical DEM (verified,
two runs, matching sha256) -- so the island is built from a recipe held in the
driver and its sha256 is checked before any tolerance is applied. A mismatch
SKIPS rather than fails, because a different GDAL or compiler is an environment
difference, not a defect. Applying a pinned tolerance to a different raster is
how a benchmark comes to assert nothing.

Measured from the driver's own run (threshold 100, no explicit outlets),
5 basins / 783 nodes:

| metric | observed | pinned |
|---|---|---|
| slope | 0.999999666 - 1.000000146 | [0.999999, 1.000001] |
| R^2 | 1.0000000000 | >= 0.999999999 |
| max abs residual | 5.428e-06 | <= 1.0e-05 chi units |
| basins | 5 | >= 5 (floor) |
| matched nodes | 783 | >= 780 (floor) |

The residual is float32 rounding in lsdtt3's CSV, not a difference in value.
Basin and node counts are FLOORS: a change that REDUCES the comparison is the
failure worth catching, and more basins is not a problem.

**A pinning mistake worth recording.** The first pin used 9 basins / 804 nodes,
taken from an earlier hand-run that supplied explicit outlets and
`n_largest_basins`. It failed immediately against the driver's own output.
Pin what the thing being checked actually produces, not what a
differently-configured run produced.

**The pin was verified to fail.** Running both sides at `a0 = 1e6` instead of
1.0 exits 1. Note WHICH bound caught it: `--a0` moves both sides together, so
the slope stayed 1.000000 and the RESIDUAL fired (5.45e-3 against the 1e-5
limit -- chi values grow ~170x and float32 precision degrades with them). The
slope bound guards the other failure mode independently: the pre-trapezoid
left-endpoint rule measured 1.00195-1.00638, far outside the pinned window.

### The CLI-vs-in-process asymmetry, and what it was hiding (2026-09-08)

The s-vector basin case compared an lsdtt3 CLI against reference routines timed
in-process. Three separate problems, each hiding the next.

**1. `lsdtt-basin-extraction` printed no `Compute time:` lines at all.** The
driver fell back to whole-process wall clock -- and labelled it
`cli_reported_compute` in the CSV anyway, which is the one thing the
timing-scope discipline exists to prevent. Fixed on both sides: the CLI now
reports its phases, and the driver reports `whole_process` honestly when it
finds none.

The scale of what that hid: flow routing measures **0.029 s** on a 256^2
raster inside a **~4.4 s process**. Over 99 % of the "lsdtt3 basin time" was
process startup, GDAL initialisation and file I/O.

**2. Not all phases are comparable, so they are timed separately.**

| phase | reference equivalent | in the comparison? |
|---|---|---|
| flow routing | all references build a flow structure | YES |
| junction network | **nothing else here builds one** | NO -- reported separately |
| basin selection | `dependencemap` / `watershed` | YES |

`seconds` is now flow + basin selection. Folding the junction network in would
charge lsdtt3 for a structure no reference builds -- the same error `accum_d8`
made with the s-vector.

**3. A 2.1 s phase that was not basin work.** The first placement of the
basin-selection timer ran past the selection into the outlet-file read and the
coordinate diagnostics, giving 2.16 s against 0.03 s of flow routing. Narrowed
to the selection proper, the phase is **0.0002-0.0006 s**.

**That leaves a real finding: the basin outlet coordinate diagnostics cost
~2.1 s on a 256^2 raster with 9 outlets.** That is reporting, not analysis, and
it is two orders of magnitude above the routing it follows. The shape matches
the PROJ defect already recorded for the cosmogenic inversion (a transform
rebuilt per point rather than reused, 73.7 us against 0.37 us). NOT confirmed
here, and not fixed -- but it should be looked at, because it is paid by every
`lsdtt-basin-extraction` run that supplies an outlet file.

**Measured after all three fixes** (256^2, depths 1/4/7 -- shape only, this is a
thermally limited laptop and none of these are results):

| engine | depth 1 | depth 4 | depth 7 |
|---|---|---|---|
| lsdtt3 total (flow + selection) | 0.018 | 0.037 | 0.015 |
| ... of which basin selection | **0.00041** | **0.00020** | **0.00057** |
| topotoolbox | 0.009 | 0.011 | 0.010 |
| whitebox | 0.316 | 0.763 | **1.344** |

lsdtt3's selection phase is FLAT with depth and sub-millisecond -- the s-vector
sub-interval read -- while Whitebox rises roughly linearly, which is the
per-nesting-level cost its `unnest_basins` design makes explicit. That is the
shape the case was built to expose. The magnitudes need a server.

### RichDEM was excluded from the surface comparison on a false premise (2026-09-08)

The case registry's comment gave "RichDEM has no slope tool" as its example of
why an engine may be absent from a case. **It has one.** RichDEM 2.3.0's
`rd.TerrainAttribute` offers `slope_riserun`, `slope_percentage`,
`slope_degrees`, `slope_radians`, `planform_curvature` and
`profile_curvature` — most of what `slope_only` and `polyfit_all` compare.

So the surface-metric comparison has been running with **three engines where
four were available**, because a code comment asserted a capability limit that
was never checked against the package. Found while auditing why TopoToolbox had
no D-infinity or MFD rows.

The TopoToolbox absences are real **but they are about the PYTHON BINDING, not
about TopoToolbox**, and getting that wrong in a paper would be caught
instantly by anyone who uses the MATLAB version:

| case | engine | why |
|---|---|---|
| `accum_dinf`, `accum_mfd` | topotoolbox | not bound in the Python package 0.0.12 — `FlowObject` takes `method='d8'` only. **TopoToolbox 3 implements both in MATLAB** (`md/usersguide_3_mfd.md`) |
| `fill_epsilon` | topotoolbox | not bound in 0.0.12 — `fillsinks(bc, hybrid)` takes no flat-increment argument |

So the correct statement is **"D-infinity and MFD are not reachable through the
TopoToolbox Python package at 0.0.12"**, not "TopoToolbox does not implement
them". An earlier draft of this section said the latter. It was wrong.

**This scopes the whole comparison and the paper must say so.** Three of the
four engines are driven through Python bindings and the fourth (lsdtt3) through
its CLI. The TopoToolbox binding in particular is at **0.0.12** — early, and
exposing a fraction of TopoToolbox 3. Coverage and performance measured here
characterise *those interfaces*, which for TopoToolbox is not the interface most
of its users have.

`SUMMARY.md` now distinguishes the two cases: `n/a` means the package cannot do
it, `-` means no comparable row was produced. They had been the same mark. The
engine also keeps its column when every one of its cells is `n/a`, which it
previously lost — dropping the column made "cannot" look like "never
considered".

**Resolved 2026-09-09.** `richdem_reference.py` gained a surface section
(`slope`, `aspect`, `curvature_profc`, `curvature_planc`, `curvature_all2`) and
RichDEM is wired into `slope_only` (`slope`) and `polyfit_all`
(`curvature_all2`). Both cases were already `report_only` and remain so.

The comparability question that made this more than a one-line addition has an
answer, and it is not "the values now agree". The four engines compute four
different things:

| engine | method | window |
|---|---|---|
| lsdtt3 | 6-term quadratic, **least-squares** fit, every metric from ONE fit | `polyfit_lambda_metres`, user-chosen |
| RichDEM | 9-term quadratic, **exact analytic** fit (Zevenbergen & Thorne 1987), differentiated; refits per attribute | fixed 3×3 |
| TopoToolbox | scipy convolutions, Evans 3×3 kernels; recomputes all five per curvature type | fixed 3×3 |
| Whitebox | 3×3 finite-difference operators in Rust; recomputes per curvature type | fixed 3×3 |

RichDEM is the **nearest of the three references to lsdtt3 in form** — both fit
a quadratic and differentiate it — and is still not value-comparable, because
an exact interpolation of a fixed 3×3 is not a least-squares fit over a chosen
window. Curvature sign conventions and aspect handedness were **not** checked
across the four and are not assumed to agree. No surface value is compared
anywhere in the suite.

The metric counts also differ within `polyfit_all`: lsdtt3 slope + aspect + 2
curvatures, TopoToolbox 4 curvatures, RichDEM and Whitebox 2 curvatures. The
case note says this. It is a cost-of-a-typical-request row, not a
same-work-per-engine row.

Two genuine RichDEM absences were recorded in `NOT_IMPLEMENTED` rather than
left as silent gaps, since an unexplained absence reads as "we did not bother":
`hillshade` (no hillshade attribute in 2.3.0) and `polyfit_many` (no
eight-metric call — recorded for all three references, as that row existed to
show lsdtt3's shared fit amortising; the case was removed 2026-09-13).

First measurement, alpine_1000, 1 thread, 3 repeats, **on the thermally limited
laptop — not quotable**, recorded only to show the rows are produced:

| case | lsdtt3 | richdem | topotoolbox | whitebox |
|---|---|---|---|---|
| `slope_only` (ms) | 59.8 | 47.8 | 66.3 | 55.9 |
| `polyfit_all` (ms) | 163.9 | 49.0 | 218.4 | 78.9 |

### TopoToolbox surface rows reported a zero standard deviation (2026-09-09)

`run_surface` in `topotoolbox_reference.py` returned `seconds`,
`seconds_median` and `seconds_all` but **not** `seconds_stdev`.
`run_benchmark_suite.py` defaults a missing field to 0.0, so every TopoToolbox
surface row in every campaign printed `+/- 0.00` — a spread that would have
looked like an implausibly stable measurement rather than a missing one. The
flow path in the same file emitted the field; only the surface path did not.
Fixed by emitting it explicitly. Verified: `--routine evansslope --repeats 3`
now reports a stdev of 7.8e-4 s where it previously reported 0.

### Consequence for recorded campaigns

`run_benchmark_suite.py`'s `accum_d8` case is `needs_filled`, so its Whitebox
rows were reading an lsdtt3-written filled DEM and were also routing NaN. The
fix is shared (both drivers call `whitebox_reference.py`), but **any campaign
recorded before 2026-09-08 is invalid for Whitebox flow rows and must be
re-run.**

---

## 10p. Hydrology test coverage — PLACEHOLDER, not yet written up

**Status: placeholder.** The tests below exist and pass. What they measure has
not been written up here, so this section is a stub naming them rather than a
record of findings. Treat the hydrology modules as **tested but not
documented**: the assertions are in the test sources, not in this file.

Context for whoever fills this in: the raster TOPMODEL stack was ported so a
dissertation student could use it, and it has **not been calibrated against any
observed catchment**. The tests are implementation checks — mass conservation,
closed-form limits, convergence — not evidence of predictive skill. Do not let
this section imply otherwise when it is written.

| Suite | Module | What it covers (to be confirmed and expanded) |
|---|---|---|
| `unit_topmodel_column` | `lsdtt_topmodel_cpu` | Column equations after TOPMODEL 95.02; exponential baseflow recession as a closed-form limit |
| `unit_topmodel_dynamic` | `lsdtt_topmodel_dynamic_cpu` | Prognostic saturated zone; the index relation should appear as its steady state |
| `unit_infiltration_excess` | `lsdtt_infiltration_excess_cpu` | Green-Ampt with exponentially declining conductivity; Green-Ampt infiltration capacity as a closed form. **Also needs:** whether `use_dimensional_ponding_test` is exercised in both settings, since the default reproduces the original FORTRAN's dimensionally inconsistent ponding test |
| `unit_canopy_interception` | `lsdtt_canopy_interception_cpu` | Rutter-type store; `C = k_S·LAI` and `p = exp(−κ·LAI)`. **Also needs:** a check that `PET_residual` divides by `evap_factor` rather than subtracting `E_canopy` directly — the header flags this as the easiest thing to get wrong |
| `unit_hydro_roughness` | `lsdtt_hydro_roughness_cpu` | The three-regime partition; `rigorous_roughness_partition` asserting that a land-cover change leaves channel `n` untouched; `jarrett` and `ferguson_vpe` depth dependence |
| `unit_hydro_routing` | `lsdtt_hydro_routing_cpu` | Linear reservoir, Muskingum-Cunge, kinematic wave |
| `unit_hydro_velocity_field` | `lsdtt_hydro_velocity_field_cpu` | Velocity, travel time, width-function IUH |
| `unit_giuh_discharge` | `lsdtt_giuh_cpu` | Width function; `v_i = v_ref (S_i/S_ref)^β` with β = 0.5 |

**Not covered anywhere, and worth deciding about:**

- **Graphflood has no entry here** beyond the DAGGER parity reference in §4.
  Its solver is stationary, explicit (forward Euler plus a flux limiter), and
  rebuilds topology from `z + h`; none of those three properties is pinned by a
  test recorded in this document.
- **No mass-conservation figure is quoted anywhere** for the distributed
  hydrology path. The design doc claims mass conservation is checked; this
  document does not say to what tolerance.
- **No cross-implementation reference.** Unlike flow routing (TopoToolbox,
  RichDEM, Whitebox) and cosmogenics (CRONUS, OCTOPUS), the TOPMODEL stack has
  no external oracle. (The S-HYPE column does: HYPE itself, §10q.) Whether one
  is worth building is an open question given
  the intended use.

---

## 10q. The S-HYPE column (`column_model: hype`, 2026-10-04)

`lsdtt-distributed-hydrology` with `column_model: hype` reimplements the land and river processes of HYPE 5.36.0 (SMHI; Lindström et al. 2010), with HYPE's default options. Unlike the TOPMODEL stack in §10p, it has an **external oracle**: HYPE itself. The equations are in `docs/design_and_algorithms/hydrological_modelling.md`. As with TOPMODEL, nothing here is evidence of predictive skill: **the model has not been calibrated against any observed catchment**.

### 10q.1 The oracle, and why CTest never touches Fortran

`tests/expected/hype/` holds output written by HYPE 5.36.0 itself, compiled from the release source (SHA-256 and git commit in `tests/expected/hype/README.md`).

- `tests/scripts/generate_hype_oracle.py` downloads, builds and runs HYPE. It is run by hand, never by CTest. `--check` regenerates the files and confirms they are byte-identical. It stops if a newer HYPE release appears, until rerun with `--accept-new-version`.
- CTest reads only the committed files. No test downloads, compiles or links third-party code.
- There are seven cases. Six each isolate one process with the others switched off: `snow`, `et`, `diversion`, `percolation`, `soil_runoff`, `routing`. The seventh, `default_all`, runs everything together over two sub-basins in series. Each case is a short daily run with a compressed seasonal cycle. Every case is defined once, in `cases/<case>.json`, and both the HYPE input files and the lsdtt3 tables are generated from that definition.
- SourceForge serves an HTML page instead of the archive to a browser user-agent. The script therefore sends a tool user-agent and checks for the gzip magic number.

### 10q.2 The tolerance has a measured basis

The released HYPE build is single precision; lsdtt3 runs in double. Each case was also run through a double-precision HYPE build (`-fdefault-real-8`), and the largest single-vs-double difference per variable is recorded in the oracle README. That difference is the **floor**: below it, a double-precision reimplementation cannot be judged against the single-precision oracle. The largest floor over all cases is **1.08e-4 mm** (layer-3 soil water in `default_all`, which holds about 440 mm).

The column, routing and model tests pass a value when |lsdtt − HYPE| ≤ 5e-4 + 2e-5·|HYPE|, about five times the floor, in mm, or in m³ s⁻¹ for flows. Each run prints the largest difference actually seen for every variable, so drift is visible before it fails. The largest measured differences are 1.08e-4 mm for the column, equal to the floor, and 2.8e-6 m³ s⁻¹ for sub-basin outflow.

### 10q.3 What each suite asserts

| Suite | Module | Checks | What it asserts |
|---|---|---|---|
| `unit_distributed_hydrology_forcing` | `lsdtt_distributed_hydrology_forcing` | 18 | Rain/PET loaders (moved out of the CLI unchanged), and the signed temperature loader |
| `unit_hype_parameters` | `lsdtt_hype_parameters` | 95 | Table readers: an absent parameter is 0, unknown names are errors. Validation: negative values, fractions out of [0, 1], `wcwp+wcfc+wcep > 1` and non-finite GeoData are rejected. Downstream ordering. The oracle cases' tables load as written |
| `unit_hype_column` | `lsdtt_hype_column_cpu` | 102 | Every class variable of the six process cases and `default_all` against HYPE within the tolerance above. Per-unit water balance on every step. The rain/snow ramp in its three regimes. Percolation never fills a receiving layer beyond its pore volume |
| `unit_hype_routing` | `lsdtt_hype_routing_cpu` | 41 | Translation shift and volume, and the damping box against the exact linear-reservoir solution. `cloc`/`cout` of the `routing` case against HYPE. The ring-buffer queue is bitwise equal to a reference copy of the original shifting queue |
| `unit_hype_model` | `lsdtt_hype_model_cpu` | 95 | `default_all` and `routing` through the assembled model against HYPE. Domain balance closure (about 3e-16 relative). Turning off balance tracking or per-unit flux storage changes no flow or state, bitwise. A unit's exception thrown inside the parallel land loop surfaces as the lowest-index failure at any thread count |
| `unit_subbasin_partition` | `lsdtt_subbasin_partition_cpu` | 113 | Outlet selection (user points, Strahler order and order range, with the base-level edge-strip merge). The incremental partition, GeoData statistics and class fractions. Flow-distance river lengths, all on synthetic DEMs whose answers are known by hand |
| `unit_hype_raster_setup` | `lsdtt_hype_raster_setup_cpu` | 18 | Class raster from land-use × soil, set-up from rasters, and the HYPE-format `GeoData.txt`/`GeoClass.txt` writer |
| `unit_hype_thread_invariance` | the model | 5 | 2,000 sub-basins × 8 classes × 365 steps: sub-basin runoff and local and total outflow on every step, the final unit states and the balance terms are **bitwise identical** at 1 and 8 threads |

**CLI tests** (`smoke_`/`failure_distributed_hydrology_hype_*`) run on a generated fixture (`fixture_hype_cli`), never on committed rasters. The fixture is a 60×60, 10 m DEM with a Y-shaped network, two land uses × two soils and 120 days of forcing with snow.

- **Smoke tests:**
  - all three outlet modes;
  - both snapping methods that move a point, with gauges placed three columns off the channel: `expanding_square` lands at (16, 21) and `downstream` at (18, 22). Each lands in a different place, so the sub-basin areas differ (45,300 against 54,200 m²);
  - the daily default for `dt_seconds`;
  - a run shorter than the forcing;
  - nodata in the soil raster.
- **The nodata case:** a 50-pixel block of nodata soil straddles a class boundary. The test pins the renormalised fractions in `GeoData.txt`, and these were checked against an independent pixel count of the rasters. The pixels stay in the sub-basin area but leave the class fractions.
- **Failure tests** (each asserts the message, not only the exit code):
  - missing temperature;
  - non-daily step;
  - `fsceff ≠ 0`;
  - unknown parameter;
  - missing class pair;
  - bad `column_model`;
  - no outlet mode;
  - snapping `none` on an off-channel gauge;
  - a sub-basin with no class pixel;
  - forcing series of different lengths;
  - `n_timesteps` past the end of the forcing;
  - an impossible and a malformed `hype_start_date`.

When the HYPE path was added to the CLI, the TOPMODEL outputs were compared before and after for four configurations. Discharge, balance, IUH and rasters were byte-identical. That was a one-off check, not a registered test.

### 10q.4 Speed: what was changed, and the measured effect

Every speed-up kept the tests above unchanged, and thread invariance had to hold bitwise. The timing harness is `lsdtt-test-hype-performance` (always built; registered as `performance_hype_model` under `LSDTT_ENABLE_PERF_TESTS`; synthetic workloads small 100×5, medium 2,000×8 and large 10,000×10 sub-basins×classes). It reports ns per unit-step.

| Change | Effect (medians, interleaved runs) |
|---|---|
| Balance tracking made optional | −35 % at medium, −13 % at large (1 thread). At large it had been 23 % of step time |
| OpenMP over units | 1.35× at 4 threads (large); no gain from 4 to 8 threads |
| Ring-buffer translation queue | Bitwise equal; not measurable (routing is about 3 % of step time) |
| PET seasonal factor computed once per step; units slimmed from 664 to 360 bytes; per-unit flux records optional | −48 % at medium (1 thread) |
| All of the above, large workload, 4 threads | **157 → 48 ns per unit-step** against the serial baseline |

Scaling stops at about 4 threads, and cost per unit-step rises with unit count (91 → 130 → 157 ns serial). Both point at memory traffic rather than compute, which is why slimming the unit paid more than threads did.

!!! warning "Laptop timings"

    These numbers come from one thermally and power-limited laptop
    (Intel i5-1145G7, 4 cores / 8 threads, WSL2, g++ 11.4, `-O3`), where
    run-to-run spread reached about 1.5×. Only differences within an
    interleaved pair of runs are meaningful. Treat them as direction and rough
    size, and re-measure on a server before quoting a speed-up. For threaded
    runs, the harness's phase split (land / aggregation / routing / balance)
    is not valid; only the total per unit-step is.

### 10q.5 Not tested

- **No real catchment:** there is no comparison against an observed hydrograph.
- **Not run in HYPE:** no lsdtt3-derived set-up has been run in HYPE, although the GeoData writer makes this possible.
- **Disabled options:** the HYPE options switched off in this version (fractional snow cover, sublimation, tile drainage, frozen soil, lakes) are rejected at start-up rather than tested.

---

## 10r. Flat resolution: the Barnes et al. (2014b) options

The method is described in `docs/design_and_algorithms/flow_routing.md`
(section "Flat resolution: Barnes et al. (2014b)"). This section records how it
was tested, what each oracle asserts now, and what the tests found. The option
is `fill_flat_resolution_choice`: `inline_gradient` (the default, unchanged),
`barnes2014b_gradient` and `barnes2014b_float_step`.

### 10r.1 Where the code is

- `src/cpu/lsdtt_raster_flat_resolution_cpu.{hpp,cpp}`: the Barnes mask, the
  float32 step, both resolvers, the rounds and the final pass, and the
  parse / name functions.
- `src/cpu/lsdtt_raster_flow_preprocessing_cpu.{hpp,cpp}`: the `fill_pits`
  overloads that take a `flat_resolution_choice`. `inline_gradient` forwards to
  the existing fill with the same arguments.
- `src/cpu/lsdtt_fill_flat_resolution_reporting_cpu.{hpp,cpp}`: the one helper
  every CLI calls for the choice check, the report lines, the WARNING and the
  resolve-only path used when `raster_is_filled` is true.
- `src/cpu/lsdtt_flowgraph_cpu.cpp`: the MFD underflow fallback (10r.7).
- `src/core/lsdtt_common_parameters.cpp`: the parameter declaration and its help
  text, which the reference pages are generated from.

### 10r.2 What RichDEM does, measured first

Before any code, RichDEM's `ResolveFlats` was run on hand-built float32 fixtures
(richdem 2.3.0, pixi `benchmark` environment) and the raise of every cell was
measured in ULPs. The implementation follows these findings:

- Interior low-edge cells are raised by 2 (the mask's low-edge value).
- Raster-border cells are never modified, even where the mask would raise them.
  The towards pass still starts from border low edges, so the interior values
  are unchanged.
- Cells next to interior nodata count as draining and are raised.
- A float32 raster is not upcast; the output is float32.
- RichDEM's "Cells inappropriately raised above surrounding terrain" counts
  ordered neighbour **pairs**, not cells. The comparison is `>=`: a valley with
  one cell raised two float steps gives 8 pairs under `>=` and 7 under `>`. The
  C++ driver pins the 8.
- On a floor at 0 m a float32 step is a denormal (1.4e-45 m). The maximum raise
  on that fixture was 6.7e-44 m, 48 denormal steps. This is why the near-zero
  guard exists.
- RichDEM's counts of NO_FLOW, low-edge and high-edge cells equal the mask's.

### 10r.3 Test inventory

| test | labels | what it covers |
|---|---|---|
| `unit_flat_resolution` | `unit;flow` | `tests/cpp/test_flat_resolution.cpp`, in-process fixtures only (groups below) |
| `unit_flat_resolution_real_terrain` | `unit;flow` | the same driver in `--real-terrain` mode on the committed Copernicus crop (30 checks) |
| `smoke_fill_flat_resolution_float_step`, `smoke_fill_flat_resolution_gradient` | `smoke;flow;cli` | C1, C2: both options run and write. Plumbing only: `tiny_dem.asc` has no flats after a `min_slope: 0` fill |
| `failure_fill_flat_resolution_unknown_choice`, `..._gradient_zero_slope`, `..._fill_spill_merge` | `failure;cli;flow` | C3, C4, C5: exact refusal messages |
| `smoke_chi_analysis_flat_resolution_gradient`, `smoke_channel_extraction_flat_resolution_gradient`, `smoke_drainage_analysis_flat_resolution_float_step` | `smoke;cli;flow` | C6: the option reaches a routing CLI. Each asserts a line only the barnes path prints |
| `smoke_chi_analysis_flat_resolution_filled_float_step`, `failure_chi_analysis_flat_resolution_filled_inline`, `failure_chi_analysis_flat_resolution_unfilled_input` | `smoke` / `failure;cli;flow` | C7: `raster_is_filled` true runs the resolve-only path. The inline control must fail with "could not assign a receiver". An unfilled input is refused with the no-outlet message |
| `rigorous_fill_flat_resolution_valley_float_step`, `rigorous_fill_flat_resolution_valley_gradient`, `rigorous_fill_flat_resolution_knob_rounds` | `rigorous;cli;flow` | exact fill-difference ranges and stdout pins on the valley fixtures, including the round count on the one-ULP knob |
| `smoke_basin_extraction_*`, `failure_basin_extraction_*`, `smoke_basin_stats_*`, `failure_basin_stats_*`, `smoke_basin_surface_metrics_*`, `smoke_channel_network_analysis_*`, `smoke_cross_sections_*`, `smoke_bankfull_width_*`, `smoke_distributed_hydrology_*`, `smoke_fluvial_surfaces_*`, `failure_fluvial_surfaces_*`, `smoke_ridgeline_*`, `smoke_valley_centreline_*` (all `flat_resolution`) | `smoke;cli;flow` or `failure;cli;flow` | one barnes run in each remaining wired CLI, asserting the "Flat resolution:" block and a value the default fill does not give. The `fill_before_routing: false` NOTE tests cover the two CLIs that skip the fill |
| `rigorous_drainage_analysis_flat_resolution_real_terrain_float_step`, `..._gradient`, `smoke_drainage_analysis_flat_resolution_real_terrain_inline_gradient` | `rigorous;cli;flow` (control: `smoke`) | the Appendix B command on the committed crop (10r.4) |
| `smoke_flat_resolution_real_terrain_routing` | `smoke;flow`, network tier | `tests/scripts/check_flat_resolution_real_terrain_routing.py` (10r.6) |
| `benchmark_flat_resolution_richdem` | `benchmark;flow`, `tier_pixi` | oracle O1 (10r.5) |
| `benchmark_flats_float_step_alpine_small`, `benchmark_flats_float_step_lakeland_small` | `benchmark;flow`, `tier_pixi` | oracle O2 (10r.5) |

The option is wired into lsdtt-raster-preprocessing and into 13 routing CLIs:
chi-analysis, channel-extraction, basin-extraction, drainage-analysis,
channel-network-analysis, basin-stats, basin-surface-metrics, cross-sections,
bankfull-width, distributed-hydrology, fluvial-surfaces, ridgeline-analysis and
valley-centreline. Deliberately not wired: lsdtt-lem and divide-migration
(periodic boundaries, an evolving surface), the three cosmogenic CLIs (locked
numerics with sha256-pinned sentinels), graphflood (its own internal fill),
internal-drainage (`fill_pits_with_sinks`), stream-burning (several internal
fills), and feature-finder and mountain-front (a hard-coded 1e-4).

The driver `unit_flat_resolution` groups its checks as follows. No group F28
exists; the numbering skips it.

- **F1** `step_up_float32` is bitwise equal to a looped `nextafter`, including
  signed zero, denormals, the crossing from negative to positive, and the
  overflow and NaN behaviour.
- **F2** The valley mask equals the frozen RichDEM raise, element by element.
  Fixtures 2 to 4 agree with RichDEM in the interior.
- **F3** Descent: on the valley, the terraces and a filled 256 × 256 fractal,
  every NO_FLOW cell has a lower same-label neighbour.
- **F4** Routability: `build_d8` succeeds, `validate()` passes, and an
  independent check finds a strictly lower neighbour for every non-exit cell.
- **F5** The valley profile at column 12, rows 1-9, is pinned to
  2 3 4 5 104 5 4 3 2. The default fill gives 24 12 12 12 12 12 12 12 24 at the
  same cells, and the test asserts that its ratio stays within 2.
- **F6** Gradient exact: every D8 step on the valley has slope at least
  `fill_min_slope` in double precision, with no tolerance.
- **F7, F7b** The rim cap: a 60 × 60 flat under a rim 0.01 m higher is capped and
  routable. A rim 3 ULPs up has a quantised cap of 0, and the result equals the
  float-step result bit for bit.
- **F8** The big flat with its rim 1 ULP above it: 701 pairs and 232 cells, as
  RichDEM.
- **F9** Near zero: with the guard, the smallest float32 D8 slope is at least
  2^-23 / (√2 L). Without it some slopes are exactly zero.
- **F10** Idempotence. **F11a, F11b** The exact messages for an unfilled input
  and for a periodic boundary. **F12** One thread and four threads give the same
  bits, including every new field. **F13** An interior nodata hole.
  **F14** `fill_pits` with `inline_gradient` equals the existing fill on the
  valley, the unfilled fractal and two fixtures with interior nodata, through
  the overloads, thread counts and queue backends. **F15** The MFD underflow fallback (10r.7). **F16** Border
  cells are never raised. **F17** Edge cases: tiny rasters, all nodata, one
  elevation everywhere, NaN nodata, and validity taken from the mask array.
  **F18** The guard's bound, method dispatch and the no-outlet message.
- **F19** Quantisation: 1 m cells at 4000 m (u below the ULP), and a raise across
  1024 m.
- **F20** A ramp whose floor is an odd multiple of q/2. Half-ULP ties leave cells
  without descent, so verification fails and the region falls back to float
  steps. Even-multiple controls verify.
- **F21** A cell or a rim a few ULPs above a flat (category a: stranded rim
  cells). **F22** Two flats a few ULPs apart (category b: stranded low edges).
  **F23** Cascades of 3, 5 and 8 terraces: no cell is ever lowered, the stranded
  front moves up each round, and the final pass changes no cell on a strictly
  descending path to an exit. **F24** The round-cap and option messages.
  **F25** The final pass keeps the near-zero guard. **F26** Gradient diagnostics
  are recounted from the per-cell flag, with an explicit allowance for cells
  whose quantised unit equals the float step. **F27** Rounds from 2 on raise
  only NO_FLOW cells. **F29** `max_rounds = 1` with the final pass off equals
  frozen FNV-1a hashes of the pre-loop output, and the defaults equal it on
  every fixture that strands nothing.

The real-terrain mode checks the committed crop: the single pass strands the
pinned cells, and the router refuses them; the defaults take the pinned rounds
and the final pass and route; the `fill_pits` dispatch equals the resolver; the
round-1 fields and the gradient recount agree; and the `inline_gradient` control
routes.

### 10r.4 The real-terrain gate, and how it was shown to catch the bug

The failure was found on real terrain, not on any synthetic fixture. One Barnes
pass left between 20 and 2475 cells without a lower neighbour on each of the
four 500 and 1000 m benchmark rasters (alpine and lakeland), for both options,
and `lsdtt-drainage-analysis` exited 1 with "could not assign a receiver". The
synthetic fixtures had rims metres above their flats, and the RichDEM comparison
could not see the failure, because RichDEM strands the same cells. The default
`inline_gradient` routed the same DEMs.

The committed crop, `tests/data/flat_resolution_lakeland_crop_dem.tif`, is a
128 × 128 float32 window (37,947 bytes, UTM 15N, Copernicus GLO-30 attribution in
`tests/README.md`). On it the single pass strands 117 cells with the float step
and 115 with the gradient, and the pre-loop CLI exits 1 at row 5, column 123.

The gate was tested by mutation, in a scratch copy with the round cap set to 1
and the final pass off. Then both rigorous CLI tests fail with exit 1 and
"could not assign a receiver", and `unit_flat_resolution_real_terrain` fails two
checks. The `inline_gradient` control passes. The gate is therefore sensitive to
the failure it was added for.

### 10r.5 Oracles O1 and O2: what they assert now

Both compare with RichDEM 2.3.0 on float32 bits. The rules changed when the
loop was added, because RichDEM has no loop: its output is the first round only.

**O1, `check_flat_resolution_vs_richdem.py`** runs on fixtures it builds (valley,
open flat edge, interior nodata, two terraces, a mesa, a big flat with a rim one
ULP up, a 2-ULP bump, floors at +0, -0, -50 and +4000 m, two 256 × 256 fractals,
and the knob, terrace and cascade fixtures used by F21-F23). It also runs the
committed crop.

- The number of differing unguarded cells, minus the number of cells that
  lsdtt3 reports changing after round 1, is exactly 0.
- Every differing cell is strictly higher in lsdtt3 than in RichDEM (ordered
  float32 keys). No lsdtt3 cell is below RichDEM.
- Cells in near-zero-guarded flats are exempt, and a fixture that mixes guarded
  and loop-changed cells fails: no fixture does.
- Cells without descent: RichDEM's count equals lsdtt3's round-1 count, and
  lsdtt3's final count is 0. On the stranding fixtures RichDEM's count is above
  0, which documents that it has the same failure. On the others both are 0.
- Pairs: lsdtt3's round-1 pair and cell counts equal the recount from the
  RichDEM-comparable output. RichDEM's printed pair count must lie within the
  range its loop allows. Its loop raises cells in place, so where two raised
  flats touch the printed count depends on processing order: on the crop it
  prints 367, while a recount from its own output gives 364.

On the fixtures without stranding, the result is one round, zero changed cells
and zero cells without descent; equality is exact.

**O2, `flats_float_step` in `run_flow_benchmark.py`** runs on the fetched
alpine and lakeland small DEMs, against `FillDepressions(epsilon=False)`
followed by `ResolveFlats`. Its tolerances are: unexplained differences 0 (the
differing cells must equal the cells changed after round 1), cells below
RichDEM 0, and lsdtt3's cells without descent 0. The count of differing
unguarded cells is no longer pinned, since both fixtures strand in round 1. Its
measured values, from the same runs as the tolerances:

| fixture | cells compared | flat cells | differing (all changed after round 1) | cells without descent (RichDEM / lsdtt3) |
|---|---|---|---|---|
| alpine_small | 65,536 | 723 | 1 | 1 / 0 |
| lakeland_small | 65,536 | 61,913 | 51 | 41 / 0 |

Negative checks, run in scratch wrappers that perturb one cell, were caught:
a 1-ULP change to an unchanged cell gives `n_unexplained_differences = 1`, and
lowering a loop-changed cell below RichDEM gives `n_below_richdem = 1`.

### 10r.6 The network-tier routing script

`tests/scripts/check_flat_resolution_real_terrain_routing.py` is stdlib-only. It
runs `lsdtt-drainage-analysis` with `inline_gradient` (control),
`barnes2014b_float_step` and `barnes2014b_gradient` on the flow-benchmark
fixtures (alpine and lakeland, small and medium, fetched when not cached) and on
the 500 and 1000 m benchmark rasters when they are already cached. Every run
must exit 0, and each barnes run must end with `Rim cells without descent: 0`.
It prints the rounds and per-round counts. The fixture registry lives in
`tests/scripts/flow_benchmark_fixtures.py`, a stdlib module, because importing
`run_flow_benchmark.py` from a stdlib script failed the test-inventory audit (it
follows local imports into the function-level numpy and osgeo imports).

Round counts recorded on the small fixtures: alpine_small 2 (1, 0),
lakeland_small 2 (41, 0). The medium fixtures needed 4 rounds each.

### 10r.7 Round count and the cap

Before the cap was chosen, the loop was run to caps of 64 and then 256 rounds.
Deterministic counts: on lakeland_500 the float-step option converged after 76
rounds and the gradient option after 74, and on lakeland_1000 after 66 rounds.
The chain length was set by the ULP relief of the lake surface, not by the number
of flats. The cost per round was similar to round 1, so the total grew with the
round count, and a 64-round cap would have left cells without descent on these
DEMs. The chosen policy is four rounds plus a final strict-rise priority-flood
pass. It gives a routable surface with a bounded number of passes, and the pass
changes no cell that already drains. The report gives the cells it changed.

### 10r.8 The MFD underflow fallback

With a large exponent the float32 weights `S^p` underflow. A cell whose summed
weight is zero silently becomes a sink. This was not specific to the float step.
An instrumented count of the branch over the existing flow, MFD and
distributed-hydrology tests never fired it (all of those run at exponent 1). Runs by hand at
exponent 8 did reach it, on the alpine and lakeland small fixtures under the
default fill as well as under the float step. The earlier code failed there with
"invalid receiver weight in receiver graph", and left two silent sinks on the
lakeland fixture.

The fix recomputes the weights of a cell from `(S_k / S_max)^p` only when the
summed weight is not a normal finite float. F15 covers it (22 checks; the old
kernel fails five of them). The outputs were compared byte for byte before and
after the change: drainage-analysis MFD rasters at exponents 1 and 4 on the
valley, alpine small and lakeland small fixtures (both fills), and at exponent 8
on the valley under the default fill; and distributed-hydrology TOPMODEL MFD
output on the valley at exponent 1 (both fills) and exponent 8 under the default
fill. All were identical. The runs that failed before now complete.

### 10r.9 Lessons

- **Synthetic fixtures and RichDEM comparisons both missed the routing failure.**
  Synthetic rims sit metres above their flats, and RichDEM strands the same cells
  on real DEMs, so a comparison with RichDEM cannot detect a failure that RichDEM
  shares. Any flat step needs a routing gate on real terrain, not only a value
  comparison.
- **The later rounds must not raise low edges.** A low edge of an already
  resolved flat, raised by two units, can strand that flat again. This was found
  by reasoning about the invariant before writing the loop, and it is pinned by
  F27.
- **An oracle count is not always a function of the output.** RichDEM's printed
  pair count depends on its loop order where two raised flats touch, so O1
  asserts a range rather than an equality.
- **Pixi exports GDAL, PROJ and CONDA variables.** A system-GDAL CLI launched
  from a pixi Python loads pixi's GDAL plugins, with errors on stderr and a
  large slowdown per call. The flat-resolution harnesses strip those variables
  from the subprocess environment. Any future pixi-driven CLI harness must do the
  same.
- **Stream state leaks between stages.** An earlier stage left `std::cout` in
  fixed four-decimal mode, so the maximum raise printed as `0.0000 m`. The report
  now sets and restores the precision.

### 10r.10 What is not measured

- No speed claim is made for the barnes options. Timings on the server have not
  been recorded, and they must be before any statement about cost. The loop adds
  rounds, so its cost depends on the terrain.
- Peak memory on a large DEM (8000 × 8000 cells and beyond) has not been measured.
  The per-cell memory of the resolver is described in its header.
- The per-CLI tests check that each CLI runs the option and prints the report.
  They do not check the routed output of every CLI against a reference.

## 10s. Channel-link bearing and junction-angle statistics (2026-10-17)

`lsdtt-channel-network-analysis` (version 1.1.0) summarises and compares the
distributions of its channel-link bearings and junction angles by group within
one run: `write_channel_bearing_statistics` and `write_junction_angle_statistics`,
grouped by Strahler order, basin or region and optionally stratified by a second
key. The new keys are marked experimental. `lsdtt-bearing-comparison`, whose only
real input was this CLI's own `_channel_link_bearings` CSV, was removed in the same
change; its four analytic smoke tests became native checks (10s.2). The Python
`lsdtt3.compare_bearings` is unchanged.

### 10s.1 Where the code is

- `src/cpu/lsdtt_channel_bearing_statistics_cpu.{hpp,cpp}`: circular mean, mean
  resultant length R, circular variance and SD, the axial (doubled-angle) form,
  Rayleigh (Zar approximation), the D8 fraction, centred rose bins, Holm, the
  linear junction-angle statistics, the category-keyed summarisers and
  comparators, the reference-CSV reader and the CSV writers. The two-sample
  numerics are the unchanged functions of `src/utils/lsdtt_statistical_tests`.
- `src/cpu/lsdtt_channel_network_grouping_cpu.{hpp,cpp}`: the grouping keys
  (order with pooling, donor order pair, outlet and order-k basins, region),
  per-pixel zone labels (polygon field or code raster), majority-with-purity over
  a pixel list, and the membership CSV.
- `src/cpu/lsdtt_channel_network_analysis_cpu.hpp`: `lsdtt_channel_link_bearing`
  carries `sampled_flats`, and `lsdtt_junction_angle_record` carries
  `donor_1_sampled_flats`, `donor_2_sampled_flats` and `receiver_sampled_flats`:
  the pixels each bearing was fitted on. They are not written to any output.

### 10s.2 What the tests assert

All driver modes are in `tests/cpp/test_channel_bearing_statistics.cpp`
(`lsdtt-test-channel-bearing-statistics`), labelled `unit;network`, native.

| test | asserts |
|---|---|
| `unit_channel_bearing_statistics_circular` | The eight bearings {10, 50, 100, 170, 200, 260, 300, 345}: against themselves W1 = 0, Kuiper V = 0, p = 1; against +20 W1 = 20; rotation-invariant W1 = 0 with optimal rotation 340; against +180 as axial W1 = 0. Point masses: 10 vs 350 W1 = 20; axial 0 vs 90 W1 = 90; axial 10 vs 170 W1 = 20; disjoint arcs {10, 20} vs {190, 200} V = 1. Tolerance 1e-9. The first four are the cases the removed `smoke_bearing_comparison_*` tests pinned; their CLI outputs at removal were `8,8,0,0,360,0,0,0,1`, `8,8,0,0,360,20,0,0.125,1`, `8,8,0,1,360,0,340,0.125,1` and `8,8,1,0,180,0,0,0,1` (n_x, n_y, axial, rotation_invariant, period, W1, optimal rotation, V, p). Also the rotation-invariant size guard (n_a n_b above the limit gives NaN and a flag, the Kuiper result is still computed). |
| `..._summaries` | {350, 10} directional: mean 0 (not 180, not 360), R = cos 10°. {170, 10} axial: mean 0, R = cos 20°. {0, 90} axial: R = 0, mean NaN. {0, 90} with weights {1, 3}: mean atan2(3, 1) = 71.565°. Ten identical values: R = 1, Rayleigh p < 1e-4, SD 0. {0, 45, 90, 135} × 25: axial R = 0 (the four-fold D8 pattern cancels on doubled angles), D8 fraction 1. Kish n_eff, axial SD halving. Rose (w = 10): 355 → bin centred 0, 5.0 → bin centred 10, axial 175 → bin centred 0; counts and lengths sum to the totals; a width that does not divide 180 is refused. Category summaries drop excluded records (category -1) and count distinct blocks; a combined key (order × 100000 + code) needs nothing special. |
| `..._comparisons` | Holm {0.01, 0.04, 0.03} → {0.03, 0.06, 0.06} in input order; NaN is not counted. Four eligible groups: pairwise 6 rows, one_vs_rest 4 rows, a reference adds 4 rows; a group below the minimum size is summarised (`meets_min_group_size` 0) and takes no part. Stratified rows carry each stratum and the pooled `all`. Reference reader: the exact message `Column 'x' not found in <path>`, non-numeric and blank cells skipped, the order filter keeps only matching rows. |
| `..._junction_angles` | {60, 90} against itself: KS D = 0, linear W1 = 0; {30} vs {90}: W1 = 60; sectors 10 and 350 are 340 apart (linear, not circular); beta_min / beta_max do not change when the donor labels are swapped. |
| `..._trellis` | The synthetic trellis below: see 10s.3. |
| `..._grouping` | On the trellis: region halves (north / south) partition the links; a link that crosses the split has purity equal to its majority pixel fraction (12 links cross; 4 of them are below 0.6 and are exactly the ones a 0.6 purity threshold excludes); order-2 axial W1 north vs south 0.6°; order pooling ("2+"), donor order pairs ("1_2"); order-2 basins (below); innermost-wins nesting (299 junctions upslope of an inner trunk junction take it, the other 301 of the valley the outer outlet, all others none). |
| `..._isotropic_null` | A diamond-square surface (seed 2, 257 × 257 at 30 m, sd 100 m, threshold 20 px): all-link axial R < 0.15 (see 10s.4). |
| `fixture_channel_bearing_statistics` | Writes the CLI fixture (`--write-fixtures`). |

CLI tests (Python harness, `LSDTT_ENABLE_PYTHON_TESTS=ON`, fixture
`channel_bearing_statistics`):
`smoke_channel_network_analysis_bearing_statistics_order` (all six CSVs, the
summary and comparison headers, the group lines),
`rigorous_..._bearing_statistics_reference` (the rotated trellis against
`tests/data/channel_bearing_statistics/ns_trunk_bearings.csv` through the order
filter: four reference bearings, order-2 axial W1 printed `90.0`),
`rigorous_..._bearing_statistics_region` (region within order: order-2 axial W1
between the halves printed `0.0`), and six `failure_..._bearing_*` tests (stratify
equal to group, region key without a file, unknown key, `donor_order_pair` for
bearings, an off-grid region raster with no CRS, a missing reference column).

The existing CNA outputs did not change. Before the change, a run of the
previous binary with every output switch on (two DEMs: a 30 m COP30 crop of the
Allegheny Front and the 128 × 128 Lakeland crop, threshold 50 px) wrote 32 CSV and
FlatGeobuf files; after it, the same parameters gave 32 byte-identical files, and
so did the run with both new switches on. The Tokunaga `--expect-sha256` of
`rigorous_channel_network_analysis_legacy_threshold` still passes.

### 10s.3 The trellis fixture, and what it measured

10 m cells, 300 rows × 160 columns, z = 0.02 × (distance north of the south edge)
+ h(x), with h made of half-cosine flanks of relief 50 m. Four N-S valleys at
columns 6, 46, 114 and 154. Each valley has one 6-pixel flank (too short to reach
the 8-pixel channel threshold) and one 34-pixel flank that sheds one tributary per
row; the short flank is on the west in the two western valleys and on the east in
the two eastern ones, so tributaries flow west in half the valleys and east in
the other half. Boundaries `nnon`. The transposed surface (valleys E-W, draining
west, `nnno`) is the rotated case.

Two features of the design were forced by D8 and are worth knowing before the
fixture is changed:
- With both flanks channelled every trunk pixel takes a tributary from each side
  and has three donors, and `calculate_junction_angles` analyses only two-donor
  junctions; one channelled flank per valley gives two-donor junctions.
- At the open south edge the last row's flank drains along row 299 into the
  outlet pixel beside the trunk, so the outlet junction is order 3 and each valley
  has two order-2 basin mouths: the trunk (row 298) and that edge channel. The
  grouping test counts the trunk basins (4, one per valley) and the edge basins
  (4) separately.

Measured (full_link sampling, min_nodes 2): 2400 junctions, 2396 links (1204 of
order 1, 1192 of order 2), 1192 two-donor junctions. Order-1: 1200 E-W
tributaries and 4 N-S head links; axial mean 89.98°, axial R 0.9923 against
(n_trib − n_head)/(n_trib + n_head) = 0.9934; directional R 0.026. Order-2:
axial mean 0, axial R 0.9933, directional mean 180. Junction alpha median
88.67°. Kuiper V (order 1 vs 2, axial) 0.993. Rotated: order-2 axial mean 90,
directional mean 270; W1 against the N-S trellis order 2 89.40° axial and 89.70°
directional. The pinned tolerances (2° on means, 5° on alpha and the order-1
mean, ±0.05 on R, ±3° on W1) are derived from the geometry: the last one or two
flank steps are diagonal (the 0.02 down-valley gradient makes the diagonal
steeper where the cross-valley gradient falls below about 0.05).

### 10s.4 The isotropic null

On a D8 network of an isotropic surface, the all-link axial R measured (257 ×
257 diamond-square at 30 m, roughness 0.55, sd 100 m, threshold 20 px, full_link)
for seeds 1-8 was 0.0751, 0.0063, 0.0523, 0.0191, 0.0457, 0.0272, 0.0369 and
0.0149. Seed 1 is 1.998 times below the 0.15 bound, under the factor of 2 asked
for, so the test uses seed 2. The directional R on the same surfaces was 0.09 to
0.36: the directional form is not a usable isotropy statistic on a D8 network,
which is the reason the axial form is the one bounded.

### 10s.5 First real run: Allegheny Front, West Virginia (2026-10-17)

COP30 (`cop30_aws`, keyless), box −79.80/−79.05/38.85/39.15, EPSG:32617 at 30 m:
2185 × 1148 cells, tiles N38W080 and N39W080. Parameters: threshold 300 px,
`analyse_all_junctions`, `only_complete_basins`, full_link sampling,
`link_sampling_min_nodes` 5, orders ≥ 4 pooled (`channel_bearing_max_order_group`
4), independence blocks = order-3 basins. 4204 junctions, 3330 kept by the
complete-basin filter. Means and R are count-weighted unless marked; axial means
are in [0, 180), grid north of UTM 17N.

| order | n | n_blocks | axial mean | axial R | axial mean (length) | axial R (length) | directional mean | directional R | D8 fraction | median alpha |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1765 | 108 | 108.7 | 0.146 | 112.3 | 0.159 | 85.7 | 0.054 | 0.049 | – |
| 2 | 791 | 106 | 117.5 | 0.142 | 118.1 | 0.158 | 83.2 | 0.067 | 0.053 | 58.0 (n 611) |
| 3 | 311 | 83 | 136.6 | 0.091 | 128.1 | 0.168 | 125.2 | 0.166 | 0.029 | 71.7 (n 262) |
| 4+ | 107 | 0 | 79.8 | 0.230 | 70.5 | 0.192 | 31.3 | 0.083 | 0.093 | 75.0 (n 80) |

Junction alpha is grouped by the order at the confluence pixel, so there is no
order-1 row. n_blocks is 0 for 4+ because a link of order above the block order
lies in no order-3 basin.

Pairwise, axial (count-based; p assumes independent links):

| a vs b | n | n_blocks | W1 (deg) | Kuiper V | p_iid | p Holm |
|---|---|---|---|---|---|---|
| 1 vs 2 | 1765 / 791 | 108 / 106 | 1.9 | 0.056 | 0.38 | 0.49 |
| 1 vs 3 | 1765 / 311 | 108 / 83 | 4.7 | 0.107 | 0.046 | 0.14 |
| 1 vs 4+ | 1765 / 107 | 108 / 0 | 7.3 | 0.227 | 0.00086 | 0.0034 |
| 2 vs 3 | 791 / 311 | 106 / 83 | 3.9 | 0.094 | 0.24 | 0.49 |
| 2 vs 4+ | 791 / 107 | 106 / 0 | 8.8 | 0.249 | 0.00025 | 0.0013 |
| 3 vs 4+ | 311 / 107 | 83 / 0 | 10.6 | 0.309 | 7.9e-06 | 4.7e-05 |

Sensitivity (axial R by order 1 / 2 / 3 / 4+, and n):
- threshold 150 px: 0.149 / 0.128 / 0.116 / 0.085 (n 3265 / 1560 / 658 / 276);
  median alpha 55.1 / 69.5 / 77.7 for orders 2 / 3 / 4+.
- threshold 600 px: 0.164 / 0.148 / 0.038 / 0.251 (n 945 / 376 / 128 / 36);
  median alpha 61.8 / 76.3 / 77.9.
- `link_sampling_min_nodes` 2: 0.135 / 0.145 / 0.075 / 0.191 (n 2001 / 860 / 344 / 125);
  the D8 fraction rose from 0.049 / 0.053 / 0.029 / 0.093 to 0.121 / 0.098 / 0.087 / 0.168.

Grid-isotropy null: a diamond-square surface (`lsdtt-synthetic-surfaces`, seed 1,
same 2185 × 1148 at 30 m, mean 783.5 m and sd 261.1 m matched to the DEM, default
synthetic georeference) with the same CNA parameters: axial R 0.063 / 0.056 /
0.074 / 0.137 for orders 1 / 2 / 3 / 4+ (n 1914 / 887 / 333 / 235); directional R
0.231 / 0.296 / 0.347 / 0.467 with directional means 229-237°; median alpha 50.9 /
59.8 / 56.9 for orders 2 / 3 / 4+.

Wall times on a thermally limited WSL2 laptop (8 cores; one run each, includes
fill, routing and all outputs): 2.3 s (threshold 300), 3.5 s (150), 2.5 s (600),
3.6 s (min_nodes 2), 2.6 s (null); peak RSS about 217-222 MB. These are laptop
timings and not a performance claim.

The province run (region within order) has not been run: the province polygon is
to be drawn by hand on a hillshade.

### 10s.6 Automatic block order and block permutation (2026-10-17)

`channel_bearing_block_order` defaults to 0 = automatic: the largest Strahler order k
whose order-k basins holding at least one analysed record (link or junction, after the
basin selection) number at least `channel_bearing_min_blocks` (default 30). With B
blocks the smallest attainable permutation p is 1 / (number of distinct label
assignments + 1); with 30 blocks even a 5-vs-25 split has C(30, 5) = 142 506
assignments, so the 999 permutations set the resolution, while with about 10 blocks an
unbalanced split can have fewer than 100. If no order reaches the minimum, no blocks
are defined: every permutation p is nan with the note "network too small for block
permutation: largest order-k basin count n < min_blocks", the descriptive statistics
are written, and a WARNING is printed. A positive value is used as given, with a note
when it has fewer blocks than the minimum. The counts per order are written to
`<wf>_channel_network_block_order_counts.csv` (`strahler_order, n_basins_with_records,
chosen`).

`channel_bearing_permutations` (default 999; 0 = off) and
`channel_bearing_permutation_seed` (default 1) give each comparison a block permutation
p = (1 + #{T* >= T}) / (1 + permutations), with T = Kuiper V (axial or directional) for
bearings and KS D for junction angles, over the records that have a block and a finite
value (the number left out is in the note). Spatial keys (basin, region, lithology):
each block takes the majority group of its records in the comparison (ties to group a)
and the block labels are permuted among blocks; `permutation_statistic` is T on those
block-majority samples. Orders (strahler_order, donor_order_pair): groups are permuted
among the records of each block. The p is nan, with a note, when a group has no blocked
record, when every block has the same majority group, or (within blocks) when no block
holds both groups. Records above the block order have no block: their summary rows say
so in `note`, and n_unblocked counts them. Each comparison row has its own
std::mt19937_64 seeded from the seed and the row's labels, so a p does not depend on
which other comparisons were run. The unit tests are listed with
`unit_lithology_channel_network_stats` in section 8c.4.

West Virginia, the run of section 10s.5 with the automatic order (order-k basins with
links or junctions, k = 1..6): wvbig 1765 / 471 / 108 / 23 / 6 / 2, chosen order 3;
the isotropic null 1914 / 525 / 136 / 33 / 9 / 2, chosen order 4. Pairwise by order,
permutation p (999) beside the iid p:

| run, type | a vs b | n | n_blocks | Kuiper V | p_iid | permutation p |
|---|---|---|---|---|---|---|
| wvbig, axial | 1 vs 2 | 1765 / 791 | 108 / 106 | 0.056 | 0.379 | 0.329 |
| wvbig, axial | 1 vs 3 | 1765 / 311 | 108 / 83 | 0.107 | 0.046 | 0.012 |
| wvbig, axial | 2 vs 3 | 791 / 311 | 106 / 83 | 0.094 | 0.244 | 0.539 |
| wvbig, axial | x vs 4+ | – / 107 | – / 0 | 0.227-0.309 | < 0.001 | nan (order 4+ has no order-3 block) |
| wvbig, directional | 1 vs 2 | 1765 / 791 | 108 / 106 | 0.041 | 0.861 | 0.140 |
| wvbig, directional | 1 vs 3 | 1765 / 311 | 108 / 83 | 0.142 | 0.0008 | 0.006 |
| wvbig, directional | 2 vs 3 | 791 / 311 | 106 / 83 | 0.137 | 0.006 | 0.031 |
| null, axial | 1 vs 2 | 1914 / 887 | 33 / 33 | 0.044 | 0.700 | 0.255 |
| null, axial | 1 vs 3 | 1914 / 333 | 33 / 32 | 0.087 | 0.198 | 0.562 |
| null, axial | 1 vs 4+ | 1914 / 235 | 33 / 21 | 0.114 | 0.079 | 0.853 |
| null, axial | 2 vs 3 | 887 / 333 | 33 / 32 | 0.080 | 0.439 | 0.241 |
| null, axial | 2 vs 4+ | 887 / 235 | 33 / 21 | 0.114 | 0.130 | 0.620 |
| null, axial | 3 vs 4+ | 333 / 235 | 32 / 21 | 0.110 | 0.371 | 0.965 |
| null, directional | 1 vs 4+ | 1914 / 235 | 33 / 21 | 0.209 | < 0.0001 | 0.014 |
| null, directional | 2 vs 4+ | 887 / 235 | 33 / 21 | 0.172 | 0.0005 | 0.289 |

Junction alpha (KS), wvbig 2 vs 3: D 0.228, p_iid < 0.0001, permutation p 0.001; null
2 vs 3: D 0.139, p_iid 0.0014, permutation p 0.037; null 2 vs 4+: D 0.126, p_iid 0.021,
permutation p 0.044. Wall times with 999 permutations (WSL2 laptop, thermally
limited, one run each, not a performance claim): wvbig order 3.3 s, null 4.4 s,
lithology by unit 8.1 s.

### 10s.7 Plateau and Valley and Ridge by map-unit age (2026-10-17)

A cross-check of the province split built from the map, not drawn: each Macrostrat
`name` in the wvbig map is assigned by its older age bound (`best_age_bottom`):
<= 358.9 Ma (Mississippian or younger) = allegheny_plateau, otherwise
valley_and_ridge. 28 names: 12 plateau (including Quaternary Alluvium, 2.6 Ma, by the
rule) and 16 Valley and Ridge; no name has conflicting bounds. The crosswalk is a
`value,code,name` CSV used through the lithology key (`lithology_field_name: name`).
By construction, deep plateau valleys that cut into Devonian rocks (Cheat near Parsons,
Blackwater Canyon) are labelled valley_and_ridge by this rule. Same parameters as
section 10s.5, automatic block order (3), pairwise; 153 links are "mixed" (83 with a
plateau majority, 70 Valley and Ridge).

| province, order | n | n_blocks | axial mean / R (count) | axial mean / R (length) | D8 fraction |
|---|---|---|---|---|---|
| plateau, 1 | 790 | 64 | 113.2 / 0.142 | 118.5 / 0.159 | 0.065 |
| V&R, 1 | 874 | 71 | 105.4 / 0.156 | 108.6 / 0.166 | 0.040 |
| plateau, 2 | 343 | 57 | 135.4 / 0.107 | 139.5 / 0.100 | 0.055 |
| V&R, 2 | 406 | 59 | 109.1 / 0.185 | 110.6 / 0.221 | 0.052 |
| plateau, 3 | 144 | 33 | 170.4 / 0.107 | 167.9 / 0.146 | 0.028 |
| V&R, 3 | 160 | 56 | 118.5 / 0.136 | 116.9 / 0.261 | 0.031 |
| plateau, 4+ | 59 | 0 | 66.0 / 0.118 | 42.7 / 0.204 | 0.136 |
| V&R, 4+ | 45 | 0 | 83.2 / 0.355 | 86.5 / 0.294 | 0.044 |
| plateau, all | 1336 | 68 | 120.2 / 0.101 | 126.5 / 0.105 | 0.061 |
| V&R, all | 1485 | 76 | 106.4 / 0.163 | 109.6 / 0.192 | 0.042 |

Plateau vs Valley and Ridge (axial / directional): order 1 V 0.057 / 0.093, p_iid
0.59 / 0.019, permutation p 0.83 / 0.70; order 2 V 0.141 / 0.111, p_iid 0.015 / 0.16,
permutation p 0.31 / 0.59; order 3 V 0.169 / 0.356, p_iid 0.18 / < 0.0001, permutation
p 0.34 / 0.032; all orders V 0.082 / 0.121, p_iid 0.003 / < 0.0001, permutation p
0.27 / 0.24; order 4+ has no block (nan). Junction alpha median (all junctions, by the
junction pixel): plateau 63.4 [60.9, 67.7] (n 451), Valley and Ridge 60.8 [58.8, 65.3]
(n 502); order 3: 77.4 [71.4, 80.9] (124) and 66.9 [59.5, 73.1] (138), KS D 0.178,
p_iid 0.027, permutation p 0.060.

### 10s.6 Angles on half-open domains in the CSVs (2026-10-19)

The writers print 10 significant digits. A link bearing of 359.99999999999579 (in
[0, 360), as `_channel_link_bearings.csv` carries it) printed as `360` in
`_channel_bearing_group_membership.csv`, and its axial value as `180`, outside the
documented [0, 360) and [0, 180). `format_channel_bearing_angle(angle, period)`
(`lsdtt_channel_bearing_statistics_cpu`) now writes a value that rounds to the period
as `0`; it is used for the membership bearing and axial bearing, the summary circular
means (count and length) and the comparison `optimal_rotation_deg`. Rose bin edges
are unchanged (an upper edge of 360 is correct there).
`unit_channel_bearing_statistics_summaries` checks the formatter (359.99999999999579
and 179.99999999999579 give `0`; 359.9999, 179.99999, 0 and NaN unchanged) and the
membership writer row for row.

Measured: re-running `wvbig_cna_order` and `wvnull_cna_order` (the parameter files of
10s.5) changed only `_channel_bearing_group_membership.csv`: 2 rows on wvbig_order
(links 1962-1910 and 1986-1962) and 3 on wvnull_order (113-51, 3619-3575, 4206-4163),
each `360,180` to `0,0`. Every other output of both runs, the summaries and comparisons
included, is byte-identical.

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
