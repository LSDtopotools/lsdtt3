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

Both tools implement both algorithm families, and both default to a *different*
one:

| | reconstruction (O(N), FIFO queue) | priority-flood (O(N log N), heap) |
|---|---|---|
| lsdtt3 | `vincent1993` | `barnes2014` ← **CLI default** |
| TopoToolbox | `hybrid=True` ← **default** | `hybrid=False` |

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
`hybrid=True`, and a new `fill_priority_flood` case pins `barnes2014` against
`hybrid=False` and RichDEM (Zhou 2016). `topotoolbox_reference.py` gained a
`fill_priority_flood` routine to expose the non-hybrid path. Both families are
now measured; neither is compared across.

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
`docs/what_is_lsdtt3/flow_routing.md`.

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
matrix; see the plan's §9.5f for the measurements and validation caveat.

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

These are outstanding work rather than findings, so they are tracked in the
benchmarking plan (uncommitted scratch) rather than here. This document records
what was measured; it is not a to-do list.

#### Sequential cases are no longer swept (changed 2026-09-01)

Depression filling (both algorithm families), least-cost breaching, and
D8/D-inf/MFD graph construction are **ordered traversals**: each step depends on
the result of the one before it. There is no parallelism in them to measure, in
lsdtt3 or in any of the references.

Sweeping threads over them measured nothing. It multiplied the campaign length
and — worse — filled the scaling table with flat rows that invite a reader to
conclude something is broken. Nine cases are now marked `sequential` in the case
registry and are **measured once, at the baseline thread count**:

`fill`, `fill_priority_flood`, `fill_epsilon`, `breach`, `accum_d8`,
`accum_dinf`, `accum_mfd`, `fill_spill_merge`, `fill_wang_liu`

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

| Suite | Cases | Status |
|---|---|---|
| `flow` | fill, fill_priority_flood, fill_epsilon, breach, breach_least_cost, accum_d8, accum_d8_no_svector | timing; values asserted separately in §4-5 |
| `mfd` | accum_dinf, accum_mfd | **timing only** -- Tarboton vs Quinn vs Freeman proportioning genuinely differ |
| `surface` | slope_only, polyfit_all, hillshade | **timing only** -- curvature definitions and aspect handedness differ |
| `depression` | fill_spill_merge, fill_wang_liu | **timing only**, and the depression hierarchy has no counterpart at all |

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
  lsdtt3's and RichDEM's Lindsay (2016). `breach` is now lsdtt3 vs RichDEM, and
  the new `breach_least_cost` compares lsdtt3's least-cost breach with Whitebox's
  at max_dist 100 cells with unresolved depressions filled. Whitebox's
  `minimize_dist` and `flat_increment` defaults are not reported by the
  package, so they are not verified equal to lsdtt3's.

Campaigns before 2026-09-13 carry the old pairings in those three rows.

`fill_spill_merge` measures lsdtt3's nested depression hierarchy (Barnes 2020)
plus Fill-Spill-Merge (Barnes 2021), lsdtt3 only. It is recorded because the
absolute cost of the hierarchy is a number users ask for. Until 2026-09-13 it
was paired with Whitebox's `sink`, which only labels depression cells, so that
ratio compared nothing and the pairing was removed. Note the hierarchy is
not separately switchable -- it is built as the first half of the
`barnesetal2021` fill, which is the only CLI path that reaches it.

Two cases reuse rows rather than re-measuring (`reuse_from` in the case
registry, `reused_from` in the CSV). `accum_d8_no_svector` copies the
reference rows of `accum_d8`, because the references run the identical call in
both. `fill_wang_liu` copies lsdtt3's `fill_priority_flood` row, because lsdtt3
has no Wang & Liu fill and its side of that case was the same barnes2014 command.

**Memory.** Every command is launched through `tests/scripts/rusage_launcher.py`.
A Linux child inherits its parent's lifetime peak RSS and `wait4` reports it as
the child's own. Measured on WSL2: after a Python parent touched and freed
1.5 GiB, `/bin/true` launched from it reported 1509.5 MiB; from a fresh process,
9.0 MiB. The 2026-09-09 forth02 campaign launched directly from a driver that
had peaked at 3.5 GB in `describe_fixture()`, and 441 of its 504 rows read
exactly 3556.98 MiB. Memory columns from campaigns before 2026-09-13 are invalid
except where a value exceeds that floor. Each campaign now records `rss_floor`
(`/bin/true` through the launcher: 11.0 MiB on WSL2, against 93.0 MiB launched
directly) and fails a tripwire above 64 MiB.

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

**Speed.** Not yet measured on a server; the before/after is part of the next
server allocation (`315e81a`, `7e81a39` and `eb65499` built side by side).

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
fixed under plan phase R2b. It is recorded because the *shape* of it recurs.

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
from the window OR fails it. The speedup was confirmed in direction on the
laptop; the quotable figure awaits the server campaign.

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
`ctest -L fluvial_surfaces -LE heavy`. What counts as a genuine regression:

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
