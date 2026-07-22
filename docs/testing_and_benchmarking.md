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
| `flow` | — | flow routing and depression-handling properties |
| `fetch` | — | native DEM fetching; offline logic by default |
| `network` | network | tests that actually download |
| `benchmark` | network + RichDEM/TopoToolbox | comparison against external references on real DEMs |
| `performance` | — | in-process timings of the flow routines |

```sh
cmake -S . -B build && cmake --build build -j
ctest --test-dir build --output-on-failure          # default: offline suite

ctest --test-dir build -L unit --output-on-failure
ctest --test-dir build -L flow --output-on-failure
```

Opt-in tiers, each behind its own CMake option so enabling one does not drag in
the others:

```sh
cmake -S . -B build \
  -DLSDTT_ENABLE_NETWORK_TESTS=ON \   # tests that download
  -DLSDTT_ENABLE_BENCHMARKS=ON \      # RichDEM comparisons (needs pixi)
  -DLSDTT_ENABLE_PERF_TESTS=ON        # performance timings
```

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

```sh
cmake -S . -B build -DLSDTT_ENABLE_PERF_TESTS=ON
ctest --test-dir build -L performance

# or directly, including a real raster
build/tests/lsdtt-test-flow-performance out_dir --full \
    --raster tests/data/dem_cache/alpine_small.tif
```

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
| **Sentinel pinning suite** | A matrix of exact concentration/ε values — every (nuclide × scheme) at several (ε, d, d_s, d_t, t, ε_old/ε_new, S_tot) points, covering CAIRN Eq. (14), knickpoint Eqs. (7)–(8), the Λ_v solve and the inversion — frozen to full double precision (17 significant digits) in `tests/expected/cosmo/numerics_sentinel.csv` and compared at 1e-12 relative tolerance. Deliberately stricter than physical accuracy requires: its job is to trip on *any* change — refactor, algebraic rearrangement, compiler flag, parameter tweak. |
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
| **vs CRONUS v3 St** | **measured 2026-07-16: +1.8 % to +6.65 %, mean +3.6 %; band frozen at 10 %** | Balco (2017) model-1A muons; grows with ε as the muon share rises |
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

---

## 8. Implicit hillslope solver validation

The backward-Euler + damped-Picard solver for the Roering (1999) nonlinear
hillslope flux law (`compute_roering_hillslope_transport_implicit`) is validated
by the standalone driver `cli/lsdtt-test-hillslope-implicit.cpp`, plus in-LEM
hillslope-strip runs driven entirely from `lsdtt-lem` parameter files. The driver
emits `plot_data/*.csv` for plotting.

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
| 5 | Linear-mode contrast (`nonlinear = false`) against the exact parabola | 2nd-order in `h`; quantifies the gap to the Roering solution, confirming the nonlinear term is live |
| 6 | Timestep independence across `dt_init = tau/{10, 40, 160}` | rel ≈ **2.4e-6** |
| 7 | Conservation (U = 0, `'nnnn'`, no sinks) — the discrete volume budget must close | max relative volume drift **2.3e-8** |
| 8 | Performance benchmark, one implicit solve, flat, U = 0 | see §8.3 |
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
to the standalone driver**, whose steady-state tests all finish on full Picard, so
every unit-level check passed — while coupled LEM steady states were biased about
**40 % high**. It was only exposed by the in-LEM integration run. A lone
semi-implicit solve now uses `omega = 1`.

> A fast path that is only taken in regimes the unit tests do not finish in is a
> fast path the unit tests cannot validate. Where a code path exists specifically
> for the cheap case, an end-to-end test that actually spends time in it is not
> optional.

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

## 12. Detailed plans

Working plans with the full derivation of the findings above live in
`docs/plan_docs/flow_routing_benchmark_plan.md`, `docs/plan_docs/fetch_raster_cli_plan.md`,
`docs/plan_docs/cosmo_cli_testing_plan.md` and `docs/plan_docs/roering_implicit_transient_plan.md`.
Those are local working documents and **may not be present in any given
checkout**; nothing here depends on them. **This document is self-contained and
is the durable record.**
