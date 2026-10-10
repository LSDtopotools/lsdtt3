# Flow Routing in lsdtt3

This page explains how flow routing is implemented in `lsdtt3`, what each
method is for, and where `lsdtt3` deliberately differs from similar software packages.

---

## Graph vocabulary used on this page

Flow routing is graph (in the mathematical sense) work, and to read this section you will need to know what some terms mean. 
They are defined once here.

- **Node and edge.** Each valid DEM cell is a **node**. An **edge** points from a
  cell to a cell it sends water to. "Receiver" means the cell downstream;
  "donor" means a cell upstream that sends water in.

- **Tree.** A graph where every node has **exactly one** outgoing edge, and
  following those edges never returns you to where you started. D8 produces a
  tree (strictly, a forest - one tree per outlet). Because each cell has one
  receiver, the set of cells draining to a point is a clean, non-overlapping
  group, which is what makes basins and drainage area exact.

- **DAG - directed acyclic graph.** "Directed" because edges have a direction;
  "acyclic" because you can never go round in a circle. A DAG relaxes the tree
  rule: a node may have **several** outgoing edges. D-infinity and MFD produce
  DAGs, because a cell splits its water between two or more receivers. Every
  tree is a DAG; not every DAG is a tree. The difference is why the D8 flow routing
  routines are numerically different from routines using multiple flow directions. 

- **Topological order.** An ordering of the nodes in which every node appears
  *after* all the nodes that feed into it. Any DAG has one; a graph with a cycle
  has none. It is what lets flow accumulation be computed in a single sweep:
  when you reach a cell, everything upstream of it is already finished. Many of the 
  big gains in computational efficiency come from this, as compare to brute force 
  searches, which were the standard method in DEM processing prior to the introduction
  of graph theory.

- **CSR - compressed sparse row.** A compact way to store a graph in which each
  node has a variable number of edges, using two flat arrays instead of a list
  per node: one array holds all the edges end to end, and an *offsets* array
  records where each node's edges begin. Node `i`'s edges are the slice from
  `offsets[i]` to `offsets[i+1]`. It is the standard sparse-matrix layout, and it
  is used here because a D-infinity or MFD cell may have one to eight receivers.

- **DFS - depth-first search.** Traversing a graph by following one branch as far
  as it goes before backtracking to the next. A **preorder** DFS records each
  node when it is first reached. Done from an outlet up a D8 tree, it lists every
  cell of a drainage basin consecutively - which is what `s_vector` exploits. More on that below.

- **Kahn's algorithm and the "ready queue".** A way of producing a topological
  order without recursion (Kahn 1962). Count how many donors each node is still
  waiting on. Any node waiting on none is **ready**, so put it in a queue. Pop a
  node, emit it, and decrement its receiver's count; if that receiver's count
  reaches zero it is now ready, so push it. Repeat until the queue empties. If
  nodes remain unemitted, the graph had a cycle - which is a useful built-in
  check on a routing graph that should not have one. This is another algorithm that has been
  around in computer science for a long time but has only recently been injected into 
  flow routing routines. 

---

## Routing objects and use-cases

`lsdtt3` has two flow-routing backends:

- **`lsdtt_flowinfo` (D8)**: one receiver per valid cell, tree/forest topology,
  optional `s_vector` index for O(1) upslope-membership tests and contiguous
  upslope slices.
- **`lsdtt_flowgraph_cpu` (D-infinity / MFD)**: multiple weighted receivers per
  cell, DAG topology (CSR storage), intended for distributed-flow metrics such
  as D-infinity/MFD area and TWI.

In practice:

- Use D8 (`lsdtt_flowinfo`) for channel-network and basin workflows.
- Use D-infinity or MFD (`lsdtt_flowgraph_cpu`) when fractional flow partition
  is the intended model.

---

## Provenance of the routing schemes

The three routing schemes in `lsdtt3` are not variants of one idea. They are
three answers to the question of how much a grid cell
is allowed to know about where its water goes. Each was a response to a specific
failure of the one before. Please indulge us in a bit of pedantic history telling. 

### D8 - O'Callaghan and Mark (1984)

The original. Each cell sends **all** of its water to exactly one of its eight
neighbours: the one of steepest descent. The network that results is a forest of
trees, every cell having one receiver and any number of donors.

The crucial detail, and the one most often got wrong, is that steepest descent is
**drop divided by distance**, not lowest elevation. A diagonal neighbour is
`sqrt(2)` cell-widths away, so it must be `sqrt(2)` times lower to win against a
cardinal one. `lsdtt3` implements the drop-over-distance rule. Implementations
that instead pick the lowest neighbour by elevation give systematically different
networks - see "Why RichDEM's D8 differs" below, which is a worked example of
exactly this.

D8's weaknesses were understood from the beginning and are geometric, not
computational:

- **Flow is quantised to 45 degrees.** A hillslope draining at 20 degrees is
  represented as alternating cardinal and diagonal steps.
- **It cannot diverge.** On a ridge, a fan, or any convex surface, real water
  spreads and D8 cannot. Every cell has one receiver, so upslope area is a set of
  disjoint trees, and adjacent hillslope cells can be assigned to different
  basins on the strength of a centimetre of elevation difference.
- **Flat areas have no defined direction at all**, which is why depression
  handling (below) is a prerequisite rather than a refinement.

What it buys in exchange is the reason it is still the default: the network is a
tree, so upslope area is exact and integer, basins partition the DEM cleanly,
channel extraction and chi analysis are well posed, and the whole thing is
linear-time (see "Why D8 is effectively O(n)").

### MFD - Freeman (1991) and Quinn et al. (1991)

Multiple-flow-direction routing attacks the divergence problem directly: a cell
gives water to **every** downslope neighbour, in proportion to

    weight_i  proportional to  (tan beta_i)^p

where `tan beta_i` is the slope to neighbour `i` and `p` is a tuning exponent
(`mfd_exponent`). Freeman and Quinn et al. arrived at the same formulation
independently in the same year, from different directions - Freeman from
catchment-area computation, Quinn et al. from TOPMODEL-style distributed
hydrology, where the D8 assumption that a hillslope drains along a single line
made the wetness index behave badly. MFD looks much more like a finite difference approach 
than a flow routing approach. 

The exponent is the whole argument. `p = 1` spreads flow broadly and is defensible
on hillslopes; large `p` concentrates flow and converges toward D8 in the limit.
There is no physically derived value, which is MFD's central weakness: the answer
depends on a parameter chosen by the user.

MFD also gives up the tree. Upslope "area" becomes a fractional, overlapping
quantity - a cell can be partly in several basins - so basin partitioning,
chi analysis, and O(1) upslope queries are no longer available in the same form.

### D-infinity - Tarboton (1997)

Tarboton's method lies between these end members. The eight neighbours define
eight triangular facets around the central cell; the steepest downslope direction
is computed **on those facets** as a continuous angle, not snapped to a
neighbour. Flow is then split between the **two** neighbours that bracket that
angle, in proportion to how close the angle is to each.

This keeps most of what MFD was after - flow is no longer quantised to 45
degrees, and divergence on a convex surface is representable - while limiting
dispersion to two receivers rather than up to eight. It also removes MFD's free
parameter: the partition is fixed by the geometry.

The price is that it is still not a tree, so the same structural limits apply as
for MFD, and the facet computation is more expensive per cell.

### How to choose

| | D8 | D-infinity | MFD |
|---|---|---|---|
| receivers per cell | 1 | 1-2 | up to 8 |
| flow direction | 45-degree steps | continuous | continuous, dispersed |
| free parameter | none | none | `mfd_exponent` |
| structure | tree | DAG | DAG |
| exact basins, chi, O(1) upslope | yes | no | no |

Use D8 when the tree structure is the point: basin delineation, channel networks,
chi analysis, most cosmogenic work. Use D-infinity or MFD when hillslope
divergence matters and you do not need a partition - hydrological wetness indices,
sediment or contaminant dispersal, debris-flow inundation.

### How the D8 graph and the multi-receiver graph differ

The difference between D8 and D-infinity/MFD is not a tuning choice; it changes
the *kind of graph* the library is holding, and most of the downstream
consequences follow from that.

| | D8 | D-infinity / MFD |
|---|---|---|
| outgoing edges per cell | exactly 1 | 1 to 8 |
| graph type | **tree** (forest) | **DAG** |
| edges carry a weight? | no - all the water goes | yes - a fraction per edge |
| storage | one receiver per cell, a fixed-size array | **CSR**, because the count varies |
| upslope set of a cell | a disjoint group of cells | overlapping, fractional |
| drainage area | integer count of cells | a weighted sum |
| basins | partition the DEM exactly | overlap; no clean partition |

Two practical consequences:

- **You can index a tree in a way you cannot index a DAG.** In a tree, the cells
  upslope of a point form one contiguous block under a preorder DFS, so
  membership is an interval test - this is `s_vector`, below. In a DAG a cell can
  be partly upslope of two different points at once, so no such interval exists.
- **Storage must be variable-width for a DAG.** One receiver per cell fits a
  plain array; one-to-eight weighted receivers per cell is exactly the case CSR
  exists for.

`lsdtt_flowgraph_cpu` builds the receiver DAG in CSR form and walks it in
topological order using a Kahn ready queue (see the vocabulary above), following
the multi-receiver linear-layout idea in Anand et al. (2020).

### Hydraulic graph routing (GraphFlood lineage)

`lsdtt_hydraulic_graph_cpu` is a graph-based hydraulic solver: rather than
recomputing routing each step, it fixes the routing graph and advances water
depth over it through time. This is the approach of **GraphFlood**
(Gailleton et al., 2024), which approximates 2-D shallow-water hydrodynamics by
solving for a steady water surface on the flow graph, at a small fraction of the
cost of a full hydrodynamic model. `lsdtt3`'s implementation is credited to those
concepts in the code comments.

The relevance to this page is that GraphFlood consumes a **multi-receiver** graph:
water spreading across a floodplain is precisely the situation D8's single
receiver cannot represent, so the DAG is not an optional refinement here but a
prerequisite. Much of the logic behind our MFD algorithms stems from the Anand paper and 
the Gailleton et al. 2024 papers. 

---

## Braun and Willett: graph theory enters the field

Braun and Willett (2013) is not a routing method, and it is worth being precise
about what it did contribute, because the vocabulary it introduced is now
everywhere in this codebase.

**The problem it solved.** Landscape evolution models integrate the stream power
law, in which a cell's incision depends on its own elevation and on the elevation
of the cell it drains to. Solved explicitly, the timestep is limited by the
Courant condition and shrinks as resolution rises. Solved implicitly, you appear
to need a large sparse linear solve every step.

**The observation.** Under D8 the drainage network is a **directed tree**: every
cell has exactly one receiver, so the dependency graph of the implicit update is
a tree, not a general sparse matrix. If you process cells in an order where every
cell comes *after* the cell it drains into, each update needs only a value already
computed. No matrix, no iteration - one sweep, `O(n)`, and unconditionally stable.

**The machinery.** Braun and Willett laid out how to build that order in linear
time, and their construction is a textbook graph algorithm dressed in
geomorphological clothes:

1. count donors per cell (`NDonorsVector`);
2. prefix-sum those counts into offsets (`DeltaVector`) and fill a flat donor
   list (`DonorStackVector`) - this is a **counting sort into CSR adjacency**,
   exactly the layout a sparse-matrix library would use;
3. depth-first traverse from each base level, pushing cells onto a stack
   (`SVector`) - a **preorder DFS**, giving both the update order and, as a
   by-product, contiguous upslope intervals.

Those four array names appear verbatim in LSDTopoTools2's `LSDFlowInfo`, which
implements the paper directly.

**What lsdtt3 kept and what it changed.** The *idea* - order the cells so that
donors precede receivers, then sweep once - is retained and is the basis of
every accumulation in the library. The *construction* is different:

- `lsdtt3` orders cells with **Kahn's algorithm** (Kahn 1962) on a donor
  *count*, not on a stored donor list. Every cell holds a single number: how many
  donors it is still waiting for. Cells waiting for none are emitted, and each
  time a cell is emitted its receiver's count drops by one; when a receiver
  reaches zero it becomes ready in turn.

  **Why change it?** Braun and Willett's construction has to know, for each cell,
  *which* cells drain into it, so it builds the full donor list (the CSR arrays).
  Kahn's algorithm only ever needs to know *how many* are left. That difference
  is worth three things:

  - **Memory.** A D8 cell has at most 8 donors, so the count fits in one byte.
    The donor list it replaces costs 12 bytes per cell (`NDonorsVector` +
    `DeltaVector` + `DonorStackVector`), and building it needs a second pass and
    a temporary. See the memory section below.
  - **Cycle detection for free.** If the queue empties while cells remain
    unemitted, the receiver graph contains a cycle. A recursive DFS would have to
    be instrumented to notice; here it falls out of the loop condition, and
    `lsdtt3` uses it to mark the object non-acyclic rather than producing a
    silently wrong ordering.
  - **No recursion.** The DFS formulation recurses to the depth of the longest
    flow path, which on a large DEM can be tens of thousands of cells - a real
    stack-overflow risk. Kahn's algorithm is a flat loop over a queue.

  The trade is that Kahn's order is breadth-like rather than depth-first, so it
  does **not** produce contiguous upslope blocks as a by-product. That is why
  `s_vector` is a separate, optional structure in `lsdtt3` rather than something
  the ordering hands you for nothing.
- The preorder-DFS product is kept, but as an **optional** structure
  (`s_vector`), built only when O(1) upslope queries are actually wanted. 
  Note that for most geomorphology operations you will want to retain the `s_vector` 
  so the "savings" gained from leaving it out are a bit notional.  

**What "owing a debt" to the paper means here.** `lsdtt3` does not run Braun and
Willett's algorithm: none of their four arrays exists in the code. However, 
the code in lsdtt3 has legacy from LSDTopoTools2 and LSDTopoTools, which 
were direct implementations of the Brain and Willet algorithm, coded very shortly 
after the lead developer (SMM) was blown away by Jean's presentation at EGU in 2013. 
What lsdtt3 takes from Braun and Willet is the *framing*: that a D8 drainage network is a directed tree, that
the useful question is what order to visit its cells in, and that the right order
turns an apparently expensive implicit calculation into one linear sweep. That
insight is now so embedded in the field that it is easy to forget it was an
argument someone had to make. Every accumulation in this library is a single pass
over a topological order because of it, even though the order is built by a
different algorithm. 

The paper is therefore cited for the idea, not the implementation - and the
shared vocabulary in LSDTopoTools2 (`SVector`, `DeltaVector` and the rest) is a
direct trace of the older code following the paper array for array.

---

## Why D8 is effectively O(n)

For `n` valid cells, D8 build and accumulation are linear in cell count (with a
small constant factor):

- **Receiver assignment**: each cell checks at most 8 neighbours.
- **Donor counting**: one pass over cells.
- **Topological ordering**: each cell enqueued/dequeued once; one receiver edge
  processed per cell.
- **Accumulation**: one forward push along the topological order.

So work scales as a bounded number of full-grid passes, not an all-to-all
operation.

---

## Depression filling: two algorithms, one result

`lsdtt3` ships two depression-filling backends. They produce **bit-identical**
output at `fill_min_slope: 0` and differ only in cost.

| `fill_algorithm_choice` | algorithm | complexity | memory access |
|---|---|---|---|
| `barnes2014` (CLI default) | Improved Priority-Flood (Barnes et al. 2014) | O(N log N) — nearly every cell through a binary heap | heap-ordered, cache-hostile |
| `vincent1993` | grayscale reconstruction by erosion (Vincent 1993), hybrid form: forward scan, backward scan, then a FIFO queue for the residue | O(N), small constant | almost entirely sequential |

Measured on a 1000×1000 Copernicus tile (Himalaya), `fill_min_slope: 0`:

| backend | compute | throughput |
|---|---|---|
| `barnes2014` | 0.1282 s | 7.80 Mcells/s |
| `vincent1993` | 0.0446 s | **22.42 Mcells/s** |

**2.87× faster, 0 differing cells out of 1,000,000.**

### Which to use

- **`fill_min_slope: 0`** → `vincent1993`. Same surface, ~3× faster.
- **`fill_min_slope > 0`** → `barnes2014`, and there is no choice.
  Reconstruction fills to exactly the spill level and has no mechanism for
  imposing a gradient across the resulting flats. Priority-flood applies the
  increment inline as it floods. Flats left by a `min_slope: 0` fill can
  instead be resolved with the opt-in Barnes et al. (2014b) options, which
  make flow converge across the flat; see
  [Flat resolution: Barnes et al. (2014b)](#flat-resolution-barnes-et-al-2014b).
- **Periodic boundaries** → `barnes2014`. The reconstruction fast path has its
  own neighbour walk and does not wrap.

Note that lsdtt3's D8 router **requires** a gradient on flats — it throws
`could not assign a receiver` on a truly flat surface — so a `vincent1993` fill
is an endpoint (or an input to something other than D8 routing), not a
precursor to routing. That is why `fill_min_slope` defaults to `1e-4`.

!!! note "The library dispatches; the CLI does not"
    `fill_pits()` selects the backend automatically from `min_slope`
    (`<= 0` → reconstruction). `lsdtt-raster-preprocessing` instead dispatches
    on `fill_algorithm_choice`, which defaults to `barnes2014`, so a CLI run
    with `fill_min_slope: 0` pays the heap cost unless you ask for
    `vincent1993` explicitly.

### This is also what TopoToolbox does

TopoToolbox's `fillsinks()` takes a `hybrid` flag defaulting to **`True`** —
the same reconstruction algorithm as `vincent1993`. Measured on the same tile:

| | reconstruction | priority-flood | ratio |
|---|---|---|---|
| lsdtt3 | 0.0446 s | 0.1282 s | 2.87× |
| TopoToolbox | 0.0541 s | 0.1067 s | 1.97× |

Both tools show the same ordering, and matched algorithm against algorithm the
two implementations are close. **Comparing across families is the mistake to
avoid** — an early version of the benchmark suite ran lsdtt3's default (heap)
against TopoToolbox's default (reconstruction) and reported lsdtt3 2.9× slower,
which measured the choice of algorithm rather than the quality of either
implementation. See `docs/developer_notes/testing_and_benchmarking.md`.

---

## Flat resolution: Barnes et al. (2014b)

The default, `fill_flat_resolution_choice: inline_gradient`, imposes
`fill_min_slope` during the Priority-Flood fill, as described above. Two opt-in
options replace the imposed ramp on flats with the flat-resolution method of
Barnes, Lehman and Mulla (2014b):

| `fill_flat_resolution_choice` | what it does on flats left by a `fill_min_slope: 0` fill |
|---|---|
| `inline_gradient` (default) | nothing separate: the Priority-Flood fill imposes `fill_min_slope` as it floods |
| `barnes2014b_gradient` | the Barnes mask scaled to a physical gradient set by `fill_min_slope` |
| `barnes2014b_float_step` | the Barnes mask applied as float32 steps: RichDEM `ResolveFlats` on float32 in its first pass |

The option is a parameter of `lsdtt-raster-preprocessing` and of the routing
CLIs that run a fill (see their CLI reference pages). With the default, the
output is unchanged, bit for bit. The implementation is in
`src/cpu/lsdtt_raster_flat_resolution_cpu.{hpp,cpp}`.

### Why the default gives parallel flow on flats

The inline ramp raises each flat cell by `fill_min_slope × cell size` (times
√2 on a diagonal) for each step it takes from the outlet. On an 8-neighbour grid
every cell of a cross-valley column is the same number of steps from the
outlet, so the ramp has no cross-valley component. Flow on the flat then runs in
parallel lines, with no convergence towards the outlet.

The valley test case is 11 × 24 cells at 10 m, with a floor at 100 m, walls at
110 m, and the outlet at row 5, column 0. D8 accumulation, in cells, across the
valley at column 12, rows 1-9:

| method | rows 1-9 |
|---|---|
| default `fill_pits`, `fill_min_slope` 1e-4 | 24 12 12 12 12 12 12 12 24 |
| `barnes2014b_float_step` | 2 3 4 5 **104** 5 4 3 2 |

With the Barnes mask, the flow from the whole flat converges on the centre line
(row 5) and runs to the outlet. The default spreads the same area evenly across
the floor. This is test F5 in `tests/cpp/test_flat_resolution.cpp`.

### The Barnes 2014b mask

The mask is computed on a depression-free DEM (the output of
`fill_pits_reconstruction`). Each valid cell is classified as follows.

* **Draining** if it has a strictly lower valid neighbour, or if it is an
  **exit**. An exit is a cell on the raster border, or a cell touching interior
  nodata. This is the same exit rule as the fill, so a flat next to an interior
  hole drains into it.
* **NO_FLOW** otherwise.
* A **low edge** is a draining cell with an equal-elevation NO_FLOW neighbour.
* A **high edge** is a NO_FLOW cell with a strictly higher neighbour.

A flat region is the set of 8-connected, equal-elevation cells connected to a
low edge. Within each region two breadth-first distances are measured over
NO_FLOW cells:

* `a`, the distance from the high edges (`a = 1` on a high edge), and `H`, the
  largest `a` in the region;
* `t`, the distance from the low edges (`t = 1` on a low edge).

The mask of a NO_FLOW cell is `2t + H - a` in a region with high edges, and
`2t` otherwise. Low edges get 2, which is the RichDEM convention. Every NO_FLOW
cell has a neighbour at least one mask unit lower, so every cell drains. The
step onto a low edge can be up to `H + 1` units.

Two rules follow from the reference implementation:

* **Border cells are never raised**, whatever their mask. RichDEM's
  `ResolveFlats` never modifies them, and the same rule is used here so that the
  two agree. Cells next to interior nodata are raised.
* The distances and maxima are taken over a fixed graph, so the mask does not
  depend on scan order, queue order or thread count.

Reference: Barnes, R., Lehman, C., Mulla, D. (2014). An efficient assignment of
drainage direction over flat surfaces in raster digital elevation models.
*Computers & Geosciences* 62, 128-135. doi:10.1016/j.cageo.2013.01.009. The code
calls it `barnes2014b`, to keep it apart from the Priority-Flood paper
(`barnes2014`).

### Option 1: `barnes2014b_float_step`

Each flat cell is raised by its mask in float32 steps: `nextafter` applied
`mask` times. The steps are computed in O(1) by `step_up_float32`, which is
bit-identical to the loop. This is RichDEM's `ResolveFlats` on a float32 raster.
On flats at |z| ≥ 1 m the first round is bit-for-bit equal to RichDEM, border
cells included, on the hand-built fixtures (oracle O1). On the fetched alpine and
lakeland small DEMs (oracle O2), the only cells that differ from RichDEM are those
that the later rounds changed, and none is below RichDEM's value.

**The slopes are not physical.** A float-step drop is one or a few ULPs of the
elevation, so the gradient across a resolved flat is about ULP(z) / cell size.
It depends on elevation and jumps by a factor of two at 512, 1024, 2048 and
4096 m. On a 30 m grid (D8 diagonal) it is about 2e-7 m/m at 100 m, 1.4e-6 m/m
at 1000 m and 2.3e-5 m/m at 8848 m. The consequences:

* lsdtt3's D8 router compares elevations in double precision, so it routes any
  positive drop, even at 0 m.
* Most physics consumers floor slope at 1e-5 to 1e-3 (channel extent, GIUH,
  sediment routing, the hydraulic graph). To them a float-step flat is still flat.
* Any consumer that forms the slope in float32, or takes its logarithm, sees a
  small number that depends on elevation. At 0 m it sees zero.

**The elevation change can be large.** The raise is `mask × ULP`. The float32
ULP is 7.6e-6 m on [64, 128) m, 1.22e-4 m on [1024, 2048) m and 4.88e-4 m on
[4096, 8192) m. For a 30 m flat the raise is:

| flat (30 m grid) | at 100 m | at 1500 m | at 3700 m | at 4700 m |
|---|---|---|---|---|
| floodplain, mask 350 (5 km × 1 km) | 2.7 mm | 4.3 cm | 8.5 cm | 17 cm |
| lake, mask 720 (10 km × 3 km) | 5.5 mm | 8.8 cm | 18 cm | 35 cm |
| salar or high lake, mask 5000-7500 | 4-6 cm | 0.6-0.9 m | 1.2-1.8 m | 2.4-3.7 m |
| lidar 1 m floodplain, mask 4150 | 3.2 cm | 0.5 m | 1.0 m | 2.0 m |

On most flats the change is millimetres to centimetres. On large high flats it
can be metres, and a rim a few ULPs above a flat can be overtopped. The rounds
described below repair those rims.

**The near-zero guard.** The float32 ULP shrinks towards zero elevation. At 0 m
it is 1.4e-45 m, a denormal number, so a flat at 0 m would get steps that
float32 consumers see as zero slope (TWI becomes +inf and the MFD drops the
receiver). For flats with |z| < 1 m each mask unit is therefore a fixed
2^-23 m (1.19e-7 m): `z' = float(z + mask × 2^-23)`. Above 1 m the native ULP
is already at least 2^-23 m, so the step is continuous at |z| = 1 m. The guard
keeps the Barnes ordering exactly and is always on; it is not a parameter. It
is bit-equal to RichDEM everywhere except on those flats. The step stays
distinct only while z' < 2 m, and the code throws beyond that bound (a mask of
about 8.4 million at |z| < 1 m). On a 30 m grid the smallest float32 slope it
leaves is 2^-23 / (√2 × 30 m) = 2.8e-9 m/m.

In `lsdtt-raster-preprocessing` and `lsdtt-basin-stats`, `min_elevation`
(default 0) sets cells at or below that elevation to nodata before the fill.
By default, flats at or below 0 m are therefore removed, and the guard acts on
flats between 0 and 1 m. Lowering `min_elevation` below 0 brings sea-level
flats into the fill, and then the guard acts on them.

### Option 2: `barnes2014b_gradient`

The same mask and the same raised cells, but one mask unit is a physical length:

    u = fill_min_slope × √2 × L

where √2 L is the router's own diagonal path length. Every raised NO_FLOW cell
has a neighbour one unit lower at most √2 L away, so every D8 step on an
uncapped flat has slope at least `fill_min_slope`, in the router's double
arithmetic. Inside a flat the steps are one to three units, so the slopes run
from `fill_min_slope` to 3√2 × `fill_min_slope` on a cardinal step. The unit
includes the √2 because a unit of `fill_min_slope × L` would guarantee only
0.71 × `fill_min_slope` on a diagonal drop.

Four measures make the gradient hold on real float32 data.

1. **Quantised unit.** Rounding `z_f + m × u` to float32 cell by cell can
   shorten a drop. Measured on synthetic flats during development: one-unit drops
   of 0.978 u at 1500-4000 m on a 30 m grid, and, on 1 m cells at 4000 m (u =
   1.4e-4 m, below the float32 ULP of 2.4e-4 m), 1262 of 3000 one-unit drops were
   zero.
   So for each flat, `q` is the float32 spacing (a power of two) at the top of its
   raise, and `u_q = ceil(u / q) × q`. Where `u < q` the gradient is coarser than
   requested (`u_q = q`), but it is never below `fill_min_slope`.
2. **Rim cap.** A raised cell that reaches or passes a formerly higher neighbour
   leaves that neighbour with no lower cell. Each flat's unit is therefore capped
   at the minimum, over raised cells `c` and higher neighbours `h`, of
   `(z_h - z_f) / (mask(c) + 1)`, quantised down to a multiple of `q`. A capped
   flat keeps the Barnes ordering but falls below `fill_min_slope`. That is
   reported, not hidden.
3. **Verification and fallback.** Each gradient-resolved flat is checked on the
   float32 output: every NO_FLOW cell needs a strictly lower final neighbour, and
   no raised cell may reach a formerly higher neighbour. A flat that fails, or
   whose capped unit quantises to zero, is raised with float steps instead. That
   fallback is counted in the report.
4. **Across rounds.** A later round can change the neighbourhood of a gradient
   cell (see the next section), so the gradient guarantee holds within the round
   that resolved the flat (see "What is guaranteed").

The raise is larger than the inline ramp's. Its upper bound per flat is
`(2T + H) × √2 × fill_min_slope × L`, where `T` and `H` are the towards and
away distances of the flat, against about `fill_min_slope × path length` for the
inline ramp. The measured maxima on the real DEMs are in the table below.

### Real terrain: one pass is not routable

On real float32 DEMs a single Barnes pass leaves cells without a lower neighbour.
They are rim cells that a flat has overtopped by a few ULPs or millimetres, and
low edges whose neighbouring lower flat was raised past them. The counts below
are the cells left without descent after one pass, from
`lsdtt-raster-preprocessing` with `fill_min_slope` 1e-4:

| fixture | gradient | float step |
|---|---|---|
| alpine_500 | 27 | 28 |
| alpine_1000 | 20 | 21 |
| lakeland_500 | 918 | 918 |
| lakeland_1000 | 2469 | 2475 |

The D8 router then fails with `could not assign a receiver`. Synthetic fixtures,
whose rims sit metres above their flats, did not show this. RichDEM's
`ResolveFlats` also leaves cells without a lower neighbour on the fetched small
DEMs (oracle O2 counts 1 on alpine_small and 41 on lakeland_small).

The resolver runs up to four rounds:

1. **Round 1** is the Barnes pass described above.
2. If any cell still has no lower neighbour, the **whole surface is refilled**
   with `fill_pits_reconstruction`. That only raises pits and the cells that
   drain into them, to their spill level. A cell that already drains is not
   changed.
3. The flats are **resolved again** on the refilled surface, with the same
   method, unit, cap, verification and guard. From round 2, **low edges are not
   raised**; only NO_FLOW cells (mask 4 or more) are. Raising a low edge that
   belongs to an already resolved flat would strand that flat again.
4. Repeat from step 2 until no cell lacks a lower neighbour, or four rounds have
   run.
5. **Final pass.** If cells remain after the last round, one priority-flood pass,
   seeded at the exits, raises only a cell that is not above the cell it is
   reached from. It raises it by one float32 step (2^-23 m below 1 m, with the
   guard). A cell with a strictly descending path to an exit is not changed.

No cell is ever lowered, by any round or by the final pass. The resolver checks
this and throws if it does not hold. Measured by `lsdtt-raster-preprocessing` on
the real DEMs (one run per row, deterministic counts):

| fixture | option | rounds (cells without descent after each round) | changed after round 1 (max extra raise) | final pass cells (max raise) | max raise | capped / float-step fallback | min effective gradient (m/m) | cells below min slope / NO_FLOW |
|---|---|---|---|---|---|---|---|---|
| alpine_500 | gradient | 4 (27, 16, 9, 2) | 60 (0.042 m) | 2 (0.98 mm) | 0.220 m | 90 / 19 | 1.15e-5 | 835 / 8260 |
| alpine_500 | float step | 4 (28, 15, 9, 2) | 62 (0.042 m) | 2 (0.98 mm) | 0.096 m | - | - | - |
| alpine_1000 | gradient | 2 (20, 0) | 20 (0.028 m) | 0 | 0.352 m | 261 / 56 | 1.15e-5 | 4784 / 24338 |
| alpine_1000 | float step | 2 (21, 0) | 25 (0.028 m) | 0 | 0.096 m | - | - | - |
| lakeland_500 | gradient | 4 (918, 388, 489, 374) | 10199 (0.020 m) | 5898 (3.1 mm) | 0.159 m | 1222 / 314 | 7.2e-7 | 21775 / 66374 |
| lakeland_500 | float step | 4 (918, 386, 488, 373) | 10199 (0.012 m) | 5912 (5.0 mm) | 0.014 m | - | - | - |
| lakeland_1000 | gradient | 4 (2469, 403, 642, 480) | 11209 (0.046 m) | 5618 (2.6 mm) | 0.192 m | 4969 / 1722 | 7.2e-7 | 67029 / 275305 |
| lakeland_1000 | float step | 4 (2475, 406, 642, 487) | 11260 (0.043 m) | 5629 (2.6 mm) | 0.044 m | - | - | - |

The final pass ran on alpine_500 and on both lakeland DEMs. On alpine_1000 the
rounds alone were enough.

### What is guaranteed

**Routable.** With the defaults (four rounds and the final pass), both barnes
options return a surface in which no valid cell, other than an exit, lacks a
strictly lower valid neighbour. The D8 router accepts it. This is tested on the
hermetic real-terrain crop and on the fetched real DEMs through
`lsdtt-drainage-analysis`, which exits 0.

**Gradient, within a round.** Every D8 step of every uncapped gradient-resolved
flat has slope at least `fill_min_slope`, in the round that resolved it. The
guarantee is not complete overall:

* capped flats fall below `fill_min_slope` by design;
* a later round, or the final pass, can raise the receiver of a gradient cell
  and reduce its slope.

So the report counts the gradient cells below `fill_min_slope` (`Cells below min
slope`) and gives the minimum gradient measured on the final surface. On the real
DEMs above, 10 % to 33 % of the NO_FLOW cells end below `fill_min_slope`, and the
minimum measured gradient is 1.15e-5 m/m on the alpine DEMs and 7.2e-7 m/m on
the lakeland DEMs. In short, the gradient is at least `fill_min_slope` where the
rim allows it, and the cells where it does not are counted.

**Float step, elevations.** Elevation changes are as in the table in Option 1:
millimetres to centimetres on most flats, metres on large high flats.

**Interior nodata.** A flat next to an interior hole drains into it, and the
hole's neighbours are never stranded. This assumes `interior_nodata_is_outlet`
is true, which is the flow-routing default and what all the CLIs set. With it
false, a flat next to a hole can be a sink.

### Relation to RichDEM

* On flats at |z| ≥ 1 m, round 1 of `barnes2014b_float_step` is bit-for-bit
  RichDEM `ResolveFlats` on float32. On guarded flats (|z| < 1 m) it differs by
  design.
* The later rounds and the final pass do work that RichDEM does not do. On the
  lakeland 500 and 1000 DEMs they change 10,199 to 11,260 cells after round 1 (the
  table above). The final output therefore equals RichDEM only where no later
  round changed a cell.
* RichDEM's "Cells inappropriately raised above surrounding terrain" counts
  ordered neighbour pairs (c, n) with z(n) > z(c) before and z'(c) ≥ z'(n) after,
  so the count uses `≥`. Where two raised flats touch, RichDEM's printed count
  depends on its loop order. On one real crop it prints 367 where a recount from
  its own output gives 364. The round-1 count of lsdtt3 is compared with the
  recount.

### Large MFD exponents (underflow fallback)

The MFD weight of each receiver is a power of its slope, computed in float32.
On gentle surfaces with a large exponent the float32 power underflows (a slope of
about 1e-7 with an exponent of 6.5 or more). The weights can then sum to zero, in
which case the cell silently becomes a sink, or to a subnormal number, or
`1 / total` overflows. This failure is not specific to the float step. With the
default fill at exponent 8 the earlier code already failed on the alpine and
lakeland small fixtures, and it left silent sinks on the lakeland fixture.

The fix applies only when the summed weight is not a normal, finite float. Then
the cell's weights are recomputed as `(S_k / S_max)^p` and normalised. Normalisation
removes the common factor, so the weights are the same. Every other cell is
computed exactly as before. Outputs are bitwise unchanged at exponents 1 and 4 on
the valley, alpine small and lakeland small fixtures, and at exponent 8 on the
valley. The runs that failed before (exponent 8 on the two real fixtures) now
complete.

---

## s_vector rationale (D8 only)

`s_vector` is the key D8 indexing feature in `lsdtt3`:

- Stores a preorder DFS layout of each drainage tree.
- Each cell `i` gets `start = s_vector_start[i]` and `size = s_vector_size[i]`.
- The upslope set of `i` is the contiguous interval
  `[start, start + size)` in `s_vector`.
- `is_upslope_of(a, b)` becomes an interval-membership check (O(1)).

Why this matters:

- Repeated upslope queries (chi, basin summaries, thresholded source finding,
  many cosmogenic operations) avoid repeated graph searches.
- Memory is paid once at build time; query-time cost is predictable.

### What it costs at build time

The S-vector is not free, and it is worth knowing how much of a D8 build it
accounts for. Measured on a 1000×1000 filled tile, single-threaded, by toggling
each `lsdtt_flowinfo_options` flag:

| configuration | build time |
|---|---|
| all on (what the CLIs request) | 186.4 ms |
| without the acyclicity check | 190.0 ms (free) |
| without the topological order | 178.3 ms |
| **without the S-vector** | **80.1 ms** |
| receivers only | 52.4 ms |

**The S-vector is ~106 ms of the 186 ms — more than everything else combined.**

Because `get_drainage_area_raster()` reads the precomputed `s_vector_size`
intervals, the drainage-area call itself then costs ~1.4 ms: the accumulation
work was already paid for during the build. Accumulating over the topological
order instead gives the identical raster (0 differing cells of 1,000,000) with
a cheaper total:

| route | build | accumulate | total |
|---|---|---|---|
| `s_vector` subtree sizes | 0.1934 s | 0.0016 s | 0.1951 s |
| topological order + `accumulate()` | 0.0941 s | 0.0400 s | **0.1341 s** |

**This is deliberately not changed.** If drainage area were the only product,
the second route would be 1.45× faster — but the S-vector is what makes
`is_upslope_of` an O(1) interval test, and chi analysis, basin summaries,
source finding and the cosmogenic routines all depend on it. Building it once
during routing is cheaper than reconstructing upslope sets on demand later.

The consequence for benchmarking is that lsdtt3's D8 accumulation timing
includes work no other tool performs. That asymmetry is recorded in the
benchmark case notes rather than engineered away.

For D-infinity/MFD, upslope sets overlap fractionally, so an equivalent
contiguous interval index is not available in the same form.

---

## 64-bit indexing policy and rationale

`lsdtt3` uses 64-bit flat indices in flow-routing objects:

- Removes structural 32-bit index ceilings for very large rasters.
- Keeps indexing semantics uniform across analyses.
- Trades some per-cell bytes for scalability beyond the practical limits seen in
  older 32-bit-indexed designs.

See [`overview.md`](../what_is_lsdtt3/overview.md) (Memory use and Very large DEMs) for measured and derived
capacity/memory discussion. The use of a 64 bit integer has a memory cost, but `lsdtt3`
is aimed at processing very large DEMs, so we don't want to run out of indices!! 

---

## Memory: what lsdtt3 changed from LSDTopoTools2

Both libraries build the same D8 topology. They store it very differently, and
the difference is mostly a matter of what is **kept** versus what is **recomputed
or never needed**.

Writing `n` for grid cells and `m` for valid (non-nodata) cells:

**LSDTopoTools2** (`LSDFlowInfo`, protected members):

| structure | type | extent | bytes |
|---|---|---|---|
| `NodeIndex` | `Array2D<int>` | n | 4n |
| `FlowDirection` | `Array2D<int>` | n | 4n |
| `FlowLengthCode` | `Array2D<int>` | n | 4n |
| `RowIndex`, `ColIndex` | `vector<int>` x2 | m | 8m |
| `NDonorsVector`, `DeltaVector`, `DonorStackVector` | `vector<int>` x3 | m | 12m |
| `ReceiverVector` | `vector<int>` | m | 4m |
| `SVector`, `SVectorIndex`, `NContributingNodes`, `BLBasinVector` | `vector<int>` x4 | m | 16m |

roughly **12n + 40m bytes**, or about **52 bytes per cell** when most cells are
valid.

**lsdtt3** (`lsdtt_flowinfo`, private members):

| structure | type | extent | bytes | when |
|---|---|---|---|---|
| `mask_` | `uint8` | n | 1n | always |
| `receiver_` | `int64` | n | 8n | always |
| `receiver_dir_` | `uint8` | n | 1n | always |
| `donor_count_` | `uint8` | n | 1n | always |
| `topological_order_` | `int64` | m | 8m | released if not requested |
| `s_vector_` + `s_vector_start_` + `s_vector_size_` | `int64`/`int64`/`uint64` | m + n + n | 8m + 16n | only if built |

roughly **11 bytes per cell** for the core graph, **19** with a retained
topological order, and **43** with the O(1) upslope index.

### Where the savings come from

- **The flat index *is* the node id.** LSDTopoTools2 maintains a `(row, col) ->
  NodeIndex` map plus `RowIndex`/`ColIndex` to invert it. `lsdtt3` addresses
  cells by flat index throughout, so `NodeIndex`, `RowIndex` and `ColIndex` -
  three arrays - simply do not exist.
- **`FlowLengthCode` is derivable.** LSDTopoTools2 stores a full `int` grid whose
  only values are 1 (cardinal) and 2 (diagonal). `lsdtt3` recovers that from the
  direction code, which it already has.
- **The donor CSR is never materialised.** `NDonorsVector` + `DeltaVector` +
  `DonorStackVector` - 12 bytes per cell - collapse to a single `uint8`
  `donor_count_`. A D8 cell can have at most 8 donors, so a byte is sufficient,
  and Kahn's algorithm needs only the count.
- **Narrow types where the range allows.** Direction codes and donor counts are
  `uint8`, not `int`.
- **Optional structures are optional.** The topological order is released after
  build if the caller did not ask for it, and `s_vector` is built only on demand.

### What it costs

64-bit indexing. `receiver_` is `int64` over the full grid - 8 bytes per cell
where LSDTopoTools2 spends 4 per valid node. That is a deliberate trade: it
removes the 32-bit ceiling on raster size (see below) at a cost the eliminations
above more than repay: the overall memory usage in `lsdtt3` is less than `LSDTopoTools2`.

### Net effect

The per-structure accounting above resolves to a single headline figure and a
caveat about how it was obtained. Rather than restate them here, both live in
[`overview.md`](../what_is_lsdtt3/overview.md) under *Memory use*, which is
where a reader looking for the comparison will go first.

For cross-implementation comparisons against RichDEM, TopoToolbox and Whitebox
— timings, fixtures, machine and the caveats that go with each — see
[`testing_and_benchmarking.md`](../developer_notes/testing_and_benchmarking.md),
which is the durable record. If you need publication-grade numbers, run the
repository benchmark harness and report the environment and fixture details
alongside them.

---

## Depression hierarchy

`lsdtt_depression_hierarchy_cpu` is lsdtt3's port of the Barnes et al. (2020)
depression hierarchy: a tree in which every closed depression is a node, leaves
are the innermost hollows, and a parent is the meta-depression two children merge
into once water rises to their shared saddle. It underpins
`lsdtt-internal-drainage`'s `depression_hierarchy` detector and the fill-spill-merge
routines in `lsdtt-raster-preprocessing`.

It is validated against the upstream reference implementation by
`unit_depression_hierarchy_oracle`, which compares depression count and the
per-depression fields (`parent`, `lchild`, `rchild`, `pit_elev`, `out_elev`,
`cell_count`, `ocean_parent`) across 25 committed test DEMs. Provenance of those
fixtures, including the upstream commits they pin, is in
`tests/expected/depression_hierarchy/PROVENANCE.md`.

### Choosing a detector: what the hierarchy adds, and what it risks

`lsdtt-internal-drainage` offers two sink detectors through
`sink_detection_algorithm_choice`. **The default is `fill_difference`.**

**`fill_difference`** takes the difference between the filled and raw DEM and
treats each connected patch that filling would raise as one sink, keeping its
deepest cell as the pit. It is the method the published Greenland moulin work
used. It is fast, robust on rough surfaces, and has no notion of nesting: a lake
containing several hollows is one sink.

**`depression_hierarchy`** builds the full tree and adds three things
`fill_difference` cannot give you:

- **Nesting.** A lake containing several hollows becomes a parent with children,
  reported through `parent_sink_id` and `hierarchy_level`. With
  `sink_leaves_only: true` you get the innermost hollows; with `false` you also
  get the enclosing depressions.
- **True depression volume**, accumulated over the tree, rather than a summed
  fill difference.
- **Spill relationships between depressions** — which depression overflows into
  which, and which reach the ocean directly.

Use it where that structure is the point: playas, karst, quarry benches, closed
tectonic basins — terrain with real depressions inside depressions.

#### The risk: it can report nothing at all

On a rough surface, **every innermost hollow is tiny**. Combining
`depression_hierarchy` with `sink_leaves_only: true` and a
`min_pixels_in_patch` threshold can therefore filter out *every* sink and
return an empty result, with no error.

This is not hypothetical. On an 18 m Greenland ablation-zone DEM (17.4 Mpx) with
the configuration that reproduces the published work
(`min_pixels_in_patch: 1000`):

| detector | sinks found | runtime |
|---|---|---|
| `fill_difference` | **143** | 5.4 s |
| `depression_hierarchy` | **0** | 24.6 s |

The hierarchy is not wrong there — it finds 468,878 leaf depressions. But every
one is a metres-across dimple in the ice, none reaches 1000 pixels, and the size
filter removes them all. The fill-difference patches, by contrast, are the
filled lakes the analysis is actually about.

So before using `depression_hierarchy` on a noisy surface, either smooth the DEM
first, or filter on depth or volume rather than on pixel count, or set
`sink_leaves_only: false` so enclosing depressions survive the threshold. If a
run reports zero sinks, this is the first thing to check.

Note that `sink_leaves_only` applies **only** to `depression_hierarchy`. Under
`fill_difference` there is no tree to take the leaves of, and setting it to
`false` prints a warning and is otherwise ignored.

### Defining where water leaves the grid

The hierarchy needs to know which cells are "ocean" — the places water leaves the
model. There are two ways to say it, and they can be combined:

- **Boundary conditions** (`boundary_conditions`, default `"oooo"`). Valid cells
  that drain off an open grid edge become ocean seeds. This is lsdtt3's native
  convention and matches the rest of the flow-routing layer.
- **An explicit ocean definition**, for DEMs where the sea is encoded in the data
  rather than at the edge. **These are library options, not `.param` keys** —
  they are reachable from the Python bindings (see
  [`drainage`](../python_api/drainage.md)) and from C++, but
  `lsdtt-internal-drainage` does not declare them, so putting them in a
  parameter file has no effect:
  - `nodata_is_ocean` — treat NODATA cells as ocean. Common where a coastline or
    ice margin has been masked out.
  - `ocean_level` — treat valid cells at exactly this elevation as ocean. Common
    where the sea is burned in as a constant.
  - `ocean_from_boundary` — set to `false` to use only the explicit rules.

Seeds are the **union** of whichever rules are enabled. If none of them selects
any cell, `build_depression_hierarchy()` throws rather than guessing.

This matters beyond convenience: a DEM whose sea is a nodata region cannot be
described by boundary conditions at all, because no valid cell touches the grid
edge. Upstream Barnes takes an ocean mask for the same reason.

### Depressions that spill at the same height

Imagine slowly raising the water level across a landscape. Each depression fills
until it reaches its lowest rim, then spills over and joins whatever is next
door. The hierarchy records the order in which those joins happen: the depression
that spills first sits lowest in the tree.

That ordering is obvious when every rim is at a different height. **It is not
obvious when two rims are at the same height.** If two depressions would spill at
exactly the same moment, which one joins first? Nothing in the landscape decides
it — the answer has to come from the software.

Equal rim heights are common rather than exotic. Any DEM stored as whole-number
elevations produces them constantly, and so does any flat saddle.

**What the original code does.** It sorts the spill points by height and takes
them in that order. When two have the same height, it leaves their relative order
to whatever the sorting routine happens to do — which in practice depends on the
order items came out of an internal lookup table. That order is not defined by
anything meaningful. It can differ between one computer and another, between
compiler versions, or even if the table happens to be a different size.

**What lsdtt3 does.** When two spill points are at the same height, lsdtt3
compares the identifying numbers of the depressions involved, and then the
position of the spill cell itself, and uses that to decide. These are arbitrary
tie-breakers, but they are *fixed*: the same DEM gives the same tree on any
machine, every time. For a tool producing published results, reproducibility is
worth more than matching another implementation's accident.

**What this costs.** On a DEM with equal-height rims, lsdtt3 and the original
build **different trees — and both are correct**. They disagree about the order
of events that genuinely happened simultaneously.

This shows up in testing. Of the 25 reference DEMs checked against the original
implementation, 8 produce a different tree. In every one of those 8:

- the **depressions themselves are identical** — same locations, same floor and
  rim elevations, same number of cells in each;
- only the way they are *linked together* differs.

In five of the eight, **every spill point in the DEM is at the same elevation**,
so there is no correct nesting to find at all; any tree is as good as any other.

Because the depressions agree even when the tree does not, the test on those 8
DEMs checks the depressions and ignores the tree, rather than skipping the DEM.
That way a genuine mistake in finding depressions would still be caught.

### Multipart footprints

A sink's cells need not be contiguous. When `min_pixels_in_patch` filters out a
small hollow, its cells are attributed to the nearest surviving ancestor, and
that ancestor's footprint can end up in two or more pieces separated by its other
children. `sink_depressions_to_polygons()` emits **one polygon record per
connected part**, all carrying the same sink attributes, so no terrain is lost.
`test_min_pixels_conserves_terrain` pins the invariant: the union of all reported
footprints must not change when the size filter is applied, because the filter
merges sinks and never removes ground.

---

## Important behavioral differences

### Why RichDEM's D8 differs from lsdtt3 and TopoToolbox

RichDEM should **not** be used as the reference for steepest-descent D8. This is
a genuine definitional difference, not an implementation bug, and it is worth
stating precisely because it has been described loosely elsewhere.

**The mechanism.** `richdem/flowmet/d8_flowdirs.hpp`, `d8_FlowDir()`:

```cpp
T   minimum_elevation = elevations(x,y);
int flowdir           = NO_FLOW;
for(int n=1;n<=8;n++)
  if( elevations(x+d8x[n],y+d8y[n]) < minimum_elevation
      || (elevations(x+d8x[n],y+d8y[n]) == minimum_elevation
            && flowdir>0 && flowdir%2==0 && n%2==1) )
  { minimum_elevation = ...; flowdir = n; }
```

The selection criterion is **lowest elevation**. Distance never enters. True D8
(O'Callaghan and Mark 1984, which both libraries cite) selects on **drop divided
by distance**, so a diagonal neighbour must be `sqrt(2)` times lower to win.
lsdtt3 and TopoToolbox both implement drop-over-distance; RichDEM does not.

**A worked case.** Cardinal south is 1.0 m lower; diagonal south-east is 1.2 m
lower on a 1 m grid:

| | drop | distance | gradient | chosen by |
|---|---|---|---|---|
| cardinal S | 1.0 | 1.0 | **1.000** | lsdtt3, TopoToolbox |
| diagonal SE | 1.2 | 1.414 | 0.849 | **RichDEM** |

RichDEM takes the diagonal because it is lower; steepest descent takes the
cardinal because it is steeper. Verified directly against the RichDEM source at
commit `af61c65`: `d8_FlowDir` returns 6 (diagonal SE) for this configuration.

**Two corrections to earlier descriptions in this repository.** Both were checked
against the source and by running `d8_FlowDir` directly:

1. The tie-break is often described as *preferring diagonals*. It does the
   opposite. In RichDEM's numbering, `n8_diag = {false,false,true,false,true,
   false,true,false,true}` - **odd `n` is cardinal**, even is diagonal - so the
   clause `flowdir%2==0 && n%2==1` switches *from* a diagonal *to* a cardinal at
   equal elevation. RichDEM's own source carries the comment
   `//TODO: What is this modulus stuff for?` beside it.
2. RichDEM is often said to **fail a planar-ramp test**, producing diagonal
   stripes. It does not. On a surface tilted purely north-south,
   `d8_FlowDir` returns 7 - cardinal south, the correct answer - because the
   tie-break above pulls it off the diagonal.

Neither correction weakens the conclusion. The divergence does not depend on ties
at all: it appears wherever the lowest neighbour is diagonal and a cardinal
neighbour is steeper, which is common on real terrain. Measured consequences on
`alpine_small`:

| | peak D8 accumulation |
|---|---|
| lsdtt3 | 41,437 |
| TopoToolbox | **41,437** |
| RichDEM | 40,446 |

with a D8 rank correlation of 0.51 against RichDEM, rising to 0.74 for D-infinity
and 0.99 for MFD - dispersion washes the single-receiver choice out.

RichDEM remains a valid reference for **depression handling**, where lsdtt3 tests
against it directly. It is only D8 flow direction that is incomparable.

### TopoToolbox as a cross-check

`lsdtt3` testing docs report that TopoToolbox's D8 behavior matches the
drop-over-distance steepest-descent rule and is used as an asserted comparison
for D8 accumulation checks.

### Barnes / Lindsay preprocessing papers are upstream of routing

- **Barnes et al. (2014)** and **Lindsay (2016/2020)** in `lsdtt3` apply to
  depression handling (fill/breach), i.e. conditioning the DEM before routing.
- They are not alternative D8 routing equations.
- In low-relief settings, flat-resolution choices (`fill_min_slope`, breach
  policy) can strongly influence subsequent D8 network realization.

---

## Primary references

- O'Callaghan, J.F., and Mark, D.M. (1984). The extraction of drainage networks
  from digital elevation data. *Computer Vision, Graphics, and Image Processing*
  28(3), 323-344.
- Tarboton, D.G. (1997). A new method for the determination of flow directions
  and upslope areas in grid digital elevation models. *Water Resources
  Research* 33(2), 309-319.
- Quinn, P., Beven, K., Chevallier, P., and Planchon, O. (1991). The prediction
  of hillslope flow paths for distributed hydrological modelling using digital
  terrain models. *Hydrological Processes* 5(1), 59-79.
- Freeman, T.G. (1991). Calculating catchment area with divergent flow based on
  a regular grid. *Computers & Geosciences* 17(3), 413-422.
- Barnes, R., Lehman, C., and Mulla, D. (2014). Priority-flood: An optimal
  depression-filling and watershed-labeling algorithm for digital elevation
  models. *Computers & Geosciences* 62, 117-127.
- Barnes, R., Lehman, C., and Mulla, D. (2014). An efficient assignment of
  drainage direction over flat surfaces in raster digital elevation models.
  *Computers & Geosciences* 62, 128-135. doi:10.1016/j.cageo.2013.01.009. The
  flat-resolution mask behind `barnes2014b_gradient` and `barnes2014b_float_step`.
- Lindsay, J.B. (2016). Efficient hybrid breaching-filling sink removal methods
  for flow path enforcement in digital elevation models. *Hydrological
  Processes* 30(6), 846-857.
- Lindsay, J.B. (2020). Pit-centric depression removal methods.
  *Geomorphometry 2020 Conference Proceedings*.
- Braun, J., and Willett, S.D. (2013). A very efficient O(n), implicit and
  parallel method to solve the stream power equation governing fluvial
  incision and landscape evolution. *Geomorphology* 180-181, 170-179.
- Anand, S.K., Hooshyar, M., and Porporato, A. (2020). Linear layout of
  multiple flow-direction networks for landscape-evolution simulations.
  *Environmental Modelling & Software* 133, 104804.
- Kahn, A.B. (1962). Topological sorting of large networks.
  *Communications of the ACM* 5(11), 558-562. The ready-queue topological sort
  used to order both the D8 tree and the multi-receiver DAG.
- Gailleton, B., Steer, P., Davy, P., Schwanghart, W., and Bernard, T. (2024).
  GraphFlood 1.0: an efficient algorithm to approximate 2D hydrodynamics for
  landscape evolution models. *Earth Surface Dynamics* 12, 1295-1313. The
  graph-based hydraulic approach behind `lsdtt_hydraulic_graph_cpu`.

For implementation specifics, also see:

- `docs/developer_notes/using_lsdtt_flowinfo.md`
- `docs/cli_reference/cli_reference_lsdtt_drainage_analysis.md`
- `docs/developer_notes/testing_and_benchmarking.md`
