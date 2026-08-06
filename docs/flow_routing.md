# Flow Routing in lsdtt3

This page explains how flow routing is implemented in `lsdtt3`, what each
method is for, and where `lsdtt3` deliberately differs from common references.

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

## Methods and references

### D8

- Single-flow-direction routing after O'Callaghan and Mark (1984): each cell
  routes to one steepest-descending neighbour.
- `lsdtt3` uses drop divided by distance (cardinal distance = 1, diagonal =
  sqrt(2)), consistent with the D8 definition.

### D-infinity

- Tarboton (1997) triangular-facet method.
- Each cell routes to up to two receivers with angle-based partitioning.

### MFD

- Quinn et al. (1991) / Freeman (1991) multiple-flow-direction routing.
- Flow is partitioned across all downslope neighbours, weighted by
  `slope^exponent` (`mfd_exponent` in CLI/config).

### DAG traversal for multi-receiver routing

- `lsdtt_flowgraph_cpu` builds a receiver DAG and traverses it in topological
  order (Kahn-style ready queue), following the multi-receiver layout idea used
  in Anand et al. (2020).

### Hydraulic graph routing (Gailleton lineage)

- `lsdtt_hydraulic_graph_cpu` is a graph-based hydraulic solver that uses a
  fixed routing graph and updates water state over time.
- In `lsdtt3` documentation and code comments this is explicitly credited as
  based on Gailleton et al. (2024) concepts.

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

For D-infinity/MFD, upslope sets overlap fractionally, so an equivalent
contiguous interval index is not available in the same form.

---

## 64-bit indexing policy and rationale

`lsdtt3` uses 64-bit flat indices in flow-routing objects:

- Removes structural 32-bit index ceilings for very large rasters.
- Keeps indexing semantics uniform across analyses.
- Trades some per-cell bytes for scalability beyond the practical limits seen in
  older 32-bit-indexed designs.

See `docs/overview.md` (Memory use and Very large DEMs) for measured and derived
capacity/memory discussion.

---

## Memory and performance comparisons (what is currently evidenced)

From the benchmark and overview documentation:

- **`lsdtt3` vs LSDTopoTools2 flow object**: about 20% lower memory for
  equivalent D8 topology in the documented 10k x 10k example, with a larger
  reduction when the O(1) upslope index is not needed.
- **`lsdtt3` vs RichDEM / TopoToolbox (routing behavior)**: comparison is
  strongest on algorithmic agreement/correctness framing (below), with timings
  reported in `docs/testing_and_benchmarking.md`.
- **`lsdtt3` vs TopoToolbox (surface-metric benchmark context)**: memory wins
  are documented there too, but those numbers are for surface workflows rather
  than a pure flow-object byte accounting.

If you need publication-grade memory comparisons between routing engines, use
the repository benchmark harness and report environment + fixture details.

---

## Important behavioral differences

### RichDEM as a D8 reference

`lsdtt3` docs explicitly warn that RichDEM should not be used as the definitive
D8 correctness reference for steepest-descent D8. The documented reason is that
RichDEM's D8 implementation (as source-checked in `docs/testing_and_benchmarking.md`)
selects the lowest neighbour with tie behavior that can differ from
drop-over-distance steepest descent.

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

## Braun and Willett context

Braun and Willett (2013) appears in `lsdtt3` in the implicit stream-power
incision solver context (landscape-evolution modeling), not as a replacement D8
router. In other words:

- D8/D-infinity/MFD define how flow is routed.
- Braun-Willett-style implicit updates use that routing/topological structure in
  fluvial incision calculations.

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

For implementation specifics, also see:

- `docs/developer_notes/using_lsdtt_flowinfo.md`
- `docs/cli_reference/cli_reference_lsdtt_drainage_analysis.md`
- `docs/testing_and_benchmarking.md`
