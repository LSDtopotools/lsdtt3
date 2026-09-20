# Using lsdtt_flowinfo

`lsdtt_flowinfo` is the D8 flow routing object in lsdtt3. It builds and stores
the complete receiver graph for a DEM, and provides fast queries for flow
relationships between pixels.

This document explains how to construct and use the object. It is written for
someone adding a new analysis routine, not for someone working on the
infrastructure layer.

---

## The pixel index system

Every pixel in the DEM is identified by a **flat row-major index** — a single
integer that encodes both row and column. This replaces the old NodeIndex.

```
flat_index = row * n_cols + col
```

Flat indices are `std::int64_t` throughout, so flow indexing is 64-bit. There
is deliberately no `flat_index_type` typedef: the concrete type is used
directly so that the width is visible at every call site.

### Converting between flat index and row/col

```cpp
// row, col → flat index
std::int64_t flat = fi.flat_index_unchecked(row, col);

// flat index → row and col
int row = fi.row_from_flat_unchecked(flat);
int col = fi.col_from_flat_unchecked(flat);
```

Use the `_unchecked` variants inside loops where you know the inputs are valid.
The checked variants (`flat_index`, `row_from_flat`, `col_from_flat`) perform
bounds checking and throw on invalid input — use them at the boundary of your
function where you are accepting untrusted input.

### Sentinel values

Two flat index values have special meaning:

| Constant | Value | Meaning |
|---|---|---|
| `lsdtt_flowinfo::nodata_receiver` | -1 | Cell is nodata or invalid |
| `lsdtt_flowinfo::outlet_receiver` | -2 | Cell drains off-grid (baselevel node) |

---

## Construction

Before constructing a `lsdtt_flowinfo`, the DEM must be pit-free and
flat-resolved. Use `fill_pits` or `breach_depressions` from
`lsdtt_raster_flow_preprocessing_cpu.hpp` first.

```cpp
#include "lsdtt_flowinfo_cpu.hpp"

// Options — defaults are fine for most use cases
lsdtt::lsdtt_flowinfo_options opts;
opts.boundary_conditions       = "oooo";  // open on all four sides
opts.interior_nodata_is_outlet = true;    // nodata islands are outlets
opts.require_acyclic_graph     = true;    // validate; throw on a cycle
opts.build_topological_order   = true;    // needed for ordered traversal
opts.build_s_vector            = true;    // needed for upslope queries
opts.n_threads                 = 0;       // 0 = use the OpenMP default
opts.allow_pits_as_outlets     = false;   // unfilled pits are an error

// Construct from a CPU raster
lsdtt::lsdtt_flowinfo fi(
    dem.get_data_array(),   // lsdtt_array2d_f
    dem.get_mask_array(),   // lsdtt_array2d_u8
    dem.get_meta(),         // lsdtt_raster_meta
    opts
);
```

Setting `build_topological_order` or `build_s_vector` to `false` skips that
part of construction; the corresponding queries then throw, and
`has_topological_order()` / `has_s_vector()` report which are available.
Turning either off must not change routing — only which queries are served.

### Forced outlets (internally drained basins)

`opts.forced_outlet_flats` nominates cells as unconditional flow terminations,
given as flat indices. Each listed cell gets `outlet_receiver` and
`no_direction` in preference to any downhill neighbour. This is for genuinely
internally drained basins — moulins, swallow holes, playas, quarries — where
flow leaves the surface at the pit rather than spilling over the depression's
lowest saddle.

**The DEM must have been prepared with `fill_pits_with_sinks()` (or
`breach_depressions_with_sinks()`) over the same cell set.** A plain
`fill_pits()` raises the depression until it spills, so by the time the flow
object sees it the floor is flat and drains outwards — the forced outlet then
becomes a terminus with nothing draining into it. Seeding the fill from the
same cells keeps a monotone descent to each sink. This pairing is not checked
at construction and gets you a silently wrong answer if you break it.

Indices are validated: out-of-range indices, and indices landing on masked or
nodata cells, throw `std::invalid_argument`. Duplicates are ignored. Left
empty (the default), routing is bit-for-bit identical to before the option
existed.

The constructor fills the receiver graph, computes donor counts, finds outlets,
builds the topological order, and builds the s_vector. Construction is
proportional to the number of valid pixels. The object is then read-only — all
query methods are `const`.

### Boundary condition codes

The `boundary_conditions` string has four characters in order: north, east,
south, west.

| Code | Meaning |
|---|---|
| `o` | Open — water exits the domain at this edge |
| `n` | No-flux — water cannot leave through this edge |
| `p` | Periodic — wraps to the opposite edge (north/south and east/west pairs must both be `p`) |
| `b` | Forced baselevel — every cell on this edge is an outlet regardless of elevation |

---

## Querying flow relationships

### Receiver (downstream pixel)

```cpp
std::int64_t rec = fi.receiver_of_unchecked(flat);

if (rec == lsdtt::lsdtt_flowinfo::outlet_receiver)
{
    // this pixel drains off-grid — it is a baselevel node
}
else if (rec == lsdtt::lsdtt_flowinfo::nodata_receiver)
{
    // this pixel is nodata
}
else
{
    int rec_row = fi.row_from_flat_unchecked(rec);
    int rec_col = fi.col_from_flat_unchecked(rec);
}
```

The flow direction code (0–7, using the compass convention below) is also
available if you need the direction rather than the destination pixel:

```
7  0  1
6  x  2
5  4  3
```

```cpp
uint8_t dir = fi.receiver_dir_of_unchecked(flat);
// dir == lsdtt_flowinfo::no_direction (255) for outlet cells
```

### Donor pixels (immediate upstream neighbours)

The number of donors is stored per cell:

```cpp
uint8_t n_donors = fi.donor_count_of_unchecked(flat);
// 0 = headwater cell; max 8
```

To get the actual donor pixels, scan the 8 neighbours and check which ones
point to this cell:

```cpp
constexpr int d_row[8] = { -1, -1,  0,  1,  1,  1,  0, -1 };
constexpr int d_col[8] = {  0,  1,  1,  1,  0, -1, -1, -1 };

int row = fi.row_from_flat_unchecked(flat);
int col = fi.col_from_flat_unchecked(flat);

for (int dir = 0; dir < 8; ++dir)
{
    int nb_row = row + d_row[dir];
    int nb_col = col + d_col[dir];

    if (!fi.in_bounds(nb_row, nb_col)) continue;

    auto nb_flat = fi.flat_index_unchecked(nb_row, nb_col);

    if (fi.is_valid_cell(nb_flat) &&
        fi.receiver_of_unchecked(nb_flat) == flat)
    {
        // nb_flat is a direct donor of flat
    }
}
```

### All upslope pixels

`upslope_nodes` returns a lightweight view (`lsdtt_flat_index_slice`) into the
s_vector. It is not a copy — it is a pointer and size into data already owned
by the `lsdtt_flowinfo` object. It is valid for the lifetime of `fi`.

The slice includes the queried pixel itself as its first element.

```cpp
lsdtt::lsdtt_flat_index_slice upslope = fi.upslope_nodes_unchecked(flat);

for (const auto* it = upslope.begin(); it != upslope.end(); ++it)
{
    std::int64_t up_flat = *it;
    int up_row = fi.row_from_flat_unchecked(up_flat);
    int up_col = fi.col_from_flat_unchecked(up_flat);
    // process upslope pixel...
}
```

Range-based for also works:

```cpp
for (auto up_flat : fi.upslope_nodes_unchecked(flat))
{
    // process up_flat...
}
```

Or use the iteration helper, which accepts a lambda:

```cpp
fi.for_each_upslope_node_unchecked(flat, [&](auto up_flat)
{
    // process up_flat...
});
```

The upslope set is returned in **preorder DFS** order — each node appears before
its donors. This means the queried pixel is always index 0 in the slice, and
you can process nodes from headwaters down to the query point by iterating in
reverse.

### Number of upslope pixels (drainage area proxy)

```cpp
uint32_t n_up = fi.s_vector_size_unchecked(flat);
// includes the pixel itself; proportional to drainage area
```

Multiply by `fi.cell_size() * fi.cell_size()` for drainage area in m².

### Testing whether one pixel is upslope of another

```cpp
bool result = fi.is_upslope_of_unchecked(candidate_flat, target_flat);
```

This is O(1) — it checks whether `candidate_flat`'s position in the s_vector
falls within `target_flat`'s upslope interval, without scanning any pixels.

---

## Baselevel nodes (outlets)

Baselevel nodes are pixels with no receiver — they drain off-grid or to a
boundary. The full list is available directly:

```cpp
const std::vector<std::int64_t>& outlets = fi.outlets();

for (auto outlet_flat : outlets)
{
    int row = fi.row_from_flat_unchecked(outlet_flat);
    int col = fi.col_from_flat_unchecked(outlet_flat);
    // process baselevel node...
}
```

To test a specific pixel:

```cpp
bool is_base = fi.is_outlet_unchecked(flat);
```

The outlet list is in row-major order (top-left to bottom-right). There is
always at least one outlet in a valid flow graph.

---

## Iterating over all valid pixels

To iterate over all valid (non-nodata) pixels in topological order (headwaters
first, outlets last):

```cpp
const std::vector<std::int64_t>& order =
    fi.topological_order();

for (auto flat : order)
{
    // pixels are visited from headwaters toward outlets
    // receiver of flat is guaranteed to appear later in the loop
}
```

This ordering is useful for accumulation algorithms: process each pixel, then
add its contribution to its receiver. Because every donor appears before its
receiver, the receiver always has the correct accumulated value by the time it
is processed.

---

## Object state and diagnostics

```cpp
fi.empty()                // true if not yet constructed
fi.is_acyclic()           // true if the flow graph has no cycles (should always be true after a clean DEM)
fi.n_valid_cells()        // number of non-nodata pixels
fi.n_outlets()            // number of baselevel nodes
fi.n_rows()               // DEM dimensions
fi.n_cols()
fi.cell_size()            // cell size in metres
fi.has_s_vector()         // true if upslope queries are available
fi.has_topological_order() // true if topological_order() is available
```

---

## Checked vs unchecked methods

Most query methods come in two variants:

| Variant | When to use |
|---|---|
| `receiver_of(flat)` | At the boundary of a function, accepting external input. Throws `std::runtime_error` on invalid cell. |
| `receiver_of_unchecked(flat)` | Inside loops where the cell is known to be valid. No bounds check — faster. |

The same pattern applies to `row_from_flat`, `col_from_flat`, `flat_index`,
`upslope_nodes`, `s_vector_size`, `s_vector_start`, `is_outlet`,
`receiver_dir_of`, `donor_count_of`, `is_upslope_of`.

A typical analysis function validates its inputs at entry (using the checked
variants or explicit `is_valid_cell` tests) and then uses `_unchecked` variants
in its inner loops.

---

## Complete example: accumulating in topological order

**In production, do not write this function.** `lsdtt_flowinfo` already
provides `get_drainage_area_raster()`, and
`get_drainage_area_raster_with_injection()` for the case where upstream cells
carry an imposed value. For accumulating an arbitrary per-cell quantity, use
`accumulate()` (returns `std::vector<double>`) or `accumulate_raster()`
(returns `lsdtt_raster_result`).

The example below is kept because it shows the traversal pattern — seed every
valid cell, then push each cell's value into its receiver walking the
topological order from headwaters to outlets — which is the shape of most
flow-routed computations you will write.

```cpp
#include "lsdtt_flowinfo_cpu.hpp"
#include "../utils/lsdtt_array2d.hpp"

lsdtt::lsdtt_array2d_f compute_drainage_area(
    const lsdtt::lsdtt_flowinfo& fi)
{
    const int    nrows     = fi.n_rows();
    const int    ncols     = fi.n_cols();
    const double cell_area = fi.cell_size() * fi.cell_size();
    const float  nodata    = static_cast<float>(fi.meta().no_data_value);

    lsdtt::lsdtt_array2d_f area(nrows, ncols, nodata);

    // Seed each valid cell with one pixel worth of area
    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            if (fi.is_valid_cell(r, c))
            {
                area(r, c) = static_cast<float>(cell_area);
            }
        }
    }

    // Accumulate in topological order (headwaters → outlets)
    for (auto flat : fi.topological_order())
    {
        const auto rec = fi.receiver_of_unchecked(flat);

        if (rec >= 0)   // rec >= 0 means a valid downstream pixel (not outlet, not nodata)
        {
            const int r_rec = fi.row_from_flat_unchecked(rec);
            const int c_rec = fi.col_from_flat_unchecked(rec);
            area(r_rec, c_rec) += area(fi.row_from_flat_unchecked(flat),
                                       fi.col_from_flat_unchecked(flat));
        }
    }

    return area;
}
```

---

## Where to put new analysis functions

New analysis functions that use flow relationships follow the same pattern as
the CPU analysis free functions in `src/cpu/`:

- Take `const lsdtt::lsdtt_flowinfo& fi` as the flow argument
- Take `const lsdtt::lsdtt_raster_cpu& dem` if elevation data is also needed
- Return `lsdtt::lsdtt_raster_cpu` or `lsdtt::lsdtt_array2d_f`
- Live in `src/cpu/` with a `_cpu` suffix on the filename

The function does not need to know about `lsdtt_flow_base` — that is
infrastructure. Include
`lsdtt_flowinfo_cpu.hpp` and work with `lsdtt_flowinfo` directly.
