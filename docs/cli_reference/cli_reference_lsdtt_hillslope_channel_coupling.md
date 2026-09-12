# CLI reference: lsdtt-hillslope-channel-coupling

Measures hillslopes against the channels they drain to. For every pixel on a
drainage divide it records the hilltop curvature, traces the flow path down to
the channel network, and reports the hillslope length and relief — together with
the steepness of the channel that path reached.

The output is one row per hilltop pixel carrying the non-dimensional pair
**E\*** and **R\*** of Roering et al. (2007), which is what lets a measured
landscape be compared against nonlinear hillslope transport theory.

This is a port of the refactored path of `lsdtt-hillslope-channel-coupling` in
LSDTopoTools2. The legacy alternatives (`use_legacy_HFR`,
`complete_hilltop_masking`, the five `HilltopFlowRouting*` variants) are **not**
ported.

For the streamline tracer underneath it see [Flow routing](../flow_routing.md);
for the channel network that supplies the sinks see
[lsdtt-channel-extraction](cli_reference_lsdtt_channel_extraction.md).

## Usage

```
lsdtt-hillslope-channel-coupling <path> <parameter_file>
lsdtt-hillslope-channel-coupling <parameter_file>   # path defaults to ./
lsdtt-hillslope-channel-coupling -h                 # print help and exit
lsdtt-hillslope-channel-coupling -v                 # print version and exit
```

## Parameter file format

```
# Lines beginning with # are comments.
key: value
```

---

## Before you start: check the vertical units

**Every quantity this tool reports scales linearly with elevation.** If a DEM's
elevations are in feet on a metre grid — the usual US state lidar delivery — then
curvature, relief, gradient, E\* and R\* are all too large by 3.28, and **nothing
looks broken**. The landscape merely appears steeper than it is. No raster format
records a vertical unit, so no reader can detect this.

It is not a hypothetical. On a real 1 m dataset the error presented not as a bug
but as a *physical finding*: three quarters of hilltops breached the hard bound
`R* < 1`, which reads as evidence that the critical gradient `S_c` is wrong for
the terrain. It is not — it is evidence that the elevations were in feet.

**Establish the vertical unit before interpreting any R\* population.** Run the
DEM through [lsdtt-raster-preprocessing](cli_reference_lsdtt_raster_preprocessing.md)
first: it warns when the median slope is physically implausible and, in US 3DEP
coverage, queries USGS to settle metres against feet. This tool accepts
`vertical_units` too, but the warning and the check live in preprocessing,
because that is where a corrected raster gets written.

---

## What the tool computes

### Divides and hilltop pixels

A drainage divide runs **between** DEM pixels, not along them. `lsdtt_ridgeline`
is an ordered sequence of pixel *corners* tracing that divide. The **hilltop
pixels** are the DEM pixels the ridgeline separates — the ones immediately either
side of it.

Each hilltop pixel therefore carries where it sits along its divide
(`ridgeline_id`, `ridgeline_node_index`, `distance_along_ridgeline_m`) and which
of the two opposing hillslopes it belongs to (`ridge_side`: 0 left of travel,
1 right).

`basin_selection_choice: sibling` is recommended: it makes an interfluve
structural, since a divide's two flanking basins share a parent by construction.

### Tracing

Each hilltop pixel is traced downslope by Lea (1992) aspect-driven kinematic
routing with Tarboton (1997) edge handling, at sub-pixel resolution, until it
reaches a channel cell. `trace_surface_choice` selects the raw DEM (default, so
paths follow true topography) or the fitted polynomial surface.

### The metrics

```
Lh   hillslope length   path length from hilltop to channel
R    relief             z_hilltop - z_channel
S    mean path gradient R / Lh
E*   = 2 |Cht| Lh / S_c
R*   = R / (Lh S_c)
```

`R* < 1` is a hard bound of the framework, not a fitted quantity.

!!! danger "Cht is ALWAYS the Laplacian — this is not configurable"
    `E* = 2|Cht|Lh/S_c` is defined with the Laplacian `∇²z` (Roering et al.
    2007; Hurst et al. 2012), so that is what the tool uses. There is no
    parameter to change it.

    A `curvature_metric_choice` parameter used to exist and defaulted to
    `mean_curvature` — lsdtt3's **differential-geometry** mean curvature,

    ```
    H = ((1+q²)r − 2pqs + (1+p²)t) / (2(1+p²+q²)^{3/2})
    ```

    which is a different quantity. It equals `∇²z/2` only where the gradient is
    zero; at finite gradient the metric terms make the ratio slope-dependent.
    Using it biased `E*` by about a factor of two, and biased gentle and steep
    hilltops by different amounts.

    **Nothing about the wrong answer looked wrong**: `R*` stayed below 1, the
    cloud stayed tight, no warning was printed — the points simply sat in the
    wrong place on the curve. On a synthetic landscape with known forcing
    (`U/D` = 0.0080 m⁻¹), `mean_curvature` recovered 0.0035 (0.43×) and the
    Laplacian 0.0071 (0.89×). A parameter file that still sets the old key gets
    a note saying it was ignored.

!!! warning "Curvature sign"
    lsdtt3 uses the strict mathematical convention, so a **convex-up hilltop has
    NEGATIVE curvature**. E\* uses `|Cht|`. A positive `Cht` is a concave pixel,
    which `remove_positive_curvature_hilltops` discards by default.

!!! note "What `Lh` measures"
    `Lh` is the path length to where the trace **enters** the channel cell, not
    to the cell centre — short by up to half a diagonal. `mean_path_gradient`
    inherits a small high bias from this (~2 % on a synthetic test).

---

## Trace outcomes

**Every hilltop pixel gets a row, including the ones that were never traced.**
`trace_status` says what happened, and **only status 0 has a meaningful `Lh`**:

| status | meaning |
|---:|---|
| **0** | **reached a channel** — the only rows with valid `Lh`, `R`, E\*, R\* |
| 1 | hit the edge of the DEM |
| 2 | **loop detected** — revisited a cell where it could still have descended: a genuine routing failure |
| 3 | no aspect at the source pixel |
| 4 | never offered to the tracer (sits on the DEM's outer ring) |
| 5 | **stalled on flat** — descended correctly, then ran out of descent with no channel to terminate on |

### Why 2 and 5 are different, and why it matters

They look like the same failure and are not.

A **stall (5)** is a statement about the **channel network's coverage**, not
about the routing. The path is correct: it descends monotonically from ridge to
valley floor and then oscillates on the flat because there is no channel cell
there to stop it. On a real 1 m dataset such traces ended a median of 116 m from
the nearest channel, while successful traces ended *on* one 100 % of the time.
**74 % of what was formerly reported as looping was this.**

A **loop (2)** is a genuine routing failure — the trace revisited a cell on
ground it could still have descended. On the same dataset this was **1.8 %** of
hilltops within one hillslope length of a channel.

The split is decided by `flat_terminus_gradient` (default 0.02): the steepest
descent available at the terminating cell.

**A high stall count means extend the channel network, or read the `Lh`
distribution as covering only the part of the DEM the network reaches.** A high
loop count means the `Lh` distribution is not reliable. The tool prints different
advice for each.

---

## Ingesting a channel network

`channel_network_fname` takes the output of `lsdtt-channel-extraction` directly
(FlatGeobuf, GeoPackage or CSV), so both tools share exactly one definition of
where the channels are. Channel heads are derived as channel cells with no
channel donor.

**This matters more than it looks: the channel head position is the single thing
that sets hillslope length.** A curvature-based detector (`wiener`) places it on
evidence; `threshold_contributing_pixels` places it on a chosen number. Prefer the
former, and ingest it rather than re-deriving it.

The network must come from **the same DEM** this tool reads. Points falling on
nodata are dropped with a count; if all of them do, the run stops and says so.

---

## Outputs

Line outputs follow `write_vector_format`; hilltop points may override it with
`hilltop_points_write_vector_format`, which also accepts `csv`. Rasters follow
`write_output_format`.

| file | contents |
|---|---|
| `*_hilltop_points` | **The main product.** One point per hilltop pixel, with curvature, `Lh`, relief, E\*, R\*, `trace_status`, the ridgeline ordering, and the joined channel properties. |
| `*_ridgelines` | The ordered divides. |
| `*_hillslope_traces` | The traced flow paths, thinned by `trace_thinning` on write. |
| `*_ridge_mask`, `*_Lh`, `*_channel_mask` | Optional rasters. |
| `*_README.csv` / `.html` | Every parameter used, with its help text. |

!!! warning "Undefined values are EMPTY"
    A blank field means undefined. It is never `0` and never `-9999`. A blank
    `E_star` means E\* could not be computed, not that it is zero. **Filter to
    `trace_status == 0` before computing any distribution.**

### Channel properties joined to each hilltop

The channel a hilltop drains to is known from the trace, so its properties are
attached in process: `channel_k_sn`, `channel_chi`, `channel_elevation`,
`channel_flow_distance_m`, `channel_drainage_area_m2`, `channel_basin_index`,
`channel_segment_id`. A post-hoc CSV join on floating-point coordinates would
drop rows silently.

Rows whose channel pixel was not part of a fitted network keep an empty
`channel_k_sn` — they are not zero and are not dropped.

---

## Worked example

Southern Ohio, 1 m lidar, 5.33 km². Full write-up including the corrected DEM
and the diagnostics in `../data/IndianCreek_1m/hillslope_channel_coupling_readme.md`.

```bash
# 1. correct the DEM first (this one is in FEET, with nodata stored as zero)
lsdtt-raster-preprocessing ./ ft_to_m.param

# 2. channels, from curvature rather than an area threshold
lsdtt-channel-extraction ./ IC_wiener_final.param

# 3. couple
lsdtt-hillslope-channel-coupling ./ hcc_final.param
```

The settings that matter (the full file adds output-format switches and the
optional rasters):

```
# hcc_final.param
read_path: ./
read_fname: indian_creek_metres.tif
write_path: ./
write_fname: IndianCreek_hcc

boundary_conditions: oooo
fill_before_routing: true

channel_network_fname: ./IndianCreek_channel_network_points.fgb

basin_selection_choice: sibling
min_basin_pixels: 500

ridge_extraction_algorithm_choice: ridgeline
stream_network_padding_pixels: 3
threshold_hilltop_gradient: 0.8
remove_positive_curvature_hilltops: true

surface_fitting_window_radius_metres: 6.0
critical_hillslope_slope: 0.8

run_chi_analysis: true
m_over_n: 0.45

hilltop_points_write_vector_format: flatgeobuf
write_ridgelines: true
write_hillslope_traces: true
trace_thinning: 50
write_output_format: gtiff
n_threads: 4
```

Result: 62,365 hilltop pixels, 55 % reaching a channel, median `Lh` 60 m, median
R\* 0.427, and 0.2 % breaching `R* < 1` at `S_c = 0.8`.

### Figures

`lsdviztools3` renders the output set:

```bash
pixi run python -m lsdviztools3.cli.plot_hilltop_coupling \
    --prefix   ./IndianCreek_hcc \
    --dem      ./indian_creek_metres.tif \
    --channels ./IndianCreek_hcc_channel_mask.tif \
    --outdir   ./hillslope_coupling_figures
```

`traces_on_hillshade`, `trace_outcomes`, `looped_trace_profiles`,
`ridgelines_and_sides`, `length_distributions`, `ridgeline_profiles`,
`estar_rstar`.

---

## Parameters

Generated from the tool's own help registry; `lsdtt-hillslope-channel-coupling -h`
prints the same information.

### Input

| parameter | type | default | what it does |
|---|---|---|---|
| `boundary_conditions` | string | `nnnn` | Four-character N,E,S,W boundary code for the flow router. |
| `fill_before_routing` | bool | `true` | Pit-fill the DEM before building the flow network. |
| `min_elevation` | float | `-9999` | Pixels at or below this elevation become nodata. The disabled sentinel preserves below-sea-level study areas; set 0 explicitly to remove sea or nodata-as-zero pixels. |
| `vertical_units` | string | `metres` | Vertical unit of the input elevations: metres, feet or us_survey_feet. |

### Basin selection

| parameter | type | default | what it does |
|---|---|---|---|
| `basin_selection_choice` | string | `sibling` | How basins are chosen; ridgelines come from their outlines. |
| `min_basin_pixels` | int | `100` | Skip basins smaller than this. |

### Ridge extraction

| parameter | type | default | what it does |
|---|---|---|---|
| `diagonal_step_policy` | string | `skip` | What to do with a non-4-connected step in a ridgeline. |
| `divide_merge_search_steps` | int | `5` | Merge tolerance K for the divide_pairs fallback route. |
| `remove_positive_curvature_hilltops` | bool | `true` | Drop hilltop pixels whose curvature is positive. |
| `ridge_extraction_algorithm_choice` | string | `ridgeline` | How hilltop pixels are found. |
| `stream_network_padding_pixels` | int | `2` | Drop hilltop pixels within this many pixels of a channel. |
| `threshold_hilltop_gradient` | float | `0.4` | Drop hilltop pixels steeper than this. <= 0 disables. |

### Channel network

| parameter | type | default | what it does |
|---|---|---|---|
| `channel_network_fname` | string | `` | Optional channel network to ingest instead of extracting one. |
| `threshold_contributing_pixels` | int | `1000` | Contributing pixels defining a channel head. |

### Floodplain

| parameter | type | default | what it does |
|---|---|---|---|
| `find_floodplain` | bool | `false` | Detect a floodplain with a simple gradient + HAND rule. |
| `floodplain_gradient_threshold` | float | `0.1` | Gradient below which a pixel may be floodplain. |
| `floodplain_hand_threshold` | float | `5.0` | HAND (m) below which a pixel may be floodplain. |
| `floodplain_mask_fname` | string | `` | Optional binary floodplain raster to ingest. |

### Surface metrics

| parameter | type | default | what it does |
|---|---|---|---|
| `surface_fitting_window_radius_metres` | float | `0` | Polynomial surface-fitting window radius in metres. 0 = 3 * cell size. |

### Tracing

| parameter | type | default | what it does |
|---|---|---|---|
| `flat_terminus_gradient` | float | `0.02` | Gradient below which a trace that revisits a cell is reported as stalled_on_flat (5) rather than loop_detected (2). |
| `hilltop_thinning` | int | `1` | Analyse every Nth hilltop pixel. 1 = all of them. |
| `plan_curvature_threshold` | float | `0.001` | abs(plan curvature) above this counts a path cell as curved. |
| `trace_max_steps` | int | `100000` | Hard cap on cell steps per trace. |
| `trace_surface_choice` | string | `raw` | Which surface the flow paths are traced on: raw or fitted. |
| `trace_thinning` | int | `20` | Write every Nth trace geometry. Affects the FILE only. |

### Chi

| parameter | type | default | what it does |
|---|---|---|---|
| `chi_A0` | float | `1.0` | Reference drainage area A0 for the chi integration (m²). **Leave at 1.0** unless you know why: a chi-elevation gradient equals `k_sn` only when `A0 = 1 m²`; any other value silently rescales it by `A0^(m/n)`. *(Renamed from `chi_A_0`; the old key still works for one release, with a warning.)* |
| `chi_critical_divergence_metres` | float | `5.0` | Segment-fitting tolerance: max divergence (m) before a new segment. |
| `m_over_n` | float | `0.45` | Concavity index m/n used for the chi integration. |
| `run_chi_analysis` | bool | `false` | Fit the chi network and attach channel steepness to every hilltop row. |

### Non-dimensional

| parameter | type | default | what it does |
|---|---|---|---|
| `critical_hillslope_slope` | float | `0.8` | Critical gradient S_c used to non-dimensionalise E* and R*. |

### Output

| parameter | type | default | what it does |
|---|---|---|---|
| `write_vector_format` | string | `flatgeobuf` | Format for traced flow paths and ridgeline polylines. |
| `hilltop_points_write_vector_format` | string | `same_as_write_vector_format` | Optional hilltop point-output override; also accepts `csv`. |
| `write_channel_mask_raster` | bool | `false` | Write the channel mask actually used as sinks. |
| `write_floodplain_mask_raster` | bool | `false` | Write the floodplain mask actually used. |
| `write_hand_raster` | bool | `false` | Write the HAND raster used by the floodplain detector. |
| `write_hillslope_length_raster` | bool | `false` | Write Lh at the hilltop pixels. |
| `write_hillslope_traces` | bool | `false` | Write the traced flow paths as lines, thinned by trace_thinning. |
| `write_hilltop_curvature_raster` | bool | `false` | Write Cht at the hilltop pixels. |
| `write_hilltop_points` | bool | `true` | Write the hilltop table (one row per hilltop pixel). |
| `write_output_format` | string | `cog` | Raster output format: cog, gtiff or png. |
| `write_relief_raster` | bool | `false` | Write hillslope relief at the hilltop pixels. |
| `write_ridge_mask_raster` | bool | `false` | Optional coupling diagnostic: hilltop pixels coded by side. Unlike feature-finder's primary mask output, it is off by default. |
| `write_ridgelines` | bool | `false` | Write the ordered ridgelines as lines. |

### Performance

| parameter | type | default | what it does |
|---|---|---|---|
| `n_threads` | int | `0` | Thread count. 0 = hardware default. |

---

## References

Lea, N.L. (1992), An aspect-driven kinematic routing algorithm, in *Overland
Flow: Hydraulics and Erosion Mechanics*.

Roering, J.J., Perron, J.T. and Kirchner, J.W. (2007), Functional relationships
between denudation rate and hilltop curvature. *EPSL* 264:245-258.

Tarboton, D.G. (1997), A new method for the determination of flow directions and
upslope areas in grid digital elevation models. *Water Resources Research*
33:309-319.

Hurst, M.D., Mudd, S.M., Walcott, R., Attal, M. and Yoo, K. (2012), Using
hilltop curvature to derive the spatial distribution of erosion rates. *JGR Earth
Surface* 117:F02017.
