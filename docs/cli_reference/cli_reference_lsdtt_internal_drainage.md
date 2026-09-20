# CLI reference: lsdtt-internal-drainage

Finds closed depressions that never drain to the edge of a DEM — moulins,
playas, sabkhas, karst swallow holes, quarries, tailings ponds — and routes flow
into them instead of over an invented spill point.

Ordinary flow routing destroys these features: filling or breaching raises or
carves the depression until it overtops its lowest saddle, and the water is then
routed away through a saddle that does not exist in the landscape. This tool
treats the depression as the terminus it really is.

For a complete worked example on a real dataset, including the human review step,
The examples below show the complete review and routing workflow.
For the sink detectors and the graph machinery behind them, see
[Flow routing](../design_and_algorithms/flow_routing.md).

## Usage

```
lsdtt-internal-drainage <path> <parameter_file>
lsdtt-internal-drainage <parameter_file>       # path defaults to ./
lsdtt-internal-drainage -h                     # print help and exit
lsdtt-internal-drainage -v                     # print version and exit
```

## Parameter file format

```
# Lines beginning with # are comments.
key: value
```

---

## The three stages

`internal_drainage_stage` selects what a run does. The intended sequence is
**detect → (human review) → route**, because deciding which candidate
depressions are real sinks is a judgement a person makes over imagery, not
something the software can decide.

| stage | what it does | needs |
|---|---|---|
| `detect` | Finds candidate depressions and writes them as points, a CSV and polygons for review. | a DEM |
| `route` | Consumes a reviewed candidate file, routes flow into each accepted sink, tags basins and optionally extracts channels. | a DEM + `sink_candidate_fname` |
| `both` | Runs detection and routing in one pass, accepting every candidate without review. | a DEM |

Use `both` only when you already trust the detector on your terrain — for
example on synthetic surfaces, or a landscape where every closed depression is
genuinely a sink. On real DEMs the review step is the point.

---

## Choosing a detector

`sink_detection_algorithm_choice` has two settings and the difference matters
more than most parameters.

**`fill_difference`** (the default) takes each connected patch that filling would
raise and treats it as one sink, keeping its deepest cell as the pit. It is fast,
robust on rough surfaces, and has no concept of nesting: a lake containing
several hollows is a single sink.

**`depression_hierarchy`** builds the full nesting tree, so a lake containing
hollows becomes a parent with children. It additionally reports
`parent_sink_id`, `hierarchy_level` and a true accumulated depression volume. Use
it where that structure is the point — playas, karst, quarry benches.

**The risk, stated plainly:** on a rough surface every innermost hollow is tiny.
Combining `depression_hierarchy` with `sink_leaves_only: true` and a
`min_pixels_in_patch` threshold can filter out *every* sink and report
nothing at all, with no error. On an 18 m Greenland ice-surface DEM at
`min_pixels_in_patch: 1000`, `fill_difference` finds 143 sinks and
`depression_hierarchy` finds **zero** — not because it failed, but because all
468,878 leaf hollows it found are smaller than the threshold.

If a run reports no sinks, check this first. Remedies: smooth the DEM, filter on
`min_sink_depth` or `min_sink_volume` rather than pixel count, or set
`sink_leaves_only: false`.

### What `sink_leaves_only` means

A large lake often contains several smaller closed hollows inside it. The
innermost hollow — the one with nothing nested inside it — is a **leaf**.

```
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        rim of the lake
       /                                         \
      /    ___________         ___________        \
     /    /  hollow A \       /  hollow B \         \      each hollow is a LEAF
    /    /   .-----.   \     /   .-----.   \         \
   /    /   /   C   \   \   /    \       /  \         \    C nests inside A,
  /    /    \ (leaf)/    \ /      `-----`    \         \   so A is NOT a leaf
 /     \     `-----`     /  \                 /         \
/       `---------------`    `---------------`           \
```

- `sink_leaves_only: true` (default) returns **C and B** — the innermost hollows,
  one sink each.
- `sink_leaves_only: false` also returns **A** and the lake itself, so you get the
  enclosing depressions as well. Those may share pit cells with their children
  and must be reviewed before routing, which is why
  `internal_drainage_stage: both` requires `sink_leaves_only: true`.

This applies **only** to `depression_hierarchy`. Under `fill_difference` there is
no tree to take the leaves of; the whole lake is one sink, and setting
`sink_leaves_only: false` prints a warning and is otherwise ignored.

---

## Reviewing candidates in QGIS

This is the round trip the vector output exists for, and it needs no knowledge of
`lsdtt3` internals.

1. Run `detect`. It writes `<write_fname>_sink_candidates.fgb` (points),
   `..._sink_candidates.csv` (the same records as a table) and
   `..._sink_depressions.fgb` (the depression footprints as polygons).
2. Open the `.fgb` in QGIS over satellite or aerial imagery. Add the depression
   polygons too — the footprint is usually what tells you whether a candidate is
   a real closed basin or DEM noise.
3. **Delete the features that are not real sinks.** You may also **add** features
   by hand where you can see a sink the DEM did not resolve: place a point and
   leave the other attributes blank.
4. Save the layer.
5. Pass the saved file as `sink_candidate_fname` to a `route` run.

Keep the saved layer. It is the one part of the workflow that cannot be
regenerated, and it is data rather than output.

### Candidate file formats

`sink_candidate_fname` accepts FlatGeobuf, GeoPackage and CSV.

- **Vector formats (`.fgb`, `.gpkg`) are simplest** — position comes from the
  geometry, so no column naming is involved. FlatGeobuf is the default.
- **CSV** must carry either `latitude` and `longitude` (read as EPSG:4326) or
  `easting` and `northing` — `x` and `y` also work — read as being **in the
  raster's own coordinate system**. If a file has both, latitude and longitude
  win.

Features are matched on their durable `sink_id`, never on file position, so a
FlatGeobuf spatial index that reorders features on write is harmless.

> **GeoPackage on WSL2**: writing `.gpkg` to a mounted Windows drive
> (`/mnt/c/...`) can fail or corrupt, because SQLite needs file locking that the
> mount does not provide reliably. Use FlatGeobuf, or write to the Linux
> filesystem.

---

## Routing modes

`sink_mode` decides how an accepted sink is made to swallow water.

| mode | what it does | when to use it |
|---|---|---|
| `forced_outlet` (default) | Leaves the DEM untouched and forces the routing so the sink cells act as outlets. | Whenever you want elevations preserved — the non-destructive choice. |
| `nodata_hole` | Punches the sink cells out of the DEM as nodata, so water leaves there. | Reproducing analyses that worked this way, or when you want the hole in the output DEM. |

Both produce the same basin tagging on their common domain.

`breach_before_fill` is currently rejected: breaching with sinks is not yet
implemented, and asking for it raises an error rather than silently filling.

---

## Outputs

| file | written when |
|---|---|
| `_sink_candidates.fgb` / `.csv` | `detect`, `both` |
| `_sink_depressions.fgb` / `.csv` | `detect`, `both` |
| `_fill_diff.tif` | `fill_difference` detection |
| `_accepted_sinks.fgb` / `.csv` | `route`, `both` |
| `_accepted_depressions.fgb` | `route`, `both` |
| `_sink_id.tif` | the basin draining to each sink |
| `_sink_drainage_area.tif` | contributing area per cell |
| `_sink_stats.csv` | per-sink area, depth, volume, elevation and slope percentiles, channel length by stream order |
| `_sink_basins.fgb` | basin polygons, with `write_sink_basin_polygons` |
| `_sink_channels.fgb`, `_sink_channel_nodes.csv`, `_sink_junctions.fgb` | with `extract_sink_channel_networks` |
| `_sink_chi.tif` | with `calculate_sink_chi` |
| `_punched.tif` | `nodata_hole` mode, with `write_sink_punched_raster` |

A sink whose footprint is split into several pieces — which happens when
`min_pixels_in_patch` merges a small hollow into a larger ancestor — is
written as one polygon record **per piece**, all carrying the same sink
attributes. No terrain is lost.

---

## Parameters

Generated from the tool's own help registry; `lsdtt-internal-drainage -h`
prints the same information.

### Run control

| parameter | type | default | what it does |
|---|---|---|---|
| `n_threads` | int | `0` | Thread count. 0 uses available threads. |

### Workflow

| parameter | type | default | what it does |
|---|---|---|---|
| `internal_drainage_stage` | string | `detect` | Workflow stage. detect writes candidates; route consumes reviewed candidates; both runs without review. |

### Detection

| parameter | type | default | what it does |
|---|---|---|---|
| `min_pixels_in_patch` | int | `1` | Minimum sink pixels. Legacy-compatible minimum depression size. |
| `min_sink_area` | float | `0.0` | Minimum depression area. Square metres. |
| `min_sink_depth` | float | `0.0` | Minimum sink depth. Reject shallower candidates. |
| `min_sink_volume` | float | `0.0` | Minimum depression volume. Cubic metres. |
| `sink_detection_algorithm_choice` | string | `fill_difference` | Sink detector. fill_difference (default) takes each connected patch that filling would raise as one sink: robust, fast, and what the published Greenland work used. depression_hierarchy additionally resolves NESTING, so a lake containing several hollows becomes a tree rather than one patch, and it reports parent_sink_id, hierarchy_level and true depression volume. Use it for genuinely nested terrain such as playas, karst or quarry benches. RISK: on rough surfaces every innermost hollow is tiny, so combining it with sink_leaves_only and a min_pixels_in_patch threshold can filter out EVERY sink and report nothing at all. |
| `sink_leaves_only` | bool | `true` | Emit one sink per innermost closed hollow. A large lake can contain several smaller closed hollows. true returns each innermost hollow, meaning one with nothing nested inside it; false also returns enclosing depressions, which may share pit cells and must be reviewed before routing. |

### Review

| parameter | type | default | what it does |
|---|---|---|---|
| `sink_candidate_fname` | string | `(none)` | Reviewed candidate file. Required for route; CSV, FlatGeobuf, or GeoPackage. |
| `sink_candidate_input_format` | string | `same_as_extension` | Candidate format. same_as_extension, csv, flatgeobuf, or gpkg. |
| `sink_candidate_layer_name` | string | `(none)` | Candidate layer name. Optional layer name for multi-layer GeoPackages. |

### Routing

| parameter | type | default | what it does |
|---|---|---|---|
| `boundary_conditions` | string | `oooo` | Flow boundaries. Phase 4 supports open boundaries only for sink-seeded filling. |
| `breach_before_fill` | bool | `false` | Breach before filling (currently rejected). |
| `fill_min_slope` | float | `0.0001` | Fill slope. Minimum gradient imposed while preserving retained sinks. |
| `sink_mode` | string | `forced_outlet` | Sink routing mode. forced_outlet is non-destructive; nodata_hole reproduces legacy routing. |
| `upslope_tagging_distance` | float | `1e12` | Tagging distance. Maximum upstream flow distance from each sink in metres. |

### Channels

| parameter | type | default | what it does |
|---|---|---|---|
| `extract_sink_channel_networks` | bool | `false` | Extract channels in sink basins. |
| `threshold_contributing_pixels` | int | `1000` | Contributing-pixel threshold defining channel sources. |
| `write_sink_channel_lines` | bool | `true` | Write sink channel lines. |
| `write_sink_channel_nodes_csv` | bool | `true` | Write sink channel nodes. |
| `write_sink_junctions` | bool | `false` | Write sink junction points. |

### Chi

| parameter | type | default | what it does |
|---|---|---|---|
| `chi_A0` | float | `1.0` | Reference drainage area A0 for the chi integration (m²), positive. **Leave at 1.0** unless you know why: a chi-elevation gradient equals `k_sn` only when `A0 = 1 m²`; any other value silently rescales it by `A0^(m/n)`. *(Renamed from `A_0`; the old key still works for one release, with a warning.)* |
| `calculate_sink_chi` | bool | `false` | Calculate full-raster sink-aware chi. |
| `m_over_n` | float | `0.45` | Chi theta. Chi is descriptive for multi-elevation moulins, not an erosion-rate proxy. |
| `write_chi_raster` | bool | `true` | Write sink-aware chi raster. |

### Output

| parameter | type | default | what it does |
|---|---|---|---|
| `write_vector_format` | string | `flatgeobuf` | Vector output format. flatgeobuf is recommended; gpkg is supported with a filesystem warning. |
| `write_accepted_depression_polygons` | bool | `true` | Write accepted depression polygons. |
| `write_accepted_sink_csv` | bool | `true` | Write accepted CSV. |
| `write_accepted_sink_points` | bool | `true` | Write accepted point vector. |
| `write_flow_direction_raster` | bool | `false` | Write ArcMap D8 directions. |
| `write_output_format` | string | `cog` | Raster output format. cog, gtiff, or same_as_input. PNG is rejected because it rescales analytical labels and values. |
| `write_sink_basin_polygons` | bool | `false` | Write tagged basin polygons. |
| `write_sink_candidate_csv` | bool | `true` | Write candidate CSV. |
| `write_sink_candidate_points` | bool | `true` | Write candidate point vector. |
| `write_sink_depression_csv` | bool | `false` | Write depression attribute CSV. |
| `write_sink_depression_polygons` | bool | `true` | Write candidate depression polygons. |
| `write_sink_drainage_area_raster` | bool | `true` | Write per-cell contributing area. |
| `write_sink_id_raster` | bool | `true` | Write routing sink labels. |
| `write_sink_punched_raster` | bool | `false` | Write punched DEM in hole mode. |
| `write_sink_stats_csv` | bool | `true` | Write joined sink statistics. |
| `write_vector_spatial_index` | bool | `false` | Build vector spatial indices. |

<!-- BEGIN GENERATED CLI PARAMETERS -->
## Parameters

This section is generated from the executable's `-h json` registry. Do not edit it manually.

### File I/O

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `read_path` | `string` | `./` | `io` | - | Directory from which input files are read. | Defaults to the directory holding the parameter file. |
| `read_fname` | `string` | - | `io` | - | Input filename, with or without extension (e.g. my_dem.tif). | Read from read_path. lsdtt3 uses GDAL, so any GDAL-readable raster format works. If empty, defaults to the parameter filename stem. |
| `write_path` | `string` | - | `io` | - | Directory to which output files are written. | Defaults to read_path when not set. |
| `write_fname` | `string` | - | `io` | - | Output filename prefix (no path, no extension). | All outputs are named &lt;write_fname&gt;_&lt;descriptor&gt;.&lt;ext&gt;. Defaults to the stem of read_fname when not set. |

### Channel network

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `threshold_contributing_pixels` | `int` | `1000` | `param` | - | Contributing-pixel threshold defining a channel head. | The total number of upslope contributing pixels from the flow-routing object, not the number of direct donors to a pixel. Reduce this value to extract a denser network. This is a PIXEL COUNT and therefore depends on the grid resolution: the equivalent drainage area is threshold_contributing_pixels x cell area. Where a genuine area is wanted instead, the convention is &lt;what&gt;_area_threshold_m2. |

### Workflow

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `internal_drainage_stage` | `string` | `detect` | `switch` | - | Workflow stage. | detect writes candidates; route consumes reviewed candidates; both runs without review. |

### Detection

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `sink_detection_algorithm_choice` | `string` | `fill_difference` | `switch` | - | Sink detector. | fill_difference (default) takes each connected patch that filling would raise as one sink: robust, fast, and what the published Greenland work used. depression_hierarchy additionally resolves NESTING, so a lake containing several hollows becomes a tree rather than one patch, and it reports parent_sink_id, hierarchy_level and true depression volume. Use it for genuinely nested terrain such as playas, karst or quarry benches. RISK: on rough surfaces every innermost hollow is tiny, so combining it with sink_leaves_only and a min_pixels_in_patch threshold can filter out EVERY sink and report nothing at all. |
| `min_pixels_in_patch` | `int` | `1` | `param` | - | Minimum sink pixels. | Legacy-compatible minimum depression size. |
| `min_sink_depth` | `float` | `0.0` | `param` | - | Minimum sink depth. | Reject shallower candidates. |
| `min_sink_area` | `float` | `0.0` | `param` | - | Minimum depression area. | Square metres. |
| `min_sink_volume` | `float` | `0.0` | `param` | - | Minimum depression volume. | Cubic metres. |
| `sink_leaves_only` | `bool` | `true` | `switch` | - | Emit one sink per innermost closed hollow. | A large lake can contain several smaller closed hollows. true returns each innermost hollow, meaning one with nothing nested inside it; false also returns enclosing depressions, which may share pit cells and must be reviewed before routing. |

### Routing

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `sink_mode` | `string` | `forced_outlet` | `switch` | - | Sink routing mode. | forced_outlet is non-destructive; nodata_hole reproduces legacy routing. |
| `boundary_conditions` | `string` | `oooo` | `param` | - | Flow boundaries. | Phase 4 supports open boundaries only for sink-seeded filling. |
| `upslope_tagging_distance` | `float` | `1e12` | `param` | - | Tagging distance. | Maximum upstream flow distance from each sink in metres. |
| `fill_min_slope` | `float` | `0.0001` | `param` | - | Fill slope. | Minimum gradient imposed while preserving retained sinks. |
| `priority_flood_queue` | `string` | `binary` | `param` | - | Priority-Flood queue. | binary (default) or radix. Used for retained-sink and punched-sink routing. |
| `breach_before_fill` | `bool` | `false` | `switch` | - | Breach before filling (currently rejected). | Breach before filling (currently rejected). |

### Review

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `sink_candidate_fname` | `string` | - | `input` | `internal_drainage_stage` | Reviewed candidate file. | Required for route; CSV, FlatGeobuf, or GeoPackage. |
| `sink_candidate_input_format` | `string` | `same_as_extension` | `param` | - | Candidate format. | same_as_extension, csv, flatgeobuf, or gpkg. |
| `sink_candidate_layer_name` | `string` | - | `param` | - | Candidate layer name. | Optional layer name for multi-layer GeoPackages. |

### Output

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_output_format` | `string` | `cog` | `param` | - | Raster output format. | cog, gtiff, or same_as_input. PNG is rejected because it rescales analytical labels and values. |
| `write_vector_format` | `string` | `flatgeobuf` | `param` | - | Vector output format. | flatgeobuf is recommended; gpkg is supported with a filesystem warning. |
| `write_sink_candidate_points` | `bool` | `true` | `switch` | - | Write candidate point vector. | Write candidate point vector. |
| `write_sink_candidate_csv` | `bool` | `true` | `switch` | - | Write candidate CSV. | Write candidate CSV. |
| `write_sink_depression_polygons` | `bool` | `true` | `switch` | - | Write candidate depression polygons. | Write candidate depression polygons. |
| `write_sink_depression_csv` | `bool` | `false` | `switch` | - | Write depression attribute CSV. | Write depression attribute CSV. |
| `write_accepted_sink_points` | `bool` | `true` | `switch` | - | Write accepted point vector. | Write accepted point vector. |
| `write_accepted_sink_csv` | `bool` | `true` | `switch` | - | Write accepted CSV. | Write accepted CSV. |
| `write_accepted_depression_polygons` | `bool` | `true` | `switch` | - | Write accepted depression polygons. | Write accepted depression polygons. |
| `write_sink_basin_polygons` | `bool` | `false` | `switch` | - | Write tagged basin polygons. | Write tagged basin polygons. |
| `write_sink_punched_raster` | `bool` | `false` | `switch` | - | Write punched DEM in hole mode. | Write punched DEM in hole mode. |
| `write_sink_id_raster` | `bool` | `true` | `switch` | - | Write routing sink labels. | Write routing sink labels. |
| `write_sink_drainage_area_raster` | `bool` | `true` | `switch` | - | Write per-cell contributing area. | Write per-cell contributing area. |
| `write_flow_direction_raster` | `bool` | `false` | `switch` | - | Write ArcMap D8 directions. | Write ArcMap D8 directions. |
| `write_sink_stats_csv` | `bool` | `true` | `switch` | - | Write joined sink statistics. | Write joined sink statistics. |
| `write_vector_spatial_index` | `bool` | `false` | `switch` | - | Build vector spatial indices. | Build vector spatial indices. |

### Run control

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `n_threads` | `int` | `0` | `io` | - | Thread count. | 0 uses available threads. |

### Chi

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `m_over_n` | `float` | `0.45` | `param` | `calculate_sink_chi` | Chi theta. | Chi is descriptive for multi-elevation moulins, not an erosion-rate proxy. |
| `chi_A0` | `float` | `1.0` | `param` | `calculate_sink_chi` | Reference drainage area A0 for the chi integration (m2). | Positive reference area in square metres. LEAVE THIS AT 1.0 UNLESS YOU KNOW WHY YOU ARE CHANGING IT: a chi-elevation gradient equals the channel steepness index k_sn only when A0 = 1 m2, and any other value silently rescales it by A0^(m/n). |
| `calculate_sink_chi` | `bool` | `false` | `switch` | - | Calculate full-raster sink-aware chi. | Calculate full-raster sink-aware chi. |
| `write_chi_raster` | `bool` | `true` | `switch` | - | Write sink-aware chi raster. | Write sink-aware chi raster. |

### Channels

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `extract_sink_channel_networks` | `bool` | `false` | `switch` | - | Extract channels in sink basins. | Extract channels in sink basins. |
| `write_sink_channel_lines` | `bool` | `true` | `switch` | - | Write sink channel lines. | Write sink channel lines. |
| `write_sink_channel_nodes_csv` | `bool` | `true` | `switch` | - | Write sink channel nodes. | Write sink channel nodes. |
| `write_sink_junctions` | `bool` | `false` | `switch` | - | Write sink junction points. | Write sink junction points. |

<!-- END GENERATED CLI PARAMETERS -->
