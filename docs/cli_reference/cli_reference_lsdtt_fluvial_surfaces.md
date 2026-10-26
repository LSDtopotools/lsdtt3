# lsdtt-fluvial-surfaces

Sorts flat ground into fluvial levels — floodplain, terrace treads, fan lobes —
by banding a height reference and labelling connected components within each
band.

---

## Where this sits

Three tools, each doing one job:

| tool | question | needs |
|---|---|---|
| `lsdtt-flat-surface-finder` | which ground is flat? | a DEM |
| **`lsdtt-fluvial-surfaces`** | **which fluvial level is each flat surface on?** | a channel, or a supplied height field |
| `lsdtt-valley-centreline` | trace a centreline through this valley | a valley mask |

The split is at the point where a channel network becomes necessary.
`lsdtt-flat-surface-finder` is pure morphometry: no flow routing, no channel, no
height reference. This tool does everything that needs one, and writes the
valley mask `lsdtt-valley-centreline` consumes.

The normal two-step workflow:

```
lsdtt-flat-surface-finder  →  <stem>_combined_flat_mask.tif
                                        ↓  flat_mask_fname
lsdtt-fluvial-surfaces     →  <stem>_valley_mask.tif, _band_index.tif, …
                                        ↓  valley_mask_fname
lsdtt-valley-centreline    →  <stem>_valley_centreline.fgb
```

You can run this tool on its own — leave `flat_mask_fname` empty and it falls
back to a built-in relief or slope admission. That is for a first look, not for
a result you intend to publish; the texture metrics in
`lsdtt-flat-surface-finder` are much better at deciding what is flat.

---

## The idea

### Band first, then label

The obvious way to segment flat surfaces is to label connected components and
then work out how high each one is. That forces the riser network between
terraces to be **continuous**: a road, a gully or a stock track cutting through
a riser lets two treads at different heights merge into a single component, and
no amount of repair afterwards can separate them again.

This tool inverts the order. It bands the height reference first, then labels
components **within each band**. Two pixels can only share a label if they share
a band, and two treads at different heights are by construction in different
bands — so the merge cannot happen. The riser never has to be found, and never
has to be continuous.

That is the whole reason the tool exists. Everything else is in service of it.

### There is no valley mode and no terrace mode

Finding a valley floor and segmenting a terrace flight are the same computation
at different band counts. The valley floor is the lowest band of a two-band cut;
a terrace flight is an n-band ladder on the same population. So there is one
parameter, `ladder_level_count`:

- `2` — the valley cut. Band 0 is the valley floor.
- `n` — an n-level ladder.
- `0` (default) — automatic: the largest level count the DEM can actually
  resolve, see below.

### Which reference field

**`hand`** (default) — height above the in-scope channel network. Terraces along
a reach sit at a roughly constant height above the channel, so HAND makes a
terrace flight a step function. This is the field for fluvial terraces and
floodplains.

**`fitted_surface`** — the signed residual to a surface fitted to the flat
ground. A fan is *not* a step function in HAND: it climbs continuously from toe
to apex, so a HAND ladder cuts one surface into arcs banded by distance from the
apex rather than by age. What *is* roughly constant on a fan is the depth of an
inset lobe below the projection of the older surface — which is what the
residual measures.

Choose the fitted surface with `fitted_surface_algorithm_choice`:

| value | when |
|---|---|
| `cone` (default) | **a fan or a debris cone.** z falls linearly with distance from an apex. |
| `plane_ransac` | a piedmont, a terrace tread, a broadly planar surface. RANSAC because incised channels and inset lobes are large-amplitude outliers that would tilt an OLS fit. |
| `poly2` | a broadly curved surface with no single centre. |

**Use `cone` for a fan.** A fan is radial, and a cone is not a low-order
polynomial in x and y. Measured on a synthetic fan (2 km radius, 0.03 gradient,
lobes inset 4 m), the spread of the older surface's own residual — which should
be flat, so its spread is pure misfit — was about 34 m under a RANSAC plane and
10 m under a poly2, against a 4 m signal. Under a cone it was 0.3 m.

The cone's apex is found automatically as the centroid of the highest 0.1 % of
fit-mask cells, or you can state it with `fitted_surface_apex_x` /
`_apex_y`. It is supplied rather than solved for because it is a strongly
non-linear parameter and you can usually point at a fan head more reliably than
a solver can find it.

### How many levels, and when to believe them

A level count is **admissible** when the tightest gap between consecutive band
centres clears `ladder_resolvability_kappa * vertical_accuracy_metres`. The
answer is the largest admissible count.

`vertical_accuracy_metres` is metadata you must supply — it cannot be derived
from the DEM. Leaving it at 0 disables the check, and then **every window
reports `resolved` and the ladder carries no instrument guarantee at all**. The
tool warns when you do. (Note that a DEM whose vertical units are feet on a
metre grid has an effective accuracy 3.28 times larger than it looks.)

Where no count qualifies, the window reports `unresolved` rather than returning
a number. **Read `<stem>_ladder_windows.csv` before describing your levels**:
whether they may be called resolved levels at all depends on which branch the
ladder took, and that is not visible in any raster.

An explicit `ladder_level_count` is honoured and reported honestly — if that
partition is not admissible the window still records the count you asked for,
with status `unresolved`, rather than silently substituting a different one.

---

## Composability: `hand_fname`

HAND is the expensive stage, and the level count is what you iterate on. So:

```
# once
write_hand_raster: true

# then, as often as you like
hand_fname: <stem>_hand.tif
ladder_level_count: 5
```

The supplied-HAND path produces **bit-identical** bands and labels to the
computed path, which is pinned by a test.

One limit worth knowing: reach-relative windowing needs each pixel's HAND
*donor node*, and a raster cannot carry those. A supplied HAND is therefore
always banded globally. Iterate on the level count with a supplied HAND, then do
the final run computing HAND if you want the windowed ladder.

---

## Outputs

| file | what |
|---|---|
| `_band_index.tif` | Int16, nodata −1. Which level each pixel is on. Authoritative. |
| `_banded_labels.tif` | UInt32. Patch label, a global running counter 1..N; 0 is background. |
| `_valley_mask.tif` | Band 0. The input to `lsdtt-valley-centreline`. |
| `_surface_class.tif` | floodplain / terrace class. `hand` field only. |
| `_surface_patches.csv` | One row per patch: label, band, class, source_id, n_pixels, area, median field value, aspect. **The join key for everything downstream.** |
| `_ladder_windows.csv` | The resolvability report. Read it first. |
| `_reference_field.tif` | The field that was actually banded, whichever it was. |
| `_hand.tif` | The HAND raster, for feeding back as `hand_fname`. |

The band is **not** packed into the label. The label is a plain running counter
and the band travels in `_band_index.tif` and in the patch table, so nothing has
to divide a label to learn its level.

`_surface_class.tif` and `_surface_class.csv` are written only for the `hand`
field. Floodplain and terrace are statements about height above a channel; they
mean nothing on a fitted-surface residual, so no class is invented for one.

---

## Things that will bite you

**`antimode` cannot band a fitted surface.** The antimode break finder builds its
histogram on `log10` of the field, which suits HAND (non-negative by
construction) and destroys a signed one: every negative residual is floored into
the lowest bin, so the entire below-surface population — the inset lobes, the
thing you are looking for — collapses into a single spike. Selecting
`fitted_surface` switches the break finder to `jenks` automatically if you did
not name one, and rejects an explicit `antimode` with that explanation.

**A terrace flight often wants `jenks`, not the default `antimode`.** A flight is
spaced roughly linearly in height above the channel, and the log10 binning
compresses the upper treads until they merge.

**`min_band_component_pixels` is per band**: a floor on the size of one surface
on one level, not the whole-mask filter that `min_component_pixels` is in
`lsdtt-flat-surface-finder` - which is why the two names differ.

**Cleaning the valley mask can make it bigger.** `clean_valley_mask` erodes,
drops patches that do not reach the channel, then regrows — and the regrow fills
enclosed holes. It is not a filter.

**The tool needs flow routing unless you supply `hand_fname`**, so the
single-command case is not cheap. That is what the two-step workflow above is
for.

---

## Usage

```bash
lsdtt-fluvial-surfaces <path> <parameter_file>
lsdtt-fluvial-surfaces <parameter_file>      # path defaults to ./
lsdtt-fluvial-surfaces -h                    # parameter help, then exit
lsdtt-fluvial-surfaces -h <parameter_name>   # full guidance on one parameter
lsdtt-fluvial-surfaces -h html               # write an HTML/CSV reference
lsdtt-fluvial-surfaces --help-json           # machine-readable parameter list
lsdtt-fluvial-surfaces -v                    # version
lsdtt-fluvial-surfaces -c                    # citation
```

Parameter file lines are `key: value`; `#` starts a comment. Every run also
writes `<write_fname>_README.csv` and `.html` alongside its outputs, recording
exactly what it was given.

---

## Parameters

### File I/O and run control

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | `./` | Directory from which input files are read. |
| `read_fname` | string | | Input DEM, with or without extension. |
| `write_path` | string | same as `read_path` | Output directory. |
| `write_fname` | string | stem of `read_fname` | Output prefix for all products. |
| `n_threads` | int | `0` | Threads. `0` = all available. |
| `write_output_format` | string | `cog` | Raster format: `cog`, `gtiff`, `same_as_input`. |
| `write_vector_format` | string | `flatgeobuf` | Vector format: `flatgeobuf`, `gpkg`, `csv`. |
| `flat_mask_fname` | string | | Flat mask from `lsdtt-flat-surface-finder`, in `read_path`. Empty falls back to the built-in admission below. |
| `hand_fname` | string | | Precomputed height-above-channel raster, in `read_path`. Skips the HAND computation. |
| `fitted_surface_mask_fname` | string | | Optional raster restricting where the surface is fitted, in `read_path`. |
| `drainage_injection_points_fname` | string | | Vector of points injecting extra drainage area, in `read_path`. |
| `source_points_fname` | string | | Channel source points (CSV or vector), in `read_path`. |
| `source_points_input_format` | string | `auto_detect` | `auto_detect`, `csv`, `ogr_vector_layer`. |
| `source_points_layer_name` | string | | Layer name within a source point vector file. |
| `source_points_latitude_field_name` | string | `latitude` | Latitude column in a source point CSV. |
| `source_points_longitude_field_name` | string | `longitude` | Longitude column in a source point CSV. |
| `outlet_points_fname` | string | | Channel outlet points (CSV or vector), in `read_path`. |
| `outlet_points_input_format` | string | `auto_detect` | As the source format. |
| `outlet_points_layer_name` | string | | Layer name within an outlet point vector file. |
| `outlet_points_latitude_field_name` | string | `latitude` | Latitude column in an outlet point CSV. |
| `outlet_points_longitude_field_name` | string | `longitude` | Longitude column in an outlet point CSV. |

### Admission

| Parameter | Type | Default | Description |
|---|---|---|---|
| `admission_algorithm_choice` | string | `flat_mask` | `flat_mask` reads `flat_mask_fname`; `relief` and `gradient` are a coarse built-in admission for a quick look. With `flat_mask_fname` empty, `flat_mask` falls back to `relief`. |
| `relief_window_radius_metres` | float | `0` | **Radius** (m) of the square window for the built-in relief admission, converted to whole cells. `0` = 5 cells, an 11×11 window. The same name and meaning as in `lsdtt-flat-surface-finder` and `lsdtt-mountain-front`. |
| `max_relief_metres` | float | `1.0` | Local relief below which a cell is admitted as flat [m]. |
| `max_admission_gradient` | float | `0.05` | Gradient below which a cell is admitted as flat [m/m]. |
| `min_band_component_pixels` | int | `50` | Minimum patch size in pixels, **applied per band**. |

### Reference field

| Parameter | Type | Default | Description |
|---|---|---|---|
| `reference_field_algorithm_choice` | string | `hand` | Field to band: `hand` or `fitted_surface`. |
| `fitted_surface_algorithm_choice` | string | `cone` | Surface to fit: `cone`, `plane_ransac`, `poly2`. |
| `fitted_surface_ransac_tolerance_metres` | float | `5.0` | RANSAC inlier tolerance for `plane_ransac` [m]. Set it to the roughness of the surface being fitted, **not** to the depth of the features being isolated. |
| `fitted_surface_ransac_iterations` | int | `1000` | Number of RANSAC hypotheses. |
| `fitted_surface_seed` | int | `1` | RANSAC random seed, for reproducibility. |
| `fitted_surface_apex_x` | float | `0.0` | Fan apex easting in the DEM's projected CRS [m]; `cone` only. |
| `fitted_surface_apex_y` | float | `0.0` | Fan apex northing [m]; `cone` only. Leave both at 0 to take the centroid of the highest 0.1 % of fit-mask cells. |

### Ladder: level count and break finding

| Parameter | Type | Default | Description |
|---|---|---|---|
| `ladder_level_count` | int | `0` | Number of levels: `0` = auto, `2` = the valley cut, `n` = an n-level ladder. |
| `ladder_max_level_count` | int | `12` | Ceiling on the level count considered under auto selection. |
| `ladder_break_algorithm_choice` | string | `antimode` | Break finder: `gap`, `antimode` or `jenks`. |
| `ladder_gap_break_position` | string | `midpoint` | Where in the empty span the break goes: `lower`, `midpoint`, `upper`. |
| `ladder_gap_rank` | string | `lowest` | Order the gap finder offers candidates in: `lowest` or `width`. |
| `ladder_gap_trough_max_ratio` | float | `0.6` | Offer a local minimum as a candidate when its floor is below this fraction of the smaller flanking peak. `0` = empty spans only. |
| `ladder_gap_empty_fraction` | float | `0.02` | A histogram bin below this fraction of the tallest bin counts as empty. |
| `ladder_gap_min_width_metres` | float | `0.5` | Narrowest empty span the gap finder will call a separation [m]. |
| `ladder_gap_hist_bins` | int | `80` | Linear histogram bins the gap finder searches for empty runs. |
| `ladder_max_separation_ratio` | float | `0.0` | Reject a partition whose worst break sits where the density is above this fraction of the populations either side. `0` disables. |
| `ladder_min_break_metres` | float | `0.0` | Reject any break below this height [m]. `0` = unbounded. |
| `ladder_max_break_metres` | float | `0.0` | Reject any break above this height [m]. `0` = unbounded. |
| `ladder_hist_bins` | int | `64` | Histogram bin count for the antimode break finder. |
| `ladder_max_valley_fraction` | float | `0.6` | Antimode acceptance: dip depth relative to the smaller flanking peak. |
| `ladder_min_field_metres` | float | `0.05` | Log floor for the antimode histogram [m]. |
| `ladder_k_bootstrap` | int | `0` | Resample each window and re-run the break search this many times, to report how stable the chosen level count is. `0` = off. |
| `ladder_k_bootstrap_seed` | int | `1` | Seed for the k bootstrap. |

### Ladder: resolvability and windowing

| Parameter | Type | Default | Description |
|---|---|---|---|
| `vertical_accuracy_metres` | float | `0.0` | Vertical accuracy of the DEM, `E_v` [m]. **`0` disables the floor** and every window then reports `resolved` with no instrument guarantee. |
| `ladder_resolvability_kappa` | float | `3.0` | Multiplier on `E_v`. Required band-centre gap = `kappa * E_v`. |
| `vertical_accuracy_algorithm_choice` | string | `fixed` | How `E_v` is arrived at: `fixed` or `slope_dependent`. |
| `vertical_accuracy_gentle_metres` | float | `0.0` | `E_v` [m] for windows at or below the break gradient. `slope_dependent` only. |
| `vertical_accuracy_steep_metres` | float | `0.0` | `E_v` [m] for windows above the break gradient. `slope_dependent` only. |
| `vertical_accuracy_break_gradient` | float | `0.2` | Gradient [m/m] separating gentle from steep `E_v`; 0.2 is the 20 percent grade the TanDEM-X specification splits on. |
| `vertical_accuracy_gradient_window_radius_metres` | float | `0` | Window radius [m] of the polynomial fit the window gradient is measured over. `0` = five cells. |
| `ladder_reach_window_metres` | float | `0.0` | Along-channel window length for reach-relative breaks [m]. `0` = auto, negative = one global ladder. Ignored for `fitted_surface`. |
| `ladder_min_window_pixels` | int | `1000` | Minimum admitted pixels for a window to find its own breaks. |

### Classification

| Parameter | Type | Default | Description |
|---|---|---|---|
| `max_floodplain_hand_metres` | float | `1.0` | Upper HAND bound of the floodplain class [m]. An optional hard cap; the ladder already decides where a level ends. The same name as in `lsdtt-hillslope-channel-coupling`. |
| `min_terrace_hand_metres` | float | `0.0` | Lower HAND bound of the terrace class [m]. |
| `channel_adjacency_radius_metres` | float | `0.0` | Distance within which a surface counts as channel-adjacent [m]. |
| `valley_corridor_half_width_metres` | float | `5000.0` | Valley-corridor half-width [m]. `0` or negative = no limit. |
| `clip_corridor_per_pixel` | bool | `true` | Apply the corridor limit per pixel as well as per surface. |
| `max_band_spread_metres` | float | `50.0` | QC tripwire on the field spread within one surface [m]. `0` disables. Under band-then-label it should almost never fire. |

### Flow routing and channel scope

| Parameter | Type | Default | Description |
|---|---|---|---|
| `fill_before_routing` | bool | `true` | Choose whether to fill pits before building flow routing. `false` routes on the unfilled DEM, so paths may terminate at sinks. |
| `fill_min_slope_for_routing` | float | `0.0001` | Minimum imposed slope when filling [m/m]. |
| `boundary_conditions` | string | `oooo` | Edge boundary codes, north east south west. |
| `threshold_contributing_pixels` | int | `1000` | Contributing-pixel threshold defining a channel. |
| `break_flats_with_channels` | bool | `false` | Remove in-scope channel pixels from the admitted set before banding. |
| `connect_surfaces_to_channels` | bool | `false` | Record which in-scope channel each surface is linked to. |
| `auto_detect_outlet_points` | bool | `true` | Infer source and outlet roles from a single combined point file. |
| `snap_sources_to_channel_threshold` | bool | `true` | Snap supplied source points onto the threshold channel network. |
| `source_snap_max_search_radius_pixels` | int | `20` | Search radius when snapping a source point [pixels]. |
| `outlet_snap_max_distance_metres` | float | `500.0` | Maximum distance an outlet may be moved onto a traced channel [m]. |

### Valley floor

| Parameter | Type | Default | Description |
|---|---|---|---|
| `force_valley_hand_metres` | float | `0.0` | Define the valley floor directly as height above channel below this [m], instead of using the ladder's break. `0` = use the ladder. |
| `force_valley_max_distance_metres` | float | `0.0` | Euclidean distance to an in-scope channel bounding the forced valley floor [m]. `0` = no limit. |

### Valley-mask cleaning

All of these affect `<stem>_valley_mask.tif` only. The band index, the labels and
the patch table keep their holes and their speckle.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `fill_valley_mask_holes` | bool | `true` | Fill enclosed holes before writing the mask. A hole in the valley floor is a hole in the domain the centreline is traced through. |
| `clean_valley_mask` | bool | `true` | Erode, drop patches that do not reach the channel, then regrow. **The regrow fills enclosed holes, so cleaning can make the mask bigger.** |
| `valley_mask_closing_width_fraction` | float | `0.0` | Set the closing radius as a fraction of the *measured* valley width. `0` = use the radius in metres instead. About `0.2` is a reasonable start. |
| `valley_width_section_spacing_metres` | float | `200.0` | Spacing of the cross sections used to measure the valley width [m]. |
| `valley_width_max_half_length_metres` | float | `3000.0` | How far each width section looks either side before giving up [m]. |
| `valley_mask_closing_radius_metres` | float | `0.0` | Seal narrow gaps with a morphological closing [m]. `0` = off. |
| `valley_dilation_pixels` | int | `0` | Erosion / regrow radius when cleaning [pixels]. `0` = no erosion. |
| `valley_max_detour_ratio` | float | `0.0` | Maximum along-ground path to the channel, as a multiple of the straight-line distance. |
| `restrict_valley_to_channel_connected` | bool | `true` | Drop valley patches with no pixel on or beside the channel. |
| `restrict_valley_to_line_of_sight` | bool | `false` | Require a straight line from each valley pixel to the channel to stay inside the mask. **A strong filter** — on the Royalton 10 m DEM it removed 82 % of band 0 on its own. Use it where two valleys meet around a col and connectivity alone merges them. |

### Output switches

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_reference_field_raster` | bool | `true` | Write the field the ladder actually ran on. |
| `write_hand_raster` | bool | `false` | Write the HAND raster, for feeding back as `hand_fname`. |
| `write_band_index_raster` | bool | `true` | Write the per-pixel band index. |
| `write_banded_labels_raster` | bool | `true` | Write the per-pixel patch label. |
| `write_valley_mask_raster` | bool | `true` | Write band 0 as a valley mask. |
| `write_surface_class_raster` | bool | `true` | Write the floodplain / terrace class raster. `hand` field only. |
| `write_floodplain_mask_raster` | bool | `false` | Write the floodplain class on its own. |
| `write_surface_patches_csv` | bool | `true` | One row per labelled patch. The join key for everything downstream. |
| `write_surface_class_csv` | bool | `true` | One row per classification record. |
| `write_ladder_windows_csv` | bool | `true` | The per-window resolvability report. **Read it first.** |
| `write_ladder_break_metric_csv` | bool | `false` | Write the break *decision*, not just the break. |
| `write_ladder_histogram_csv` | bool | `false` | Write the reference-field density each window's breaks were chosen from. |
| `ladder_hist_plot_bins` | int | `80` | Bin count for that per-window density. |
| `write_ladder_window_id_raster` | bool | `false` | Which along-channel window each pixel's breaks came from. |
| `write_reference_field_profile_csv` | bool | `false` | Per-pixel table of the reference field along the channel. |
| `write_nearest_source_id_raster` | bool | `false` | Which channel each pixel's height is measured to. |
| `write_patch_terrain_aspect` | bool | `false` | Add terrain aspect and its resultant length to the patch table. |
| `patch_aspect_window_radius_metres` | float | `0` | Window radius [m] of the polynomial fit the aspect is measured over. `0` = five cells. |
| `write_channel_scope_raster` | bool | `false` | The in-scope channel network as a raster. |
| `write_channel_scope_nodes_csv` | bool | `false` | The in-scope channel nodes as a table. |
| `write_point_snap_csv` | bool | `false` | Where each supplied point snapped to. |
| `write_channel_network` | bool | `false` | Write the channel network as vector layers. |
| `write_channel_network_lines` | bool | `true` | Compact default representation of the channel network used to fit the surfaces; per-node points remain opt-in. |
| `write_channel_network_points` | bool | `false` | Channel nodes as points. |
| `write_junction_network_points` | bool | `false` | Junctions as points. |

---

## See also

- [Floodplains, terraces and valleys](../design_and_algorithms/fluvial_surfaces_and_valleys.md) — the
  band-then-label ordering, the three break finders, the resolvability floor,
  and where the methods come from
- [CLI reference: lsdtt-flat-surface-finder](cli_reference_lsdtt_flat_surface_finder.md)
- [CLI reference: lsdtt-valley-centreline](cli_reference_lsdtt_valley_centreline.md)

### Renamed parameters (2026-10-26)

Old keys still run, with a warning naming the replacement; converted values are
written to `_ingested_params.param`.

| Old key | Now |
|---|---|
| `fallback_admission_algorithm_choice` | `admission_algorithm_choice` — still only applies when `flat_mask_fname` is empty; `slope` is now `gradient` |
| `fallback_relief_window_radius_pixels`, `fallback_relief_window_pixels` | `relief_window_radius_metres` — converted with the cell size |
| `fallback_relief_threshold_metres` | `max_relief_metres` |
| `fallback_slope_threshold` | `max_admission_gradient` |
| `min_component_pixels` | `min_band_component_pixels` |
| `vertical_accuracy_mode` | `vertical_accuracy_algorithm_choice` |
| `vertical_accuracy_break_degrees` | `vertical_accuracy_break_gradient` — converted with gradient = tan(angle) |
| `vertical_accuracy_slope_window_radius_metres` | `vertical_accuracy_gradient_window_radius_metres` |
| `floodplain_max_height_metres` | `max_floodplain_hand_metres` |
| `terrace_min_height_metres` | `min_terrace_hand_metres` |

The `_ladder_windows.csv` column `median_slope_degrees` is now `median_gradient`
(m/m).

<!-- BEGIN GENERATED CLI PARAMETERS -->
## Parameters

This section is generated from the executable's `-h json` registry. Do not edit it manually.

**Status: experimental.** This CLI is not yet fully tested.

### File I/O

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `read_path` | `string` | `./` | `io` | - | Directory from which input files are read. | Defaults to the directory holding the parameter file. |
| `read_fname` | `string` | - | `io` | - | Input filename, with or without extension (e.g. my_dem.tif). | Read from read_path. lsdtt3 uses GDAL, so any GDAL-readable raster format works. If empty, defaults to the parameter filename stem. |
| `write_path` | `string` | - | `io` | - | Directory to which output files are written. | Defaults to read_path when not set. |
| `write_fname` | `string` | - | `io` | - | Output filename prefix (no path, no extension). | All outputs are named &lt;write_fname&gt;_&lt;descriptor&gt;.&lt;ext&gt;. Defaults to the stem of read_fname when not set. |

### General

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `n_threads` | `int` | `0` | `io` | - | Number of threads. 0 = all available. | OpenMP thread count for every parallel stage. |
| `write_output_format` | `string` | `cog` | `io` | - | Raster output format: cog, gtiff or same_as_input. | cog is a Cloud Optimised GeoTIFF and is the project default. |
| `write_vector_format` | `string` | `flatgeobuf` | `io` | - | Vector output format: flatgeobuf, gpkg or csv. | Governs every vector layer this tool writes, including the channel network. FlatGeobuf is preferred: it is compact, carries its own georeferencing, and does not have GeoPackage's problems writing across drives (a common WSL2 situation). |

### Method

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `fluvial_surface_algorithm_choice` | `string` | `ladder` | `switch` | - | Method preset: ladder or clubb2017. | Fills every stage parameter left at `preset`. ladder (default): the admitted ground (flat_mask_fname, or the built-in relief admission) is banded by a data-driven ladder on HAND, components labelled within each band. clubb2017: Clubb et al. (2017, Earth Surf. Dynam. 5:369-385) - floodplain where HAND above a Strahler-order channel and gradient are both below bounds read off quantile-quantile plots, floodplain patches must contain a channel pixel, terraces above the floodplain bound measured from the nearest channel (swath baseline). Any stage named in the parameter file overrides the preset, so e.g. clubb2017 with level_algorithm_choice: ladder bands the Clubb admission with the ladder. Stage values (ladder / clubb2017): admission_algorithm_choice flat_mask / floodplain_thresholds; floodplain_threshold_algorithm_choice fixed / qq_gaussian; channel_definition_algorithm_choice contributing_pixels / stream_order; level_algorithm_choice ladder / fixed_thresholds; terrace_reference_field_algorithm_choice same / swath_baseline. |
| `level_algorithm_choice` | `string` | `preset` | `param` | `fluvial_surface_algorithm_choice` | How surfaces are put on levels: preset, ladder or fixed_thresholds. | ladder: data-driven breaks in the reference field, as many levels as the DEM can resolve (the ladder_* parameters). fixed_thresholds: two levels from two bounds - floodplain below max_floodplain_hand_metres (and touching the channel), terrace above it by min_terrace_height_above_floodplain_metres - with the bounds chosen by floodplain_threshold_algorithm_choice. fixed_thresholds needs reference_field_algorithm_choice hand. |
| `channel_definition_algorithm_choice` | `string` | `preset` | `param` | `fluvial_surface_algorithm_choice` | What counts as a channel: preset, contributing_pixels or stream_order. | contributing_pixels: the whole network down from heads at threshold_contributing_pixels. stream_order: only its pixels of Strahler order &gt;= min_channel_stream_order; HAND, the floodplain contact test and the swath baseline then all use that channel. The threshold still decides which first-order channels exist. |
| `min_channel_stream_order` | `int` | `3` | `param` | `channel_definition_algorithm_choice` | Lowest Strahler order counted as channel (stream_order only). | Clubb et al. (2017) used third order at all eight of their sites. LSDTopoTools2 called this threshold_SO. |
| `terrace_reference_field_algorithm_choice` | `string` | `preset` | `param` | `fluvial_surface_algorithm_choice` | What terrace height is measured from: preset, same or swath_baseline. | same: the main reference field (HAND). swath_baseline: elevation above the NEAREST channel pixel within swath_half_width_metres, which does not follow drainage, so a terrace above a confluence is referenced to the nearest channel rather than the one its flow path reaches (Clubb et al. 2017). fixed_thresholds only. |
| `swath_half_width_metres` | `float` | `500` | `param` | `terrace_reference_field_algorithm_choice` | Half-width [m] of the swath baseline: terraces further than this from the channel are not measured. | A HALF-width, as half_width_metres in lsdtt-swath-profile. LSDTopoTools2's swath_width (default 1000) was a full width. |

### Admission

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `flat_mask_fname` | `string` | - | `io` | - | Flat mask raster from lsdtt-flat-surface-finder, in read_path. | Cells &gt; 0.5 are admitted. This is normally &lt;stem&gt;_combined_flat_mask.tif from lsdtt-flat-surface-finder, which is the intended two-step workflow: detect flat ground there, sort it into levels here. Leave empty to use the built-in relief or gradient admission (admission_algorithm_choice) instead, which is coarser and is meant for a quick look rather than for a result you intend to publish. |
| `admission_algorithm_choice` | `string` | `preset` | `param` | - | Which ground is admitted: preset, flat_mask, relief, gradient or floodplain_thresholds. | floodplain_thresholds (the clubb2017 preset): cells whose gradient, read over floodplain_gradient_window_radius_metres, is below max_floodplain_gradient as chosen by floodplain_threshold_algorithm_choice, among cells with a defined HAND. flat_mask (the ladder preset) reads flat_mask_fname, normally &lt;stem&gt;_combined_flat_mask.tif from lsdtt-flat-surface-finder. relief and gradient are a coarse built-in admission for a quick look when there is no flat mask: relief admits cells whose local relief over relief_window_radius_metres is below max_relief_metres; gradient admits cells whose gradient is below max_admission_gradient. flat_mask with flat_mask_fname empty falls back to relief, as earlier versions did. Replaces fallback_admission_algorithm_choice (whose `slope` is now `gradient`). |
| `relief_window_radius_metres` | `float` | `0` | `param` | `admission_algorithm_choice` | Radius (m) of the square window local relief (max - min) is measured over. 0 = 5 cells. | Converted to the nearest whole number of cells, at least one; the run log states the radius actually used. A radius r cells is a (2r+1) x (2r+1) window, so the default is a 11 x 11 window. Replaces relief_window_side_pixels (a full side) and fallback_relief_window_radius_pixels (a radius in cells), which spelled the same idea two ways. |
| `max_relief_metres` | `float` | `1` | `param` | `admission_algorithm_choice` | Local relief (m) below which a pixel is admitted as flat. | Read a value off the relief CDF the tool writes. Where the tool has relief_threshold_algorithm_choice, this is the value used by `fixed`. |
| `max_admission_gradient` | `float` | `0.05` | `param` | `admission_algorithm_choice` | Gradient (m/m) below which a pixel is admitted as flat. | A gradient, not an angle: 0.087 is 5 degrees, 0.176 is 10, 0.268 is 15 (gradient = tan(angle)). Measured from a fitted polynomial surface rather than a 3x3 difference. |
| `min_band_component_pixels` | `int` | `50` | `param` | - | Minimum patch size in pixels, applied PER BAND. | Components are formed inside each band, so this is a floor on the size of one surface on one level - not the whole-mask filter that min_component_pixels is in lsdtt-flat-surface-finder, which is why the name differs. Replaces min_component_pixels here. |

### Reference field

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `reference_field_algorithm_choice` | `string` | `hand` | `switch` | - | Field to band: hand or fitted_surface. | hand = height above the in-scope channel network; the right field for a terrace flight, where treads sit at a constant height above the channel. fitted_surface = the signed residual to a surface fitted to the flat ground; the right field for a fan, which climbs continuously in HAND and would be cut into meaningless arcs by a HAND ladder. |
| `hand_fname` | `string` | - | `io` | - | Precomputed height-above-channel raster, in read_path. | THE COMPOSABILITY LEVER. HAND is the expensive stage of this tool and the level count is what you iterate on, so compute HAND once, write it with write_hand_raster, then re-run with hand_fname set and try different ladder settings for the cost of the ladder alone. Flow routing is still built when a channel scope or network output is requested, but the HAND computation is skipped. |
| `fitted_surface_algorithm_choice` | `string` | `cone` | `param` | `reference_field_algorithm_choice` | Surface to fit: cone, plane_ransac or poly2. | cone is the fan option and the default: z falling linearly with distance from an apex. A fan is radial, and a cone is not a low-order polynomial in x and y - on the synthetic fan fixture the older surface's own residual spread is about 34 m under a RANSAC plane and 10 m under a poly2, against a 4 m lobe inset, so the misfit swamps the signal; under a cone it is about 0.3 m. Use plane_ransac for a piedmont or a terrace tread and poly2 for a broadly curved surface with no single centre. |
| `fitted_surface_mask_fname` | `string` | - | `io` | `reference_field_algorithm_choice` | Optional raster restricting where the surface is fitted, in read_path. | Cells &gt; 0.5 take part in the fit. Use it to fit on the older fan surface alone. Optional: a cone fitted to everything, inset lobes included, still separates the populations - the lobes drag the whole fit down rather than blurring the boundary, and the ladder does not care about the absolute offset. Empty = fit over the admitted mask. |
| `fitted_surface_ransac_tolerance_metres` | `float` | `5.0` | `param` | `reference_field_algorithm_choice` | RANSAC inlier tolerance for plane_ransac [m]. | Set this to roughly the roughness of the surface being fitted, NOT to the depth of the features being isolated: a tolerance as large as a lobe inset lets the lobes join the fit and flattens the signal. |
| `fitted_surface_ransac_iterations` | `int` | `1000` | `param` | `reference_field_algorithm_choice` | Number of RANSAC hypotheses for plane_ransac. | Increase when the inlier fraction is expected to be below about 30%. |
| `fitted_surface_seed` | `int` | `1` | `param` | `reference_field_algorithm_choice` | RANSAC random seed, for reproducibility. | Same seed and same input give the same plane. |
| `fitted_surface_apex_x` | `float` | `0.0` | `param` | `reference_field_algorithm_choice` | Fan apex easting in the DEM's projected CRS [m]; cone only. | Leave both apex coordinates at 0 to locate the apex automatically, as the centroid of the highest 0.1% of fit-mask cells. The apex is supplied rather than solved for because it is a strongly non-linear parameter and a user can point at a fan head more reliably than a solver can find it. |
| `fitted_surface_apex_y` | `float` | `0.0` | `param` | `reference_field_algorithm_choice` | Fan apex northing in the DEM's projected CRS [m]; cone only. | See fitted_surface_apex_x. |

### Ladder

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `ladder_level_count` | `int` | `0` | `param` | - | Number of levels: 0 = auto, 2 = the valley cut, n = an n-level ladder. | 0 lets the resolvability floor choose the largest level count the DEM can actually carry. 2 gives the valley/wall separation that lsdtt-valley-centreline consumes. An explicit count is honoured and reported honestly: if that partition is not resolvable the window still records the count you asked for, with status unresolved, rather than silently substituting a different one. |
| `ladder_max_level_count` | `int` | `12` | `param` | - | Ceiling on the level count considered under auto selection. | A ceiling on what one reach can plausibly carry, not a storage limit. A reach resolving more than a dozen distinct fluvial levels is more likely over-fitting noise than recording twelve terraces. Raise it freely for a genuine long Quaternary staircase; there is no numeric cliff behind it. |
| `ladder_break_algorithm_choice` | `string` | `antimode` | `param` | - | Break finder: gap, antimode or jenks. | antimode finds dips in a smoothed log10 density and is shared with the valley engine. jenks is Fisher-Jenks on the linear field; it takes the level count directly and is the source of the reported goodness-of-variance fit. NOTE: antimode bins on log10 and therefore cannot band a field with negative values, which every fitted_surface residual has. This tool switches to jenks automatically when you select fitted_surface without naming a break algorithm, and rejects the combination if you ask for it explicitly. |
| `ladder_gap_break_position` | `string` | `midpoint` | `param` | `ladder_break_algorithm_choice` | Where in the empty span between clusters the break goes: lower, midpoint or upper. | The span is nearly empty by construction, so any position separates the same two clusters - but the sparse pixels INSIDE it (the bank, the riser, the floodplain margin) go to whichever side the break leaves them on. Use `upper` when the low band IS the valley floor and you want the whole margin in it; `lower` for the strictest reading of the active floodplain; `midpoint` splits the difference. |
| `ladder_gap_rank` | `string` | `lowest` | `param` | `ladder_break_algorithm_choice` | Order the gap finder offers its candidate breaks in: lowest or width. | gap algorithm only, and it decides which break wins when more than one is admissible. `lowest` walks the ladder from the bottom, so the valley cut is taken before any terrace or hillslope cut. `width` prefers the widest empty span, which is only a good guide where the whole field range is a few metres: in a mountain catchment the widest empty span sits in the sparse tail a hundred metres above the river, and taking it bands the entire terrace staircase as valley floor. |
| `ladder_gap_trough_max_ratio` | `float` | `0.6` | `param` | `ladder_break_algorithm_choice` | Offer a local minimum as a candidate break when its floor is below this fraction of the smaller flanking peak. 0 = empty spans only. | gap algorithm only. `ladder_gap_empty_fraction` measures emptiness against the TALLEST bin, so where one level dominates the population a real dip above it never reads as empty and the finder has no candidate anywhere near the valley floor. This offers those dips as candidates; `ladder_max_separation_ratio` still has to accept them before any becomes a break, which is why this is the looser of the two numbers. |
| `ladder_gap_empty_fraction` | `float` | `0.02` | `param` | `ladder_break_algorithm_choice` | A histogram bin below this fraction of the tallest bin counts as empty. | gap algorithm only. Raise it to treat sparsely-populated ground as a gap; lower it to require a near-total absence of pixels before calling a space between clusters. |
| `ladder_gap_min_width_metres` | `float` | `0.5` | `param` | `ladder_break_algorithm_choice` | Narrowest empty span the gap finder will call a separation [m]. | gap algorithm only. Below this a run of empty bins is noise, not the space between two surfaces. |
| `ladder_gap_hist_bins` | `int` | `80` | `param` | `ladder_break_algorithm_choice` | Linear histogram bins the gap finder searches for empty runs. | Too few and a real gap is bridged by a single populated bin; too many and ordinary sampling noise fragments a cluster into false gaps. |
| `ladder_max_separation_ratio` | `float` | `0.0` | `param` | - | Reject a partition whose worst break sits where the density is above this fraction of the populations either side. 0 disables. | A SECOND CONSTRAINT alongside the resolvability floor, never blended with it into a score. The floor asks whether the DEM can MEASURE a separation; this asks whether there IS one. Both are necessary and neither implies the other: two lidar terraces 30 cm apart pass this and need a small vertical_accuracy_metres to pass the floor, while a 26 m cut through one broad mode passes any floor and fails this outright. 0.6 matches the antimode finder's own acceptance limit. Off by default because turning it on changes which partitions are admissible. |
| `ladder_min_break_metres` | `float` | `0.0` | `param` | - | Reject any break below this height [m]. 0 = unbounded. | A hard bound applied to every candidate before selection, for when the lowest plausible level boundary is known independently. |
| `ladder_max_break_metres` | `float` | `0.0` | `param` | - | Reject any break above this height [m]. 0 = unbounded. | The companion to ladder_min_break_metres, and the direct lever against a break landing far above the valley in a landscape whose flat uplands dominate the population. |
| `ladder_k_bootstrap` | `int` | `0` | `param` | - | Resample each window and re-run the break search this many times, to report how stable the chosen level count is. 0 = off. | `k_selected` on its own says nothing about whether a different draw of the same landscape would have given the same answer, and a ladder that returns k = 7 every time is not the same claim as one that returns it once in fifteen. This costs B times the break search per window, and the break search is where the cost sits, not the pixel count - every finder bins before it searches. Measured on a 3.58 M cell DEM with 41 windows, 15 replicates added 4.9 s under `gap` (16.4 to 21.3 s). The same DEM's whole ladder is 6.2 s under `antimode` and 37.7 s under `jenks`, so scale accordingly: affordable on gap, think first on jenks. Results go to _ladder_windows.csv as k_bootstrap_mode, k_stability, k_ci_lo, k_ci_hi and n_bootstrap. |
| `ladder_k_bootstrap_seed` | `int` | `1` | `param` | `ladder_k_bootstrap` | Seed for the k bootstrap, so a stability figure is reproducible. | Each window is seeded from this plus its own source and index, so a window's stability does not shift because a neighbouring reach changed. The generator is specified in the source rather than taken from the standard library, so the number reproduces across platforms. |
| `vertical_accuracy_algorithm_choice` | `string` | `fixed` | `param` | - | How E_v is arrived at: fixed or slope_dependent. | Photogrammetric and InSAR DEMs degrade on slopes - TanDEM-X relative LE90 is roughly 2 m on gentle ground and 4 m on steep - so a single number for a raster spanning a valley floor and a 40-degree wall states an accuracy that is wrong in both places. Under slope_dependent each along-channel WINDOW reads one E_v off the median gradient of its own admitted cells. Per window, not per pixel: the resolvability floor constrains the separation between band CENTRES, and a band centre is a property of a partition rather than of a place, so it has no slope of its own. The value each window used is written to _ladder_windows.csv. |
| `vertical_accuracy_gentle_metres` | `float` | `0.0` | `param` | - | E_v [m] for windows at or below the break gradient. slope_dependent only. | For TanDEM-X this is the quoted relative LE90 on gentle ground, about 2 m. |
| `vertical_accuracy_steep_metres` | `float` | `0.0` | `param` | - | E_v [m] for windows above the break gradient. slope_dependent only. | About 4 m for TanDEM-X. It must be the LARGER of the two or the mode is doing nothing useful. |
| `vertical_accuracy_break_gradient` | `float` | `0.2` | `param` | - | Gradient [m/m] separating the gentle and steep E_v. | 0.2 is a 20 percent grade, which is where the TanDEM-X accuracy specification splits. Replaces vertical_accuracy_break_degrees; an old value is converted with gradient = tan(angle). |
| `vertical_accuracy_gradient_window_radius_metres` | `float` | `0` | `param` | - | Window radius [m] of the polynomial fit the window gradient is measured over. 0 = five cells. slope_dependent only. | Stated in cells by default because a length tuned on lidar is sub-pixel on 30 m COP30. This is the same measurement lsdtt-flat-surface-finder makes for max_admission_gradient, and it is worth using the same window in both so the gradient a surface was admitted on is the gradient its E_v was read from. |
| `vertical_accuracy_metres` | `float` | `0.0` | `param` | - | Vertical accuracy of the DEM, E_v [m]. 0 disables the floor. | Two levels are only claimed as distinct when their band centres are at least ladder_resolvability_kappa * E_v apart. This is metadata you must supply: it cannot be derived from the DEM. 0 disables the check and every window then reports resolved, carrying no instrument guarantee at all. Note that a DEM whose vertical units are feet on a metre grid has an effective E_v 3.28 times larger than you think. |
| `ladder_resolvability_kappa` | `float` | `3.0` | `param` | - | Multiplier on E_v in the resolvability floor. | Required band-centre gap = kappa * vertical_accuracy_metres. |
| `ladder_reach_window_metres` | `float` | `0.0` | `param` | - | Along-channel window length for reach-relative breaks [m]. | 0 = auto: the larger of 2 km and an eighth of each channel's length. Negative = one global partition over the whole population. Windowing lets a break level track downstream changes in incision instead of being averaged into a single basin-wide ladder. Ignored for the fitted_surface field, which has no along-channel coordinate. |
| `ladder_min_window_pixels` | `int` | `1000` | `param` | - | Minimum admitted pixels for a window to find its own breaks. | A sparser window falls back to the global partition rather than inventing structure from too few pixels. |
| `ladder_hist_bins` | `int` | `64` | `param` | - | Histogram bin count for the antimode break finder. | Clamped to at least 16. |
| `ladder_max_valley_fraction` | `float` | `0.6` | `param` | - | Antimode acceptance: dip depth relative to the smaller flanking peak. | A dip is accepted as a break when its smoothed density is below this fraction of the smaller of its two flanking modal peaks. |
| `ladder_min_field_metres` | `float` | `0.05` | `param` | - | Log floor for the antimode histogram [m]. | Field values below this are clamped before taking log10. |

### Outputs

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_patch_terrain_aspect` | `bool` | `false` | `output` | - | Add the terrain aspect and its resultant length to the patch table. | The existing `aspect_degrees` column is the direction the BANDED FIELD climbs, fitted as a plane through it - on a HAND field that is roughly `away from the channel`, not the direction the ground faces. This adds the topographic aspect as a CIRCULAR mean (the arithmetic mean of 350 and 10 degrees is 180, the opposite way) together with `aspect_r_bar`, the mean resultant length in [0, 1]. R_bar is the number that says whether the mean direction is worth quoting: 1 means every cell in the patch faces the same way, 0 means they face every way equally and the aspect is meaningless. Costs one polynomial fit. |
| `patch_aspect_window_radius_metres` | `float` | `0` | `output` | `write_patch_terrain_aspect` | Window radius [m] of the polynomial fit the terrain aspect is measured over. 0 = five cells. | Stated in cells by default because a length tuned on lidar is sub-pixel on 30 m COP30. |
| `write_floodplain_threshold_csv` | `bool` | `true` | `output` | `floodplain_threshold_algorithm_choice` | Write the table each data-chosen floodplain bound was read from. | Output files: &lt;stem&gt;_floodplain_hand_qq.csv and &lt;stem&gt;_floodplain_gradient_qq.csv (qq_gaussian), or _adaptive_cdf.csv (cdf_divergence). A threshold is only interpretable against the distribution it was cut from. |
| `write_terrace_field_raster` | `bool` | `false` | `output` | `terrace_reference_field_algorithm_choice` | Write the terrace reference field (fixed_thresholds). | Output file: &lt;stem&gt;_terrace_field.tif - swath-baseline relief, or HAND under same. |
| `write_floodplain_gradient_raster` | `bool` | `false` | `output` | `admission_algorithm_choice` | Write the gradient raster (m/m) the floodplain bound was applied to. | Output file: &lt;stem&gt;_floodplain_gradient.tif. |
| `write_reference_field_raster` | `bool` | `true` | `switch` | - | Write the field the ladder actually ran on. | Output file: &lt;stem&gt;_reference_field.tif. Whichever field was selected, so a run is self-documenting. |
| `write_hand_raster` | `bool` | `false` | `switch` | - | Write the height-above-channel raster. | Output file: &lt;stem&gt;_hand.tif. Feed it back as hand_fname to iterate on the ladder without recomputing HAND. |
| `write_band_index_raster` | `bool` | `true` | `switch` | - | Write the per-pixel band index. | Output file: &lt;stem&gt;_band_index.tif (Int16, nodata -1). The authoritative statement of which level each pixel is on; the band is NOT packed into the patch label. |
| `write_banded_labels_raster` | `bool` | `true` | `switch` | - | Write the per-pixel patch label. | Output file: &lt;stem&gt;_banded_labels.tif (UInt32). A global running counter 1..N assigned band by band; 0 is background. Join it to &lt;stem&gt;_surface_patches.csv on the label. |
| `write_valley_mask_raster` | `bool` | `true` | `switch` | - | Write band 0 as a valley mask. | Output file: &lt;stem&gt;_valley_mask.tif. This is the input lsdtt-valley-centreline expects as valley_mask_fname. With ladder_level_count 2 it is exactly the valley/wall cut. |
| `valley_mask_closing_width_fraction` | `float` | `0.0` | `param` | `clean_valley_mask` | Set the closing radius as a fraction of the measured valley width. 0 = use valley_mask_closing_radius_metres instead. | A radius in metres does not travel: a value that seals the slivers in a 400 m floodplain erases a 100 m one, so it has to be re-guessed for every DEM. This measures the width instead - cross sections cast perpendicular to each channel, the CUMULATIVE length of valley cells along each (so a bar or an island does not count as floor), median across sections - and takes the radius as this fraction of it. About 0.2 is a reasonable starting point. Overrides valley_mask_closing_radius_metres when set. |
| `valley_width_section_spacing_metres` | `float` | `200.0` | `param` | `valley_mask_closing_width_fraction` | Spacing of the cross sections used to measure the valley width [m]. | Only used when valley_mask_closing_width_fraction is set. Closer spacing gives a better-sampled median at more cost. |
| `valley_width_max_half_length_metres` | `float` | `3000.0` | `param` | `valley_mask_closing_width_fraction` | How far each cross section looks either side before giving up [m]. | A section reaching this on BOTH sides is discarded rather than counted, because it is measuring the search limit and not the valley. Set it comfortably wider than the widest valley expected. |
| `valley_mask_closing_radius_metres` | `float` | `0.0` | `param` | `clean_valley_mask` | Seal narrow gaps in the valley mask with a morphological closing [m]. 0 = off. | A dilation followed by an erosion of the same radius. It bridges thin OPEN-ENDED gaps - the parallel slivers left by scroll bars, levees, abandoned channels and field boundaries lying a little above the surrounding floor - which hole filling cannot touch, because a finger open at one end is not an enclosed hole. Because the erosion undoes the dilation, the outer extent is left about where it was rather than grown. Set it to roughly the width of the widest sliver you want closed; too large and it will bridge across a genuine terrace riser. |
| `fill_valley_mask_holes` | `bool` | `true` | `switch` | `write_valley_mask_raster` | Fill enclosed holes in the valley mask before writing it. | A hole in the valley floor - a rise the ladder put in a higher band, a patch the size filter dropped - is a hole in the DOMAIN the centreline is traced through, and the least-cost walk goes around it. The centreline then hugs one wall instead of running down the middle. Filling enclosed holes is what makes the mask a valley rather than a valley with islands. Only affects &lt;stem&gt;_valley_mask.tif; the band index and the patch table keep the holes, because there the rise is a real observation. |
| `clean_valley_mask` | `bool` | `true` | `switch` | `write_valley_mask_raster` | Tidy band 0 before writing it as the valley mask. | Band 0 straight off the ladder is speckled and contains flat patches that have nothing to do with the channel, so it is not a drop-in for lsdtt-valley-centreline until it has been cleaned the way that tool used to clean its own. Erode, drop patches that do not reach the channel, then regrow. Only affects &lt;stem&gt;_valley_mask.tif; the band index, the labels and the patch table are untouched. |
| `valley_dilation_pixels` | `int` | `0` | `param` | `clean_valley_mask` | Erosion / regrow radius used when cleaning the valley mask [pixels]. 0 = no erosion. | DEFAULT 0, and erosion is usually the wrong tool here. It exists to SEVER thin threads of flat ground joining one valley floor to a neighbour's across a col, so that connected-component filtering can tell them apart - it is a means of making connectivity discriminating, not a way to shrink the mask. The cost is that a reach narrower than 2 * this + 1 pixels erodes to nothing and cannot regrow: on the Royalton 10 m DEM a value of 4 (64 m) destroyed 30% of a valley floor that connectivity alone had got right. Raise it only when two valleys are genuinely merging through a col. |
| `valley_max_detour_ratio` | `float` | `0.0` | `param` | `clean_valley_mask` | Maximum along-ground path to the channel, as a multiple of the straight-line distance. | A flat that sits NEAR the channel but is only reachable the long way round is a different landform, not this channel's valley floor. 0 or negative disables the test, which is the DEFAULT: on a meandering floodplain the along-ground path around a bend legitimately exceeds the straight line several times over, and the test removes real valley floor. Turn it on where flats from separate drainages are connected by a long corridor. |
| `restrict_valley_to_line_of_sight` | `bool` | `false` | `param` | `clean_valley_mask` | Require a straight line from each valley pixel to the channel to stay inside the mask. | A STRONG filter, written for separating one valley floor from a neighbouring one across a broad col, where connectivity alone is not enough. Default FALSE here: on a meandering floodplain a great deal of genuine valley floor has no straight line to the channel that stays inside the mask. On the Royalton 10 m DEM it removed 82% of band 0 on its own, with every other filter disabled - which made the valley mask SMALLER than the floodplain class computed from the same field. Turn it on where two valleys meet around a col and connectivity alone merges them. |
| `restrict_valley_to_channel_connected` | `bool` | `true` | `param` | `clean_valley_mask` | Drop valley patches with no pixel on or beside the channel. | The valley floor of a channel touches that channel. A patch that touches none is flat ground somewhere else. |
| `write_surface_class_raster` | `bool` | `true` | `switch` | - | Write the floodplain / terrace class raster. | Output file: &lt;stem&gt;_surface_class.tif. Only written for the hand field: floodplain and terrace are statements about height above a channel, and mean nothing on a fitted-surface residual. |
| `write_floodplain_mask_raster` | `bool` | `false` | `switch` | - | Write the floodplain class on its own. | Output file: &lt;stem&gt;_floodplain_mask.tif. hand field only. |
| `write_surface_patches_csv` | `bool` | `true` | `switch` | - | Write one row per labelled patch. | Output file: &lt;stem&gt;_surface_patches.csv with label, band, class, source_id, n_pixels, area_m2, median_field_m, aspect_degrees and window_centre_m. The join key for everything downstream. |
| `write_surface_class_csv` | `bool` | `true` | `switch` | - | Write one row per classification record. | Output file: &lt;stem&gt;_surface_class.csv. hand field only. |
| `write_reference_field_profile_csv` | `bool` | `false` | `switch` | - | Write a per-pixel table of the reference field along the channel. | Output file: &lt;stem&gt;_reference_field_profile.csv with x, y, field_m, channel_flow_distance_m, source_id, class, label and band. One row per labelled pixel: the input for a height-against-distance profile figure, and the table that proves the raster labels join to the patch and class tables. Requires the hand field with a computed HAND, since the along-channel distance comes from the pixel's HAND donor node. |
| `write_nearest_source_id_raster` | `bool` | `false` | `switch` | - | Write which channel each pixel's height is measured to. | Output file: &lt;stem&gt;_nearest_source_id.tif (Int32, -1 off-network). The value is the source id of the pixel's HAND DONOR - the channel node it actually drains to - not the nearest channel by straight-line distance, so it partitions the DEM the same way the height reference does. Use it to look at one channel's results in isolation. Needs a computed HAND: a supplied hand_fname carries no donor nodes. |
| `write_ladder_break_metric_csv` | `bool` | `false` | `switch` | - | Write the break DECISION, not just the break. | Output file: &lt;stem&gt;_ladder_break_metric.csv, one row per histogram bin per window: the raw count, the smoothed density the separation ratio is measured on, which detected gap the bin falls in, and whether a candidate or the chosen break lands there. The chosen break on a raw histogram looks arbitrary, and whether it beat several candidates or was the only one on offer cannot be recovered from the answer alone. |
| `write_ladder_window_id_raster` | `bool` | `false` | `switch` | - | Write which along-channel window each pixel's breaks came from. | Output file: &lt;stem&gt;_ladder_window_id.tif (Int16, -1 where the global partition was used or the pixel is not banded). The value is the window's ordinal within its channel, counting downstream - a reach number. Use it to trace a break that looks wrong back to the ground that produced it. |
| `write_ladder_histogram_csv` | `bool` | `false` | `switch` | - | Write the reference-field density each window's breaks were chosen from. | Output file: &lt;stem&gt;_ladder_histogram.csv, one row per histogram bin per window, with the bin containing each break flagged. THIS is the diagnostic that says whether a break separates two populations or merely cuts through one - a deep, wide dip between two modes means two surfaces; a shallow scoop through one broad mode means one surface and an arbitrary cut. Neither the resolvability floor nor the chosen break value can tell you which you have. Plot it with lsdviz3-ladder-density. |
| `ladder_hist_plot_bins` | `int` | `80` | `param` | `write_ladder_histogram_csv` | Bin count for the per-window density in write_ladder_histogram_csv. | Linear in field units, not log10, so it can be read directly against the breaks. Distinct from ladder_hist_bins, which is the log10 histogram the antimode finder searches for dips. |
| `write_ladder_windows_csv` | `bool` | `true` | `switch` | - | Write the per-window resolvability report. | Output file: &lt;stem&gt;_ladder_windows.csv with the level count chosen, the status, the breaks, the band centres and the gap that had to be cleared. READ THIS FIRST: what may be claimed about a ladder depends on whether its windows resolved, and that is not visible in any raster. |
| `write_channel_scope_raster` | `bool` | `false` | `switch` | - | Write the in-scope channel network as a raster. | Output file: &lt;stem&gt;_channel_scope.tif (UInt32, source_id + 1). |
| `write_channel_scope_nodes_csv` | `bool` | `false` | `switch` | - | Write the in-scope channel nodes as a table. | Output file: &lt;stem&gt;_channel_scope_nodes.csv. |
| `write_point_snap_csv` | `bool` | `false` | `switch` | - | Write where each supplied point snapped to. | Outputs: &lt;stem&gt;_point_snaps.csv plus raw and snapped FlatGeobuf point layers, so both can be overlaid in one QC figure. |
| `write_channel_network` | `bool` | `false` | `switch` | - | Write the channel network as vector layers. | Master switch for the three layers below. The network this tool measured its surfaces against is written here so you do not have to chain a second CLI to see it. |
| `write_channel_network_lines` | `bool` | `true` | `switch` | `write_channel_network` | Write channel links as lines. | Lines are the compact default representation of the network used to fit the surfaces; per-node points remain opt-in. Output file: &lt;stem&gt;_channel_network_lines.fgb. |
| `write_channel_network_points` | `bool` | `false` | `switch` | `write_channel_network` | Write channel nodes as points. | Output file: &lt;stem&gt;_channel_network_points.fgb. |
| `write_junction_network_points` | `bool` | `false` | `switch` | `write_channel_network` | Write junctions as points. | Output file: &lt;stem&gt;_junction_network_points.fgb. |

### Floodplain

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `floodplain_threshold_algorithm_choice` | `string` | `preset` | `param` | - | How the floodplain bounds are chosen: preset, fixed, qq_gaussian or cdf_divergence. | fixed uses max_floodplain_hand_metres and max_floodplain_gradient as given. qq_gaussian is Clubb et al. (2017): for each quantity, the lowest value at which the data rejoin a Gaussian drawn through two percentiles of the quantile-quantile plot (floodplain_qq_*), and stay on it. cdf_divergence fits a family of distributions and takes the value at which the empirical CDF most exceeds the best fit. The chosen values are printed and written to the ingested-parameters file, so a run can be repeated with fixed. |
| `max_floodplain_gradient` | `float` | `0.1` | `param` | `floodplain_threshold_algorithm_choice` | Gradient (m/m) above which ground is not floodplain. | A gradient, not an angle (gradient = tan(angle)). The gradient bound of admission floodplain_thresholds, used as given under floodplain_threshold_algorithm_choice fixed. |
| `floodplain_gradient_window_radius_metres` | `float` | `0` | `param` | `admission_algorithm_choice` | Window radius [m] of the polynomial fit the floodplain gradient is read over. 0 = three cells. | LSDTopoTools2's surface_fitting_radius (30 m by default, three cells on 10 m data). |
| `min_terrace_height_above_floodplain_metres` | `float` | `0` | `param` | `level_algorithm_choice` | How far above the floodplain bound a terrace starts [m] (fixed_thresholds). Must be &gt;= 0. | Terraces are ground whose terrace field exceeds max_floodplain_hand_metres plus this, so the two never overlap. LSDTopoTools2's minimum_terrace_height (whose window could overlap the floodplain). Different from min_terrace_hand_metres, the ladder's absolute bound. |
| `max_terrace_height_metres` | `float` | `0` | `param` | `level_algorithm_choice` | Highest a terrace may stand on the terrace field [m] (fixed_thresholds). 0 = no ceiling. | With no ceiling the swath half-width and the gradient bound are the only limits on how high terraces reach. |
| `fill_floodplain_holes` | `bool` | `false` | `param` | `level_algorithm_choice` | Fill enclosed holes in the floodplain (fixed_thresholds). | LSDTopoTools2's fill_floodplain. |
| `max_floodplain_hand_metres` | `float` | `1` | `param` | - | HAND (m) above which ground is not floodplain. | Here an optional absolute override on the lowest band's upper bound. Under band-then-label the ladder already decides where one level ends and the next begins; this is a hard cap for cases where you know the floodplain height independently. Only meaningful for the hand field. Replaces floodplain_max_height_metres. |

### Threshold selection

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `floodplain_qq_lower_percentile` | `float` | `25` | `param` | `floodplain_threshold_algorithm_choice` | Lower percentile of the QQ cloud the Gaussian reference line is drawn through (qq_gaussian). | 25 and 75 (the IQR) are the values Clubb et al. (2017) used. |
| `floodplain_qq_upper_percentile` | `float` | `75` | `param` | `floodplain_threshold_algorithm_choice` | Upper percentile of the QQ cloud the Gaussian reference line is drawn through (qq_gaussian). | Must be above the lower percentile. |
| `floodplain_qq_threshold_fraction` | `float` | `0.01` | `param` | `floodplain_threshold_algorithm_choice` | Fraction of the data range within which a quantile counts as on the Gaussian line (qq_gaussian). | 0.01 (1%) is the paper's value. The threshold is the lowest quantile that is on the line and stays on it for the next 200 quantiles of 10000. Larger values admit more ground. |

### Classification

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `min_terrace_hand_metres` | `float` | `0.0` | `param` | - | Lower HAND bound of the terrace class [m]. | As max_floodplain_hand_metres, at the other end. Only meaningful for the hand field. Replaces terrace_min_height_metres. |
| `channel_adjacency_radius_metres` | `float` | `0.0` | `param` | - | Distance within which a surface counts as channel-adjacent [m]. | 0 = touching only. |
| `valley_corridor_half_width_metres` | `float` | `5000.0` | `param` | - | Valley-corridor half-width [m]. 0 or negative = no limit. | A terrace must have at least one pixel within this Euclidean distance of an in-scope channel. This is what keeps the regional upland out of the terrace class: an interfluve is also flat and also elevated above the trunk channel, but it is not a fluvial terrace. |
| `clip_corridor_per_pixel` | `bool` | `true` | `param` | `valley_corridor_half_width_metres` | Apply the corridor limit per pixel as well as per surface. | Stops a surface that merely touches the valley edge from carrying its class across the whole interfluve. |
| `max_band_spread_metres` | `float` | `50.0` | `param` | - | QC tripwire on the field spread within one surface [m]. 0 disables. | A surface whose robust field spread still exceeds this cannot take the terrace class. Under band-then-label a band's spread is bounded by its own breaks, so this should almost never fire; it is retained as a tripwire, not as a mechanism. |

### Pit filling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `fill_before_routing` | `bool` | `true` | `switch` | - | Choose whether to pit-fill the DEM before building the flow network. | ASKS: should the DEM be filled at all? Filling is OPTIONAL here and this is a genuine analytical choice. True (the default) routes over a pit-free surface. False builds the flow network on the UNFILLED surface, so flow paths terminate at sinks instead of being carried through them. Do not confuse this with raster_is_filled, which is an assertion that the input is already filled. |
| `fill_flat_resolution_choice` | `string` | `inline_gradient` | `param` | - | How flats are made drainable: inline_gradient, barnes2014b_gradient or barnes2014b_float_step. | inline_gradient (default, unchanged): the Priority-Flood fill imposes fill_min_slope_for_routing during the flood. The barnes2014b options make flow converge across flats (Barnes et al. 2014, doi:10.1016/j.cageo.2013.01.009). barnes2014b_gradient uses fill_min_slope_for_routing (m/m, must be &gt; 0) as the minimum gradient; flats are raised about 3-4x more, capped below their rims, and capped flats are reported. barnes2014b_float_step is equivalent to RichDEM ResolveFlats on float32 in its first pass: it ignores fill_min_slope_for_routing and its slopes depend on elevation (of order ULP(z) / cell size). With either barnes2014b option the result is always routable: cells the first pass leaves without a lower neighbour are refilled and resolved again (up to 4 rounds), and any still left are raised by a final priority-flood pass by one float32 step each; the report gives the rounds and the cells changed after the first pass. |

### Flow routing

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `fill_min_slope_for_routing` | `float` | `0.0001` | `param` | `fill_before_routing` | Minimum imposed slope when filling for routing [m/m]. | Used only when fill_before_routing is true. inline_gradient (default): imposed during the Barnes et al. 2014 Priority-Flood. barnes2014b_gradient: the minimum gradient across flats. barnes2014b_float_step: ignored. |
| `boundary_conditions` | `string` | `oooo` | `param` | - | Edge boundary codes, north east south west. | o = open, n = no-flux, p = periodic, b = forced outlet. |

### Channel network

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `threshold_contributing_pixels` | `int` | `1000` | `param` | - | Contributing-pixel threshold defining a channel head. | The total number of upslope contributing pixels from the flow-routing object, not the number of direct donors to a pixel. Reduce this value to extract a denser network. This is a PIXEL COUNT and therefore depends on the grid resolution: the equivalent drainage area is threshold_contributing_pixels x cell area. Where a genuine area is wanted instead, the convention is &lt;what&gt;_area_threshold_m2. |

### Channel scope

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `break_flats_with_channels` | `bool` | `false` | `switch` | - | Remove in-scope channel pixels from the admitted set before banding. | Default FALSE, changed from the old flat-surface-finder behaviour. It used to be needed to stop one surface spanning a channel; banding by height does that on its own, since the two banks of an incised channel are not in the same band as the channel. All it does now is keep the channel itself out of a surface, which matters for area statistics and little else. |
| `connect_surfaces_to_channels` | `bool` | `false` | `switch` | - | Record which in-scope channel each surface is linked to. | Adds the connection diagnostics used to attribute a surface to a reach. The patch table already carries source_id from the height reference's donor, so this is for the cases where the linkage itself is the question. |
| `drainage_injection_points_fname` | `string` | - | `io` | - | Vector of points injecting extra drainage area, in read_path. | For a channel entering the DEM from outside it, where the contributing-pixel threshold alone would not start a channel at the boundary. Empty = no injection. |
| `source_points_fname` | `string` | - | `io` | - | Channel source points (CSV or vector), in read_path. | Empty = every channel above the threshold is in scope. Supplying sources restricts the analysis to those channels, which is how you measure terraces along one trunk rather than every rill in the DEM. |
| `source_points_input_format` | `string` | `auto_detect` | `io` | `source_points_fname` | Source point format: auto_detect, csv or ogr_vector_layer. | auto_detect infers from the file extension. |
| `source_points_layer_name` | `string` | - | `io` | `source_points_fname` | Layer name within a source point vector file. | Empty = the first layer. |
| `source_points_latitude_field_name` | `string` | `latitude` | `io` | `source_points_fname` | Latitude column name in a source point CSV. | CSV points are read as EPSG:4326 and transformed to the DEM CRS. |
| `source_points_longitude_field_name` | `string` | `longitude` | `io` | `source_points_fname` | Longitude column name in a source point CSV. | See source_points_latitude_field_name. |
| `outlet_points_fname` | `string` | - | `io` | - | Channel outlet points (CSV or vector), in read_path. | Truncates each traced channel at its matched outlet, so the analysis covers a reach rather than everything down to the DEM edge. |
| `outlet_points_input_format` | `string` | `auto_detect` | `io` | `outlet_points_fname` | Outlet point format: auto_detect, csv or ogr_vector_layer. | As source_points_input_format. |
| `outlet_points_layer_name` | `string` | - | `io` | `outlet_points_fname` | Layer name within an outlet point vector file. | Empty = the first layer. |
| `outlet_points_latitude_field_name` | `string` | `latitude` | `io` | `outlet_points_fname` | Latitude column name in an outlet point CSV. | See source_points_latitude_field_name. |
| `outlet_points_longitude_field_name` | `string` | `longitude` | `io` | `outlet_points_fname` | Longitude column name in an outlet point CSV. | See source_points_latitude_field_name. |
| `auto_detect_outlet_points` | `bool` | `true` | `switch` | `source_points_fname` | Infer source and outlet roles from a single combined point file. | A supplied point that snaps strictly downstream of another is treated as that one's outlet. Ignored when outlet_points_fname is given, because then the roles are stated directly. |
| `snap_sources_to_channel_threshold` | `bool` | `true` | `switch` | `source_points_fname` | Snap supplied source points onto the threshold channel network. | A hand-digitised point rarely lands exactly on the routed channel. |
| `source_snap_max_search_radius_pixels` | `int` | `20` | `param` | `snap_sources_to_channel_threshold` | Search radius when snapping a source point [pixels]. | A point that finds no threshold channel within this radius is reported as unsnapped in the point snap diagnostics. |
| `outlet_snap_max_distance_metres` | `float` | `500.0` | `param` | `outlet_points_fname` | Maximum distance an outlet may be moved onto a traced channel [m]. | Outlets snap onto the traced flow paths, not onto the channel network, so a slightly misplaced outlet still truncates the intended channel. |

### Valley floor

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `force_valley_hand_metres` | `float` | `0.0` | `param` | - | Define the valley floor directly as height above channel below this [m], instead of using the ladder's break. 0 = use the ladder. | The ladder finds a break from the data; sometimes you already know the answer - from field mapping, a flood stage, a published terrace height - and want it applied rather than rediscovered. Combined with force_valley_max_distance_metres this defines the valley floor outright: HAND below the threshold AND within the distance. The band index, labels and patch table are UNAFFECTED and still come from the ladder; only &lt;stem&gt;_valley_mask.tif is forced, so the two can be compared. |
| `force_valley_max_distance_metres` | `float` | `0.0` | `param` | - | Euclidean distance to an in-scope channel bounding the forced valley floor [m]. 0 = no distance limit. | The companion to force_valley_hand_metres, and the reason the pair is useful: a HAND threshold alone admits every low patch in the DEM, including ground belonging to other drainages entirely. The distance limit is what confines it to this channel's valley. Usable on its own as well, if the ladder is finding the break but reaching too far. |

### Interior nodata

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `interior_nodata_method` | `string` | `sink` | `param` | - | How nodata pixels inside the DEM are treated: sink or fill. | Interior nodata are nodata pixels not connected to the raster edge through other nodata (isolated voids, moulins, masked lakes). sink (the default, unchanged behaviour): they are outlets, so flow that reaches them ends there, as at the DEM edge; use this for real internal drainage (e.g. moulins on an ice sheet). fill: they are interpolated with nodata_fill_method_choice before any filling or routing, so flow passes across them; use this for data voids. Nodata connected to the edge (e.g. sea removed by min_elevation) is never filled and always acts as an open boundary. |
| `nodata_fill_method_choice` | `string` | `laplacian` | `param` | `interior_nodata_method=fill` | Method for filling interior nodata pixels. | laplacian: smooth Dirichlet interpolation (recommended). iterative_weighted_average: fast, well-suited to thin voids. Used when interior_nodata_method is fill. |
| `nodata_fill_max_iter` | `int` | `10000` | `param` | `interior_nodata_method=fill` | Maximum iterations for iterative nodata fill methods. | Used when interior_nodata_method is fill. |

### Orientations

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `bearings_relative_to_true_north` | `bool` | `false` | `param` | - | Write bearings, aspects and orientations clockwise from true north. | Off (default): they are clockwise from grid north (grid up), so in an oblique frame with the feature along +x, 90 degrees means along the feature. On: each value is converted with the grid convergence at its own location (each cell of a raster, each row of a table), which also removes the up to about 3 degree difference between UTM grid north and true north. Angles between two directions (junction angles, deflections, contrasts) are unchanged, and hillshade illumination always stays relative to grid up. Needs a raster with a CRS. |

<!-- END GENERATED CLI PARAMETERS -->
