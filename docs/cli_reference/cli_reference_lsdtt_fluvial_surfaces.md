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

**`min_component_pixels` is per band here**, unlike in
`lsdtt-flat-surface-finder` where the same name filters components of the whole
mask. Here it is a floor on the size of one surface on one level.

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
| `fallback_admission_algorithm_choice` | string | `relief` | Built-in admission when `flat_mask_fname` is empty: `relief` or `slope`. |
| `fallback_relief_window_radius_pixels` | int | `5` | Moving-window **radius** in pixels for the fallback relief admission. Not a window side: `5` gives an 11×11 window. `lsdtt-flat-surface-finder` and `lsdtt-mountain-front` spell their equivalent `relief_window_side_pixels`, which is a full **side length**, so the same number means a different window there. The old name `fallback_relief_window_pixels` is accepted as a deprecated alias. |
| `fallback_relief_threshold_metres` | float | `1.0` | Local relief below which a cell is admitted as flat [m]. |
| `fallback_slope_threshold` | float | `0.05` | Gradient magnitude below which a cell is admitted as flat [m/m]. |
| `min_component_pixels` | int | `50` | Minimum patch size in pixels, **applied per band**. |

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
| `vertical_accuracy_mode` | string | `fixed` | How `E_v` is arrived at: `fixed` or `slope_dependent`. |
| `vertical_accuracy_gentle_metres` | float | `0.0` | `E_v` [m] for windows at or below the break gradient. `slope_dependent` only. |
| `vertical_accuracy_steep_metres` | float | `0.0` | `E_v` [m] for windows above the break gradient. `slope_dependent` only. |
| `vertical_accuracy_break_degrees` | float | `11.31` | Gradient [degrees] separating gentle from steep `E_v`. |
| `vertical_accuracy_slope_window_metres` | float | `0` | Window diameter [m] of the polynomial fit the window gradient is measured over. `0` = five cells. |
| `ladder_reach_window_metres` | float | `0.0` | Along-channel window length for reach-relative breaks [m]. `0` = auto, negative = one global ladder. Ignored for `fitted_surface`. |
| `ladder_min_window_pixels` | int | `1000` | Minimum admitted pixels for a window to find its own breaks. |

### Classification

| Parameter | Type | Default | Description |
|---|---|---|---|
| `floodplain_max_height_metres` | float | `1.0` | Upper HAND bound of the floodplain class [m]. An optional hard cap; the ladder already decides where a level ends. |
| `terrace_min_height_metres` | float | `0.0` | Lower HAND bound of the terrace class [m]. |
| `channel_adjacency_radius_metres` | float | `0.0` | Distance within which a surface counts as channel-adjacent [m]. |
| `valley_corridor_half_width_metres` | float | `5000.0` | Valley-corridor half-width [m]. `0` or negative = no limit. |
| `clip_corridor_per_pixel` | bool | `true` | Apply the corridor limit per pixel as well as per surface. |
| `max_band_spread_metres` | float | `50.0` | QC tripwire on the field spread within one surface [m]. `0` disables. Under band-then-label it should almost never fire. |

### Flow routing and channel scope

| Parameter | Type | Default | Description |
|---|---|---|---|
| `fill_dem_for_flow_routing` | bool | `true` | Fill pits before building flow routing. |
| `fill_min_slope_for_routing` | float | `0.0001` | Minimum imposed slope when filling [m/m]. |
| `boundary_conditions` | string | `oooo` | Edge boundary codes, north east south west. |
| `channel_threshold_contributing_pixels` | int | `1000` | Contributing-pixel threshold defining a channel. |
| `break_flats_with_channels` | bool | `false` | Remove in-scope channel pixels from the admitted set before banding. |
| `connect_surfaces_to_channels` | bool | `false` | Record which in-scope channel each surface is linked to. |
| `auto_detect_outlet_points` | bool | `true` | Infer source and outlet roles from a single combined point file. |
| `snap_sources_to_channel_threshold` | bool | `true` | Snap supplied source points onto the threshold channel network. |
| `source_snap_max_radius_pixels` | int | `20` | Search radius when snapping a source point [pixels]. |
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
| `patch_aspect_window_metres` | float | `0` | Window diameter [m] of the polynomial fit the aspect is measured over. `0` = five cells. |
| `write_channel_scope_raster` | bool | `false` | The in-scope channel network as a raster. |
| `write_channel_scope_nodes_csv` | bool | `false` | The in-scope channel nodes as a table. |
| `write_point_snap_csv` | bool | `false` | Where each supplied point snapped to. |
| `write_channel_network` | bool | `false` | Write the channel network as vector layers. |
| `write_channel_network_lines` | bool | `true` | Channel links as lines. |
| `write_channel_network_points` | bool | `false` | Channel nodes as points. |
| `write_junction_network_points` | bool | `false` | Junctions as points. |

---

## See also

- [Floodplains, terraces and valleys](../fluvial_surfaces_and_valleys.md) — the
  band-then-label ordering, the three break finders, the resolvability floor,
  and where the methods come from
- [CLI reference: lsdtt-flat-surface-finder](cli_reference_lsdtt_flat_surface_finder.md)
- [CLI reference: lsdtt-valley-centreline](cli_reference_lsdtt_valley_centreline.md)
