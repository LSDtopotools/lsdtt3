# CLI reference: lsdtt-valley-centreline

**Trace a centreline through a supplied valley, and measure it.** Takes a valley
mask, carves a synthetic trough through it, routes the centreline down that
trough, fits a pinned monotonic thalweg baseline, and casts cross sections
orthogonal to the centreline to measure valley width and cross-section
hydraulics.

**It does not identify the valley.** That moved to
[`lsdtt-fluvial-surfaces`](cli_reference_lsdtt_fluvial_surfaces.md), where a
valley is band 0 of a reference-field ladder — a formulation that also works on a
fan, which a corridor HAND break does not. `valley_mask_fname` is **required**;
running without it is a hard error naming that tool. For the algorithms and their
provenance in Clubb et al. (2022), see
[Floodplains, terraces and valleys](../design_and_algorithms/fluvial_surfaces_and_valleys.md).

## Usage

```bash
lsdtt-valley-centreline <path> <parameter_file>
lsdtt-valley-centreline <parameter_file>      # path defaults to ./
lsdtt-valley-centreline -h                    # parameter help, then exit
lsdtt-valley-centreline -h <parameter_name>   # full guidance on one parameter
lsdtt-valley-centreline -h html               # write an HTML/CSV reference
lsdtt-valley-centreline --help-json           # machine-readable parameter list
lsdtt-valley-centreline -v                    # version
lsdtt-valley-centreline -c                    # citation
```

Parameter file lines are `key: value`; `#` starts a comment.

---

## Pipeline

| Step | What happens |
|---|---|
| 1 | Load DEM, pit-fill for routing (`fill_min_slope_for_routing`) |
| 2 | Build the D8 flow object (`boundary_conditions`) |
| 3 | Load source and outlet points, map them to pixels |
| 4 | Resolve channel scope: snapping, optional outlet truncation, optional auto role detection |
| 5 | Load `valley_mask_fname`; each channel's valley is that mask intersected with its corridor |
| 6 | Carve the synthetic trough and route the centreline down it |
| 7 | Fit the pinned monotonic thalweg baseline |
| 8 | Optionally refine the supplied mask against that baseline |
| 9 | Cast cross sections, measure width, detect features, derive hydraulics |
| 10 | Write rasters, tables and vectors |

**One centreline per resolved source.** There is no threshold-only all-channel
mode; `source_points_fname` is required.

---

## File I/O and run control

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | `./` | Directory from which input files are read. |
| `read_fname` | string | | Input DEM, with or without extension. |
| `write_path` | string | same as `read_path` | Output directory. |
| `write_fname` | string | stem of `read_fname` | Output prefix for all products. |
| `valley_mask_fname` | string | | **Required.** Valley mask raster, in `read_path`. Cells > 0.5 are valley. Normally `<stem>_valley_mask.tif` from `lsdtt-fluvial-surfaces`. |
| `n_threads` | int | `0` | OpenMP threads. `0` = all available. |
| `write_output_format` | string | `cog` | Raster format: `cog`, `gtiff`, `same_as_input`. Existing parameter files that explicitly select `gtiff` continue to write plain GeoTIFFs. |

---

## Source and outlet points

Supply either a source-only file, a source file plus an explicit outlet file, or
a single combined file with `auto_detect_outlet_points: true` (the default),
which classifies the downstream points as outlets. Point files may be CSV
(`latitude` / `longitude` in EPSG:4326) or any OGR vector (`.gpkg`, `.fgb`, …).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `source_points_fname` | string | | **Required.** Vector or CSV file with source points, in `read_path` (an absolute path is used as given). |
| `source_points_input_format` | string | `auto_detect` | `auto_detect`, `csv_latitude_longitude_epsg4326`, `ogr_vector_layer`. |
| `source_points_layer_name` | string | | Optional layer name for vector input. |
| `source_points_latitude_field_name` | string | `latitude` | CSV latitude column. |
| `source_points_longitude_field_name` | string | `longitude` | CSV longitude column. |
| `outlet_points_fname` | string | | Optional outlet points file, to truncate channels. |
| `outlet_points_input_format` | string | `auto_detect` | As the source format. |
| `outlet_points_layer_name` | string | | Optional layer name. |
| `outlet_points_latitude_field_name` | string | `latitude` | CSV latitude column. |
| `outlet_points_longitude_field_name` | string | `longitude` | CSV longitude column. |

> **Renamed 2026-10-26.** These were `source_points_path` /
> `outlet_points_path`, taken relative to the working directory. They are now
> `source_points_fname` / `outlet_points_fname`, resolved against `read_path`
> as in `lsdtt-fluvial-surfaces`. The old keys still run, with a warning; an
> absolute path behaves as before.

## Channel scope

| Parameter | Type | Default | Description |
|---|---|---|---|
| `threshold_contributing_pixels` | int | `1000` | Contributing-pixel threshold defining channel cells for source snapping. |
| `auto_detect_outlet_points` | bool | `true` | Classify points in `source_points_fname` into sources and outlets. Disabled when an explicit outlet file is given. |
| `snap_sources_to_channel_threshold` | bool | `true` | Snap source points to the nearest cell meeting the threshold. |
| `source_snap_max_search_radius_pixels` | int | `25` | Maximum expanding-square search radius for snapping. |
| `outlet_snap_max_distance_metres` | float | `0` | Maximum distance between a raw outlet point and its traced-path node. `0` = unlimited. |
| `boundary_conditions` | string | `oooo` | Flow boundary conditions (N, E, S, W). |
| `fill_min_slope_for_routing` | float | `0.0001` | Minimum slope imposed when pit-filling for routing. |

---

## Valley refinement

| Parameter | Type | Default | Description |
|---|---|---|---|
| `refine_valley_with_baseline` | bool | `true` | Re-test the supplied mask against the centreline baseline after routing, removing the stale-reference wedges a nearest-channel-pixel test grows onto downstream-sloping terraces where the channel veers away. |
| `refine_max_height_metres` | float | `0.0` | Height above the pinned baseline below which a pixel is valley [m]. **Required by the refinement, which does nothing while this is 0.** |

The threshold must be stated rather than inherited: this tool no longer derives a
HAND break of its own, and a refinement run against a fabricated threshold would
look like a result and be one only by accident. The refined mask is written for
you and is **never fed back** to `lsdtt-fluvial-surfaces`.

---

## Trough carving and routing

The valley is punched out of the DEM, each valley pixel is set to the minimum
bank elevation nearby and then depressed in proportion to its distance from the
wall, and the result is breached, filled and re-carved. The D8 path down the
resulting V-shaped trough is the centreline. This follows Clubb et al. (2022).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `trough_scaling_factor` | float | `0.1` | Trough depth per metre of distance to the valley wall (α in Clubb et al. 2022, Eq. 2). Doubled automatically on routing failure, up to 2 retries. |
| `trough_along_channel_slope` | float | `0.02` | Synthetic downstream tilt of the trough (m per m along channel). Stops a near-flat or through valley routing its centreline out of the wrong tile edge. |
| `centreline_loops` | int | `5` | Trough carve / breach / fill iterations. |
| `min_bank_elevation_window_radius_metres` | float | `200` | Radius of the local-minimum filter on the nearest-wall elevation field. |
| `trough_nodata_buffer_pixels` | int | `2` | The valley hole is dilated by this many pixels before bank elevations are sampled. |
| `trough_fill_min_slope` | float | `0.0001` | Minimum slope for pit filling of the trough surface. |
| `min_centreline_nodes` | int | `10` | A routed centreline shorter than this counts as a failure. |

---

## Cross sections

Stations are placed along the centreline and a transect is cast across it,
orthogonal to a direction smoothed over `cross_section_tangent_window_nodes`
nodes each side. Valley width at a station is the span of the contiguous
in-valley run through the centreline; a side not found before the transect runs
out marks the station not-found and excludes it from the per-channel statistic.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_valley_cross_sections` | bool | `true` | Write regularly spaced cross sections orthogonal to each centreline. |
| `cross_section_spacing_metres` | float | `500` | Along-centreline spacing between stations (m). |
| `cross_section_tangent_window_nodes` | int | `8` | Nodes each side used to smooth the local direction the transect is drawn orthogonal to. |
| `cross_section_length_multiplier` | float | `1.5` | Transect total length as a multiple of the widest valley width for the channel. |
| `cross_section_min_half_length_metres` | float | `0` | Minimum transect half-length (m), so a narrow channel's sections still show wall and upland context. |
| `cross_section_max_half_length_metres` | float | `5000` | Hard cap (m) on the outward march when finding a valley side. |
| `clip_overlapping_sections` | bool | `true` | Clip each transect side where the ray becomes closer to a neighbouring station centre than to its own (local Voronoi bisector clipping). |
| `neighbour_clip_window` | int | `2` | Neighbouring stations each side considered for overlap clipping. |
| `min_section_reach_fraction` | float | `0.35` | Drop a station whose clipped transect falls below this fraction of the nominal length. |

### Curved cross sections

Straight rays are correct in straight reaches and cross one another in tight
bends. Curved sections are traced as isocontours of a smooth along-valley
coordinate field — straight where the valley is straight, curving gradually
through bends like the grid lines of a boundary-fitted mesh.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `curved_cross_sections` | bool | `false` | Trace sections as streamwise isocontours instead of straight rays. |
| `curved_wall_clearance_fraction` | float | `0.10` | A shortened side never clips below (1 + this) × the distance to the valley wall. |
| `curved_streamwise_smoothing_passes` | int | `40` | Anchored Laplacian smoothing passes on the streamwise coordinate field. Higher = smoother sections, fewer creases at bends. |
| `curved_coarse_spacing_metres` | float | `700` | Sections are placed coarse-to-fine: a sparse backbone at this spacing, then midpoints down to `cross_section_spacing_metres`. Keeps coverage even instead of clustering next to gaps. |

### Features and hydraulics

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_cross_section_features` | bool | `true` | Detect and write, along each section: D8 channel crossings, left/right valley walls, hydraulic-step boundaries, maximum-curvature peaks and multidirectional-hillshade peaks. |
| `cross_section_curvature_radius_metres` | float | `30` | Polynomial-fit window radius (m) for the maximum-curvature raster used to find slope-break peaks. |
| `cross_section_feature_prominence_fraction` | float | `0.15` | Minimum peak prominence, as a fraction of that profile's range over the transect. |
| `write_cross_section_hydraulics` | bool | `false` | Flood each profile stage by stage and record flow width, wetted area, wetted perimeter and hydraulic depth (A/W). Breakpoints in hydraulic depth against elevation above thalweg mark bankfull and inset flat-surface levels (Delchiaro et al. 2025, *Computers & Geosciences*, doi:10.1016/j.cageo.2025.106001). |
| `write_cross_section_informed_valley_mask_raster` | bool | `false` | Build a valley-floor mask from the sections' hydraulic breaks: each corridor pixel's height above the nearest node's thalweg baseline is compared with the channel's **median** lowest hydraulic break, so one anomalous section cannot streak onto the plateau. Requires `write_cross_section_features`. |

---

## Other outputs

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_valley_mask_raster` | bool | `true` | Write the valley mask actually used (0 outside; `source_id + 1` in each channel's valley). |
| `write_valley_centreline` | bool | `true` | Write the centrelines as vector lines. |
| `write_valley_centreline_nodes` | bool | `true` | Write the per-node centreline table. |
| `write_hillshade_raster` | bool | `false` | Standard single-azimuth hillshade of the DEM — a terrain backdrop for the section figures. |
| `write_multidirectional_hillshade_raster` | bool | `false` | Mean of hillshades from several evenly spaced azimuths. Reduces directional bias and picks out valley and terrace margins. |
| `multidirectional_hillshade_n_directions` | int | `8` | Number of azimuths for the multidirectional hillshade. |

### Files written

| File | Contents |
|---|---|
| `<write_fname>_valley_centreline.fgb` | Centrelines as polylines |
| `<write_fname>_valley_centreline_nodes.csv` | Per-node table: position, flow distance, baseline elevation |
| `<write_fname>_valley_cross_sections.*` | Section geometries |
| `<write_fname>_valley_cross_section_samples.csv` | Per-sample elevation along each section |
| `<write_fname>_valley_cross_section_features.*` | Detected walls, channel crossings, curvature and hillshade peaks |
| `<write_fname>_valley_cross_section_sides.*`, `_boundaries.*` | Per-side extents and clipping outcome |
| `<write_fname>_valley_cross_section_hydraulics.csv` | Stage / width / area / wetted perimeter / hydraulic depth |
| `<write_fname>_valley_cross_section_hydraulic_breaks.csv` | Breakpoints in the hydraulic-depth curve |
| `<write_fname>_valley_mask.tif` | The mask used, per channel |
| `<write_fname>_cross_section_informed_valley_mask.tif` | Hydraulic-break valley mask, when enabled |
| `<write_fname>_hillshade.tif`, `_multidirectional_hillshade.tif` | Backdrops |
| `<write_fname>_README.csv` / `.html` | The resolved parameters for that run |

---

## Retired parameters

Valley identification moved out. These keys are **accepted for one release and
ignored**, with a message naming where the equivalent now lives; they do not stop
the run, because a missing `valley_mask_fname` stops it first with a message that
says the same thing once.

| Ignored key | Now |
|---|---|
| `valley_fallback_hand_metres` | **no successor.** There is no fallback break: where no level count is admissible the ladder in `lsdtt-fluvial-surfaces` reports `unresolved` rather than substituting a number |
| `valley_force_hand_metres` | `force_valley_hand_metres` in `lsdtt-fluvial-surfaces` |
| `valley_min_hand_metres`, `valley_max_hand_metres` | `ladder_min_break_metres`, `ladder_max_break_metres` there |
| `valley_wall_peak_fraction`, `valley_break_fraction`, `valley_break_fraction_min`, `valley_break_fraction_max`, `valley_break_method`, `valley_break_relief_window_metres` | the ladder's break-finder parameters there |
| `valley_break_reach_window_metres` | `ladder_reach_window_metres` there |
| `valley_corridor_width_metres` | renamed to `valley_corridor_half_width_metres`, then moved there |
| `valley_dilation_pixels`, `valley_max_detour_ratio` | same names, in `lsdtt-fluvial-surfaces` — mask cleaning happens where the mask is made |
| `write_valley_hand_histogram`, `write_valley_break_metric`, `write_valley_reach_id_raster` | the ladder diagnostics there |

---

## See also

- [Floodplains, terraces and valleys](../design_and_algorithms/fluvial_surfaces_and_valleys.md) — the
  trough construction, its origin in Clubb et al. (2022), and where lsdtt3
  departs from it
- [CLI reference: lsdtt-fluvial-surfaces](cli_reference_lsdtt_fluvial_surfaces.md)
- [CLI reference: lsdtt-flat-surface-finder](cli_reference_lsdtt_flat_surface_finder.md)
- [CLI reference: lsdtt-swath-profile](cli_reference_lsdtt_swath_profile.md)

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

### Channel network

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `threshold_contributing_pixels` | `int` | `1000` | `param` | - | Contributing-pixel threshold defining a channel head. | The total number of upslope contributing pixels from the flow-routing object, not the number of direct donors to a pixel. Reduce this value to extract a denser network. This is a PIXEL COUNT and therefore depends on the grid resolution: the equivalent drainage area is threshold_contributing_pixels x cell area. Where a genuine area is wanted instead, the convention is &lt;what&gt;_area_threshold_m2. |

### General

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `source_points_fname` | `string` | - | - | - | Vector or CSV file with source points, in read_path (an absolute path is used as given). | REQUIRED. Each source defines one channel (traced downstream); one valley + centreline is produced per channel. A combined file may also contain outlets (see auto_detect_outlet_points). |
| `source_points_input_format` | `string` | `auto_detect` | - | - | Point input format. | auto_detect, csv_latitude_longitude_epsg4326, or ogr_vector_layer. |
| `source_points_layer_name` | `string` | - | - | - | Optional layer name for source-point vector input. | - |
| `source_points_latitude_field_name` | `string` | `latitude` | - | - | Latitude column for CSV source points. | - |
| `source_points_longitude_field_name` | `string` | `longitude` | - | - | Longitude column for CSV source points. | - |
| `outlet_points_fname` | `string` | - | - | - | Vector or CSV file with outlet points, in read_path (an absolute path is used as given). | Optional. Bounds each traced channel at the matched outlet. |
| `outlet_points_input_format` | `string` | `auto_detect` | - | - | Point input format. | auto_detect, csv_latitude_longitude_epsg4326, or ogr_vector_layer. |
| `outlet_points_layer_name` | `string` | - | - | - | Optional layer name for outlet-point vector input. | - |
| `outlet_points_latitude_field_name` | `string` | `latitude` | - | - | Latitude column for CSV outlet points. | - |
| `outlet_points_longitude_field_name` | `string` | `longitude` | - | - | Longitude column for CSV outlet points. | - |

### Channel scope

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `auto_detect_outlet_points` | `bool` | `true` | `param` | - | Auto-classify points in source_points_fname into sources and outlets. | A point downstream of another supplied point becomes an outlet. Disabled when an explicit outlet_points_fname is given. |

### Channel scope (mirrors lsdtt-flat-surface-finder)

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `snap_sources_to_channel_threshold` | `bool` | `true` | `param` | - | Snap source points to the nearest cell meeting the channel threshold. | - |
| `source_snap_max_search_radius_pixels` | `int` | `25` | `param` | - | Maximum expanding-square search radius for source snapping. | - |
| `outlet_snap_max_distance_metres` | `float` | `0` | `param` | - | Maximum distance between a raw outlet point and its traced-path node. | 0 = unlimited. |
| `boundary_conditions` | `string` | `oooo` | `param` | - | Flow boundary conditions (N,E,S,W). | - |
| `fill_min_slope_for_routing` | `float` | `0.0001` | `param` | - | Minimum slope (m/m) used when pit-filling the DEM for flow routing. | inline_gradient (default): imposed during the Priority-Flood fill. barnes2014b_gradient: the minimum gradient across flats. barnes2014b_float_step: ignored. Separate from trough_fill_min_slope. |

### Valley width network

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `run_valley_width_network` | `bool` | `false` | `switch` | - | Measure valley widths on every channel from the channel threshold instead of on supplied source points. | The network mode of LSDTopoTools2's lsdtt-valley-metrics (extract_all_channels). Channel heads are traced downstream longest-first and each stops where it meets an earlier trace, so the trunk of each basin is its longest flow path and every tributary segment ends at its parent's channel: each valley is measured once. source_points_fname must be empty. Writes _valley_segments.csv, _valley_segment_channels.fgb (each segment's D8 channel), _valley_widths.csv (one row per station, with drainage area at the station's channel pixel, matched in downstream order so it never decreases), and adds segment_id, parent_segment_id, nesting_depth and excluded_reason to the cross-section outputs. |
| `network_centreline_algorithm_choice` | `string` | `medial_path` | `param` | `run_valley_width_network` | Line each segment's sections hang from: medial_path or smoothed_channel. | Not the raw D8 channel: sections cast from a meandering channel swing with every bend and cross one another. medial_path: the least-cost walk along the valley's medial axis (cost = largest distance-to-wall minus distance-to-wall), head to confluence, confined to the valley mask plus the segment's own channel pixels. smoothed_channel: the D8 channel smoothed by a Gaussian-smoothed cubic spline whose sigma scales with the valley width; it may not leave the valley mask, so where it would, it runs along the mask edge until it comes back in. |
| `network_centreline_window_margin_multiple` | `float` | `3` | `param` | `network_centreline_algorithm_choice` | medial_path: each walk runs in the segment's bounding box padded by this multiple of its largest distance-to-wall (plus 10 cells). | Bounds the work per segment. The medial axis lies within a valley width of the channel, so 3 leaves room on both sides. |
| `network_centreline_smoothing_width_multiple` | `float` | `1` | `param` | `network_centreline_algorithm_choice` | smoothed_channel: Gaussian sigma (along-channel metres) as a multiple of the segment's valley width. | Valley width = 2 x the 90th percentile of distance-to-wall over the segment's channel pixels; sigma is never below two cells. Larger values remove longer meanders but cut more corners. |
| `run_tributary_stop_at_parent` | `bool` | `true` | `switch` | `run_valley_width_network` | Stop each tributary's sections where they reach the parent valley. | Going downstream, a tributary's stations stop at the first one whose centre is within tributary_stop_width_fraction x the parent's local half-width of the parent's centreline. The local half-width is judged per tributary station: the median of the parent's kept sections' offset on the side that station lies, over the parent stations within tributary_parent_width_window_multiple x the tributary's valley width of the one nearest it. A station whose width line crosses the parent's centreline also stops it. Stopped stations are kept with excluded_reason near_parent_channel. |
| `tributary_stop_width_fraction` | `float` | `1` | `param` | `run_tributary_stop_at_parent` | Fraction of the parent's local half-width at which a tributary's sections stop. | 1 = at the parent's valley wall. |
| `tributary_parent_width_window_multiple` | `float` | `2` | `param` | `run_tributary_stop_at_parent` | Half-length of the along-parent window whose median offset is the parent's local half-width, as a multiple of the tributary's valley width. | The parent sections at the tributary's mouth are widened by it; the window is kept longer than the mouth (and never fewer than the nearest parent station and its two neighbours) so they stay a minority of the median at any station spacing. Valley width = 2 x the 90th percentile of distance-to-wall over the tributary's channel pixels. |
| `run_tributary_crossing_removal` | `bool` | `true` | `switch` | `run_valley_width_network` | Exclude parent sections that cross a tributary's kept sections. | A parent section whose width line (left side to right side) crosses the width line of a kept section of one of its direct tributaries has run up that tributary's valley, which makes it anomalously wide. Applied to every parent against its own tributaries, so it nests down the tree. Excluded sections are kept with excluded_reason tributary_crossing. |
| `channel_definition_algorithm_choice` | `string` | `contributing_pixels` | `param` | `run_valley_width_network` | What counts as a channel: contributing_pixels or stream_order. | contributing_pixels: the whole network down from heads at threshold_contributing_pixels. stream_order: only its pixels of Strahler order &gt;= min_channel_stream_order (Clubb et al. 2017 used third order). Only used with run_valley_width_network. |
| `min_channel_stream_order` | `int` | `3` | `param` | `channel_definition_algorithm_choice` | Lowest Strahler order counted as channel (stream_order only). | LSDTopoTools2 called this threshold_SO. |

### Pit filling

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `fill_flat_resolution_choice` | `string` | `inline_gradient` | `param` | - | How flats are made drainable: inline_gradient, barnes2014b_gradient or barnes2014b_float_step. | inline_gradient (default, unchanged): the Priority-Flood fill imposes fill_min_slope_for_routing during the flood. The barnes2014b options make flow converge across flats (Barnes et al. 2014, doi:10.1016/j.cageo.2013.01.009). barnes2014b_gradient uses fill_min_slope_for_routing (m/m, must be &gt; 0) as the minimum gradient; flats are raised about 3-4x more, capped below their rims, and capped flats are reported. barnes2014b_float_step is equivalent to RichDEM ResolveFlats on float32 in its first pass: it ignores fill_min_slope_for_routing and its slopes depend on elevation (of order ULP(z) / cell size). With either barnes2014b option the result is always routable: cells the first pass leaves without a lower neighbour are refilled and resolved again (up to 4 rounds), and any still left are raised by a final priority-flood pass by one float32 step each; the report gives the rounds and the cells changed after the first pass. |

### Valley mask

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `valley_mask_fname` | `string` | - | `io` | - | Valley mask raster from lsdtt-fluvial-surfaces, in read_path. | Cells &gt; 0.5 are valley. Normally &lt;stem&gt;_valley_mask.tif from lsdtt-fluvial-surfaces, which is band 0 of its reference-field ladder. When supplied, this tool does NOT identify the valley: it skips the corridor HAND break search entirely and goes straight to the medial-axis centreline, the trough carve and the pinned baseline, which are the parts that are genuinely its own. Every valley_* break parameter is then ignored, and each channel reports no break rather than a fabricated one. Leave empty to derive the valley here as before (deprecated: that path moves out entirely at the next step). |
| `valley_mask_closing_radius_metres` | `float` | `0.0` | `param` | - | Seal narrow gaps in the valley mask with a morphological closing [m]. 0 = off. | A dilation followed by an erosion of the same radius. It removes thin non-valley features that reach into the valley - a road, a riser, a levee, scroll bars, a field boundary a little above the floor - which hole filling cannot touch, because a finger open at one end is not an enclosed hole. Because the erosion undoes the dilation, the outer extent stays about where it was. Set it to roughly half the width of the widest intrusion to close; too large and it bridges a genuine terrace riser. Applied before the opening and the hole filling. |
| `valley_mask_opening_radius_metres` | `float` | `0.0` | `param` | - | Remove narrow tendrils of the valley mask with a morphological opening [m]. 0 = off. | An erosion followed by a dilation of the same radius: valley that is narrower than about twice the radius - a thread running up a gully or along a ditch - disappears and the rest is restored to its extent. The converse of the closing, which removes thin NON-valley features. Applied after the closing and before the hole filling. A reach of real valley narrower than twice the radius is removed too. |
| `fill_valley_mask_holes` | `bool` | `false` | `switch` | - | Fill enclosed holes in the valley mask. | A background region not connected to the outside of the valid DEM is a hole. Applied last, so holes the closing has sealed off are filled too. |

### Relative elevation

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_relative_elevation_raster` | `bool` | `false` | `switch` | - | Write the relative elevation raster (RER): elevation above the valley centreline's baseline, inside the valley mask. | Output file: &lt;stem&gt;_relative_elevation.tif. Like HAND, but measured against the centreline's baseline (centreline_baseline_algorithm_choice) rather than the channel the pixel's flow reaches, so a meandering or wall-hugging channel does not imprint itself on the field (LSDTopoTools2 lsdtt-basic-metrics get_relative_elevation_raster). The larger channel wins: each centreline node owns the disc of radius distance-to-the-valley-wall round it (inside the valley by construction); a pixel covered by several centrelines' discs takes the one whose nearest covering node has the largest drainage area, then the nearest node of it. Nodes inside a larger centreline's discs are dropped, so a tributary's reference stops where it enters its parent's valley. Pixels no disc covers take the node nearest along the ground through the mask. Pixel elevations are the DEM's; the baseline's are the routing (filled) DEM's channel elevations. Works with source points and with run_valley_width_network. |
| `write_relative_elevation_owner_raster` | `bool` | `false` | `switch` | `write_relative_elevation_raster` | Write which centreline each relative-elevation pixel is measured against. | Output file: &lt;stem&gt;_relative_elevation_owner.tif (source_id / segment_id; nodata outside the valley). For checking where the larger channel takes over. |
| `centreline_baseline_algorithm_choice` | `string` | `projected_monotonic` | `param` | - | How a centreline's baseline is pinned to its channel: projected_monotonic or channel_crossings. | projected_monotonic: every channel pixel pins the centreline node nearest it (the lowest elevation winning on a shared node), then each pin is lowered to the lowest pin upstream of it and the baseline is linear between pins - it never rises downstream and is never above the channel. channel_crossings (LSDTopoTools2): pins only where the channel crosses the centreline in plan view, at the channel elevation interpolated at the crossing (and at centreline nodes that are channel pixels), linear between crossings and flat beyond the first and last, with no other constraint. Used by the baselines in the centreline outputs, the baseline refinement and the relative elevation raster. |

### Valley identification

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `refine_max_height_metres` | `float` | `0.0` | `param` | `refine_valley_with_baseline` | Height above the pinned baseline below which a pixel is valley [m]. | Required by refine_valley_with_baseline, which does nothing while this is 0. It used to reuse the HAND break this tool derived itself; identification has moved to lsdtt-fluvial-surfaces, so the threshold is now stated rather than inherited. Refusing to run without one is deliberate: a refinement against a fabricated threshold would look like a result and be one only by accident. The refined mask is written for you and is never fed back to lsdtt-fluvial-surfaces. |
| `refine_valley_with_baseline` | `bool` | `true` | `switch` | - | Re-test the valley mask against the centreline baseline (relative elevation) after routing. | Removes the stale-reference wedges the nearest-channel-pixel test grows onto downstream-sloping terraces where the channel veers away. |

### Trough carving / routing

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `trough_scaling_factor` | `float` | `0.1` | `param` | - | Trough depth per metre of distance to the valley wall. | Doubled automatically on routing failure (up to 2 retries). |
| `trough_along_channel_slope` | `float` | `0.02` | `param` | - | Synthetic downstream tilt of the trough (m per m along channel). | Stops near-flat / through valleys from routing the centreline out the wrong (upstream) tile edge. |
| `centreline_loops` | `int` | `5` | `param` | - | Trough carve/breach/fill iterations. | - |
| `min_bank_elevation_window_radius_metres` | `float` | `200` | `param` | - | Radius of the local-minimum filter on the nearest-wall elevation field (minimum bank elevation). | - |
| `trough_nodata_buffer_pixels` | `int` | `2` | `param` | - | Valley hole is dilated by this many pixels before bank elevations are sampled. | - |
| `trough_fill_min_slope` | `float` | `0.0001` | `param` | - | Minimum slope for pit filling of the trough surface. | - |
| `min_centreline_nodes` | `int` | `10` | `param` | - | A routed centreline shorter than this counts as a failure. | - |
| `n_threads` | `int` | `0` | `io` | - | Number of OpenMP threads (0 = all available). | - |

### Outputs

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_valley_mask_raster` | `bool` | `true` | `switch` | - | Write the valley mask raster (0 outside; source_id + 1 in each channel's valley). | Output suffix: _valley_mask_channels.tif (UInt32: source_id + 1, 0 = not valley). This is the SUPPLIED mask partitioned by which channel claimed each cell, not a copy of it - and deliberately not called _valley_mask.tif, which is lsdtt-fluvial-surfaces' output and would be overwritten by this one under a shared write_fname. |
| `write_hillshade_raster` | `bool` | `false` | `switch` | - | Write a standard single-azimuth hillshade of the DEM (a natural-looking terrain backdrop for the cross-section / feature figures). | Output suffix: _hillshade.tif. |
| `write_multidirectional_hillshade_raster` | `bool` | `false` | `switch` | - | Write a multidirectional hillshade of the DEM (mean of hillshades from several evenly spaced azimuths). It reduces directional bias and highlights slope breaks (valley/terrace margins) well. | Output suffix: _multidirectional_hillshade.tif. |
| `multidirectional_hillshade_n_directions` | `int` | `8` | `param` | - | Number of evenly spaced azimuths for the multidirectional hillshade. | - |
| `write_valley_centreline` | `bool` | `true` | `switch` | - | Write the valley centrelines as vector lines. | Output suffix: _valley_centreline.fgb (one line per channel with source_id). |
| `write_valley_centreline_nodes` | `bool` | `true` | `switch` | - | Write per-node centreline table. | Output suffix: _valley_centreline_nodes.csv with source_id, node_seq, x, y, elevation, distance_m. |
| `write_valley_cross_sections` | `bool` | `true` | `switch` | - | Write regularly-spaced cross sections orthogonal to each valley centreline. | Outputs _valley_cross_sections.fgb (transect lines), _valley_cross_section_sides.fgb (labelled left/right valley-side points), _valley_cross_section_boundaries.fgb (valley sides + hydraulic step-midpoint boundary points), and _valley_cross_section_samples.csv (per-sample elevation profiles). Each channel's transects share one length = cross_section_length_multiplier x the widest valley on that channel. |
| `cross_section_spacing_metres` | `float` | `500` | `param` | - | Along-centreline spacing between cross-section stations (m). | - |
| `cross_section_length_multiplier` | `float` | `1.5` | `param` | - | Cross-section total length as a multiple of the widest valley width for the channel. | 1.5 spans 1.5x the widest valley so every section shows both walls with margin. |
| `cross_section_tangent_window_nodes` | `int` | `8` | `param` | - | Centreline nodes on each side used to smooth the local direction the transect is drawn orthogonal to. | - |
| `cross_section_min_half_length_metres` | `float` | `0` | `param` | - | Minimum transect half-length (m). Floors the transect length so a narrow channel's sections still extend a sensible distance beyond the valley to show wall/upland context. | 0 = length purely from valley width; set e.g. 800-1500 for narrow tributaries. Still capped by cross_section_max_half_length_metres. |
| `cross_section_max_half_length_metres` | `float` | `5000` | `param` | - | Hard cap (m) on the outward march when finding a valley side, and hence on any transect half-length. | Where the valley mask connects to a reservoir / coastal flat / the sea the in-valley run through the centre can span the whole tile; a station still in-valley at this cap has that side marked not found and is excluded from the widest-valley statistic. Raise for genuinely very wide valleys. |
| `clip_overlapping_sections` | `bool` | `true` | `switch` | - | Phase-1 non-overlap: clip each transect side where that ray becomes closer to a neighbouring station centre than its own centre (local Voronoi bisector clipping). | Removes concave-side overlap while keeping straight orthogonal transects. |
| `curved_cross_sections` | `bool` | `false` | `param` | - | Trace each cross section as an isocontour of a smooth along-valley (streamwise) coordinate field instead of a straight ray: straight in straight reaches, curving gradually through bends like the constant-streamwise grid lines of a boundary-fitted mesh. Sections keep a uniform length and a greedy pass decimates stations in tight bends so survivors stay long. | Boundary-fitted valley-floor refinement. Sides are clipped at the Voronoi bisector of the nearest kept neighbours, never below the valley wall plus curved_wall_clearance_fraction. |
| `curved_wall_clearance_fraction` | `float` | `0.10` | `param` | - | For curved cross sections: a shortened side never clips below (1 + this) x the distance to the valley wall, so a section always clears the pre-existing valley by this margin. | Only used when curved_cross_sections is true. |
| `curved_streamwise_smoothing_passes` | `int` | `40` | `param` | - | For curved cross sections: number of anchored Laplacian smoothing passes applied to the streamwise coordinate field. Higher = smoother sections with fewer creases at bends. | Only used when curved_cross_sections is true. |
| `curved_coarse_spacing_metres` | `float` | `700` | `param` | - | For curved cross sections: sections are placed coarse-to-fine - a sparse backbone at this spacing first, then midpoints between accepted sections (this halved, quartered, ... down to cross_section_spacing_metres). Keeps coverage even instead of dense clusters next to gaps. | Rounded down to a power-of-two multiple of the base spacing. Set &lt;= the base spacing to disable the hierarchy. Only used when curved_cross_sections is true. |
| `neighbour_clip_window` | `int` | `2` | `param` | - | Number of neighbouring stations on each side considered for overlap clipping. | 2 checks s+/-1 and s+/-2; larger is more conservative. |
| `min_section_reach_fraction` | `float` | `0.35` | `param` | - | Drop a station if the clipped transect length is below this fraction of the nominal full length. | Avoids tiny degenerate sections in tight hairpins. |
| `write_cross_section_features` | `bool` | `true` | `switch` | - | Along each cross section, detect and write features: D8 channel crossings, left/right valley walls, hydraulic-step boundaries, maximum-curvature peaks and multidirectional-hillshade peaks. | Outputs _valley_cross_section_features.csv (one row per feature) and adds curvature / hillshade / is_channel columns to the per-sample CSV. Needs the DEM to derive curvature and hillshade; the channel mask comes from the resolved channel scope. |
| `write_cross_section_hydraulics` | `bool` | `false` | `switch` | - | For every cross section, write the BankfullMapper-style hydraulic depth function: the profile is flooded stage by stage and, at each stage, the flow width, wetted area, wetted perimeter and hydraulic depth (area-to-width ratio A/W) are recorded. Breakpoints in the hydraulic depth vs elevation-above-thalweg curve mark bankfull / inset flat-surface levels (Delchiaro et al. 2025, Computers &amp; Geosciences, doi:10.1016/j.cageo.2025.106001). | Outputs _valley_cross_section_hydraulics.csv (one row per stage per station: stage_m, depth_above_thalweg_m, width_m, area_m2, wetted_perimeter_m, hydraulic_depth_m). |
| `write_cross_section_levels_csv` | `bool` | `false` | `output` | - | Cluster every section's hydraulic breaks into levels per channel and write them. | Output file: &lt;stem&gt;_cross_section_levels.csv (source_id, level_hand_m, n_sections, spread_m, status, ...). The breaks of all sections on a channel are pooled and split by the ladder's own rule: the largest number of levels whose medians are at least ladder_resolvability_kappa * vertical_accuracy_metres apart. Feed the file to lsdtt-fluvial-surfaces as supplied_levels_fname (fluvial_surface_algorithm_choice cross_section_levels), with the same source points and threshold so the source ids match. Needs write_cross_section_features (the breaks come from that path). |
| `cross_section_level_min_section_fraction` | `float` | `0.5` | `param` | `write_cross_section_levels_csv` | Fraction of a channel's sections a level must appear on to be kept. | A tread is crossed by most sections; a step only a few sections see - one running up a riser gap's ramp, or onto a valley wall - is not a level of the channel. Dropped levels are listed in the run log. |
| `write_cross_section_informed_valley_mask_raster` | `bool` | `false` | `switch` | - | Build a valley-floor mask from the cross sections' HYDRAULIC BREAKS instead of the corridor HAND histogram: each corridor pixel's height above the nearest centreline node's thalweg baseline is compared with the channel's MEDIAN lowest hydraulic break (a single robust floodplain level, so an anomalous section cannot streak or leak onto the plateau). Needs write_cross_section_features (the breaks come from that path). | Outputs _cross_section_informed_valley_mask.tif (Byte; 0 = above the break, 1 = at/below the median lowest break [floodplain]) plus per-channel _cross_section_informed_valley_mask_ch&lt;sid&gt;.tif. |
| `cross_section_curvature_radius_metres` | `float` | `30` | `param` | - | Polynomial-fit window radius (m) for the maximum-curvature raster used to find slope-break peaks along the cross sections. | Larger smooths over fine breaks (fewer, coarser peaks); smaller resolves finer risers but is noisier. |
| `cross_section_feature_prominence_fraction` | `float` | `0.15` | `param` | - | Minimum prominence of a curvature / hillshade peak as a fraction of that profile's value range over the transect, for it to be reported. | Raise to keep only the strongest breaks; lower to catch subtler ones. |
| `write_output_format` | `string` | `cog` | `io` | - | Raster output format. | Allowed values: cog (default), gtiff, same_as_input. Set gtiff explicitly to retain plain GeoTIFF output. |

### Ladder

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `ladder_max_level_count` | `int` | `12` | `param` | - | Ceiling on the level count considered under auto selection. | A ceiling on what one reach can plausibly carry, not a storage limit. A reach resolving more than a dozen distinct fluvial levels is more likely over-fitting noise than recording twelve terraces. Raise it freely for a genuine long Quaternary staircase; there is no numeric cliff behind it. |
| `ladder_break_algorithm_choice` | `string` | `antimode` | `param` | - | Break finder: gap, antimode or jenks. | antimode finds dips in a smoothed log10 density. jenks is Fisher-Jenks on the linear values; it takes the level count directly and is the source of the reported goodness-of-variance fit. gap finds empty spans between clusters. |
| `vertical_accuracy_metres` | `float` | `0.0` | `param` | - | Vertical accuracy of the DEM, E_v [m]. 0 disables the floor. | Two levels are only claimed as distinct when their band centres are at least ladder_resolvability_kappa * E_v apart. This is metadata you must supply: it cannot be derived from the DEM. 0 disables the check and every window then reports resolved, carrying no instrument guarantee at all. Note that a DEM whose vertical units are feet on a metre grid has an effective E_v 3.28 times larger than you think. |
| `ladder_resolvability_kappa` | `float` | `3.0` | `param` | - | Multiplier on E_v in the resolvability floor. | Required band-centre gap = kappa * vertical_accuracy_metres. |

### Interior nodata

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `interior_nodata_method` | `string` | `sink` | `param` | - | How nodata pixels inside the DEM are treated: sink or fill. | Interior nodata are nodata pixels not connected to the raster edge through other nodata (isolated voids, moulins, masked lakes). sink (the default, unchanged behaviour): they are outlets, so flow that reaches them ends there, as at the DEM edge; use this for real internal drainage (e.g. moulins on an ice sheet). fill: they are interpolated with nodata_fill_method_choice before any filling or routing, so flow passes across them; use this for data voids. Nodata connected to the edge (e.g. sea removed by min_elevation) is never filled and always acts as an open boundary. |
| `nodata_fill_method_choice` | `string` | `laplacian` | `param` | `interior_nodata_method=fill` | Method for filling interior nodata pixels. | laplacian: smooth Dirichlet interpolation (recommended). iterative_weighted_average: fast, well-suited to thin voids. Used when interior_nodata_method is fill. |
| `nodata_fill_max_iter` | `int` | `10000` | `param` | `interior_nodata_method=fill` | Maximum iterations for iterative nodata fill methods. | Used when interior_nodata_method is fill. |

<!-- END GENERATED CLI PARAMETERS -->
