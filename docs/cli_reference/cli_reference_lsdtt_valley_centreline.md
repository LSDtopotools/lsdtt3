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
[Floodplains, terraces and valleys](../fluvial_surfaces_and_valleys.md).

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
mode; `source_points_path` is required.

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
| `write_output_format` | string | `gtiff` | Raster format: `cog`, `gtiff`, `same_as_input`. |

---

## Source and outlet points

Supply either a source-only file, a source file plus an explicit outlet file, or
a single combined file with `auto_detect_outlet_points: true` (the default),
which classifies the downstream points as outlets. Point files may be CSV
(`latitude` / `longitude` in EPSG:4326) or any OGR vector (`.gpkg`, `.fgb`, …).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `source_points_path` | string | | **Required.** Vector or CSV file with source points. |
| `source_points_input_format` | string | `auto_detect` | `auto_detect`, `csv_latitude_longitude_epsg4326`, `ogr_vector_layer`. |
| `source_points_layer_name` | string | | Optional layer name for vector input. |
| `source_points_latitude_field_name` | string | `latitude` | CSV latitude column. |
| `source_points_longitude_field_name` | string | `longitude` | CSV longitude column. |
| `outlet_points_path` | string | | Optional outlet points file, to truncate channels. |
| `outlet_points_input_format` | string | `auto_detect` | As the source format. |
| `outlet_points_layer_name` | string | | Optional layer name. |
| `outlet_points_latitude_field_name` | string | `latitude` | CSV latitude column. |
| `outlet_points_longitude_field_name` | string | `longitude` | CSV longitude column. |

> **Naming note.** This tool still spells these `source_points_path` /
> `outlet_points_path`, while `lsdtt-fluvial-surfaces` spells the same concept
> `source_points_fname` / `outlet_points_fname`. The two are not interchangeable;
> use the spelling of the tool you are running.

## Channel scope

| Parameter | Type | Default | Description |
|---|---|---|---|
| `threshold_contributing_pixels` | int | `1000` | Contributing-pixel threshold defining channel cells for source snapping. |
| `auto_detect_outlet_points` | bool | `true` | Classify points in `source_points_path` into sources and outlets. Disabled when an explicit outlet file is given. |
| `snap_sources_to_channel_threshold` | bool | `true` | Snap source points to the nearest cell meeting the threshold. |
| `source_snap_max_radius_pixels` | int | `25` | Maximum expanding-square search radius for snapping. |
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
| `cross_section_curvature_window_metres` | float | `60` | Polynomial-fit window diameter (m) for the maximum-curvature raster used to find slope-break peaks. |
| `cross_section_feature_prominence_fraction` | float | `0.15` | Minimum peak prominence, as a fraction of that profile's range over the transect. |
| `write_cross_section_hydraulics` | bool | `false` | Flood each profile stage by stage and record flow width, wetted area, wetted perimeter and hydraulic depth (A/W). Breakpoints in hydraulic depth against elevation above thalweg mark bankfull and inset flat-surface levels (Delchiaro et al. 2025, *Computers & Geosciences*, doi:10.1016/j.cageo.2025.106001). |
| `write_cross_section_informed_valley_mask` | bool | `false` | Build a valley-floor mask from the sections' hydraulic breaks: each corridor pixel's height above the nearest node's thalweg baseline is compared with the channel's **median** lowest hydraulic break, so one anomalous section cannot streak onto the plateau. Requires `write_cross_section_features`. |

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

- [Floodplains, terraces and valleys](../fluvial_surfaces_and_valleys.md) — the
  trough construction, its origin in Clubb et al. (2022), and where lsdtt3
  departs from it
- [CLI reference: lsdtt-fluvial-surfaces](cli_reference_lsdtt_fluvial_surfaces.md)
- [CLI reference: lsdtt-flat-surface-finder](cli_reference_lsdtt_flat_surface_finder.md)
- [CLI reference: lsdtt-swath-profile](cli_reference_lsdtt_swath_profile.md)
