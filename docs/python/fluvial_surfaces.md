# Floodplains, terraces and valleys

The Python bindings expose the steps behind four command-line tools:
`lsdtt-flat-surface-finder` (flat ground), `lsdtt-fluvial-surfaces` (floodplain
and terrace classes), `lsdtt-valley-centreline` (valley centrelines, relative
elevation and terrace profiles) and the riser search of `lsdtt-feature-finder`.
Each function calls the C++ routine that the tool calls. Where
`lsdtt-flat-surface-finder` and `lsdtt-fluvial-surfaces` spell a keyword
differently, the bindings use the `lsdtt-fluvial-surfaces` spelling, and the
older names are accepted (see the table under the floodplain classes).

The method presets of `lsdtt-fluvial-surfaces` (`ladder`, `clubb2017`) are
not single calls here. A script chooses the steps itself, as the examples
below do.

The examples assume a filled DEM and its D8 flow:

```python
import lsdtt3

dem = lsdtt3.Raster.read("dem.tif")
filled = lsdtt3.fill_pits(dem)
flow = lsdtt3.FlowInfo(filled)
scope = flow.channel_scope(threshold_contributing_pixels=500)
```

## Flat ground

`find_flat_surfaces` labels ground that is smooth and nearly level:

```python
flat = lsdtt3.find_flat_surfaces(dem, relief_window_radius_metres=50.0,
                                 max_relief_metres=1.0)
flat.labelled_flat_mask    # Raster: one label per flat
```

`relief_window_radius_metres` is the radius of the local-relief window, in
metres, rounded to whole cells. The default, 0, is a 5 by 5 cell window.
`max_relief_metres` is the local relief (maximum minus minimum) at or above
which a cell is not flat. `lsdtt-flat-surface-finder` still reads
`relief_window_side_pixels` (a side in cells) and `relief_threshold_metres`
for these two settings; both older names are accepted here too (see the
table under the floodplain classes).

## Floodplain and terrace classes

`classify_flat_surfaces_with_scope` splits the flats into floodplain and
terrace by height above the channel (HAND, in metres):

```python
classes = lsdtt3.classify_flat_surfaces_with_scope(
    flat.labelled_flat_mask, filled, scope,
    max_floodplain_hand_metres=2.0,
    min_terrace_hand_metres=2.0,
    valley_corridor_half_width_metres=5000.0,
)
classes.surface_class      # Raster: 0 none, 1 floodplain, 2 terrace
classes.table              # one row per flat label
classes.banded_labels      # Raster: patch labels (1..N, 0 = background)
classes.patches            # one row per patch: label, band, n_pixels, area_m2, median_field
classes.ladder_windows     # one row per reach window: breaks, status, ...
```

The threshold values in the example are placeholders; choose them for the
landscape. A channel-adjacent flat is floodplain when its median HAND is at
most `max_floodplain_hand_metres`. A flat linked to a channel is terrace when
its median HAND is at least `min_terrace_hand_metres` and it is not
floodplain. A terrace must have a cell within `valley_corridor_half_width_metres`
of a channel, which excludes upland flats; 0 means no limit.

These keywords replace older names. The old names still work, with a
`DeprecationWarning`:

| Old name | Name now |
|---|---|
| `floodplain_max_height_metres` | `max_floodplain_hand_metres` |
| `terrace_min_height_metres` | `min_terrace_hand_metres` |
| `terrace_max_distance_to_channel_metres` | `valley_corridor_half_width_metres` |
| `relief_threshold_metres` | `max_relief_metres` |
| `relief_window_side_pixels` | `relief_window_radius_metres` |

`relief_window_side_pixels` gave the window as a full side in cells, which
must be odd and at least 3; the older `relief_window_pixels` means the same
side. The radius is `(side - 1) / 2` cells, converted to metres with the DEM
cell size.

Three defaults differ from the command line, so set them to the CLI values to
reproduce a CLI run:

| Keyword | Python default | CLI default |
|---|---|---|
| `min_band_component_pixels` | 1 | 50 |
| `ladder_unresolved_window_algorithm_choice` | `single_band` | `nearest_resolved` |
| `valley_corridor_half_width_metres` | 0 (no limit) | 5000 |

The ladder keywords (`ladder_*`, `vertical_accuracy_*`, `clip_corridor_per_pixel`,
`max_band_spread_metres`, `segmentation_algorithm_choice`) keep the CLI's names
and meanings. `help(lsdtt3.classify_flat_surfaces_with_scope)` lists them all.

Two further options:

- **Training polygons.** `training_surfaces_fname` names a vector file of
  polygons classed `floodplain`, `terrace` or `terraceN`. They set the height
  breaks in each reach window. Their statistics are in
  `classes.training_polygons`.
- **Whole-mask segmentation.** `segmentation_algorithm_choice="whole_mask"`
  segments the admitted ground as one mask, with `riser_mask` cells as
  barriers, for example `find_risers(...).centreline` (see below).

## Floodplain thresholds

`select_threshold` picks a threshold from a sample of values. A Raster sample
is its valid cells, so the example takes every HAND cell; mask the HAND raster
to the cells of interest first if the threshold should come from them only.
`classify_floodplain_terraces_by_thresholds` then applies thresholds in the
Clubb et al. (2017) style. It takes HAND and the channel mask of the scope,
and the admitted ground (here, the flats):

```python
hand = flow.height_above_channel(filled, scope.mask)
admitted = flat.labelled_flat_mask
bound = lsdtt3.select_threshold(hand, "qq_gaussian")
if bound["found"]:
    thresholds = lsdtt3.classify_floodplain_terraces_by_thresholds(
        admitted, hand, scope.mask,
        max_floodplain_hand_metres=bound["threshold"])
    thresholds.floodplain_mask
    thresholds.terrace_mask
```

`threshold_algorithm_choice` accepts `fixed`, `qq_gaussian`, `cdf_divergence`
and `qq_upper_departure`. `found` is False when no threshold could be chosen;
`threshold` is then not meaningful.

## Level tests

`cluster_surface_levels` tests whether a set of heights sits on levels.
Heights are split wherever neighbours are at least `required_gap_metres`
apart. The result gives the clusters and a p-value against uniformly drawn
heights:

```python
levels = lsdtt3.cluster_surface_levels(classes.patches["median_field"],
                                       required_gap_metres=0.5)
levels["clusters"]     # DataFrame: lower, upper and centre of each level
levels["p_value"]
```

`required_gap_metres` is typically `ladder_resolvability_kappa` times
`vertical_accuracy_metres`: the gap that two levels must have to be resolved.

## Comparing and scoring classes

Surface-class rasters from several methods, or from a mapped survey, can be
compared:

```python
table = lsdtt3.compare_surface_classes({"ladder": a, "clubb": b})
bits = lsdtt3.surface_class_agreement({"ladder": a, "clubb": b}, "floodplain")
score = lsdtt3.score_surface_classes(
    a, "mapped_surfaces.gpkg", mapped_surfaces_class_field_name="class",
    mapped_surfaces_extent_fname="survey_extent.gpkg")
score["floodplain_iou"]
```

`compare_surface_classes` gives, for each pair of methods, the intersection
over union of each class and, of the cells either method classed, the fraction
both gave the same class.
`surface_class_agreement` gives a bit field: bit k is set where method k gave
the class. `score_surface_classes` scores against mapped polygons; cells
outside `mapped_surfaces_extent_fname` are not counted.

`relief_above_nearest_baseline` gives the height of each cell above its
nearest baseline cell, for example the channel mask. It is the swath-baseline
terrace reference of `lsdtt-fluvial-surfaces`.

## Risers

`find_risers` finds short, steep steps between gentler surfaces, such as a
terrace edge, and returns them as line segments:

```python
risers = lsdtt3.find_risers(dem, search_mask=valley_mask)
risers.centreline          # Raster: 1 on riser cells
risers.segments            # GeoDataFrame: length, step height, face gradient, geometry
```

`search_mask` limits the search to cells greater than 0.5, for example a
valley mask. `min_step_height_metres` (default 1.0 m) is the smallest step
kept. `score_risers_against_mapped(risers.centreline, "mapped_edges.gpkg")`
scores the detected cells against mapped terrace-edge lines, within
`mapped_risers_tolerance_metres` (default 20 m).

## Valley centrelines and terraces

Valley centrelines run along a valley floor, one per channel segment:

```python
valley_mask = lsdtt3.Raster.read("valley_mask.tif")   # cells > 0.5 are valley
network = lsdtt3.build_valley_width_network(filled, scope, valley_mask)
network.segments           # one row per channel segment
centrelines = network.centrelines

rel_elev, owner = centrelines.relative_elevation(filled)
bins, patch_table = centrelines.terrace_profile(
    filled, classes.banded_labels, terrace_profile_bin_metres=100.0)
```

`relative_elevation` returns the height of each cell above its centreline
baseline, and the `source_id` of the centreline that owns it. Nodata is
outside the domain, which is the valley floor unless a `domain_mask` is given.
`terrace_profile` bins the patches along each centreline. Both need the
centrelines to come from a `ChannelScope` or from `build_valley_width_network`,
which rank the channels by drainage area.

Centrelines from a single scope use `extract_valley_centrelines(filled, scope,
valley_mask=valley_mask)` instead.

Cross sections across the centrelines give surface levels for the
classification:

```python
sections = lsdtt3.extract_valley_cross_sections(filled, centrelines,
                                                spacing_metres=500.0)
levels = sections.levels()
supplied = {int(sid): list(group["level_hand_metres"])
            for sid, group in levels.groupby("source_id")}

classes = lsdtt3.classify_flat_surfaces_with_scope(
    flat.labelled_flat_mask, filled, scope, supplied_levels=supplied)
```

`levels()` pools the breaks of each channel's cross sections and splits them
into levels, as `write_cross_section_levels_csv` does. `supplied_levels` takes
the heights per channel `source_id`; the classification bands at the midpoints
between them instead of finding levels itself.

See the [Floodplains, terraces and valleys API](../python_api/fluvial.md) for
the return objects and every keyword.
