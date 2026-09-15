# lsdtt-outlet-matching

## Intent

`lsdtt-outlet-matching` transfers outlet locations from a reference DEM to a
target or resampled DEM. It matches D8 contributing-basin areas, not channel
networks or a channel-area threshold. The matched CSV is suitable as an outlet
input for `lsdtt-ridgeline-analysis` and `lsdtt-basin-extraction`; it contains
the required `latitude` and `longitude` columns.

## Pipeline

1. Read both DEMs and the reference latitude/longitude CSV.
2. Build independent D8 routing graphs with topological order and s-vectors.
3. Transform reference outlet coordinates and validate approximate DEM extent overlap.
4. Compute each reference outlet's basin area as its reference s-vector size
   multiplied by reference cell area.
5. Build target contributing-pixel counts and edge-influence values with D8
   accumulation.
6. Search valid target cells within `local_search_radius_metres`. Select the
   non-edge-incomplete candidate with the smallest absolute area difference;
   distance is used only to break area ties.

## Usage

```text
lsdtt-outlet-matching <path> <parameter_file>
lsdtt-outlet-matching <parameter_file>
lsdtt-outlet-matching -h
```

The input CSV must contain `latitude` and `longitude` in EPSG:4326. DEMs must
have projected CRS metadata with compatible linear units, normally metres.

## Parameters

| Parameter | Meaning |
|---|---|
| `reference_dem_fname` | DEM defining the source drainage masks |
| `target_dem_fname` | DEM on which outlets are selected |
| `reference_outlet_fname` | Reference outlet CSV |
| `local_search_radius_metres` | Candidate search radius |
| `raster_is_filled` | Whether both DEMs already satisfy acyclic routing requirements |
| `matched_outlet_fname`, `diagnostics_fname` | Optional output names |

Output defaults are `<write_fname>_matched_outlets.csv` and
`<write_fname>_outlet_matching_diagnostics.csv`.

## Matching and outputs

The reference area is `reference_flow->s_vector_size(reference outlet flat)`
multiplied by reference cell area. Each target area is its D8 contributing-pixel
count multiplied by target cell area. Candidate target cells must be valid, lie
within the local search radius, and have no upstream influence from the target
raster edge. Among those candidates, the minimum absolute target-minus-reference
area difference wins, with distance as the tie-break.

The diagnostics CSV contains `point_id`, candidate row/column/flat,
`distance_metres`, target and reference areas, signed area difference,
`area_difference_fraction`, `edge_incomplete`, and `accepted`. The matched CSV
contains the selected target coordinates and indices, both areas, signed and
fractional area differences, search distance, and status. No mask rasterization
or precision, recall, F1, or IoU threshold is used.

## Example

```text
read_path: ./data
write_path: ./out
read_fname: ignored.tif
write_fname: outlet_match
reference_dem_fname: reference.tif
target_dem_fname: target.tif
reference_outlet_fname: reference_outlets.csv
local_search_radius_metres: 750
```
