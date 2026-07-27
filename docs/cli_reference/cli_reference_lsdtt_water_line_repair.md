# CLI reference: lsdtt-water-line-repair

Repairs fragmented water-line polylines (for example OS VectorMap Local channels truncated at road crossings; an example parameter file for the OS Vector Map Local is provided below). The tool is pure vector: it reads a water line layer and a road-centreline layer, finds dangling water tips, bridges plausible pairs, and writes a repaired line layer containing originals plus inserted bridge segments.

## Usage

```
lsdtt-water-line-repair <path> <parameter_file>
lsdtt-water-line-repair <parameter_file>         # path defaults to ./
lsdtt-water-line-repair -h                       # print help and exit
lsdtt-water-line-repair -v                       # print version and exit
lsdtt-water-line-repair -c                       # print citation and exit
```

## Pipeline overview

| Stage | Description |
|---|---|
| 1 | Load water and road line layers (must share CRS) |
| 2 | Build endpoint graph; keep degree-1 (dangling) water tips |
| 3 | Estimate outward tangent at each tip |
| 4 | Build tip-pair candidates with distance / facing / heading gates |
| 5 | Mark candidate as road-associated if connector intersects or lies near a road |
| 6 | Select one-to-one accepted bridges (optionally mutual-best only) |
| 7 | Write repaired layer (originals + bridges) |
| 8 | Optionally write tip diagnostics |

---

## Input / output parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Base path for relative input paths. |
| `read_fname` | string | parameter file stem | Not used directly by this CLI; present for lsdtt3 consistency. |
| `write_path` | string | same as `read_path` | Output directory. |
| `write_fname` | string | same as `read_fname` | Output filename stem. |

---

## Required layer inputs

| Parameter | Type | Default | Description |
|---|---|---|---|
| `input_vector` | string | _(empty)_ | Dataset containing the water layer (and optionally roads). Required. Resolved against `read_path` if relative. |
| `water_layer` | string | _(empty)_ | Name of the water LINESTRING layer to repair. Required. |
| `road_layer` | string | _(empty)_ | Name of the road-centreline LINESTRING layer used for road-crossing association. Required. |
| `road_vector` | string | _(empty)_ | Optional separate dataset holding the road layer. If empty, roads are read from `input_vector`. |

Both layers must be in the same CRS. A projected metre CRS is strongly recommended.

---

## Bridge candidate geometry and scoring

| Parameter | Type | Default | Description |
|---|---|---|---|
| `snap_tolerance` | float | `0.5` | Endpoint coincidence tolerance (m) used for degree counting. |
| `tangent_fit_length` | float | `15.0` | Inward arc length (m) used to estimate each tip's outward tangent direction. |
| `max_gap_length` | float | `25.0` | Maximum gap length (m) for road-associated bridges. |
| `max_nonroad_gap` | float | `30.0` | Maximum gap length (m) for non-road bridges. Set `0` to disable non-road bridging. |
| `max_facing_angle` | float | `45.0` | Strict (non-road) gate: each tip must face its partner within this angle (deg). |
| `max_heading_misalignment` | float | `50.0` | Strict (non-road) gate: tangents should be near anti-parallel (`180 deg`); this sets max misalignment (deg). |
| `road_max_facing_angle` | float | `60.0` | Road tier gate: at least one tip must face within this angle (deg). No heading constraint in road tier. |
| `road_intersection_tolerance` | float | `10.0` | Distance (m) used to classify a connector as road-associated when there is no strict intersection. |
| `require_mutual_best` | bool | `true` | Accept a bridge only if each tip is the other's best-scoring partner. Strongly reduces false joins. |

Selection is one-to-one: once a tip is used in an accepted bridge, it cannot be reused.

---

## Output writing options

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_output_format` | string | `flatgeobuf` | Output format for repaired lines. Supported here: `flatgeobuf`/`fgb`, `gpkg`. |
| `calculate_spatial_index` | bool | `true` | Build spatial index on the repaired output layer. |
| `road_name_field` | string | `roadName` | Road attribute copied to `crossed_road_name` for road bridges. Falls back internally to `roadNumber` when missing/empty. |
| `merge_chains` | bool | `false` | Placeholder switch for dissolving connected chains. **Not implemented**; leave false. |

---

## Optional diagnostics

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_tip_diagnostics` | bool | `false` | Write all dangling tips as points to `<write_fname>_repair_tips.fgb`, with bridge flag and nearest-road distance. |

---

## Output files

| File | Written when | Description |
|---|---|---|
| `<write_fname>_water_line_repaired.<ext>` | always | Repaired line layer (original features plus inserted bridge segments). `<ext>` from `write_output_format`. |
| `<write_fname>_repair_tips.fgb` | `write_tip_diagnostics: true` | Point diagnostics for dangling tips. |

---

## Repaired layer schema

The repaired line output includes fields:

- `repair_type` (`original`, `bridge_road`, `bridge_gap`)
- `bridge_length_m` (NULL for originals)
- `crossed_road_name` (NULL for non-road bridges/originals)
- `src_fid` (copied from input `Fid` where available)
- `src_featcode` (copied from input `FeatCode` where available)

---

## Tip diagnostics schema

`<write_fname>_repair_tips.fgb` fields:

- `bridged` (1 if that tip was used in an accepted bridge)
- `road_dist_m` (distance to nearest road segment)
- `line_index`
- `at_start` (1 for line start tip, 0 for line end tip)
- `tx`, `ty` (outward unit tangent components)

---

## Example parameter file

```
read_path: ./
write_path: ./lsdtt3_outputs/
write_fname: nr73nw_vml

input_vector: vml-nr73nw.gpkg
water_layer: vml-nr73nw_Water_Line
road_layer: vml-nr73nw_Road_Centreline

write_output_format: gpkg
calculate_spatial_index: true

snap_tolerance: 0.5
max_gap_length: 25.0
max_nonroad_gap: 30.0
max_facing_angle: 45.0
max_heading_misalignment: 50.0
road_max_facing_angle: 60.0
tangent_fit_length: 15.0
road_intersection_tolerance: 10.0
require_mutual_best: true

write_tip_diagnostics: true
```
