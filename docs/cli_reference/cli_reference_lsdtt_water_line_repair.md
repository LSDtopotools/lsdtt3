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

!!! note "Multipart lines"
    `input_vector_fname`, `road_vector_fname` accept multipart geometry (`MultiLineString` and friends). Each feature is split into its parts on read, with the attributes copied onto each, and the split is reported in the log.

    Shapefiles are the commonest source of multipart data, and **versions before August 2026 discarded it silently** — see [Vector data](../what_is_lsdtt3/vector_data.md).

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
| `input_vector_fname` | string | _(empty)_ | Dataset containing the water layer (and optionally roads). Required. Resolved against `read_path` if relative. |
| `water_layer` | string | _(empty)_ | Name of the water LINESTRING layer to repair. Required. |
| `road_layer` | string | _(empty)_ | Name of the road-centreline LINESTRING layer used for road-crossing association. Required. |
| `road_vector_fname` | string | _(empty)_ | Optional separate dataset holding the road layer. If empty, roads are read from `input_vector_fname`. |

Both layers must be in the same CRS. A projected metre CRS is strongly recommended.

---

## Bridge candidate geometry and scoring

| Parameter | Type | Default | Description |
|---|---|---|---|
| `snap_tolerance_metres` | float | `0.5` | Endpoint coincidence tolerance (m) used for degree counting. |
| `tangent_fit_length_metres` | float | `15.0` | Inward arc length (m) used to estimate each tip's outward tangent direction. |
| `max_gap_length_metres` | float | `25.0` | Maximum gap length (m) for road-associated bridges. |
| `max_nonroad_gap_metres` | float | `30.0` | Maximum gap length (m) for non-road bridges. Set `0` to disable non-road bridging. |
| `max_facing_angle_degrees` | float | `45.0` | Strict (non-road) gate: each tip must face its partner within this angle (deg). |
| `max_heading_misalignment_degrees` | float | `50.0` | Strict (non-road) gate: tangents should be near anti-parallel (`180 deg`); this sets max misalignment (deg). |
| `road_max_facing_angle_degrees` | float | `60.0` | Road tier gate: at least one tip must face within this angle (deg). No heading constraint in road tier. |
| `road_intersection_tolerance_metres` | float | `10.0` | Distance (m) used to classify a connector as road-associated when there is no strict intersection. |
| `require_mutual_best` | bool | `true` | Accept a bridge only if each tip is the other's best-scoring partner. Strongly reduces false joins. |

Selection is one-to-one: once a tip is used in an accepted bridge, it cannot be reused.

---

## Output writing options

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_vector_format` | string | `flatgeobuf` | Output **vector** format for repaired lines: `flatgeobuf`/`fgb`, `gpkg`. This is the vector key, not the raster key `write_output_format` used by the raster tools. |
| `calculate_spatial_index` | bool | `true` | Build spatial index on the repaired output layer. |
| `road_name_field` | string | `roadName` | Road attribute copied to `crossed_road_name` for road bridges. Falls back internally to `roadNumber` when missing/empty. |
| `merge_chains` | bool | `false` | Placeholder switch for dissolving connected chains. **Not implemented**; leave false. |

`write_vector_format` was called `write_output_format` in earlier releases. The
old name still works as a deprecated alias: the value is carried across and a
warning naming the replacement is printed. Update your parameter files, because
`write_output_format` means a *raster* format in every other CLI.

---

## Optional diagnostics

| Parameter | Type | Default | Description |
|---|---|---|---|
| `write_tip_diagnostics` | bool | `false` | Write all dangling tips as points to `<write_fname>_repair_tips.fgb`, with bridge flag and nearest-road distance. |

---

## Output files

| File | Written when | Description |
|---|---|---|
| `<write_fname>_water_line_repaired.<ext>` | always | Repaired line layer (original features plus inserted bridge segments). `<ext>` from `write_vector_format`. |
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

input_vector_fname: vml-nr73nw.gpkg
water_layer: vml-nr73nw_Water_Line
road_layer: vml-nr73nw_Road_Centreline

write_vector_format: gpkg
calculate_spatial_index: true

snap_tolerance_metres: 0.5
max_gap_length_metres: 25.0
max_nonroad_gap_metres: 30.0
max_facing_angle_degrees: 45.0
max_heading_misalignment_degrees: 50.0
road_max_facing_angle_degrees: 60.0
tangent_fit_length_metres: 15.0
road_intersection_tolerance_metres: 10.0
require_mutual_best: true

write_tip_diagnostics: true
```

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

### Input vectors

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `input_vector_fname` | `string` | - | `io` | - | Input vector dataset (e.g. a GeoPackage) holding both layers. | Required. A bare filename is read from read_path. Both layers must share the same CRS (projected metres recommended). |
| `water_layer` | `string` | - | `param` | - | Name of the water LINESTRING layer to repair. | Required. e.g. vml-nr73sw_Water_Line. |
| `road_layer` | `string` | - | `param` | - | Name of the road centreline LINESTRING layer used as the crossing test. | Required. e.g. vml-nr73sw_Road_Centreline. |
| `road_vector_fname` | `string` | - | `io` | - | Dataset holding the road layer, if different from input_vector_fname. | Use when the water network is in one file (e.g. the centreline .fgb) and roads are in another (e.g. the source .gpkg). Empty = same as input_vector_fname. A bare filename is read from read_path. |
| `road_name_field` | `string` | `roadName` | `param` | - | Road attribute copied into the bridge crossed_road_name field. | Falls back to 'roadNumber' when this field is empty/absent. |

### Outputs

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `write_vector_format` | `string` | `flatgeobuf` | `io` | - | Output VECTOR format: flatgeobuf or gpkg. | Governs every vector layer this tool writes. Distinct from the raster key write_output_format used by the raster tools: this CLI writes vectors only. FlatGeobuf is preferred: it is compact, single-layer, carries its own georeferencing, and does not have GeoPackage's problems writing across drives (a common WSL2 situation). GeoPackage is an SQLite database, and SQLite cannot take the file locks it needs on a drive mounted into WSL2, so a gpkg write there fails. |
| `calculate_spatial_index` | `bool` | `true` | `io` | - | Build a spatial index when writing the output layer. | Speeds GIS queries; buffers the layer in memory before writing. |
| `write_tip_diagnostics` | `bool` | `false` | `switch` | - | Write a diagnostic point layer of every dangling tip. | Output: &lt;write_fname&gt;_repair_tips.fgb. Attributes: bridged (1 if the tip was joined), road_dist_m (distance to nearest road), line_index. Use it to see which tips near roads stayed unbridged and why. |

### Gap bridging

| Parameter | Type | Default | Role | Parent | Description | Guidance |
|---|---|---|---|---|---|---|
| `snap_tolerance_metres` | `float` | `0.5` | `param` | - | Endpoint coincidence tolerance (m) for degree counting. | Endpoints within this distance are treated as the same node, so shared confluences are not mistaken for dangling tips. |
| `max_gap_length_metres` | `float` | `25.0` | `param` | - | Maximum ROAD-crossing bridge length (m); longer road gaps are rejected. | Roads are narrow, but raise for wide / dual-carriageway crossings or wide verges where a true crossing gap exceeds the default. |
| `max_nonroad_gap_metres` | `float` | `30.0` | `param` | - | Maximum bridge length (m) for gaps NOT associated with a road. | Fills digitising breaks regardless of roads. The STRICT facing/heading gates still apply, so genuine headwaters and parallel streams are not joined. Set 0 to disable non-road bridging. |
| `max_facing_angle_degrees` | `float` | `45.0` | `param` | - | Strict tier: max angle (deg) EACH tip's tangent may deviate from the line to its partner. Used for non-road gaps. | Both tips must point toward each other across the gap. |
| `max_heading_misalignment_degrees` | `float` | `50.0` | `param` | - | Strict tier: max deviation from anti-parallel (deg) between the two tips' tangents. Used for non-road gaps. | A true continuation has the two outward tangents ~opposite (180 deg). Suppresses bridges between parallel unrelated streams. |
| `road_max_facing_angle_degrees` | `float` | `60.0` | `param` | - | Road tier: a road-associated bridge only needs AT LEAST ONE tip to face within this angle (no heading constraint). | The channel kinks at a road crossing so the dangling tip's own tangent is often unreliable; the partner across the road faces well. |
| `tangent_fit_length_metres` | `float` | `15.0` | `param` | - | Inward arc length (m) used to estimate each tip's outward tangent. | Longer values smooth over digitising wiggle near the tip. |
| `road_intersection_tolerance_metres` | `float` | `10.0` | `param` | - | Road-association distance (m): a bridge is 'road-associated' if its connector comes within this distance of a road centreline. | Covers staggered severances where the stub ends sit a few metres off the road. Uses segment-to-segment distance, not just a midpoint test. |
| `require_mutual_best` | `bool` | `true` | `param` | - | Require each tip to be the other's best partner before bridging. | Strongly suppresses spurious bridges; disable only to be more inclusive. |
| `merge_chains` | `bool` | `false` | `param` | - | Dissolve connected chains into single features (NOT yet implemented). | Leave false: keep originals + add bridge lines. Coincident endpoints already make the network topologically continuous downstream. |

<!-- END GENERATED CLI PARAMETERS -->
