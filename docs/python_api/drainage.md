# Drainage And Depressions API

This page is generated from the built Python runtime manifest. Do not edit signatures here.

## `DepressionHierarchyResult`

Complete Barnes depression hierarchy with raster labels and a node table.

Parameters
----------
*args
    Constructor arguments described by the generated binding signature.
**kwargs
    Constructor keyword arguments described by the generated binding
    signature.

Attributes
----------
See the generated binding signature for available properties and methods.

Notes
-----
Complete Barnes depression hierarchy with raster labels and a node table.

Authors
-------
B. Gailleton

### Members

#### `flow_direction`

(self) -> lsdtt3.Raster

#### `labels`

(self) -> lsdtt3.Raster

#### `n_nodes`

(self) -> int

#### `nodes`

(self) -> object

---

## `DrainResult`

No additional binding notes are available.

Parameters
----------
*args
    Constructor arguments described by the generated binding signature.
**kwargs
    Constructor keyword arguments described by the generated binding
    signature.

Attributes
----------
See the generated binding signature for available properties and methods.

Notes
-----
No additional binding notes are available.

Authors
-------
B. Gailleton

### Members

#### `accumulation_curvature`

(self) -> lsdtt3.Raster

#### `blurred_log10_abs_accumulation_curvature`

(self) -> lsdtt3.Raster

#### `drain_mask`

(self) -> lsdtt3.Raster

#### `frangi_best_scale_metres`

(self) -> lsdtt3.Raster

#### `frangi_response`

(self) -> lsdtt3.Raster

#### `log10_abs_accumulation_curvature`

(self) -> lsdtt3.Raster

#### `masked_frangi_scale`

(self) -> lsdtt3.Raster

#### `threshold_mask`

(self) -> lsdtt3.Raster

---

## `InternalDrainageResult`

Closed-depression candidates and, optionally, their GIS-ready footprints.

Parameters
----------
*args
    Constructor arguments described by the generated binding signature.
**kwargs
    Constructor keyword arguments described by the generated binding
    signature.

Attributes
----------
See the generated binding signature for available properties and methods.

Notes
-----
Closed-depression candidates and, optionally, their GIS-ready footprints.

Candidates are point GeoDataFrames with hierarchy, depth, area, and storage
metadata. Footprints are polygons keyed by ``sink_id``.

Authors
-------
B. Gailleton

### Members

#### `candidates`

(self) -> object

#### `n_candidates`

(self) -> int

#### `polygons`

(self) -> object

---

## `STICategoricalRaster`

Rasterized categorical STI field and code lookup.

Parameters
----------
*args
    Constructor arguments described by the generated binding signature.
**kwargs
    Constructor keyword arguments described by the generated binding
    signature.

Attributes
----------
See the generated binding signature for available properties and methods.

Notes
-----
Rasterized categorical STI field and code lookup.

Authors
-------
B. Gailleton

### Members

#### `codes`

(self) -> lsdtt3.Raster

#### `lookup`

(self) -> object

---

## `STILnkdTable`

Parsed soil hydraulic-conductivity lookup table.

Parameters
----------
*args
    Constructor arguments described by the generated binding signature.
**kwargs
    Constructor keyword arguments described by the generated binding
    signature.

Attributes
----------
See the generated binding signature for available properties and methods.

Notes
-----
Parsed soil hydraulic-conductivity lookup table.

Authors
-------
B. Gailleton

### Members

#### `summary`

(self) -> dict

---

## `STISprTable`

Parsed soil-profile runoff lookup table.

Parameters
----------
*args
    Constructor arguments described by the generated binding signature.
**kwargs
    Constructor keyword arguments described by the generated binding
    signature.

Attributes
----------
See the generated binding signature for available properties and methods.

Notes
-----
Parsed soil-profile runoff lookup table.

Authors
-------
B. Gailleton

### Members

#### `summary`

(self) -> dict

---

## `breach_depressions`

breach_depressions(raster: lsdtt3.Raster, n_threads: int = 0) -> lsdtt3.Raster

Parameters
----------
*args
    Positional arguments described by the generated binding signature.
**kwargs
    Keyword arguments described by the generated binding signature.

Returns
-------
object
    Result described by the generated binding signature.

Notes
-----
breach_depressions(raster: lsdtt3.Raster, n_threads: int = 0) -> lsdtt3.Raster

Breach depressions by carving drainage paths.

Authors
-------
B. Gailleton

```python
def breach_depressions(raster: lsdtt3.Raster, n_threads: int = 0) -> lsdtt3.Raster
```

---

## `depression_hierarchy`

depression_hierarchy(dem: lsdtt3.Raster, boundary_conditions: str = 'oooo', ocean_from_boundary: bool = True, nodata_is_ocean: bool = False, ocean_level: object | None = None) -> lsdtt3.DepressionHierarchyResult

Parameters
----------
*args
    Positional arguments described by the generated binding signature.
**kwargs
    Keyword arguments described by the generated binding signature.

Returns
-------
object
    Result described by the generated binding signature.

Notes
-----
depression_hierarchy(dem: lsdtt3.Raster, boundary_conditions: str = 'oooo', ocean_from_boundary: bool = True, nodata_is_ocean: bool = False, ocean_level: object | None = None) -> lsdtt3.DepressionHierarchyResult

Build the complete nested depression hierarchy for a DEM.

Parameters
----------
dem : Raster
    Projected elevation raster.
boundary_conditions : str, optional
    Open edges seed the ocean when ``ocean_from_boundary=True``.
ocean_from_boundary, nodata_is_ocean : bool, optional
    Define ocean seeds from open edges and/or nodata cells.
ocean_level : float, optional
    Treat valid cells at this elevation as ocean seeds.

Returns
-------
DepressionHierarchyResult
    Leaf labels, D8 spill-routing directions, and leaf/meta-depression nodes.

Authors
-------
B. Gailleton

```python
def depression_hierarchy(dem: lsdtt3.Raster, boundary_conditions: str = 'oooo', ocean_from_boundary: bool = True, nodata_is_ocean: bool = False, ocean_level: object | None = None) -> lsdtt3.DepressionHierarchyResult
```

---

## `find_drains`

find_drains(raster: lsdtt3.Raster, sigma_min_metres: float = 1.0, sigma_max_metres: float = 5.0, n_scales: int = 4, response_threshold: float = 0.15, min_component_pixels: int = 6, min_major_axis_pixels: float = 5.0, min_elongation: float = 3.0, accumulation_curvature_lambda: float = 4.0, max_log10_abs_accumulation_curvature: float = -5.0, max_frangi_scale_metres: float = 1.8, n_threads: int = 0) -> lsdtt3.DrainResult

Parameters
----------
*args
    Positional arguments described by the generated binding signature.
**kwargs
    Keyword arguments described by the generated binding signature.

Returns
-------
object
    Result described by the generated binding signature.

Notes
-----
find_drains(raster: lsdtt3.Raster, sigma_min_metres: float = 1.0, sigma_max_metres: float = 5.0, n_scales: int = 4, response_threshold: float = 0.15, min_component_pixels: int = 6, min_major_axis_pixels: float = 5.0, min_elongation: float = 3.0, accumulation_curvature_lambda: float = 4.0, max_log10_abs_accumulation_curvature: float = -5.0, max_frangi_scale_metres: float = 1.8, n_threads: int = 0) -> lsdtt3.DrainResult

Find narrow drainage ditches from Frangi and curvature evidence.

Authors
-------
B. Gailleton

```python
def find_drains(raster: lsdtt3.Raster, sigma_min_metres: float = 1.0, sigma_max_metres: float = 5.0, n_scales: int = 4, response_threshold: float = 0.15, min_component_pixels: int = 6, min_major_axis_pixels: float = 5.0, min_elongation: float = 3.0, accumulation_curvature_lambda: float = 4.0, max_log10_abs_accumulation_curvature: float = -5.0, max_frangi_scale_metres: float = 1.8, n_threads: int = 0) -> lsdtt3.DrainResult
```

---

## `find_internal_drainage`

find_internal_drainage(dem: lsdtt3.Raster, method: str = 'hierarchy', boundary_conditions: str = 'oooo', min_pixels: int = 1, min_depth_metres: float = 0.0, min_area_m2: float = 0.0, min_volume_m3: float = 0.0, leaves_only: bool = True, include_polygons: bool = True, n_threads: int = 0) -> lsdtt3.InternalDrainageResult

Parameters
----------
*args
    Positional arguments described by the generated binding signature.
**kwargs
    Keyword arguments described by the generated binding signature.

Returns
-------
object
    Result described by the generated binding signature.

Notes
-----
find_internal_drainage(dem: lsdtt3.Raster, method: str = 'hierarchy', boundary_conditions: str = 'oooo', min_pixels: int = 1, min_depth_metres: float = 0.0, min_area_m2: float = 0.0, min_volume_m3: float = 0.0, leaves_only: bool = True, include_polygons: bool = True, n_threads: int = 0) -> lsdtt3.InternalDrainageResult

Detect internally drained topographic depressions in a projected DEM.

Parameters
----------
dem : Raster
    Projected elevation raster.
method : {"hierarchy", "fill_difference"}, optional
    Nested Barnes depression hierarchy, or the legacy filled-minus-raw method.
boundary_conditions : str, optional
    Four-character north/east/south/west boundary condition string.
min_pixels : int, optional
    Minimum depression footprint size in raster cells.
 min_depth_metres, min_area_m2, min_volume_m3 : float, optional
    Post-detection filters on depression properties.
leaves_only : bool, optional
    For ``method="hierarchy"``, return only innermost depressions.
include_polygons : bool, optional
    Derive footprint polygons as well as candidate points.
n_threads : int, optional
    OpenMP thread count; zero uses the runtime default.

Returns
-------
InternalDrainageResult
    Point candidates and optional polygon footprints, both in the DEM CRS.

Notes
-----
Use the candidate ``pit_row`` and ``pit_col`` columns with
``fill_pits_with_sinks`` and ``FlowInfo(..., forced_outlets=...)`` when
modelling selected depressions as drainage exits.

Authors
-------
B. Gailleton
S. M. Mudd

**Deprecated compatibility names:**
- `min_area_metres2` -> `min_area_m2`; earliest removal: following minor.
- `min_volume_metres3` -> `min_volume_m3`; earliest removal: following minor.

```python
def find_internal_drainage(dem: lsdtt3.Raster, method: str = 'hierarchy', boundary_conditions: str = 'oooo', min_pixels: int = 1, min_depth_metres: float = 0.0, min_area_m2: float = 0.0, min_volume_m3: float = 0.0, leaves_only: bool = True, include_polygons: bool = True, n_threads: int = 0) -> lsdtt3.InternalDrainageResult
```

---

## `load_sti_lnkd_table`

load_sti_lnkd_table(csv_path: str, depth_cap_cm: float = 100.0, min_valid_depth_cm: float = 30.0, ks_unit_scale_to_m_per_day: float = 1.0, eps: float = 1e-12) -> lsdtt3.STILnkdTable

Parameters
----------
*args
    Positional arguments described by the generated binding signature.
**kwargs
    Keyword arguments described by the generated binding signature.

Returns
-------
object
    Result described by the generated binding signature.

Notes
-----
load_sti_lnkd_table(csv_path: str, depth_cap_cm: float = 100.0, min_valid_depth_cm: float = 30.0, ks_unit_scale_to_m_per_day: float = 1.0, eps: float = 1e-12) -> lsdtt3.STILnkdTable

Load the KSD STI lookup table from a soil-profile CSV.

Authors
-------
B. Gailleton

```python
def load_sti_lnkd_table(csv_path: str, depth_cap_cm: float = 100.0, min_valid_depth_cm: float = 30.0, ks_unit_scale_to_m_per_day: float = 1.0, eps: float = 1e-12) -> lsdtt3.STILnkdTable
```

---

## `load_sti_spr_table`

load_sti_spr_table(csv_path: str, eps: float = 1e-12) -> lsdtt3.STISprTable

Parameters
----------
*args
    Positional arguments described by the generated binding signature.
**kwargs
    Keyword arguments described by the generated binding signature.

Returns
-------
object
    Result described by the generated binding signature.

Notes
-----
load_sti_spr_table(csv_path: str, eps: float = 1e-12) -> lsdtt3.STISprTable

Load the SPR STI lookup table from a CSV.

Authors
-------
B. Gailleton

```python
def load_sti_spr_table(csv_path: str, eps: float = 1e-12) -> lsdtt3.STISprTable
```

---

## `rasterize_sti_categories`

rasterize_sti_categories(reference: lsdtt3.Raster, vector_path: str, field_name: str, layer_name: str, oversample: int = 8) -> lsdtt3.STICategoricalRaster

Parameters
----------
*args
    Positional arguments described by the generated binding signature.
**kwargs
    Keyword arguments described by the generated binding signature.

Returns
-------
object
    Result described by the generated binding signature.

Notes
-----
rasterize_sti_categories(reference: lsdtt3.Raster, vector_path: str, field_name: str, layer_name: str, oversample: int = 8) -> lsdtt3.STICategoricalRaster

Rasterize a string vector attribute and retain its categorical lookup.

Authors
-------
B. Gailleton

```python
def rasterize_sti_categories(reference: lsdtt3.Raster, vector_path: str, field_name: str, layer_name: str, oversample: int = 8) -> lsdtt3.STICategoricalRaster
```

---

## `soil_topographic_index`

soil_topographic_index(twi: lsdtt3.Raster, series_codes: lsdtt3.STICategoricalRaster, land_use: lsdtt3.STICategoricalRaster, lnkd_table: lsdtt3.STILnkdTable, spr_table: lsdtt3.STISprTable, mode: str = 'ksd', spr_transform: str = 'ln_1_minus_spr', eps: float = 1e-12) -> lsdtt3.SoilTopographicIndexResult

Parameters
----------
*args
    Positional arguments described by the generated binding signature.
**kwargs
    Keyword arguments described by the generated binding signature.

Returns
-------
object
    Result described by the generated binding signature.

Notes
-----
soil_topographic_index(twi: lsdtt3.Raster, series_codes: lsdtt3.STICategoricalRaster, land_use: lsdtt3.STICategoricalRaster, lnkd_table: lsdtt3.STILnkdTable, spr_table: lsdtt3.STISprTable, mode: str = 'ksd', spr_transform: str = 'ln_1_minus_spr', eps: float = 1e-12) -> lsdtt3.SoilTopographicIndexResult

Compute STI from TWI, categorical rasters, and lookup tables.

Authors
-------
B. Gailleton

```python
def soil_topographic_index(twi: lsdtt3.Raster, series_codes: lsdtt3.STICategoricalRaster, land_use: lsdtt3.STICategoricalRaster, lnkd_table: lsdtt3.STILnkdTable, spr_table: lsdtt3.STISprTable, mode: str = 'ksd', spr_transform: str = 'ln_1_minus_spr', eps: float = 1e-12) -> lsdtt3.SoilTopographicIndexResult
```
