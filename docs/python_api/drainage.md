# Drainage And Depressions API

This page is generated from the built Python runtime manifest. Do not edit signatures here.

## `DepressionHierarchyResult`

Complete Barnes depression hierarchy with raster labels and a node table.

Returned by ``depression_hierarchy()``; not constructible directly.

Attributes
----------
flow_direction : Raster
labels : Raster
n_nodes : int
nodes : object

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

Drain (ditch) detection output from :func:`find_drains`.

``drain_mask`` is the final 0/1 mask; the other rasters are the Frangi and
accumulation-curvature evidence used to build it.

Returned by ``find_drains()``; not constructible directly.

Attributes
----------
accumulation_curvature : Raster
blurred_log10_abs_accumulation_curvature : Raster
drain_mask : Raster
frangi_best_scale_metres : Raster
frangi_response : Raster
log10_abs_accumulation_curvature : Raster
masked_frangi_scale : Raster
threshold_mask : Raster

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

Candidates are point GeoDataFrames with hierarchy, depth, area, and storage
metadata. Footprints are polygons keyed by ``sink_id``.

Returned by ``find_internal_drainage()``; not constructible directly.

Attributes
----------
candidates : geopandas.GeoDataFrame
    Depression candidates as points with depth, area and storage attributes.
n_candidates : int
polygons : geopandas.GeoDataFrame | None
    Depression footprints keyed by sink_id, or None if not requested.

Authors
-------
B. Gailleton

### Members

#### `candidates`

Depression candidates as points with depth, area and storage attributes.

#### `n_candidates`

(self) -> int

#### `polygons`

Depression footprints keyed by sink_id, or None if not requested.

---

## `STICategoricalRaster`

Rasterized categorical STI field and code lookup.

Returned by ``rasterize_sti_categories()``; not constructible directly.

Attributes
----------
codes : Raster
lookup : object

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

Returned by ``load_sti_lnkd_table()``; not constructible directly.

Attributes
----------
summary : dict

Authors
-------
B. Gailleton

### Members

#### `summary`

(self) -> dict

---

## `STISprTable`

Parsed soil-profile runoff lookup table.

Returned by ``load_sti_spr_table()``; not constructible directly.

Attributes
----------
summary : dict

Authors
-------
B. Gailleton

### Members

#### `summary`

(self) -> dict

---

## `breach_depressions`

Breach depressions by carving drainage paths.

Hybrid breaching (Lindsay 2016): each pit is drained by lowering a path out
of it rather than raising the depression, which changes less of the DEM
than filling in incised terrain. See also ``breach_least_cost``.

Parameters
----------
raster : Raster
    Elevation raster (m).
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
Raster
    Breached DEM (m).

Authors
-------
B. Gailleton

```python
def breach_depressions(raster: lsdtt3.Raster, n_threads: int = 0) -> lsdtt3.Raster
```

---

## `depression_hierarchy`

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

Find narrow drainage ditches from Frangi and curvature evidence.

Frangi valley enhancement is thresholded, small or blob-like components are
removed, and pixels with strong accumulation curvature (natural channels)
or a wide Frangi scale are masked out.

Parameters
----------
raster : Raster
    Elevation raster (m), usually a high-resolution lidar DEM.
sigma_min_metres, sigma_max_metres : float, default=1.0, 5.0
    Frangi smoothing widths (m).
n_scales : int, default=4
    Number of Frangi scales.
response_threshold : float, default=0.15
    Frangi response at or above which a pixel is a candidate.
min_component_pixels : int, default=6
    Smallest connected component kept (pixels).
min_major_axis_pixels : float, default=5.0
    Smallest component major-axis length kept (pixels).
min_elongation : float, default=3.0
    Smallest component elongation (major / minor axis) kept.
accumulation_curvature_lambda : float, default=4.0
    Polyfit length scale (m) for accumulation curvature.
max_log10_abs_accumulation_curvature : float, default=-5.0
    Candidates whose blurred log10 |accumulation curvature| (1/m^2) exceeds
    this are removed.
max_frangi_scale_metres : float, default=1.8
    Candidates whose best Frangi scale (m) exceeds this are removed.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
DrainResult

Authors
-------
B. Gailleton

```python
def find_drains(raster: lsdtt3.Raster, sigma_min_metres: float = 1.0, sigma_max_metres: float = 5.0, n_scales: int = 4, response_threshold: float = 0.15, min_component_pixels: int = 6, min_major_axis_pixels: float = 5.0, min_elongation: float = 3.0, accumulation_curvature_lambda: float = 4.0, max_log10_abs_accumulation_curvature: float = -5.0, max_frangi_scale_metres: float = 1.8, n_threads: int = 0) -> lsdtt3.DrainResult
```

---

## `find_internal_drainage`

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

Load the KSD STI lookup table from a soil-profile CSV.

Parameters
----------
csv_path : str
    Soil-profile CSV with saturated hydraulic conductivity by horizon.
depth_cap_cm : float, default=100.0
    Depth (cm) to which profiles are integrated.
min_valid_depth_cm : float, default=30.0
    Profiles shallower than this (cm) are rejected.
ks_unit_scale_to_m_per_day : float, default=1.0
    Factor converting the CSV conductivity to m/day.
eps : float, default=1e-12
    Small positive floor that avoids log(0).

Returns
-------
STILnkdTable

Authors
-------
B. Gailleton

```python
def load_sti_lnkd_table(csv_path: str, depth_cap_cm: float = 100.0, min_valid_depth_cm: float = 30.0, ks_unit_scale_to_m_per_day: float = 1.0, eps: float = 1e-12) -> lsdtt3.STILnkdTable
```

---

## `load_sti_spr_table`

Load the SPR STI lookup table from a CSV.

Parameters
----------
csv_path : str
    CSV of standard percentage runoff (SPR) values.
eps : float, default=1e-12
    Small positive floor that avoids log(0).

Returns
-------
STISprTable

Authors
-------
B. Gailleton

```python
def load_sti_spr_table(csv_path: str, eps: float = 1e-12) -> lsdtt3.STISprTable
```

---

## `rasterize_sti_categories`

Rasterize a string vector attribute and retain its categorical lookup.

Parameters
----------
reference : Raster
    Raster whose grid is used.
vector_path : str
    Polygon vector file.
field_name : str
    Attribute holding the category (for example a soil series code).
layer_name : str
    Layer in the vector file.
oversample : int, default=8
    Sub-pixel oversampling used to pick each pixel's majority category.

Returns
-------
STICategoricalRaster
    Integer codes plus the code-to-category lookup.

Authors
-------
B. Gailleton

```python
def rasterize_sti_categories(reference: lsdtt3.Raster, vector_path: str, field_name: str, layer_name: str, oversample: int = 8) -> lsdtt3.STICategoricalRaster
```

---

## `soil_topographic_index`

Compute STI from TWI, categorical rasters, and lookup tables.

Parameters
----------
twi : Raster
    Topographic wetness index (``topographic_wetness_index``).
series_codes : STICategoricalRaster
    Soil series categories (``rasterize_sti_categories``).
land_use : STICategoricalRaster
    Land-use categories on the same grid.
lnkd_table : STILnkdTable
    Conductivity lookup (``load_sti_lnkd_table``).
spr_table : STISprTable
    Runoff lookup (``load_sti_spr_table``).
mode : {"ksd", "spr"}, default="ksd"
    Soil term from conductivity x depth or from standard percentage runoff.
spr_transform : {"ln_1_minus_spr", "ln_spr"}, default="ln_1_minus_spr"
    Transform of the SPR term.
eps : float, default=1e-12
    Small positive floor that avoids log(0).

Returns
-------
SoilTopographicIndexResult

Authors
-------
B. Gailleton

```python
def soil_topographic_index(twi: lsdtt3.Raster, series_codes: lsdtt3.STICategoricalRaster, land_use: lsdtt3.STICategoricalRaster, lnkd_table: lsdtt3.STILnkdTable, spr_table: lsdtt3.STISprTable, mode: str = 'ksd', spr_transform: str = 'ln_1_minus_spr', eps: float = 1e-12) -> lsdtt3.SoilTopographicIndexResult
```
