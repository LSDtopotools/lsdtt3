# Raster And Terrain API

This page is generated from the built Python runtime manifest. Do not edit signatures here.

## `FlatSurfaceClassification`

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

#### `floodplain_mask`

(self) -> lsdtt3.Raster

#### `global_breaks_metres`

(self) -> list[float]

#### `hand_band`

(self) -> lsdtt3.Raster

#### `n_floodplain_labels`

(self) -> int

#### `n_terrace_labels`

(self) -> int

#### `surface_class`

(self) -> lsdtt3.Raster

#### `table`

(self) -> object

---

## `FlatSurfaceResult`

Results from texture- and relief-based flat-surface detection.

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
Results from texture- and relief-based flat-surface detection.

Authors
-------
B. Gailleton

### Members

#### `combined_flat_mask`

(self) -> lsdtt3.Raster

#### `healed_flat_mask`

(self) -> lsdtt3.Raster

#### `iqr_flat_mask`

(self) -> lsdtt3.Raster

#### `iqr_residuals`

(self) -> lsdtt3.Raster

#### `iqr_texture`

(self) -> lsdtt3.Raster

#### `labelled_flat_mask`

(self) -> lsdtt3.Raster

#### `madk2_flat_mask`

(self) -> lsdtt3.Raster

#### `madk2_roughness`

(self) -> lsdtt3.Raster

#### `nvd_flat_mask`

(self) -> lsdtt3.Raster

#### `nvd_s3`

(self) -> lsdtt3.Raster

#### `relief`

(self) -> lsdtt3.Raster

#### `relief_flat_mask`

(self) -> lsdtt3.Raster

---

## `FrangiResult`

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

#### `best_scale_metres`

(self) -> lsdtt3.Raster

#### `response`

(self) -> lsdtt3.Raster

---

## `Raster`

An in-memory float raster with an explicit validity mask.

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
    An in-memory float raster with an explicit validity mask.

Instances are value objects: functions in this module return new rasters and
do not modify their inputs.

    Authors
    -------
    B. Gailleton


### Members

#### `from_numpy`

from_numpy(data: object, mask: object, metadata: lsdtt3.RasterMeta) -> lsdtt3.Raster

Create a raster from copied two-dimensional NumPy-compatible arrays.

Parameters
----------
data : array-like of float
    Values with shape ``metadata.n_rows, metadata.n_cols``.
mask : array-like of bool
    True for valid pixels. False pixels are stored as ``nodata_value``.
metadata : RasterMeta
    Grid metadata matching both array dimensions.

```python
def from_numpy(data: object, mask: object, metadata: lsdtt3.RasterMeta) -> lsdtt3.Raster
```

#### `mask_numpy`

mask_numpy(self) -> object

Return a copied Boolean array, True where raster values are valid.

```python
def mask_numpy(self) -> object
```

#### `metadata`

(self) -> lsdtt3.RasterMeta

#### `read`

read(path: str) -> lsdtt3.Raster

Read a single-band GDAL-compatible raster.

Parameters
----------
path : path-like
    Raster file to load.

Returns
-------
Raster
    Raster values, metadata, and nodata mask.

```python
def read(path: str) -> lsdtt3.Raster
```

#### `shape`

(self) -> tuple

#### `to_numpy`

to_numpy(self) -> object

Return a copied ``float32`` array of raster values.

Invalid cells contain the raster nodata value; use ``mask_numpy`` to identify
them.

```python
def to_numpy(self) -> object
```

#### `valid_pixel_count`

(self) -> int

#### `write`

write(self, path: str, format: str = 'cog') -> None

Write this raster using an LSDTopoTools GDAL output format.

Parameters
----------
path : path-like
    Destination file.
format : str, default="cog"
    Output format, for example ``"cog"`` or ``"gtiff"``.

```python
def write(self, path: str, format: str = 'cog') -> None
```

---

## `RasterHierarchy`

Resolution hierarchy for one in-memory float Raster.

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
    Resolution hierarchy for one in-memory float Raster.

Level 0 is finest. Standard levels are exact factor-two grids; explicit
integer factors can be added over the same extent.
The hierarchy stores no solver state and each call returning a Raster returns
a copy, consistent with the rest of the Python Raster API.

    Authors
    -------
    B. Gailleton


```python
def __init__(self, finest: lsdtt3.Raster, storage: str = 'retain_allocated', keep_finest: bool = True, levels: int = 0) -> None
```

### Members

#### `__init__`

__init__(self, finest: lsdtt3.Raster, storage: str = 'retain_allocated', keep_finest: bool = True, levels: int = 0) -> None

```python
def __init__(self, finest: lsdtt3.Raster, storage: str = 'retain_allocated', keep_finest: bool = True, levels: int = 0) -> None
```

#### `active_level`

(self) -> int

#### `add_level_for_factor`

add_level_for_factor(self, factor: int) -> int

```python
def add_level_for_factor(self, factor: int) -> int
```

#### `allocate_to`

allocate_to(self, level: int) -> None

```python
def allocate_to(self, level: int) -> None
```

#### `clear`

clear(self) -> None

```python
def clear(self) -> None
```

#### `clear_except`

clear_except(self, level: int) -> None

```python
def clear_except(self, level: int) -> None
```

#### `has_level`

has_level(self, level: int) -> bool

```python
def has_level(self, level: int) -> bool
```

#### `level`

level(self, level: int) -> lsdtt3.Raster

```python
def level(self, level: int) -> lsdtt3.Raster
```

#### `level_factor`

level_factor(self, level: int) -> int

Integer coarsening factor relative to level zero.

```python
def level_factor(self, level: int) -> int
```

#### `level_shape`

level_shape(self, level: int) -> tuple

Grid shape as (n_rows, n_cols).

```python
def level_shape(self, level: int) -> tuple
```

#### `move_to`

move_to(self, level: int, method: str) -> lsdtt3.Raster

```python
def move_to(self, level: int, method: str) -> lsdtt3.Raster
```

#### `n_levels`

(self) -> int

#### `prolong_to`

prolong_to(self, level: int, method: str = 'linear') -> lsdtt3.Raster

```python
def prolong_to(self, level: int, method: str = 'linear') -> lsdtt3.Raster
```

#### `restrict_to`

restrict_to(self, level: int, method: str = 'conservative_mean') -> lsdtt3.Raster

```python
def restrict_to(self, level: int, method: str = 'conservative_mean') -> lsdtt3.Raster
```

#### `set_active`

set_active(self, raster: lsdtt3.Raster) -> None

```python
def set_active(self, raster: lsdtt3.Raster) -> None
```

---

## `RasterMeta`

Metadata describing a north-up, square-cell raster grid.

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
    Metadata describing a north-up, square-cell raster grid.

Use ``RasterMeta.from_shape`` to construct an instance.

    Authors
    -------
    B. Gailleton


### Members

#### `cell_size`

(self) -> float

#### `crs_wkt`

(self) -> str

#### `epsg_code`

(self) -> int

#### `from_shape`

from_shape(n_rows: int, n_cols: int, x_min: float = 0.0, y_min: float = 0.0, cell_size: float = 1.0, nodata_value: float = -9999.0, crs_wkt: str = '', epsg_code: int = 0) -> lsdtt3.RasterMeta

Create raster metadata from a grid shape and lower-left origin.

Parameters
----------
n_rows, n_cols : int
    Positive grid dimensions.
x_min, y_min, cell_size : float
    Lower-left origin and square cell size.
nodata_value : float, optional
    Value stored for invalid pixels.
crs_wkt : str, optional
    Coordinate reference system in WKT form.
epsg_code : int, optional
    EPSG identifier when known.

```python
def from_shape(n_rows: int, n_cols: int, x_min: float = 0.0, y_min: float = 0.0, cell_size: float = 1.0, nodata_value: float = -9999.0, crs_wkt: str = '', epsg_code: int = 0) -> lsdtt3.RasterMeta
```

#### `is_valid`

(self) -> bool

#### `n_cols`

(self) -> int

#### `n_pixels`

(self) -> int

#### `n_rows`

(self) -> int

#### `nodata_value`

(self) -> float

#### `pixel_area`

(self) -> float

#### `x_max`

(self) -> float

#### `x_min`

(self) -> float

#### `y_max`

(self) -> float

#### `y_min`

(self) -> float

---

## `SoilTopographicIndexResult`

Soil-topographic-index outputs and diagnostics.

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
Soil-topographic-index outputs and diagnostics.

Authors
-------
B. Gailleton

### Members

#### `lookup`

(self) -> object

#### `stats`

(self) -> dict

#### `sti`

(self) -> lsdtt3.Raster

#### `term`

(self) -> lsdtt3.Raster

#### `valid_mask`

(self) -> lsdtt3.Raster

---

## `SpectralResult`

Results from a full DEM spectral analysis.

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
    Results from a full DEM spectral analysis.

Spectral arrays are returned as copied NumPy arrays. Frequencies are in
cycles per metre and wavelengths are in metres.

    Authors
    -------
    B. Gailleton


### Members

#### `background_psd`

(self) -> object

#### `background_psd_2d`

(self) -> object

#### `binned`

(self) -> dict

#### `frequency`

(self) -> object

#### `frequency_x`

(self) -> object

#### `frequency_y`

(self) -> object

#### `normalised_psd`

(self) -> object

#### `padded_shape`

(self) -> tuple

#### `psd`

(self) -> object

#### `psd_2d`

(self) -> object

#### `rollover_beta`

(self) -> float

#### `rollover_frequency`

(self) -> float

#### `rollover_wavelength`

(self) -> float

#### `sub_rollover_beta`

(self) -> float

#### `wavelength`

(self) -> object

---

## `SwathProfile`

Results from a corridor-based longitudinal swath profile.

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
    Results from a corridor-based longitudinal swath profile.

The distance rasters are Raster objects. ``table`` is a copied pandas
DataFrame with one row per along-axis bin.

    Authors
    -------
    B. Gailleton


### Members

#### `along_axis`

(self) -> lsdtt3.Raster

#### `perpendicular_distance`

(self) -> lsdtt3.Raster

#### `quantile_levels`

(self) -> object

#### `signed_perpendicular_distance`

(self) -> lsdtt3.Raster

#### `table`

(self) -> object

---

## `add`

add(a: lsdtt3.Raster, b: lsdtt3.Raster) -> lsdtt3.Raster

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
    add(a: lsdtt3.Raster, b: lsdtt3.Raster) -> lsdtt3.Raster

Add two grid-compatible rasters.

    Authors
    -------
    B. Gailleton


```python
def add(a: lsdtt3.Raster, b: lsdtt3.Raster) -> lsdtt3.Raster
```

---

## `add_uniform_noise`

add_uniform_noise(raster: lsdtt3.Raster, seed: int, magnitude: float) -> lsdtt3.Raster

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
    add_uniform_noise(raster: lsdtt3.Raster, seed: int, magnitude: float) -> lsdtt3.Raster

Add reproducible uniform noise to valid pixels.

    Authors
    -------
    B. Gailleton


```python
def add_uniform_noise(raster: lsdtt3.Raster, seed: int, magnitude: float) -> lsdtt3.Raster
```

---

## `aspect`

aspect(raster: lsdtt3.Raster, lambda_metres: float = 0.0, n_threads: int = 0, flat_slope_tolerance: float = 1e-12) -> lsdtt3.Raster

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
    aspect(raster: lsdtt3.Raster, lambda_metres: float = 0.0, n_threads: int = 0, flat_slope_tolerance: float = 1e-12) -> lsdtt3.Raster

Compute polyfit aspect in degrees; lambda_metres defaults to three cells.

    Authors
    -------
    B. Gailleton


```python
def aspect(raster: lsdtt3.Raster, lambda_metres: float = 0.0, n_threads: int = 0, flat_slope_tolerance: float = 1e-12) -> lsdtt3.Raster
```

---

## `bandpass_filter`

bandpass_filter(raster: lsdtt3.Raster, short_metres: float, long_metres: float, window: str = 'hann', clip_to_valid_rectangle: bool = False, n_threads: int = 0) -> lsdtt3.Raster

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
    bandpass_filter(raster: lsdtt3.Raster, short_metres: float, long_metres: float, window: str = 'hann', clip_to_valid_rectangle: bool = False, n_threads: int = 0) -> lsdtt3.Raster

Return the DEM component between two wavelengths.

    Authors
    -------
    B. Gailleton


```python
def bandpass_filter(raster: lsdtt3.Raster, short_metres: float, long_metres: float, window: str = 'hann', clip_to_valid_rectangle: bool = False, n_threads: int = 0) -> lsdtt3.Raster
```

---

## `breach_least_cost`

breach_least_cost(raster: lsdtt3.Raster, max_distance: int = 1000, max_cost: float = 3.4028234663852886e+38, prefer_shortest_path: bool = False, flat_increment: float = 0.0, fill_remainder: bool = True, n_threads: int = 0, edge_policy: str = 'outlet', fill_min_slope: float = 9.999999747378752e-05) -> lsdtt3.Raster

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
    breach_least_cost(raster: lsdtt3.Raster, max_distance: int = 1000, max_cost: float = 3.4028234663852886e+38, prefer_shortest_path: bool = False, flat_increment: float = 0.0, fill_remainder: bool = True, n_threads: int = 0, edge_policy: str = 'outlet', fill_min_slope: float = 9.999999747378752e-05) -> lsdtt3.Raster

Breach depressions along least-cost paths (Lindsay 2020).

A port of WhiteboxTools BreachDepressionsLeastCost (whitebox_workflows 2.0.6).

Parameters
----------
max_distance : int, default=1000
    Maximum breach path length in cells, counted in steps from the pit (a
    diagonal step counts as one). 1 to 32767.
max_cost : float, optional
    A pit whose cheapest path costs more is left unbreached. The cost is the
    depth each path cell is lowered, summed: metres, or metres x cells when
    ``prefer_shortest_path`` is true. Default unlimited.
prefer_shortest_path : bool, default=False
    Weight each step's cost by its length in cells, favouring shorter paths.
flat_increment : float, default=0
    Elevation drop (m) per cell along a breach channel. 0 or less selects twice
    the float32 spacing at the largest absolute elevation.
fill_remainder : bool, default=True
    Fill what the breach leaves (Priority-Flood with ``fill_min_slope``) so
    every interior cell has a lower neighbour.
n_threads : int, default=0
    Ignored; the algorithm is serial.
edge_policy : {"outlet", "whitebox"}, default="outlet"
    ``outlet``: the raster edge is an outlet, like a nodata cell.
    ``whitebox``: as in WhiteboxTools, a breach must reach an edge cell that
    is already low enough.
fill_min_slope : float, default=0.0001
    Gradient (m/m) imposed by the remainder fill.

    Authors
    -------
    B. Gailleton


```python
def breach_least_cost(raster: lsdtt3.Raster, max_distance: int = 1000, max_cost: float = 3.4028234663852886e+38, prefer_shortest_path: bool = False, flat_increment: float = 0.0, fill_remainder: bool = True, n_threads: int = 0, edge_policy: str = 'outlet', fill_min_slope: float = 9.999999747378752e-05) -> lsdtt3.Raster
```

---

## `carve_depth`

carve_depth(raster: lsdtt3.Raster, n_threads: int = 0) -> lsdtt3.Raster

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
    carve_depth(raster: lsdtt3.Raster, n_threads: int = 0) -> lsdtt3.Raster

Return the non-negative elevation removed by depression breaching.

    Authors
    -------
    B. Gailleton


```python
def carve_depth(raster: lsdtt3.Raster, n_threads: int = 0) -> lsdtt3.Raster
```

---

## `classify_flat_surfaces`

classify_flat_surfaces(labelled_flat_mask: lsdtt3.Raster, height_above_channel: lsdtt3.Raster, channel_source_ids: lsdtt3.Raster, floodplain_max_height_metres: float = 1.0, terrace_min_height_metres: float = 0.0, channel_adjacency_radius_metres: float = 0.0, terrace_max_distance_to_channel_metres: float = 0.0, enable_leak_split: bool = True, enable_hand_banding: bool = False, donor_flat: object | None = None, channel_flow_distance: object | None = None, n_threads: int = 0) -> lsdtt3.FlatSurfaceClassification

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
    classify_flat_surfaces(labelled_flat_mask: lsdtt3.Raster, height_above_channel: lsdtt3.Raster, channel_source_ids: lsdtt3.Raster, floodplain_max_height_metres: float = 1.0, terrace_min_height_metres: float = 0.0, channel_adjacency_radius_metres: float = 0.0, terrace_max_distance_to_channel_metres: float = 0.0, enable_leak_split: bool = True, enable_hand_banding: bool = False, donor_flat: object | None = None, channel_flow_distance: object | None = None, n_threads: int = 0) -> lsdtt3.FlatSurfaceClassification

Classify labelled flats as floodplain or terrace from HAND and channel IDs.

``channel_source_ids`` is a co-registered Raster containing non-negative
source IDs on in-scope channel pixels and negative values elsewhere.

    Authors
    -------
    B. Gailleton


```python
def classify_flat_surfaces(labelled_flat_mask: lsdtt3.Raster, height_above_channel: lsdtt3.Raster, channel_source_ids: lsdtt3.Raster, floodplain_max_height_metres: float = 1.0, terrace_min_height_metres: float = 0.0, channel_adjacency_radius_metres: float = 0.0, terrace_max_distance_to_channel_metres: float = 0.0, enable_leak_split: bool = True, enable_hand_banding: bool = False, donor_flat: object | None = None, channel_flow_distance: object | None = None, n_threads: int = 0) -> lsdtt3.FlatSurfaceClassification
```

---

## `classify_flat_surfaces_with_scope`

classify_flat_surfaces_with_scope(labelled_flat_mask: lsdtt3.Raster, dem: lsdtt3.Raster, scope: lsdtt3.ChannelScope, floodplain_max_height_metres: float = 1.0, terrace_min_height_metres: float = 0.0, channel_adjacency_radius_metres: float = 0.0, terrace_max_distance_to_channel_metres: float = 0.0, enable_leak_split: bool = True, enable_hand_banding: bool = False, n_threads: int = 0) -> lsdtt3.FlatSurfaceClassification

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
    classify_flat_surfaces_with_scope(labelled_flat_mask: lsdtt3.Raster, dem: lsdtt3.Raster, scope: lsdtt3.ChannelScope, floodplain_max_height_metres: float = 1.0, terrace_min_height_metres: float = 0.0, channel_adjacency_radius_metres: float = 0.0, terrace_max_distance_to_channel_metres: float = 0.0, enable_leak_split: bool = True, enable_hand_banding: bool = False, n_threads: int = 0) -> lsdtt3.FlatSurfaceClassification

Classify flat surfaces using a ``ChannelScope``.

HAND, its nearest-channel donors, flow distance, channel mask, and source IDs
are derived from the scope's retained D8 flow topology.

    Authors
    -------
    B. Gailleton


```python
def classify_flat_surfaces_with_scope(labelled_flat_mask: lsdtt3.Raster, dem: lsdtt3.Raster, scope: lsdtt3.ChannelScope, floodplain_max_height_metres: float = 1.0, terrace_min_height_metres: float = 0.0, channel_adjacency_radius_metres: float = 0.0, terrace_max_distance_to_channel_metres: float = 0.0, enable_leak_split: bool = True, enable_hand_banding: bool = False, n_threads: int = 0) -> lsdtt3.FlatSurfaceClassification
```

---

## `classify_surface_shape`

classify_surface_shape(raster: lsdtt3.Raster, lambda_metres: float = 0.0, z_factor: float = 1.0, n_threads: int = 0, kt: float = 0.0) -> dict

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
    classify_surface_shape(raster: lsdtt3.Raster, lambda_metres: float = 0.0, z_factor: float = 1.0, n_threads: int = 0, kt: float = 0.0) -> dict
classify_surface_shape(maximum_curvature: lsdtt3.Raster, minimum_curvature: lsdtt3.Raster, kt: float = 0.0) -> dict

Overloaded function.

1. ``classify_surface_shape(raster: lsdtt3.Raster, lambda_metres: float = 0.0, z_factor: float = 1.0, n_threads: int = 0, kt: float = 0.0) -> dict``

Compute polyfit principal curvatures, then classify surface shape.

2. ``classify_surface_shape(maximum_curvature: lsdtt3.Raster, minimum_curvature: lsdtt3.Raster, kt: float = 0.0) -> dict``

Classify surface shape from existing principal-curvature rasters.

    Authors
    -------
    B. Gailleton


```python
def classify_surface_shape(raster: lsdtt3.Raster, lambda_metres: float = 0.0, z_factor: float = 1.0, n_threads: int = 0, kt: float = 0.0) -> dict
def classify_surface_shape(maximum_curvature: lsdtt3.Raster, minimum_curvature: lsdtt3.Raster, kt: float = 0.0) -> dict
```

---

## `clip_to_maximum_rectangle`

clip_to_maximum_rectangle(raster: lsdtt3.Raster) -> lsdtt3.Raster

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
    clip_to_maximum_rectangle(raster: lsdtt3.Raster) -> lsdtt3.Raster

Crop to the largest all-valid axis-aligned rectangle.

    Authors
    -------
    B. Gailleton


```python
def clip_to_maximum_rectangle(raster: lsdtt3.Raster) -> lsdtt3.Raster
```

---

## `clip_to_polygons`

clip_to_polygons(raster: lsdtt3.Raster, polygons: object, coverage: str = 'majority', crop: bool = True) -> lsdtt3.Raster

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
    clip_to_polygons(raster: lsdtt3.Raster, polygons: object, coverage: str = 'majority', crop: bool = True) -> lsdtt3.Raster

Crop to polygon bounds and set pixels outside polygons to nodata.

Parameters
----------
polygons : geopandas.GeoDataFrame
    Polygon geometries with a CRS exactly matching raster.
coverage : {"any", "majority", "full"}, default="majority"
    Pixel-inclusion rule.
crop : bool, default=True
    Crop to polygon bounds before masking. False retains the source extent.

    Authors
    -------
    B. Gailleton


```python
def clip_to_polygons(raster: lsdtt3.Raster, polygons: object, coverage: str = 'majority', crop: bool = True) -> lsdtt3.Raster
```

---

## `detect_directional_periodicity`

detect_directional_periodicity(raster: lsdtt3.Raster, window_size_pixels: int, step_pixels: int, min_spacing_metres: float, max_spacing_metres: float, min_valid_fraction: float = 0.6, n_threads: int = 0) -> dict

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
    detect_directional_periodicity(raster: lsdtt3.Raster, window_size_pixels: int, step_pixels: int, min_spacing_metres: float, max_spacing_metres: float, min_valid_fraction: float = 0.6, n_threads: int = 0) -> dict

Detect local regularly spaced parallel features with sliding-window FFTs.

    Authors
    -------
    B. Gailleton


```python
def detect_directional_periodicity(raster: lsdtt3.Raster, window_size_pixels: int, step_pixels: int, min_spacing_metres: float, max_spacing_metres: float, min_valid_fraction: float = 0.6, n_threads: int = 0) -> dict
```

---

## `double_raster`

double_raster(source: lsdtt3.Raster, method: str = 'linear') -> lsdtt3.Raster

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
    double_raster(source: lsdtt3.Raster, method: str = 'linear') -> lsdtt3.Raster

Double both raster dimensions over the same extent using the optimized path.

    Authors
    -------
    B. Gailleton


```python
def double_raster(source: lsdtt3.Raster, method: str = 'linear') -> lsdtt3.Raster
```

---

## `fill_depth`

fill_depth(raster: lsdtt3.Raster, min_slope: float = 9.999999747378752e-05, n_threads: int = 0) -> lsdtt3.Raster

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
    fill_depth(raster: lsdtt3.Raster, min_slope: float = 9.999999747378752e-05, n_threads: int = 0) -> lsdtt3.Raster

Return the non-negative elevation added by pit filling.

    Authors
    -------
    B. Gailleton


```python
def fill_depth(raster: lsdtt3.Raster, min_slope: float = 9.999999747378752e-05, n_threads: int = 0) -> lsdtt3.Raster
```

---

## `fill_interior_nodata`

fill_interior_nodata(raster: lsdtt3.Raster, method: str = 'laplacian', max_iter: int = 10000, tolerance: float = 9.999999974752427e-07, n_threads: int = 0) -> lsdtt3.Raster

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
    fill_interior_nodata(raster: lsdtt3.Raster, method: str = 'laplacian', max_iter: int = 10000, tolerance: float = 9.999999974752427e-07, n_threads: int = 0) -> lsdtt3.Raster

Fill nodata regions that are not connected to the raster edge.

    Authors
    -------
    B. Gailleton


```python
def fill_interior_nodata(raster: lsdtt3.Raster, method: str = 'laplacian', max_iter: int = 10000, tolerance: float = 9.999999974752427e-07, n_threads: int = 0) -> lsdtt3.Raster
```

---

## `fill_pits`

fill_pits(raster: lsdtt3.Raster, min_slope: float = 9.999999747378752e-05, n_threads: int = 0, priority_flood_backend: str = 'binary') -> lsdtt3.Raster

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
    fill_pits(raster: lsdtt3.Raster, min_slope: float = 9.999999747378752e-05, n_threads: int = 0, priority_flood_backend: str = 'binary') -> lsdtt3.Raster

Fill depressions while imposing an optional minimum drainage slope.

    Authors
    -------
    B. Gailleton


```python
def fill_pits(raster: lsdtt3.Raster, min_slope: float = 9.999999747378752e-05, n_threads: int = 0, priority_flood_backend: str = 'binary') -> lsdtt3.Raster
```

---

## `fill_pits_boundary_aware`

fill_pits_boundary_aware(raster: lsdtt3.Raster, boundary_conditions: str = 'oooo', min_slope: float = 9.999999747378752e-05, n_threads: int = 0, priority_flood_backend: str = 'binary') -> lsdtt3.Raster

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
    fill_pits_boundary_aware(raster: lsdtt3.Raster, boundary_conditions: str = 'oooo', min_slope: float = 9.999999747378752e-05, n_threads: int = 0, priority_flood_backend: str = 'binary') -> lsdtt3.Raster

Fill pits with NESW open, closed, periodic, or base-level boundaries.

    Authors
    -------
    B. Gailleton


```python
def fill_pits_boundary_aware(raster: lsdtt3.Raster, boundary_conditions: str = 'oooo', min_slope: float = 9.999999747378752e-05, n_threads: int = 0, priority_flood_backend: str = 'binary') -> lsdtt3.Raster
```

---

## `fill_pits_priority_flood`

fill_pits_priority_flood(raster: lsdtt3.Raster, min_slope: float = 9.999999747378752e-05, n_threads: int = 0, use_radix_heap: bool = False) -> lsdtt3.Raster

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
    fill_pits_priority_flood(raster: lsdtt3.Raster, min_slope: float = 9.999999747378752e-05, n_threads: int = 0, use_radix_heap: bool = False) -> lsdtt3.Raster

Explicit Improved Priority-Flood; optionally use its radix-heap queue.

    Authors
    -------
    B. Gailleton


```python
def fill_pits_priority_flood(raster: lsdtt3.Raster, min_slope: float = 9.999999747378752e-05, n_threads: int = 0, use_radix_heap: bool = False) -> lsdtt3.Raster
```

---

## `fill_pits_with_outlets`

fill_pits_with_outlets(raster: lsdtt3.Raster, outlet_mask: lsdtt3.Raster, min_slope: float = 9.999999747378752e-05, n_threads: int = 0, priority_flood_backend: str = 'binary') -> lsdtt3.Raster

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
    fill_pits_with_outlets(raster: lsdtt3.Raster, outlet_mask: lsdtt3.Raster, min_slope: float = 9.999999747378752e-05, n_threads: int = 0, priority_flood_backend: str = 'binary') -> lsdtt3.Raster

Fill pits using positive cells in outlet_mask as drainage outlets.

    Authors
    -------
    B. Gailleton


```python
def fill_pits_with_outlets(raster: lsdtt3.Raster, outlet_mask: lsdtt3.Raster, min_slope: float = 9.999999747378752e-05, n_threads: int = 0, priority_flood_backend: str = 'binary') -> lsdtt3.Raster
```

---

## `fill_pits_with_sinks`

fill_pits_with_sinks(raster: lsdtt3.Raster, sinks: lsdtt3.Raster, min_slope: float = 9.999999747378752e-05, n_threads: int = 0, priority_flood_backend: str = 'binary') -> lsdtt3.Raster

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
    fill_pits_with_sinks(raster: lsdtt3.Raster, sinks: lsdtt3.Raster, min_slope: float = 9.999999747378752e-05, n_threads: int = 0, priority_flood_backend: str = 'binary') -> lsdtt3.Raster
fill_pits_with_sinks(raster: lsdtt3.Raster, sinks: object, min_slope: float = 9.999999747378752e-05, n_threads: int = 0, priority_flood_backend: str = 'binary') -> lsdtt3.Raster

Overloaded function.

1. ``fill_pits_with_sinks(raster: lsdtt3.Raster, sinks: lsdtt3.Raster, min_slope: float = 9.999999747378752e-05, n_threads: int = 0, priority_flood_backend: str = 'binary') -> lsdtt3.Raster``

Fill pits while retaining positive cells in sinks as declared sinks.

2. ``fill_pits_with_sinks(raster: lsdtt3.Raster, sinks: object, min_slope: float = 9.999999747378752e-05, n_threads: int = 0, priority_flood_backend: str = 'binary') -> lsdtt3.Raster``

Fill pits while retaining zero-based (row, col) sinks.

    Authors
    -------
    B. Gailleton


```python
def fill_pits_with_sinks(raster: lsdtt3.Raster, sinks: lsdtt3.Raster, min_slope: float = 9.999999747378752e-05, n_threads: int = 0, priority_flood_backend: str = 'binary') -> lsdtt3.Raster
def fill_pits_with_sinks(raster: lsdtt3.Raster, sinks: object, min_slope: float = 9.999999747378752e-05, n_threads: int = 0, priority_flood_backend: str = 'binary') -> lsdtt3.Raster
```

---

## `find_flat_surfaces`

find_flat_surfaces(dem: lsdtt3.Raster, use_madk2: bool = True, madk2_window_radius_metres: float = 100.0, madk2_threshold_percentile: float = 25.0, use_residual_iqr: bool = True, iqr_method: str = 'fftw_highpass', iqr_fit_lambda_metres: float = 500.0, iqr_detrend_radius_metres: float = 500.0, iqr_window_radius_metres: float = 50.0, iqr_window_shape: str = 'circular', iqr_threshold_percentile: float = 25.0, use_nvd: bool = False, nvd_polyfit_lambda: float = 10.0, nvd_window_radius_metres: float = 100.0, nvd_window_shape: str = 'circular', nvd_threshold_percentile: float = 25.0, require_all_metrics_pass: bool = True, directional_close_line_length_pixels: int = 0, directional_close_n_angles: int = 18, min_component_pixels: int = 100, relief_window_side_pixels: int = 5, relief_threshold_metres: float = 1.0, remove_risers: bool = False, n_threads: int = 0) -> lsdtt3.FlatSurfaceResult

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
    find_flat_surfaces(dem: lsdtt3.Raster, use_madk2: bool = True, madk2_window_radius_metres: float = 100.0, madk2_threshold_percentile: float = 25.0, use_residual_iqr: bool = True, iqr_method: str = 'fftw_highpass', iqr_fit_lambda_metres: float = 500.0, iqr_detrend_radius_metres: float = 500.0, iqr_window_radius_metres: float = 50.0, iqr_window_shape: str = 'circular', iqr_threshold_percentile: float = 25.0, use_nvd: bool = False, nvd_polyfit_lambda: float = 10.0, nvd_window_radius_metres: float = 100.0, nvd_window_shape: str = 'circular', nvd_threshold_percentile: float = 25.0, require_all_metrics_pass: bool = True, directional_close_line_length_pixels: int = 0, directional_close_n_angles: int = 18, min_component_pixels: int = 100, relief_window_side_pixels: int = 5, relief_threshold_metres: float = 1.0, remove_risers: bool = False, n_threads: int = 0) -> lsdtt3.FlatSurfaceResult

Detect and label flat surfaces from texture and local relief.

    Authors
    -------
    B. Gailleton


**Deprecated compatibility names:**
- `relief_window_pixels` -> `relief_window_side_pixels`; earliest removal: following minor.

```python
def find_flat_surfaces(dem: lsdtt3.Raster, use_madk2: bool = True, madk2_window_radius_metres: float = 100.0, madk2_threshold_percentile: float = 25.0, use_residual_iqr: bool = True, iqr_method: str = 'fftw_highpass', iqr_fit_lambda_metres: float = 500.0, iqr_detrend_radius_metres: float = 500.0, iqr_window_radius_metres: float = 50.0, iqr_window_shape: str = 'circular', iqr_threshold_percentile: float = 25.0, use_nvd: bool = False, nvd_polyfit_lambda: float = 10.0, nvd_window_radius_metres: float = 100.0, nvd_window_shape: str = 'circular', nvd_threshold_percentile: float = 25.0, require_all_metrics_pass: bool = True, directional_close_line_length_pixels: int = 0, directional_close_n_angles: int = 18, min_component_pixels: int = 100, relief_window_side_pixels: int = 5, relief_threshold_metres: float = 1.0, remove_risers: bool = False, n_threads: int = 0) -> lsdtt3.FlatSurfaceResult
```

---

## `frangi_filter`

frangi_filter(raster: lsdtt3.Raster, sigma_min_metres: float = 1.0, sigma_max_metres: float = 5.0, n_scales: int = 4, beta: float = 0.5, c: float = 1.0, target: str = 'ridges', n_threads: int = 0) -> lsdtt3.FrangiResult

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
    frangi_filter(raster: lsdtt3.Raster, sigma_min_metres: float = 1.0, sigma_max_metres: float = 5.0, n_scales: int = 4, beta: float = 0.5, c: float = 1.0, target: str = 'ridges', n_threads: int = 0) -> lsdtt3.FrangiResult

Run a multiscale Frangi ridge or valley filter.

    Authors
    -------
    B. Gailleton


```python
def frangi_filter(raster: lsdtt3.Raster, sigma_min_metres: float = 1.0, sigma_max_metres: float = 5.0, n_scales: int = 4, beta: float = 0.5, c: float = 1.0, target: str = 'ridges', n_threads: int = 0) -> lsdtt3.FrangiResult
```

---

## `fundamental_form_principal_curvatures`

fundamental_form_principal_curvatures(raster: lsdtt3.Raster, n_threads: int = 0) -> dict

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
    fundamental_form_principal_curvatures(raster: lsdtt3.Raster, n_threads: int = 0) -> dict

Compute differential-geometry maximum and minimum principal curvatures.

    Authors
    -------
    B. Gailleton


```python
def fundamental_form_principal_curvatures(raster: lsdtt3.Raster, n_threads: int = 0) -> dict
```

---

## `halve_raster`

halve_raster(source: lsdtt3.Raster, method: str = 'conservative_mean') -> lsdtt3.Raster

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
    halve_raster(source: lsdtt3.Raster, method: str = 'conservative_mean') -> lsdtt3.Raster

Halve both even raster dimensions over the same extent using the optimized path.

    Authors
    -------
    B. Gailleton


```python
def halve_raster(source: lsdtt3.Raster, method: str = 'conservative_mean') -> lsdtt3.Raster
```

---

## `highpass_filter`

highpass_filter(raster: lsdtt3.Raster, cutoff_metres: float, window: str = 'hann', clip_to_valid_rectangle: bool = False, n_threads: int = 0) -> lsdtt3.Raster

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
    highpass_filter(raster: lsdtt3.Raster, cutoff_metres: float, window: str = 'hann', clip_to_valid_rectangle: bool = False, n_threads: int = 0) -> lsdtt3.Raster

Return the short-wavelength DEM component.

    Authors
    -------
    B. Gailleton


```python
def highpass_filter(raster: lsdtt3.Raster, cutoff_metres: float, window: str = 'hann', clip_to_valid_rectangle: bool = False, n_threads: int = 0) -> lsdtt3.Raster
```

---

## `hillshade`

hillshade(raster: lsdtt3.Raster, azimuth: float = 315.0, altitude: float = 45.0, z_factor: float = 1.0) -> lsdtt3.Raster

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
    hillshade(raster: lsdtt3.Raster, azimuth: float = 315.0, altitude: float = 45.0, z_factor: float = 1.0) -> lsdtt3.Raster

Compute a single-illumination hillshade in the range 0 to 255.

    Authors
    -------
    B. Gailleton


```python
def hillshade(raster: lsdtt3.Raster, azimuth: float = 315.0, altitude: float = 45.0, z_factor: float = 1.0) -> lsdtt3.Raster
```

---

## `interpolate_to_nearest_power_of_two`

interpolate_to_nearest_power_of_two(source: lsdtt3.Raster, direction: str = 'nearest', method: str = 'linear') -> lsdtt3.Raster

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
    interpolate_to_nearest_power_of_two(source: lsdtt3.Raster, direction: str = 'nearest', method: str = 'linear') -> lsdtt3.Raster

Resize to lower, higher, or nearest power-of-two dimensions.

The full extent is preserved. Rasters whose two dimension scale factors differ
cannot be represented with square cells and are rejected.

    Authors
    -------
    B. Gailleton


```python
def interpolate_to_nearest_power_of_two(source: lsdtt3.Raster, direction: str = 'nearest', method: str = 'linear') -> lsdtt3.Raster
```

---

## `invert_binary_mask`

invert_binary_mask(raster: lsdtt3.Raster) -> lsdtt3.Raster

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
    invert_binary_mask(raster: lsdtt3.Raster) -> lsdtt3.Raster

Invert a binary raster mask while preserving nodata.

    Authors
    -------
    B. Gailleton


```python
def invert_binary_mask(raster: lsdtt3.Raster) -> lsdtt3.Raster
```

---

## `least_cost_carve_depth`

least_cost_carve_depth(raster: lsdtt3.Raster, max_distance: int = 1000, max_cost: float = 3.4028234663852886e+38, prefer_shortest_path: bool = False, flat_increment: float = 0.0, fill_remainder: bool = True, n_threads: int = 0, edge_policy: str = 'outlet', fill_min_slope: float = 9.999999747378752e-05) -> lsdtt3.Raster

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
    least_cost_carve_depth(raster: lsdtt3.Raster, max_distance: int = 1000, max_cost: float = 3.4028234663852886e+38, prefer_shortest_path: bool = False, flat_increment: float = 0.0, fill_remainder: bool = True, n_threads: int = 0, edge_policy: str = 'outlet', fill_min_slope: float = 9.999999747378752e-05) -> lsdtt3.Raster

Return the non-negative elevation removed by least-cost breaching (parameters as breach_least_cost).

    Authors
    -------
    B. Gailleton


```python
def least_cost_carve_depth(raster: lsdtt3.Raster, max_distance: int = 1000, max_cost: float = 3.4028234663852886e+38, prefer_shortest_path: bool = False, flat_increment: float = 0.0, fill_remainder: bool = True, n_threads: int = 0, edge_policy: str = 'outlet', fill_min_slope: float = 9.999999747378752e-05) -> lsdtt3.Raster
```

---

## `lowpass_filter`

lowpass_filter(raster: lsdtt3.Raster, cutoff_metres: float, window: str = 'hann', clip_to_valid_rectangle: bool = False, n_threads: int = 0) -> lsdtt3.Raster

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
    lowpass_filter(raster: lsdtt3.Raster, cutoff_metres: float, window: str = 'hann', clip_to_valid_rectangle: bool = False, n_threads: int = 0) -> lsdtt3.Raster

Return the long-wavelength DEM component.

    Authors
    -------
    B. Gailleton


```python
def lowpass_filter(raster: lsdtt3.Raster, cutoff_metres: float, window: str = 'hann', clip_to_valid_rectangle: bool = False, n_threads: int = 0) -> lsdtt3.Raster
```

---

## `make_boundary_conditions`

make_boundary_conditions(north: str = 'open', east: str = 'open', south: str = 'open', west: str = 'open') -> str

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
    make_boundary_conditions(north: str = 'open', east: str = 'open', south: str = 'open', west: str = 'open') -> str

Create an LSDTopoTools NESW boundary-condition string.

Parameters
----------
north, east, south, west : {"open", "closed", "periodic", "base_level"}
    Boundary treatment. The returned string has NESW order.

    Authors
    -------
    B. Gailleton


```python
def make_boundary_conditions(north: str = 'open', east: str = 'open', south: str = 'open', west: str = 'open') -> str
```

---

## `make_outlet_mask`

make_outlet_mask(raster: lsdtt3.Raster, north: bool = False, east: bool = False, south: bool = False, west: bool = False) -> lsdtt3.Raster

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
    make_outlet_mask(raster: lsdtt3.Raster, north: bool = False, east: bool = False, south: bool = False, west: bool = False) -> lsdtt3.Raster

Create a valid binary mask with selected raster edges set to one.

    Authors
    -------
    B. Gailleton


```python
def make_outlet_mask(raster: lsdtt3.Raster, north: bool = False, east: bool = False, south: bool = False, west: bool = False) -> lsdtt3.Raster
```

---

## `mask_by_binary`

mask_by_binary(values: lsdtt3.Raster, keep_mask: lsdtt3.Raster) -> lsdtt3.Raster

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
    mask_by_binary(values: lsdtt3.Raster, keep_mask: lsdtt3.Raster) -> lsdtt3.Raster

Set values to nodata where keep_mask is false or nodata.

    Authors
    -------
    B. Gailleton


```python
def mask_by_binary(values: lsdtt3.Raster, keep_mask: lsdtt3.Raster) -> lsdtt3.Raster
```

---

## `mask_from_points`

mask_from_points(raster: lsdtt3.Raster, points: object) -> lsdtt3.Raster

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
    mask_from_points(raster: lsdtt3.Raster, points: object) -> lsdtt3.Raster

Create a binary mask from zero-based ``(row, col)`` point coordinates.

Parameters
----------
points : array-like, shape (n, 2)
    Row-major NumPy coordinates. Points must be valid raster cells.

    Authors
    -------
    B. Gailleton


```python
def mask_from_points(raster: lsdtt3.Raster, points: object) -> lsdtt3.Raster
```

---

## `mask_to_nodata`

mask_to_nodata(raster: lsdtt3.Raster, threshold: float, below: bool = True) -> lsdtt3.Raster

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
    mask_to_nodata(raster: lsdtt3.Raster, threshold: float, below: bool = True) -> lsdtt3.Raster

Set valid pixels below (or above) a threshold to nodata.

    Authors
    -------
    B. Gailleton


```python
def mask_to_nodata(raster: lsdtt3.Raster, threshold: float, below: bool = True) -> lsdtt3.Raster
```

---

## `masked_surface_recovery`

masked_surface_recovery(raster: lsdtt3.Raster, mask: lsdtt3.Raster, cutoff_metres: float, window: str = 'hann', n_threads: int = 0) -> lsdtt3.Raster

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
    masked_surface_recovery(raster: lsdtt3.Raster, mask: lsdtt3.Raster, cutoff_metres: float, window: str = 'hann', n_threads: int = 0) -> lsdtt3.Raster

Recover a smooth surface from selected pixels.

    Authors
    -------
    B. Gailleton


```python
def masked_surface_recovery(raster: lsdtt3.Raster, mask: lsdtt3.Raster, cutoff_metres: float, window: str = 'hann', n_threads: int = 0) -> lsdtt3.Raster
```

---

## `mosaic`

mosaic(rasters: collections.abc.Sequence[lsdtt3.Raster]) -> lsdtt3.Raster

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
    mosaic(rasters: collections.abc.Sequence[lsdtt3.Raster]) -> lsdtt3.Raster

Mosaic loaded rasters into their union extent.

Parameters
----------
rasters : sequence of Raster
    Rasters with matching CRS, cell size, and pixel alignment. Later valid
    pixels take precedence in overlaps.

Returns
-------
Raster
    Mosaic covering the union extent.

    Authors
    -------
    B. Gailleton


```python
def mosaic(rasters: collections.abc.Sequence[lsdtt3.Raster]) -> lsdtt3.Raster
```

---

## `multidirectional_hillshade`

multidirectional_hillshade(raster: lsdtt3.Raster, n_directions: int = 8, altitude: float = 45.0, z_factor: float = 1.0, start_azimuth: float = 0.0) -> lsdtt3.Raster

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
    multidirectional_hillshade(raster: lsdtt3.Raster, n_directions: int = 8, altitude: float = 45.0, z_factor: float = 1.0, start_azimuth: float = 0.0) -> lsdtt3.Raster

Compute an averaged multi-illumination hillshade in the range 0 to 255.

    Authors
    -------
    B. Gailleton


```python
def multidirectional_hillshade(raster: lsdtt3.Raster, n_directions: int = 8, altitude: float = 45.0, z_factor: float = 1.0, start_azimuth: float = 0.0) -> lsdtt3.Raster
```

---

## `multiply`

multiply(a: lsdtt3.Raster, b: lsdtt3.Raster) -> lsdtt3.Raster

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
    multiply(a: lsdtt3.Raster, b: lsdtt3.Raster) -> lsdtt3.Raster

Multiply two grid-compatible rasters.

    Authors
    -------
    B. Gailleton


```python
def multiply(a: lsdtt3.Raster, b: lsdtt3.Raster) -> lsdtt3.Raster
```

---

## `normalize_unit_interval`

normalize_unit_interval(raster: lsdtt3.Raster) -> lsdtt3.Raster

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
    normalize_unit_interval(raster: lsdtt3.Raster) -> lsdtt3.Raster

Linearly rescale valid values to the closed interval [0, 1].

    Authors
    -------
    B. Gailleton


```python
def normalize_unit_interval(raster: lsdtt3.Raster) -> lsdtt3.Raster
```

---

## `opentopography_datasets`

opentopography_datasets() -> list

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
    opentopography_datasets() -> list

Return OpenTopography DEM catalogue records as dictionaries.

    Authors
    -------
    B. Gailleton


```python
def opentopography_datasets() -> list
```

---

## `pad_to_nearest_power_of_two`

pad_to_nearest_power_of_two(source: lsdtt3.Raster) -> lsdtt3.Raster

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
    pad_to_nearest_power_of_two(source: lsdtt3.Raster) -> lsdtt3.Raster

Pad the north and east sides with nodata to the next power-of-two dimensions.

    Authors
    -------
    B. Gailleton


```python
def pad_to_nearest_power_of_two(source: lsdtt3.Raster) -> lsdtt3.Raster
```

---

## `plan_curvature`

plan_curvature(raster: lsdtt3.Raster) -> lsdtt3.Raster

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
    plan_curvature(raster: lsdtt3.Raster) -> lsdtt3.Raster

Compute 3-by-3 plan curvature; concave-up is positive.

    Authors
    -------
    B. Gailleton


```python
def plan_curvature(raster: lsdtt3.Raster) -> lsdtt3.Raster
```

---

## `polyfit_metrics`

polyfit_metrics(raster: lsdtt3.Raster, metrics: collections.abc.Sequence[str], lambda_metres: float = 0.0, z_factor: float = 1.0, n_threads: int = 0, flat_slope_tolerance: float = 1e-12) -> dict

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
    polyfit_metrics(raster: lsdtt3.Raster, metrics: collections.abc.Sequence[str], lambda_metres: float = 0.0, z_factor: float = 1.0, n_threads: int = 0, flat_slope_tolerance: float = 1e-12) -> dict

Compute several quadratic-fit surface metrics in one shared fit.

Parameters
----------
raster : Raster
    Input elevation raster.
metrics : list[str]
    Requested metric names. See the Python documentation for valid names.
lambda_metres : float, default=0
    Fit scale in metres. Zero uses three raster cells.
z_factor : float, default=1
    Elevation scale factor applied before metric calculation.
n_threads : int, default=0
    OpenMP thread count; zero selects the runtime default.
flat_slope_tolerance : float, default=1e-12
    Directional metrics are nodata where fitted slope is below this value.

Returns
-------
dict[str, Raster]
    One newly allocated raster for every requested metric.

    Authors
    -------
    B. Gailleton


```python
def polyfit_metrics(raster: lsdtt3.Raster, metrics: collections.abc.Sequence[str], lambda_metres: float = 0.0, z_factor: float = 1.0, n_threads: int = 0, flat_slope_tolerance: float = 1e-12) -> dict
```

---

## `rasterize_polygons`

rasterize_polygons(raster: lsdtt3.Raster, polygons: object, coverage: str = 'majority') -> lsdtt3.Raster

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
    rasterize_polygons(raster: lsdtt3.Raster, polygons: object, coverage: str = 'majority') -> lsdtt3.Raster

Rasterize GeoPandas polygon geometries onto a raster grid.

Parameters
----------
polygons : geopandas.GeoDataFrame
    Polygon geometries with a CRS exactly matching raster.
coverage : {"any", "majority", "full"}, default="majority"
    Pixel-inclusion rule. ``"majority"`` tests the pixel centre.

Returns
-------
Raster
    Binary mask raster with valid pixels throughout the reference grid.

Raises
------
ValueError
    If the GeoDataFrame CRS does not match the raster CRS.

    Authors
    -------
    B. Gailleton


```python
def rasterize_polygons(raster: lsdtt3.Raster, polygons: object, coverage: str = 'majority') -> lsdtt3.Raster
```

---

## `reproject`

reproject(raster: lsdtt3.Raster, epsg: int, cell_size: float = 0.0, resampling: str = 'bilinear') -> lsdtt3.Raster

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
    reproject(raster: lsdtt3.Raster, epsg: int, cell_size: float = 0.0, resampling: str = 'bilinear') -> lsdtt3.Raster

Reproject a raster to an EPSG coordinate reference system.

Parameters
----------
epsg : int
    Target EPSG identifier.
cell_size : float, default=0
    Target cell size. A non-positive value lets GDAL choose it.
resampling : str, default="bilinear"
    GDAL resampling method.

    Authors
    -------
    B. Gailleton


```python
def reproject(raster: lsdtt3.Raster, epsg: int, cell_size: float = 0.0, resampling: str = 'bilinear') -> lsdtt3.Raster
```

---

## `reproject_to_grid`

reproject_to_grid(raster: lsdtt3.Raster, target_metadata: lsdtt3.RasterMeta, resampling: str = 'bilinear') -> lsdtt3.Raster

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
    reproject_to_grid(raster: lsdtt3.Raster, target_metadata: lsdtt3.RasterMeta, resampling: str = 'bilinear') -> lsdtt3.Raster

Reproject and resample a raster onto an explicit target grid.

Parameters
----------
raster : Raster
target_metadata : RasterMeta
    Target extent, resolution, dimensions, and CRS.
resampling : str, default="bilinear"
    GDAL resampling method.

Returns
-------
Raster
    New raster on target_metadata's grid.

    Authors
    -------
    B. Gailleton


```python
def reproject_to_grid(raster: lsdtt3.Raster, target_metadata: lsdtt3.RasterMeta, resampling: str = 'bilinear') -> lsdtt3.Raster
```

---

## `resample_raster`

resample_raster(source: lsdtt3.Raster, target_metadata: lsdtt3.RasterMeta, method: str = 'linear') -> lsdtt3.Raster

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
    resample_raster(source: lsdtt3.Raster, target_metadata: lsdtt3.RasterMeta, method: str = 'linear') -> lsdtt3.Raster

Resample an in-memory raster onto an explicit target grid without GDAL.

The source and target CRS must match. Sampling methods are ``nearest``,
``linear`` and ``cubic``. Aggregation methods are ``mean``, ``sum``, ``min``,
``max``, ``mode`` and the overlap-area-weighted ``conservative_mean``.

    Authors
    -------
    B. Gailleton


```python
def resample_raster(source: lsdtt3.Raster, target_metadata: lsdtt3.RasterMeta, method: str = 'linear') -> lsdtt3.Raster
```

---

## `sample_points`

sample_points(raster: lsdtt3.Raster, points: object, column: str = 'raster_value') -> object

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
    sample_points(raster: lsdtt3.Raster, points: object, column: str = 'raster_value') -> object

Sample a raster at GeoPandas point geometries.

Parameters
----------
points : geopandas.GeoDataFrame
    Point geometries with a CRS exactly matching raster.
column : str, default="raster_value"
    Name of the output value column.

Returns
-------
geopandas.GeoDataFrame
    Copy of points with sampled values. Outside or nodata samples are None.

    Authors
    -------
    B. Gailleton


```python
def sample_points(raster: lsdtt3.Raster, points: object, column: str = 'raster_value') -> object
```

---

## `scale`

scale(raster: lsdtt3.Raster, factor: float) -> lsdtt3.Raster

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
    scale(raster: lsdtt3.Raster, factor: float) -> lsdtt3.Raster

Multiply valid raster values by a scalar factor.

    Authors
    -------
    B. Gailleton


```python
def scale(raster: lsdtt3.Raster, factor: float) -> lsdtt3.Raster
```

---

## `set_nodata_rectangle`

set_nodata_rectangle(raster: lsdtt3.Raster, row_start: int, row_stop: int, col_start: int, col_stop: int) -> lsdtt3.Raster

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
    set_nodata_rectangle(raster: lsdtt3.Raster, row_start: int, row_stop: int, col_start: int, col_stop: int) -> lsdtt3.Raster

Set a half-open row/column rectangle to nodata.

    Authors
    -------
    B. Gailleton


```python
def set_nodata_rectangle(raster: lsdtt3.Raster, row_start: int, row_stop: int, col_start: int, col_stop: int) -> lsdtt3.Raster
```

---

## `slope`

slope(raster: lsdtt3.Raster) -> lsdtt3.Raster

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
    slope(raster: lsdtt3.Raster) -> lsdtt3.Raster

Compute 3-by-3 gradient-magnitude slope in rise/run units.

    Authors
    -------
    B. Gailleton


```python
def slope(raster: lsdtt3.Raster) -> lsdtt3.Raster
```

---

## `spectral_analysis`

spectral_analysis(raster: lsdtt3.Raster, log_bin_width: float = 0.05, n_background_iterations: int = 20, seed: int = -1, window: str = 'hann', clip_to_valid_rectangle: bool = False, n_threads: int = 0) -> lsdtt3.SpectralResult

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
    spectral_analysis(raster: lsdtt3.Raster, log_bin_width: float = 0.05, n_background_iterations: int = 20, seed: int = -1, window: str = 'hann', clip_to_valid_rectangle: bool = False, n_threads: int = 0) -> lsdtt3.SpectralResult

Run full FFT-based spectral analysis and return a SpectralResult.

    Authors
    -------
    B. Gailleton


```python
def spectral_analysis(raster: lsdtt3.Raster, log_bin_width: float = 0.05, n_background_iterations: int = 20, seed: int = -1, window: str = 'hann', clip_to_valid_rectangle: bool = False, n_threads: int = 0) -> lsdtt3.SpectralResult
```

---

## `spectral_available`

spectral_available() -> bool

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
    spectral_available() -> bool

Return whether this build includes FFTW spectral support.

    Authors
    -------
    B. Gailleton


```python
def spectral_available() -> bool
```

---

## `subtract`

subtract(a: lsdtt3.Raster, b: lsdtt3.Raster) -> lsdtt3.Raster

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
    subtract(a: lsdtt3.Raster, b: lsdtt3.Raster) -> lsdtt3.Raster

Subtract two grid-compatible rasters.

    Authors
    -------
    B. Gailleton


```python
def subtract(a: lsdtt3.Raster, b: lsdtt3.Raster) -> lsdtt3.Raster
```

---

## `swath_distances`

swath_distances(raster: lsdtt3.Raster, baseline: object, n_threads: int = 0) -> dict

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
    swath_distances(raster: lsdtt3.Raster, baseline: object, n_threads: int = 0) -> dict

Compute pixel distances relative to a GeoPandas baseline.

Parameters
----------
raster : Raster
baseline : geopandas.GeoDataFrame
    Exactly one LineString in the same projected metre CRS as raster.
n_threads : int, default=0
    OpenMP thread count; zero uses the available threads.

Returns
-------
dict[str, Raster]
    ``along_axis``, ``perpendicular``, and ``signed_perpendicular`` rasters.

    Authors
    -------
    B. Gailleton


```python
def swath_distances(raster: lsdtt3.Raster, baseline: object, n_threads: int = 0) -> dict
```

---

## `swath_profile`

swath_profile(reference: lsdtt3.Raster, baseline: object, values: object | None = None, half_width_metres: float = 0.0, bin_width_metres: float = 1000.0, n_threads: int = 0) -> lsdtt3.SwathProfile

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
    swath_profile(reference: lsdtt3.Raster, baseline: object, values: object | None = None, half_width_metres: float = 0.0, bin_width_metres: float = 1000.0, n_threads: int = 0) -> lsdtt3.SwathProfile

Summarise raster values in a corridor along a GeoPandas LineString.

Parameters
----------
reference : Raster
    Defines the swath grid and valid-pixel domain.
baseline : geopandas.GeoDataFrame
    Exactly one LineString in the same projected metre CRS as reference.
values : Raster, optional
    Co-registered raster to summarise. Defaults to reference.
half_width_metres : float, default=0
    Corridor half-width. Zero includes all valid pixels.
bin_width_metres : float, default=1000
    Along-axis bin width; must be positive.
n_threads : int, default=0
    OpenMP thread count; zero uses the available threads.

Returns
-------
SwathProfile
    Distance rasters and per-bin mean and percentile statistics.

    Authors
    -------
    B. Gailleton


```python
def swath_profile(reference: lsdtt3.Raster, baseline: object, values: object | None = None, half_width_metres: float = 0.0, bin_width_metres: float = 1000.0, n_threads: int = 0) -> lsdtt3.SwathProfile
```

---

## `threshold_above`

threshold_above(raster: lsdtt3.Raster, threshold: float) -> lsdtt3.Raster

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
    threshold_above(raster: lsdtt3.Raster, threshold: float) -> lsdtt3.Raster

Return a binary mask that is True where valid values exceed threshold.

    Authors
    -------
    B. Gailleton


```python
def threshold_above(raster: lsdtt3.Raster, threshold: float) -> lsdtt3.Raster
```

---

## `threshold_below`

threshold_below(raster: lsdtt3.Raster, threshold: float) -> lsdtt3.Raster

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
    threshold_below(raster: lsdtt3.Raster, threshold: float) -> lsdtt3.Raster

Return a binary mask that is True where valid values are below threshold.

    Authors
    -------
    B. Gailleton


```python
def threshold_below(raster: lsdtt3.Raster, threshold: float) -> lsdtt3.Raster
```

---

## `topographic_wetness_index`

topographic_wetness_index(drainage_area: lsdtt3.Raster, slope: lsdtt3.Raster) -> lsdtt3.Raster

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
    topographic_wetness_index(drainage_area: lsdtt3.Raster, slope: lsdtt3.Raster) -> lsdtt3.Raster

Compute ``log((A / cell_area) / slope)`` from drainage area and slope rasters.

    Authors
    -------
    B. Gailleton


```python
def topographic_wetness_index(drainage_area: lsdtt3.Raster, slope: lsdtt3.Raster) -> lsdtt3.Raster
```

---

## `total_curvature`

total_curvature(raster: lsdtt3.Raster, radius_metres: float) -> lsdtt3.Raster

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
    total_curvature(raster: lsdtt3.Raster, radius_metres: float) -> lsdtt3.Raster

Compute moving-window total curvature at radius_metres.

    Authors
    -------
    B. Gailleton


```python
def total_curvature(raster: lsdtt3.Raster, radius_metres: float) -> lsdtt3.Raster
```

---

## `wiener_filter`

wiener_filter(raster: lsdtt3.Raster, signal_wavelength_long_metres: float = 1000.0, signal_wavelength_short_metres: float = 100.0, window: str = 'hann', clip_to_valid_rectangle: bool = False, n_threads: int = 0) -> lsdtt3.Raster

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
    wiener_filter(raster: lsdtt3.Raster, signal_wavelength_long_metres: float = 1000.0, signal_wavelength_short_metres: float = 100.0, window: str = 'hann', clip_to_valid_rectangle: bool = False, n_threads: int = 0) -> lsdtt3.Raster

Attenuate high-frequency white noise with an adaptive Wiener filter.

    Authors
    -------
    B. Gailleton


**Deprecated compatibility names:**
- `signal_wavelength_long_m` -> `signal_wavelength_long_metres`; earliest removal: following minor.
- `signal_wavelength_short_m` -> `signal_wavelength_short_metres`; earliest removal: following minor.

```python
def wiener_filter(raster: lsdtt3.Raster, signal_wavelength_long_metres: float = 1000.0, signal_wavelength_short_metres: float = 100.0, window: str = 'hann', clip_to_valid_rectangle: bool = False, n_threads: int = 0) -> lsdtt3.Raster
```
