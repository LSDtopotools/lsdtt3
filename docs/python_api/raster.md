# Raster And Terrain API

This page is generated from the built Python runtime manifest. Do not edit signatures here.

## `FlatSurfaceClassification`

Floodplain and terrace classification of labelled flats.

Made by :func:`classify_flat_surfaces` or
:func:`classify_flat_surfaces_with_scope`. ``table`` has one row per label:
``label``, ``surface_class`` (``"floodplain"``, ``"terrace"`` or
``"none"``), ``source_id``, ``n_pixels``, ``median_height``,
``mean_height`` (HAND, m), ``min_distance_to_channel_metres``,
``adjacent_to_channel``.

Returned by ``classify_flat_surfaces()``, ``classify_flat_surfaces_with_scope()``; not constructible directly.

Attributes
----------
floodplain_mask : Raster
global_breaks_metres : list[float]
hand_band : Raster
n_floodplain_labels : int
n_terrace_labels : int
surface_class : Raster
table : pandas.DataFrame
    One row per flat label; columns as in the class documentation.

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

One row per flat label; columns as in the class documentation.

---

## `FlatSurfaceResult`

Results from texture- and relief-based flat-surface detection.

Returned by ``find_flat_surfaces()``; not constructible directly.

Attributes
----------
combined_flat_mask : Raster
healed_flat_mask : Raster
iqr_flat_mask : Raster
iqr_residuals : Raster
iqr_texture : Raster
labelled_flat_mask : Raster
madk2_flat_mask : Raster
madk2_roughness : Raster
nvd_flat_mask : Raster
nvd_s3 : Raster
relief : Raster
relief_flat_mask : Raster

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

Multiscale Frangi filter output from :func:`frangi_filter`.

``response`` is the largest Frangi line-likeness response over the scales and
``best_scale_metres`` the scale (m) where it occurred.

Returned by ``frangi_filter()``; not constructible directly.

Attributes
----------
best_scale_metres : Raster
response : Raster

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

An in-memory single-band float32 raster with an explicit validity mask.

Instances are value objects: functions in this module return new rasters and
do not modify their inputs. Values are float32; each pixel is either valid or
nodata (``mask_numpy()`` is True where valid). Array row 0 is the northern
edge (see ``RasterMeta``). Typical round trip::

    dem = lsdtt3.Raster.read("dem.tif")
    z = dem.to_numpy()                 # float32, shape (n_rows, n_cols)
    out = lsdtt3.Raster.from_numpy(z * 2, dem.mask_numpy(), dem.metadata)
    out.write("doubled.tif")

Not constructible directly; create one with ``Raster.from_numpy()`` or ``Raster.read()``. It is also returned by ``add()``, ``add_uniform_noise()``, ``aspect()``, ``breach_depressions()`` and 243 more.

Attributes
----------
metadata : RasterMeta
    Grid metadata (RasterMeta): shape, lower-left corner, cell size, nodata value, CRS.
shape : tuple[int, int]
    Grid shape as (n_rows, n_cols), the shape of to_numpy().
valid_pixel_count : int
    Number of valid (non-nodata) pixels.

Methods
-------
Raster.from_numpy(data: object, mask: object, metadata: RasterMeta) -> Raster
    Create a raster from copied two-dimensional NumPy-compatible arrays.
mask_numpy() -> object
    Return a copied Boolean array, True where raster values are valid.
Raster.read(path: str) -> Raster
    Read the first band of a GDAL-readable raster file.
to_numpy() -> object
    Return a copied ``float32`` array of raster values.
write(path: str, format: str = 'cog') -> None
    Write this raster to a file with GDAL.

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
    Values with shape ``(metadata.n_rows, metadata.n_cols)``; row 0 is the
    northern edge. Copied and converted to float32.
mask : array-like of bool
    Same shape; True for valid pixels. False pixels are stored as
    ``metadata.nodata_value``. Use ``np.isfinite(data)`` or another raster's
    ``mask_numpy()``.
metadata : RasterMeta
    Grid metadata matching both array dimensions, usually an existing
    raster's ``.metadata`` or ``RasterMeta.from_shape(...)``.

Returns
-------
Raster
    New in-memory raster. Raises ``ValueError`` if a shape does not match.

```python
def from_numpy(data: object, mask: object, metadata: lsdtt3.RasterMeta) -> lsdtt3.Raster
```

#### `mask_numpy`

mask_numpy(self) -> object

Return a copied Boolean array, True where raster values are valid.

Returns
-------
numpy.ndarray
    bool array of shape ``(n_rows, n_cols)``.

```python
def mask_numpy(self) -> object
```

#### `metadata`

Grid metadata (RasterMeta): shape, lower-left corner, cell size, nodata value, CRS.

#### `read`

read(path: str) -> lsdtt3.Raster

Read the first band of a GDAL-readable raster file.

The raster must be north-up with square pixels. Pixels equal to the file's
nodata value become invalid.

Parameters
----------
path : str
    Raster file to load (GeoTIFF, COG, or another GDAL-readable format).

Returns
-------
Raster
    Raster values, metadata, and nodata mask.

```python
def read(path: str) -> lsdtt3.Raster
```

#### `shape`

Grid shape as (n_rows, n_cols), the shape of to_numpy().

#### `to_numpy`

to_numpy(self) -> object

Return a copied ``float32`` array of raster values.

Invalid cells contain the raster nodata value (not NaN); use ``mask_numpy``
to identify them, for example
``np.where(r.mask_numpy(), r.to_numpy(), np.nan)``.

Returns
-------
numpy.ndarray
    float32 array of shape ``(n_rows, n_cols)``; row 0 is the northern edge.

```python
def to_numpy(self) -> object
```

#### `valid_pixel_count`

Number of valid (non-nodata) pixels.

#### `write`

write(self, path: str, format: str = 'cog') -> None

Write this raster to a file with GDAL.

Parameters
----------
path : str
    Destination file (for example ``"slope.tif"``).
format : str, default="cog"
    ``"cog"`` (Cloud-Optimised GeoTIFF), ``"gtiff"`` (plain GeoTIFF) or
    ``"png"``. An unknown name prints a warning and writes a COG.

Returns
-------
None

```python
def write(self, path: str, format: str = 'cog') -> None
```

---

## `RasterHierarchy`

Resolution hierarchy for one in-memory float Raster.

Level 0 is finest. Standard levels are exact factor-two grids; explicit
integer factors can be added over the same extent.
The hierarchy stores no solver state and each call returning a Raster returns
a copy, consistent with the rest of the Python Raster API.

Construct it as

    RasterHierarchy(finest: Raster, storage: str = 'retain_allocated', keep_finest: bool = True, levels: int = 0)

Parameters
----------
finest : Raster
    Finest-resolution raster (level 0).
storage : {"active_only", "retain_allocated", "preallocated"}, default="retain_allocated"
    Which level rasters are kept in memory: only the active level, every
    level once allocated, or every level from the start.
keep_finest : bool, default=True
    Always keep level 0 in memory.
levels : int, default=0
    Number of standard factor-two levels to lay out; 0 uses every exact
    factor-two level the grid allows.

Attributes
----------
active_level : int
n_levels : int

Methods
-------
add_level_for_factor(factor: int) -> int
allocate_to(level: int) -> None
clear() -> None
clear_except(level: int) -> None
has_level(level: int) -> bool
level(level: int) -> Raster
level_factor(level: int) -> int
    Integer coarsening factor relative to level zero.
level_shape(level: int) -> tuple
    Grid shape as (n_rows, n_cols).
move_to(level: int, method: str) -> Raster
prolong_to(level: int, method: str = 'linear') -> Raster
restrict_to(level: int, method: str = 'conservative_mean') -> Raster
set_active(raster: Raster) -> None

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

``x_min`` and ``y_min`` are the map coordinates of the lower-left (south-west)
corner of the grid, but array row 0 is the TOP (north) row, as in GDAL and
NumPy images: pixel ``(row, col)`` has its centre at
``x = x_min + (col + 0.5) * cell_size`` and
``y = y_max - (row + 0.5) * cell_size``. Coordinates and ``cell_size`` are in
the units of the CRS (metres for the projected CRSs lsdtt3 analyses need).
Every ``Raster`` and ``FlowInfo`` carries one as ``.metadata``; create one with
``RasterMeta.from_shape`` to build a raster from NumPy arrays.

Not constructible directly; create one with ``RasterMeta.from_shape()``. It is also returned by ``Raster.metadata``, ``FlowInfo.metadata``, ``FlowGraph.metadata``.

Attributes
----------
cell_size : float
    Square cell width and height, in CRS units (metres when projected).
crs_wkt : str
    Coordinate reference system as WKT; empty if none.
epsg_code : int
    EPSG code of the CRS if known, else 0.
is_valid : bool
    True if n_rows, n_cols and cell_size are all positive.
n_cols : int
    Number of columns (west to east).
n_pixels : int
    Total number of pixels, n_rows * n_cols (valid or not).
n_rows : int
    Number of rows (north to south).
nodata_value : float
    Value stored in invalid pixels.
pixel_area : float
    Area of one pixel, cell_size squared (m^2 for a metre CRS).
x_max : float
    Easting of the eastern edge: x_min + n_cols * cell_size.
x_min : float
    Easting of the western edge of the grid, in CRS units.
y_max : float
    Northing of the northern edge (the top of row 0): y_min + n_rows * cell_size.
y_min : float
    Northing of the southern edge of the grid (the bottom of the last row).

Methods
-------
RasterMeta.from_shape(n_rows: int, n_cols: int, x_min: float = 0.0, y_min: float = 0.0, cell_size: float = 1.0, nodata_value: float = -9999.0, crs_wkt: str = '', epsg_code: int = 0) -> RasterMeta
    Create raster metadata from a grid shape and the grid's lower-left corner.
lonlat_to_row_col(longitude: object, latitude: object, crs: str = 'EPSG:4326') -> tuple
    Return the ``(row, col)`` of the pixels containing longitude/latitude points.
row_col_to_xy(row: object, col: object) -> tuple
    Return the map coordinates ``(x, y)`` of pixel centres.
xy_to_row_col(x: object, y: object) -> tuple
    Return the ``(row, col)`` of the pixels containing map points ``(x, y)``.

Authors
-------
B. Gailleton
S. M. Mudd

### Members

#### `cell_size`

Square cell width and height, in CRS units (metres when projected).

#### `crs_wkt`

Coordinate reference system as WKT; empty if none.

#### `epsg_code`

EPSG code of the CRS if known, else 0.

#### `from_shape`

from_shape(n_rows: int, n_cols: int, x_min: float = 0.0, y_min: float = 0.0, cell_size: float = 1.0, nodata_value: float = -9999.0, crs_wkt: str = '', epsg_code: int = 0) -> lsdtt3.RasterMeta

Create raster metadata from a grid shape and the grid's lower-left corner.

``x_min, y_min`` locate the lower-left (south-west) corner of the grid; the
first array row (row 0) is still the northern edge, at
``y_max = y_min + n_rows * cell_size``.

Parameters
----------
n_rows, n_cols : int
    Positive grid dimensions (rows run north to south, columns west to east).
x_min, y_min : float
    Easting and northing of the lower-left corner of the grid, in CRS units
    (metres for a projected CRS).
cell_size : float
    Square cell width and height in CRS units; must be positive.
nodata_value : float, optional
    Value stored in invalid pixels and written as the nodata value.
crs_wkt : str, optional
    Coordinate reference system in WKT form; empty means no CRS.
epsg_code : int, optional
    EPSG identifier when known (for example 32633 for UTM zone 33N); 0 means
    unknown. This only records the code: pass ``crs_wkt`` as well if the
    raster will be written or analysed in map units.

Returns
-------
RasterMeta
    New metadata. Raises ``ValueError`` if ``n_rows``, ``n_cols`` or
    ``cell_size`` is not positive.

```python
def from_shape(n_rows: int, n_cols: int, x_min: float = 0.0, y_min: float = 0.0, cell_size: float = 1.0, nodata_value: float = -9999.0, crs_wkt: str = '', epsg_code: int = 0) -> lsdtt3.RasterMeta
```

#### `is_valid`

True if n_rows, n_cols and cell_size are all positive.

#### `lonlat_to_row_col`

lonlat_to_row_col(self, longitude: object, latitude: object, crs: str = 'EPSG:4326') -> tuple

Return the ``(row, col)`` of the pixels containing longitude/latitude points.

The points are transformed from ``crs`` into the grid's CRS (``crs_wkt``)
with GDAL/PROJ, then located as in ``xy_to_row_col``. Longitude comes first,
as x. To move a point onto a channel use ``FlowInfo.snap_points`` instead.

Parameters
----------
longitude, latitude : float or array-like of float
    Coordinates in decimal degrees (in ``crs`` units for a projected
    ``crs``). Arrays are broadcast against each other.
crs : str, default='EPSG:4326'
    CRS of the input points: an ``"EPSG:<code>"`` string, WKT or PROJ
    string; the default is WGS84 longitude/latitude.

Returns
-------
tuple of (int, int) or tuple of (numpy.ndarray, numpy.ndarray)
    ``(row, col)`` as in ``xy_to_row_col``; -1 for points off the grid or
    that cannot be projected. Raises ``ValueError`` if the grid has no CRS
    or ``crs`` is not understood.

Examples
--------
>>> dem = lsdtt3.fetch_dem(15.63, 37.88, 16.21, 38.31)
>>> row, col = dem.metadata.lonlat_to_row_col(15.79008, 37.94234)

```python
def lonlat_to_row_col(self, longitude: object, latitude: object, crs: str = 'EPSG:4326') -> tuple
```

#### `n_cols`

Number of columns (west to east).

#### `n_pixels`

Total number of pixels, n_rows * n_cols (valid or not).

#### `n_rows`

Number of rows (north to south).

#### `nodata_value`

Value stored in invalid pixels.

#### `pixel_area`

Area of one pixel, cell_size squared (m^2 for a metre CRS).

#### `row_col_to_xy`

row_col_to_xy(self, row: object, col: object) -> tuple

Return the map coordinates ``(x, y)`` of pixel centres.

``x = x_min + (col + 0.5) * cell_size`` and
``y = y_max - (row + 0.5) * cell_size``: row 0 is the northern (top) edge, so
``y`` decreases as ``row`` increases.

Parameters
----------
row, col : int or array-like of int
    Zero-based pixel positions (rows north to south, columns west to east).
    Arrays are broadcast against each other. Positions are not checked
    against the grid: an off-grid position gives the centre the pixel would
    have.

Returns
-------
tuple of (float, float) or tuple of (numpy.ndarray, numpy.ndarray)
    ``(x, y)`` in the grid's CRS units (metres for a projected CRS): Python
    floats for scalar input, otherwise float64 arrays of the broadcast
    input shape.

Examples
--------
>>> meta = lsdtt3.RasterMeta.from_shape(100, 200, x_min=500000.0, y_min=4200000.0, cell_size=30.0)
>>> meta.row_col_to_xy(0, 1)
(500045.0, 4202985.0)

```python
def row_col_to_xy(self, row: object, col: object) -> tuple
```

#### `x_max`

Easting of the eastern edge: x_min + n_cols * cell_size.

#### `x_min`

Easting of the western edge of the grid, in CRS units.

#### `xy_to_row_col`

xy_to_row_col(self, x: object, y: object) -> tuple

Return the ``(row, col)`` of the pixels containing map points ``(x, y)``.

A pixel covers ``x_min + col * cell_size <= x < x_min + (col + 1) * cell_size``
and ``y_max - (row + 1) * cell_size < y <= y_max - row * cell_size``, so row 0
is the northern (top) edge: ``col = floor((x - x_min) / cell_size)`` and
``row = floor((y_max - y) / cell_size)``. The inverse of ``row_col_to_xy`` for
pixel centres.

Parameters
----------
x, y : float or array-like of float
    Easting and northing in the grid's own CRS (``crs_wkt``; metres for a
    projected CRS). Arrays are broadcast against each other.

Returns
-------
tuple of (int, int) or tuple of (numpy.ndarray, numpy.ndarray)
    ``(row, col)``: Python ints for scalar input, otherwise int64 arrays of
    the broadcast input shape. Points outside the grid (including points on
    its eastern or southern boundary line) and non-finite coordinates give
    row = col = -1; they do not raise.

Examples
--------
>>> meta = lsdtt3.RasterMeta.from_shape(100, 200, x_min=500000.0, y_min=4200000.0, cell_size=30.0)
>>> meta.xy_to_row_col(500045.0, 4202985.0)
(0, 1)

```python
def xy_to_row_col(self, x: object, y: object) -> tuple
```

#### `y_max`

Northing of the northern edge (the top of row 0): y_min + n_rows * cell_size.

#### `y_min`

Northing of the southern edge of the grid (the bottom of the last row).

---

## `SoilTopographicIndexResult`

Soil-topographic-index outputs and diagnostics.

Returned by ``soil_topographic_index()``; not constructible directly.

Attributes
----------
lookup : object
stats : dict
sti : Raster
term : Raster
valid_mask : Raster

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

Spectral arrays are returned as copied NumPy arrays. Frequencies are in
cycles per metre and wavelengths are in metres.

Returned by ``spectral_analysis()``; not constructible directly.

Attributes
----------
background_psd : numpy.ndarray
background_psd_2d : numpy.ndarray
binned : dict
frequency : numpy.ndarray
frequency_x : numpy.ndarray
frequency_y : numpy.ndarray
normalised_psd : numpy.ndarray
padded_shape : tuple
psd : numpy.ndarray
psd_2d : numpy.ndarray
rollover_beta : float
rollover_frequency : float
rollover_wavelength : float
sub_rollover_beta : float
wavelength : numpy.ndarray

Authors
-------
B. Gailleton

### Members

#### `background_psd`

background_psd(self) -> numpy.ndarray

#### `background_psd_2d`

background_psd_2d(self) -> numpy.ndarray

#### `binned`

(self) -> dict

#### `frequency`

frequency(self) -> numpy.ndarray

#### `frequency_x`

frequency_x(self) -> numpy.ndarray

#### `frequency_y`

frequency_y(self) -> numpy.ndarray

#### `normalised_psd`

normalised_psd(self) -> numpy.ndarray

#### `padded_shape`

(self) -> tuple

#### `psd`

psd(self) -> numpy.ndarray

#### `psd_2d`

psd_2d(self) -> numpy.ndarray

#### `rollover_beta`

(self) -> float

#### `rollover_frequency`

(self) -> float

#### `rollover_wavelength`

(self) -> float

#### `sub_rollover_beta`

(self) -> float

#### `wavelength`

wavelength(self) -> numpy.ndarray

---

## `SwathProfile`

Results from a corridor-based longitudinal swath profile.

The distance rasters are Raster objects. ``table`` is a copied pandas
DataFrame with one row per along-axis bin.

Returned by ``swath_profile()``; not constructible directly.

Attributes
----------
along_axis : Raster
    Distance (m) of each pixel along the baseline.
perpendicular_distance : Raster
    Unsigned distance (m) of each pixel from the baseline.
quantile_levels : numpy.ndarray
    Percentile levels of the table's p<k> columns (float64 array).
signed_perpendicular_distance : Raster
    Signed distance (m) from the baseline (sign gives the side).
table : pandas.DataFrame
    One row per along-axis bin: along_axis_centre (m), n_pixels, mean, and percentile columns p0, p5, ..., p100 (5% steps).

Authors
-------
B. Gailleton

### Members

#### `along_axis`

Distance (m) of each pixel along the baseline.

#### `perpendicular_distance`

Unsigned distance (m) of each pixel from the baseline.

#### `quantile_levels`

Percentile levels of the table's p<k> columns (float64 array).

#### `signed_perpendicular_distance`

Signed distance (m) from the baseline (sign gives the side).

#### `table`

One row per along-axis bin: along_axis_centre (m), n_pixels, mean, and percentile columns p0, p5, ..., p100 (5% steps).

---

## `absolute`

Absolute value of every valid pixel (for example abs(subtract(a, b)) for a DEM of absolute difference).

Parameters
----------
raster : Raster
    Input raster.

Returns
-------
Raster

Authors
-------
S. M. Mudd
B. Gailleton

```python
def absolute(raster: lsdtt3.Raster) -> lsdtt3.Raster
```

---

## `add`

Add two grid-compatible rasters.

Parameters
----------
a, b : Raster
    Rasters on the same grid (shape, cell size and georeferencing).

Returns
-------
Raster
    ``a + b``; nodata where either input is nodata.

Authors
-------
B. Gailleton

```python
def add(a: lsdtt3.Raster, b: lsdtt3.Raster) -> lsdtt3.Raster
```

---

## `add_uniform_noise`

Add reproducible uniform noise to valid pixels.

Parameters
----------
raster : Raster
    Input raster.
seed : int
    Random seed; the same seed gives the same noise.
magnitude : float
    Noise is drawn uniformly from ``[-magnitude, magnitude]``, in the units
    of the raster (m for a DEM); non-negative.

Returns
-------
Raster
    New in-memory raster; nodata pixels are unchanged.

Authors
-------
B. Gailleton

```python
def add_uniform_noise(raster: lsdtt3.Raster, seed: int, magnitude: float) -> lsdtt3.Raster
```

---

## `aspect`

Compute polyfit aspect in degrees; lambda_metres defaults to three cells.

Aspect is the downslope direction the surface faces, in compass degrees
(0 = north, 90 = east, clockwise), from the polyfit surface.

Parameters
----------
raster : Raster
    Elevation raster (m) in a projected metre CRS.
lambda_metres : float, default=0.0
    Polyfit length scale (m), as in ``polyfit_metrics``; 0 uses three cells.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.
flat_slope_tolerance : float, default=1e-12
    Pixels whose fitted gradient (m/m) is below this are treated as flat
    (no defined aspect).

Returns
-------
Raster
    Aspect in degrees, 0 to 360.

Authors
-------
B. Gailleton

```python
def aspect(raster: lsdtt3.Raster, lambda_metres: float = 0.0, n_threads: int = 0, flat_slope_tolerance: float = 1e-12) -> lsdtt3.Raster
```

---

## `bandpass_filter`

Return the DEM component between two wavelengths.

Parameters
----------
raster : Raster
    Elevation raster (m) in a projected metre CRS.
short_metres, long_metres : float
    Shortest and longest wavelength (m) kept; ``short_metres`` must be less
    than ``long_metres``.
window : {"hann", "hamming", "none"}, default="hann"
    Taper applied before the FFT.
clip_to_valid_rectangle : bool, default=False
    First crop to the largest all-valid rectangle (the result is then
    smaller than the input).
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
Raster

Authors
-------
B. Gailleton

```python
def bandpass_filter(raster: lsdtt3.Raster, short_metres: float, long_metres: float, window: str = 'hann', clip_to_valid_rectangle: bool = False, n_threads: int = 0) -> lsdtt3.Raster
```

---

## `breach_least_cost`

Breach depressions along least-cost paths (Lindsay 2020).

A port of WhiteboxTools BreachDepressionsLeastCost (whitebox_workflows 2.0.6).

Parameters
----------
raster : Raster
    Elevation raster (m).
max_distance : int or None, default=None
    Maximum breach path length in cells, counted in steps from the pit (a
    diagonal step counts as one). 1 to 32767. None selects the edge policy's
    default: 1000 for ``outlet`` (long enough to bridge large water bodies)
    and 100 for ``whitebox`` (the WhiteboxTools default).
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

Returns
-------
Raster
    New in-memory raster; save it with ``.write(path)``.

Authors
-------
B. Gailleton
S. M. Mudd

```python
def breach_least_cost(raster: lsdtt3.Raster, max_distance: int | None = None, max_cost: float = 3.4028234663852886e+38, prefer_shortest_path: bool = False, flat_increment: float = 0.0, fill_remainder: bool = True, n_threads: int = 0, edge_policy: str = 'outlet', fill_min_slope: float = 9.999999747378752e-05) -> lsdtt3.Raster
```

---

## `carve_depth`

Return the non-negative elevation removed by depression breaching.

Parameters
----------
raster : Raster
    Elevation raster (m).
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
Raster
    ``raster - breach_depressions(raster)`` in metres.

Authors
-------
B. Gailleton

```python
def carve_depth(raster: lsdtt3.Raster, n_threads: int = 0) -> lsdtt3.Raster
```

---

## `classify_flat_surfaces`

Classify labelled flats as floodplain or terrace from HAND and channel IDs.

``channel_source_ids`` is a co-registered Raster containing non-negative
source IDs on in-scope channel pixels and negative values elsewhere.

Parameters
----------
labelled_flat_mask : Raster
    Flat labels (``FlatSurfaceResult.labelled_flat_mask``).
height_above_channel : Raster
    HAND (m) on the same grid (``FlowInfo.height_above_channel``).
channel_source_ids : Raster
    Source ID of each in-scope channel pixel, negative elsewhere.
floodplain_max_height_metres : float, default=1.0
    A channel-adjacent label is floodplain when its median HAND (m) is at
    most this.
terrace_min_height_metres : float, default=0.0
    A label linked to a channel is terrace when its median HAND (m) is at
    least this (and it is not floodplain).
channel_adjacency_radius_metres : float, default=0.0
    A label is channel-adjacent when a pixel lies within this distance (m)
    of an in-scope channel; 0 means touching (8-neighbour).
terrace_max_distance_to_channel_metres : float, default=0.0
    A terrace must have a pixel within this distance (m) of a channel,
    excluding upland flats; 0 means no limit.
enable_leak_split : bool, default=True
    Split a label whose HAND distribution is bimodal (a floodplain leaking
    into a terrace across a broken riser) at the antimode.
enable_hand_banding : bool, default=False
    Split labels into reach-relative HAND bands along the channel before
    classifying.
donor_flat : array-like of int or None, default=None
    HAND banding only: flat index of the channel pixel each pixel drains to
    (``FlowInfo.height_above_channel_with_donor()["channel_donor"]``).
channel_flow_distance : Raster or None, default=None
    HAND banding only: flow distance (m) used to place donors along the
    channel.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
FlatSurfaceClassification

Authors
-------
B. Gailleton

```python
def classify_flat_surfaces(labelled_flat_mask: lsdtt3.Raster, height_above_channel: lsdtt3.Raster, channel_source_ids: lsdtt3.Raster, floodplain_max_height_metres: float = 1.0, terrace_min_height_metres: float = 0.0, channel_adjacency_radius_metres: float = 0.0, terrace_max_distance_to_channel_metres: float = 0.0, enable_leak_split: bool = True, enable_hand_banding: bool = False, donor_flat: object | None = None, channel_flow_distance: object | None = None, n_threads: int = 0) -> lsdtt3.FlatSurfaceClassification
```

---

## `classify_flat_surfaces_with_scope`

Classify flat surfaces using a ``ChannelScope``.

HAND, its nearest-channel donors, flow distance, channel mask, and source IDs
are derived from the scope's retained D8 flow topology.

Parameters
----------
labelled_flat_mask : Raster
    Flat labels (``FlatSurfaceResult.labelled_flat_mask``).
dem : Raster
    Elevation (m) on the scope's flow grid.
scope : ChannelScope
    In-scope channels (``FlowInfo.channel_scope``).
floodplain_max_height_metres : float, default=1.0
    A channel-adjacent label is floodplain when its median HAND (m) is at
    most this.
terrace_min_height_metres : float, default=0.0
    A label linked to a channel is terrace when its median HAND (m) is at
    least this (and it is not floodplain).
channel_adjacency_radius_metres : float, default=0.0
    A label is channel-adjacent when a pixel lies within this distance (m)
    of an in-scope channel; 0 means touching (8-neighbour).
terrace_max_distance_to_channel_metres : float, default=0.0
    A terrace must have a pixel within this distance (m) of a channel,
    excluding upland flats; 0 means no limit.
enable_leak_split : bool, default=True
    Split a label whose HAND distribution is bimodal (a floodplain leaking
    into a terrace across a broken riser) at the antimode.
enable_hand_banding : bool, default=False
    Split labels into reach-relative HAND bands along the channel before
    classifying.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
FlatSurfaceClassification

Authors
-------
B. Gailleton

```python
def classify_flat_surfaces_with_scope(labelled_flat_mask: lsdtt3.Raster, dem: lsdtt3.Raster, scope: lsdtt3.ChannelScope, floodplain_max_height_metres: float = 1.0, terrace_min_height_metres: float = 0.0, channel_adjacency_radius_metres: float = 0.0, terrace_max_distance_to_channel_metres: float = 0.0, enable_leak_split: bool = True, enable_hand_banding: bool = False, n_threads: int = 0) -> lsdtt3.FlatSurfaceClassification
```

---

## `classify_surface_shape`

Compute polyfit principal curvatures, then classify surface shape.

This function has 2 call forms:

1. ``classify_surface_shape(raster: Raster, lambda_metres: float = 0.0, z_factor: float = 1.0, n_threads: int = 0, kt: float = 0.0) -> dict``

   Compute polyfit principal curvatures, then classify surface shape.

2. ``classify_surface_shape(maximum_curvature: Raster, minimum_curvature: Raster, kt: float = 0.0) -> dict``

   Classify surface shape from existing principal-curvature rasters.

Parameters
----------
raster : Raster
    Elevation raster (m) in a projected metre CRS.
lambda_metres : float, default=0.0
    Polyfit length scale (m), as in ``polyfit_metrics``; 0 uses three cells.
z_factor : float, default=1.0
    Elevation scale factor.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.
kt : float, default=0.0
    Call form 1: Curvature threshold (1/m): principal curvatures with magnitude <= kt
    count as zero (flat in that direction).
    Call form 2: Curvature threshold (1/m) below which a curvature counts as zero.
maximum_curvature, minimum_curvature : Raster
    Principal curvatures k_max and k_min (1/m), for example from
    ``polyfit_metrics`` or ``fundamental_form_principal_curvatures``.

Returns
-------
dict[str, Raster]
    Call form 1: ``"shape_class"``: SMAP class codes -4 perfect saddle, -3 basin
    (bowl), -2 synformal saddle, -1 synform (valley), 0 plane, 1 antiform
    (ridge), 2 antiformal saddle, 3 peak (dome); ``"gaussian_curvature"``
    (1/m^2) and ``"mean_curvature"`` (1/m) after the kt threshold.
dict[str, Raster]
    Call form 2: ``"shape_class"`` (codes as in the first call form),
    ``"gaussian_curvature"`` (1/m^2) and ``"mean_curvature"`` (1/m).

Authors
-------
B. Gailleton

```python
def classify_surface_shape(raster: lsdtt3.Raster, lambda_metres: float = 0.0, z_factor: float = 1.0, n_threads: int = 0, kt: float = 0.0) -> dict
def classify_surface_shape(maximum_curvature: lsdtt3.Raster, minimum_curvature: lsdtt3.Raster, kt: float = 0.0) -> dict
```

---

## `clip_to_maximum_rectangle`

Crop to the largest all-valid axis-aligned rectangle.

Parameters
----------
raster : Raster
    Input raster, for example a reprojected DEM with nodata corners.

Returns
-------
Raster
    New, smaller raster with no nodata pixels and updated georeferencing.

Authors
-------
B. Gailleton

```python
def clip_to_maximum_rectangle(raster: lsdtt3.Raster) -> lsdtt3.Raster
```

---

## `clip_to_polygons`

Crop to polygon bounds and set pixels outside polygons to nodata.

Parameters
----------
raster : Raster
    Raster to clip.
polygons : geopandas.GeoDataFrame
    Geometries with a CRS equivalent to the raster's (compared with GDAL
    ``IsSame``); otherwise ``ValueError``. Reproject first with
    ``gdf.to_crs(raster.metadata.crs_wkt)``.
coverage : {"any", "majority", "full"}, default="majority"
    Pixel-inclusion rule, as in ``rasterize_polygons``.
crop : bool, default=True
    Crop to polygon bounds before masking. False retains the source extent.

Returns
-------
Raster
    Raster with pixels outside the polygons set to nodata.

Authors
-------
B. Gailleton

```python
def clip_to_polygons(raster: lsdtt3.Raster, polygons: object, coverage: str = 'majority', crop: bool = True) -> lsdtt3.Raster
```

---

## `coherence_diffusion`

Coherence-enhancing diffusion: smoothing along linear features.

Smooths along the local feature orientation from the structure tensor and
little across it, bridging and denoising ridges, valleys and scarps.
time_step must be in (0, 0.25]. A radius of 0 means one cell
(derivative_radius_metres) or three cells (integration_radius_metres).
contrast sets how selective the along-feature smoothing is; min_diffusivity
is the across-feature floor, in (0, 1). The same smoothing as
smooth_algorithm_choice: coherence_diffusion in lsdtt-raster-preprocessing.

Parameters
----------
raster : Raster
    Elevation raster (m).
n_iterations : int, default=10
    Number of diffusion iterations.
time_step : float, default=0.15
    Explicit time step, in (0, 0.25].
derivative_radius_metres : float, default=0.0
    Smoothing radius (m) of the gradients; 0 means one cell.
integration_radius_metres : float, default=0.0
    Radius (m) over which the structure tensor is averaged; 0 means three
    cells.
contrast : float, default=0.2
    How selective the along-feature smoothing is.
min_diffusivity : float, default=0.001
    Across-feature diffusivity floor, in (0, 1).
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
Raster
    Smoothed DEM (m).

Authors
-------
S. M. Mudd
B. Gailleton

```python
def coherence_diffusion(raster: lsdtt3.Raster, n_iterations: int = 10, time_step: float = 0.15, derivative_radius_metres: float = 0.0, integration_radius_metres: float = 0.0, contrast: float = 0.2, min_diffusivity: float = 0.001, n_threads: int = 0) -> lsdtt3.Raster
```

---

## `coregister_dem`

Coregister a DEM to a reference with the Nuth and Kaab (2011) method.

Iteratively fits the horizontal and vertical shift that best explains the
elevation differences as a function of aspect, and applies it to moving.
Both rasters must be on the same grid (use reproject_to_grid first). The same
coregistration as run_dem_coregistration in lsdtt-raster-calculator, with
its defaults.

Parameters
----------
reference : Raster
    Reference DEM (m).
moving : Raster
    DEM (m) to shift, on the same grid as ``reference``.
max_iterations : int, default=8
    Largest number of Nuth-Kaab iterations.
offset_threshold_pixels : float, default=0.001
    Stop when an iteration's horizontal step is below this (pixels).
aspect_bin_count : int, default=72
    Number of aspect bins in the fit.
min_samples : int, default=500
    Smallest number of valid samples required.
subsample_count : int, default=500000
    Largest number of pixels sampled per iteration.
apply_vertical_shift : bool, default=True
    Also remove the vertical offset dz.
max_horizontal_step_metres : float, default=0.0
    Cap (m) on each iteration's horizontal step; 0 means none.
initial_dx_metres, initial_dy_metres, initial_dz_metres : float, default=0.0
    Starting shift (m) in x, y and z.

Returns
-------
dict
    ``"raster"`` (moving, shifted), ``"dx_metres"``, ``"dy_metres"``,
    ``"dz_metres"`` (the total shift, m), ``"iterations"``, ``"converged"``
    and ``"diagnostics"`` (a pandas DataFrame with one row per iteration).

Authors
-------
S. M. Mudd
B. Gailleton

```python
def coregister_dem(reference: lsdtt3.Raster, moving: lsdtt3.Raster, max_iterations: int = 8, offset_threshold_pixels: float = 0.001, aspect_bin_count: int = 72, min_samples: int = 500, subsample_count: int = 500000, apply_vertical_shift: bool = True, max_horizontal_step_metres: float = 0.0, initial_dx_metres: float = 0.0, initial_dy_metres: float = 0.0, initial_dz_metres: float = 0.0) -> dict
```

---

## `detect_directional_periodicity`

Detect local regularly spaced parallel features with sliding-window FFTs.

Families of evenly spaced parallel lines (ditches, tracks, cultivation
ridges) give a sharp peak in the 2D power spectrum at 1 / spacing. A
window slides over the raster and the strongest peak in the spacing band
is recorded; each pixel keeps its best window.

Parameters
----------
raster : Raster
    Input raster, for example a DEM or a curvature response.
window_size_pixels : int
    Window side (pixels, >= 8); should span several of the largest spacings.
step_pixels : int
    Window stride (pixels, >= 1), typically half the window.
min_spacing_metres, max_spacing_metres : float
    Range of spacings (m) searched; the minimum must exceed 2 cell sizes and
    the maximum be below the window extent.
min_valid_fraction : float, default=0.6
    Windows with a smaller fraction of valid pixels are skipped.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
dict
    ``"strength"`` (peak prominence), ``"spacing_metres"`` and
    ``"orientation_degrees"`` (line orientation) Rasters.

Authors
-------
B. Gailleton

```python
def detect_directional_periodicity(raster: lsdtt3.Raster, window_size_pixels: int, step_pixels: int, min_spacing_metres: float, max_spacing_metres: float, min_valid_fraction: float = 0.6, n_threads: int = 0) -> dict
```

---

## `double_raster`

Double both raster dimensions over the same extent using the optimized path.

Parameters
----------
source : Raster
    Raster to refine; the cell size is halved.
method : str, default="linear"
    Sampling method (``"nearest"``, ``"linear"``, ``"cubic"``).

Returns
-------
Raster

Authors
-------
B. Gailleton

```python
def double_raster(source: lsdtt3.Raster, method: str = 'linear') -> lsdtt3.Raster
```

---

## `fill_depth`

Return the non-negative elevation added by pit filling.

Parameters
----------
raster : Raster
    Elevation raster (m).
min_slope : float, default=0.0001
    Gradient (m/m) imposed across filled flats: each filled pixel is raised
    at least ``min_slope * cell_size`` (``* sqrt(2)`` on a diagonal step)
    above the pixel it drains to, so every pixel has a downhill route. 0
    fills to a flat surface.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
Raster
    ``fill_pits(raster) - raster`` in metres (0 outside depressions).

Authors
-------
B. Gailleton

```python
def fill_depth(raster: lsdtt3.Raster, min_slope: float = 9.999999747378752e-05, n_threads: int = 0) -> lsdtt3.Raster
```

---

## `fill_interior_nodata`

Fill nodata regions that are not connected to the raster edge.

Holes enclosed by valid pixels are filled by interpolation; nodata regions
that touch the raster edge are left as nodata.

Parameters
----------
raster : Raster
    Input raster.
method : {"laplacian", "iterative"}, default="laplacian"
    ``"laplacian"`` solves Laplace's equation across each hole (smooth
    surface through its rim); ``"iterative"`` repeatedly averages valid
    neighbours.
max_iter : int, default=10000
    Maximum number of iterations (positive).
tolerance : float, default=1e-6
    Convergence tolerance on the change per iteration, in raster units
    (positive).
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
Raster
    New in-memory raster.

Authors
-------
B. Gailleton

```python
def fill_interior_nodata(raster: lsdtt3.Raster, method: str = 'laplacian', max_iter: int = 10000, tolerance: float = 9.999999974752427e-07, n_threads: int = 0) -> lsdtt3.Raster
```

---

## `fill_pits`

Fill depressions while imposing an optional minimum drainage slope.

Priority-Flood filling (Barnes et al. 2014): every pit and closed depression
is raised to its spill level, so all water can drain to the raster edge or
to nodata. The raster edge and pixels next to nodata act as outlets. Use the
result to build a ``FlowInfo``.

Parameters
----------
raster : Raster
    Elevation raster (m).
min_slope : float, default=0.0001
    Gradient (m/m) imposed across filled flats: each filled pixel is raised
    at least ``min_slope * cell_size`` (``* sqrt(2)`` on a diagonal step)
    above the pixel it drains to, so every pixel has a downhill route. 0
    fills to a flat surface.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.
priority_flood_backend : {"binary", "radix"}, default="binary"
    Priority-queue implementation: ``"binary"`` heap or ``"radix"``
    (hybrid radix heap). Both give the same result; only speed differs.

Returns
-------
Raster
    Filled DEM (m).

Authors
-------
B. Gailleton

```python
def fill_pits(raster: lsdtt3.Raster, min_slope: float = 9.999999747378752e-05, n_threads: int = 0, priority_flood_backend: str = 'binary') -> lsdtt3.Raster
```

---

## `fill_pits_boundary_aware`

Fill pits with NESW open, closed, periodic, or base-level boundaries.

Closed (``n``) edges are not outlets; periodic (``p``) edges wrap, so
depressions spilling across the seam are filled to their true spill level.
Fill with the same ``boundary_conditions`` you will give ``FlowInfo``.

Parameters
----------
raster : Raster
    Elevation raster (m).
boundary_conditions : str, default='oooo'
    Four NESW codes ``o``, ``n``, ``p``, ``b`` (see
    ``make_boundary_conditions``). Both axes cannot be periodic; with a
    periodic axis the other edges must be open.
min_slope : float, default=0.0001
    Gradient (m/m) imposed across filled flats: each filled pixel is raised
    at least ``min_slope * cell_size`` (``* sqrt(2)`` on a diagonal step)
    above the pixel it drains to, so every pixel has a downhill route. 0
    fills to a flat surface.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.
priority_flood_backend : {"binary", "radix"}, default="binary"
    Priority-queue implementation: ``"binary"`` heap or ``"radix"``
    (hybrid radix heap). Both give the same result; only speed differs.

Returns
-------
Raster
    Filled DEM (m).

Authors
-------
B. Gailleton

```python
def fill_pits_boundary_aware(raster: lsdtt3.Raster, boundary_conditions: str = 'oooo', min_slope: float = 9.999999747378752e-05, n_threads: int = 0, priority_flood_backend: str = 'binary') -> lsdtt3.Raster
```

---

## `fill_pits_priority_flood`

Explicit Improved Priority-Flood; optionally use its radix-heap queue.

As ``fill_pits`` but always uses the Priority-Flood algorithm, even when
``min_slope`` is 0.

Parameters
----------
raster : Raster
    Elevation raster (m).
min_slope : float, default=0.0001
    Gradient (m/m) imposed across filled flats: each filled pixel is raised
    at least ``min_slope * cell_size`` (``* sqrt(2)`` on a diagonal step)
    above the pixel it drains to, so every pixel has a downhill route. 0
    fills to a flat surface.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.
use_radix_heap : bool, default=False
    Use the radix-heap priority queue instead of a binary heap.

Returns
-------
Raster
    Filled DEM (m).

Authors
-------
B. Gailleton

```python
def fill_pits_priority_flood(raster: lsdtt3.Raster, min_slope: float = 9.999999747378752e-05, n_threads: int = 0, use_radix_heap: bool = False) -> lsdtt3.Raster
```

---

## `fill_pits_with_outlets`

Fill pits using positive cells in outlet_mask as drainage outlets.

Only the selected pixels act as outlets (the raster edge does not, unless
selected), so a DEM can be filled towards chosen edges.

Parameters
----------
raster : Raster
    Elevation raster (m).
outlet_mask : Raster
    0/1 raster on the same grid; pixels > 0.5 are outlets (for example from
    ``make_outlet_mask``). Must select at least one valid pixel.
min_slope : float, default=0.0001
    Gradient (m/m) imposed across filled flats: each filled pixel is raised
    at least ``min_slope * cell_size`` (``* sqrt(2)`` on a diagonal step)
    above the pixel it drains to, so every pixel has a downhill route. 0
    fills to a flat surface.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.
priority_flood_backend : {"binary", "radix"}, default="binary"
    Priority-queue implementation: ``"binary"`` heap or ``"radix"``
    (hybrid radix heap). Both give the same result; only speed differs.

Returns
-------
Raster
    Filled DEM (m).

Authors
-------
B. Gailleton

```python
def fill_pits_with_outlets(raster: lsdtt3.Raster, outlet_mask: lsdtt3.Raster, min_slope: float = 9.999999747378752e-05, n_threads: int = 0, priority_flood_backend: str = 'binary') -> lsdtt3.Raster
```

---

## `fill_pits_with_sinks`

Fill pits while retaining positive cells in sinks as declared sinks.

This function has 2 call forms:

1. ``fill_pits_with_sinks(raster: Raster, sinks: Raster, min_slope: float = 0.0001, n_threads: int = 0, priority_flood_backend: str = 'binary') -> Raster``

   Fill pits while retaining positive cells in sinks as declared sinks.

   Depressions draining to a declared sink (moulin, swallow hole, playa) are
   not filled; use the same pixels as ``forced_outlets`` of ``FlowInfo``.

2. ``fill_pits_with_sinks(raster: Raster, sinks: object, min_slope: float = 0.0001, n_threads: int = 0, priority_flood_backend: str = 'binary') -> Raster``

   Fill pits while retaining zero-based (row, col) sinks.

Parameters
----------
raster : Raster
    Elevation raster (m).
sinks : Raster or array-like of int
    Call form 1: 0/1 raster on the same grid; pixels > 0.5 are sinks.
    Call form 2: N x 2 array of ``(row, col)`` sink pixels (valid pixels).
min_slope : float, default=0.0001
    Gradient (m/m) imposed across filled flats: each filled pixel is raised
    at least ``min_slope * cell_size`` (``* sqrt(2)`` on a diagonal step)
    above the pixel it drains to, so every pixel has a downhill route. 0
    fills to a flat surface.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.
priority_flood_backend : {"binary", "radix"}, default="binary"
    Priority-queue implementation: ``"binary"`` heap or ``"radix"``
    (hybrid radix heap). Both give the same result; only speed differs.

Returns
-------
Raster
    Filled DEM (m).

Authors
-------
B. Gailleton

```python
def fill_pits_with_sinks(raster: lsdtt3.Raster, sinks: lsdtt3.Raster, min_slope: float = 9.999999747378752e-05, n_threads: int = 0, priority_flood_backend: str = 'binary') -> lsdtt3.Raster
def fill_pits_with_sinks(raster: lsdtt3.Raster, sinks: object, min_slope: float = 9.999999747378752e-05, n_threads: int = 0, priority_flood_backend: str = 'binary') -> lsdtt3.Raster
```

---

## `find_flat_surfaces`

Detect and label flat surfaces from texture and local relief.

Each enabled roughness metric is computed in a moving window and thresholded
at a percentile of its own distribution (low roughness = flat); the metric
masks are combined, local relief must also be below
``relief_threshold_metres``, and connected flats are labelled.

Parameters
----------
dem : Raster
    Elevation raster (m).
use_madk2 : bool, default=True
    Use MADk2 roughness (median absolute deviation of second differences;
    insensitive to uniform slope).
madk2_window_radius_metres : float, default=100.0
    MADk2 window radius (m).
madk2_threshold_percentile : float, default=25.0
    Pixels below this percentile of MADk2 count as flat.
use_residual_iqr : bool, default=True
    Use the interquartile range of residuals to a smooth reference surface.
iqr_method : {"fftw_highpass", "wavelet_quadratic", "moving_window_mean"}, default="fftw_highpass"
    How the reference surface is removed (the first two need FFTW).
iqr_fit_lambda_metres : float, default=500.0
    Reference-surface wavelength (m) for ``fftw_highpass`` and
    ``wavelet_quadratic``.
iqr_detrend_radius_metres : float, default=500.0
    Moving-window radius (m) for ``moving_window_mean``.
iqr_window_radius_metres : float, default=50.0
    Radius (m) of the window in which the IQR is computed.
iqr_window_shape : {"circular", "square"}, default="circular"
    Shape of that window.
iqr_threshold_percentile : float, default=25.0
    Pixels below this percentile of the IQR count as flat.
use_nvd : bool, default=False
    Use normal-vector dispersion (needs FFTW).
nvd_polyfit_lambda : float, default=10.0
    Wavelength (m) of the surface fit giving the normals.
nvd_window_radius_metres : float, default=100.0
    NVD window radius (m).
nvd_window_shape : {"circular", "square"}, default="circular"
    NVD window shape.
nvd_threshold_percentile : float, default=25.0
    Pixels below this percentile of NVD count as flat.
require_all_metrics_pass : bool, default=True
    True: a pixel must pass every enabled metric (AND); False: any (OR).
directional_close_line_length_pixels : int, default=0
    Length (pixels) of the line used to close small gaps in narrow non-flat
    corridors (channels, scarps) before labelling; below 3 disables it.
directional_close_n_angles : int, default=18
    Number of line orientations tried (18 = 10-degree steps).
min_component_pixels : int, default=100
    Smallest labelled flat (pixels).
relief_window_side_pixels : int, default=5
    Side (pixels, odd, >= 3) of the square local-relief window.
relief_threshold_metres : float, default=1.0
    Pixels with local relief (max - min, m) at or above this are not flat.
remove_risers : bool, default=False
    Fit a plane to each flat and remove pixels with high residual relief
    (risers between terrace levels), then relabel.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
FlatSurfaceResult
    Roughness rasters, per-metric and combined masks, and labelled flats.

Authors
-------
B. Gailleton
S. M. Mudd

**Deprecated compatibility names:**
- `relief_window_pixels` -> `relief_window_side_pixels`; earliest removal: following minor.

```python
def find_flat_surfaces(dem: lsdtt3.Raster, use_madk2: bool = True, madk2_window_radius_metres: float = 100.0, madk2_threshold_percentile: float = 25.0, use_residual_iqr: bool = True, iqr_method: str = 'fftw_highpass', iqr_fit_lambda_metres: float = 500.0, iqr_detrend_radius_metres: float = 500.0, iqr_window_radius_metres: float = 50.0, iqr_window_shape: str = 'circular', iqr_threshold_percentile: float = 25.0, use_nvd: bool = False, nvd_polyfit_lambda: float = 10.0, nvd_window_radius_metres: float = 100.0, nvd_window_shape: str = 'circular', nvd_threshold_percentile: float = 25.0, require_all_metrics_pass: bool = True, directional_close_line_length_pixels: int = 0, directional_close_n_angles: int = 18, min_component_pixels: int = 100, relief_window_side_pixels: int = 5, relief_threshold_metres: float = 1.0, remove_risers: bool = False, n_threads: int = 0) -> lsdtt3.FlatSurfaceResult
```

---

## `frangi_filter`

Run a multiscale Frangi ridge or valley filter.

The Hessian of the Gaussian-smoothed raster is evaluated at ``n_scales``
smoothing widths and the Frangi line-likeness response is kept at its
maximum over the scales.

Parameters
----------
raster : Raster
    Input raster, usually elevation (m).
sigma_min_metres, sigma_max_metres : float, default=1.0, 5.0
    Smallest and largest Gaussian smoothing width (m).
n_scales : int, default=4
    Number of scales.
beta : float, default=0.5
    Frangi blob-suppression parameter (positive); smaller values penalise
    blob-like (non-linear) features more.
c : float, default=1.0
    Frangi structure-sensitivity parameter, in the units of the second
    derivatives (1/m for a DEM).
target : {"ridges", "valleys"}, default="ridges"
    Enhance ridges (convex lines) or valleys (concave lines).
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
FrangiResult

Authors
-------
B. Gailleton
S. M. Mudd

```python
def frangi_filter(raster: lsdtt3.Raster, sigma_min_metres: float = 1.0, sigma_max_metres: float = 5.0, n_scales: int = 4, beta: float = 0.5, c: float = 1.0, target: str = 'ridges', n_threads: int = 0) -> lsdtt3.FrangiResult
```

---

## `fundamental_form_principal_curvatures`

Compute differential-geometry maximum and minimum principal curvatures.

The DEM is treated as the surface (x, y, z) and the curvatures come from
its first and second fundamental forms, built from central differences.

Parameters
----------
raster : Raster
    Elevation raster (m) in a projected metre CRS.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
dict[str, Raster]
    ``"maximum_curvature"`` and ``"minimum_curvature"`` (1/m).

Authors
-------
B. Gailleton

```python
def fundamental_form_principal_curvatures(raster: lsdtt3.Raster, n_threads: int = 0) -> dict
```

---

## `halve_raster`

Halve both even raster dimensions over the same extent using the optimized path.

Parameters
----------
source : Raster
    Raster with even dimensions; the cell size is doubled.
method : str, default="conservative_mean"
    Resampling method, as in ``resample_raster``.

Returns
-------
Raster

Authors
-------
B. Gailleton

```python
def halve_raster(source: lsdtt3.Raster, method: str = 'conservative_mean') -> lsdtt3.Raster
```

---

## `highpass_filter`

Return the short-wavelength DEM component.

Parameters
----------
raster : Raster
    Elevation raster (m) in a projected metre CRS.
cutoff_metres : float
    Wavelengths (m) shorter than this are kept.
window : {"hann", "hamming", "none"}, default="hann"
    Taper applied before the FFT.
clip_to_valid_rectangle : bool, default=False
    First crop to the largest all-valid rectangle (the result is then
    smaller than the input).
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
Raster

Authors
-------
B. Gailleton

```python
def highpass_filter(raster: lsdtt3.Raster, cutoff_metres: float, window: str = 'hann', clip_to_valid_rectangle: bool = False, n_threads: int = 0) -> lsdtt3.Raster
```

---

## `hillshade`

Compute a single-illumination hillshade in the range 0 to 255.

Lambertian shading from Zevenbergen and Thorne (1987) derivatives. Edge
pixels and pixels next to nodata are nodata.

Parameters
----------
raster : Raster
    Elevation raster (m) in a projected metre CRS.
azimuth : float, default=315.0
    Direction the light comes from, in compass degrees (0 = north,
    clockwise); 315 is north-west.
altitude : float, default=45.0
    Sun elevation above the horizon, in degrees (0-90).
z_factor : float, default=1.0
    Vertical exaggeration multiplied into elevations; use it to exaggerate
    relief or to convert feet to metres (0.3048).

Returns
-------
Raster
    Shading values 0 (dark) to 255 (bright).

Authors
-------
B. Gailleton

```python
def hillshade(raster: lsdtt3.Raster, azimuth: float = 315.0, altitude: float = 45.0, z_factor: float = 1.0) -> lsdtt3.Raster
```

---

## `interpolate_to_nearest_power_of_two`

Resize to lower, higher, or nearest power-of-two dimensions.

The full extent is preserved. Rasters whose two dimension scale factors differ
cannot be represented with square cells and are rejected.

Parameters
----------
source : Raster
    Raster to resize.
direction : {"lower", "higher", "nearest"}, default="nearest"
    Which power of two to resize each dimension to.
method : str, default="linear"
    Sampling or aggregation method, as in ``resample_raster``.

Returns
-------
Raster

Authors
-------
B. Gailleton

```python
def interpolate_to_nearest_power_of_two(source: lsdtt3.Raster, direction: str = 'nearest', method: str = 'linear') -> lsdtt3.Raster
```

---

## `invert_binary_mask`

Invert a binary raster mask while preserving nodata.

Parameters
----------
raster : Raster
    0/1 raster; values > 0.5 count as 1.

Returns
-------
Raster
    0 where the input is > 0.5, 1 elsewhere; nodata stays nodata.

Authors
-------
B. Gailleton

```python
def invert_binary_mask(raster: lsdtt3.Raster) -> lsdtt3.Raster
```

---

## `least_cost_carve_depth`

Return the non-negative elevation removed by least-cost breaching (parameters as breach_least_cost).

Parameters
----------
raster : Raster
    Elevation raster (m).
max_distance, max_cost, prefer_shortest_path, flat_increment, fill_remainder, n_threads, edge_policy, fill_min_slope
    As in ``breach_least_cost``.

Returns
-------
Raster
    ``raster - breach_least_cost(raster, ...)`` in metres.

Authors
-------
B. Gailleton
S. M. Mudd

```python
def least_cost_carve_depth(raster: lsdtt3.Raster, max_distance: int | None = None, max_cost: float = 3.4028234663852886e+38, prefer_shortest_path: bool = False, flat_increment: float = 0.0, fill_remainder: bool = True, n_threads: int = 0, edge_policy: str = 'outlet', fill_min_slope: float = 9.999999747378752e-05) -> lsdtt3.Raster
```

---

## `lowpass_filter`

Return the long-wavelength DEM component.

Parameters
----------
raster : Raster
    Elevation raster (m) in a projected metre CRS.
cutoff_metres : float
    Wavelengths (m) longer than this are kept.
window : {"hann", "hamming", "none"}, default="hann"
    Taper applied before the FFT.
clip_to_valid_rectangle : bool, default=False
    First crop to the largest all-valid rectangle (the result is then
    smaller than the input).
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
Raster

Authors
-------
B. Gailleton

```python
def lowpass_filter(raster: lsdtt3.Raster, cutoff_metres: float, window: str = 'hann', clip_to_valid_rectangle: bool = False, n_threads: int = 0) -> lsdtt3.Raster
```

---

## `make_boundary_conditions`

Create an LSDTopoTools NESW boundary-condition string.

Parameters
----------
north, east, south, west : {"open", "closed", "periodic", "base_level"}
    Boundary treatment of each edge (``"no_flux"`` and the one-letter codes
    ``o``, ``n``, ``p``, ``b`` are also accepted). ``open``: flow may leave
    across the edge; ``closed``: no flow leaves; ``periodic``: the edge
    wraps to the opposite one (give both opposite edges); ``base_level``:
    every edge pixel is an outlet.

Returns
-------
str
    Four-character NESW code such as ``"onon"``, for ``FlowInfo``,
    ``fill_pits_boundary_aware`` or ``run_graphflood``.

Authors
-------
B. Gailleton

```python
def make_boundary_conditions(north: str = 'open', east: str = 'open', south: str = 'open', west: str = 'open') -> str
```

---

## `make_outlet_mask`

Create a valid binary mask with selected raster edges set to one.

Parameters
----------
raster : Raster
    Raster whose grid the mask uses.
north, east, south, west : bool, default=False
    Set the pixels of that edge (row 0, last column, last row, column 0) to 1.

Returns
-------
Raster
    0/1 raster with every pixel valid, for ``fill_pits_with_outlets`` or as
    an outlet mask for ``run_graphflood``.

Authors
-------
B. Gailleton

```python
def make_outlet_mask(raster: lsdtt3.Raster, north: bool = False, east: bool = False, south: bool = False, west: bool = False) -> lsdtt3.Raster
```

---

## `mask_by_binary`

Set values to nodata where keep_mask is false or nodata.

Parameters
----------
values : Raster
    Raster whose values are kept.
keep_mask : Raster
    0/1 raster on the same grid; values are kept where it is > 0.5, for
    example ``BasinCollection.mask(0)`` or ``threshold_above(...)``.

Returns
-------
Raster
    ``values`` where ``keep_mask`` > 0.5, nodata elsewhere.

Authors
-------
B. Gailleton

```python
def mask_by_binary(values: lsdtt3.Raster, keep_mask: lsdtt3.Raster) -> lsdtt3.Raster
```

---

## `mask_from_points`

Create a binary mask from zero-based ``(row, col)`` point coordinates.

Parameters
----------
raster : Raster
    Raster whose grid the mask uses.
points : array-like, shape (n, 2)
    Integer ``(row, col)`` pixels (row 0 is the northern edge); each must be
    a valid pixel of ``raster``.

Returns
-------
Raster
    0/1 raster, 1 at the points, every pixel valid.

Authors
-------
B. Gailleton

```python
def mask_from_points(raster: lsdtt3.Raster, points: object) -> lsdtt3.Raster
```

---

## `mask_to_nodata`

Set valid pixels below (or above) a threshold to nodata.

With ``below=True`` pixels with value ``<= threshold`` become nodata; with
``below=False`` pixels with value ``> threshold`` become nodata. A common
use is masking the sea in a Copernicus DEM, which stores sea as exactly
0 m: ``mask_to_nodata(dem, 0.0)``.

Parameters
----------
raster : Raster
    Input raster.
threshold : float
    Threshold in the units of the raster (m for a DEM).
below : bool, default=True
    True masks values ``<= threshold``; False masks values ``> threshold``.

Returns
-------
Raster
    New in-memory raster with the selected pixels set to nodata.

Authors
-------
B. Gailleton

```python
def mask_to_nodata(raster: lsdtt3.Raster, threshold: float, below: bool = True) -> lsdtt3.Raster
```

---

## `masked_surface_recovery`

Recover a smooth surface from selected pixels.

The DEM is detrended, zeroed outside ``mask`` and low-pass filtered, so a
smooth long-wavelength surface is interpolated under the masked-out area
(for example the piedmont surface beneath a mountain range).

Parameters
----------
raster : Raster
    Elevation raster (m).
mask : Raster
    0/1 raster; pixels > 0.5 are used.
cutoff_metres : float
    Low-pass cut-off wavelength (m).
window : {"hann", "hamming", "none"}, default="hann"
    Taper applied before the FFT.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
Raster
    Smooth surface (m), valid wherever the DEM is valid.

Authors
-------
B. Gailleton

```python
def masked_surface_recovery(raster: lsdtt3.Raster, mask: lsdtt3.Raster, cutoff_metres: float, window: str = 'hann', n_threads: int = 0) -> lsdtt3.Raster
```

---

## `mosaic`

Mosaic loaded rasters into their union extent.

Parameters
----------
rasters : sequence of Raster
    Rasters with identical CRS, equal cell size and aligned pixel edges
    (otherwise ``ValueError``). Later valid pixels take precedence in
    overlaps. Use ``reproject_to_grid`` first to align rasters.

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

Compute an averaged multi-illumination hillshade in the range 0 to 255.

The mean of ``n_directions`` hillshades lit from evenly spaced azimuths
``start_azimuth + k * 360 / n_directions``; this reduces the directional
bias of a single light source.

Parameters
----------
raster : Raster
    Elevation raster (m) in a projected metre CRS.
n_directions : int, default=8
    Number of light directions.
altitude : float, default=45.0
    Sun elevation above the horizon, in degrees.
z_factor : float, default=1.0
    Vertical exaggeration multiplied into elevations.
start_azimuth : float, default=0.0
    First light direction in compass degrees (0 = north, clockwise).

Returns
-------
Raster
    Mean shading values 0 to 255.

Authors
-------
B. Gailleton

```python
def multidirectional_hillshade(raster: lsdtt3.Raster, n_directions: int = 8, altitude: float = 45.0, z_factor: float = 1.0, start_azimuth: float = 0.0) -> lsdtt3.Raster
```

---

## `multiply`

Multiply two grid-compatible rasters.

Parameters
----------
a, b : Raster
    Rasters on the same grid.

Returns
-------
Raster
    ``a * b``; nodata where either input is nodata.

Authors
-------
B. Gailleton

```python
def multiply(a: lsdtt3.Raster, b: lsdtt3.Raster) -> lsdtt3.Raster
```

---

## `normalize_unit_interval`

Linearly rescale valid values to the closed interval [0, 1].

Parameters
----------
raster : Raster
    Input raster; its minimum maps to 0 and maximum to 1.

Returns
-------
Raster

Authors
-------
B. Gailleton

```python
def normalize_unit_interval(raster: lsdtt3.Raster) -> lsdtt3.Raster
```

---

## `opentopography_datasets`

Return OpenTopography DEM catalogue records as dictionaries.

Parameters
----------
None

Returns
-------
list of dict
    One dict per dataset with keys ``demtype`` (the ``ot_dem_type`` value
    for ``fetch_dem``), ``name``, ``resolution``, ``coverage``,
    ``description``, ``citation`` and ``max_area_km2``.

Authors
-------
B. Gailleton

```python
def opentopography_datasets() -> list
```

---

## `pad_to_nearest_power_of_two`

Pad the north and east sides with nodata to the next power-of-two dimensions.

Parameters
----------
source : Raster
    Raster to pad.

Returns
-------
Raster

Authors
-------
B. Gailleton

```python
def pad_to_nearest_power_of_two(source: lsdtt3.Raster) -> lsdtt3.Raster
```

---

## `perona_malik_smoothing`

Edge-preserving Perona-Malik smoothing of a DEM.

Explicit anisotropic diffusion: a neighbour difference below the edge
threshold K is smoothed, one above it is preserved. K is the k_quantile of the
DEM's own neighbour-difference distribution (stencil_difference_magnitude),
or k_metres_per_cell when k_quantile is 0. time_step must be in (0, 0.25];
conductance is "exponential" or "quadratic". Nodata cells are untouched and
never diffused into. The same smoothing as smooth_algorithm_choice:
perona_malik in lsdtt-raster-preprocessing.

Parameters
----------
raster : Raster
    Elevation raster (m).
n_iterations : int, default=10
    Number of diffusion iterations.
time_step : float, default=0.2
    Explicit time step, in (0, 0.25].
conductance : {"exponential", "quadratic"}, default="exponential"
    Perona-Malik conductance function.
k_quantile : float, default=0.9
    Quantile (0-1) of the neighbour-difference distribution used as K; 0
    uses ``k_metres_per_cell``.
k_metres_per_cell : float, default=0.0
    Edge threshold K (m of elevation difference per cell) when
    ``k_quantile`` is 0.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
dict
    ``"raster"`` (smoothed DEM), ``"k_metres_per_cell"`` (the K used) and
    ``"n_iterations"``.

Authors
-------
S. M. Mudd
B. Gailleton

```python
def perona_malik_smoothing(raster: lsdtt3.Raster, n_iterations: int = 10, time_step: float = 0.2, conductance: str = 'exponential', k_quantile: float = 0.9, k_metres_per_cell: float = 0.0, n_threads: int = 0) -> dict
```

---

## `plan_curvature`

Compute 3-by-3 plan curvature; concave-up is positive.

Plan (contour) curvature from the Zevenbergen-Thorne 3 x 3 stencil:
positive where contours are concave (flow converges, as in hollows and
channels), negative where convex (flow diverges, as on ridge noses).

Parameters
----------
raster : Raster
    Elevation raster (m) in a projected metre CRS.

Returns
-------
Raster
    Plan curvature in 1/m; 0 on flat pixels; nodata at edges and next to
    nodata.

Authors
-------
B. Gailleton

```python
def plan_curvature(raster: lsdtt3.Raster) -> lsdtt3.Raster
```

---

## `polyfit_metrics`

Compute several quadratic-fit surface metrics in one shared fit.

A quadratic surface ``z = a x^2 + b y^2 + c x y + d x + e y + f`` is fitted
in a square window around every pixel and the metrics are derived from its
coefficients (Florinsky 2017 names). x is east and y is north. Pixels whose
window is incomplete (near edges or nodata) are nodata.

Valid metric names and units:

- ``"fitted_surface"``: fitted elevation (m)
- ``"dzdx"``, ``"dzdy"``: gradient components east and north (m/m)
- ``"d2zdx2"``, ``"d2zdxdy"``, ``"d2zdy2"``: second derivatives (1/m)
- ``"slope"``: gradient magnitude ``sqrt(dzdx**2 + dzdy**2)`` (m/m)
- ``"slope_degrees"``: the same slope as an angle (degrees)
- ``"aspect_degrees"``, ``"aspect_radians"``: downslope direction,
  clockwise from north
- ``"laplacian"``: ``d2zdx2 + d2zdy2`` (1/m)
- ``"hessian_determinant"`` (1/m^2)
- ``"mean_curvature"`` (1/m), ``"gaussian_curvature"`` (1/m^2)
- ``"maximum_curvature"``, ``"minimum_curvature"``: principal curvatures
  k_max and k_min (1/m)
- ``"curvedness"`` (1/m), ``"shape_index"`` (dimensionless)
- ``"profile_curvature"`` (vertical, k_v), ``"plan_curvature"``,
  ``"tangential_curvature"`` (horizontal, k_h) (1/m)
- ``"difference_curvature"`` (0.5 (k_v - k_h)),
  ``"horizontal_excess_curvature"`` (k_h - k_min),
  ``"vertical_excess_curvature"`` (k_v - k_min) (1/m)
- ``"accumulation_curvature"`` (k_h k_v), ``"ring_curvature"`` (1/m^2)
- ``"rotor"`` (1/m)

Parameters
----------
raster : Raster
    Input elevation raster (m) in a projected metre CRS.
metrics : list[str]
    Requested metric names from the list above; an unknown name raises
    ``ValueError``.
lambda_metres : float, default=0
    Fitting length scale in metres, roughly the window width: the window
    is ``2W + 1`` pixels square with ``W = max(1, round(lambda_metres /
    (2 * cell_size)))``. Zero uses ``3 * cell_size``. Larger values smooth
    more. Must not be less than the cell size.
z_factor : float, default=1
    Elevation scale factor applied before metric calculation (for example
    0.3048 for elevations in feet).
n_threads : int, default=0
    OpenMP thread count; zero selects the runtime default.
flat_slope_tolerance : float, default=1e-12
    Where the fitted gradient (m/m) is below this, the pixel is treated as
    too flat for direction-dependent metrics (aspect, plan, tangential and
    profile curvature, rotor).

Returns
-------
dict[str, Raster]
    One newly allocated raster for every requested metric, keyed by name.

Authors
-------
B. Gailleton

```python
def polyfit_metrics(raster: lsdtt3.Raster, metrics: collections.abc.Sequence[str], lambda_metres: float = 0.0, z_factor: float = 1.0, n_threads: int = 0, flat_slope_tolerance: float = 1e-12) -> dict
```

---

## `rasterize_polygons`

Rasterize GeoPandas polygon geometries onto a raster grid.

Parameters
----------
raster : Raster
    Raster whose grid is used.
polygons : geopandas.GeoDataFrame
    Geometries with a CRS equivalent to the raster's (compared with GDAL
    ``IsSame``); otherwise ``ValueError``. Reproject first with
    ``gdf.to_crs(raster.metadata.crs_wkt)``.
coverage : {"any", "majority", "full"}, default="majority"
    Pixel-inclusion rule: ``"any"`` a pixel touched by a polygon;
    ``"majority"`` the pixel centre lies inside (about 50% cover);
    ``"full"`` the pixel lies wholly inside.

Returns
-------
Raster
    0/1 mask (1 inside the polygons) with valid pixels throughout the grid.

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

Reproject a raster to an EPSG coordinate reference system.

Parameters
----------
raster : Raster
    Raster to warp.
epsg : int
    Target EPSG code, for example a UTM zone from ``utm_epsg_for_point``
    (32633 is UTM 33N).
cell_size : float, default=0
    Target cell size in target CRS units (m for UTM). A non-positive value
    lets GDAL choose it.
resampling : str, default="bilinear"
    GDAL resampling method (``"near"``, ``"bilinear"``, ``"cubic"``, ...).

Returns
-------
Raster
    New raster in the target CRS; pixels outside the source are nodata.

Authors
-------
B. Gailleton

```python
def reproject(raster: lsdtt3.Raster, epsg: int, cell_size: float = 0.0, resampling: str = 'bilinear') -> lsdtt3.Raster
```

---

## `reproject_to_grid`

Reproject and resample a raster onto an explicit target grid.

Use it to put a second raster exactly on a DEM's grid, for example
``reproject_to_grid(other, dem.metadata)``.

Parameters
----------
raster : Raster
    Raster to warp (any CRS).
target_metadata : RasterMeta
    Target extent, cell size, dimensions, nodata value and CRS (must have a
    CRS), usually another raster's ``.metadata``.
resampling : str, default="bilinear"
    GDAL resampling method: ``"near"``, ``"bilinear"``, ``"cubic"``,
    ``"cubicspline"``, ``"lanczos"``, ``"average"``, ``"mode"``, ... Use
    ``"near"`` for categorical rasters.

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

Resample an in-memory raster onto an explicit target grid without GDAL.

The source and target CRS must match. Sampling methods are ``nearest``,
``linear`` and ``cubic``. Aggregation methods are ``mean``, ``sum``, ``min``,
``max``, ``mode`` and the overlap-area-weighted ``conservative_mean``.

Parameters
----------
source : Raster
    Raster to resample.
target_metadata : RasterMeta
    Target grid in the same CRS.
method : str, default="linear"
    One of the methods above.

Returns
-------
Raster

Authors
-------
B. Gailleton

```python
def resample_raster(source: lsdtt3.Raster, target_metadata: lsdtt3.RasterMeta, method: str = 'linear') -> lsdtt3.Raster
```

---

## `sample_points`

Sample a raster at GeoPandas point geometries.

Each point takes the value of the pixel it falls in (no interpolation).

Parameters
----------
raster : Raster
    Raster to sample.
points : geopandas.GeoDataFrame
    Point geometries with a CRS exactly matching raster. The CRS is
    compared with GDAL ``IsSame``; reproject first with
    ``gdf.to_crs(raster.metadata.crs_wkt)``.
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

Multiply valid raster values by a scalar factor.

Parameters
----------
raster : Raster
    Input raster.
factor : float
    Multiplier, for example 0.3048 to convert feet to metres.

Returns
-------
Raster

Authors
-------
B. Gailleton

```python
def scale(raster: lsdtt3.Raster, factor: float) -> lsdtt3.Raster
```

---

## `set_nodata_rectangle`

Set a half-open row/column rectangle to nodata.

Pixels with ``row_start <= row < row_stop`` and ``col_start <= col <
col_stop`` (as in NumPy slicing; row 0 is the northern edge) become nodata.

Parameters
----------
raster : Raster
    Input raster.
row_start, row_stop : int
    Row range; non-empty and within the raster.
col_start, col_stop : int
    Column range; non-empty and within the raster.

Returns
-------
Raster
    New in-memory raster.

Authors
-------
B. Gailleton

```python
def set_nodata_rectangle(raster: lsdtt3.Raster, row_start: int, row_stop: int, col_start: int, col_stop: int) -> lsdtt3.Raster
```

---

## `slope`

Compute 3-by-3 gradient-magnitude slope in rise/run units.

Centred differences on the four cardinal neighbours:
``sqrt(dz/dx**2 + dz/dy**2)``. For a smoother, scale-dependent slope use
``polyfit_metrics(raster, ["slope"], lambda_metres=...)``.

Parameters
----------
raster : Raster
    Elevation raster (m) in a projected metre CRS.

Returns
-------
Raster
    Gradient in m/m (dimensionless rise over run; 1.0 is a 45-degree
    slope). Edge pixels and pixels next to nodata are nodata.

Authors
-------
B. Gailleton

```python
def slope(raster: lsdtt3.Raster) -> lsdtt3.Raster
```

---

## `spectral_analysis`

Run full FFT-based spectral analysis and return a SpectralResult.

The detrended, tapered DEM is Fourier transformed and its power spectral
density is compared with a Monte Carlo background of red-noise surfaces
(Perron et al. 2008). Needs FFTW (``spectral_available()``).

Parameters
----------
raster : Raster
    Elevation raster (m) in a projected metre CRS.
log_bin_width : float, default=0.05
    Width of the log10-frequency bins of the binned spectrum.
n_background_iterations : int, default=100
    Number of background noise surfaces averaged.
seed : int, default=-1
    Base random seed of the background; negative uses a fixed built-in
    sequence (still reproducible).
window : {"hann", "hamming", "none"}, default="hann"
    Taper applied before the FFT.
clip_to_valid_rectangle : bool, default=False
    First crop to the largest all-valid rectangle (the result is then
    smaller than the input).
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
SpectralResult

Authors
-------
B. Gailleton
S. M. Mudd

```python
def spectral_analysis(raster: lsdtt3.Raster, log_bin_width: float = 0.05, n_background_iterations: int = 100, seed: int = -1, window: str = 'hann', clip_to_valid_rectangle: bool = False, n_threads: int = 0) -> lsdtt3.SpectralResult
```

---

## `spectral_available`

Return whether this build includes FFTW spectral support.

Parameters
----------
None

Returns
-------
bool
    True if the FFT-based functions (spectral analysis, filters, Wiener
    channel heads) are available.

Authors
-------
B. Gailleton

```python
def spectral_available() -> bool
```

---

## `stencil_difference_magnitude`

Largest absolute elevation difference from each cell to a 4-neighbour.

In metres; nodata where the cell is nodata or has no valid neighbour. This is
the quantity Perona-Malik smoothing evaluates its conductance on, so its
distribution is what k_quantile is a quantile of.

Parameters
----------
raster : Raster
    Elevation raster (m).

Returns
-------
Raster
    Largest neighbour difference (m).

Authors
-------
S. M. Mudd
B. Gailleton

```python
def stencil_difference_magnitude(raster: lsdtt3.Raster) -> lsdtt3.Raster
```

---

## `subtract`

Subtract two grid-compatible rasters.

Parameters
----------
a, b : Raster
    Rasters on the same grid.

Returns
-------
Raster
    ``a - b`` (for example a DEM of difference); nodata where either input
    is nodata.

Authors
-------
B. Gailleton

```python
def subtract(a: lsdtt3.Raster, b: lsdtt3.Raster) -> lsdtt3.Raster
```

---

## `swath_distances`

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

Return a binary mask that is True where valid values exceed threshold.

Pixels with value ``>= threshold`` are 1 and others 0 (the threshold
itself counts as above).

Parameters
----------
raster : Raster
    Input raster.
threshold : float
    Threshold in the units of the raster.

Returns
-------
Raster
    0/1 raster; nodata where the input is nodata or not finite.

Authors
-------
B. Gailleton

```python
def threshold_above(raster: lsdtt3.Raster, threshold: float) -> lsdtt3.Raster
```

---

## `threshold_below`

Return a binary mask that is True where valid values are below threshold.

Pixels with value ``<= threshold`` are 1 and others 0 (the threshold
itself counts as below).

Parameters
----------
raster : Raster
    Input raster.
threshold : float
    Threshold in the units of the raster.

Returns
-------
Raster
    0/1 raster; nodata where the input is nodata or not finite.

Authors
-------
B. Gailleton

```python
def threshold_below(raster: lsdtt3.Raster, threshold: float) -> lsdtt3.Raster
```

---

## `topographic_wetness_index`

Compute ``log((A / cell_area) / slope)`` from drainage area and slope rasters.

Parameters
----------
drainage_area : Raster
    Drainage area (m^2), for example ``FlowInfo.drainage_area()`` or
    ``FlowGraph.drainage_area()``.
slope : Raster
    Gradient (m/m) on the same grid, for example ``slope(dem)``.

Returns
-------
Raster
    Natural-log wetness index; nodata where area or slope is not positive.

Authors
-------
B. Gailleton

```python
def topographic_wetness_index(drainage_area: lsdtt3.Raster, slope: lsdtt3.Raster) -> lsdtt3.Raster
```

---

## `total_curvature`

Compute moving-window total curvature at radius_metres.

The Laplacian ``z_xx + z_yy`` of a quadratic surface fitted by least
squares to the valid pixels in a square window.

Parameters
----------
raster : Raster
    Elevation raster (m) in a projected metre CRS.
radius_metres : float
    Window half-width in metres.

Returns
-------
Raster
    Total curvature in 1/m: positive for concave-up (valleys, hollows),
    negative for convex (ridges, hilltops).

Authors
-------
B. Gailleton

```python
def total_curvature(raster: lsdtt3.Raster, radius_metres: float) -> lsdtt3.Raster
```

---

## `vertical_units_to_metres`

Convert DEM elevations to metres.

vertical_units is "metres" (no change), "feet" (international foot, 0.3048 m)
or "us_survey_feet" (1200/3937 m). Never applied automatically: establish the
vertical unit from the data provider's metadata.

Parameters
----------
dem : Raster
    Elevation raster.
vertical_units : {"metres", "feet", "us_survey_feet"}
    Unit the elevations are stored in.

Returns
-------
Raster
    Elevations in metres.

Authors
-------
S. M. Mudd
B. Gailleton

```python
def vertical_units_to_metres(dem: lsdtt3.Raster, vertical_units: str) -> lsdtt3.Raster
```

---

## `wiener_filter`

Attenuate high-frequency white noise with an adaptive Wiener filter.

A power law fitted to the spectrum between the two signal wavelengths gives
the signal; the noise level comes from the near-Nyquist spectrum; each
frequency is weighted by signal / (signal + noise).

Parameters
----------
raster : Raster
    Elevation raster (m) in a projected metre CRS.
signal_wavelength_long_metres : float, default=1000.0
    Longest wavelength (m) of the band used to fit the signal.
signal_wavelength_short_metres : float, default=100.0
    Shortest wavelength (m) of that band.
window : {"hann", "hamming", "none"}, default="hann"
    Taper applied before the FFT.
clip_to_valid_rectangle : bool, default=False
    First crop to the largest all-valid rectangle (the result is then
    smaller than the input).
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
Raster
    Filtered DEM (m).

Authors
-------
B. Gailleton
S. M. Mudd

**Deprecated compatibility names:**
- `signal_wavelength_long_m` -> `signal_wavelength_long_metres`; earliest removal: following minor.
- `signal_wavelength_short_m` -> `signal_wavelength_short_metres`; earliest removal: following minor.

```python
def wiener_filter(raster: lsdtt3.Raster, signal_wavelength_long_metres: float = 1000.0, signal_wavelength_short_metres: float = 100.0, window: str = 'hann', clip_to_valid_rectangle: bool = False, n_threads: int = 0) -> lsdtt3.Raster
```
