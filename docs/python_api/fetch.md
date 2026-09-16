# Fetching And Geospatial Data API

This page is generated from the built Python runtime manifest. Do not edit signatures here.

## `OutletFetchResult`

In-memory DEM groups fetched around outlet points.

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
    In-memory DEM groups fetched around outlet points.

``rasters`` and ``groups`` share order: raster ``i`` is the DEM for group
``groups.iloc[i]``. ``outlets`` records the initial WGS84 envelopes used to
build those groups. ``basin_labels`` and ``snapped_outlets`` describe the
verified selected watersheds and their snapped channel locations.

    Authors
    -------
    B. Gailleton


### Members

#### `basin_labels`

(self) -> list

#### `groups`

(self) -> object

#### `outlets`

(self) -> object

#### `rasters`

(self) -> list

#### `snapped_outlets`

(self) -> object

---

## `VectorBounds`

WGS84 bounds derived from a vector layer.

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
    WGS84 bounds derived from a vector layer.

The input may contain polygons, or exactly two points defining opposite box
corners. Geometry is transformed to WGS84 before its bounds are calculated.

    Authors
    -------
    B. Gailleton


### Members

#### `bbox`

(self) -> tuple

#### `east`

(self) -> float

#### `geometry_kind`

(self) -> str

#### `layer_name`

(self) -> str

#### `message`

(self) -> str

#### `n_features`

(self) -> int

#### `north`

(self) -> float

#### `ok`

(self) -> bool

#### `source_crs`

(self) -> str

#### `south`

(self) -> float

#### `warning`

(self) -> str

#### `west`

(self) -> float

---

## `bounds_from_vector`

bounds_from_vector(path: str, layer: str = '') -> lsdtt3.VectorBounds

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
    bounds_from_vector(path: str, layer: str = '') -> lsdtt3.VectorBounds

Derive WGS84 ``(west, south, east, north)`` bounds from a vector file.

The layer must contain polygons or exactly two points. Use the result's
``bbox`` directly as keyword arguments to :func:`fetch_dem`.

    Authors
    -------
    B. Gailleton


```python
def bounds_from_vector(path: str, layer: str = '') -> lsdtt3.VectorBounds
```

---

## `cop30_tiles_for_bbox`

cop30_tiles_for_bbox(west: float, south: float, east: float, north: float) -> list[str]

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
cop30_tiles_for_bbox(west: float, south: float, east: float, north: float) -> list[str]

Authors
-------
B. Gailleton

```python
def cop30_tiles_for_bbox(west: float, south: float, east: float, north: float) -> list[str]
```

---

## `fetch_dem`

fetch_dem(west: float, south: float, east: float, north: float, source: str = 'cop30_aws', target_epsg: int = 0, grid_spacing: float = 30.0, resampling: str = 'cubic', nodata_value: float = -9999.0, size_pixels: int = 0, cache_dir: str = '', ot_dem_type: str = 'COP30', pc_collection: str = 'cop-dem-glo-30', pc_asset: str = 'data', http_timeout_seconds: int = 60, http_connect_timeout_seconds: int = 15, http_max_retry: int = 3) -> lsdtt3.Raster

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
    fetch_dem(west: float, south: float, east: float, north: float, source: str = 'cop30_aws', target_epsg: int = 0, grid_spacing: float = 30.0, resampling: str = 'cubic', nodata_value: float = -9999.0, size_pixels: int = 0, cache_dir: str = '', ot_dem_type: str = 'COP30', pc_collection: str = 'cop-dem-glo-30', pc_asset: str = 'data', http_timeout_seconds: int = 60, http_connect_timeout_seconds: int = 15, http_max_retry: int = 3) -> lsdtt3.Raster

Fetch a DEM for a WGS84 bounding box as an in-memory projected raster.

Parameters
----------
west, south, east, north : float
    WGS84 bounds in decimal degrees. Antimeridian-crossing requests are not
    supported: ``west`` must be less than ``east``.
source : {"cop30_aws", "opentopography", "planetary_computer"}, optional
    ``cop30_aws`` is a keyless global 30 m default. OpenTopography reads its
    key only from ``OT_API_KEY`` or ``OT_API_KEY_FILE``.
target_epsg : int, optional
    Output EPSG code. Zero selects UTM from the bounding-box centre.
grid_spacing : float, optional
    Output grid spacing in projected units.

cache_dir : str, optional
    Persistent cache directory. This is the only optional disk output; the
    returned raster itself is kept in memory.

Returns
-------
Raster
    A projected raster. Network, GDAL, and invalid-request failures raise
    ``RuntimeError``.

    Authors
    -------
    B. Gailleton


**Deprecated compatibility names:**
- `size_px` -> `size_pixels`; earliest removal: following minor.

```python
def fetch_dem(west: float, south: float, east: float, north: float, source: str = 'cop30_aws', target_epsg: int = 0, grid_spacing: float = 30.0, resampling: str = 'cubic', nodata_value: float = -9999.0, size_pixels: int = 0, cache_dir: str = '', ot_dem_type: str = 'COP30', pc_collection: str = 'cop-dem-glo-30', pc_asset: str = 'data', http_timeout_seconds: int = 60, http_connect_timeout_seconds: int = 15, http_max_retry: int = 3) -> lsdtt3.Raster
```

---

## `fetch_dem_for_outlets`

fetch_dem_for_outlets(outlets: object, source: str = 'cop30_aws', target_epsg: int = 0, grid_spacing: float = 30.0, initial_radius_km: float = 100.0, auto_split_disconnected_basins: bool = True, split_cluster_gap_km: float = 25.0, max_final_dem_area_km2: float = 2500.0, split_overlap_km: float = 2.0, max_basins_per_dem: int = 0, cache_dir: str = '', resampling: str = 'cubic', nodata_value: float = -9999.0, snapping_method: str = 'downstream', contributing_pixels_threshold: int = 1000, edge_buffer_pixels: int = 3, nodata_buffer_pixels: int = 1, tight_bbox_buffer_pixels: int = 3, verify_expand_factor: float = 1.25, verify_max_iters: int = 3, tight_bbox_refetch: bool = True, min_search_radius_pixels: int = 0, max_search_radius_pixels: int = 25, mask_outside_basins: bool = False) -> lsdtt3.OutletFetchResult

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
    fetch_dem_for_outlets(outlets: object, source: str = 'cop30_aws', target_epsg: int = 0, grid_spacing: float = 30.0, initial_radius_km: float = 100.0, auto_split_disconnected_basins: bool = True, split_cluster_gap_km: float = 25.0, max_final_dem_area_km2: float = 2500.0, split_overlap_km: float = 2.0, max_basins_per_dem: int = 0, cache_dir: str = '', resampling: str = 'cubic', nodata_value: float = -9999.0, snapping_method: str = 'downstream', contributing_pixels_threshold: int = 1000, edge_buffer_pixels: int = 3, nodata_buffer_pixels: int = 1, tight_bbox_buffer_pixels: int = 3, verify_expand_factor: float = 1.25, verify_max_iters: int = 3, tight_bbox_refetch: bool = True, min_search_radius_pixels: int = 0, max_search_radius_pixels: int = 25, mask_outside_basins: bool = False) -> lsdtt3.OutletFetchResult

Fetch one or more in-memory DEMs around WGS84 outlet points.

The initial radius around each point is converted to WGS84 bounds, nearby
outlets are grouped, and large groups are split before each group is fetched.
The input GeoDataFrame is reprojected to WGS84 internally. Each group is
verified with D8 basin delineation, expanded if a basin reaches an edge or
nodata, and optionally refetched to its tight basin envelope.

Parameters
----------
outlets : geopandas.GeoDataFrame
    Point geometries defining outlet locations, in any declared CRS.
initial_radius_km : float, optional
    Initial half-width around each outlet before grouping.
max_final_dem_area_km2 : float, optional
    Largest allowed grouped initial extent before splitting it.
snapping_method : {"none", "expanding_square", "downstream"}, optional
    Outlet snapping used before each basin delineation.
max_search_radius_pixels : int, optional
    Largest expanding-square snap radius when ``snapping_method`` is
    ``"expanding_square"``.
edge_buffer_pixels, nodata_buffer_pixels : int, optional
    Required basin clearance from the DEM edge and invalid cells.
tight_bbox_refetch : bool, optional
    Refetch a smaller final DEM around the verified basin footprint.
mask_outside_basins : bool, optional
    Set cells outside the selected basin union to nodata in returned DEMs.

Returns
-------
OutletFetchResult
    Grouped raster collection and copied outlet/group diagnostic tables.

    Authors
    -------
    B. Gailleton


```python
def fetch_dem_for_outlets(outlets: object, source: str = 'cop30_aws', target_epsg: int = 0, grid_spacing: float = 30.0, initial_radius_km: float = 100.0, auto_split_disconnected_basins: bool = True, split_cluster_gap_km: float = 25.0, max_final_dem_area_km2: float = 2500.0, split_overlap_km: float = 2.0, max_basins_per_dem: int = 0, cache_dir: str = '', resampling: str = 'cubic', nodata_value: float = -9999.0, snapping_method: str = 'downstream', contributing_pixels_threshold: int = 1000, edge_buffer_pixels: int = 3, nodata_buffer_pixels: int = 1, tight_bbox_buffer_pixels: int = 3, verify_expand_factor: float = 1.25, verify_max_iters: int = 3, tight_bbox_refetch: bool = True, min_search_radius_pixels: int = 0, max_search_radius_pixels: int = 25, mask_outside_basins: bool = False) -> lsdtt3.OutletFetchResult
```

---

## `gdal_data_description`

gdal_data_description() -> str

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
    gdal_data_description() -> str

As proj_data_description(), for GDAL's own data directory.

    Authors
    -------
    B. Gailleton


```python
def gdal_data_description() -> str
```

---

## `gdal_version`

gdal_version() -> str

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
    gdal_version() -> str

GDAL release name of the library actually loaded, e.g. "3.8.4".

    Authors
    -------
    B. Gailleton


```python
def gdal_version() -> str
```

---

## `http_available`

http_available() -> bool

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
    http_available() -> bool

Return whether the bundled GDAL has HTTP support.

    Authors
    -------
    B. Gailleton


```python
def http_available() -> bool
```

---

## `proj_data_description`

proj_data_description() -> str

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
    proj_data_description() -> str

Where PROJ is getting proj.db, as a human-readable string.

One of ``"bundled: <dir>"`` (lsdtt3 pinned it), ``"PROJ_DATA=<dir>"``,
``"PROJ_LIB=<dir>"``, or ``"PROJ built-in search"``.

Put this in any bug report about projections or coordinates that look subtly
wrong -- a mismatched proj.db usually returns plausible values rather than
failing.

    Authors
    -------
    B. Gailleton


```python
def proj_data_description() -> str
```

---

## `proj_version`

proj_version() -> str

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
    proj_version() -> str

PROJ version actually in use at runtime, e.g. "9.3.1".

Reported by the loaded library rather than by the headers this was compiled
against -- that distinction is the whole point.

    Authors
    -------
    B. Gailleton


```python
def proj_version() -> str
```

---

## `utm_epsg_for_point`

utm_epsg_for_point(longitude: float, latitude: float) -> int

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
    utm_epsg_for_point(longitude: float, latitude: float) -> int

Return the UTM EPSG code containing a WGS84 point.

    Authors
    -------
    B. Gailleton


```python
def utm_epsg_for_point(longitude: float, latitude: float) -> int
```
