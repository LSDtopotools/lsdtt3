# Fetching And Geospatial Data API

This page is generated from the built Python runtime manifest. Do not edit signatures here.

## `OutletFetchResult`

In-memory DEM groups fetched around outlet points.

``rasters`` and ``groups`` share order: raster ``i`` is the DEM for group
``groups.iloc[i]``. ``outlets`` records the initial WGS84 envelopes used to
build those groups. ``basin_labels`` and ``snapped_outlets`` describe the
verified selected watersheds and their snapped channel locations.

Returned by ``fetch_dem_for_outlets()``; not constructible directly.

Attributes
----------
basin_labels : list
groups : pandas.DataFrame
    One row per DEM: group_id, west, south, east, north (WGS84), area_km2, outlet_ids, basins_verified, verification_message.
outlets : pandas.DataFrame
    Initial outlet envelopes: outlet_id, longitude, latitude, status, iterations, west, south, east, north, message.
rasters : list
snapped_outlets : pandas.DataFrame
    Snapped outlets: group_id, outlet_id, row, col, x, y (pixel centre in the group raster's CRS).

Authors
-------
B. Gailleton

### Members

#### `basin_labels`

(self) -> list

#### `groups`

One row per DEM: group_id, west, south, east, north (WGS84), area_km2, outlet_ids, basins_verified, verification_message.

#### `outlets`

Initial outlet envelopes: outlet_id, longitude, latitude, status, iterations, west, south, east, north, message.

#### `rasters`

(self) -> list

#### `snapped_outlets`

Snapped outlets: group_id, outlet_id, row, col, x, y (pixel centre in the group raster's CRS).

---

## `VectorBounds`

WGS84 bounds derived from a vector layer.

The input may contain polygons, or exactly two points defining opposite box
corners. Geometry is transformed to WGS84 before its bounds are calculated.

Returned by ``bounds_from_vector()``; not constructible directly.

Attributes
----------
bbox : tuple
east : float
geometry_kind : str
layer_name : str
message : str
n_features : int
north : float
ok : bool
source_crs : str
south : float
warning : str
west : float

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

Derive WGS84 ``(west, south, east, north)`` bounds from a vector file.

The layer must contain polygons or exactly two points (opposite corners of
a box). Geometry is transformed to WGS84 first. ``bbox`` is a tuple
``(west, south, east, north)``, so ``fetch_dem(*b.bbox)`` fetches the area.

Parameters
----------
path : str
    Vector file readable by GDAL/OGR (GeoPackage, shapefile, GeoJSON, ...).
layer : str, default=""
    Layer name; empty uses the first layer.

Returns
-------
VectorBounds
    Bounds in decimal degrees plus diagnostics (``ok``, ``message``).

Authors
-------
B. Gailleton
S. M. Mudd

```python
def bounds_from_vector(path: str, layer: str = '') -> lsdtt3.VectorBounds
```

---

## `cop30_tiles_for_bbox`

Return the Copernicus 30 m DEM tile names covering a WGS84 box.

Parameters
----------
west, south, east, north : float
    WGS84 bounds in decimal degrees.

Returns
-------
list of str
    Tile identifiers of the 1 x 1 degree tiles intersecting the box.

Authors
-------
B. Gailleton

```python
def cop30_tiles_for_bbox(west: float, south: float, east: float, north: float) -> list[str]
```

---

## `fetch_dem`

Fetch a DEM for a WGS84 bounding box as an in-memory projected raster.

The source tiles are downloaded, mosaicked and warped to a projected grid
(UTM by default). The Copernicus DEM (the default ``cop30_aws`` source)
stores the sea as exactly 0 m, not as nodata, so mask it before flow
routing: ``dem = mask_to_nodata(dem, 0.0)``. Cells outside the source
footprint are nodata.

Parameters
----------
west, south, east, north : float
    WGS84 bounds in decimal degrees (longitude west/east, latitude
    south/north). Antimeridian-crossing requests are not supported:
    ``west`` must be less than ``east``. ``VectorBounds.bbox`` from
    ``bounds_from_vector`` can be unpacked here: ``fetch_dem(*b.bbox)``.
source : {"cop30_aws", "opentopography", "planetary_computer"}, optional
    ``cop30_aws`` is a keyless global 30 m Copernicus DEM default.
    OpenTopography reads its key only from the ``OT_API_KEY`` or
    ``OT_API_KEY_FILE`` environment variables (never pass a key as an
    argument). ``planetary_computer`` is keyless.
target_epsg : int, optional
    Output EPSG code. Zero selects the UTM zone of the bounding-box centre.
grid_spacing : float, optional
    Output cell size in target CRS units (m for UTM); default 30.
resampling : str, optional
    GDAL warp resampling method (``"cubic"`` default, ``"bilinear"``,
    ``"near"``, ...).
nodata_value : float, optional
    Nodata value of the output (default -9999).
size_pixels : int, optional
    If > 0, crop an exact ``size_pixels`` x ``size_pixels`` square from the
    centre of the warped grid; 0 keeps the full extent.
cache_dir : str, optional
    Persistent cache directory. This is the only optional disk output; the
    returned raster itself is kept in memory. Empty means no caching.
ot_dem_type : str, optional
    OpenTopography dataset, for example ``"COP30"``, ``"SRTMGL1"``,
    ``"NASADEM"`` (see ``opentopography_datasets()``); other sources ignore it.
pc_collection, pc_asset : str, optional
    Planetary Computer STAC collection and asset key (defaults
    ``"cop-dem-glo-30"`` and ``"data"``); other sources ignore them.
http_timeout_seconds : int, optional
    Total time limit per HTTP request in seconds.
http_connect_timeout_seconds : int, optional
    Connection time limit per HTTP request in seconds.
http_max_retry : int, optional
    Number of retries of a failed HTTP request.

Returns
-------
Raster
    A projected raster of elevation (m). Network, GDAL, and invalid-request
    failures raise ``RuntimeError``.

Authors
-------
B. Gailleton
S. M. Mudd

**Deprecated compatibility names:**
- `size_px` -> `size_pixels`; earliest removal: following minor.

```python
def fetch_dem(west: float, south: float, east: float, north: float, source: str = 'cop30_aws', target_epsg: int = 0, grid_spacing: float = 30.0, resampling: str = 'cubic', nodata_value: float = -9999.0, size_pixels: int = 0, cache_dir: str = '', ot_dem_type: str = 'COP30', pc_collection: str = 'cop-dem-glo-30', pc_asset: str = 'data', http_timeout_seconds: int = 60, http_connect_timeout_seconds: int = 15, http_max_retry: int = 3) -> lsdtt3.Raster
```

---

## `fetch_dem_for_outlets`

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
source : {"cop30_aws", "opentopography", "planetary_computer"}, optional
    DEM source, as in ``fetch_dem``.
target_epsg : int, optional
    Output EPSG code; 0 selects UTM from each group's centre.
grid_spacing : float, optional
    Output cell size in target CRS units (m for UTM).
initial_radius_km : float, optional
    Initial half-width (km) of the box around each outlet before grouping.
auto_split_disconnected_basins : bool, optional
    Split groups whose outlets form separate clusters.
split_cluster_gap_km : float, optional
    Gap (km) between outlet clusters above which a group is split.
max_final_dem_area_km2 : float, optional
    Largest allowed grouped initial extent (km^2) before splitting it.
split_overlap_km : float, optional
    Overlap (km) between split extents.
max_basins_per_dem : int, optional
    Largest number of outlets per DEM; 0 means no limit.
cache_dir : str, optional
    Persistent download cache directory; empty means none.
resampling : str, optional
    GDAL warp resampling method.
nodata_value : float, optional
    Nodata value of the output DEMs.
snapping_method : {"none", "expanding_square", "downstream", "largest_area"}, optional
    Outlet snapping used before each basin delineation; the same rules as
    ``FlowInfo.snap_points`` and the CLI ``basin_outlet_snapping_method``.
contributing_pixels_threshold : int, optional
    Channel threshold (contributing pixels) that outlets are snapped to.
edge_buffer_pixels, nodata_buffer_pixels : int, optional
    Required basin clearance (pixels) from the DEM edge and invalid cells.
tight_bbox_buffer_pixels : int, optional
    Margin (pixels) added around the verified basins for the tight refetch.
verify_expand_factor : float, optional
    Factor (> 1) by which an extent grows when a basin reaches its edge.
verify_max_iters : int, optional
    Largest number of fetch-and-verify attempts per group.
tight_bbox_refetch : bool, optional
    Refetch a smaller final DEM around the verified basin footprint.
min_search_radius_pixels : int, optional
    Smallest expanding-square snap radius (pixels).
max_search_radius_pixels : int, optional
    Largest expanding-square snap radius when ``snapping_method`` is
    ``"expanding_square"``, and the half-width of the square searched for
    ``"largest_area"``.
mask_outside_basins : bool, optional
    Set cells outside the selected basin union to nodata in returned DEMs.

Returns
-------
OutletFetchResult
    Grouped raster collection and copied outlet/group diagnostic tables.

Authors
-------
B. Gailleton
S. M. Mudd

```python
def fetch_dem_for_outlets(outlets: object, source: str = 'cop30_aws', target_epsg: int = 0, grid_spacing: float = 30.0, initial_radius_km: float = 100.0, auto_split_disconnected_basins: bool = True, split_cluster_gap_km: float = 25.0, max_final_dem_area_km2: float = 2500.0, split_overlap_km: float = 2.0, max_basins_per_dem: int = 0, cache_dir: str = '', resampling: str = 'cubic', nodata_value: float = -9999.0, snapping_method: str = 'downstream', contributing_pixels_threshold: int = 1000, edge_buffer_pixels: int = 3, nodata_buffer_pixels: int = 1, tight_bbox_buffer_pixels: int = 3, verify_expand_factor: float = 1.25, verify_max_iters: int = 3, tight_bbox_refetch: bool = True, min_search_radius_pixels: int = 0, max_search_radius_pixels: int = 25, mask_outside_basins: bool = False) -> lsdtt3.OutletFetchResult
```

---

## `gdal_data_description`

As proj_data_description(), for GDAL's own data directory.

Parameters
----------
None

Returns
-------
str

Authors
-------
S. M. Mudd
B. Gailleton

```python
def gdal_data_description() -> str
```

---

## `gdal_version`

GDAL release name of the library actually loaded, e.g. "3.8.4".

Parameters
----------
None

Returns
-------
str

Authors
-------
S. M. Mudd
B. Gailleton

```python
def gdal_version() -> str
```

---

## `http_available`

Return whether the bundled GDAL has HTTP support.

Parameters
----------
None

Returns
-------
bool
    True if remote downloads (``fetch_dem``) can work.

Authors
-------
B. Gailleton

```python
def http_available() -> bool
```

---

## `proj_data_description`

Where PROJ is getting proj.db, as a human-readable string.

One of ``"bundled: <dir>"`` (lsdtt3 pinned it), ``"PROJ_DATA=<dir>"``,
``"PROJ_LIB=<dir>"``, or ``"PROJ built-in search"``.

Put this in any bug report about projections or coordinates that look subtly
wrong -- a mismatched proj.db usually returns plausible values rather than
failing.

Parameters
----------
None

Returns
-------
str

Authors
-------
S. M. Mudd
B. Gailleton

```python
def proj_data_description() -> str
```

---

## `proj_version`

PROJ version actually in use at runtime, e.g. "9.3.1".

Reported by the loaded library rather than by the headers this was compiled
against -- that distinction is the whole point.

Parameters
----------
None

Returns
-------
str

Authors
-------
S. M. Mudd
B. Gailleton

```python
def proj_version() -> str
```

---

## `utm_epsg_for_point`

Return the UTM EPSG code containing a WGS84 point.

Parameters
----------
longitude, latitude : float
    WGS84 position in decimal degrees (longitude first).

Returns
-------
int
    EPSG code of the WGS84 UTM zone: 326zz in the northern hemisphere,
    327zz in the southern (for example 32633 for UTM 33N).

Authors
-------
B. Gailleton

```python
def utm_epsg_for_point(longitude: float, latitude: float) -> int
```
