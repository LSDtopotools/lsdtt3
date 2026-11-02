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

## `RotationFrame`

An oblique analysis frame aligned with a feature, and its rotate-back target.

A frame is a Hotine Oblique Mercator (variant B) CRS whose central line runs
along a fault, river, coastline or elongated range, with the feature along
+x and the frame centre at (0, 0). A raster in the frame is an ordinary
north-up raster in that CRS, so every lsdtt3 analysis runs on it unchanged,
on a grid that holds the strip rather than its bounding box. Bearings from
analyses are then clockwise from grid up (90 = along the feature) unless
``bearings_relative_to_true_north=True`` is passed.

The frame also carries the lattice that ``rotate_back`` writes products onto:
the source DEM's own grid for a frame built with ``source=``, else the UTM
(or polar stereographic) CRS of the frame centre on the lattice a fetch uses.
Bring back products (basins, channels, IDs), never the DEM.

Get one with ``RotationFrame.from_feature`` (a line, two points or polygon
file), ``RotationFrame.from_valid_data`` (the footprint of a raster's valid
cells), ``RotationFrame.read`` (a ``*_rotation_frame.txt`` written by a CLI
or by ``write``), or ``RotationFrame.from_raster`` (rebuilt from a frame
raster's own CRS when the frame file is lost). Typical use::

    frame = lsdtt3.RotationFrame.from_feature("fault.gpkg", width=10000.0,
                                              grid_spacing=30.0)
    dem = lsdtt3.fetch_dem(frame=frame, grid_spacing=30.0)
    frame.write("fault_rotation_frame.txt")
    # ... analysis on dem ...
    basins_utm = lsdtt3.rotate_back(basins, frame, resampling="near")

Property names spell out the units of the frame-file keys: ``centre_lon`` is
``centre_longitude``, ``azimuth_deg`` is ``azimuth_degrees``, ``width_m`` is
``width_metres``, ``target_crs_epsg`` is ``target_epsg``, and so on.

Not constructible directly; create one with ``RotationFrame.from_feature()`` or ``RotationFrame.from_raster()`` or ``RotationFrame.from_valid_data()`` or ``RotationFrame.read()``.

Attributes
----------
axis_length_metres : float
    Geodesic length of the axis between its end points (without the end buffer), m.
azimuth_degrees : float
    Azimuth of the axis at the centre, clockwise from true north, in (-90, 90].
base_crs_wkt : str
    WKT of the geographic CRS (datum) the frame is built on.
buffer_all_sides : bool
    True when the end buffer is also on both long sides (polygon frames).
cell_size : float
    Frame cell size, m.
centre_latitude : float
    Latitude of the frame centre, degrees in the base geographic CRS.
centre_longitude : float
    Longitude of the frame centre, degrees in the base geographic CRS.
crs_wkt : str
    The frame CRS as WKT.
end_buffer_metres : float
    Length added beyond each axis end (and each long side for a polygon), m.
end_points : tuple
    The axis end points ((lon, lat), (lon, lat)), west one first.
extent : tuple
    Frame extent (x_min, x_max, y_min, y_max) in frame metres.
feature_fname : str
    File name of the feature the frame was built from, or empty.
feature_layer : str
    Layer of the feature file the frame was built from, or empty.
frame_source : str
    How the frame was defined: line, two_points, polygon, data_footprint or crs.
n_cols : int
    Number of columns of a raster filling the frame.
n_rows : int
    Number of rows of a raster filling the frame.
proj_string : str
    The frame CRS as a PROJ string.
source_dem_fname : str
    File name of the source DEM, when the frame file records one.
target_cell_size : float
    Cell size of the rotate-back target grid, in its CRS units.
target_crs_wkt : str
    WKT of the rotate-back target CRS when it has no EPSG code, else empty.
target_epsg : int
    EPSG code of the rotate-back target CRS (0 when it is given as WKT).
target_origin : tuple
    A node (x, y) of the rotate-back lattice: the source DEM's lower-left corner, or (0, 0).
warnings : list[str]
    Notes from building or reading the frame (also issued as warnings).
width_metres : float
    Full width of the strip across the axis, m.

Methods
-------
RotationFrame.from_feature(path: str, width: float = 0.0, end_buffer: float = -1.0, grid_spacing: float = 0.0, layer: str = '', source: object | None = None, target_epsg: int = 0) -> RotationFrame
    Build the frame of a line, two-point or polygon feature file.
RotationFrame.from_raster(raster: Raster) -> RotationFrame
    Rebuild the frame of a raster that is already in an lsdtt3 oblique frame.
RotationFrame.from_valid_data(source: Raster, end_buffer: float = -1.0, grid_spacing: float = 0.0) -> RotationFrame
    Build the frame of a raster's valid cells (its minimum-area rectangle).
grid_convergence_degrees(x: float, y: float) -> float
    Angle of true north, clockwise from grid up, at a point of the frame.
RotationFrame.read(path: str) -> RotationFrame
    Read a rotation frame file (``*_rotation_frame.txt``).
write(path: str) -> None
    Write the frame and its rotate-back target as a rotation frame file.

Authors
-------
S. M. Mudd
B. Gailleton

### Members

#### `axis_length_metres`

Geodesic length of the axis between its end points (without the end buffer), m.

#### `azimuth_degrees`

Azimuth of the axis at the centre, clockwise from true north, in (-90, 90].

#### `base_crs_wkt`

WKT of the geographic CRS (datum) the frame is built on.

#### `buffer_all_sides`

True when the end buffer is also on both long sides (polygon frames).

#### `cell_size`

Frame cell size, m.

#### `centre_latitude`

Latitude of the frame centre, degrees in the base geographic CRS.

#### `centre_longitude`

Longitude of the frame centre, degrees in the base geographic CRS.

#### `crs_wkt`

The frame CRS as WKT.

#### `end_buffer_metres`

Length added beyond each axis end (and each long side for a polygon), m.

#### `end_points`

The axis end points ((lon, lat), (lon, lat)), west one first.

#### `extent`

Frame extent (x_min, x_max, y_min, y_max) in frame metres.

#### `feature_fname`

File name of the feature the frame was built from, or empty.

#### `feature_layer`

Layer of the feature file the frame was built from, or empty.

#### `frame_source`

How the frame was defined: line, two_points, polygon, data_footprint or crs.

#### `from_feature`

from_feature(path: str, width: float = 0.0, end_buffer: float = -1.0, grid_spacing: float = 0.0, layer: str = '', source: object | None = None, target_epsg: int = 0) -> lsdtt3.RotationFrame

Build the frame of a line, two-point or polygon feature file.

A line (or a multi-part line) uses its two furthest end points as the axis;
two points are the axis ends; a polygon uses its minimum-area enclosing
rectangle, longest side as the axis, and takes its width from the polygon.
The same routine as ``fetch_from_linear_feature`` in lsdtt-fetch-raster
and ``rotate_to_linear_feature`` in lsdtt-raster-preprocessing.

Parameters
----------
path : str
    Feature file (any OGR vector: GeoPackage, GeoJSON, shapefile, ...), in
    any CRS.
width : float, default=0.0
    Full width of the strip across the axis, in metres. Required (> 0) for
    a line or two points; ignored, with a warning, for a polygon.
end_buffer : float, default=-1.0
    Length added beyond each axis end, in metres (on all four sides for a
    polygon). Negative means 10 cells.
grid_spacing : float, default=0.0
    Frame cell size in metres. Without ``source``, 0 means 30 m (the
    ``fetch_dem`` default). With ``source``, 0 means the source's own cell
    size when it is projected in metres, else the size reprojection would
    choose, as lsdtt-raster-preprocessing does.
layer : str, default=""
    Layer to read; empty reads the first layer (warned if there are several).
source : Raster or None, default=None
    The DEM this frame will rotate (``rotate_to_frame``). The frame is built
    on its geographic datum, and the rotate-back target is its own grid
    (CRS, cell size and lower-left corner) when it is projected in metres.
    None builds a frame for ``fetch_dem`` on WGS84.
target_epsg : int, default=0
    Without ``source`` only: EPSG code of the rotate-back target. 0 picks the
    UTM (or polar stereographic beyond 75 degrees) CRS of the frame centre.

Returns
-------
RotationFrame
    The frame. Reader and builder notes (several layers, an ill-defined
    polygon axis) are issued as ``UserWarning``; a file that cannot be read
    or used raises ``ValueError``.

```python
def from_feature(path: str, width: float = 0.0, end_buffer: float = -1.0, grid_spacing: float = 0.0, layer: str = '', source: object | None = None, target_epsg: int = 0) -> lsdtt3.RotationFrame
```

#### `from_raster`

from_raster(raster: lsdtt3.Raster) -> lsdtt3.RotationFrame

Rebuild the frame of a raster that is already in an lsdtt3 oblique frame.

The fallback when the frame file is lost: the frame comes from the raster's
own CRS and extent, and the rotate-back target is the UTM (or polar
stereographic) CRS of the frame centre on the lattice anchored at (0, 0). A
``UserWarning`` says that a source DEM origin off that lattice cannot be
recovered without the frame file.

Parameters
----------
raster : Raster
    A raster in an lsdtt3 oblique Mercator frame.

Returns
-------
RotationFrame
    The rebuilt frame; ``ValueError`` for a raster in any other CRS.

```python
def from_raster(raster: lsdtt3.Raster) -> lsdtt3.RotationFrame
```

#### `from_valid_data`

from_valid_data(source: lsdtt3.Raster, end_buffer: float = -1.0, grid_spacing: float = 0.0) -> lsdtt3.RotationFrame

Build the frame of a raster's valid cells (its minimum-area rectangle).

The polygon of the outer corners of the first and last valid cell of every
row, framed as a polygon feature would be; the frame of
``rotate_to_data_footprint`` in lsdtt-raster-preprocessing. A footprint with
no clear long axis (a square DEM) is warned.

Parameters
----------
source : Raster
    The DEM; it must have a CRS. The rotate-back target is its own grid
    when it is projected in metres.
end_buffer : float, default=-1.0
    Length added on all four sides, in metres; negative means 10 cells.
grid_spacing : float, default=0.0
    Frame cell size in metres; 0 chooses it from ``source`` as
    ``from_feature`` does.

Returns
-------
RotationFrame
    The frame; ``ValueError`` for a raster with no CRS or no valid cells.

```python
def from_valid_data(source: lsdtt3.Raster, end_buffer: float = -1.0, grid_spacing: float = 0.0) -> lsdtt3.RotationFrame
```

#### `grid_convergence_degrees`

grid_convergence_degrees(self, x: float, y: float) -> float

Angle of true north, clockwise from grid up, at a point of the frame.

A grid bearing ``g`` at ``(x, y)`` is the true bearing ``g - c``. At the
frame centre ``c`` is ``90 - azimuth_degrees``.

Parameters
----------
x, y : float
    Frame coordinates in metres.

Returns
-------
float
    The convergence in degrees.

```python
def grid_convergence_degrees(self, x: float, y: float) -> float
```

#### `n_cols`

Number of columns of a raster filling the frame.

#### `n_rows`

Number of rows of a raster filling the frame.

#### `proj_string`

The frame CRS as a PROJ string.

#### `read`

read(path: str) -> lsdtt3.RotationFrame

Read a rotation frame file (``*_rotation_frame.txt``).

The files written by lsdtt-fetch-raster, lsdtt-raster-preprocessing and
``RotationFrame.write`` are flat YAML mappings; this reader checks the
stored CRS against the stored centre and azimuth and runs the frame's axis
self-check.

Parameters
----------
path : str
    The frame file.

Returns
-------
RotationFrame
    The frame and its rotate-back target; ``ValueError`` if the file is
    missing, newer than this lsdtt3, or inconsistent.

```python
def read(path: str) -> lsdtt3.RotationFrame
```

#### `source_dem_fname`

File name of the source DEM, when the frame file records one.

#### `target_cell_size`

Cell size of the rotate-back target grid, in its CRS units.

#### `target_crs_wkt`

WKT of the rotate-back target CRS when it has no EPSG code, else empty.

#### `target_epsg`

EPSG code of the rotate-back target CRS (0 when it is given as WKT).

#### `target_origin`

A node (x, y) of the rotate-back lattice: the source DEM's lower-left corner, or (0, 0).

#### `warnings`

Notes from building or reading the frame (also issued as warnings).

#### `width_metres`

Full width of the strip across the axis, m.

#### `write`

write(self, path: str) -> None

Write the frame and its rotate-back target as a rotation frame file.

The format is the one the CLIs write and read (``rotation_frame_fname`` in
lsdtt-raster-preprocessing), so a frame made in Python can rotate products
back from the command line and the other way round.

Parameters
----------
path : str
    Output file, conventionally ``<stem>_rotation_frame.txt``.

Returns
-------
None
    Nothing; ``RuntimeError`` if the file cannot be written.

```python
def write(self, path: str) -> None
```

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

## `dem_sources`

Return the ``fetch_dem`` DEM source catalogue as dictionaries.

The same rows ``lsdtt-fetch-raster -h dem_source`` describes.

Parameters
----------
None

Returns
-------
list of dict
    One dict per source with keys ``dem_source`` (the ``source`` value for
    ``fetch_dem``), ``title``, ``provider``, ``description``, ``coverage``,
    ``native_resolutions_metres``, ``default_resolution_metres``,
    ``surface_models`` (``"dtm"`` and/or ``"dsm"``),
    ``default_surface_model``, ``vertical_datum``, ``collections`` (the
    ``source_collection_choice`` labels; empty for most sources),
    ``attribution`` (to cite when you use the data), ``accessed`` (when the
    access details were last checked) and ``requires_lerc``.

Authors
-------
B. Gailleton

```python
def dem_sources() -> list
```

---

## `fetch_dem`

Fetch a DEM for a WGS84 bounding box as an in-memory projected raster.

The source tiles are downloaded, mosaicked and warped to a projected grid
(UTM by default), with the same code and the same parameter names as
``lsdtt-fetch-raster``. The Copernicus DEM (the default ``cop30_aws``
source) stores the sea as exactly 0 m, not as nodata, so mask it before
flow routing: ``dem = mask_to_nodata(dem, 0.0)``. Cells outside the source
footprint are nodata.

Besides the global 30 m sources there are keyless national and polar
high-resolution DEMs (``us_3dep``, ``arctic_arcticdem``,
``antarctic_rema``, ``no_dtm``, ``nl_ahn``, ``es_mdt``, ``pl_nmt``,
``ca_hrdem``, ``ch_swissalti3d``, ``nz_lidar``, ``gb_sct_lidar``,
``jp_gsi``); ``dem_sources()`` lists each one's resolutions, surface models
(DSM or DTM), vertical datum, coverage and attribution. Warnings from the
fetch (partial coverage, overviews read, a LERC-compressed source this GDAL
cannot read, older Scottish collections covering the area) are issued as
``UserWarning``.

Or, with ``frame`` (a ``RotationFrame``) or ``fetch_from_linear_feature``,
fetch a strip along a fault, river or range warped straight from the source
tiles into an oblique frame (one resampling pass), with the feature along
+x. The return type is always ``Raster``; keep the ``RotationFrame`` to
rotate products back with ``rotate_back``.

Parameters
----------
west, south, east, north : float, optional
    WGS84 bounds in decimal degrees (longitude west/east, latitude
    south/north). Antimeridian-crossing requests are not supported:
    ``west`` must be less than ``east``. ``VectorBounds.bbox`` from
    ``bounds_from_vector`` can be unpacked here: ``fetch_dem(*b.bbox)``.
    Required unless ``frame`` or ``fetch_from_linear_feature`` gives the
    area; refused with them.
source : str, optional
    ``cop30_aws`` (default) is a keyless global 30 m Copernicus DEM.
    OpenTopography reads its key only from the ``OT_API_KEY`` or
    ``OT_API_KEY_FILE`` environment variables (never pass a key as an
    argument). ``planetary_computer`` is keyless. The national and polar
    sources are listed by ``dem_sources()``.
target_epsg : int, optional
    Output EPSG code. Zero selects the UTM zone of the bounding-box centre.
grid_spacing : float or None, optional
    Output cell size in target CRS units (m for UTM). None (default) uses
    the frame's cell size with ``frame``, the source's native resolution
    (``source_resolution_metres`` if set) for a national or polar source,
    and 30 otherwise. Output cells are always square.
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
frame : RotationFrame or None, optional
    Fetch the frame's rectangle, warped into the frame. ``grid_spacing``,
    if given, must equal ``frame.cell_size``; ``target_epsg`` and ``size_pixels``
    are refused (the frame carries its rotate-back target).
fetch_from_linear_feature : bool, optional
    Build the frame from ``fetch_feature_fname`` here, as lsdtt-fetch-raster
    does; the same as ``frame=RotationFrame.from_feature(fetch_feature_fname,
    linear_feature_fetch_width, linear_feature_fetch_end_buffer,
    grid_spacing, fetch_feature_layer)``, which is the frame to keep for
    ``rotate_back``.
fetch_feature_fname : str, optional
    Line, two-point or polygon file (any OGR format and CRS).
fetch_feature_layer : str, optional
    Layer of the feature file; empty reads the first.
linear_feature_fetch_width : float, optional
    Full strip width in metres; required for a line or two points.
linear_feature_fetch_end_buffer : float, optional
    Length beyond each axis end in metres (all sides for a polygon);
    negative means 10 cells.
rotate_after_fetch : bool, optional
    With ``fetch_from_linear_feature``: True (default) returns the strip in
    the oblique frame; False returns the frame's WGS84 envelope in UTM, as a
    plain fetch.
source_resolution_metres : float, optional
    Native resolution (m) to fetch from a source that offers several (for
    example 1, 2 or 10 for ``arctic_arcticdem``); 0 (default) uses the
    source's default. Refused for a source that does not offer it.
surface_model_choice : {"source_default", "dtm", "dsm"}, optional
    Terrain (DTM) or surface (DSM) model. ``source_default`` takes the DTM
    where the source has one. Asking for a model the source lacks raises
    ``RuntimeError``; the Copernicus sources are DSMs only.
source_collection_choice : str, optional
    ``gb_sct_lidar`` only: ``latest`` (default) uses the newest collection
    with tiles in the area (a warning names older ones that also cover it);
    a collection label (``phase-1`` ... ``phase-6``,
    ``national-lidar-programme``) uses that one. Collections are never mixed.
max_output_pixels : int, optional
    Largest output grid (cells) allowed, checked before any download; for
    the WCS sources the native cells requested are counted too.
max_download_megabytes : float, optional
    Largest download allowed (MB): an estimate checked before the fetch,
    and a running count that stops the transfer where lsdtt3 does the HTTP.
max_tiles : int, optional
    Largest number of tiles, source files or service requests allowed.
request_delay_milliseconds : int, optional
    Pause between requests to tile services (ms).

Returns
-------
Raster
    A projected raster of elevation (m): UTM (or ``target_epsg``), or the
    oblique frame. Network, GDAL, and invalid-request failures (including
    an area outside the source's coverage or over a size limit) raise
    ``RuntimeError``; an inconsistent combination of arguments raises
    ``ValueError``.

Authors
-------
B. Gailleton
S. M. Mudd

**Deprecated compatibility names:**
- `size_px` -> `size_pixels`; earliest removal: following minor.

```python
def fetch_dem(west: float = 'nan', south: float = 'nan', east: float = 'nan', north: float = 'nan', source: str = 'cop30_aws', target_epsg: int = 0, grid_spacing: float | None = None, resampling: str = 'cubic', nodata_value: float = -9999.0, size_pixels: int = 0, cache_dir: str = '', ot_dem_type: str = 'COP30', pc_collection: str = 'cop-dem-glo-30', pc_asset: str = 'data', http_timeout_seconds: int = 60, http_connect_timeout_seconds: int = 15, http_max_retry: int = 3, frame: object | None = None, fetch_from_linear_feature: bool = False, fetch_feature_fname: str = '', fetch_feature_layer: str = '', linear_feature_fetch_width: float = 0.0, linear_feature_fetch_end_buffer: float = -1.0, rotate_after_fetch: bool = True, source_resolution_metres: float = 0.0, surface_model_choice: str = 'source_default', source_collection_choice: str = 'latest', max_output_pixels: int = 250000000, max_download_megabytes: float = 2048.0, max_tiles: int = 256, request_delay_milliseconds: int = 200) -> lsdtt3.Raster
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
    DEM source, as in ``fetch_dem``. The national and polar sources are
    refused here, as in ``lsdtt-fetch-raster`` outlets mode: find the groups
    with ``cop30_aws``, then ``fetch_dem`` each group's area.
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
max_output_pixels, max_download_megabytes, max_tiles, request_delay_milliseconds
    Size limits and request pacing for each group's fetch, as in
    ``fetch_dem``.

Returns
-------
OutletFetchResult
    Grouped raster collection and copied outlet/group diagnostic tables.

Authors
-------
B. Gailleton
S. M. Mudd

```python
def fetch_dem_for_outlets(outlets: object, source: str = 'cop30_aws', target_epsg: int = 0, grid_spacing: float = 30.0, initial_radius_km: float = 100.0, auto_split_disconnected_basins: bool = True, split_cluster_gap_km: float = 25.0, max_final_dem_area_km2: float = 2500.0, split_overlap_km: float = 2.0, max_basins_per_dem: int = 0, cache_dir: str = '', resampling: str = 'cubic', nodata_value: float = -9999.0, snapping_method: str = 'downstream', contributing_pixels_threshold: int = 1000, edge_buffer_pixels: int = 3, nodata_buffer_pixels: int = 1, tight_bbox_buffer_pixels: int = 3, verify_expand_factor: float = 1.25, verify_max_iters: int = 3, tight_bbox_refetch: bool = True, min_search_radius_pixels: int = 0, max_search_radius_pixels: int = 25, mask_outside_basins: bool = False, max_output_pixels: int = 250000000, max_download_megabytes: float = 2048.0, max_tiles: int = 256, request_delay_milliseconds: int = 200) -> lsdtt3.OutletFetchResult
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

## `rotate_back`

Warp a product in an oblique frame back onto the frame's rotate-back grid.

The routine of ``run_rotate_back`` in lsdtt-raster-preprocessing. A raster at
the frame's cell size goes to the target cell size, any other keeps its own;
the output lattice is anchored at the target origin, so after
``rotate_to_frame`` with a frame built from a projected DEM the product lands
on that DEM's own pixels. Rotate back products (basins, channel masks, IDs,
metrics), never the DEM itself: a round trip resamples twice.

Parameters
----------
raster : Raster
    A product in ``frame.crs_wkt``; any other CRS raises ``ValueError``.
frame : RotationFrame
    The frame the product was made in.
resampling : str, default="bilinear"
    GDAL warp method (as ``rotate_to_frame``); use ``"near"`` for
    categorical rasters (basin IDs, masks, junctions, lithology).
target_epsg : int, default=0
    EPSG code overriding the frame's target CRS (``rotate_back_epsg``); a
    CRS other than the frame's target is gridded on the lattice anchored at
    (0, 0), with a warning. 0 uses the frame's target.

Returns
-------
Raster
    The product on the target grid.

Authors
-------
S. M. Mudd
B. Gailleton

```python
def rotate_back(raster: lsdtt3.Raster, frame: lsdtt3.RotationFrame, resampling: str = 'bilinear', target_epsg: int = 0) -> lsdtt3.Raster
```

---

## `rotate_to_frame`

Warp a raster into an oblique frame, in one resampling pass from its own CRS.

The output is a north-up raster in ``frame.crs_wkt`` covering the frame
(``frame.n_cols`` x ``frame.n_rows`` cells of ``frame.cell_size``), with the
feature along +x; cells outside the source are nodata. The routine of
``rotate_to_linear_feature`` / ``rotate_to_data_footprint`` in
lsdtt-raster-preprocessing, so the same frame and resampling give the same
raster.

Parameters
----------
raster : Raster
    The DEM (or any raster) to rotate; it must have a CRS.
frame : RotationFrame
    The frame, normally built with ``source=raster`` so that rotate-back
    lands on the raster's own pixels.
resampling : str, default="cubic"
    GDAL warp method: near, bilinear, cubic, cubicspline, lanczos, average,
    mode, max, min, med, q1, q3, sum (``resampling_method`` in
    lsdtt-raster-preprocessing).

Returns
-------
Raster
    The raster in the frame.

Authors
-------
S. M. Mudd
B. Gailleton

```python
def rotate_to_frame(raster: lsdtt3.Raster, frame: lsdtt3.RotationFrame, resampling: str = 'cubic') -> lsdtt3.Raster
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
