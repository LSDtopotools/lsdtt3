# Flow, Channels, And Basins API

This page is generated from the built Python runtime manifest. Do not edit signatures here.

## `BasinCollection`

Selected D8 drainage basins (possibly nested) on a ``FlowInfo`` topology.

Made by ``FlowInfo.basins(outlets)``. Basin ``basin_id`` (0-based) is the
whole upslope area of its outlet pixel; basins may nest (a basin whose
outlet lies inside another basin). Each basin is stored as a view of the
flow topology, so collections are cheap. ``table`` describes each basin,
``polygons`` gives outlines, ``labels`` and ``mask()`` give rasters, and
``statistics(dem)`` summarises elevation, slope and shape.

Returned by ``line_crossing_basins()``, ``FlowInfo.basins()``, ``SideTributaryBasins.basins``; not constructible directly.

Attributes
----------
labels : Raster
    Raster of basin membership: ``basin_id + 1``, 0 outside all basins.
polygons : geopandas.GeoDataFrame
    Basin outline polygons as a GeoDataFrame.
size : int
    Number of basins.
table : pandas.DataFrame
    One row per basin (pandas.DataFrame).

Methods
-------
mask(basin_id: int) -> Raster
    Return a 0/1 Raster of one basin.
membership(row: int, col: int) -> object
    Return all selected basin IDs containing this pixel.
outline_ridgelines(dem: Raster, split_method: str = 'longest_channel') -> object
    Split every selected basin outline into two ridgelines.
pixels(basin_id: int) -> object
    Return the pixels of one basin.
statistics(dem: Raster, include_slope: bool = True, include_flow_distance: bool = True, include_shape_metrics: bool = False, perimeter_simplify_cells: float = 1.0, bearings_relative_to_true_north: bool = False) -> object
    Return per-basin elevation, slope, flow-distance, and optional shape statistics.

Authors
-------
B. Gailleton
S. M. Mudd

### Members

#### `labels`

Raster of basin membership: ``basin_id + 1``, 0 outside all basins.

Basins are painted in ``basin_id`` order, so a pixel inside several nested
basins carries the largest ``basin_id + 1`` among them. Nodata where the DEM
is nodata.

#### `mask`

mask(self, basin_id: int) -> lsdtt3.Raster

Return a 0/1 Raster of one basin.

Parameters
----------
basin_id : int
    Zero-based basin index (a ``table.basin_id`` value).

Returns
-------
Raster
    1 inside the basin, 0 on other valid pixels.

```python
def mask(self, basin_id: int) -> lsdtt3.Raster
```

#### `membership`

membership(self, row: int, col: int) -> object

Return all selected basin IDs containing this pixel.

Parameters
----------
row, col : int
    Zero-based pixel position of a valid pixel.

Returns
-------
numpy.ndarray
    Int32 array of zero-based basin IDs (empty if the pixel is in none).

```python
def membership(self, row: int, col: int) -> object
```

#### `outline_ridgelines`

outline_ridgelines(self, dem: lsdtt3.Raster, split_method: str = 'longest_channel') -> object

Split every selected basin outline into two ridgelines.

Parameters
----------
dem : Raster
    Elevation (m) on the flow grid.
split_method : {"longest_channel", "maximum_elevation", "equal_length"}, default="longest_channel"
    How the outline is split in two: at the
    head of the longest flow path, at the highest outline point, or halfway
    round the perimeter.

Returns
-------
geopandas.GeoDataFrame
    Left/right ridge lines with columns ``basin_id``, ``ridgeline_id``,
    ``outlet_flat``, ``n_nodes``, ``length_metres``, ``elevation_min``,
    ``elevation_max`` (m), ``geometry``.

```python
def outline_ridgelines(self, dem: lsdtt3.Raster, split_method: str = 'longest_channel') -> object
```

#### `pixels`

pixels(self, basin_id: int) -> object

Return the pixels of one basin.

Parameters
----------
basin_id : int
    Zero-based basin index.

Returns
-------
numpy.ndarray
    N x 2 int64 array of ``(row, col)``; the first row is the outlet.

```python
def pixels(self, basin_id: int) -> object
```

#### `polygons`

Basin outline polygons as a GeoDataFrame.

Columns ``basin_id``, ``outlet_flat``, ``geometry`` (polygon traced along
pixel edges), in the DEM CRS.

#### `size`

Number of basins.

#### `statistics`

statistics(self, dem: lsdtt3.Raster, include_slope: bool = True, include_flow_distance: bool = True, include_shape_metrics: bool = False, perimeter_simplify_cells: float = 1.0, bearings_relative_to_true_north: bool = False) -> object

Return per-basin elevation, slope, flow-distance, and optional shape statistics.

Parameters
----------
dem : Raster
    Elevation (m) on the flow grid.
include_slope : bool, default=True
    Add ``slope_25th``, ``slope_50th``, ``slope_75th``: percentiles of the
    3 x 3 gradient (m/m, as ``slope()``).
include_flow_distance : bool, default=True
    Add ``flow_distance_25th/50th/75th``: percentiles of D8 flow distance
    to the outlet (m).
include_shape_metrics : bool, default=False
    Add planform shape metrics: ``centroid_x``, ``centroid_y``,
    ``basin_length_metres``, ``elongation_ratio``, ``form_factor``,
    ``circularity_ratio``, ``gravelius_compactness``, ``solidity``,
    ``eccentricity``, ``orientation_degrees`` (principal axis, compass,
    0-180, from grid north unless ``bearings_relative_to_true_north``),
    ``hypsometric_integral``.
perimeter_simplify_cells : float, default=1.0
    Simplification tolerance of the outline, in cells, used for the
    perimeter in the shape metrics.
bearings_relative_to_true_north : bool, default=False
    Convert ``orientation_degrees`` to true north with the grid convergence
    at the basin centroid, as the CLI's ``bearings_relative_to_true_north``;
    needs a flow grid with a CRS. False keeps it relative to grid north.

Returns
-------
pandas.DataFrame
    One row per basin: ``basin_id``, ``n_pixels``, ``area_metres2``,
    ``elevation_min``, ``elevation_max``, ``elevation_mean``,
    ``elevation_median``, ``elevation_stddev`` (m), plus the optional
    columns above.

```python
def statistics(self, dem: lsdtt3.Raster, include_slope: bool = True, include_flow_distance: bool = True, include_shape_metrics: bool = False, perimeter_simplify_cells: float = 1.0, bearings_relative_to_true_north: bool = False) -> object
```

#### `table`

One row per basin (pandas.DataFrame).

Columns: ``basin_id`` (0-based), ``outlet_flat``, ``outlet_row``,
``outlet_col``, ``n_pixels``, ``area_metres2`` (m^2), ``complete`` (False if
any basin pixel, including the outlet, is next to a nodata pixel or to an
open or base-level DEM edge, so the basin may be cut off), ``outlet_on_dem_edge``,
``parent_id`` (the smallest other selected basin containing this outlet,
-1 if none), ``nesting_depth`` (0 for a basin with no parent).

---

## `ChannelCandidatesResult`

Channel-bottom candidate score, mask and diagnostics from :func:`find_channel_candidates`.

``score`` is a weighted 0-1 channelness score and ``mask`` its 0/1 threshold;
the other rasters are the intermediate evidence (per-scale curvature, scale
persistence, local depth, slope, Hessian elongation).

Returned by ``find_channel_candidates()``; not constructible directly.

Attributes
----------
best_scale_metres : Raster
curvature_raw_by_scale : list[Raster]
curvature_signal_by_scale : list[Raster]
depth_signal : Raster
hessian_elongation : Raster
mask : Raster
persistence_by_adjacent_scale_pair : list[Raster]
persistent : Raster
scales_metres : list[float]
score : Raster
slope_signal : Raster

Authors
-------
B. Gailleton

### Members

#### `best_scale_metres`

(self) -> lsdtt3.Raster

#### `curvature_raw_by_scale`

(self) -> list[lsdtt3.Raster]

#### `curvature_signal_by_scale`

(self) -> list[lsdtt3.Raster]

#### `depth_signal`

(self) -> lsdtt3.Raster

#### `hessian_elongation`

(self) -> lsdtt3.Raster

#### `mask`

(self) -> lsdtt3.Raster

#### `persistence_by_adjacent_scale_pair`

(self) -> list[lsdtt3.Raster]

#### `persistent`

(self) -> lsdtt3.Raster

#### `scales_metres`

(self) -> list[float]

#### `score`

(self) -> lsdtt3.Raster

#### `slope_signal`

(self) -> lsdtt3.Raster

---

## `ChannelHeads`

Channel-head locations resolved from a D8 flow topology.

Made by :func:`extract_channel_heads`. Use ``pixels`` as explicit
``FlowInfo.channel_scope(..., sources=heads.pixels)`` input, or pass the
object to ``FlowInfo.scope_from_heads``. The diagnostic rasters are None
unless ``method="wiener"`` and ``diagnostics=True``.

Returned by ``extract_channel_heads()``; not constructible directly.

Attributes
----------
channel_mask : Raster | None
    0/1 candidate channel mask after pruning, or None.
filtered_dem : Raster | None
    Wiener-filtered DEM (m), or None.
method : str
    Method used: 'threshold' or 'wiener'.
n_heads : int
    Number of channel heads.
pixels : numpy.ndarray
    Channel heads as an N x 2 int64 array of (row, col).
qq_table : pandas.DataFrame
    Quantile-quantile table of positive curvature: standard_normal_variate, data_quantile, normal_reference (empty unless diagnostics).
skeleton : Raster | None
    0/1 thinned skeleton of the channel mask, or None.
table : pandas.DataFrame
    One row per head: head_id, flat, row, col, x, y (pixel centre, CRS units), drainage_area (m^2), flow_distance_to_outlet (m).
tangential_curvature : Raster | None
    Tangential curvature (1/m) of the filtered DEM, or None.

Authors
-------
B. Gailleton

### Members

#### `channel_mask`

0/1 candidate channel mask after pruning, or None.

#### `filtered_dem`

Wiener-filtered DEM (m), or None.

#### `method`

Method used: 'threshold' or 'wiener'.

#### `n_heads`

Number of channel heads.

#### `pixels`

Channel heads as an N x 2 int64 array of (row, col).

#### `qq_table`

Quantile-quantile table of positive curvature: standard_normal_variate, data_quantile, normal_reference (empty unless diagnostics).

#### `skeleton`

0/1 thinned skeleton of the channel mask, or None.

#### `table`

One row per head: head_id, flat, row, col, x, y (pixel centre, CRS units), drainage_area (m^2), flow_distance_to_outlet (m).

#### `tangential_curvature`

Tangential curvature (1/m) of the filtered DEM, or None.

---

## `ChannelNetwork`

A D8 channel network as a graph of junctions and links.

Made by ``FlowInfo.channel_network(scope)``: every resolved source of the
``ChannelScope`` is traced downstream to its D8 outlet (outlet truncation in
the scope does not cut the network). A junction is a source, a confluence or
an outlet pixel; a link is the channel between a junction and the next
junction downstream. Stream orders are Strahler orders. Use
``lines_geodataframe()`` / ``points_geodataframe()`` for GIS layers,
``channel_mask`` as the ``channels`` raster of HAND, and ``chi()`` for chi
analysis.

Returned by ``FlowInfo.channel_network()``, ``DivideMigrationResult.network``; not constructible directly.

Attributes
----------
channel_mask : Raster
    0/1 Raster: 1 on network channel pixels, 0 on other valid pixels.
junction_id : numpy.ndarray
    Int32 array (n_rows, n_cols) of the junction ID at junction pixels, -1 elsewhere.
n_channel_cells : int
    Number of channel pixels in the network.
n_junctions : int
    Number of junctions (sources, confluences and outlets).

Methods
-------
basin_com_bearings(bearings_relative_to_true_north: bool = False) -> object
    Return outlet-to-basin-centre-of-mass bearings for each junction basin.
chi(outlets: object | None = None, m_over_n: float = 0.45, reference_area: float = 1.0, denominator: object | None = None) -> ChiCollection
    Extract outlet-seeded chi coordinates on this channel network.
chi_sweep(dem: Raster, m_over_n_values: object, outlets: object | None = None, reference_area: float = 1.0, critical_divergence_metres: float = 5.0, denominator: object | None = None) -> object
    Fit each tested m/n and return per-basin disorder.
donors(junction_id: int) -> list[int]
    Return the IDs of the junctions that drain directly into a junction.
junction_angles(junctions: object | None = None, mode: str = 'full_link', fixed_length_metres: float = 0.0, fixed_n_pixels: int = 0, min_nodes: int = 2) -> object
    Measure two-donor confluence angles in compass degrees.
junctions() -> object
    Return junction topology as a pandas DataFrame.
junctions_geodataframe() -> object
    Return channel junctions as a GeoPandas point layer.
lines_geodataframe() -> object
    Return compact non-outlet channel links as a GeoPandas line layer.
link_bearings(junctions: object | None = None, mode: str = 'full_link', fixed_length_metres: float = 0.0, fixed_n_pixels: int = 0, anchor: str = 'upstream', min_nodes: int = 2, bearings_relative_to_true_north: bool = False) -> object
    Fit compass bearings to compact channel links.
link_profiles(values: Raster, junctions: object | None = None) -> object
    Sample a raster along compact downstream channel links.
links() -> object
    Return junction-to-junction links as a pandas DataFrame.
points_geodataframe() -> object
    Return channel pixels as a GeoPandas point layer.

Authors
-------
B. Gailleton
S. M. Mudd

### Members

#### `basin_com_bearings`

basin_com_bearings(self, bearings_relative_to_true_north: bool = False) -> object

Return outlet-to-basin-centre-of-mass bearings for each junction basin.

Parameters
----------
bearings_relative_to_true_north : bool, default=False
    False: clockwise from grid north (grid up). True: converted to true north
    with the grid convergence at the midpoint of the outlet and the centre of mass, as the CLI's
    ``bearings_relative_to_true_north``; needs a flow grid with a CRS.

Returns
-------
pandas.DataFrame
    Columns ``junction_id``, ``stream_order``, ``outlet_x``, ``outlet_y``,
    ``centre_of_mass_x``, ``centre_of_mass_y``, ``bearing_deg`` (compass
    degrees from the outlet to the centre of mass, from grid north or true
    north as above), ``distance_metres``,
    ``n_pixels``.

```python
def basin_com_bearings(self, bearings_relative_to_true_north: bool = False) -> object
```

#### `channel_mask`

0/1 Raster: 1 on network channel pixels, 0 on other valid pixels.

#### `chi`

chi(self, outlets: object | None = None, m_over_n: float = 0.45, reference_area: float = 1.0, denominator: object | None = None) -> lsdtt3.ChiCollection

Extract outlet-seeded chi coordinates on this channel network.

chi is the integral of ``(A0 / A)**m_over_n`` along the channel from the
outlet upstream (trapezoidal rule on D8 steps), so it has units of metres
and is 0 at each outlet. Plotting elevation against chi gives straight
lines where channel steepness is uniform.

Parameters
----------
outlets : array-like of int or None, default=None
    N x 2 array of ``(row, col)`` outlet pixels on the network, one basin
    each. None uses every outlet of the network.
m_over_n : float, default=0.45
    Concavity index theta = m/n (dimensionless, positive).
reference_area : float, default=1.0
    Reference drainage area A0 in m^2 (positive). With A0 = 1 m^2 the
    fitted chi-elevation slope equals the steepness index k_sn.
denominator : Raster or None, default=None
    Replacement for drainage area A in the integrand (for example a
    discharge raster, in the same units as A0); must be positive on the
    channels. None uses D8 drainage area (m^2).

Returns
-------
ChiCollection
    Chi values for each outlet's channel pixels.

```python
def chi(self, outlets: object | None = None, m_over_n: float = 0.45, reference_area: float = 1.0, denominator: object | None = None) -> lsdtt3.ChiCollection
```

#### `chi_sweep`

chi_sweep(self, dem: lsdtt3.Raster, m_over_n_values: object, outlets: object | None = None, reference_area: float = 1.0, critical_divergence_metres: float = 5.0, denominator: object | None = None) -> object

Fit each tested m/n and return per-basin disorder.

For each value the chi profiles are fitted (``ChiCollection.fit`` with
``mainstem="max_flow_length"``) and the chi disorder is computed; the m/n
with the lowest disorder is the best-collapsing concavity.

Parameters
----------
dem : Raster
    Elevation (m) on the network grid.
m_over_n_values : array-like of float
    Concavity values to test (non-empty).
outlets, reference_area, denominator
    As in ``chi``.
critical_divergence_metres : float, default=5.0
    As in ``ChiCollection.fit``.

Returns
-------
pandas.DataFrame
    Columns ``m_over_n``, ``basin_id``, ``outlet_flat``, ``disorder``.

```python
def chi_sweep(self, dem: lsdtt3.Raster, m_over_n_values: object, outlets: object | None = None, reference_area: float = 1.0, critical_divergence_metres: float = 5.0, denominator: object | None = None) -> object
```

#### `donors`

donors(self, junction_id: int) -> list[int]

Return the IDs of the junctions that drain directly into a junction.

Parameters
----------
junction_id : int
    Junction ID (see ``junctions()``).

Returns
-------
list of int
    Donor junction IDs; empty for a source.

```python
def donors(self, junction_id: int) -> list[int]
```

#### `junction_angles`

junction_angles(self, junctions: object | None = None, mode: str = 'full_link', fixed_length_metres: float = 0.0, fixed_n_pixels: int = 0, min_nodes: int = 2) -> object

Measure two-donor confluence angles in compass degrees.

``donor_opening_angle_deg`` is the angle between the two donor channels
(the "junction angle" of Strong and Mudd 2023); ``donor_1/2_receiver_angle_deg``
are the angles between each donor and the receiver.

Parameters
----------
junctions : array-like of int or None, default=None
    Junction IDs to measure; None measures every two-donor junction.
mode, fixed_length_metres, fixed_n_pixels, min_nodes
    How much of each link is fitted, as in ``link_bearings``.

Returns
-------
pandas.DataFrame
    Columns ``junction_id``, ``donor_1``, ``donor_2``, ``receiver``,
    ``stream_order``, ``donor_1_area_m2``, ``donor_2_area_m2``,
    ``junction_area_m2``, ``donor_opening_angle_deg``,
    ``donor_1_receiver_angle_deg``, ``donor_2_receiver_angle_deg``,
    ``bending_angle_signed_deg``, ``bending_angle_abs_deg``.

```python
def junction_angles(self, junctions: object | None = None, mode: str = 'full_link', fixed_length_metres: float = 0.0, fixed_n_pixels: int = 0, min_nodes: int = 2) -> object
```

#### `junction_id`

Int32 array (n_rows, n_cols) of the junction ID at junction pixels, -1 elsewhere.

#### `junctions`

junctions(self) -> object

Return junction topology as a pandas DataFrame.

Returns
-------
pandas.DataFrame
    One row per junction: ``junction_id``, ``flat``, ``row``, ``col``,
    ``receiver_junction`` (next junction downstream, -1 at an outlet),
    ``donor_count`` (junctions draining directly into it), ``stream_order``
    (Strahler), ``is_outlet``.

```python
def junctions(self) -> object
```

#### `junctions_geodataframe`

junctions_geodataframe(self) -> object

Return channel junctions as a GeoPandas point layer.

Returns
-------
geopandas.GeoDataFrame
    Columns as ``junctions()`` plus ``geometry`` (pixel centres, DEM CRS).

```python
def junctions_geodataframe(self) -> object
```

#### `lines_geodataframe`

lines_geodataframe(self) -> object

Return compact non-outlet channel links as a GeoPandas line layer.

Returns
-------
geopandas.GeoDataFrame
    One LineString per link (pixel centres from a junction down to the next
    junction), in the DEM CRS. Columns ``upstream_junction``,
    ``downstream_junction``, ``stream_order`` (Strahler), ``length_metres``,
    ``drainage_area`` (m^2, at the upstream junction), ``geometry``.

```python
def lines_geodataframe(self) -> object
```

#### `link_bearings`

link_bearings(self, junctions: object | None = None, mode: str = 'full_link', fixed_length_metres: float = 0.0, fixed_n_pixels: int = 0, anchor: str = 'upstream', min_nodes: int = 2, bearings_relative_to_true_north: bool = False) -> object

Fit compass bearings to compact channel links.

Parameters
----------
junctions : array-like of int or None, default=None
    Upstream junction IDs of the links to fit; None fits every link.
mode : {"full_link", "fixed_length", "fixed_n_pixels"}, default="full_link"
    Fit the whole link, or only a sub-chain of the link.
fixed_length_metres : float, default=0
    Sub-chain path length (m) for ``mode="fixed_length"``; must be positive
    in that mode.
fixed_n_pixels : int, default=0
    Sub-chain length in pixels for ``mode="fixed_n_pixels"`` (at least 2).
anchor : {"upstream", "downstream"}, default="upstream"
    End of the link the sub-chain starts from.
min_nodes : int, default=2
    Minimum number of pixels a fit needs (at least 2).
bearings_relative_to_true_north : bool, default=False
    False: clockwise from grid north (grid up). True: converted to true north
    with the grid convergence at each link's representative point (``representative_x``,
    ``representative_y``); ``direction_x`` and ``direction_y`` stay grid
    components, as the CLI's
    ``bearings_relative_to_true_north``; needs a flow grid with a CRS.

Returns
-------
pandas.DataFrame
    Columns ``upstream_junction``, ``downstream_junction``,
    ``stream_order``, ``upstream_flat``, ``downstream_flat``,
    ``drainage_area_m2``, ``n_nodes_total``, ``n_nodes_used``,
    ``total_length_metres``, ``sampled_length_metres``, ``direction_x``,
    ``direction_y``, ``bearing_deg`` (compass degrees, 0 = grid north,
    or true north as above, clockwise), ``fit_linearity``, ``representative_flat``,
    ``representative_x``, ``representative_y``.

```python
def link_bearings(self, junctions: object | None = None, mode: str = 'full_link', fixed_length_metres: float = 0.0, fixed_n_pixels: int = 0, anchor: str = 'upstream', min_nodes: int = 2, bearings_relative_to_true_north: bool = False) -> object
```

#### `link_profiles`

link_profiles(self, values: lsdtt3.Raster, junctions: object | None = None) -> object

Sample a raster along compact downstream channel links.

Parameters
----------
values : Raster
    Raster on the network grid to sample (for example the DEM).
junctions : array-like of int or None, default=None
    Upstream junction IDs of the links to sample; None samples every link.

Returns
-------
pandas.DataFrame
    One row per link pixel: ``upstream_junction``, ``downstream_junction``,
    ``node_index``, ``flat``, ``row``, ``col``,
    ``distance_from_upstream_metres``, ``value`` (NaN where nodata).

```python
def link_profiles(self, values: lsdtt3.Raster, junctions: object | None = None) -> object
```

#### `links`

links(self) -> object

Return junction-to-junction links as a pandas DataFrame.

Returns
-------
pandas.DataFrame
    One row per link: ``upstream_junction``, ``downstream_junction``,
    ``stream_order`` (Strahler order of the upstream junction),
    ``length_metres`` (D8 path length).

```python
def links(self) -> object
```

#### `n_channel_cells`

Number of channel pixels in the network.

#### `n_junctions`

Number of junctions (sources, confluences and outlets).

#### `points_geodataframe`

points_geodataframe(self) -> object

Return channel pixels as a GeoPandas point layer.

Returns
-------
geopandas.GeoDataFrame
    One point per channel pixel at the pixel centre, in the DEM CRS. Columns
    ``flat``, ``row``, ``col``, ``junction_id`` (junction of the link the
    pixel belongs to), ``is_junction``, ``drainage_area`` (m^2),
    ``flow_distance_to_outlet`` (m), ``geometry``.

```python
def points_geodataframe(self) -> object
```

---

## `ChannelScope`

The set of channel pixels selected for analysis, and the source of each.

Made by ``FlowInfo.channel_scope`` (or ``FlowInfo.scope_from_heads``). The
mode depends on what was supplied: ``"all_channels_threshold"`` (no
sources: every pixel with at least the threshold contributing pixels),
``"user_sources"`` (each source traced downstream to its outlet) or
``"user_sources_bounded"`` (each traced channel stops at its matched
outlet). Pass it to ``FlowInfo.channel_network`` to build a junction graph.

Returned by ``FlowInfo.channel_scope()``, ``FlowInfo.scope_from_heads()``; not constructible directly.

Attributes
----------
mask : Raster
    0/1 Raster of in-scope channel pixels.
mode : str
    Scope mode: 'all_channels_threshold', 'user_sources' or 'user_sources_bounded'.
outlets : numpy.ndarray
    Per resolved source, the (row, col) where its channel is cut, or (-1, -1) if it runs to the DEM outlet.
source_id : numpy.ndarray
    Int32 array (n_rows, n_cols): index of the source that first claims each channel pixel, -1 elsewhere.
sources : numpy.ndarray
    Resolved (snapped) source pixels as an N x 2 int64 (row, col) array; unsnappable inputs are dropped.

Methods
-------
nodes(source_id: int) -> object
    Return the ordered channel pixels of one source, upstream to downstream.
table() -> object
    Return ordered channel nodes as a pandas DataFrame.

Authors
-------
B. Gailleton

### Members

#### `mask`

0/1 Raster of in-scope channel pixels.

#### `mode`

Scope mode: 'all_channels_threshold', 'user_sources' or 'user_sources_bounded'.

#### `nodes`

nodes(self, source_id: int) -> object

Return the ordered channel pixels of one source, upstream to downstream.

Parameters
----------
source_id : int
    Zero-based index into ``sources``.

Returns
-------
numpy.ndarray
    N x 2 int64 array of ``(row, col)``. In ``"all_channels_threshold"``
    mode a list stops where it joins a channel claimed by an earlier source.

```python
def nodes(self, source_id: int) -> object
```

#### `outlets`

Per resolved source, the (row, col) where its channel is cut, or (-1, -1) if it runs to the DEM outlet.

#### `source_id`

Int32 array (n_rows, n_cols): index of the source that first claims each channel pixel, -1 elsewhere.

#### `sources`

Resolved (snapped) source pixels as an N x 2 int64 (row, col) array; unsnappable inputs are dropped.

#### `table`

table(self) -> object

Return ordered channel nodes as a pandas DataFrame.

Returns
-------
pandas.DataFrame
    Columns ``source_id`` (index into ``sources``), ``node_index`` (0 at the
    source, increasing downstream), ``flat``, ``row``, ``col``.

```python
def table(self) -> object
```

---

## `ConnectivityResult`

Sediment-connectivity rasters following Cavalli et al. (2013).

Returned by ``sediment_connectivity_index()``; not constructible directly.

Attributes
----------
downstream : Raster
index : Raster
upstream : Raster
weight : Raster

Authors
-------
B. Gailleton

### Members

#### `downstream`

(self) -> lsdtt3.Raster

#### `index`

(self) -> lsdtt3.Raster

#### `upstream`

(self) -> lsdtt3.Raster

#### `weight`

(self) -> lsdtt3.Raster

---

## `DivideMigrationResult`

Across-divide comparison made by :func:`divide_migration`.

``network`` (ChannelNetwork), ``chi`` (ChiCollection) and ``basin_labels``
(Raster) are the intermediate products; ``outlets``, ``units`` and
``pair_summary`` are pandas DataFrames and ``divides`` and
``channel_heads`` are GeoDataFrames.

Returned by ``divide_migration()``; not constructible directly.

Attributes
----------
basin_labels : Raster
channel_heads : geopandas.GeoDataFrame
    Channel heads near the divides with their metrics, as points.
chi : ChiCollection
divides : geopandas.GeoDataFrame
    Divide segments between adjacent units, as lines.
network : ChannelNetwork
outlets : pandas.DataFrame
    One row per input outlet with its snapped position.
pair_summary : pandas.DataFrame
    One row per pair of adjacent units: across-divide differences and their uncertainties.
units : pandas.DataFrame
    One row per unit (group of outlets).

Authors
-------
B. Gailleton

### Members

#### `basin_labels`

(self) -> lsdtt3.Raster

#### `channel_heads`

Channel heads near the divides with their metrics, as points.

#### `chi`

(self) -> lsdtt3.ChiCollection

#### `divides`

Divide segments between adjacent units, as lines.

#### `network`

(self) -> lsdtt3.ChannelNetwork

#### `outlets`

One row per input outlet with its snapped position.

#### `pair_summary`

One row per pair of adjacent units: across-divide differences and their uncertainties.

#### `units`

One row per unit (group of outlets).

---

## `FlowGraph`

Fractional D-infinity or MFD flow-routing graph.

Use ``FlowInfo`` for single-receiver D8 channel trees. Use ``FlowGraph`` for
distributed-flow drainage area and later TWI/connectivity calculations.

Construct it as

    FlowGraph(dem: Raster, method: str = 'mfd', boundary_conditions: str = 'oooo', interior_nodata_is_outlet: bool = True, mfd_exponent: float = 1.0, n_threads: int = 0)

Parameters
----------
dem : Raster
    Depression-conditioned elevation raster (m).
method : {"mfd", "dinf"}, default="mfd"
    ``mfd``: flow split among all downslope neighbours weighted by
    slope**mfd_exponent (Quinn et al. 1991; Freeman 1991). ``dinf``:
    D-infinity (Tarboton 1997), flow split between the two neighbours of
    the steepest triangular facet.
boundary_conditions : str, default='oooo'
    Four NESW edge codes, as for ``FlowInfo``.
interior_nodata_is_outlet : bool, default=True
    Pixels next to nodata may drain into it.
mfd_exponent : float, default=1.0
    Exponent on slope in the MFD split (1.0 Quinn, 1.1 Freeman).
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Attributes
----------
donor_count : numpy.ndarray
    Uint8 array: number of cells draining into each cell.
metadata : RasterMeta
method : str
n_valid_cells : int
outlets : numpy.ndarray
    Outlet pixels as an N x 2 (row, col) array.
shape : tuple
valid_mask : numpy.ndarray
    Bool array (n_rows, n_cols), True for valid pixels.

Methods
-------
accumulate(weights: object | None = None) -> Raster
    Accumulate weights through fractional receivers; default weights are cell areas.
direction() -> Raster
drainage_area() -> Raster
    Return the distributed-flow drainage area (m^2) as a Raster.
is_outlet(row: int, col: int) -> bool
    Return True if the pixel (row, col) is an outlet.
outlet_mask() -> Raster
    Return a 0/1 Raster of outlet pixels.
receivers_of(row: int, col: int) -> dict
    Return fractional receiver pixels and weights.
slope() -> Raster

Authors
-------
B. Gailleton

```python
def __init__(self, dem: lsdtt3.Raster, method: str = 'mfd', boundary_conditions: str = 'oooo', interior_nodata_is_outlet: bool = True, mfd_exponent: float = 1.0, n_threads: int = 0) -> None
```

### Members

#### `__init__`

__init__(self, dem: lsdtt3.Raster, method: str = 'mfd', boundary_conditions: str = 'oooo', interior_nodata_is_outlet: bool = True, mfd_exponent: float = 1.0, n_threads: int = 0) -> None

```python
def __init__(self, dem: lsdtt3.Raster, method: str = 'mfd', boundary_conditions: str = 'oooo', interior_nodata_is_outlet: bool = True, mfd_exponent: float = 1.0, n_threads: int = 0) -> None
```

#### `accumulate`

accumulate(self, weights: object | None = None) -> lsdtt3.Raster

Accumulate weights through fractional receivers; default weights are cell areas.

```python
def accumulate(self, weights: object | None = None) -> lsdtt3.Raster
```

#### `direction`

direction(self) -> lsdtt3.Raster

```python
def direction(self) -> lsdtt3.Raster
```

#### `donor_count`

Uint8 array: number of cells draining into each cell.

#### `drainage_area`

drainage_area(self) -> lsdtt3.Raster

Return the distributed-flow drainage area (m^2) as a Raster.

```python
def drainage_area(self) -> lsdtt3.Raster
```

#### `is_outlet`

is_outlet(self, row: int, col: int) -> bool

Return True if the pixel (row, col) is an outlet.

```python
def is_outlet(self, row: int, col: int) -> bool
```

#### `metadata`

(self) -> lsdtt3.RasterMeta

#### `method`

(self) -> str

#### `n_valid_cells`

(self) -> int

#### `outlet_mask`

outlet_mask(self) -> lsdtt3.Raster

Return a 0/1 Raster of outlet pixels.

```python
def outlet_mask(self) -> lsdtt3.Raster
```

#### `outlets`

Outlet pixels as an N x 2 (row, col) array.

#### `receivers_of`

receivers_of(self, row: int, col: int) -> dict

Return fractional receiver pixels and weights.

```python
def receivers_of(self, row: int, col: int) -> dict
```

#### `shape`

(self) -> tuple

#### `slope`

slope(self) -> lsdtt3.Raster

```python
def slope(self) -> lsdtt3.Raster
```

#### `valid_mask`

Bool array (n_rows, n_cols), True for valid pixels.

---

## `FlowInfo`

D8 drainage topology (steepest-descent receivers) of a conditioned DEM.

Each valid pixel drains to the one of its eight neighbours with the steepest
downhill gradient, or is an outlet. ``FlowInfo`` stores this routing (and
each pixel's upslope pixels), not the elevations, so later calls that need
elevations take the DEM again. Fill or breach depressions first (for example
``fill_pits``); a pixel with no downhill neighbour and no allowed exit raises
an error unless ``allow_pits_as_outlets=True``. Pixels are addressed by
``(row, col)`` (row 0 is the northern edge) or by flat index
``row * n_cols + col``.

Typical workflow::

    F = lsdtt3.FlowInfo(lsdtt3.fill_pits(dem))
    scope = F.channel_scope(threshold_contributing_pixels=1000)
    N = F.channel_network(scope)
    B = F.basins(outlets, min_upslope_pixels=0)

Construct it as

    FlowInfo(dem: Raster, boundary_conditions: str = 'oooo', interior_nodata_is_outlet: bool = True, forced_outlets: object | None = None, allow_pits_as_outlets: bool = False, n_threads: int = 0)

Parameters
----------
dem : Raster
    Depression-conditioned elevation raster (m) in a projected metre CRS.
boundary_conditions : str, default='oooo'
    Four characters for the north, east, south and west edges: ``o`` open
    (an edge pixel may drain out of the grid when it has no lower in-grid
    neighbour), ``n`` no-flux (flow cannot leave), ``p`` periodic (wraps to
    the opposite edge; must be given for both opposite edges), ``b``
    base level (every pixel on the edge is an outlet). See
    ``make_boundary_conditions``.
interior_nodata_is_outlet : bool, default=True
    If True, a pixel next to a nodata pixel may drain into it (become an
    outlet) when it has no lower valid neighbour, as at a coastline or the
    edge of a clipped DEM. If False such pixels must drain over valid pixels.
forced_outlets : array-like of int or None, default=None
    N x 2 array of ``(row, col)`` pixels that are made outlets regardless of
    their neighbours (moulins, swallow holes, playas). Prepare the DEM with
    ``fill_pits_with_sinks`` using the same pixels. Out-of-grid or nodata
    pixels raise ``ValueError``; duplicates are ignored.
allow_pits_as_outlets : bool, default=False
    If True, a pixel that has no downhill neighbour and no allowed exit
    becomes an outlet instead of raising an error (useful in landscape
    evolution runs where transient pits form).
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Attributes
----------
boundary_conditions : str
    The four-character NESW boundary-condition string used.
direction : numpy.ndarray
    Uint8 array (n_rows, n_cols) of internal D8 direction codes; see flow_direction() for a raster.
donor_count : numpy.ndarray
    Uint8 array (n_rows, n_cols): number of neighbours draining directly into each pixel.
metadata : RasterMeta
    Grid metadata (RasterMeta) of the DEM the topology was built from.
n_outlets : int
    Number of outlet pixels (pixels with no receiver).
n_valid_cells : int
    Number of valid (routed) pixels.
outlets : numpy.ndarray
    Outlet pixels as an N x 2 int64 array of (row, col), in flat-index order.
receiver : numpy.ndarray
    Int64 array (n_rows, n_cols) of each pixel's receiver flat index; -2 at outlets, -1 at nodata.
shape : tuple[int, int]
    Grid shape as (n_rows, n_cols).
valid_mask : numpy.ndarray
    Bool array (n_rows, n_cols), True for valid pixels.

Methods
-------
accumulate(weights: object | None = None) -> Raster
    Accumulate per-pixel weights downstream; default weights are cell areas.
basins(outlets: object | None = None, min_upslope_pixels: int = 0) -> BasinCollection
    Create a ``BasinCollection`` from explicit outlet pixels or all qualifying D8 outlets.
channel_mask(threshold_contributing_pixels: int) -> Raster
    Return a 0/1 raster where contributing pixels meet the threshold.
channel_network(scope: ChannelScope) -> ChannelNetwork
    Build a compact junction graph from a ChannelScope.
channel_point_attributes(filled_dem: Raster, points: geopandas.GeoDataFrame, *, threshold_contributing_pixels: int, slope_window_half_length_metres: float = 25.0, channel_points_snap_max_search_radius_pixels: int = 10) -> geopandas.GeoDataFrame
    Snap points (e.g. field sites) to the nearest channel cell and add its drainage area and along-channel slope.
channel_scope(threshold_contributing_pixels: int, sources: object | None = None, outlets: object | None = None, snap_sources: bool = True, source_snap_max_radius_pixels: int = 25, outlet_snap_max_distance_metres: float = 0.0, auto_detect_roles: bool = False, channel_definition_algorithm_choice: str = 'contributing_pixels', min_channel_stream_order: int = 3) -> ChannelScope
    Resolve threshold or supplied source paths into a ChannelScope.
chi_from_outlets(outlets: object | None = None, m_over_n: float = 0.45, reference_area: float = 1.0, denominator: object | None = None) -> Raster
    Integrate zero-based chi across complete upstream raster domains.
chi_from_outlets_inherited(outlets: object, seed_chi: Raster, m_over_n: float = 0.45, reference_area: float = 1.0, denominator: object | None = None) -> Raster
    Integrate chi from supplied outlets using seed-raster base values.
contributing_pixels() -> Raster
    Return the number of pixels draining through each pixel, including itself.
distance_to_channel(threshold_contributing_pixels: int) -> Raster
    Return the D8 flow-path length from each pixel to the nearest downstream channel.
distance_to_farthest_ridge() -> Raster
    Return the longest D8 flow-path length from any upslope pixel to this pixel.
downstream_sum(values: Raster, sinks: Raster) -> Raster
    Sum a raster along each pixel's D8 path down to the first sink pixel.
drainage_area() -> Raster
    Return the D8 drainage area of every pixel, including the pixel itself.
flat_index(row: int, col: int) -> int
    Return the flat index ``row * n_cols + col`` of a pixel.
flow_direction() -> Raster
    Return the D8 receiver direction of every pixel as a Raster of codes.
flow_distance_to_outlet() -> Raster
    Return the D8 flow-path length from each pixel to its outlet.
flow_lines_geodataframe() -> object
    Return all D8 receiver edges as a GeoPandas line layer.
height_above_channel(dem: Raster, channels: Raster) -> Raster
    Compute HAND (height above nearest drainage) from an explicit channel raster.
height_above_channel_with_donor(dem: Raster, channels: Raster) -> dict
    Compute HAND and the channel pixel ("donor") each pixel drains to.
height_to_farthest_ridge(dem: Raster) -> Raster
    Return the elevation difference to the farthest upslope ridge pixel.
is_outlet(row: int, col: int) -> bool
    Return True if the pixel is an outlet (has no downstream receiver).
is_upslope(candidate_row: int, candidate_col: int, target_row: int, target_col: int) -> bool
    Return True if the candidate pixel drains through the target pixel.
receiver_of(row: int, col: int) -> object
    Return the pixel this pixel drains to.
row_col(flat_index: int) -> tuple
    Return the ``(row, col)`` of a flat index.
scope_from_heads(heads: ChannelHeads, outlets: object | None = None, outlet_snap_max_distance_metres: float = 0.0) -> ChannelScope
    Resolve ChannelHeads into explicit D8 source paths for channel and chi analysis.
snap_points(x: object | None = None, y: object | None = None, *, longitude: object | None = None, latitude: object | None = None, points: object | None = None, crs: str | None = None, method: str = 'downstream', threshold_contributing_pixels: int | None = None, threshold_drainage_area_m2: float | None = None, min_search_radius_pixels: int = 0, max_search_radius_pixels: int = 25, max_search_radius_metres: float | None = None) -> object
    Snap points (for example basin outlets or sample sites) onto the flow network.
sources(threshold_contributing_pixels: int) -> object
    Return channel sources (channel heads) as an N by 2 (row, col) array.
upslope_mean(values: Raster, include_self: bool = True) -> Raster
    Average a raster over the upslope area of every pixel.
upslope_pixels(row: int, col: int) -> object
    Return every pixel that drains through this pixel, including itself.
upslope_sum(values: Raster, include_self: bool = True) -> Raster
    Sum a raster over the upslope area of every pixel.

Authors
-------
B. Gailleton
S. M. Mudd

**Deprecated compatibility names:**
- `threshold_pixels` -> `threshold_contributing_pixels`; earliest removal: following minor.

```python
def __init__(self, dem: lsdtt3.Raster, boundary_conditions: str = 'oooo', interior_nodata_is_outlet: bool = True, forced_outlets: object | None = None, allow_pits_as_outlets: bool = False, n_threads: int = 0) -> None
```

### Members

#### `__init__`

__init__(self, dem: lsdtt3.Raster, boundary_conditions: str = 'oooo', interior_nodata_is_outlet: bool = True, forced_outlets: object | None = None, allow_pits_as_outlets: bool = False, n_threads: int = 0) -> None

Build D8 flow topology from a depression-conditioned DEM.

``boundary_conditions`` is a north/east/south/west string using open (``o``),
no-flux (``n``), periodic (``p``), or base-level (``b``) edges. Forced outlets
are an N by 2 array of ``(row, col)`` surface-flow terminations.

```python
def __init__(self, dem: lsdtt3.Raster, boundary_conditions: str = 'oooo', interior_nodata_is_outlet: bool = True, forced_outlets: object | None = None, allow_pits_as_outlets: bool = False, n_threads: int = 0) -> None
```

#### `accumulate`

accumulate(self, weights: object | None = None) -> lsdtt3.Raster

Accumulate per-pixel weights downstream; default weights are cell areas.

Parameters
----------
weights : Raster or None, default=None
    Per-pixel values on the flow grid (same grid and validity mask) added to
    every downstream pixel, for example runoff depth x area (m^3) or
    sediment supply. Nodata weights count as 0. None uses each pixel's area
    (``cell_size**2``, m^2), which gives drainage area.

Returns
-------
Raster
    Sum of the weights of each pixel and everything upslope of it, in the
    units of ``weights`` (m^2 by default).

```python
def accumulate(self, weights: object | None = None) -> lsdtt3.Raster
```

#### `basins`

basins(self, outlets: object | None = None, min_upslope_pixels: int = 0) -> lsdtt3.BasinCollection

Create a ``BasinCollection`` from explicit outlet pixels or all qualifying D8 outlets.

Explicit outlets may be any valid pixels, including channel pixels inside
another basin, which creates nested basins.

Parameters
----------
outlets : array-like of int or None, default=None
    N x 2 integer array of ``(row, col)`` outlet pixels; basin ``i`` is the
    upslope area of the ``i``-th distinct outlet. Duplicate outlets are
    silently dropped (so IDs follow first appearance). Out-of-grid or
    nodata pixels raise ``ValueError``. None selects every D8 outlet of the
    flow topology (pixels that drain off the grid or into nodata).
min_upslope_pixels : int, default=0
    With ``outlets=None``, keep only outlets with at least this many
    contributing pixels (including the outlet). Ignored when ``outlets`` is
    given. Non-negative.

Returns
-------
BasinCollection

```python
def basins(self, outlets: object | None = None, min_upslope_pixels: int = 0) -> lsdtt3.BasinCollection
```

#### `boundary_conditions`

The four-character NESW boundary-condition string used.

#### `channel_mask`

channel_mask(self, threshold_contributing_pixels: int) -> lsdtt3.Raster

Return a 0/1 raster where contributing pixels meet the threshold.

Parameters
----------
threshold_contributing_pixels : int
    A pixel is a channel (1) when its contributing-pixel count, including
    itself, is >= this value; must be positive. Area threshold in m^2 is
    this times ``metadata.pixel_area``.

Returns
-------
Raster
    1 for channel pixels, 0 elsewhere; nodata where the DEM is nodata.

```python
def channel_mask(self, threshold_contributing_pixels: int) -> lsdtt3.Raster
```

#### `channel_network`

channel_network(self, scope: lsdtt3.ChannelScope) -> lsdtt3.ChannelNetwork

Build a compact junction graph from a ChannelScope.

Parameters
----------
scope : ChannelScope
    Made by this ``FlowInfo``'s ``channel_scope`` (or ``scope_from_heads``);
    must have at least one resolved source. Each resolved source is traced
    to its D8 outlet.

Returns
-------
ChannelNetwork

```python
def channel_network(self, scope: lsdtt3.ChannelScope) -> lsdtt3.ChannelNetwork
```

#### `channel_point_attributes`

channel_point_attributes(self, filled_dem: Raster, points: geopandas.GeoDataFrame, *, threshold_contributing_pixels: int, slope_window_half_length_metres: float = 25.0, channel_points_snap_max_search_radius_pixels: int = 10) -> geopandas.GeoDataFrame

Snap points (e.g. field sites) to the nearest channel cell and add its drainage area and along-channel slope.

lsdtt-channel-extraction ``run_channel_point_attributes``, with the same
C++ routine and parameter names. Each point moves to the channel cell
(at least ``threshold_contributing_pixels`` contributing pixels) whose
centre is nearest in straight-line distance, within
``channel_points_snap_max_search_radius_pixels``; it does not follow the
D8 path. The slope is the least-squares slope of ``filled_dem`` against
distance along the channel over +/- ``slope_window_half_length_metres``;
upstream the window follows the largest-area tributary.

Parameters
----------
filled_dem : Raster
    The DEM this FlowInfo was built on (pit-filled); elevations and the
    slope are taken from it.
points : geopandas.GeoDataFrame
    Point geometries with a CRS (reprojected to the grid's). Every column
    is kept.
threshold_contributing_pixels : int, keyword-only
    Channel cells have at least this many contributing pixels; also the
    least a tributary must carry for the slope window to follow it.
slope_window_half_length_metres : float, default=25.0
    Half-length (m) of the along-channel slope window; at least two cell
    sizes.
channel_points_snap_max_search_radius_pixels : int, default=10
    Furthest a point may move, in pixels (Chebyshev rings).

Returns
-------
geopandas.GeoDataFrame
    ``points`` in the grid CRS with the geometry moved to the snapped cell
    centre, and the columns original_x, original_y (grid CRS),
    snap_distance_m, snap_status, snapped_flat, contributing_pixels,
    drainage_area_m2, flow_distance_m, elevation_m, smoothed_elevation_m,
    slope_m_per_m, slope_valid (0/1), slope_window_length_m. A point that
    does not snap keeps its position, with its snap_status and empty (NaN /
    <NA>) values. An added column whose name is already in ``points`` is
    renamed with the suffix ``_lsdtt`` (with a warning), as in the CLI.

Examples
--------
>>> filled = lsdtt3.fill_pits(dem)
>>> F = lsdtt3.FlowInfo(filled)
>>> sites = geopandas.read_file("field_sites.gpkg")
>>> out = F.channel_point_attributes(filled, sites, threshold_contributing_pixels=500,
...                                  slope_window_half_length_metres=150.0)

```python
def channel_point_attributes(self, filled_dem: Raster, points: geopandas.GeoDataFrame, *, threshold_contributing_pixels: int, slope_window_half_length_metres: float = 25.0, channel_points_snap_max_search_radius_pixels: int = 10) -> geopandas.GeoDataFrame
```

#### `channel_scope`

channel_scope(self, threshold_contributing_pixels: int, sources: object | None = None, outlets: object | None = None, snap_sources: bool = True, source_snap_max_radius_pixels: int = 25, outlet_snap_max_distance_metres: float = 0.0, auto_detect_roles: bool = False, channel_definition_algorithm_choice: str = 'contributing_pixels', min_channel_stream_order: int = 3) -> lsdtt3.ChannelScope

Resolve threshold or supplied source paths into a ChannelScope.

With no ``sources``, every pixel whose contributing-pixel count reaches the
threshold is a channel. With ``sources``, only the channels traced
downstream from those points are used, optionally stopped at ``outlets``.

Parameters
----------
threshold_contributing_pixels : int
    Channel threshold: a pixel is a channel when the number of pixels
    draining through it (including itself) is >= this value; positive.
    Threshold area in m^2 is this times ``metadata.pixel_area`` (for
    example 1000 pixels on a 30 m grid is 0.9 km^2).
sources : array-like of int or None, default=None
    N x 2 array of ``(row, col)`` channel-head pixels (valid pixels). Each
    is traced downstream to its D8 outlet. None selects all channels.
outlets : array-like of int or None, default=None
    N x 2 array of ``(row, col)`` points where traced channels stop; needs
    ``sources``. Each outlet is matched to the nearest pixel on a traced
    path, and that path is cut there.
snap_sources : bool, default=True
    Move each source to the nearest pixel that meets the threshold
    (expanding-square search); sources that cannot be snapped are dropped.
    If False sources are used as given.
source_snap_max_radius_pixels : int, default=25
    Largest search radius, in pixels, for source snapping.
outlet_snap_max_distance_metres : float, default=0.0
    Largest distance (m) from an outlet point to the traced path it cuts;
    outlets farther away are ignored. 0 or negative means unlimited.
auto_detect_roles : bool, default=False
    Pass all points as ``sources`` and let roles be inferred: a point that
    lies downstream of another supplied point becomes an outlet that cuts
    that channel. Requires ``outlets=None``; sources are always snapped.
channel_definition_algorithm_choice : {"contributing_pixels", "stream_order"}, default="contributing_pixels"
    What counts as a channel pixel, as in lsdtt-fluvial-surfaces:
    ``contributing_pixels`` keeps every traced pixel; ``stream_order``
    keeps only traced pixels of Strahler order >= ``min_channel_stream_order``
    (Clubb et al. 2017 use 3). The threshold still decides which
    first-order channels exist, and so the order every pixel reaches.
min_channel_stream_order : int, default=3
    Lowest Strahler order kept under ``stream_order``.

Returns
-------
ChannelScope
    Selected channel pixels, resolved sources and outlets.

```python
def channel_scope(self, threshold_contributing_pixels: int, sources: object | None = None, outlets: object | None = None, snap_sources: bool = True, source_snap_max_radius_pixels: int = 25, outlet_snap_max_distance_metres: float = 0.0, auto_detect_roles: bool = False, channel_definition_algorithm_choice: str = 'contributing_pixels', min_channel_stream_order: int = 3) -> lsdtt3.ChannelScope
```

#### `chi_from_outlets`

chi_from_outlets(self, outlets: object | None = None, m_over_n: float = 0.45, reference_area: float = 1.0, denominator: object | None = None) -> lsdtt3.Raster

Integrate zero-based chi across complete upstream raster domains.

Unlike ``ChannelNetwork.chi``, chi is computed on every pixel upslope of
each outlet, not only on channels.

Parameters
----------
outlets : array-like of int or None, default=None
    N x 2 array of ``(row, col)`` outlet pixels; None uses every outlet.
m_over_n : float, default=0.45
    Concavity index theta (positive).
reference_area : float, default=1.0
    Reference drainage area A0 in m^2 (positive).
denominator : Raster or None, default=None
    Replacement for drainage area A (m^2) in the integrand; None uses D8
    drainage area.

Returns
-------
Raster
    Chi (m), 0 at each outlet; nodata outside the selected basins.

```python
def chi_from_outlets(self, outlets: object | None = None, m_over_n: float = 0.45, reference_area: float = 1.0, denominator: object | None = None) -> lsdtt3.Raster
```

#### `chi_from_outlets_inherited`

chi_from_outlets_inherited(self, outlets: object, seed_chi: lsdtt3.Raster, m_over_n: float = 0.45, reference_area: float = 1.0, denominator: object | None = None) -> lsdtt3.Raster

Integrate chi from supplied outlets using seed-raster base values.

Each outlet starts at the value of ``seed_chi`` at that pixel instead of 0,
so chi can be continued from a larger (for example coarser) domain.

Parameters
----------
outlets : array-like of int
    N x 2 array of ``(row, col)`` outlet pixels (non-empty).
seed_chi : Raster
    Chi (m) on the flow grid giving each outlet's starting value.
m_over_n, reference_area, denominator
    As in ``chi_from_outlets``.

Returns
-------
Raster
    Chi (m) upslope of the outlets.

```python
def chi_from_outlets_inherited(self, outlets: object, seed_chi: lsdtt3.Raster, m_over_n: float = 0.45, reference_area: float = 1.0, denominator: object | None = None) -> lsdtt3.Raster
```

#### `contributing_pixels`

contributing_pixels(self) -> lsdtt3.Raster

Return the number of pixels draining through each pixel, including itself.

Returns
-------
Raster
    Pixel counts (1 at a ridge or headwater pixel). Multiply by
    ``metadata.pixel_area`` for area in m^2, or use ``drainage_area()``.

```python
def contributing_pixels(self) -> lsdtt3.Raster
```

#### `direction`

Uint8 array (n_rows, n_cols) of internal D8 direction codes; see flow_direction() for a raster.

#### `distance_to_channel`

distance_to_channel(self, threshold_contributing_pixels: int) -> lsdtt3.Raster

Return the D8 flow-path length from each pixel to the nearest downstream channel.

Parameters
----------
threshold_contributing_pixels : int
    Channel threshold in contributing pixels (including self).

Returns
-------
Raster
    Distance in metres; 0 on channel pixels; nodata where no channel is
    reached downstream.

```python
def distance_to_channel(self, threshold_contributing_pixels: int) -> lsdtt3.Raster
```

#### `distance_to_farthest_ridge`

distance_to_farthest_ridge(self) -> lsdtt3.Raster

Return the longest D8 flow-path length from any upslope pixel to this pixel.

Returns
-------
Raster
    Distance in metres; 0 at headwater pixels.

```python
def distance_to_farthest_ridge(self) -> lsdtt3.Raster
```

#### `donor_count`

Uint8 array (n_rows, n_cols): number of neighbours draining directly into each pixel.

#### `downstream_sum`

downstream_sum(self, values: lsdtt3.Raster, sinks: lsdtt3.Raster) -> lsdtt3.Raster

Sum a raster along each pixel's D8 path down to the first sink pixel.

Parameters
----------
values : Raster
    Values on the flow grid summed along the path (the pixel itself
    included, the sink excluded).
sinks : Raster
    0/1 raster on the flow grid; pixels > 0.5 stop the sum (for example a
    channel mask).

Returns
-------
Raster
    Path sums in the units of ``values``; nodata on sink pixels.

```python
def downstream_sum(self, values: lsdtt3.Raster, sinks: lsdtt3.Raster) -> lsdtt3.Raster
```

#### `drainage_area`

drainage_area(self) -> lsdtt3.Raster

Return the D8 drainage area of every pixel, including the pixel itself.

Returns
-------
Raster
    Drainage area in m^2: ``contributing_pixels * cell_size**2``.

```python
def drainage_area(self) -> lsdtt3.Raster
```

#### `flat_index`

flat_index(self, row: int, col: int) -> int

Return the flat index ``row * n_cols + col`` of a pixel.

Parameters
----------
row, col : int
    Zero-based pixel position; row 0 is the northern edge.

Returns
-------
int
    Flat index, as used in the ``flat`` columns of lsdtt3 tables. Raises
    ``RuntimeError`` if the pixel is outside the grid.

```python
def flat_index(self, row: int, col: int) -> int
```

#### `flow_direction`

flow_direction(self) -> lsdtt3.Raster

Return the D8 receiver direction of every pixel as a Raster of codes.

Codes are ``0`` north, ``1`` north-east, ``2`` east, ``3`` south-east,
``4`` south, ``5`` south-west, ``6`` west, ``7`` north-west (clockwise from
north). Outlet pixels carry the code of the edge they drain across.

Returns
-------
Raster
    Direction codes 0-7; nodata where the DEM is nodata.

```python
def flow_direction(self) -> lsdtt3.Raster
```

#### `flow_distance_to_outlet`

flow_distance_to_outlet(self) -> lsdtt3.Raster

Return the D8 flow-path length from each pixel to its outlet.

Returns
-------
Raster
    Distance in metres along the D8 path (cardinal steps ``cell_size``,
    diagonal steps ``cell_size * sqrt(2)``); 0 at outlets.

```python
def flow_distance_to_outlet(self) -> lsdtt3.Raster
```

#### `flow_lines_geodataframe`

flow_lines_geodataframe(self) -> object

Return all D8 receiver edges as a GeoPandas line layer.

Returns
-------
geopandas.GeoDataFrame
    One two-point LineString (pixel centre to receiver centre) per
    non-outlet pixel, in the DEM CRS. Columns ``flat``, ``receiver_flat``,
    ``row``, ``col``, ``receiver_row``, ``receiver_col``, ``length_metres``,
    ``geometry``. Large for big DEMs.

```python
def flow_lines_geodataframe(self) -> object
```

#### `height_above_channel`

height_above_channel(self, dem: lsdtt3.Raster, channels: lsdtt3.Raster) -> lsdtt3.Raster

Compute HAND (height above nearest drainage) from an explicit channel raster.

HAND (Nobre et al. 2011) is the elevation of each pixel minus the elevation
of the first channel pixel reached by following D8 flow downstream.

Parameters
----------
dem : Raster
    Elevation (m) on the flow grid (same grid and validity mask).
channels : Raster
    0/1 raster on the flow grid marking channel pixels (values > 0.5), for
    example ``ChannelNetwork.channel_mask`` or ``channel_mask(threshold)``.

Returns
-------
Raster
    Height above the channel in metres; 0 on channel pixels; nodata where
    the flow path reaches an outlet without meeting a channel.

```python
def height_above_channel(self, dem: lsdtt3.Raster, channels: lsdtt3.Raster) -> lsdtt3.Raster
```

#### `height_above_channel_with_donor`

height_above_channel_with_donor(self, dem: lsdtt3.Raster, channels: lsdtt3.Raster) -> dict

Compute HAND and the channel pixel ("donor") each pixel drains to.

Parameters
----------
dem : Raster
    Elevation (m) on the flow grid.
channels : Raster
    0/1 channel raster on the flow grid (values > 0.5 are channel).

Returns
-------
dict
    ``"height"``: HAND Raster in metres (as ``height_above_channel``);
    ``"channel_donor"``: int64 array (n_rows, n_cols) of the flat index of
    the channel pixel that sets each pixel's HAND, -1 where undefined.

```python
def height_above_channel_with_donor(self, dem: lsdtt3.Raster, channels: lsdtt3.Raster) -> dict
```

#### `height_to_farthest_ridge`

height_to_farthest_ridge(self, dem: lsdtt3.Raster) -> lsdtt3.Raster

Return the elevation difference to the farthest upslope ridge pixel.

The farthest ridge pixel is the pixel in the upslope area with the longest
D8 flow path to this pixel.

Parameters
----------
dem : Raster
    Elevation (m) on the flow grid.

Returns
-------
Raster
    Ridge elevation minus pixel elevation, in metres.

```python
def height_to_farthest_ridge(self, dem: lsdtt3.Raster) -> lsdtt3.Raster
```

#### `is_outlet`

is_outlet(self, row: int, col: int) -> bool

Return True if the pixel is an outlet (has no downstream receiver).

Parameters
----------
row, col : int
    Zero-based pixel position of a valid pixel.

Returns
-------
bool
    Raises ``RuntimeError`` for an out-of-grid or nodata pixel.

```python
def is_outlet(self, row: int, col: int) -> bool
```

#### `is_upslope`

is_upslope(self, candidate_row: int, candidate_col: int, target_row: int, target_col: int) -> bool

Return True if the candidate pixel drains through the target pixel.

Parameters
----------
candidate_row, candidate_col : int
    Pixel to test.
target_row, target_col : int
    Pixel whose upslope area is searched. A pixel is upslope of itself.

Returns
-------
bool

```python
def is_upslope(self, candidate_row: int, candidate_col: int, target_row: int, target_col: int) -> bool
```

#### `metadata`

Grid metadata (RasterMeta) of the DEM the topology was built from.

#### `n_outlets`

Number of outlet pixels (pixels with no receiver).

#### `n_valid_cells`

Number of valid (routed) pixels.

#### `outlets`

Outlet pixels as an N x 2 int64 array of (row, col), in flat-index order.

#### `receiver`

Int64 array (n_rows, n_cols) of each pixel's receiver flat index; -2 at outlets, -1 at nodata.

#### `receiver_of`

receiver_of(self, row: int, col: int) -> object

Return the pixel this pixel drains to.

Parameters
----------
row, col : int
    Zero-based pixel position of a valid pixel.

Returns
-------
tuple of (int, int) or None
    ``(row, col)`` of the D8 downstream neighbour, or None if the pixel is an
    outlet. Raises ``RuntimeError`` for an out-of-grid or nodata pixel.

```python
def receiver_of(self, row: int, col: int) -> object
```

#### `row_col`

row_col(self, flat_index: int) -> tuple

Return the ``(row, col)`` of a flat index.

Parameters
----------
flat_index : int
    Flat index ``row * n_cols + col`` (for example a ``flat`` or
    ``outlet_flat`` table value).

Returns
-------
tuple of (int, int)
    ``(row, col)``; equivalent to ``divmod(flat_index, n_cols)``. Raises
    ``RuntimeError`` if the index is out of range.

```python
def row_col(self, flat_index: int) -> tuple
```

#### `scope_from_heads`

scope_from_heads(self, heads: lsdtt3.ChannelHeads, outlets: object | None = None, outlet_snap_max_distance_metres: float = 0.0) -> lsdtt3.ChannelScope

Resolve ChannelHeads into explicit D8 source paths for channel and chi analysis.

Equivalent to ``channel_scope`` with ``sources=heads.pixels`` and no source
snapping.

Parameters
----------
heads : ChannelHeads
    Heads extracted with this ``FlowInfo``.
outlets : array-like of int or None, default=None
    N x 2 array of ``(row, col)`` points where traced channels stop.
outlet_snap_max_distance_metres : float, default=0.0
    Largest distance (m) from an outlet to the path it cuts; 0 = unlimited.

Returns
-------
ChannelScope

```python
def scope_from_heads(self, heads: lsdtt3.ChannelHeads, outlets: object | None = None, outlet_snap_max_distance_metres: float = 0.0) -> lsdtt3.ChannelScope
```

#### `shape`

Grid shape as (n_rows, n_cols).

#### `snap_points`

snap_points(self, x: object | None = None, y: object | None = None, *, longitude: object | None = None, latitude: object | None = None, points: object | None = None, crs: str | None = None, method: str = 'downstream', threshold_contributing_pixels: int | None = None, threshold_drainage_area_m2: float | None = None, min_search_radius_pixels: int = 0, max_search_radius_pixels: int = 25, max_search_radius_metres: float | None = None) -> object

Snap points (for example basin outlets or sample sites) onto the flow network.

Each point is placed in the pixel that contains it and then moved by
``method``. The result is the ``(row, col)`` array that
``FlowInfo.basins(outlets)``, ``ChannelNetwork.chi(outlets=...)`` and other
``outlets`` arguments take. Give the points as exactly one of ``x`` and
``y``, ``longitude`` and ``latitude``, or ``points``.

The methods, their names and their rules are those of the command-line
programs' ``basin_outlet_snapping_method`` (the same C++ routines): a
``method``, ``threshold_contributing_pixels`` and the two radii here give
the same pixel as ``basin_outlet_snapping_method``,
``basin_outlet_contributing_pixels_threshold`` and
``basin_outlet_snap_min/max_search_radius_pixels`` with the same values. A
pixel "qualifies" when its contributing-pixel count (the pixel itself and
all pixels draining through it) is at least the threshold (``>=``).

* ``"downstream"``: follow the D8 receivers from the point's pixel to the
  first qualifying pixel. The point must be on a valid pixel and upstream
  of a qualifying pixel. The radii are not used.
* ``"expanding_square"``: search square rings of half-width
  ``min_search_radius_pixels`` to ``max_search_radius_pixels`` around the
  point's pixel, one ring at a time, and take the nearest qualifying pixel
  (straight-line distance) in the first ring that has one; ties go to the
  larger drainage area, then the smaller flat index.
* ``"largest_area"``: the pixel with the largest drainage area in the
  square of half-width ``max_search_radius_pixels`` (ties: the nearer
  pixel, then the smaller flat index), provided it qualifies. This is the
  GIS "snap pour point" rule. With a wide square it can jump onto a bigger
  neighbouring river, so keep the square to a few pixels.
  ``min_search_radius_pixels`` is not used.

Error policy: if any point is outside the grid or cannot be snapped,
``ValueError`` is raised naming each such point (by its zero-based position
in the input) and the snapping status (``input_out_of_bounds``,
``input_invalid_flow_cell``, ``threshold_not_found_in_search_window`` or
``outlet_reached_before_threshold``); no partial result is returned.

Parameters
----------
x, y : float or array-like of float, optional
    Easting and northing of the points, in the flow grid's CRS
    (``F.metadata.crs_wkt``, metres) unless ``crs`` is given.
longitude, latitude : float or array-like of float, optional, keyword-only
    Point coordinates in decimal degrees, in ``crs`` (default WGS84,
    ``"EPSG:4326"``). Longitude is the x coordinate: pass both by keyword.
points : geopandas.GeoDataFrame or geopandas.GeoSeries, optional, keyword-only
    Point geometries with a CRS; they are transformed to the grid's CRS.
crs : str or None, default=None, keyword-only
    CRS of ``x``/``y`` or ``longitude``/``latitude`` as an ``"EPSG:<code>"``
    string, WKT or PROJ string. None means the grid's CRS for ``x``/``y`` and
    ``"EPSG:4326"`` for ``longitude``/``latitude``. Not used with ``points``.
method : {"downstream", "expanding_square", "largest_area"}, default="downstream", keyword-only
    Snapping rule, as described above.
threshold_contributing_pixels : int or None, default=None, keyword-only
    Contributing pixels (>= 1) a pixel needs to qualify; the CLIs'
    ``basin_outlet_contributing_pixels_threshold``. Give this or
    ``threshold_drainage_area_m2``, not both. With neither, every method
    uses that key's CLI default, 1000 contributing pixels (as a CLI run
    without the key does); pass ``threshold_contributing_pixels=1`` to
    ``"largest_area"`` for no minimum.
threshold_drainage_area_m2 : float or None, default=None, keyword-only
    The same threshold as a drainage area (m^2, > 0). It is converted to the
    smallest whole number of contributing pixels N with
    N * ``F.metadata.pixel_area`` >= ``threshold_drainage_area_m2``, so a
    pixel qualifies exactly when its drainage area (contributing pixels *
    pixel area) is >= ``threshold_drainage_area_m2``.
min_search_radius_pixels : int, default=0, keyword-only
    Smallest ring half-width (pixels, >= 0) for ``"expanding_square"``; 0
    tests the point's own pixel first. As
    ``basin_outlet_snap_min_search_radius_pixels``.
max_search_radius_pixels : int, default=25, keyword-only
    Largest ring half-width for ``"expanding_square"`` and the square's
    half-width for ``"largest_area"`` (pixels, >= 0). As
    ``basin_outlet_snap_max_search_radius_pixels``.
max_search_radius_metres : float or None, default=None, keyword-only
    ``max_search_radius_pixels`` in metres; if given it replaces it and is
    rounded up to whole pixels.

Returns
-------
numpy.ndarray
    N x 2 int64 array of snapped ``(row, col)`` pixels, one row per input
    point and in input order (row 0 is the north edge). Points may snap to
    the same pixel; ``FlowInfo.basins`` then keeps only the first.

Examples
--------
>>> F = lsdtt3.FlowInfo(lsdtt3.fill_pits(dem))
>>> outlets = F.snap_points(longitude=[15.79008, 16.06060], latitude=[37.94234, 38.06143],
...                         method="downstream", threshold_contributing_pixels=500)
>>> B = F.basins(outlets)

```python
def snap_points(self, x: object | None = None, y: object | None = None, *, longitude: object | None = None, latitude: object | None = None, points: object | None = None, crs: str | None = None, method: str = 'downstream', threshold_contributing_pixels: int | None = None, threshold_drainage_area_m2: float | None = None, min_search_radius_pixels: int = 0, max_search_radius_pixels: int = 25, max_search_radius_metres: float | None = None) -> object
```

#### `sources`

sources(self, threshold_contributing_pixels: int) -> object

Return channel sources (channel heads) as an N by 2 (row, col) array.

A source is a channel pixel (contributing pixels >= the threshold) that
has no channel pixel draining into it.

Parameters
----------
threshold_contributing_pixels : int
    Channel threshold in contributing pixels (including self); positive.

Returns
-------
numpy.ndarray
    N x 2 int64 array of ``(row, col)``; usable as ``channel_scope(sources=...)``.

```python
def sources(self, threshold_contributing_pixels: int) -> object
```

#### `upslope_mean`

upslope_mean(self, values: lsdtt3.Raster, include_self: bool = True) -> lsdtt3.Raster

Average a raster over the upslope area of every pixel.

Parameters
----------
values : Raster
    Values on the flow grid; nodata pixels are excluded from the mean.
include_self : bool, default=True
    Include the pixel's own value.

Returns
-------
Raster
    Upslope means, in the units of ``values``; nodata where no valid value
    is upslope.

```python
def upslope_mean(self, values: lsdtt3.Raster, include_self: bool = True) -> lsdtt3.Raster
```

#### `upslope_pixels`

upslope_pixels(self, row: int, col: int) -> object

Return every pixel that drains through this pixel, including itself.

Parameters
----------
row, col : int
    Zero-based pixel position of a valid pixel.

Returns
-------
numpy.ndarray
    N x 2 int64 array of ``(row, col)``; the first row is the pixel itself.
    N is the pixel's contributing-pixel count.

```python
def upslope_pixels(self, row: int, col: int) -> object
```

#### `upslope_sum`

upslope_sum(self, values: lsdtt3.Raster, include_self: bool = True) -> lsdtt3.Raster

Sum a raster over the upslope area of every pixel.

Parameters
----------
values : Raster
    Values on the flow grid; nodata pixels contribute 0.
include_self : bool, default=True
    Include the pixel's own value; if False only strictly upslope pixels
    count (headwaters give 0).

Returns
-------
Raster
    Upslope sums, in the units of ``values``.

```python
def upslope_sum(self, values: lsdtt3.Raster, include_self: bool = True) -> lsdtt3.Raster
```

#### `valid_mask`

Bool array (n_rows, n_cols), True for valid pixels.

---

## `JunctionBasinCollection`

Junction-aware D8 basins selected from a ``ChannelNetwork``.

Each basin drains to one channel junction and carries that junction's
Strahler stream order. The ``filter_*`` methods remove basins in place and
return the collection, so they can be chained:
``junction_basins(N).filter_incomplete().filter_stream_order(2, 2)``.
Created with :func:`junction_basins`.

Returned by ``junction_basins()``, ``JunctionBasinCollection.filter_area()``, ``JunctionBasinCollection.filter_incomplete()``, ``JunctionBasinCollection.filter_largest()`` and 3 more; not constructible directly.

Attributes
----------
labels : Raster
    Raster of each basin's junction ID; nested basins overwrite their parents; nodata outside all basins.
size : int
    Number of selected basins.
table : pandas.DataFrame
    One row per basin (pandas.DataFrame).

Methods
-------
filter_area(min_pixels: int = -1, max_pixels: int = -1) -> JunctionBasinCollection
    Keep basins whose pixel count is within [min_pixels, max_pixels].
filter_incomplete() -> JunctionBasinCollection
    Remove basins that touch nodata or an open DEM edge; returns this collection.
filter_largest() -> JunctionBasinCollection
    Keep only the basin with the most pixels; returns this collection.
filter_nested() -> JunctionBasinCollection
    Remove basins whose outlet lies upslope of another selected basin's outlet; returns this collection.
filter_nested_by_stream_order() -> JunctionBasinCollection
    Remove nested basins only where they share the enclosing basin's stream order; returns this collection.
filter_stream_order(min_order: int = -1, max_order: int = -1) -> JunctionBasinCollection
    Keep basins whose Strahler order is within [min_order, max_order].
mask(basin_id: int) -> Raster
    Return a 0/1 Raster of one basin.

Authors
-------
B. Gailleton

### Members

#### `filter_area`

filter_area(self, min_pixels: int = -1, max_pixels: int = -1) -> lsdtt3.JunctionBasinCollection

Keep basins whose pixel count is within [min_pixels, max_pixels].

Parameters
----------
min_pixels, max_pixels : int, default=-1
    Bounds in pixels (multiply by ``metadata.pixel_area`` for m^2); -1
    leaves that end open.

Returns
-------
JunctionBasinCollection
    This collection, filtered in place.

```python
def filter_area(self, min_pixels: int = -1, max_pixels: int = -1) -> lsdtt3.JunctionBasinCollection
```

#### `filter_incomplete`

filter_incomplete(self) -> lsdtt3.JunctionBasinCollection

Remove basins that touch nodata or an open DEM edge; returns this collection.

```python
def filter_incomplete(self) -> lsdtt3.JunctionBasinCollection
```

#### `filter_largest`

filter_largest(self) -> lsdtt3.JunctionBasinCollection

Keep only the basin with the most pixels; returns this collection.

```python
def filter_largest(self) -> lsdtt3.JunctionBasinCollection
```

#### `filter_nested`

filter_nested(self) -> lsdtt3.JunctionBasinCollection

Remove basins whose outlet lies upslope of another selected basin's outlet; returns this collection.

```python
def filter_nested(self) -> lsdtt3.JunctionBasinCollection
```

#### `filter_nested_by_stream_order`

filter_nested_by_stream_order(self) -> lsdtt3.JunctionBasinCollection

Remove nested basins only where they share the enclosing basin's stream order; returns this collection.

```python
def filter_nested_by_stream_order(self) -> lsdtt3.JunctionBasinCollection
```

#### `filter_stream_order`

filter_stream_order(self, min_order: int = -1, max_order: int = -1) -> lsdtt3.JunctionBasinCollection

Keep basins whose Strahler order is within [min_order, max_order].

Parameters
----------
min_order, max_order : int, default=-1
    Bounds on stream order; -1 leaves that end open.

Returns
-------
JunctionBasinCollection
    This collection, filtered in place.

```python
def filter_stream_order(self, min_order: int = -1, max_order: int = -1) -> lsdtt3.JunctionBasinCollection
```

#### `labels`

Raster of each basin's junction ID; nested basins overwrite their parents; nodata outside all basins.

#### `mask`

mask(self, basin_id: int) -> lsdtt3.Raster

Return a 0/1 Raster of one basin.

Parameters
----------
basin_id : int
    Zero-based basin index (a ``table.basin_id`` value).

Returns
-------
Raster

```python
def mask(self, basin_id: int) -> lsdtt3.Raster
```

#### `size`

Number of selected basins.

#### `table`

One row per basin (pandas.DataFrame).

Columns ``basin_id`` (0-based), ``junction_id``, ``stream_order``,
``outlet_flat``, ``outlet_row``, ``outlet_col``, ``n_pixels``,
``area_metres2`` (m^2), ``complete`` (False if the basin touches nodata or
an open edge).

---

## `SideTributaryBasins`

Side basins draining into a traced D8 channel.

Made by :func:`side_tributary_basins`.

Returned by ``side_tributary_basins()``; not constructible directly.

Attributes
----------
basins : BasinCollection
    The accepted side basins as a BasinCollection.
channel_mask : Raster
    0/1 Raster of the traced trunk channel.
table : pandas.DataFrame
    One row per accepted side basin: outlet_flat, receiver_channel_flat, outlet_row, outlet_col, receiver_channel_row, receiver_channel_col, n_pixels, drainage_area_metres2, complete.

Authors
-------
B. Gailleton

### Members

#### `basins`

The accepted side basins as a BasinCollection.

#### `channel_mask`

0/1 Raster of the traced trunk channel.

#### `table`

One row per accepted side basin: outlet_flat, receiver_channel_flat, outlet_row, outlet_col, receiver_channel_row, receiver_channel_col, n_pixels, drainage_area_metres2, complete.

---

## `basin_channel_lengths`

Return total in-basin D8 channel length by Strahler stream order.

``basins`` and ``network`` must have been created from the same ``FlowInfo``.

Parameters
----------
basins : BasinCollection
    Basins to measure.
network : ChannelNetwork
    Channel network from the same ``FlowInfo``.
max_stream_order : int, default=6
    Highest Strahler order reported (positive).

Returns
-------
pandas.DataFrame
    One row per basin: ``basin_id`` and ``channel_length_order_<k>_metres``
    for k = 1 .. ``max_stream_order`` (D8 channel length in m).

Authors
-------
B. Gailleton

```python
def basin_channel_lengths(basins: lsdtt3.BasinCollection, network: lsdtt3.ChannelNetwork, max_stream_order: int = 6) -> object
```

---

## `basin_surface_metrics`

Summarise polyfit surface metrics inside each basin at several fitting scales.

Parameters
----------
basins : BasinCollection
    Basins to summarise.
dem : Raster
    Elevation (m) on the basins' flow grid.
scales : list of float
    Polyfit length scales ``lambda_metres`` (m, positive), as in
    ``polyfit_metrics``.
metrics : list of str
    Any of ``"slope"`` (m/m), ``"slope_degrees"``, ``"aspect_degrees"``,
    ``"maximum_curvature"``, ``"minimum_curvature"``, ``"mean_curvature"``,
    ``"gaussian_curvature"``, ``"profile_curvature"``, ``"plan_curvature"``,
    ``"tangential_curvature"``.
include_outline : bool, default=True
    Also sample each metric along the basin outline.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.
bearings_relative_to_true_north : bool, default=False
    Convert ``aspect_degrees`` to true north at each pixel before it is
    summarised, as the CLI's ``bearings_relative_to_true_north``; needs a
    DEM with a CRS. False keeps it relative to grid north.

Returns
-------
dict of str to pandas.DataFrame
    ``"interior"`` (basin_id, scale_metres, metric, n_values, min, mean,
    median, max), ``"outline"`` (basin_id, scale_metres, metric,
    distance_metres, value), ``"curvature_flips"`` (sign changes of
    maximum_curvature along the outline), ``"slope_area"`` (median slope in
    0.1 log10 drainage-area bins) and ``"slope_area_summary"``.

Authors
-------
B. Gailleton
S. M. Mudd

```python
def basin_surface_metrics(basins: lsdtt3.BasinCollection, dem: lsdtt3.Raster, scales: collections.abc.Sequence[float], metrics: collections.abc.Sequence[str], include_outline: bool = True, n_threads: int = 0, bearings_relative_to_true_north: bool = False) -> dict
```

---

## `divide`

Divide two grid-compatible rasters (a / b); nodata where b is 0 or either input is nodata.

Parameters
----------
a, b : Raster
    Numerator and denominator rasters on the same grid.

Returns
-------
Raster

Authors
-------
B. Gailleton
S. M. Mudd

```python
def divide(a: lsdtt3.Raster, b: lsdtt3.Raster) -> lsdtt3.Raster
```

---

## `divide_migration`

Compare channel-head elevation, chi, relief, and gradient across adjacent outlet-defined basin divides.

Basins are defined by the supplied outlets (grouped into units by
``clusters``); for every pair of adjacent units the channel heads on either
side of their shared divide are compared, as evidence for the direction of
divide migration.

Parameters
----------
flow : FlowInfo
    Flow topology of the conditioned DEM.
dem : Raster
    The conditioned DEM (m) on the flow grid.
outlets : array-like of int
    N x 2 array of ``(row, col)`` outlet pixels (at least two).
clusters : array-like or None, default=None
    One label per outlet grouping outlets into units (outlets with the same
    label form one side); None makes every outlet its own unit.
source_pixels : int, default=1000
    Channel threshold (contributing pixels) for sources, snapping and chi.
snapping : {"downstream", "none"}, default="downstream"
    Snap outlets downstream onto the channel network, or use them as given.
outlet_control : {"none", "complete_only", "elevation", "max_out_elevation"}, default="complete_only"
    ``complete_only``: drop basins touching nodata or an open edge;
    ``elevation``: drop outlets below ``min_outlet_elevation``;
    ``max_out_elevation``: put every basin on a common base level (the
    highest outlet elevation) before comparing chi; ``none``: no control.
min_outlet_elevation : float, default=-inf
    Lowest outlet elevation (m) kept by ``outlet_control="elevation"``.
m_over_n : float, default=0.45
    Concavity index for chi (positive).
reference_area : float, default=1.0
    Chi reference area A0 (m^2, positive).
uncertainty : {"std_err", "std_dev", "bootstrap", "ttest"}, default="std_err"
    How the uncertainty of each across-divide difference is reported.
bootstrap_replicates : int, default=1000
    Replicates for ``uncertainty="bootstrap"``.
seed : int, default=5489
    Random seed for the bootstrap.

Returns
-------
DivideMigrationResult
    Channel network, chi, basin labels, and the outlet, unit, divide,
    channel-head and pair-summary tables.

Authors
-------
B. Gailleton

```python
def divide_migration(flow: lsdtt3.FlowInfo, dem: lsdtt3.Raster, outlets: object, clusters: object | None = None, source_pixels: int = 1000, snapping: str = 'downstream', outlet_control: str = 'complete_only', min_outlet_elevation: float = '-inf', m_over_n: float = 0.45, reference_area: float = 1.0, uncertainty: str = 'std_err', bootstrap_replicates: int = 1000, seed: int = 5489) -> lsdtt3.DivideMigrationResult
```

---

## `extract_channel_heads`

Extract D8 channel heads by drainage threshold or Wiener curvature.

``method='threshold'`` returns the upstream ends of all pixels with at least
``threshold_contributing_pixels`` contributing pixels (as
``FlowInfo.sources``). ``method='wiener'`` Wiener-filters the DEM, computes
tangential curvature, keeps pixels whose positive curvature exceeds a
quantile-quantile threshold, prunes by drainage area and component size,
thins the mask to a skeleton and uses its end points as heads; it requires
FFTW and a projected-metre DEM. With ``diagnostics=True``, the result
retains the filtered DEM, tangential curvature, candidate channel mask,
skeleton, and QQ table.

Parameters
----------
dem : Raster
    Elevation (m) on the flow grid (same grid and validity mask).
flow : FlowInfo
    Flow topology built from the (conditioned) DEM.
method : {"threshold", "wiener"}, default="threshold"
    Head-detection method.
threshold_contributing_pixels : int, default=1000
    ``"threshold"`` only: channel threshold in contributing pixels
    (including self); positive.
pruning_drainage_area_m2 : float, default=1000.0
    ``"wiener"`` only: candidate channel pixels with a smaller D8 drainage
    area (m^2) are removed.
surface_fitting_radius_metres : float, default=6.0
    ``"wiener"`` only: polyfit length scale (m) used for curvature, as
    ``lambda_metres`` in ``polyfit_metrics``.
connected_components_threshold_pixels : int, default=100
    ``"wiener"`` only: connected channel patches with fewer pixels (D8
    connectivity) are removed.
diagnostics : bool, default=False
    ``"wiener"`` only: keep the intermediate rasters and QQ table.

Returns
-------
ChannelHeads
    Head pixels and (optionally) diagnostics.

Authors
-------
B. Gailleton
S. M. Mudd

**Deprecated compatibility names:**
- `surface_fitting_radius_m` -> `surface_fitting_radius_metres`; earliest removal: following minor.
- `threshold_pixels` -> `threshold_contributing_pixels`; earliest removal: following minor.

```python
def extract_channel_heads(dem: lsdtt3.Raster, flow: lsdtt3.FlowInfo, method: str = 'threshold', threshold_contributing_pixels: int = 1000, pruning_drainage_area_m2: float = 1000.0, surface_fitting_radius_metres: float = 6.0, connected_components_threshold_pixels: int = 100, diagnostics: bool = False) -> lsdtt3.ChannelHeads
```

---

## `find_channel_candidates`

Find multiscale channel-bottom candidates and return score and diagnostics.

Curvature is computed at ``n_scales`` polyfit scales between
``lambda_min_metres`` and ``lambda_max_metres``; pixels whose curvature
persists across adjacent scales, that are lower than their surroundings
and (optionally) line-like are scored as channel bottoms.

Parameters
----------
raster : Raster
    Elevation raster (m), typically a high-resolution lidar DEM.
lambda_min_metres, lambda_max_metres : float, default=3.0, 20.0
    Smallest and largest polyfit length scale (m).
n_scales : int, default=5
    Number of scales between them.
curvature_metric : str, default="plan_curvature"
    ``"plan_curvature"``, ``"tangential_curvature"``,
    ``"profile_curvature"``, ``"mean_curvature"``, ``"laplacian"``,
    ``"maximum_curvature"`` or ``"minimum_curvature"``.
positive_curvature_is_good : bool, default=True
    Treat positive curvature as channel-like (False inverts the sign).
depth_window_radius_metres : float, default=12.0
    Radius (m) of the window whose mean elevation minus the pixel's gives
    the local depth signal.
depth_window_shape : {"circular", "square"}, default="circular"
    Shape of that window.
weight_persistent, weight_depth, weight_slope : float, default=0.7, 0.3, 0.0
    Weights of the scale-persistence, depth and slope signals in the score.
score_threshold : float, default=0.5
    Score at or above which a pixel enters ``mask``.
use_hessian_filter : bool, default=True
    Also require line-like (elongated) curvature.
hessian_elongation_threshold : float, default=0.6
    Minimum Hessian elongation (0-1) when the filter is on.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
ChannelCandidatesResult

Authors
-------
B. Gailleton

```python
def find_channel_candidates(raster: lsdtt3.Raster, lambda_min_metres: float = 3.0, lambda_max_metres: float = 20.0, n_scales: int = 5, curvature_metric: str = 'plan_curvature', positive_curvature_is_good: bool = True, depth_window_radius_metres: float = 12.0, depth_window_shape: str = 'circular', weight_persistent: float = 0.7, weight_depth: float = 0.3, weight_slope: float = 0.0, score_threshold: float = 0.5, use_hessian_filter: bool = True, hessian_elongation_threshold: float = 0.6, n_threads: int = 0) -> lsdtt3.ChannelCandidatesResult
```

---

## `junction_basins`

Select junction-aware basins from a channel network.

Parameters
----------
network : ChannelNetwork
    Channel network whose junctions define the basins.
selection : {"all", "outlets", "nested"}, default="all"
    ``"all"``: one basin per junction, draining to the pixel just upstream
    of the next junction downstream. ``"outlets"``: one basin per supplied
    outlet. ``"nested"``: each supplied outlet's basin plus every junction
    basin nested inside it.
outlets : array-like of int or None, default=None
    N x 2 array of ``(row, col)`` outlet pixels (valid pixels); required
    for ``"outlets"`` and ``"nested"``, ignored for ``"all"``.

Returns
-------
JunctionBasinCollection
    Basins that can be filtered by area, stream order and nesting.

Authors
-------
B. Gailleton

```python
def junction_basins(network: lsdtt3.ChannelNetwork, selection: str = 'all', outlets: object | None = None) -> lsdtt3.JunctionBasinCollection
```

---

## `line_crossing_basins`

Create basins whose outlets are channel pixels crossed by a line.

Every network channel pixel within 0.71 cell sizes of the line, with at least
``minimum_contributing_pixels`` contributing pixels, becomes an outlet.

Parameters
----------
network : ChannelNetwork
    Channel network on the line's grid.
line : shapely geometry or geopandas.GeoDataFrame
    Line (for example a mountain front or a road) in the DEM CRS; a
    GeoDataFrame is merged into one geometry.
minimum_contributing_pixels : int, default=500
    Smallest channel crossing kept, in contributing pixels (positive).

Returns
-------
BasinCollection

Authors
-------
B. Gailleton

```python
def line_crossing_basins(network: lsdtt3.ChannelNetwork, line: object, minimum_contributing_pixels: int = 500) -> lsdtt3.BasinCollection
```

---

## `side_tributary_basins`

Extract the side basins that drain into channels traced from input points.

Parameters
----------
flow : FlowInfo
    Flow topology.
points_path : str
    Vector file of the points from which trunk channels are traced D8
    downstream.
minimum_contributing_pixels : int, default=1000
    Smallest side basin kept, in contributing pixels (positive).
reject_incomplete : bool, default=True
    Drop side basins that touch nodata or an open DEM edge.
single_channel : bool, default=False
    Trace a single main channel: the file holds one point (traced to the
    DEM edge) or two points (traced from the upstream one to where the
    other's flow path joins it). If False every point is a source.
include_stopping_point : bool, default=False
    Treat the last point in the file as a stopping point: each trace ends
    when it first enters a square window around it.
stopping_point_window_pixels : int, default=5
    Half-width (pixels) of the stopping-point window.

Returns
-------
SideTributaryBasins

Authors
-------
B. Gailleton

```python
def side_tributary_basins(flow: lsdtt3.FlowInfo, points_path: str, minimum_contributing_pixels: int = 1000, reject_incomplete: bool = True, single_channel: bool = False, include_stopping_point: bool = False, stopping_point_window_pixels: int = 5) -> lsdtt3.SideTributaryBasins
```
