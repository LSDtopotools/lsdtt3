# Flow, Channels, And Basins API

This page is generated from the built Python runtime manifest. Do not edit signatures here.

## `BasinCollection`

Selected nested D8 drainage basins retaining their parent ``FlowInfo`` topology.

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
    Selected nested D8 drainage basins retaining their parent ``FlowInfo`` topology.

Each basin is defined by an outlet pixel and stores an efficient native
upslope slice; no raster-sized copy is retained per basin.

    Authors
    -------
    B. Gailleton


### Members

#### `labels`

Basin IDs (one-based); nested basins overwrite parents.

#### `mask`

mask(self, basin_id: int) -> lsdtt3.Raster

```python
def mask(self, basin_id: int) -> lsdtt3.Raster
```

#### `membership`

membership(self, row: int, col: int) -> object

Return all selected basin IDs containing this pixel.

```python
def membership(self, row: int, col: int) -> object
```

#### `outline_ridgelines`

outline_ridgelines(self, dem: lsdtt3.Raster, split_method: str = 'longest_channel') -> object

Split every selected basin outline into two ridgelines.

Returns a GeoDataFrame of left/right ridge lines. ``split_method`` is
``"longest_channel"``, ``"maximum_elevation"``, or ``"equal_length"``.

```python
def outline_ridgelines(self, dem: lsdtt3.Raster, split_method: str = 'longest_channel') -> object
```

#### `pixels`

pixels(self, basin_id: int) -> object

```python
def pixels(self, basin_id: int) -> object
```

#### `polygons`

Basin outline polygons as a GeoDataFrame.

#### `size`

(self) -> int

#### `statistics`

statistics(self, dem: lsdtt3.Raster, include_slope: bool = True, include_flow_distance: bool = True, include_shape_metrics: bool = False, perimeter_simplify_cells: float = 1.0) -> object

Return per-basin elevation, slope, flow-distance, and optional shape statistics.

```python
def statistics(self, dem: lsdtt3.Raster, include_slope: bool = True, include_flow_distance: bool = True, include_shape_metrics: bool = False, perimeter_simplify_cells: float = 1.0) -> object
```

#### `table`

(self) -> object

---

## `ChannelCandidatesResult`

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
    Channel-head locations resolved from a D8 flow topology.

Use ``pixels`` as explicit ``FlowInfo.channel_scope(..., sources=...)`` input.

    Authors
    -------
    B. Gailleton


### Members

#### `channel_mask`

(self) -> object

#### `filtered_dem`

(self) -> object

#### `method`

(self) -> str

#### `n_heads`

(self) -> int

#### `pixels`

(self) -> object

#### `qq_table`

(self) -> object

#### `skeleton`

(self) -> object

#### `table`

(self) -> object

#### `tangential_curvature`

(self) -> object

---

## `ChannelNetwork`

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

#### `basin_com_bearings`

basin_com_bearings(self) -> object

Return outlet-to-basin-centre-of-mass bearings for each junction basin.

```python
def basin_com_bearings(self) -> object
```

#### `channel_mask`

(self) -> lsdtt3.Raster

#### `chi`

chi(self, outlets: object | None = None, m_over_n: float = 0.45, reference_area: float = 1.0, denominator: object | None = None) -> lsdtt3.ChiCollection

Extract outlet-seeded chi coordinates on this channel network.

```python
def chi(self, outlets: object | None = None, m_over_n: float = 0.45, reference_area: float = 1.0, denominator: object | None = None) -> lsdtt3.ChiCollection
```

#### `chi_sweep`

chi_sweep(self, dem: lsdtt3.Raster, m_over_n_values: object, outlets: object | None = None, reference_area: float = 1.0, critical_divergence_metres: float = 5.0, denominator: object | None = None) -> object

Fit each tested m/n and return per-basin disorder.

```python
def chi_sweep(self, dem: lsdtt3.Raster, m_over_n_values: object, outlets: object | None = None, reference_area: float = 1.0, critical_divergence_metres: float = 5.0, denominator: object | None = None) -> object
```

#### `donors`

donors(self, junction_id: int) -> list[int]

```python
def donors(self, junction_id: int) -> list[int]
```

#### `junction_angles`

junction_angles(self, junctions: object | None = None, mode: str = 'full_link', fixed_length_metres: float = 0.0, fixed_n_pixels: int = 0, min_nodes: int = 2) -> object

Measure two-donor confluence angles in compass degrees.

```python
def junction_angles(self, junctions: object | None = None, mode: str = 'full_link', fixed_length_metres: float = 0.0, fixed_n_pixels: int = 0, min_nodes: int = 2) -> object
```

#### `junction_id`

(self) -> object

#### `junctions`

junctions(self) -> object

Return junction topology as a pandas DataFrame.

```python
def junctions(self) -> object
```

#### `junctions_geodataframe`

junctions_geodataframe(self) -> object

Return channel junctions as a GeoPandas point layer.

```python
def junctions_geodataframe(self) -> object
```

#### `lines_geodataframe`

lines_geodataframe(self) -> object

Return compact non-outlet channel links as a GeoPandas line layer.

```python
def lines_geodataframe(self) -> object
```

#### `link_bearings`

link_bearings(self, junctions: object | None = None, mode: str = 'full_link', fixed_length_metres: float = 0.0, fixed_n_pixels: int = 0, anchor: str = 'upstream', min_nodes: int = 2) -> object

Fit compass bearings to compact channel links.

```python
def link_bearings(self, junctions: object | None = None, mode: str = 'full_link', fixed_length_metres: float = 0.0, fixed_n_pixels: int = 0, anchor: str = 'upstream', min_nodes: int = 2) -> object
```

#### `link_profiles`

link_profiles(self, values: lsdtt3.Raster, junctions: object | None = None) -> object

Sample a raster along compact downstream channel links.

```python
def link_profiles(self, values: lsdtt3.Raster, junctions: object | None = None) -> object
```

#### `links`

links(self) -> object

Return junction-to-junction links as a pandas DataFrame.

```python
def links(self) -> object
```

#### `n_channel_cells`

(self) -> int

#### `n_junctions`

(self) -> int

#### `points_geodataframe`

points_geodataframe(self) -> object

Return channel pixels as a GeoPandas point layer.

```python
def points_geodataframe(self) -> object
```

---

## `ChannelScope`

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

#### `mask`

(self) -> lsdtt3.Raster

#### `mode`

(self) -> str

#### `nodes`

nodes(self, source_id: int) -> object

```python
def nodes(self, source_id: int) -> object
```

#### `outlets`

(self) -> object

#### `source_id`

(self) -> object

#### `sources`

(self) -> object

#### `table`

table(self) -> object

Return ordered channel nodes as a pandas DataFrame.

```python
def table(self) -> object
```

---

## `ConnectivityResult`

Sediment-connectivity rasters following Cavalli et al. (2013).

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
Sediment-connectivity rasters following Cavalli et al. (2013).

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

#### `basin_labels`

(self) -> lsdtt3.Raster

#### `channel_heads`

(self) -> object

#### `chi`

(self) -> lsdtt3.ChiCollection

#### `divides`

(self) -> object

#### `network`

(self) -> lsdtt3.ChannelNetwork

#### `outlets`

(self) -> object

#### `pair_summary`

(self) -> object

#### `units`

(self) -> object

---

## `FlowGraph`

Fractional D-infinity or MFD flow-routing graph.

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
    Fractional D-infinity or MFD flow-routing graph.

Use ``FlowInfo`` for single-receiver D8 channel trees. Use ``FlowGraph`` for
distributed-flow drainage area and later TWI/connectivity calculations.

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

(self) -> object

#### `drainage_area`

drainage_area(self) -> lsdtt3.Raster

```python
def drainage_area(self) -> lsdtt3.Raster
```

#### `is_outlet`

is_outlet(self, row: int, col: int) -> bool

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

```python
def outlet_mask(self) -> lsdtt3.Raster
```

#### `outlets`

(self) -> object

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

(self) -> object

---

## `FlowInfo`

Immutable D8 drainage topology derived from a conditioned raster.

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
    Immutable D8 drainage topology derived from a conditioned raster.

``FlowInfo`` owns routing topology, not elevations. Prepare depressions before
construction; later operations explicitly receive rasters where needed.

    Authors
    -------
    B. Gailleton


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

Accumulate weights downstream; default weights are cell areas.

```python
def accumulate(self, weights: object | None = None) -> lsdtt3.Raster
```

#### `basins`

basins(self, outlets: object | None = None, min_upslope_pixels: int = 0) -> lsdtt3.BasinCollection

Create a ``BasinCollection`` from explicit outlet pixels or all qualifying D8 outlets.

Explicit outlets may be internal channel pixels, which creates nested basins.

```python
def basins(self, outlets: object | None = None, min_upslope_pixels: int = 0) -> lsdtt3.BasinCollection
```

#### `boundary_conditions`

(self) -> str

#### `channel_mask`

channel_mask(self, threshold_pixels: int) -> lsdtt3.Raster

Return a 0/1 raster where contributing pixels meet the threshold.

```python
def channel_mask(self, threshold_pixels: int) -> lsdtt3.Raster
```

#### `channel_network`

channel_network(self, scope: lsdtt3.ChannelScope) -> lsdtt3.ChannelNetwork

Build a compact junction graph from a ChannelScope.

```python
def channel_network(self, scope: lsdtt3.ChannelScope) -> lsdtt3.ChannelNetwork
```

#### `channel_scope`

channel_scope(self, threshold_pixels: int, sources: object | None = None, outlets: object | None = None, snap_sources: bool = True, source_snap_max_radius_pixels: int = 25, outlet_snap_max_distance_metres: float = 0.0, auto_detect_roles: bool = False) -> lsdtt3.ChannelScope

Resolve threshold or supplied source paths into a ChannelScope.

```python
def channel_scope(self, threshold_pixels: int, sources: object | None = None, outlets: object | None = None, snap_sources: bool = True, source_snap_max_radius_pixels: int = 25, outlet_snap_max_distance_metres: float = 0.0, auto_detect_roles: bool = False) -> lsdtt3.ChannelScope
```

#### `chi_from_outlets`

chi_from_outlets(self, outlets: object | None = None, m_over_n: float = 0.45, reference_area: float = 1.0, denominator: object | None = None) -> lsdtt3.Raster

Integrate zero-based chi across complete upstream raster domains.

```python
def chi_from_outlets(self, outlets: object | None = None, m_over_n: float = 0.45, reference_area: float = 1.0, denominator: object | None = None) -> lsdtt3.Raster
```

#### `chi_from_outlets_inherited`

chi_from_outlets_inherited(self, outlets: object, seed_chi: lsdtt3.Raster, m_over_n: float = 0.45, reference_area: float = 1.0, denominator: object | None = None) -> lsdtt3.Raster

Integrate chi from supplied outlets using seed-raster base values.

```python
def chi_from_outlets_inherited(self, outlets: object, seed_chi: lsdtt3.Raster, m_over_n: float = 0.45, reference_area: float = 1.0, denominator: object | None = None) -> lsdtt3.Raster
```

#### `contributing_pixels`

contributing_pixels(self) -> lsdtt3.Raster

```python
def contributing_pixels(self) -> lsdtt3.Raster
```

#### `direction`

(self) -> object

#### `distance_to_channel`

distance_to_channel(self, threshold_pixels: int) -> lsdtt3.Raster

```python
def distance_to_channel(self, threshold_pixels: int) -> lsdtt3.Raster
```

#### `distance_to_farthest_ridge`

distance_to_farthest_ridge(self) -> lsdtt3.Raster

```python
def distance_to_farthest_ridge(self) -> lsdtt3.Raster
```

#### `donor_count`

(self) -> object

#### `downstream_sum`

downstream_sum(self, values: lsdtt3.Raster, sinks: lsdtt3.Raster) -> lsdtt3.Raster

```python
def downstream_sum(self, values: lsdtt3.Raster, sinks: lsdtt3.Raster) -> lsdtt3.Raster
```

#### `drainage_area`

drainage_area(self) -> lsdtt3.Raster

```python
def drainage_area(self) -> lsdtt3.Raster
```

#### `flat_index`

flat_index(self, row: int, col: int) -> int

```python
def flat_index(self, row: int, col: int) -> int
```

#### `flow_direction`

flow_direction(self) -> lsdtt3.Raster

```python
def flow_direction(self) -> lsdtt3.Raster
```

#### `flow_distance_to_outlet`

flow_distance_to_outlet(self) -> lsdtt3.Raster

```python
def flow_distance_to_outlet(self) -> lsdtt3.Raster
```

#### `flow_lines_geodataframe`

flow_lines_geodataframe(self) -> object

Return all D8 receiver edges as a GeoPandas line layer.

```python
def flow_lines_geodataframe(self) -> object
```

#### `height_above_channel`

height_above_channel(self, dem: lsdtt3.Raster, channels: lsdtt3.Raster) -> lsdtt3.Raster

Compute HAND from an explicit 0/1 channel raster.

```python
def height_above_channel(self, dem: lsdtt3.Raster, channels: lsdtt3.Raster) -> lsdtt3.Raster
```

#### `height_above_channel_with_donor`

height_above_channel_with_donor(self, dem: lsdtt3.Raster, channels: lsdtt3.Raster) -> dict

```python
def height_above_channel_with_donor(self, dem: lsdtt3.Raster, channels: lsdtt3.Raster) -> dict
```

#### `height_to_farthest_ridge`

height_to_farthest_ridge(self, dem: lsdtt3.Raster) -> lsdtt3.Raster

```python
def height_to_farthest_ridge(self, dem: lsdtt3.Raster) -> lsdtt3.Raster
```

#### `is_outlet`

is_outlet(self, row: int, col: int) -> bool

```python
def is_outlet(self, row: int, col: int) -> bool
```

#### `is_upslope`

is_upslope(self, candidate_row: int, candidate_col: int, target_row: int, target_col: int) -> bool

```python
def is_upslope(self, candidate_row: int, candidate_col: int, target_row: int, target_col: int) -> bool
```

#### `metadata`

(self) -> lsdtt3.RasterMeta

#### `n_outlets`

(self) -> int

#### `n_valid_cells`

(self) -> int

#### `outlets`

(self) -> object

#### `receiver`

(self) -> object

#### `receiver_of`

receiver_of(self, row: int, col: int) -> object

```python
def receiver_of(self, row: int, col: int) -> object
```

#### `row_col`

row_col(self, flat_index: int) -> tuple

```python
def row_col(self, flat_index: int) -> tuple
```

#### `scope_from_heads`

scope_from_heads(self, heads: lsdtt3.ChannelHeads, outlets: object | None = None, outlet_snap_max_distance_metres: float = 0.0) -> lsdtt3.ChannelScope

Resolve ChannelHeads into explicit D8 source paths for channel and chi analysis.

```python
def scope_from_heads(self, heads: lsdtt3.ChannelHeads, outlets: object | None = None, outlet_snap_max_distance_metres: float = 0.0) -> lsdtt3.ChannelScope
```

#### `shape`

(self) -> tuple

#### `sources`

sources(self, threshold_pixels: int) -> object

Return channel sources as an N by 2 (row, col) array.

```python
def sources(self, threshold_pixels: int) -> object
```

#### `upslope_mean`

upslope_mean(self, values: lsdtt3.Raster, include_self: bool = True) -> lsdtt3.Raster

```python
def upslope_mean(self, values: lsdtt3.Raster, include_self: bool = True) -> lsdtt3.Raster
```

#### `upslope_pixels`

upslope_pixels(self, row: int, col: int) -> object

```python
def upslope_pixels(self, row: int, col: int) -> object
```

#### `upslope_sum`

upslope_sum(self, values: lsdtt3.Raster, include_self: bool = True) -> lsdtt3.Raster

```python
def upslope_sum(self, values: lsdtt3.Raster, include_self: bool = True) -> lsdtt3.Raster
```

#### `valid_mask`

(self) -> object

---

## `JunctionBasinCollection`

Junction-aware D8 basins selected from a ``ChannelNetwork``.

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
    Junction-aware D8 basins selected from a ``ChannelNetwork``.

The collection retains the channel-network stream order and supports in-place
selection filters. It is created with :func:`junction_basins`.

    Authors
    -------
    B. Gailleton


### Members

#### `filter_area`

filter_area(self, min_pixels: int = -1, max_pixels: int = -1) -> lsdtt3.JunctionBasinCollection

```python
def filter_area(self, min_pixels: int = -1, max_pixels: int = -1) -> lsdtt3.JunctionBasinCollection
```

#### `filter_incomplete`

filter_incomplete(self) -> lsdtt3.JunctionBasinCollection

```python
def filter_incomplete(self) -> lsdtt3.JunctionBasinCollection
```

#### `filter_largest`

filter_largest(self) -> lsdtt3.JunctionBasinCollection

```python
def filter_largest(self) -> lsdtt3.JunctionBasinCollection
```

#### `filter_nested`

filter_nested(self) -> lsdtt3.JunctionBasinCollection

```python
def filter_nested(self) -> lsdtt3.JunctionBasinCollection
```

#### `filter_nested_by_stream_order`

filter_nested_by_stream_order(self) -> lsdtt3.JunctionBasinCollection

```python
def filter_nested_by_stream_order(self) -> lsdtt3.JunctionBasinCollection
```

#### `filter_stream_order`

filter_stream_order(self, min_order: int = -1, max_order: int = -1) -> lsdtt3.JunctionBasinCollection

```python
def filter_stream_order(self, min_order: int = -1, max_order: int = -1) -> lsdtt3.JunctionBasinCollection
```

#### `labels`

(self) -> lsdtt3.Raster

#### `mask`

mask(self, basin_id: int) -> lsdtt3.Raster

```python
def mask(self, basin_id: int) -> lsdtt3.Raster
```

#### `size`

(self) -> int

#### `table`

(self) -> object

---

## `SideTributaryBasins`

Side basins draining into a traced D8 channel.

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
Side basins draining into a traced D8 channel.

Authors
-------
B. Gailleton

### Members

#### `basins`

(self) -> lsdtt3.BasinCollection

#### `channel_mask`

(self) -> lsdtt3.Raster

#### `table`

(self) -> object

---

## `basin_channel_lengths`

basin_channel_lengths(basins: lsdtt3.BasinCollection, network: lsdtt3.ChannelNetwork, max_stream_order: int = 6) -> object

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
    basin_channel_lengths(basins: lsdtt3.BasinCollection, network: lsdtt3.ChannelNetwork, max_stream_order: int = 6) -> object

Return total in-basin D8 channel length by Strahler stream order.

``basins`` and ``network`` must have been created from the same ``FlowInfo``.

    Authors
    -------
    B. Gailleton


```python
def basin_channel_lengths(basins: lsdtt3.BasinCollection, network: lsdtt3.ChannelNetwork, max_stream_order: int = 6) -> object
```

---

## `basin_surface_metrics`

basin_surface_metrics(basins: lsdtt3.BasinCollection, dem: lsdtt3.Raster, scales: collections.abc.Sequence[float], metrics: collections.abc.Sequence[str], include_outline: bool = True, n_threads: int = 0) -> dict

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
basin_surface_metrics(basins: lsdtt3.BasinCollection, dem: lsdtt3.Raster, scales: collections.abc.Sequence[float], metrics: collections.abc.Sequence[str], include_outline: bool = True, n_threads: int = 0) -> dict

Authors
-------
B. Gailleton

```python
def basin_surface_metrics(basins: lsdtt3.BasinCollection, dem: lsdtt3.Raster, scales: collections.abc.Sequence[float], metrics: collections.abc.Sequence[str], include_outline: bool = True, n_threads: int = 0) -> dict
```

---

## `divide`

divide(a: lsdtt3.Raster, b: lsdtt3.Raster) -> lsdtt3.Raster

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
    divide(a: lsdtt3.Raster, b: lsdtt3.Raster) -> lsdtt3.Raster

Divide two grid-compatible rasters.

    Authors
    -------
    B. Gailleton


```python
def divide(a: lsdtt3.Raster, b: lsdtt3.Raster) -> lsdtt3.Raster
```

---

## `divide_migration`

divide_migration(flow: lsdtt3.FlowInfo, dem: lsdtt3.Raster, outlets: object, clusters: object | None = None, source_pixels: int = 1000, snapping: str = 'downstream', outlet_control: str = 'complete_only', min_outlet_elevation: float = -inf, m_over_n: float = 0.45, reference_area: float = 1.0, uncertainty: str = 'std_err', bootstrap_replicates: int = 1000, seed: int = 5489) -> lsdtt3.DivideMigrationResult

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
    divide_migration(flow: lsdtt3.FlowInfo, dem: lsdtt3.Raster, outlets: object, clusters: object | None = None, source_pixels: int = 1000, snapping: str = 'downstream', outlet_control: str = 'complete_only', min_outlet_elevation: float = -inf, m_over_n: float = 0.45, reference_area: float = 1.0, uncertainty: str = 'std_err', bootstrap_replicates: int = 1000, seed: int = 5489) -> lsdtt3.DivideMigrationResult

Compare channel-head elevation, chi, relief, and gradient across adjacent outlet-defined basin divides.

    Authors
    -------
    B. Gailleton


```python
def divide_migration(flow: lsdtt3.FlowInfo, dem: lsdtt3.Raster, outlets: object, clusters: object | None = None, source_pixels: int = 1000, snapping: str = 'downstream', outlet_control: str = 'complete_only', min_outlet_elevation: float = '-inf', m_over_n: float = 0.45, reference_area: float = 1.0, uncertainty: str = 'std_err', bootstrap_replicates: int = 1000, seed: int = 5489) -> lsdtt3.DivideMigrationResult
```

---

## `extract_channel_heads`

extract_channel_heads(dem: lsdtt3.Raster, flow: lsdtt3.FlowInfo, method: str = 'threshold', threshold_pixels: int = 1000, pruning_drainage_area_m2: float = 1000.0, surface_fitting_radius_m: float = 6.0, connected_components_threshold_pixels: int = 100, diagnostics: bool = False) -> lsdtt3.ChannelHeads

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
    extract_channel_heads(dem: lsdtt3.Raster, flow: lsdtt3.FlowInfo, method: str = 'threshold', threshold_pixels: int = 1000, pruning_drainage_area_m2: float = 1000.0, surface_fitting_radius_m: float = 6.0, connected_components_threshold_pixels: int = 100, diagnostics: bool = False) -> lsdtt3.ChannelHeads

Extract D8 channel heads by drainage threshold or Wiener curvature.

``method='wiener'`` requires FFTW and a projected-metre DEM. With
``diagnostics=True``, the result retains the filtered DEM, tangential
curvature, candidate channel mask, skeleton, and QQ table.

    Authors
    -------
    B. Gailleton


```python
def extract_channel_heads(dem: lsdtt3.Raster, flow: lsdtt3.FlowInfo, method: str = 'threshold', threshold_pixels: int = 1000, pruning_drainage_area_m2: float = 1000.0, surface_fitting_radius_m: float = 6.0, connected_components_threshold_pixels: int = 100, diagnostics: bool = False) -> lsdtt3.ChannelHeads
```

---

## `find_channel_candidates`

find_channel_candidates(raster: lsdtt3.Raster, lambda_min_metres: float = 3.0, lambda_max_metres: float = 20.0, n_scales: int = 5, curvature_metric: str = 'plan_curvature', positive_curvature_is_good: bool = True, depth_window_radius_metres: float = 12.0, depth_window_shape: str = 'circular', weight_persistent: float = 0.7, weight_depth: float = 0.3, weight_slope: float = 0.0, score_threshold: float = 0.5, use_hessian_filter: bool = True, hessian_elongation_threshold: float = 0.6, n_threads: int = 0) -> lsdtt3.ChannelCandidatesResult

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
    find_channel_candidates(raster: lsdtt3.Raster, lambda_min_metres: float = 3.0, lambda_max_metres: float = 20.0, n_scales: int = 5, curvature_metric: str = 'plan_curvature', positive_curvature_is_good: bool = True, depth_window_radius_metres: float = 12.0, depth_window_shape: str = 'circular', weight_persistent: float = 0.7, weight_depth: float = 0.3, weight_slope: float = 0.0, score_threshold: float = 0.5, use_hessian_filter: bool = True, hessian_elongation_threshold: float = 0.6, n_threads: int = 0) -> lsdtt3.ChannelCandidatesResult

Find multiscale channel-bottom candidates and return score and diagnostics.

    Authors
    -------
    B. Gailleton


```python
def find_channel_candidates(raster: lsdtt3.Raster, lambda_min_metres: float = 3.0, lambda_max_metres: float = 20.0, n_scales: int = 5, curvature_metric: str = 'plan_curvature', positive_curvature_is_good: bool = True, depth_window_radius_metres: float = 12.0, depth_window_shape: str = 'circular', weight_persistent: float = 0.7, weight_depth: float = 0.3, weight_slope: float = 0.0, score_threshold: float = 0.5, use_hessian_filter: bool = True, hessian_elongation_threshold: float = 0.6, n_threads: int = 0) -> lsdtt3.ChannelCandidatesResult
```

---

## `junction_basins`

junction_basins(network: lsdtt3.ChannelNetwork, selection: str = 'all', outlets: object | None = None) -> lsdtt3.JunctionBasinCollection

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
    junction_basins(network: lsdtt3.ChannelNetwork, selection: str = 'all', outlets: object | None = None) -> lsdtt3.JunctionBasinCollection

Select junction-aware basins from a channel network.

``selection`` is ``"all"``, ``"outlets"``, or ``"nested"``. The latter two
require outlet pixels as an N-by-2 ``(row, col)`` array.

    Authors
    -------
    B. Gailleton


```python
def junction_basins(network: lsdtt3.ChannelNetwork, selection: str = 'all', outlets: object | None = None) -> lsdtt3.JunctionBasinCollection
```

---

## `line_crossing_basins`

line_crossing_basins(network: lsdtt3.ChannelNetwork, line: object, minimum_contributing_pixels: int = 500) -> lsdtt3.BasinCollection

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
line_crossing_basins(network: lsdtt3.ChannelNetwork, line: object, minimum_contributing_pixels: int = 500) -> lsdtt3.BasinCollection

Authors
-------
B. Gailleton

```python
def line_crossing_basins(network: lsdtt3.ChannelNetwork, line: object, minimum_contributing_pixels: int = 500) -> lsdtt3.BasinCollection
```

---

## `side_tributary_basins`

side_tributary_basins(flow: lsdtt3.FlowInfo, points_path: str, minimum_contributing_pixels: int = 1000, reject_incomplete: bool = True, single_channel: bool = False, include_stopping_point: bool = False, stopping_point_window_pixels: int = 5) -> lsdtt3.SideTributaryBasins

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
side_tributary_basins(flow: lsdtt3.FlowInfo, points_path: str, minimum_contributing_pixels: int = 1000, reject_incomplete: bool = True, single_channel: bool = False, include_stopping_point: bool = False, stopping_point_window_pixels: int = 5) -> lsdtt3.SideTributaryBasins

Authors
-------
B. Gailleton

```python
def side_tributary_basins(flow: lsdtt3.FlowInfo, points_path: str, minimum_contributing_pixels: int = 1000, reject_incomplete: bool = True, single_channel: bool = False, include_stopping_point: bool = False, stopping_point_window_pixels: int = 5) -> lsdtt3.SideTributaryBasins
```
