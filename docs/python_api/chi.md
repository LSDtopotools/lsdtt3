# Chi And Ksn API

This page is generated from the built Python runtime manifest. Do not edit signatures here.

## `ChiCollection`

Chi coordinates of the channel pixels of one or more basins.

Made by ``ChannelNetwork.chi()``. Basin ``basin_id`` (0-based) is the
``basin_id``-th outlet. Chi (m) is 0 at each outlet and increases upstream.
``nodes()`` gives the values as a table, ``raster()`` as a Raster, and
``fit(dem)`` fits piecewise-linear segments in chi-elevation space to give
channel steepness (k_sn).

Returned by ``ChannelNetwork.chi()``, ``DivideMigrationResult.chi``; not constructible directly.

Attributes
----------
n_basins : int
    Number of basins (outlets).

Methods
-------
fit(dem: Raster, critical_divergence_metres: float = 5.0, mainstem: str = 'max_flow_length') -> ChiFit
    Fit piecewise-linear channel profiles in chi-elevation space.
nodes() -> object
    Return extracted channel chi nodes as a pandas DataFrame.
raster() -> Raster
    Return chi values as a masked Raster.

Authors
-------
B. Gailleton

### Members

#### `fit`

fit(self, dem: lsdtt3.Raster, critical_divergence_metres: float = 5.0, mainstem: str = 'max_flow_length') -> lsdtt3.ChiFit

Fit piecewise-linear channel profiles in chi-elevation space.

Each basin's main stem and then each tributary is fitted with straight
segments between its ends; a segment is split at its worst-fitting pixel
while the elevation misfit there exceeds ``critical_divergence_metres``.
Each segment's slope dz/dchi is its k_sn.

Parameters
----------
dem : Raster
    Elevation (m) on the network grid, normally the conditioned DEM.
critical_divergence_metres : float, default=5.0
    Largest allowed vertical misfit (m) between the observed profile and a
    segment before the segment is split (non-negative). Larger values give
    fewer, longer segments; noisier DEMs need larger values (for example
    10 m on a 30 m DEM).
mainstem : {"max_flow_length", "max_chi"}, default="max_flow_length"
    How the main stem of each basin is chosen: the longest D8 flow path,
    or the path reaching the largest chi.

Returns
-------
ChiFit
    Fitted nodes, segments and per-basin k_sn summaries.

```python
def fit(self, dem: lsdtt3.Raster, critical_divergence_metres: float = 5.0, mainstem: str = 'max_flow_length') -> lsdtt3.ChiFit
```

#### `n_basins`

Number of basins (outlets).

#### `nodes`

nodes(self) -> object

Return extracted channel chi nodes as a pandas DataFrame.

Returns
-------
pandas.DataFrame
    One row per channel pixel: ``basin_id`` (0-based), ``outlet_flat``,
    ``flat``, ``receiver_flat`` (downstream pixel, negative at the outlet),
    ``row``, ``col``, ``chi`` (m).

```python
def nodes(self) -> object
```

#### `raster`

raster(self) -> lsdtt3.Raster

Return chi values as a masked Raster.

Returns
-------
Raster
    Chi (m) on the channel pixels of every basin; nodata elsewhere.

```python
def raster(self) -> lsdtt3.Raster
```

---

## `ChiFit`

Piecewise-linear fits of channel profiles in chi-elevation space.

Made by ``ChiCollection.fit(dem)``. Each segment is a straight line in
(chi, elevation); its slope ``ksn = dz/dchi`` is the channel steepness
index. With ``reference_area = 1`` m^2 this is k_sn in m^(2 theta), where
theta = m_over_n; with another A0 it is k_sn / A0**theta (dimensionless).
All three tables are pandas DataFrames; ``basin_id`` is 0-based and
``flat`` values are flat indices (row = flat // n_cols, col = flat % n_cols;
or ``FlowInfo.row_col(flat)``). lsdviztools3 chi plots expect the column
names ``basin_index`` and ``k_sn``: rename with
``df.rename(columns={"basin_id": "basin_index", "ksn": "k_sn"})``.

Returned by ``ChiCollection.fit()``; not constructible directly.

Attributes
----------
nodes : pandas.DataFrame
    One row per channel pixel: ``basin_id``, ``outlet_flat``, ``flat``,
    ``receiver_flat``, ``chi`` (m), ``observed_elevation`` (m),
    ``fitted_elevation`` (m), ``residual_elevation`` (observed - fitted, m),
    ``abs_residual_elevation`` (m), ``branch_id``, ``is_main_stem``,
    ``is_source``, ``is_junction``, ``is_split_point`` (pixel where segments
    meet). There are no row/col columns: use ``flat`` as above.
segments : pandas.DataFrame
    One row per segment: ``basin_id``, ``outlet_flat``, ``segment_id``,
    ``parent_segment_id``, ``branch_id``, ``downstream_flat``,
    ``upstream_flat``, ``downstream_chi``, ``upstream_chi`` (m),
    ``downstream_elevation``, ``upstream_elevation`` (m), ``ksn``,
    ``max_abs_residual`` (m), ``n_nodes``.
summary : pandas.DataFrame
    One row per basin: ``basin_id``, ``outlet_flat``, ``n_segments``,
    ``n_main_stem_segments``, ``total_chi_length`` (m), ``mean_ksn``,
    ``median_ksn``, ``min_ksn``, ``max_ksn``, ``chi_weighted_mean_ksn``,
    ``main_stem_mean_ksn``, ``main_stem_median_ksn``.

Methods
-------
disorder() -> object
    Return full-network chi-disorder per basin.
disorder_bootstrap(replicates: int = 1000, seed: int = 5489, sample_full_tributary_count: bool = True, n_side_branches: int = 3, with_replacement: bool = True) -> object
    Bootstrap disorder over side branches.
disorder_jackknife() -> object
    Leave-one-side-branch-out disorder samples.
disorder_subset_sensitivity(n_side_branches: int = 3) -> object
    Evaluate disorder across mainstem-plus-tributary subsets.
invert_ksn(m_over_n_values: object, reference_area: float = 1.0, coarse_cell_size_metres: float = 2000.0, roughness_weight: float = 10000.0, diagonal_jitter: float = 1e-08, use_mae: bool = True) -> ChiInversion
    Run regularised spatial ksn/u-star inversion over m/n trials.

Authors
-------
B. Gailleton

### Members

#### `disorder`

disorder(self) -> object

Return full-network chi-disorder per basin.

Disorder measures how badly tributary profiles fail to collapse onto the
main stem in chi-elevation space; lower is better.

Returns
-------
pandas.DataFrame
    Columns ``basin_id``, ``outlet_flat``, ``n_nodes``,
    ``n_main_stem_nodes``, ``n_branches``, ``n_side_branches``, ``disorder``.

```python
def disorder(self) -> object
```

#### `disorder_bootstrap`

disorder_bootstrap(self, replicates: int = 1000, seed: int = 5489, sample_full_tributary_count: bool = True, n_side_branches: int = 3, with_replacement: bool = True) -> object

Bootstrap disorder over side branches.

Parameters
----------
replicates : int, default=1000
    Number of bootstrap replicates (positive).
seed : int, default=5489
    Random seed, for reproducible samples.
sample_full_tributary_count : bool, default=True
    Draw as many side branches as the basin has; if False draw
    ``n_side_branches`` per replicate.
n_side_branches : int, default=3
    Side branches per replicate when ``sample_full_tributary_count`` is False.
with_replacement : bool, default=True
    Sample side branches with replacement.

Returns
-------
pandas.DataFrame
    Columns ``basin_id``, ``full_disorder``, ``replicate_id``, ``disorder``,
    ``n_nodes``, ``mean``, ``standard_deviation``, ``percentile_2p5``,
    ``percentile_97p5``.

```python
def disorder_bootstrap(self, replicates: int = 1000, seed: int = 5489, sample_full_tributary_count: bool = True, n_side_branches: int = 3, with_replacement: bool = True) -> object
```

#### `disorder_jackknife`

disorder_jackknife(self) -> object

Leave-one-side-branch-out disorder samples.

Returns
-------
pandas.DataFrame
    Columns ``basin_id``, ``full_disorder``, ``replicate_id``,
    ``removed_branch_id``, ``disorder``, ``n_nodes``, ``mean``,
    ``standard_error``.

```python
def disorder_jackknife(self) -> object
```

#### `disorder_subset_sensitivity`

disorder_subset_sensitivity(self, n_side_branches: int = 3) -> object

Evaluate disorder across mainstem-plus-tributary subsets.

Parameters
----------
n_side_branches : int, default=3
    Number of side branches combined with the main stem in each subset.

Returns
-------
pandas.DataFrame
    Columns ``basin_id``, ``full_disorder``, ``sample_id``, ``disorder``,
    ``n_nodes``, ``branch_ids``.

```python
def disorder_subset_sensitivity(self, n_side_branches: int = 3) -> object
```

#### `invert_ksn`

invert_ksn(self, m_over_n_values: object, reference_area: float = 1.0, coarse_cell_size_metres: float = 2000.0, roughness_weight: float = 10000.0, diagonal_jitter: float = 1e-08, use_mae: bool = True) -> lsdtt3.ChiInversion

Run regularised spatial ksn/u-star inversion over m/n trials.

For each tested m/n, solves on a coarse grid for the regularised field of
u-star (normalised uplift) that best reproduces the channel elevations, and
selects the m/n with the lowest misfit.

Parameters
----------
m_over_n_values : array-like of float
    Concavity values to test (non-empty).
reference_area : float, default=1.0
    Reference drainage area A0 (m^2) used for chi.
coarse_cell_size_metres : float, default=2000.0
    Cell size (m) of the coarse grid on which u-star is solved.
roughness_weight : float, default=1e4
    Weight of the roughness (smoothness) regularisation term.
diagonal_jitter : float, default=1e-8
    Small value added to the normal-equation diagonal for numerical
    stability.
use_mae : bool, default=True
    Select the best m/n by mean absolute error (True) or RMSE (False).

Returns
-------
ChiInversion
    Trial, summary and node tables, and u-star rasters.

```python
def invert_ksn(self, m_over_n_values: object, reference_area: float = 1.0, coarse_cell_size_metres: float = 2000.0, roughness_weight: float = 10000.0, diagonal_jitter: float = 1e-08, use_mae: bool = True) -> lsdtt3.ChiInversion
```

#### `nodes`

Per-pixel fitted profile table (pandas.DataFrame); see the class documentation for columns.

#### `segments`

Per-segment k_sn table (pandas.DataFrame); see the class documentation for columns.

#### `summary`

Per-basin k_sn summary table (pandas.DataFrame); see the class documentation for columns.

---

## `ChiInversion`

Results of ``ChiFit.invert_ksn``: per-basin m/n trials and u-star fields.

Returned by ``ChiFit.invert_ksn()``; not constructible directly.

Attributes
----------
nodes : pandas.DataFrame
    Per channel pixel of the best trial: basin_id, flat, receiver_flat, chi, observed_elevation, predicted_elevation, residual_elevation, u_star, coarse_cell_id.
summary : pandas.DataFrame
    Best trial per basin: basin_id, best_trial_id, best_m_over_n, mae, rmse.
trials : pandas.DataFrame
    One row per basin and m/n trial: basin_id, trial_id, m_over_n, mae, rmse, n_nodes, n_parameters, n_constrained_parameters.

Methods
-------
support_raster(basin_id: int = 0) -> Raster
    Return the number of channel nodes supporting each coarse cell of one basin.
u_star_raster(basin_id: int = 0) -> Raster
    Return the best-trial u-star field of one basin as a Raster.

Authors
-------
B. Gailleton

### Members

#### `nodes`

Per channel pixel of the best trial: basin_id, flat, receiver_flat, chi, observed_elevation, predicted_elevation, residual_elevation, u_star, coarse_cell_id.

#### `summary`

Best trial per basin: basin_id, best_trial_id, best_m_over_n, mae, rmse.

#### `support_raster`

support_raster(self, basin_id: int = 0) -> lsdtt3.Raster

Return the number of channel nodes supporting each coarse cell of one basin.

Parameters
----------
basin_id : int, default=0
    Zero-based basin index.

Returns
-------
Raster

```python
def support_raster(self, basin_id: int = 0) -> lsdtt3.Raster
```

#### `trials`

One row per basin and m/n trial: basin_id, trial_id, m_over_n, mae, rmse, n_nodes, n_parameters, n_constrained_parameters.

#### `u_star_raster`

u_star_raster(self, basin_id: int = 0) -> lsdtt3.Raster

Return the best-trial u-star field of one basin as a Raster.

Parameters
----------
basin_id : int, default=0
    Zero-based basin index.

Returns
-------
Raster

```python
def u_star_raster(self, basin_id: int = 0) -> lsdtt3.Raster
```
