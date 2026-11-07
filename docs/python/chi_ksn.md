# Chi and ksn

Chi analysis builds on the D8 channel network. The normal sequence is DEM,
`FlowInfo`, a channel scope/network, chi coordinates, and a piecewise-linear
fit.

```python
import lsdtt3

dem = lsdtt3.fill_pits(lsdtt3.Raster.read("dem.tif"))
flow = lsdtt3.FlowInfo(dem)
network = flow.channel_network(flow.channel_scope(threshold_contributing_pixels=500))

chi = network.chi(m_over_n=0.45)
fit = chi.fit(dem, critical_divergence_metres=10.0)

nodes = fit.nodes        # pandas DataFrame: one row per channel node
segments = fit.segments  # pandas DataFrame: one row per fitted segment, with ksn
summary = fit.summary    # pandas DataFrame: one row per basin
```

Without `outlets`, chi is computed upstream of every outlet of the network.
To analyse chosen basins, pass their outlet pixels as an N by 2 array of
`(row, col)`: `network.chi(outlets=outlet_pixels, m_over_n=0.45)`. The same
array gives the basins themselves with `flow.basins(outlet_pixels)` (see
[Flow, channels, and basins](flow_networks.md)).

`critical_divergence_metres` is how far, in metres of elevation, the data may
move from a straight line in chi-elevation space before a segment is split.
The default is 5 m; 10 m suits 30 m DEMs. Much smaller values (the 0.1 m in
the lidar example below) split the profile into many short segments.

The tables are plain pandas DataFrames without geometry. `nodes` has the
node's flat pixel index (`flat`; row = `flat // n_cols`, col =
`flat % n_cols`), `chi`, `observed_elevation`, `fitted_elevation` and
`branch_id`; `segments` has `downstream_chi`, `upstream_chi` and `ksn` for
each segment. Use the segment table for mapped `ksn`; join or assign its values
to the node table for point-based plotting. See the repository example
`python_bindings/examples/chi_ksn.py` for a map plus chi-elevation panel.

The basin column is `basin_id` and the steepness column `ksn`.
`lsdviztools3`'s chi plots expect the command-line programs' names,
`basin_index` and `k_sn`, so rename them first:

```python
from lsdviztools3.io.chi import assign_node_ksn
from lsdviztools3.charts import plot_chi_stacked

chi_nodes = assign_node_ksn(
    fit.nodes.rename(columns={"basin_id": "basin_index"}),
    fit.segments.rename(columns={"basin_id": "basin_index", "ksn": "k_sn"}),
)
fig, ax = plot_chi_stacked(chi_nodes)
```

## Tributaries and knickpoints

The fit also gives the tables that `lsdtt-chi-analysis` writes with
`write_chi_tributary_table` and `write_chi_knickpoint_table`, made by the same
C++ routines:

```python
tributaries = fit.tributaries               # one row per tributary (fit branch)
knickpoints = fit.knickpoints               # one row per node where two segments meet
thresholds = fit.knickpoint_thresholds      # suggested thresholds for a major knickpoint
```

A tributary is a branch of the fit, from its head down to the node where it
joins its parent; branch 0 of each basin is the main stem.

- **`tributaries`**:
  - `junction_drainage_area_m2`, the drainage area just before it joins, and
    `area_fraction`, that area over the basin's;
  - `area_rank_in_basin`, with 1 the largest tributary;
  - `ksn_min`, `ksn_median`, `ksn_max`;
  - `parent_branch_id`, and `n_descendants`, the tributaries draining into it.
- **`knickpoints`**: `ksn_upstream`, `ksn_downstream`, `delta_ksn` (upstream
  minus downstream, positive where the channel steepens upstream),
  `log_ratio_ksn`, `x`, `y`, `chi`, `elevation`, and `at_junction` (1 at a
  tributary mouth, where a change in k_sn is not a knickpoint).
- **`thresholds`**: for |delta_ksn| and |log_ratio_ksn|, the 90th and 95th
  percentiles, the median plus three MADs and Otsu's split, all over the rows
  with `at_junction` 0. Also an outlier value for segment k_sn
  (`q3_plus_3iqr`).

`source_key` is the channel key of the chi data maps. It is on every row of
`fit.nodes`, `fit.segments`, `fit.tributaries` and `fit.knickpoints`, so a
tributary can be followed across the tables. It is -1 in a basin whose upslope
area touches nodata.

To choose the major knickpoints by the Otsu threshold, for example:

```python
t = thresholds.query("measure == 'abs_delta_ksn' and method == 'otsu'").threshold.iloc[0]
major = knickpoints[(knickpoints.at_junction == 0) & (knickpoints.delta_ksn.abs() > t)]
```

lsdviztools3 draws these tables (`select_tributaries`, `lsdviz3-knickpoints`)
after the same `basin_id` to `basin_index` rename as above.

See [Chi and ksn API](../python_api/chi.md).
