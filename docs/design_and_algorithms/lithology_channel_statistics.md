# Lithology and channel steepness statistics

This page describes the statistics behind two lithology features:

- **k_sn by lithologic unit** in `lsdtt-chi-analysis` (`lithology_fname`): the
  chi-elevation slope of each fitted channel, summarised per unit. Usage is in the
  [chi analysis reference](../cli_reference/cli_reference_lsdtt_chi_analysis.md#k_sn-by-lithologic-unit).
- **Lithology grouping** of channel-link bearings and junction angles in
  `lsdtt-channel-network-analysis`. Those statistics are on the
  [channel-link bearing statistics](channel_bearing_statistics.md) page.

The code is in `src/cpu/lsdtt_lithology_map_cpu.{hpp,cpp}` (reading a map and
assigning units), `src/cpu/lsdtt_lithology_channel_stats_cpu.{hpp,cpp}` (nodes,
segments, reaches, summaries and the bootstrap) and `src/utils/lsdtt_statistics`
(the weighted quantile and the median). The tests and the measured tolerances are in
the developer testing notes, sections 8c.3 and 8c.4.

## Units on the grid

A lithology map is read once, on the DEM grid. A polygon layer is rasterised; a
code raster is used as it is, or warped onto the grid with nearest neighbour when
its grid differs. Each pixel gets a unit code, or `-1` (unmapped) when a valid
pixel is touched by no polygon or is nodata in the code raster. A map that labels no
valid DEM pixel is refused, and a partial cover is reported with its unlabelled
fraction (a warning above 10 %).

Each channel node takes the unit of its pixel. The distance from each pixel to the
nearest unit boundary is computed once, so a node's distance to its nearest contact
is known. Nodes closer than `lithology_contact_buffer_metres` to a contact are set
aside from the per-unit statistics; they stay in the node table, flagged.

## Segments

The chi fit is piecewise linear, and each segment owns a set of nodes. A segment
has a length-weighted majority unit and a purity: the share of its length (weights
are each node's step length to its receiver) in that unit. Ties go to the lower code.
Segments are **not split at contacts**, so a segment that crosses a contact gives all
of its nodes the segment's k_sn. The node-level k_sn is therefore the segment value,
and the reach-level k_sn (below) is the one that respects the unit boundary.

## Lithologic reaches and the OLS k_sn

A reach is a run of consecutive nodes in one unit on one branch. It breaks when the
unit changes or a buffer is entered. For a reach with nodes indexed by k, the
elevation z_k is regressed on the chi coordinate chi_k by ordinary least squares:

    slope = sum((chi_k - chi_mean)(z_k - z_mean)) / sum((chi_k - chi_mean)^2)

The slope is `k_sn_ols`, and it is the reach's k_sn (the chi-elevation slope, as in
the fitted segments). The standard error is `k_sn_ols_se = sqrt(SSR / (n - 2) / Sxx)`,
where SSR is the residual sum of squares, n the node count and Sxx the chi sum of
squares. `r_squared = 1 - SSR / SYY`. The SE is NaN below three nodes.

A second estimate, `k_sn_fit_chi_weighted`, is the mean of the nodes' segment k_sn
values weighted by the chi interval each node spans (only intervals with positive
chi span count). It is written for comparison with `k_sn_ols`.

A reach **qualifies** when it has at least `lithology_min_reach_pixels` nodes (the
default 10; the CLI refuses fewer than 2, since a slope needs two points), a positive
chi span and a finite slope.

The reach end types are `source`, `junction`, `outlet`, `contact` and `buffer`. They
describe why the reach stops, and they are written so that a reach can be filtered
by how it ends (for example, reaches that end at a contact can be excluded).

## Summaries

Per unit, the node statistics use the nodes not set aside by the buffer: the median,
the 25th and 75th percentiles (the code's interpolated quantiles), the mean, the
minimum and maximum, and the length-weighted mean and quartiles. The reach statistics
use the qualifying reaches: their median and quartiles (interpolated), and the
length-weighted median, weighted by reach length. `low_sample` is set when fewer than
five reaches qualify, and in that case no bootstrap is run.

### Median convention (D2)

The unweighted median is the middle value for an odd count, and the mean of the two
middle values for an even count (`median_in_place_impl` in `lsdtt_statistics`). The
interpolated quartiles use the shared `quantile_interpolated_in_place` definition.

### Weighted quantile

The length-weighted quantiles use `compute_weighted_quantile` in
`src/utils/lsdtt_statistics.cpp`. It is a **lower** weighted quantile:

1. Sort the values ascending (a stable sort, so equal values keep their input order).
2. Accumulate the weights in that order. The total is the sum over all values.
3. Return the first value whose cumulative weight is at least `q * total`.

So `q = 1` always returns the largest value. A zero-weight value is returned only if it is
the first point at which the cumulative weight reaches the target (for example, at q = 0). Weights must be finite and not
negative, and a total weight of zero is an error. With equal weights the function
agrees with the floor quantile (`compute_quantile_floor_in_place`) at q = 0, 0.5 and 1
for sizes 1 to 12; the unit test pins this, and the fixed-vector cases at q = 0, 0.25,
0.5 and 1.

Because the quantile is a value from the data (not an interpolation), a length-weighted
median is always one of the observed k_sn values.

## Bootstrap interval

The interval on the per-unit median reach k_sn is a percentile bootstrap:

1. The qualifying reaches of the unit are the sample. With `lithology_bootstrap_cluster`
   `reach`, each reach is a cluster. With `basin`, the reaches of one basin form one
   cluster (clusters are ordered by basin index before the draw).
2. Draw as many clusters as there are, with replacement, and pool their reaches.
3. Take the median of the pooled reach k_sn values.
4. Repeat `lithology_bootstrap_replicates` times (default 1000).
5. Report the 2.5 and 97.5 percentiles of the replicate medians (interpolated).

The random stream is a `std::mt19937_64`. In `lsdtt-chi-analysis` it is seeded per
m/n value with `std::seed_seq{seed, round(m/n * 10000)}` (`lithology_bootstrap_engine`).
A value's interval therefore depends only on the seed and that m/n value, and the
per-m/n output of a run with several m/n values is byte-identical to the run with one.
In `lsdtt-channel-network-analysis` the junction-angle bootstrap uses one stream per
run.

A reach-level interval treats reaches as independent. Reaches on one trunk channel are
correlated, so the reach interval is a **lower bound** on the uncertainty. A basin
cluster is the coarser check, with fewer clusters and so a wider interval.

## Measured tolerances

The test set in section 8c.3 of the testing notes is a steady-state known answer.
The DEM is 120 by 160 cells at 30 m, with two erodibilities K1 = 1e-5 and K2 = 4e-5,
uplift 2e-4 and a contact in one of two layouts. For n = 1 (m = 0.45) the expected
medians are (U/K)^(1/n), 20 and 5; for n = 2 they are sqrt(20) and sqrt(5). The
measured relative errors (identical for n = 1 and 2) were:

| Statistic | Layout | Unit K1 | Unit K2 | Ratio of medians |
|---|---|---|---|---|
| node median | A (contact parallel to drainage) | +1.14 % | +1.25 % | 0.108 % |
| reach median (OLS, 10 nodes, buffer 0) | A | +1.759 % | +1.310 % | 0.443 % |
| reach median | B (contact across drainage) | +1.192 % | +1.581 % | 0.383 % |

The medians sit high because the steady state integrates chi at each cell's own area
and the chi integral is evaluated trapezoidally. Their tests use a 2 % cap, and the
ratio tests a 1 % cap set at about twice the measured error. Through the CLI, which
re-routes flow on the steady DEM, the reach medians are 20.403 and 5.163, a ratio of
3.9515 against the true 4 (1.212 % from 4; the test is pinned at 2.42 %).

Reach counts are pinned exactly: layout A has 46 (K1) and 42 (K2) reaches (16 and 14
qualifying); layout B has 57 and 38 (18 and 20 qualifying). Layout B has 20
contact-ended reaches against 8 for A.

Two smaller checks are in the same test set: `lithology_min_reach_pixels` is honoured,
and the OLS on an exact line z = 3 + 7 chi returns k_sn 7 with R^2 1 and SE 0.

## Field example

The West Virginia Allegheny Front run (live data, not a CTest) is the measured
per-unit table in the chi reference page. It is not a validation of lithology
erodibility: the units are map units, and the caveats in that page (uniform uplift,
knickzones, strike-parallel sampling, alluvium, map-scale position) apply.
