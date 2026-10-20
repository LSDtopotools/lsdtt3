# Channel-link bearing and junction-angle statistics

This page describes the statistics that `lsdtt-channel-network-analysis` computes on
its channel-link bearings and junction angles when `write_channel_bearing_statistics`
or `write_junction_angle_statistics` is on. Usage is in the
[channel network analysis reference](../cli_reference/cli_reference_lsdtt_channel_network_analysis.md#step-7b-bearing-and-junction-angle-statistics-by-group).

The code is in `src/cpu/lsdtt_channel_bearing_statistics_cpu.{hpp,cpp}` (circular
summaries, rose bins, Holm, the comparisons and the CSV writers),
`src/cpu/lsdtt_channel_network_grouping_cpu.{hpp,cpp}` (the grouping keys and the
block labels), and the two-sample numerics in `src/utils/lsdtt_statistical_tests`,
which are unchanged by this work. The measured record is in the developer testing
notes, sections 10s.1 to 10s.6.

## Circular summaries

A bearing is an angle theta in degrees with a period of 360 (directional) or 180
(axial). Each value becomes a unit vector, and the summaries are taken on those
vectors.

Every periodic angle the CSVs carry (link bearings and axial bearings, circular means,
the optimal rotation) lies on the half-open domain [0, 360) or [0, 180). The writers
print 10 significant digits, at which a value just below the period (for example
359.99999999999579) would round to the period; such a value is written as 0.

- **Axial doubling.** For the axial form, each angle is doubled before the vectors are
  formed, so theta and theta + 180 map to the same point. The mean direction is then
  halved, which puts it in [0, 180). A four-fold D8 pattern (0, 45, 90, 135 degrees)
  cancels to zero on doubled angles; this is why the axial mean resultant length is
  robust to D8 quantisation.
- **Mean direction and mean resultant length R-bar.** R-bar is the length of the mean of
  the unit vectors: 1 if all values agree, 0 if they cancel. Below
  `k_channel_bearing_undefined_resultant` (1e-10) R-bar is treated as exactly zero and
  the mean direction is NaN, since the cancellation is round-off (four D8 directions
  cancel to about 1e-17, not 0). A weighted summary uses the weights in the same way.
- **Circular variance** is 1 - R-bar. The **axial SD** is half the doubled-angle SD, so
  it is on the same scale as the axial angle.
- **Kish effective sample size** for weights w: n_eff = (sum w)^2 / (sum w^2). It is
  given for the length-weighted summary (`n_eff_length`), because a few long links can
  dominate a sum of weights.
- **Rayleigh test.** The statistic is Z = n R-bar^2. The p-value uses Zar's approximation

      p = exp( sqrt(1 + 4n + 4(n^2 - (n R-bar)^2)) - (1 + 2n) )

  clamped to [0, 1]. It assumes independent records, so it is named
  `rayleigh_p_iid_approx`.
- **D8 fraction.** The fraction of values within 0.5 degrees of a multiple of 45
  degrees (`fraction_on_d8_multiples`). A D8 network gives values exactly on these
  multiples for short links, so this fraction shows how much of a group is quantised.

## Rose bins

The rose bins are centred on 0 degrees: bin k covers [k w - w/2, k w + w/2), wrapped,
for bin width w (`channel_bearing_rose_bin_width_degrees`, default 10). Because 0, 45 and
90 degrees are then bin centres, a D8 direction is counted in its own bin rather than
split across two. The width must divide 180 exactly. The axial rose folds the angles
first, so 175 degrees is counted in the bin centred on 0.

## Two-sample comparisons

For bearings the two groups are compared with:

- **Circular Wasserstein distance** (W1, in degrees, the cost is the shorter arc). With
  `channel_bearing_rotation_invariant`, one global rotation is also optimised; the
  cost grows with n_a times n_b, and pairs above 250000 are skipped.
- **Kuiper's V** on the circle, with `kuiper_p_iid_approx`. The p-value assumes
  independent records.

For junction angles the sectors (alpha, beta_min, beta_max) are **angular magnitudes**,
not directions. A sector of 10 degrees and one of 350 degrees are 340 degrees apart
in every sense that matters for a sector, so the junction comparisons use linear
statistics: the two-sample Kolmogorov-Smirnov statistic D (`ks_D`, `ks_p_iid_approx`)
and the linear Wasserstein distance (`linear_wasserstein_deg`). Using a circular
distance here would make the comparison wrap around where the physics does not.

When a family has several comparisons, the iid p-values are adjusted with **Holm's
step-down** method (`kuiper_p_holm`, `ks_p_holm`), per family: per stratum and per
angular type. NaN p-values are not counted. A reference comparison (against another run's
bearings) is a family of its own.

## Lithology grouping, purity and descriptive-only groups

A link takes the length-weighted majority unit of the nodes its bearing was fitted on,
with its purity (the majority's share of that length). Links below the purity threshold
form `mixed (<X% one unit)`. A junction takes its pixel's unit (`lithology`), and,
for `lithology_interior`, its unit only when its three sampled arms agree; otherwise
it is `arms disagree`.

These two groups are marked **descriptive only**. They are summarised and roses are
written, but they are excluded from every comparison, from the one-vs-rest "rest"
group, and from comparisons with a reference. A mixed link has no single lithology, so
comparing it with a unit would put a contact effect into the contrast. The excluded
count is written to stdout.

## Block order and block permutation

### Why blocks

Links that share a valley are not independent samples. A group that is spatially
clustered can look significant under an iid test simply because its records share a
few basins. The blocks are the units of quasi-independence, and the reported `n_blocks`
counts them beside `n`.

### The block rule

A block is an order-k basin (the Strahler basin of the chosen order k), and a record
belongs to the basin that holds it. The **automatic** rule (`channel_bearing_block_order`
= 0) chooses the largest order k whose order-k basins, holding at least one analysed
record, number at least `channel_bearing_min_blocks` (default 30). If no order reaches
the minimum, no blocks are defined. A positive value is used as given.

Records above the chosen order have no order-k basin. They are summarised and enter
the iid comparisons, but take no part in the block permutations, and the counts are
given in `n_unblocked`.

### Why 30 blocks

The smallest p-value a permutation test can give is 1 / (number of distinct label
assignments + 1). With B blocks and two groups in a split of s and B - s, the number of
label assignments is C(B, s). With 30 blocks even a 5-vs-25 split has C(30, 5) =
142506 assignments, so the 999 permutations (not the block count) set the resolution.
With about 10 blocks an unbalanced split can have fewer than 100 assignments, and the
p-value is then coarse. The default of 30 follows from that argument: with fewer blocks the
attainable p-values are coarse, and with 30 the number of permutations, not the number of
blocks, limits the resolution. `channel_bearing_min_blocks` can change it.

### The block permutation p-value

For a comparison with observed statistic T (Kuiper V for bearings, KS D for junction
angles), the block permutation p is

    p = (1 + #{T* >= T}) / (1 + P)

where P is `channel_bearing_permutations` (default 999) and T* is the statistic
recomputed under each permutation. What is permuted depends on the key:

- **Spatial keys** (basin, region, lithology). Each block takes the majority group of
  its records in the comparison (ties go to group a). The block labels are then
  permuted among the blocks, so a whole valley changes group together. T* is computed
  on these block-majority samples. The reported `permutation_statistic` is T on the same
  samples, and it is therefore not the record-level statistic.
- **Order keys** (strahler_order, donor_order_pair). The group labels are permuted
  among the records within each block. Block membership is kept, so records of a
  basin remain in that basin.

A p is NaN, with a note, when a group has no blocked record, when every block has the
same majority group, or, for an order key, when no block holds both groups.

Each comparison row draws from its own `std::mt19937_64`, seeded from
`channel_bearing_permutation_seed` and the row's own labels. A row's p therefore does
not change when other comparisons are added or removed from the run.

### Reading the p-values

The file carries the iid and the permutation p-values side by side, always beside the
block counts. The labels say what each is: `kuiper_p_iid_approx` or `ks_p_iid_approx`
assume independent records; `permutation_p` is the block permutation above. The
permutation p does not reach a p below 1 / (P + 1), so with 999 permutations the
smallest value is 0.001.

## Measured behaviour

The measurements are in the testing notes, section 10s; the figures here are copies of
those records. They are measurements, and the page does not draw a conclusion from them.

**Isotropic null (section 10s.4).** On a D8 network of a diamond-square surface (257 by
257 at 30 m, roughness 0.55, sd 100 m, threshold 20 px), the all-link axial R-bar for
seeds 1 to 8 was 0.0751, 0.0063, 0.0523, 0.0191, 0.0457, 0.0272, 0.0369 and 0.0149.
The CTest checks a bound of 0.15 on seed 2. Directional R-bar on the same surfaces was
0.09 to 0.36. The axial form is the one that is bounded, because the directional form
is not a usable isotropy statistic on a D8 network.

**Allegheny Front (section 10s.5 and 10s.6).** With the automatic block order, the
network chose order 3 (108 basins with records; order 4 had 23). The pairwise
comparisons by order, axial, with the iid and permutation p-values (999 permutations):

| run, type | a vs b | n | n_blocks | Kuiper V | kuiper_p_iid_approx | permutation_p |
|---|---|---|---|---|---|---|
| wvbig, axial | 1 vs 2 | 1765 / 791 | 108 / 106 | 0.056 | 0.379 | 0.329 |
| wvbig, axial | 1 vs 3 | 1765 / 311 | 108 / 83 | 0.107 | 0.046 | 0.012 |
| wvbig, axial | 2 vs 3 | 791 / 311 | 106 / 83 | 0.094 | 0.244 | 0.539 |
| wvbig, axial | x vs 4+ | - / 107 | - / 0 | 0.227 to 0.309 | < 0.001 | nan |
| wvbig, directional | 1 vs 3 | 1765 / 311 | 108 / 83 | 0.142 | 0.0008 | 0.006 |

Order 4+ has no order-3 block (n_blocks 0), so its permutation p is nan. The isotropic
null with the same parameters chose order 4 (33 basins with records). Its axial
comparisons have p values ranging from 0.241 to 0.965 under permutation, and the
null's junction-angle KS comparison of order 2 with order 4+ (D 0.126) has a permutation
p of 0.044 (iid p 0.021).

**Threshold and sampling sensitivity (section 10s.5).** Varying the threshold (150, 300
and 600 px) and the minimum nodes per sampled link (2 against 5) changed the axial R-bar by
order (the threshold 600 run gave 0.038 for order 3 against 0.148 for order 2), and the
D8 fraction for order 1 rose from 0.049 to 0.121 when the minimum fell to 2 nodes. The
Allegheny Front runs used 5 nodes.

## Tests

The unit tests in `tests/cpp/test_channel_bearing_statistics.cpp` (labelled `unit;network`)
pin the circular and linear numerics on hand-built cases, including Holm, the
comparison families, the grouping and the isotropic null. The lithology grouping is
tested in `unit_lithology_channel_network_stats` (see section 8c.4 of the testing notes).
The CLI smoke and rigorous tests run the statistics on a synthetic trellis and on a
reference file. The block order and the block permutation are tested with the lithology
tests (section 8c.4 of the testing notes).
