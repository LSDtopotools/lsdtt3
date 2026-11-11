# Chi analysis

This page explains how `lsdtt-chi-analysis` computes the chi coordinate,
chooses the basins it analyses, divides chi profiles into segments, finds
knickpoints and estimates the concavity index m/n with the disorder metric.
Every parameter is listed in the
[chi analysis CLI reference](../cli_reference/cli_reference_lsdtt_chi_analysis.md);
the same analysis from Python is in [Chi and ksn](../python/chi_ksn.md).

The code is in `src/cpu/lsdtt_chi_extraction_cpu.{hpp,cpp}` (chi and basin
selection), `src/cpu/lsdtt_chi_network_analysis_cpu.{hpp,cpp}` (segment
fitting), `src/cpu/lsdtt_chi_knickpoints_cpu.{hpp,cpp}` (tributaries and
knickpoints) and `src/cpu/lsdtt_chi_disorder_cpu.{hpp,cpp}` (disorder).

---

## The chi coordinate

Under the stream power incision model, a channel at steady state has

```
S = k_sn * A^(-m/n)
```

where *S* is the channel gradient, *A* the drainage area, *m/n* the concavity
index and *k<sub>sn</sub>* the normalised channel steepness. Perron and Royden
(2013) integrate this upstream from a base level at *x<sub>b</sub>*:

```
z(x) = z(x_b) + k_sn * chi,     chi = integral from x_b to x of (A_0 / A(x'))^(m/n) dx'
```

χ (chi) is a transformed distance along the channel. Plotted against chi, a
channel with uniform steepness is a straight line through its outlet, with
gradient *k<sub>sn</sub>*. Changes in steepness are changes in gradient, so a
chi-elevation profile can be read without measuring local slopes, which are
noisy on a DEM.

**How lsdtt3 integrates it.** Chi is accumulated upstream along the D8 flow
path, from each basin outlet (where chi = 0) to every channel pixel. Each step
from a pixel to its downstream receiver adds

```
d_chi = dx * 0.5 * [ (A_0 / A_here)^(m/n) + (A_0 / A_receiver)^(m/n) ]
```

the trapezoidal rule, where `dx` is the distance between the two pixel centres
(the grid spacing, or √2 times it on a diagonal) and *A* is the D8 drainage
area. The trapezoidal rule is second-order accurate. A rule that evaluates the
integrand only at the upstream pixel overestimates chi on every step, because
the area always grows downstream; on a 30 m grid that bias is about 0.4 %.

**Keep A<sub>0</sub> = 1 m²** (`chi_A0`, the default). Only then does the
gradient of the chi-elevation profile equal *k<sub>sn</sub>*. Any other value
scales every gradient by *A<sub>0</sub><sup>m/n</sup>*: at *A<sub>0</sub>* =
10⁶ m² and *m/n* = 0.45 that is a factor of about 500.

**Channels.** A channel starts where the number of contributing pixels reaches
`threshold_contributing_pixels`. This is a pixel count, so the drainage area it
represents changes with the grid spacing. With
`write_chi_coordinate_raster: true` chi is also integrated over the hillslope
pixels of each basin, which gives the chi map used for comparing chi across
drainage divides.

---

## Choosing the basins

Basin outlets come from `basin_outlet_fname`, or, if none is given, the tool
detects the outlets of the drainage basins itself. A detected outlet that is
not on the edge of the DEM drains to nodata inside it, such as a coastline with
the sea masked. It is moved `n_pixels_move_outlet_upstream` pixels up its main
stem, so the basin starts clear of the ragged coastal fringe.

The basins then pass three filters, in this order. All the chi outputs (data
maps, fits, disorder, inversion) use the same filtered set, so a basin keeps
the same key in every file.

1. **Short channels.** A basin whose channel network has fewer than
   `min_channel_pixels` pixels is dropped: it is too short to divide into
   segments.
2. **Incomplete basins.** A basin is incomplete if any of its pixels borders
   an open edge of the DEM or a nodata pixel. Part of its drainage area may
   then lie outside the data, so its drainage areas, and its chi values, would
   be too small. Incomplete basins are always dropped.
3. **Size.** Of the complete basins left, only the `n_largest_basins` with the
   most pixels are kept (5 by default; 0 keeps them all). The kept basins stay
   in their original order.

To analyse a basin that fails the completeness test, enlarge the DEM so the
whole basin is inside it, or fill the interior nodata that it touches
(`interior_nodata_method: fill`).

---

## Dividing chi profiles into segments

With `fit_chi_networks: true` the channel network of each basin is fitted with
a piecewise-linear function of chi: a set of straight segments in
chi-elevation space, each with its own *k<sub>sn</sub>*. The fit is
continuous: segments meet end to end, and tributaries meet their parent
channel at the junction.

**Branches.** The network is first divided into branches. The main stem
(branch 0) runs from the outlet up the longest flow path. At each pixel the
next pixel upstream is the donor with the longest flow path above it. Every
other channel joining the main stem starts a branch, from its junction up its
own longest path, and so on up the network. The main stem is fitted
first. The tributaries of a branch are then fitted in descending order of flow
length, each together with its own tributaries before the next.

**Fitting one branch.** Each branch is fitted by recursive splitting:

1. Draw a straight line in chi-elevation space from the downstream end of the
   branch to its upstream end. The downstream end is fixed at its fitted
   elevation: the outlet's observed elevation for the main stem, or the
   fitted elevation at the junction for a tributary. The upstream end takes
   its observed elevation, raised if needed so the line does not slope
   downward upstream.
2. Find the pixel between them whose observed elevation is furthest from the
   line.
3. If that distance is more than `fit_critical_divergence_metres` (5 m by
   default), split the branch there. The split point is fixed at its observed
   elevation, and steps 1–3 repeat on each half.
4. Otherwise accept the line as one segment, with
   *k<sub>sn</sub>* = Δ*z* / Δχ.

The result is that no pixel lies more than `fit_critical_divergence_metres`
from the fit. A smaller value gives more, shorter segments that follow the DEM
more closely, including its noise; a larger value gives fewer segments and
keeps only the large changes in steepness. Choose it from the vertical
accuracy of the DEM: a few metres for a 30 m global DEM, less for lidar.

This is not the segmentation of Mudd et al. (2014), which LSDTopoTools2 used:
that method chooses the number and position of segments by statistical model
selection (AICc) over many randomly sampled subsets of the channel pixels.
lsdtt3's recursive split is deterministic and controlled by one tolerance in
metres.

The fit runs once for each m/n value tested. Its outputs, per m/n value, are
`<write_fname>_chi_fit_nodes_mn_*.csv` (every pixel with its observed and
fitted elevation, branch and flags), `_chi_fit_segments_mn_*.csv` (one row per
segment, with its *k<sub>sn</sub>*), `_chi_ksn_basin_summary_mn_*.csv`, and
optionally the segments as lines (`write_chi_fit_segment_lines`).

---

## Tributaries and knickpoints

A **tributary** is a branch of the fit: a channel from its head down to where it
joins its parent branch. `<write_fname>_chi_tributaries_mn_*.csv` lists them,
with their order (0 for the main stem, 1 for its tributaries, and so on) and
the drainage area where they join.

A **knickpoint candidate** is any pixel where two fitted segments meet:

- within a branch, where one segment meets the next segment downstream;
- at a tributary mouth (flagged `at_junction = 1`), where the tributary's
  lowest segment meets the segment of its parent branch at the junction.

Each candidate in `<write_fname>_chi_knickpoints_mn_*.csv` carries the
*k<sub>sn</sub>* above and below it, their difference `delta_ksn` (upstream
minus downstream: positive where the channel steepens upstream) and
`log_ratio_ksn` (the base-10 logarithm of their ratio), with its position,
chi, elevation, flow distance and drainage area.

Every segment boundary is a candidate, so most are small. To pick out the
major ones, `<write_fname>_chi_knickpoint_thresholds_mn_*.csv` suggests
thresholds from the distribution of |`delta_ksn`| and |`log_ratio_ksn`|
over the within-branch candidates:

- the 90th and 95th percentiles;
- the median plus three times the median absolute deviation (scaled by
  1.4826, so it estimates a standard deviation for normally distributed
  values);
- Otsu's threshold, the value that best separates the distribution into two
  classes.

It also gives an outlier threshold for segment *k<sub>sn</sub>* itself: the
upper quartile plus three times the interquartile range.

The tool does not choose a threshold for you: filter the knickpoint table with
the one that suits your question. The candidates depend on the segmentation,
so `fit_critical_divergence_metres` also sets the smallest knickpoint that can
be found: a step in the profile smaller than the tolerance does not split a
segment.

---

## Estimating m/n with the disorder metric

If m/n is right, every channel in a basin at steady state with uniform
erodibility and uplift collapses onto one curve in chi-elevation space,
because elevation is then a single function of chi. With the wrong m/n,
tributaries sit above or below the main stem at the same chi. The tool tests
a range of m/n values (`mn_min` to `mn_max` in steps of `mn_step`, 0.3 to 0.6
in steps of 0.02 by default; or the list in `m_over_n`), recomputing chi and
the fit for each.

**The disorder metric** measures how far the channels are from collapsing
onto one curve. lsdtt3 implements the disorder method of Mudd et al. (2018),
which follows Hergarten et al. (2016). Sort all the channel pixels of a basin
by their observed elevation. If chi rises with elevation without reversing,
the summed steps in chi equal the total range of chi. Each reversal adds to
the sum. lsdtt3 computes

```
D = ( sum_i |chi[i+1] - chi[i]|  -  (chi_max - chi_min) ) / (chi_max - chi_min)
```

which is 0 when the channels are perfectly collapsed and grows as they
scatter. Mudd et al. (2018) normalise by the maximum chi; the network includes
the outlet, where chi = 0, so the range and the maximum are the same. The m/n value with the lowest disorder is the best estimate for that
basin. Ties in elevation are broken by chi, so the result does not depend on
the order the pixels are stored in.

`<write_fname>_chi_disorder_sweep.csv` has one row per basin and m/n value,
with the disorder of the whole network. Look for the m/n value of each
basin's minimum.

**Uncertainty.** A basin's disorder depends on which tributaries it happens to
have. Three options (each needs `fit_chi_networks: true`) resample the side
branches and recompute disorder at every m/n value, always keeping the main
stem:

- `run_subset_sensitivity`: every combination of
  `subset_n_side_branches_per_sample` side branches (3 by default) together
  with the main stem, as in Mudd et al. (2018). The m/n value with the lowest
  disorder for each combination gives a population of best-fit values; Mudd
  et al. (2018) report its median and interquartile range. It is a
  sensitivity test rather than a formal uncertainty;
- `run_jackknife`: leave out one side branch at a time, giving a jackknife
  mean and standard error;
- `run_bootstrap`: `bootstrap_replicates` resamples of the side branches
  (with replacement by default), giving a mean, standard deviation and a 95 %
  interval (2.5th and 97.5th percentiles).

An m/n value whose interval overlaps the minimum's is not distinguishable from
it. The detailed samples are written per m/n value as
`_chi_disorder_subset_samples_mn_*.csv`, `_chi_disorder_jackknife_mn_*.csv`
and `_chi_disorder_bootstrap_mn_*.csv`; the summary statistics are in the
sweep table.

The method assumes uniform erodibility and uplift and a network near steady
state. Where lithology or uplift changes within a basin, the channels do not
collapse onto one curve at any m/n, and the minimum disorder is not a
concavity estimate. The [Lithology and channel steepness statistics](lithology_channel_statistics.md)
page covers k<sub>sn</sub> by rock type. `run_inversion` is a separate
estimate of m/n from a spatial inversion of k<sub>sn</sub>; see the CLI
reference.

---

## References

- Hergarten, S., Robl, J., and Stüwe, K. (2016). Tectonic geomorphology at
  small catchment sizes – extensions of the stream-power approach and the χ
  method. *Earth Surface Dynamics*, 4(1), 1–9.
  <https://doi.org/10.5194/esurf-4-1-2016>
- Mudd, S. M., Attal, M., Milodowski, D. T., Grieve, S. W. D., and Valters,
  D. A. (2014). A statistical framework to quantify spatial variation in
  channel gradients using the integral method of channel profile analysis.
  *Journal of Geophysical Research: Earth Surface*, 119(2), 138–152.
  <https://doi.org/10.1002/2013JF002981>. The segment-fitting algorithm of
  LSDTopoTools2; not used by lsdtt3.
- Mudd, S. M., Clubb, F. J., Gailleton, B., and Hurst, M. D. (2018). How
  concave are river channels? *Earth Surface Dynamics*, 6, 505–523.
  <https://doi.org/10.5194/esurf-6-505-2018>. The disorder method for
  estimating m/n, implemented by lsdtt3.
- Perron, J. T., and Royden, L. (2013). An integral approach to bedrock river
  profile analysis. *Earth Surface Processes and Landforms*, 38(6), 570–576.
  <https://doi.org/10.1002/esp.3302>. The chi coordinate.
