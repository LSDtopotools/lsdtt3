# Floodplains, terraces and valleys in lsdtt3

This page explains the design choices behind lsdtt3's fluvial surface tooling: what
`lsdtt-flat-surface-finder`, `lsdtt-fluvial-surfaces` and `lsdtt-valley-centreline`
each do, why the work is split across three programs rather than one, which
algorithms are used at each step, and where those algorithms come from. It is not a
CLI reference — for parameter names, types and defaults see the reference pages for
[lsdtt-fluvial-surfaces](../cli_reference/cli_reference_lsdtt_fluvial_surfaces.md) and
[lsdtt-valley-centreline](../cli_reference/cli_reference_lsdtt_valley_centreline.md),
or run any of the three with `-h`.

The methods here descend from two papers, **Clubb et al. (2017)** on geomorphometric
delineation of floodplains and terraces, and **Clubb et al. (2022)** on continuous
measurement of valley floor width. Both were implemented in LSDTopoTools2 as the
single driver `lsdtt-valley-metrics`. lsdtt3 keeps the parts of that lineage that
have held up, replaces the parts that had accumulated repairs, and says explicitly
where it departs. Those departures are the substance of this page — if you are
porting a `lsdtt-valley-metrics` parameter file, read
[From lsdtt-valley-metrics](#from-lsdtt-valley-metrics) first.

---

## Overview: three tools, one job each

| tool | question it answers | what it needs |
|---|---|---|
| `lsdtt-flat-surface-finder` | **which ground is flat?** | a DEM |
| `lsdtt-fluvial-surfaces` | **which fluvial level is each flat surface on?** | a channel network, or a supplied height field |
| `lsdtt-valley-centreline` | **where is the centreline of this valley, and how wide is it?** | a valley mask |

The split is placed at exactly one point: **where a channel network becomes
necessary**. Deciding that a piece of ground is flat is pure morphometry and needs
nothing but the DEM. Deciding that a flat surface is a floodplain rather than a
terrace tread is a statement about its height above a channel, and needs flow
routing. Tracing a centreline needs a valley to trace it through, but does not need
to know how that valley was found.

The normal pipeline:

```
DEM
 |
 v
lsdtt-flat-surface-finder    ->  <stem>_combined_flat_mask.tif
 |                                     (flat_mask_fname)
 v
lsdtt-fluvial-surfaces       ->  <stem>_band_index.tif
                                 <stem>_banded_labels.tif
                                 <stem>_surface_class.tif
                                 <stem>_surface_patches.csv
                                 <stem>_ladder_windows.csv
                                 <stem>_valley_mask.tif
 |                                     (valley_mask_fname)
 v
lsdtt-valley-centreline      ->  <stem>_valley_centreline.fgb
                                 cross sections, widths, hydraulics
```

Each stage writes a file the next stage reads, so you can stop at any of them,
substitute your own product (a hand-digitised valley mask, a HAND raster from
elsewhere), or iterate on the expensive middle stage without recomputing the cheap
first one.

`lsdtt-fluvial-surfaces` will run without `flat_mask_fname`, falling back to a
built-in relief or slope calculation. That is for a first look. `lsdtt-valley-centreline`
will **not** run without `valley_mask_fname` — it exits with an error naming
`lsdtt-fluvial-surfaces` — because a centreline traced through a valley the tool
invented for itself is a result that cannot be audited.

---

## Where this comes from

### Clubb et al. (2017): floodplain and terrace identification

The 2017 method identifies floodplain and terrace pixels from two per-pixel
quantities:

- **local gradient `S`**, from a polynomial surface fitted in a circular window
  whose radius is chosen from breaks in the standard deviation and interquartile
  range of curvature with window size (Grieve et al. 2016);
- **relief relative to the nearest channel `R_c`**, the elevation difference between
  the pixel and the nearest channel pixel reached by steepest-descent flow routing,
  where "channel" means a link at or above a user-set Strahler order (third order is a good bet). 

A pixel is floodplain or terrace when **both** are below threshold. The paper's
contribution is that neither threshold is set by hand: both are read off a
**quantile–quantile plot** against a Gaussian fitted between the 25th and 75th
percentiles of the data, taking the threshold as the lowest value at which the
empirical and Gaussian distributions differ by less than 1 % of the data range —
the same construction Lashermes et al. (2007) and Passalacqua et al. (2010) use for
hillslope–valley transitions. Terraces are then separated from floodplain by a
single user-set height above the modern channel.

### Clubb et al. (2022): valley floor width measurements

The 2022 method takes the 2017 valley floor mask and measures its width
continuously. It needs a longitudinal line to measure orthogonal to, and using the flow routed
channel itself can result in errors where the river meanders inside its valley. So it
builds an artificial V-shaped valley and lets flow routing find the axis:

1. For every valley floor pixel, find the distance `d_w` to the nearest valley wall
   pixel, and the elevation of that wall pixel.
2. Set each valley pixel to `z_min`, the minimum wall elevation in a circular window
   of roughly the valley width.
3. Depress it by distance from the wall:

   ```
   z_t = z_min - alpha * d_w
   ```

   which is deepest along the geometric centre of the valley.
4. Breach (Lindsay 2016) and fill, and repeat steps 3–4 — five iterations in their
   analysis — until a single carved channel remains.
5. Route D8 down the trough. **That flow path is the centreline.**

They chose this over a morphological skeleton explicitly because it preserves flow
routing information on the centreline pixels. Width is then measured on rays
orthogonal to a smoothed centreline bearing, `b = 0.5 * (b_u - b_d) + b_d`, where
`b_u` and `b_d` are the bearings to points `n` pixels up- and downstream, walking
out to left and right until the ray leaves the valley mask.

### The legacy implementation

Both papers live in `lsdtt-valley-metrics` in LSDTopoTools2, driving `LSDFloodplain`
and `LSDTerrace`. lsdtt3 reproduces its trough construction closely enough that a
legacy parameter file's `centreline_loops`, `trough_scaling_factor` and
`min_bank_elevation_window_radius` transfer with the same meaning. The thresholding
and the surface classification, by contrast, are deliberately different, and the
next three sections say why.

---

## `lsdtt-flat-surface-finder`: which ground is flat

Pure morphometry. No flow routing, no channel network, no height reference — the
program does not link against the flow layer at all.

### Admission: relief, gradient, or both

`admission_algorithm_choice` selects the test that admits a pixel:

| choice | test |
|---|---|
| `relief` (default) | local relief in a window of radius `relief_window_radius_metres` below `max_relief_metres` |
| `gradient` | gradient (m/m) below `max_admission_gradient` |
| `relief_and_gradient` | both |
| `relief_or_gradient` | either |

**These are not the same test and they reject different ground.** Local relief
charges a surface for its own tilt, so a smooth terrace tread that slopes downvalley
fails a relief test while a rough but level surface passes it; a gradient cap does
the opposite. On steep terrain the two admit nearly disjoint sets. Clubb et al.
(2017) used gradient; the historical lsdtt3 behaviour used relief; offering both,
and their conjunction and disjunction, is the honest resolution. Use
`relief_or_gradient` when the point is to find every candidate surface and let the
ladder in the next tool sort them out.

The admission gradient is measured from a polynomial surface fitted over
`admission_gradient_window_radius_metres` rather than from a 3×3 difference, following
the same reasoning as the 2017 paper: a three-cell operator on a 30 m DEM measures
pixel noise, and the bound has to be read at the scale of the surface being admitted.
There is no default bound worth guessing — a gradient of 0.087 m/m is a
valley floor, 0.27 m/m admits fan and pediment — so run the census at more
than one bound and report the sensitivity rather than defending one number.

### Choosing the threshold: QQ plot or CDF divergence

`relief_threshold_algorithm_choice` takes `fixed`, `cdf_divergence` or `qq_gaussian`. `qq_gaussian` is the QQ rule of Clubb et al. (2017), described under [The clubb2017 method](#the-clubb2017-method). `cdf_divergence` picks the relief bound from the
data, by finding where the empirical CDF diverges most from a fitted parametric
distribution, and
writes the whole curve to `<stem>_relief_adaptive_cdf.csv` so you can see the break
you are relying on. `relief_threshold_multiplier` scales the result.

This is the same *idea* as the 2017 QQ threshold — let the distribution say where
the process domain changes, rather than the user — implemented on the CDF rather
than the quantile–quantile plot, and against a family of fitted distributions
(normal, exponential, Weibull, Pareto; see `src/utils/lsdtt_distributions.hpp`)
rather than a Gaussian pinned at two percentiles. The 2017 construction assumes the
low tail of the population is Gaussian and reports the departure from that; where a
landscape's flat ground is genuinely a separate population, the divergence of the
fitted CDF locates the same break without committing to normality.

It remains a **diagnostic you should look at**, not an oracle. The 2017 paper is
clear that the method is most accurate where there is real contrast in slope and
elevation between the feature and its surroundings; the same caveat applies here,
and it is the reason `lsdtt-fluvial-surfaces` can take a supplied mask from any
source.

### After admission

Morphological opening and closing (`flat_opening_radius_metres`,
`flat_closing_radius_metres`, `open_then_close`) clean the mask; components below
`min_component_pixels` are dropped; the survivors are labelled. Plane fitting to
each labelled surface is optional on top of that.

The break-of-slope (VEC) skeleton that once lived here could *break* the flat mask
before labelling, with endpoint bridging to re-close gaps in it. Both existed to
stop two terrace treads merging into one component through a dissected riser. The
band-then-label ordering in the next tool makes that merge impossible, so riser
search is now a morphometric product of its own: `lsdtt-feature-finder`
`feature_type: find_risers`, which writes open riser segments with a step height.
Its old keys stop the run (output switches) or warn (tuning keys) with the new
`riser_*` name.

The headline output is `<stem>_combined_flat_mask.tif`, which is
`lsdtt-fluvial-surfaces`' `flat_mask_fname`.

---

## `lsdtt-fluvial-surfaces`: which level is each surface on

This is where a channel becomes necessary, and where lsdtt3 departs furthest from
the legacy implementation.

### Method presets and the four stages

`lsdtt-fluvial-surfaces` classifies the ground in four stages. Each stage has its own parameter, and a method preset (`fluvial_surface_algorithm_choice`) sets all of them at once.

| Stage | Parameter | What it decides |
|---|---|---|
| 1. Admission | `admission_algorithm_choice` | Which ground is considered at all: `flat_mask`, `relief`, `gradient`, or `floodplain_thresholds` (the floodplain bounds on HAND and gradient) |
| 2. Reference field | `reference_field_algorithm_choice` | The height the ground is measured from: `hand`, `relative_elevation` or `fitted_surface` |
| 3. Levels | `level_algorithm_choice` | How the admitted ground is split into levels: `ladder` (breaks found in the data), `fixed_thresholds` (a floodplain bound and a terrace bound) or `supplied_levels` (read from a file) |
| 4. Segmentation | `segmentation_algorithm_choice` | How the admitted ground becomes patches: `within_bands` or `whole_mask` |

Two more stage parameters sit beside these. `channel_definition_algorithm_choice` says which pixels count as channel (`contributing_pixels`, or `stream_order` with `min_channel_stream_order`). `terrace_reference_field_algorithm_choice` says what a terrace's height is measured from (`same` as the reference field, or `swath_baseline`). The floodplain bounds are chosen by `floodplain_threshold_algorithm_choice` (`fixed` or `qq_gaussian`).

The three presets set these stages as follows. The reference field is not set by a preset, so it is `hand` unless `reference_field_algorithm_choice` says otherwise.

| Stage | `ladder` (default) | `clubb2017` | `cross_section_levels` |
|---|---|---|---|
| Admission | `flat_mask` | `floodplain_thresholds` | `flat_mask` |
| Floodplain bounds | `fixed` | `qq_gaussian` | `fixed` |
| Channel definition | `contributing_pixels` | `stream_order` | `contributing_pixels` |
| Levels | `ladder` | `fixed_thresholds` | `supplied_levels` |
| Terrace field | `same` | `swath_baseline` | `same` |
| Segmentation | `within_bands` | `whole_mask` | `within_bands` |

- `ladder` is the data-driven method described in the sections that follow. It needs a flat mask (`flat_mask_fname`, from `lsdtt-flat-surface-finder`) or the built-in relief admission.
- `clubb2017` is the method of Clubb et al. (2017), with the bounds chosen from quantile-quantile plots. See [The clubb2017 method](#the-clubb2017-method) below.
- `cross_section_levels` is the ladder, but with the levels read from `supplied_levels_fname`, normally the `<stem>_cross_section_levels.csv` written by `lsdtt-valley-centreline` (see below).

Any stage named explicitly in the parameter file overrides the preset, so `clubb2017` with `level_algorithm_choice: ladder` bands the Clubb admission with the ladder. Several presets can run in one go by listing them comma-separated with no spaces (`ladder,clubb2017`). Each writes under `<write_fname>_<preset>`, and the comparison outputs are written as well (see [Method records and comparison outputs](#method-records-and-comparison-outputs)).

### The reference field

**`hand`** (default) — height above the nearest drainage (Nobre et al. 2011),
computed on the flow graph against the in-scope channel network, so each pixel
carries the donor channel node it was measured from. This is the direct descendant
of the 2017 paper's `R_c`, and it is the right field for a terrace flight: treads
along a reach sit at a roughly constant height above the channel, which makes a
flight a step function in HAND.

**`fitted_surface`** — the signed residual to a surface fitted to the flat ground,
positive above and negative below. **A fan is not a step function in HAND.** It
climbs continuously from toe to apex, so a HAND ladder cuts one surface into arcs
banded by distance from the apex rather than by age. What *is* roughly constant on a
fan is the depth of an inset lobe below the projection of the older surface — which
is what the residual measures. The sign convention is deliberately the same
direction as HAND, so band 0 is the low, young surface under either field.

`fitted_surface_algorithm_choice` picks the surface:

| value | for | why |
|---|---|---|
| `cone` (default) | a fan or debris cone | z falls linearly with distance from an apex |
| `plane_ransac` | a piedmont or broadly planar tread | RANSAC, not OLS: incised distributaries and inset lobes are large-amplitude outliers covering a real fraction of the domain, and OLS would tilt the reference toward the very features the residual is meant to isolate |
| `poly2` | a broadly curved surface with no single centre | |

**Use `cone` for a fan**, because a fan is radial and a cone is not a low-order
polynomial in x and y. Measured on a synthetic fan (2 km radius, 0.03 gradient,
lobes inset 4 m), the residual spread of the older surface — which should be flat,
so its spread is pure misfit — is about 34 m under a RANSAC plane, 10 m under a
poly2, and 0.3 m under a cone, against a 4 m signal. The apex is taken as the
centroid of the highest 0.1 % of fit-mask cells unless you set
`fitted_surface_apex_x` / `_apex_y`; it is supplied rather than solved for because
it is a strongly non-linear parameter and you can usually point at a fan head more
reliably than a solver can find it.

### Band first, then label

This is the load-bearing change, and the reason the tool exists.

The obvious ordering — the one the legacy code used — labels connected components
over the whole flat mask and consults height afterwards:

```
labels = label(flat mask minus breaks)     # one pass, whole raster
band   = f(HAND, labels)                   # height consulted after
```

That ordering **forces the riser network between terraces to be continuous**. A
road, a gully or a stock track cutting through a riser lets two treads at different
heights merge into a single component, and nothing downstream can separate them
again. Four independent repair mechanisms had accumulated against that one failure:
a break-of-slope skeleton with A\* endpoint bridging to re-close the gaps, a
residual-relief re-break pass, a bimodal-HAND leak splitter, and a HAND banding pass
that split wide labels into sub-labels after the fact.

lsdtt3 inverts the order:

```
band   = ladder(reference field, admitted)   # height consulted FIRST
for b in bands:
    labels_b = label(admitted AND band == b) # one pass per band
```

Two pixels can only share a label if they share a band, and two treads at different
heights are by construction in different bands. **The riser never has to be found
and never has to be continuous** — the labelling passes for band `b` and band `b+1`
never see each other's pixels, so there is nothing for the gap to leak through.
The cost is one labelling pass per band instead of one in total, which is negligible
next to computing HAND.

A second consequence is worth knowing if you are comparing outputs with the legacy
tool: every patch is now a genuine connected component within one band. The old
`(label, band)` accumulation key could aggregate spatially disjoint patches into one
record and report a single median height for them.

Labels are a **global running counter 1..N**, assigned band by band; the band is
*not* packed into the label id. It travels in `_band_index.tif` and in
`_surface_patches.csv`, so nothing has to divide a label to learn its level. (There
is one hard limit: the raster backing store is float32, so a label above
2^24 = 16 777 216 could not round-trip exactly and every join keyed on it would
silently corrupt. The engine throws rather than crossing it. No realistic run comes
near.)

### There is no valley mode and no terrace mode

Finding a valley floor and segmenting a terrace flight are the same computation at
different band counts. The valley floor is the lowest band of a two-band cut; a
terrace flight is an n-band ladder on the same population. So there is one
parameter, `ladder_level_count`: `2` for the valley cut, `n` for an n-level ladder,
`0` (the default) for automatic.

This was already true in the legacy numerics without being true in the interface —
the valley engine and the terrace classifier called the same antimode break finder,
one taking the first qualifying dip and the other keeping all of them.

### Finding the breaks

`ladder_break_algorithm_choice` selects the break finder. The three are genuinely
different tools and the choice matters more than most parameters on this page.

**`antimode`** (default) finds dips in a smoothed log10 density. It is the incumbent,
shared with the valley engine, and it works where populations genuinely separate. It
has one hard constraint: **it bins on log10 and therefore cannot band a field with
negative values.** Every negative residual is floored into the lowest bin, so on a
fitted-surface field the entire below-surface population — the inset lobes, the
thing you are looking for — collapses into a single spike. Selecting
`fitted_surface` switches the finder to `jenks` automatically if you did not name
one, and rejects an explicit `antimode` with that explanation rather than silently
producing meaningless breaks.

**`jenks`** is Fisher–Jenks on the linear field. It takes the level count directly,
which makes it the natural choice when the count is explicit, and it is the source
of the reported goodness-of-variance fit. A terrace flight often wants `jenks`: a
flight is spaced roughly linearly in height above the channel, and log10 binning
compresses the upper treads until they merge.

**`gap`** looks for the **empty space between clusters** and puts the break in it.
This is what a reader's eye does on a HAND histogram — the active floodplain is the
cluster near zero, it stops somewhere, the next surface starts higher up, and the
break belongs in the space between. Neither of the other two looks for that space.
Jenks minimises within-class variance, so on a distribution dominated by one huge
upland mode it puts the cut *inside* that mode and misses the floodplain entirely;
antimode looks for a dip but on an axis that compresses everything above a few
metres into a handful of bins. The finder bins linearly and offers two kinds of
candidate:

- **empty runs** — contiguous bins each holding less than `ladder_gap_empty_fraction`
  of the tallest bin and spanning at least `ladder_gap_min_width_metres`. The clean
  case: nothing at all between two levels.
- **troughs** — interior local minima of the smoothed histogram. These are needed
  because "empty" is measured against the *tallest* bin, and in a mountain catchment
  the floodplain bin is enormous: a real dip between the floor and the first terrace
  can hold 5–20 % of that peak and never register as empty, while the sparse tail
  150 m up is empty everywhere. Without troughs the finder has no candidate anywhere
  near the valley floor and takes a break out in the tail instead — which bands the
  entire terrace staircase as valley floor.

`ladder_gap_rank` decides the order candidates are offered in. Ranking by width
prefers the widest empty span, which in a landscape with a long sparse tail is
always high up; ranking `lowest` first (the default) walks the ladder from the
bottom, which is the order the levels physically sit in. `ladder_gap_break_position`
places the break within the span: the span is nearly empty by construction so any
position separates the same two clusters, but the sparse pixels *inside* it — the
bank, the riser, the floodplain margin — go to whichever side the break leaves them
on. Use `upper` when the low band is the valley floor and you want the whole margin
in it, `lower` for the strictest reading of the active floodplain.

### How many levels, and when to believe them

This is the part with no counterpart in the 2017 or 2022 methods, and it exists to
stop the tool reporting levels the DEM cannot support.

A level count is **admissible** when the tightest gap between consecutive band
centres clears

```
kappa * E_v
```

where `E_v` is the vertical accuracy of the elevation product
(`vertical_accuracy_metres`) and `kappa` is `ladder_resolvability_kappa`, default 3.
The answer is the **largest admissible count**. Two levels cannot be claimed as
distinct if their centres sit closer together than the DEM can measure.

Three details in that rule are deliberate:

**Band centres are medians, not midpoints.** `m_j` is the median of the j-th band's
actual reference-field values, not the arithmetic midpoint `(b_{j-1} + b_j)/2` of
its break edges. The midpoint is purely geometric and depends only on where the
breaks fell; the median is where the surface *actually sits*. Since what a terrace
correlation cares about is the elevation of the tread, not the elevation halfway
between two break levels, the median is the honest input to a "can this DEM tell
these two surfaces apart" test.

**Every count is scanned; the search does not stop at the first hit.**
Admissibility is *not* monotone in the level count, precisely because the centres
are medians of real populations: merging two well-separated bands can produce a
merged band whose median lands close to a neighbour's, so a coarser partition can
have a *smaller* minimum gap than a finer one. With the count bounded at
`ladder_max_level_count` (default 12) the full scan is trivially cheap.

**It is a constraint with an extremum over the feasible set, never a trade-off.**
Do not be tempted to select the count by maximising an additive score such as
`w1 * GVF + w2 * gap`. Goodness-of-variance-fit is bounded on [0, 1] while the band
gap scales as `range / k` and is unbounded, so over a large vertical range the gap
term dominates by orders of magnitude and the criterion collapses to the smallest
count regardless of what the data say. GVF is reported for diagnosis and is never an
input to the selection.

`vertical_accuracy_metres` is **metadata you have to supply** — it cannot be derived
from the DEM. Leaving it at 0 disables the floor, and then every window reports
`resolved` and the ladder carries no instrument guarantee at all; the tool warns when
you do. (A DEM whose vertical units are feet on a metre grid has an effective
accuracy 3.28 times larger than it looks, and nothing in the file will tell you.)

`vertical_accuracy_algorithm_choice: slope_dependent` reads a gentle value below a
break gradient (`vertical_accuracy_break_gradient`, m/m) and a steep value above it, because photogrammetric and InSAR DEMs degrade
on slopes — TanDEM-X relative LE90 is roughly 2 m on gentle ground and 4 m on steep,
and one number for a raster spanning a valley floor and a steep valley wall is wrong in
both places. `E_v` is resolved **per along-channel window, not per pixel**: the floor
constrains the separation between band *centres*, and a band centre is a property of
a partition rather than of a place, so it has no slope of its own. Each window takes
the median gradient of its own admitted cells; the value it used is written to
`_ladder_windows.csv`.

**Where no count qualifies, the window reports `unresolved` rather than returning a
number.** That is not "one level is the answer"; it is "this DEM cannot resolve a
ladder here". An explicit `ladder_level_count` is honoured and reported honestly — if
that partition is not admissible the window records the count you asked for with
status `unresolved`, rather than silently substituting a different one. What you may
claim about your levels depends on which branch was taken, and that is not visible
in any raster, so **read `<stem>_ladder_windows.csv` before describing your levels**.

### Reach-relative windowing

Both legacy engines windowed their break search along the channel and auto-sized the
window as `max(2 km, channel length / 8)`; lsdtt3 keeps this as
`ladder_reach_window_metres` (0 = auto, negative = one global ladder). It matters
because a terrace flight is only a step function *relative to a reach*: over tens of
kilometres of channel the absolute height of a given tread changes, and a single
global set of breaks smears the flight.

Windowing is only meaningful for `hand`. **A fan has no along-channel coordinate**,
so for `fitted_surface` the ladder is global over the fit domain and the window
parameter is ignored with a note on stdout.

There is one composability cost. Reach-relative windowing needs each pixel's HAND
*donor node*, and a raster cannot carry those, so a HAND supplied through
`hand_fname` is always banded globally. Since HAND is the expensive stage and the
level count is what you iterate on, the practical recipe is: compute HAND once with
`write_hand_raster: true`, iterate on `ladder_level_count` against `hand_fname`,
then do the final run computing HAND if you want the windowed ladder. The
supplied-HAND path produces bit-identical bands and labels to the computed path,
which is pinned by a test.

### Classification, and where it sits relative to Clubb 2017

For the `hand` field only, patches are classified as floodplain or terrace and
written to `_surface_class.tif` and `_surface_class.csv`. `max_floodplain_hand_metres`
and `min_terrace_hand_metres` are the direct descendants of the 2017 paper's single
"threshold height above the modern river channel", but they are **demoted from
mechanism to optional override**: under band-then-label the ladder has already
decided where one level ends and the next begins, and these are hard caps for cases
where you know the floodplain height independently.

No class is written for `fitted_surface`. Floodplain and terrace are statements about
height above a channel; they mean nothing on a residual, so none is invented.

### The valley mask

`<stem>_valley_mask.tif` is band 0, cleaned. `clean_valley_mask` erodes, drops
patches that do not reach the channel, then regrows — note that the regrow fills
enclosed holes, so **cleaning can make the mask bigger**; it is not a filter.
`fill_valley_mask_holes` (default true) fills enclosed holes explicitly, which is
the counterpart of the OpenCV flood fill Clubb et al. (2022) added for
anthropogenically modified valley floors. It matters more than it sounds: a hole in
the valley floor is a hole in the *domain the centreline is traced through*, and the
routed path goes around it, hugging one wall instead of running down the middle.

`valley_mask_closing_width_fraction` sets the closing radius as a fraction of the
*measured* valley width rather than in metres, because a radius in metres does not
travel between DEMs: a value that seals the slivers in a 400 m floodplain erases a
100 m one. The width is measured by casting perpendicular sections and taking the
median of the cumulative length of valley cells along each — cumulative, so a bar or
an island does not count as floor. About 0.2 is a reasonable starting point.

`restrict_valley_to_line_of_sight` is a strong filter written for separating one
valley floor from a neighbour across a broad col, where connectivity alone merges
them. It defaults to **false**, because on a meandering floodplain a great deal of
genuine valley floor has no straight line to the channel that stays inside the mask:
on a 10 m DEM it removed a large share of band 0 on its own with every other
filter disabled, making the valley mask smaller than the floodplain class computed
from the same field.

Only the valley mask is affected by any of this. The band index, the labels and the
patch table keep their holes and their speckle, because there a rise in the valley
floor is a real observation.

---

### The clubb2017 method

Clubb et al. (2017) classify floodplain and terrace from two per-pixel quantities, the gradient of the ground and the height above the nearest channel. Both bounds are read from the data, not set by hand, using quantile-quantile (QQ) plots. `clubb2017` implements that method with HAND as the height measure. The steps, in order:

1. **Bounds from QQ plots (`floodplain_threshold_algorithm_choice: qq_gaussian`).** For each quantity, the data quantiles are compared with a Gaussian drawn through the 25th and 75th percentiles (`floodplain_qq_lower_percentile` and `floodplain_qq_upper_percentile`, the values Clubb et al. used). The bound is the lowest value at which the data rejoin that line and stay on it, to within `floodplain_qq_threshold_fraction` (0.01) of the data range. The chosen values are printed and written, with the QQ tables they came from (`<stem>_floodplain_hand_qq.csv` and `<stem>_floodplain_gradient_qq.csv`, from `write_floodplain_threshold_csv`).
2. **Gradient bound or not (`floodplain_bound_algorithm_choice`).** `hand_and_gradient` applies both bounds, as in the paper. `hand_only` applies the HAND bound alone. `auto` (the default) measures the fraction of the DEM's floodplain-gradient cells below `low_relief_max_gradient` (0.1 m/m). When that fraction is more than `low_relief_min_area_fraction` (0.8), `auto` uses `hand_only`; otherwise it uses `hand_and_gradient`. The switch is for low-relief DEMs, where most ground already passes a gradient test. The measured fraction and the branch taken are written to the method record.
3. **Admission.** A cell is admitted when its HAND is below the HAND bound (`max_floodplain_hand_metres`, once chosen) and, when the gradient test is in use, its gradient is below the gradient bound. The gradient is a polynomial-fit gradient over a window of radius `floodplain_gradient_window_radius_metres` (0 means three cells). See [Choosing a method](#choosing-a-method) for how to set that window.
4. **Floodplain.** Admitted cells form 8-connected components. Components smaller than `min_component_pixels` are dropped. A component is kept only if it contains a channel pixel of the in-scope channel. That channel contact is what makes the patch a floodplain rather than any low, flat patch. Holes are filled (`fill_floodplain_holes`, which is on under `clubb2017`).
5. **Terraces.** Terraces are admitted ground that is not floodplain, passes the gradient test (when in use), and has a terrace field above the HAND bound: `terrace_field > max_floodplain_hand_metres + min_terrace_height_above_floodplain_metres`. A `max_terrace_height_metres` above zero also caps it. Terraces need no channel contact. Under `clubb2017` the terrace field is `swath_baseline`: the relief above a baseline taken from the channel. Terraces further than `swath_half_width_metres` (500 m by default) from the channel are not measured. The terrace field can be written as `<stem>_terrace_field.tif` (`write_terrace_field_raster`).
6. **Optional second pass (`floodplain_second_pass_algorithm_choice: qq_upper_departure`, default `none`).** This is not part of Clubb et al. (2017). The QQ plot of HAND is taken over the first-pass floodplain cells, against the same Gaussian line, and the highest quantile at which the data still agree with it is the second bound. The floodplain is rebuilt with that bound, with the same gradient bound, channel contact and size filter, and its holes are always filled. Terraces are then measured from the second-pass bound. The first pass is kept (`<stem>_floodplain_first_pass_mask.tif` and the first-pass polygons), and the second QQ table is `<stem>_floodplain_hand_qq_second_pass.csv`. The option is aimed at wide flats that stand above the valley floor, such as fans at tributary mouths.

#### Differences from the LSDTopoTools2 legacy code

The method descends from the LSDTopoTools2 floodplain and terrace code (LSDFloodplain and LSDTerrace). lsdtt3 departs from it in two deliberate places:

- **Terrace window.** The legacy terrace window could overlap the floodplain whenever its minimum terrace height was below the floodplain bound. In lsdtt3 terraces start above the floodplain bound, so no cell is both.
- **Channel seeds.** Only the pixels of the in-scope channel seed the floodplain (the channel selected by `channel_definition_algorithm_choice`), not every channel pixel of any Strahler order.

The legacy `minimum_terrace_height` is `min_terrace_height_above_floodplain_metres` here, and the second pass (`qq_upper_departure`) is new. This page does not catalogue every other legacy behaviour.

### The ladder: floodplain, unresolved windows and clamping

The ladder (`level_algorithm_choice: ladder`) finds its breaks in the data, as described under "Finding the breaks" and "Reach-relative windowing". Three details matter when reading its output.

- **Which level is the floodplain.** The floodplain is the lowest level that touches an in-scope channel, found in each reach window (`ladder_reach_window_metres`). Terraces are the levels above it.
- **Unresolved windows.** A window where the ladder finds no admissible levels is unresolved. `ladder_unresolved_window_algorithm_choice` decides what it bands with. `nearest_resolved` (the default) borrows the breaks of the nearest resolved window on the same channel, so band numbers mean a comparable HAND range along the channel. `single_band` bands the whole window as one level, so the floodplain there can span a wide HAND range. In either case the window is still reported as unresolved in `<stem>_ladder_windows.csv`, and `breaks_from_window_m` names the window that lent its breaks. Read that file first.
- **Negative HAND.** A cell whose flow reaches a channel cell higher than itself has a negative HAND. This comes from fill and flat-resolution artefacts. Those cells are set to 0 before banding and counted (`negative_hand_cells_set_to_0` in the method record). The antimode break finder works on the logarithm of the field, so it needs non-negative values.
- **Patch clustering check.** `write_patch_level_clusters_csv` (on by default) tests whether the ladder's floodplain and terrace patches cluster in HAND. Each classified patch is one value, its median HAND, and the sorted values are split into clusters wherever neighbours are at least the required gap. The result is `<stem>_patch_level_clusters.csv`.

### Training polygons: supervised thresholds

`training_surfaces_fname` takes polygons (GeoPackage, FlatGeobuf or shapefile) of floodplain and terrace that you have drawn where you are sure of the surface, usually on a hillshade. The cells inside them are used to choose thresholds. `training_surfaces_layer_name` and `training_surfaces_class_field_name` (default `class`) select the layer and the class field.

- **Classes.** `floodplain` or `fp` is floodplain. `terrace`, `terraceN`, `tN` or `terrace N` is terrace. Case and spaces are ignored, and the terrace number is ignored for the thresholds: every terrace polygon is simply terrace.
- **Summary of each polygon.** Each polygon is summarised by its median field value and its core, the cells between its own 25th and 75th percentiles. A few cells off the surface, such as a gully crossed while drawing, therefore do not move it. The 10th and 90th percentiles are reported too.
- **Flagging.** A polygon whose interquartile range exceeds both `training_spread_factor` (3.0) times the required gap and `training_max_relative_spread` (0.25) times its median is flagged as a sloping surface. It still counts as terrace for the floodplain and terrace bounds, but takes no part in the level breaks. The relative limit keeps high terraces, whose heights vary more in metres, as single levels.
- **Thresholds.** Each threshold is the value that best separates the trimmed cores of the two populations it divides, by Youden's J (balanced accuracy). On a plateau it is the midpoint of the optimal interval. A population with fewer than `min_training_cells_per_class` (20) core cells in a window does not take part in that window. The window then falls back to the fit over all training cells, and if that is also too small, to the data-chosen threshold.

Hand-drawn polygons do not separate HAND cleanly, so treat the chosen thresholds as something to check against the polygons. The cells used for training are the same cells a score against those polygons would use, so the two are not independent.

### Mapped floodplains and terraces

`mapped_surfaces_fname` takes mapped floodplain and terrace polygons (GeoPackage, FlatGeobuf or shapefile) to score each method against. `mapped_surfaces_layer_name` and `mapped_surfaces_class_field_name` (default `class`) select the layer and field. Class values follow the training rules above. Any other value stops the run before any work. Terrace levels are not scored separately: every terrace class is one class.

The polygons are reprojected to the DEM's CRS and burned by cell centre only, with no growth to the edges. `mapped_surfaces_extent_fname` is a polygon of the area the mapper actually surveyed. A cell is scored only if it is inside that extent, the DEM is valid there, and the method used HAND.

Outputs:

- `<write_fname>_mapped_surface_class.tif`: 0 for surveyed cells that are neither class, 1 for floodplain, 2 for terrace.
- `<write_fname>_mapped_surface_scores.csv`: for each method, the 3 × 3 confusion matrix (mapped class against method class, in cells), and floodplain and terrace precision, recall and intersection-over-union. Precision is the correctly classed cells divided by the cells the method put in that class. Recall is the same count divided by the cells mapped as that class. IoU is the correct count divided by the cells in either.

### Method records and comparison outputs

`write_method_record_csv` (on by default) writes `<stem>_method_record.csv`, with the columns method, item, value and source. It lists every stage, the source of each value (the preset, the parameter file, or a measurement), each data-chosen bound, and the counts. The record is the place to look when you need to say which parameters produced a result.

When `fluvial_surface_algorithm_choice` lists several presets, each runs and writes its own outputs, and the comparison outputs are added:

- `<write_fname>_method_summary.csv`: one summary per method.
- `<write_fname>_method_comparison.csv`: for each pair of methods, the floodplain and terrace intersection-over-union, the agreement on cells that either method classed, and the 3 × 3 confusion in cells and in m².
- `<write_fname>_method_agreement_floodplain.tif` and `<write_fname>_method_agreement_terrace.tif`: UInt8 rasters, where bit k is set if method k (in list order) gave that class. At most eight methods can be compared.

### Relative elevation as the reference field

`reference_field_algorithm_choice: relative_elevation` bands the relative elevation (RER) instead of HAND. Set `relative_elevation_floor_mask_fname` to the valley-floor mask the centrelines were traced through (cells above 0.5 count). RER is computed over the whole DEM. Each cell's donor is the channel pixel of the centreline node that owns it, so the reach windows, training and bounds work as they do for HAND. How RER is built is described under `lsdtt-valley-centreline` below.

---

## `lsdtt-valley-centreline`: the trough, the axis and the width

This tool retains what is genuinely its own: the trough carve, the routed
centreline, the pinned monotonic thalweg baseline, cross sections and
cross-section hydraulics. It takes the valley as given.

### The trough carve

This follows Clubb et al. (2022) closely, and the parameter names transfer from the
legacy driver:

1. Punch the valley out of the DEM as nodata, dilating the hole by
   `trough_nodata_buffer_pixels` (default 2, as in the legacy driver) before bank
   elevations are sampled, so pixels just inside the valley edge are not read as
   bank.
2. For every valley pixel, take the distance to the nearest wall pixel and that
   wall pixel's elevation.
3. Take the minimum bank elevation in a circular window of radius
   `min_bank_elevation_window_radius_metres`.
4. Subtract `trough_scaling_factor * distance_to_wall`, giving the V-shaped trough of
   Eq. 2 in the paper.
5. Breach and fill, and repeat for `centreline_loops` iterations — default **5**,
   the value the paper settled on as the balance between an accurate centreline and
   computation time.
6. Build a flow object on the trough and route D8 from the channel head to the
   outlet. That path is the centreline.

Two departures. `trough_scaling_factor` defaults to **0.1** (matching the legacy
driver's default, rather than the α = 0.5 quoted for the paper's own analysis), and
it is **doubled automatically on routing failure**, up to two retries — a shallow
trough is preferable where it works and there is no reason to make the user find the
boundary by hand. `trough_along_channel_slope` (default 0.02 m/m) adds a synthetic
downstream tilt that has no counterpart in the paper; it stops a near-flat or
through valley from routing its centreline out of the wrong, upstream tile edge, a
failure the papers' field sites did not encounter.

The reason for carving a trough rather than skeletonising the mask is the paper's:
the routed path carries flow-routing information, which a morphological medial axis
does not, so everything downstream — flow distance, ordering, the monotonic
baseline — is available on the centreline nodes.

### Cross sections and width

Width is measured on sections cast across the centreline, with
`cross_section_spacing_metres` (default 500 m) between stations. The local direction
is smoothed over `cross_section_tangent_window_nodes` centreline nodes on each side
(default 8), which is the same construction as the paper's node spacing `n` and its
two-bearing average `b = 0.5(b_u - b_d) + b_d`, generalised from two points to a
smoothed tangent. Each section reports `valley_width_metres` as the span of the
contiguous in-valley run through the centreline; a station where a side is not found
before the section length runs out is marked not-found and excluded from the
per-channel statistic rather than contributing a truncated width.

`curved_cross_sections` is an extension beyond the 2022 method. Straight orthogonal
rays are correct in straight reaches and cross each other in tight bends, which
double-counts floor on the inside of the bend and truncates it on the outside.
Curved sections are traced as isocontours of a smooth along-valley coordinate field
— straight where the valley is straight, curving gradually through bends like the
grid lines of a boundary-fitted mesh — with sides clipped at the Voronoi bisector of
the nearest kept neighbours and a greedy decimation pass that drops stations in tight
bends so the survivors stay long.

`write_cross_section_hydraulics` goes further than either paper, deriving stage–width
and stage–area curves per section and locating the breakpoints in hydraulic depth
that mark bankfull and flat-surface boundaries.

### Refinement

`refine_valley_with_baseline` re-tests the supplied mask against the fitted
centreline baseline and removes wedges referenced to a stale part of the channel.
It requires an explicit `refine_max_height_metres`, because this tool no longer
derives a break of its own and a refinement run against a fabricated threshold would
look like a result and be one only by accident. The refined mask is written for you
and is **never fed back** to `lsdtt-fluvial-surfaces`. The dataflow between the three
tools runs one way, with one exception: the cross-section levels file
(`write_cross_section_levels_csv`) is read back by `lsdtt-fluvial-surfaces` as `supplied_levels_fname`.

---

## `lsdtt-valley-centreline`: relative elevation, terrace profiles and widths

### Relative elevation (RER)

`write_relative_elevation_raster` writes `<stem>_relative_elevation.tif`. Each cell's value is its elevation above the baseline of the centreline node that owns it. `write_relative_elevation_owner_raster` writes which centreline owns each cell. The raster differs from HAND in one respect: it is measured against the valley's own centreline baseline, not against the channel the cell's flow happens to reach. A meandering or wall-hugging channel therefore does not imprint itself on the surface.

**Ownership.** Each cell belongs to one centreline node, by one of two rules:

- On the valley floor (inside the floor mask), the larger channel wins. Each node claims a disc whose radius is its distance to the valley wall, and where discs overlap, the node with the larger drainage area takes the cell.
- Off the floor, the cell belongs to the nearest node along the ground, with the ground distance scaled by that node's floor half-width.

`relative_elevation_max_distance_metres` (0 means no limit) is applied after ownership is decided. It cuts cells that have an owner but are too far from it along the ground.

**Domain.** `relative_elevation_domain_algorithm_choice` sets where RER is computed: `valley_mask` (the default), `supplied` (a raster given in `relative_elevation_domain_fname`), or `whole_dem`.

**Baseline.** `centreline_baseline_algorithm_choice` sets how each node's baseline is pinned to the channel:

- `projected_monotonic` (the default): each channel pixel pins the nearest centreline node, with the lowest elevation winning where several land on one node. The baseline is then made monotonic downstream and is never above any pin.
- `channel_crossings`: the baseline is pinned where the centreline crosses the channel. A centreline node that is itself a channel pixel also pins.

**Mask clean-up.** The valley mask is cleaned before RER with `valley_mask_closing_radius_metres`, `valley_mask_opening_radius_metres` (which removes narrow tendrils) and `fill_valley_mask_holes`. These are the same settings as in `lsdtt-fluvial-surfaces`.

### Terrace profile

`write_terrace_profile_csv` profiles the labelled surface patches along the valley centreline they belong to. It takes two inputs from `lsdtt-fluvial-surfaces`: `surface_patch_labels_fname` (the `<stem>_banded_labels.tif`) and `surface_patches_fname` (the `<stem>_surface_patches.csv`, which gives each label's band and class). Each labelled pixel is tied to a centreline node by the RER ownership rule. Each patch is profiled against the centreline that owns most of its pixels. Its pixels are then binned every `terrace_profile_bin_metres` (100 m by default) along that centreline.

Outputs:

- `<stem>_terrace_profile.csv`: one row per patch per bin, with the centreline, patch label, band, class, bin centre, pixel count, and the minimum, 25th percentile, median, 75th percentile and maximum of both the elevation and the relative elevation.
- `<stem>_terrace_patches.csv`: one row per patch.

### Cross-section levels

`write_cross_section_levels_csv` clusters each channel's cross-section hydraulic breaks into levels and writes `<stem>_cross_section_levels.csv`. Its columns include `source_id`, `level_hand_m`, `n_sections`, `spread_m` and `status`. The breaks of all sections on one channel are pooled and split by the ladder's own rule. A level is kept only if it appears on at least `cross_section_level_min_section_fraction` (0.5) of that channel's sections. The file is written only when `write_cross_section_features` is on, which is its default.

The file is read back by `lsdtt-fluvial-surfaces` through `supplied_levels_fname`, with `level_algorithm_choice: supplied_levels` or the `cross_section_levels` preset.

### Valley width network

`run_valley_width_network` (off by default) measures valley widths on every channel, from the channel threshold, instead of on supplied source points. It requires `source_points_fname` to be empty. Channel segments are defined with `channel_definition_algorithm_choice` (`contributing_pixels` by default). Two rules keep tributaries from measuring the parent valley:

- `run_tributary_stop_at_parent` (on by default): a tributary's sections stop where they reach the parent's valley. `tributary_stop_width_fraction` (1, meaning the parent's valley wall) sets where that is.
- `run_tributary_crossing_removal` (on by default): parent sections that cross a tributary's kept sections are excluded.

Outputs: `<stem>_valley_segments.csv`, `<stem>_valley_segment_channels.fgb` (each segment's D8 channel), and `<stem>_valley_widths.csv`, with one row per station. Each station carries the drainage area at its channel pixel, which never decreases downstream. The cross-section outputs gain `segment_id`, `parent_segment_id`, `nesting_depth` and `excluded_reason`.

---

## `lsdtt-channel-extraction`: perimeter sources

`run_perimeter_sources` (off by default) finds the channel sources where valley floors enter a DEM across its edge. It is for a DEM clipped around a few large valleys, where the valleys run in from outside the data. It needs open boundaries on every edge (`boundary_conditions: oooo`, the default). The method works from the DEM alone:

1. The outer edge of the largest valid block is traced into a closed perimeter profile.
2. A gradient is computed at each perimeter cell from a least-squares plane over `perimeter_gradient_radius_metres`. Cells below `perimeter_max_gradient` (0.05 m/m) are gentle. Steep gaps shorter than `perimeter_gap_closing_metres` inside a gentle run are absorbed, which handles channel banks.
3. A gentle run is a valley crossing if it is at least `perimeter_min_valley_width_metres` long and is a low spot: the perimeter rises at least `perimeter_min_relief_metres` (5 m) above the run's lowest cell on both sides.
4. Outlets are D8 outlets on the edge. The outlet cells of one low-spot run form one outlet run, with their drainage areas summed. The `n_perimeter_outlets` runs with the largest summed area are kept (one by default).
5. A crossing is a source if a floor cell drains into the DEM, its flow path stays below `perimeter_inward_max_gradient` (0.1 m/m) for `perimeter_inward_check_distance_metres` (500 m), and the path ends in a kept outlet run.

Negative values select the built-in defaults for the gradient radius (the larger of 5 cells and 10 m) and the inward gradient window (the larger of 5 cells and 25 m).

Output: `<write_fname>_perimeter_profile.csv` (`write_perimeter_profile`, on by default), with one row per trace step. It records each step's position, elevation, gradient, whether it is gentle, and its run and outlet IDs.

Guidance:

- `perimeter_min_valley_width_metres` defaults to 300 m, which is set for 30 m global DEMs, where it is 10 cells. Lower it for a finer DEM. On lidar, set it in metres to the narrowest valley you want to count as a crossing, and check that the result is not picking up channel banks.
- To feed the sources into `lsdtt-valley-centreline`, set `auto_detect_outlet_points: false`. In testing, automatic classification turned one of the perimeter sources into an outlet.

---

## `lsdtt-feature-finder`: riser finding

`feature_type: find_risers` finds short, steep steps between gentler surfaces, such as terrace risers, shorelines and scarps. Each riser is returned as an open line segment with a step height and a face gradient. The steps are:

1. **Edge field.** `riser_edge_metric` chooses the field the lines are drawn on: `vertical_excess_curvature` (the default), `curvedness`, `maximum_curvature` or `multidirectional_hillshade`. For the polyfit metrics, `riser_polyfit_lambda_metres` (100 m) sets the half-wavelength.
2. **Adaptive threshold.** `riser_vec_adaptive_method: mean_std` (the default) marks a cell as riser when its edge value is more than `riser_vec_adaptive_std_factor` (0.5) local standard deviations above the local mean. The standard deviation is floored at 1/255 of the field's range, so the test stays defined on flat neighbourhoods. The `mean` and `gaussian` methods are also available. Their offset `riser_vec_adaptive_C` is in grey levels.
3. **Clean-up and skeleton.** The mask is closed with `riser_closing_radius_metres` and small components are removed. The mask is then thinned to a skeleton, and spurs and terminal loops up to `riser_max_spur_pixels` and `riser_max_terminal_loop_pixels` (both 10) are pruned.
4. **Gap closing (optional).** `riser_gap_closing_algorithm_choice: astar` bridges broken risers. Bridges are limited by `riser_connect_search_radius_metres` (200 m), `riser_connect_elevation_tolerance_metres` (2 m) and a corridor half-width.
5. **Step height.** The step across each segment is measured over `riser_step_window_radius_metres`. Segments with a median step below `riser_min_step_height_metres` (1 m) or fewer than `riser_min_segment_pixels` (5) pixels are dropped.

`riser_search_mask_fname` restricts where risers are looked for. It is empty by default, which means every valid cell.

**Why the threshold changed.** The earlier threshold rescaled the edge field to 8-bit grey levels between its minimum and maximum, and worked in those levels. The offset was rounded to whole grey levels, so an offset between -1 and 0 had no effect, and the same parameter behaved differently as the field's range changed. `mean_std` works in the units of the field and does not depend on its range.

**Pixel-counted parameters.** Several parameters count cells, not metres: `riser_vec_adaptive_block_size_pixels`, `riser_max_spur_pixels`, `riser_max_terminal_loop_pixels`, `riser_min_component_pixels`, `riser_binary_min_component_pixels`, `riser_min_segment_pixels` and `riser_connect_corridor_half_width_pixels`. Their meaning in metres is the count times the cell size, so they must be scaled when the grid spacing changes. A spur length of 10 cells is 10 m on a 1 m DEM and 100 m on a 10 m DEM. Choose them from the riser width and length in metres, divided by the cell size.

**Scoring against mapped risers.** `mapped_risers_fname` takes mapped terrace edges or risers (lines). `mapped_risers_extent_fname` gives the surveyed area, and `mapped_risers_tolerance_metres` (20 m) sets how close a detected cell must be to a mapped edge to count. The outputs are `<write_fname>_mapped_risers.tif` and `<write_fname>_mapped_riser_scores.csv`, which give per scope the cell counts, precision, recall and F1. The `extent` scope covers every valid cell inside the extent. The `search_mask` scope covers those cells inside the search mask too.

---

## Choosing a method

This section describes what each method does and what it is sensitive to. It does not rank the methods, because the right choice depends on the landscape and the data. Use the method record and the comparison outputs to see what each run did, and report the numbers together with the parameters that produced them.

- **Start with the ladder if you have no terrace heights to start from.** The ladder finds its breaks in the data. Check `<stem>_ladder_windows.csv` first, because it shows which reach windows were resolved. Where many windows are unresolved, the levels are borrowed from other windows or left as one band, so the classes are less specific there.
- **Use `clubb2017` for the published rule.** It returns a floodplain class and a single terrace class above the floodplain bound, with both bounds chosen from QQ plots. Its floodplain needs channel contact, so flat ground that is not connected to the channel in the admitted mask is not classed as floodplain. Check the measured gradient-bound fraction in the method record before deciding whether the gradient test is doing anything.
- **Scale the gradient window to the grid.** `floodplain_gradient_window_radius_metres` is a radius in metres, and 0 means three cells. A window that is wide compared with the valley floor averages the valley walls into the floor. The gradient bound then cuts the floor at a fixed distance from the walls rather than at the edge of the floor. Set the radius from the valley floor width and the grid spacing, and check that the floodplain edge follows the floor.
- **Expect the HAND reference to differ from RER where treads drain sideways.** HAND measures each cell against the channel its flow reaches. Where a tread drains into a different reach, HAND can give it the same height as the tread below. RER measures against the baseline of the centreline that owns the cell, so it separates such treads, but it needs floor centrelines and can leave more windows unresolved. The synthetic terrace staircase in `tests/cpp/test_fluvial_surfaces.cpp` has this structure.
- **Use `fitted_surface` for fans**, which do not form steps in HAND.
- **Scale the riser parameters to the grid.** The pixel-counted riser parameters need the cell size taken into account, as described above, on fine DEMs.
- **Scale the perimeter width to the grid.** The 300 m default is for 30 m DEMs.
- **Keep scoring polygons separate from training polygons.** A score against the polygons that trained the thresholds is not an independent test.

---

## From `lsdtt-valley-metrics`

A rough map from the LSDTopoTools2 driver to the three lsdtt3 programs. Where a
parameter has no successor, the reason is in the right-hand column.

| `lsdtt-valley-metrics` | goes to | lsdtt3 |
|---|---|---|
| `surface_fitting_radius` | 1 | `admission_gradient_window_radius_metres` |
| `slope_threshold`, `use_absolute_thresholds` | 1 | `max_admission_gradient` with `admission_algorithm_choice: gradient` |
| `QQ_threshold`, `*_lower_percentile`, `*_upper_percentile` | 1 | `relief_threshold_algorithm_choice: qq_gaussian` with `relief_qq_lower_percentile`, `relief_qq_upper_percentile` and `relief_qq_threshold_fraction` (the Clubb et al. 2017 rule). In `lsdtt-fluvial-surfaces` the floodplain bounds use `floodplain_threshold_algorithm_choice: qq_gaussian`. `cdf_divergence` is the alternative, with `relief_threshold_multiplier` |
| `minimum_patch_size` | 1, 2 | `min_component_pixels` in tool 1; `min_band_component_pixels` in tool 2, where it is **per band**, not global |
| `threshold_SO` | 2 | `channel_definition_algorithm_choice: stream_order` with `min_channel_stream_order` (a Strahler order), or the channel-scope model: `source_points_fname` / `outlet_points_fname`, or `threshold_contributing_pixels`. Check that the Strahler order matches the legacy `threshold_SO` on your data |
| `relief_threshold` | 2 | replaced by the ladder; `max_floodplain_hand_metres` is the optional hard cap |
| `minimum_terrace_height` | 2 | `min_terrace_height_above_floodplain_metres`, the height above the floodplain bound at which terraces start. `min_terrace_hand_metres` is the separate, absolute lower bound of the ladder |
| `get_terraces` | 2 | no switch — `ladder_level_count` is the only control, and terraces are bands 1..n |
| `fill_floodplain` | 2 | `fill_valley_mask_holes` (default true) |
| `get_valley_centreline` | 3 | the tool itself |
| `centreline_loops` | 3 | `centreline_loops`, same meaning and same default of 5 |
| `trough_scaling_factor` | 3 | same name, same default, now auto-doubled on routing failure |
| `minimum_bank_elevation_window_radius` | 3 | `min_bank_elevation_window_radius_metres` |
| `get_valley_widths`, `width_node_spacing` | 3 | `write_valley_cross_sections`, `cross_section_tangent_window_nodes` |
| `valley_banks_search_radius` | 3 | `cross_section_max_half_length_metres` and the per-section not-found rule |
| `swath_width`, `swath_bin_spacing`, `print_swath_raster` | — | swath profiling is `lsdtt-swath-profile` |
| `print_channels_to_csv`, `print_junctions_to_csv` | 2 | `write_channel_network_lines`, `write_junction_network_points` |
| `convert_csv_to_geojson` | — | vector outputs are written directly; `write_vector_format` selects FlatGeobuf or GeoPackage |

Two structural differences worth stating plainly:

- **The legacy driver ran the whole pipeline in one process, and so had to recompute
  everything on every change.** The three-tool split exists so the expensive stages
  can be cached: a flat mask survives a change of level count, and a HAND raster
  survives a change of break algorithm.
- **The legacy tool decided how many levels there were and did not report whether it
  could.** lsdtt3 reports a per-window status and the resolvability floor it was
  tested against, and refuses to name a level count it cannot support.

---

## A worked sequence

```
# 1. which ground is flat
#    lsdtt-flat-surface-finder
admission_algorithm_choice: relief_or_gradient
max_admission_gradient: 0.18
relief_threshold_algorithm_choice: cdf_divergence
write_combined_flat_mask_raster: true

# 2. which level is it on
#    lsdtt-fluvial-surfaces
flat_mask_fname: mydem_combined_flat_mask.tif
source_points_fname: sources.csv
reference_field_algorithm_choice: hand
ladder_break_algorithm_choice: gap
ladder_level_count: 0
vertical_accuracy_metres: 2.0
write_hand_raster: true
write_valley_mask_raster: true

# 2b. iterate cheaply on the level count
hand_fname: mydem_hand.tif
ladder_level_count: 4

# 3. trace and measure the valley
#    lsdtt-valley-centreline
valley_mask_fname: mydem_valley_mask.tif
write_valley_centreline: true
write_valley_cross_sections: true
cross_section_spacing_metres: 200
```

Read `mydem_ladder_windows.csv` after step 2 and before describing anything as a
terrace level.

---

## References

**The methods this tooling descends from**

- **Clubb, F. J., Mudd, S. M., Milodowski, D. T., Valters, D. A., Slater, L. J.,
  Hurst, M. D., & Limaye, A. B. (2017).** Geomorphometric delineation of floodplains
  and terraces from objectively defined topographic thresholds. *Earth Surface
  Dynamics*, 5, 369–385. doi:10.5194/esurf-5-369-2017 — the two-threshold
  (gradient and relief-to-nearest-channel) classification, with both thresholds
  selected statistically from quantile–quantile plots.
- **Clubb, F. J., Weir, E. F., & Mudd, S. M. (2022).** Continuous measurements of
  valley floor width in mountainous landscapes. *Earth Surface Dynamics*, 10,
  437–456. doi:10.5194/esurf-10-437-2022 — the synthetic V-shaped trough, the
  routed valley centreline, and continuous orthogonal width measurement.

**Methods those two papers build on**

- **Nobre, A. D., Cuartas, L. A., Hodnett, M., Rennó, C. D., Rodrigues, G.,
  Silveira, A., Waterloo, M., & Saleska, S. (2011).** Height Above the Nearest
  Drainage — a hydrologically relevant new terrain model. *Journal of Hydrology*,
  404(1–2), 13–29. doi:10.1016/j.jhydrol.2011.03.051 — HAND, the default reference
  field.
- **Lashermes, B., Foufoula-Georgiou, E., & Dietrich, W. E. (2007).** Channel
  network extraction from high resolution topography using wavelets. *Geophysical
  Research Letters*, 34, L23S04. doi:10.1029/2007GL031140 — the quantile–quantile
  construction for locating a process-domain transition.
- **Passalacqua, P., Do Trung, T., Foufoula-Georgiou, E., Sapiro, G., & Dietrich,
  W. E. (2010).** A geometric framework for channel network extraction from lidar:
  nonlinear diffusion and geodesic paths. *Journal of Geophysical Research*, 115,
  F01002. doi:10.1029/2009JF001254 — the Perona–Malik preprocessing and QQ
  thresholding used in Clubb et al. (2017).
- **Perona, P., & Malik, J. (1990).** Scale-space and edge detection using
  anisotropic diffusion. *IEEE Transactions on Pattern Analysis and Machine
  Intelligence*, 12(7), 629–639. doi:10.1109/34.56205
- **Grieve, S. W. D., Mudd, S. M., Hurst, M. D., & Milodowski, D. T. (2016).** A
  nondimensional framework for exploring the relief structure of landscapes. *Earth
  Surface Dynamics*, 4, 309–325. doi:10.5194/esurf-4-309-2016 — the
  curvature-versus-window-size criterion for choosing a surface-fitting radius.
- **Lindsay, J. B. (2016).** Efficient hybrid breaching-filling sink removal methods
  for flow path enforcement in digital elevation models. *Hydrological Processes*,
  30(6), 846–857. doi:10.1002/hyp.10648 — the breaching used between trough carve
  iterations.
- **Braun, J., & Willett, S. D. (2013).** A very efficient O(n), implicit and
  parallel method to solve the stream power equation governing fluvial incision and
  landscape evolution. *Geomorphology*, 180–181, 170–179.
  doi:10.1016/j.geomorph.2012.10.008 — the O(n) topologically ordered flow
  traversal underlying the routing in all three tools.

**Classification and clustering**

- **Jenks, G. F. (1967).** The data model concept in statistical mapping.
  *International Yearbook of Cartography*, 7, 186–190 — the natural-breaks
  optimisation offered as `ladder_break_algorithm_choice: jenks`.
- **Fisher, W. D. (1958).** On grouping for maximum homogeneity. *Journal of the
  American Statistical Association*, 53(284), 789–798.
  doi:10.1080/01621459.1958.10501479 — the exact dynamic-programming solution that
  Jenks' method implements.

**Related lsdtt3 pages**

- [CLI reference: lsdtt-fluvial-surfaces](../cli_reference/cli_reference_lsdtt_fluvial_surfaces.md)
- [CLI reference: lsdtt-valley-centreline](../cli_reference/cli_reference_lsdtt_valley_centreline.md)
- [CLI reference: lsdtt-flat-surface-finder](../cli_reference/cli_reference_lsdtt_flat_surface_finder.md)
- [CLI reference: lsdtt-swath-profile](../cli_reference/cli_reference_lsdtt_swath_profile.md)
- [Flow routing in lsdtt3](flow_routing.md)
- [Surface metrics](surface_metrics.md)
