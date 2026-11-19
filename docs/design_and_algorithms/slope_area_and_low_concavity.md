# Slope-area analysis and low-concavity reaches

A slope-area plot shows how gradient S changes with drainage area A. Along
fluvial channels it is often close to a power law (Hack 1957; Flint 1974):

```text
S = k_s A^(-theta)
```

Here `theta` is the concavity and `k_s` the steepness. On a log-log plot this is a
straight line with slope `-theta`. Most fluvial channels have theta of about
0.35-0.6 (e.g. Whipple 2004; Mudd et al. 2018). Other processes leave other
shapes:

- the hillslope part of the plot does not follow the power law at small areas;
- debris-flow valleys bend the plot over a range of small areas (see below).

`lsdtt-slope-area` measures this in two ways. It bins gradient by drainage
area over a whole DEM or basin. It can also mark the **channel reaches where
gradient changes little with area**. This page explains both. Every parameter
is listed in the
[CLI reference](../cli_reference/cli_reference_lsdtt_slope_area.md); the same
routines are in Python (see the [Python guide](../python/slope_area.md)).

---

## 1. The routing surface and the excluded pixels

Drainage area needs a surface that drains, so the DEM is conditioned before
routing:

1. **Carve** (`breach_before_fill`, default true). With
   `breach_algorithm_choice: least_cost` (the default), each pit is drained by
   lowering the cheapest path out of it. This is the least-cost breach of
   [lsdtt-raster-preprocessing](../cli_reference/cli_reference_lsdtt_raster_preprocessing.md),
   with the same parameter names. `lindsay2016` is the other choice.
2. **Fill** (`fill_before_routing`, default true). The fill raises what the
   carve did not solve.

Carving first avoids large filled areas. A filled depression is a flat, and its
gradient says nothing about the landscape.

Both steps change elevations, and **the gradient of a changed pixel is the
conditioning's, not the landscape's**. A carved reach has the breach's own
tiny gradient; a filled one is flat. So the CLI compares the routing surface
with the input DEM:

- a pixel above the input DEM is **filled** (code 1);
- a pixel below it is **carved** (code 2).

These are written to `_conditioned_pixels.tif`. With
`exclude_carved_and_filled_pixels: true` (the default) they are left out of
the bins, the D8 flow-path fits and the low-concavity fits. Drainage area
is still routed through them. The least-cost breach itself raises single-cell
pits to just below their lowest neighbour, so some pixels are filled even
though no separate fill step ran.

On a typical 30 m DEM the conditioning changes a small fraction of pixels. On
noisy or pitted data it can be much more. Check `_conditioned_pixels.tif`
before reading the plots.

## 2. Drainage area

`drainage_area_algorithm_choice` selects the area:

- **`dinf`** (default): D-infinity (Tarboton 1997). Flow is split between
  two downslope neighbours, so area spreads over hillslopes and grows
  smoothly down a channel.
- **`d8`**: all flow goes to the steepest neighbour. Area along a channel is
  constant between confluences and jumps at each one.

The choice is independent of the gradient method: the D8 flow-path gradient
always walks D8 receivers, and uses the chosen area only to pick its upstream
path (section 3).

## 3. Gradient

`slope_algorithm_choice` selects how S is measured.

### Polynomial fit (`polyfit`, default)

A quadratic surface is fitted to the elevations within `lambda` of each
pixel, and its gradient is the slope (the direct polyfit of
[surface metrics](surface_metrics.md)). One output set is written per value in
`polyfit_lambda_metres_list`.

- It measures the local surface gradient everywhere, hillslopes included. That
  is what the hillslope-to-channel transition on a slope-area plot needs.
- By default the fit runs on the **input DEM** (after interior-nodata
  handling, before any carve or fill), since the conditioning is for routing
  only. `polyfit_on_routing_surface: true` fits the conditioned surface instead.

### D8 flow-path gradient (`d8_flow_path`)

The gradient along the flow path, averaged over a window of total length
`W` (`d8_gradient_window_metres_list`, default 150 m) **centred** on the pixel:

- **downstream**, W/2 along the D8 receivers (or to an outlet);
- **upstream**, W/2 along the D8 donor of largest drainage area at each step
  (or to a divide). Upstream the path branches, and the largest-area branch
  is the main channel.

Elevation (on the routing surface) is fitted against flow distance by least
squares, and S is the negative slope of that line. This is the same
along-channel gradient
[`lsdtt-channel-extraction`](../cli_reference/cli_reference_lsdtt_channel_extraction.md)
uses. Carved and filled cells are walked through but not fitted, and their own
gradient is nodata.

A fit is valid if it has at least 3 cells and they span at least W/2 less one
diagonal step. So a pixel at a divide or an outlet, with a window on one side
only, still gets a gradient if that side is long enough. W must be at least
four cell sizes. The default of 150 m is five cells on 30 m data.

The D8 gradient is valid on fewer pixels than polyfit: a pixel needs a long
enough flow path around it, clear of conditioned cells. This is by design. The
D8 gradient is meant for channels, and polyfit for the
hillslope-to-channel transition.

## 4. Binning

The pixels used are those with a valid gradient and area, area > 0, inside any
basin clip, and not excluded. They are split into `n_log_bins` bins of equal
width in ln A. Each bin reports:

- the geometric centre of the bin;
- the 5th, 25th, 50th, 75th and 95th percentiles of S;
- the number of pixels.

These are written as `_lambda_XXm.csv` (polyfit) or `_d8window_XXm.csv` (D8).
The summary CSV reports, per lambda or window, the bin centre with the largest
median gradient: the area where the plot turns from rising (hillslopes) to
falling (channels).

## 5. Low-concavity reaches

With `run_low_concavity_detection: true` the CLI looks for **channel reaches
where gradient changes little with drainage area**: where theta is close to 0
rather than 0.35-0.6.

### The fit

At each channel pixel (D8 contributing pixels >=
`threshold_contributing_pixels`):

1. Take the pixel and the next pixels **downstream** along the D8 receivers,
   `low_concavity_window_pixels` in all (default 10), stopping at an outlet.
   The window goes downstream because the path upstream branches. It runs
   through confluences, and `n_junctions_in_window` records how many it
   crosses.
2. Drop pixels without a gradient and carved or filled pixels. If the pixel
   itself is carved or filled it gets no fit.
3. Fit `log10 S = a + b log10 A` by least squares over the remaining pairs
   (at least `low_concavity_min_pixels`, default 5, and not all with the same
   area). Then theta = -b. The CLI also writes:
   - the standard error of theta,
     `theta_se = sqrt(SS_res / (n - 2) / sum (x - mean x)^2)` with x = log10 A;
   - `r2`;
   - the spread of log10 A in the window, `log10_area_range`.
4. Flag the pixel:
   - **1** if |theta| < `low_concavity_max_abs_concavity` (default 0.15);
   - **0** if the fit is valid but |theta| is larger;
   - **-1** if there is no valid fit.

   The threshold is on the absolute value, so a weakly convex reach (small
   negative theta) counts too. Theta keeps its sign in the outputs.

S is always the D8 flow-path gradient of section 3, with its own window
`low_concavity_gradient_window_metres` (default 150 m), whatever
`slope_algorithm_choice` is. A is the area of section 2.

### Segments

The per-pixel flags are noisy, so they are also written as a vector layer:

- The channel network is cut into links: each link runs from a source, or the
  pixel below a confluence, down to the next confluence or outlet.
- Each link is split where the flag changes.
- Each segment carries its flag, length, median theta, median theta_se, median
  slope, area range and Strahler order.
- `low_concavity_min_run_pixels` relabels flagged runs shorter than that as 0.
  The short run is merged with its neighbours, not dropped, so the segments
  still cover the whole network.

### Reading the flags

**A flag says that gradient did not vary with area in the window. It does
not say that the window follows a clean power law with a small exponent.**
Where S scatters about a constant, the fitted theta is near 0 with a low
`r2`. Use `theta_se` and `r2` to tell the two apart.

**Area method.**

- With **D8 area**, area along a channel is constant between confluences and
  jumps at each one, while the averaged gradient does not jump. A window
  across a confluence therefore pulls theta toward 0, and flags tend to sit
  at confluences. Check `n_junctions_in_window`.
- **D-infinity area** grows along the channel by lateral inflow and does not
  have this effect.

**Where low concavity is expected.**

- **Debris-flow channels.** In steep unglaciated valleys cut by debris flows,
  the slope-area plot is curved: at small areas gradient increases only
  slowly as area decreases, and becomes a fluvial power law further down
  (Stock and Dietrich 2003). Montgomery and Foufoula-Georgiou (1993) read the
  break in the slope-area plot of mountain basins as the change from steep
  debris-flow channels to lower-gradient alluvial channels. The scaling break
  is commonly at about 0.1-1 km² (Stock and Dietrich 2003). Channel heads
  must reach below that (lower `threshold_contributing_pixels`) for the
  debris-flow reaches to be on the network.
- **Threshold hillslopes.** Where hillslopes are close to a limiting gradient,
  gradient no longer increases with relief or area (e.g. Burbank et al. 1996;
  Montgomery 2001). The finder runs only on channel pixels, so threshold
  hillslopes are flagged only when the channel threshold reaches onto them.
- **Alluvial fans.** Fans have a nearly uniform gradient down-fan, but their
  flow is distributary, and the routing here is convergent. Both D8 and
  D-infinity route a fan's flow along threads that a real fan does not have
  for long, so the area on a fan is a routing product. Treat flags on fans as
  a sign of uniform gradient rather than as a measured slope-area relation.
  This interpretation has not been checked against mapped fans.

## 6. Reproducing earlier runs

Before 2026-11-18 lsdtt-slope-area:

- filled the DEM without carving;
- kept the filled pixels;
- fitted polyfit to the filled DEM;
- named the area method `flow_algorithm`.

The [CLI reference](../cli_reference/cli_reference_lsdtt_slope_area.md#reproducing-runs-made-before-2026-11-18)
gives the parameter edits that reproduce those runs byte for byte. A
regression test pins the outputs.

## References

- Burbank, D. W., Leland, J., Fielding, E., Anderson, R. S., Brozovic, N.,
  Reid, M. R. and Duncan, C. (1996). Bedrock incision, rock uplift and
  threshold hillslopes in the northwestern Himalayas. *Nature*, 379, 505-510.
- Flint, J. J. (1974). Stream gradient as a function of order, magnitude, and
  discharge. *Water Resources Research*, 10(5), 969-973.
- Hack, J. T. (1957). Studies of longitudinal stream profiles in Virginia and
  Maryland. *USGS Professional Paper* 294-B.
- Montgomery, D. R. (2001). Slope distributions, threshold hillslopes, and
  steady-state topography. *American Journal of Science*, 301, 432-454.
- Montgomery, D. R. and Foufoula-Georgiou, E. (1993). Channel network source
  representation using digital elevation models. *Water Resources Research*,
  29(12), 3925-3934.
- Mudd, S. M., Clubb, F. J., Gailleton, B. and Hurst, M. D. (2018). How
  concave are river channels? *Earth Surface Dynamics*, 6, 505-523.
- Stock, J. and Dietrich, W. E. (2003). Valley incision by debris flows:
  evidence of a topographic signature. *Water Resources Research*, 39(4),
  1089. doi:10.1029/2001WR001057.
- Tarboton, D. G. (1997). A new method for the determination of flow
  directions and upslope areas in grid digital elevation models. *Water
  Resources Research*, 33(2), 309-319.
- Whipple, K. X. (2004). Bedrock rivers and the geomorphology of active
  orogens. *Annual Review of Earth and Planetary Sciences*, 32, 151-185.
