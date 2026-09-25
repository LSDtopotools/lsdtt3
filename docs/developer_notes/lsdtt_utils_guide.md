# lsdtt3 utils contents

This document explains what belongs in `src/utils/`, what does **not** belong in
`src/utils/`, and how the utility-layer statistics and comparison functions are
intended to be used.

---

## Purpose of `src/utils/`

The `utils` layer is the home for:

- generic mathematical routines
- generic numerical / statistical routines
- lightweight containers or helper types
- metadata-only geometry helpers
- functions that operate on already-valid flat data
- code that does **not** depend on GDAL / OGR
- code that is not tied to a specific backend implementation

In short, `utils` should contain things that are reusable across many parts of
`lsdtt3`, and which could also make sense in another numerical C++ project.

Good examples are:

- descriptive statistics
- histograms, PDFs, CDFs
- binning utilities
- distribution-comparison tests
- circular-data comparison tools
- simple array/container helpers
- pure raster-grid geometry helpers

---

## What should **not** go into `utils`

`utils` is **not** the place for:

- GDAL / OGR I/O
- CRS parsing or coordinate transformation through GDAL
- raster file loading / writing
- backend factories
- nodata masking rules tied to raster objects
- CPU-specific raster traversal
- GPU-specific kernels
- flow-routing construction
- moving-window extraction over raster objects
- basin extraction from flow-routing objects

Those belong elsewhere:

- `src/core/` for shared infrastructure, GDAL/OGR-aware helpers, I/O, and
  backend-neutral interfaces
- `src/cpu/` for concrete CPU implementations
- `src/gpu/` for concrete GPU implementations

---

## The valid-data contract

A core design rule in `lsdtt3` is:

> `utils` routines should assume the input data are already valid.

This means the utility layer should **not** do:

- nodata filtering
- mask checking
- raster flattening
- moving-window extraction
- basin extraction
- NaN screening tied to raster semantics

Instead, those tasks should happen **before** data reach `utils`.

### Why this rule exists

This keeps `utils`:

- simple
- reusable
- backend-independent
- easy to test
- free from raster-specific policy decisions

It also avoids mixing two distinct concerns:

1. **What are the valid values?**
2. **What statistics or tests should be computed on those values?**

---

## Where nodata handling should happen

Nodata checks and validity filtering belong in the code that owns the data
structure or the backend-specific extraction logic.

Typical examples:

- **whole-raster statistics**  
  The raster object or a CPU/GPU raster helper should flatten valid pixels and
  pass those values to `utils`.

- **moving-window statistics**  
  The moving-window extraction code should collect valid values inside the
  window and pass those values to `utils`.

- **basin statistics**  
  The basin extraction code should use flow-routing relationships (for example,
  via `s_vector`) to collect valid basin values and pass those values to
  `utils`.

### In practice

That means functions like these conceptually belong outside `utils`:

- `flatten_to_valid(...)`
- `window_to_valid(...)`
- `extract_valid_basin_values(...)`

These are naturally CPU/GPU-side or object-side tasks.

---

## Current contents of `src/utils/`

The utils layer contains the following headers. The sections below describe the
long-standing core modules in detail; the rest are summarised in the table.

| Header | Purpose |
|---|---|
| `lsdtt_array2d` | Lightweight 2D array over a contiguous `std::vector` |
| `lsdtt_binary_array2d_tracing` | Boundary tracing on binary images (Gonzalez, Woods & Eddins pp. 557–560); returns row/col in the input array's coordinate system |
| `lsdtt_distributions` | Parametric distribution fitting for valid flat numeric samples |
| `lsdtt_histogram` | Histogram, PDF/CDF and binning utilities |
| `lsdtt_image_filters` | Generic 2-D image-filter primitives on flat row-major arrays |
| `lsdtt_linear_algebra` | Lightweight linear algebra utilities |
| `lsdtt_natural_breaks` | Jenks natural-breaks classification (`lsdtt_jenks_options`) |
| `lsdtt_numeric_utils` | Numeric helpers, including inclusive `start, start + step, …` ranges |
| `lsdtt_polygon_geometry` | Header-only planar polygon geometry |
| `lsdtt_random` | Fully specified random numbers (splitmix64) |
| `lsdtt_raster_geometry` | Header-only raster geometry from metadata alone |
| `lsdtt_regression` | Linear, orthogonal and polynomial regression for valid paired samples |
| `lsdtt_signal1d` | 1-D profile signal processing: moving average, LOWESS, derivative |
| `lsdtt_splines` | Parametric spline fitting for valid ordered planar data |
| `lsdtt_statistical_tests` | Nonparametric two-sample tests and distribution distances |
| `lsdtt_statistics` | Basic 1-D descriptive statistics |
| `lsdtt_string_utils` | Header-only string helpers shared by CLIs and library |
| `lsdtt_thread_utils` | Thread-count resolution |

`lsdtt_image_filters` is worth singling out: it holds the native Otsu,
Canny, Hough and adaptive-threshold primitives that replaced OpenCV. Where a
primitive stands in for an OpenCV routine it reproduces that routine's rounding
rules, and names the routine it mirrors in its comment, so raster-level
operations built on it match the previous output exactly. Its inputs are dense
images — every element is data, with no nodata or mask handling; the raster
layer in `src/cpu/` decides what nodata becomes before calling in.

### `lsdtt_array2d.hpp`
A lightweight 2D array wrapper over a contiguous `std::vector`.

This is a good `utils` component because it is:

- generic
- reusable
- independent of GDAL
- not backend-specific

### `lsdtt_raster_geometry.hpp`
Header-only raster geometry helpers that work only with raster metadata and grid
geometry.

These functions are appropriate for `utils` because they are:

- metadata-only
- pure geometry helpers
- free of file I/O
- free of reprojection logic

### `lsdtt_statistics.hpp/.cpp`
Descriptive statistics for valid flat numeric samples.

Typical outputs include:

- mean
- median
- percentiles
- variance
- standard deviation
- mean absolute deviation
- median absolute deviation
- interquartile range

These functions assume the caller has already supplied valid values.

### `lsdtt_histogram.hpp/.cpp`
Histogram, PDF, CDF, and binning utilities for valid flat numeric samples.

This module includes:

- equal-width histograms
- log-binned histograms
- natural-breaks binning
- PDF and CDF generation
- binning raw values
- descriptive statistics within bins

### `lsdtt_statistical_tests.hpp/.cpp`
Nonparametric distribution-comparison tools for valid flat numeric samples.

This module is intended for comparing two datasets without assuming normality.

---

## What kind of functions belong in `utils`

As a rule of thumb, a function belongs in `utils` if all of the following are
true:

- it works on plain numeric values, vectors, iterator ranges, or simple helper
  containers
- it does not need GDAL / OGR
- it does not need raster file handles
- it does not need backend selection
- it does not need to know about nodata policy
- it does not care whether the caller is CPU or GPU code

Examples:

- `compute_descriptive_stats(...)`
- `compute_histogram(...)`
- `compute_pdf(...)`
- `compute_cdf(...)`
- `compute_kolmogorov_smirnov_two_sample_test(...)`
- `compute_wasserstein_distance(...)`

---

## What kind of functions should stay outside `utils`

A function should probably **not** go into `utils` if it needs to know about:

- raster masks
- nodata sentinels
- CRS
- GDAL vector layers
- moving-window boundary handling on a raster object
- flow-routing object internals
- device memory or CUDA kernels

Examples that should stay outside `utils`:

- “flatten this raster to valid values”
- “extract valid values inside this moving window”
- “extract basin values from this outlet pixel”
- “transform longitude/latitude to raster CRS”
- “load this GeoTIFF and write output COG”

---

## Statistical philosophy in `utils`

Most `lsdtt3` terrain datasets are not normally distributed, so the statistics
and comparison tools in `utils` should be designed with that in mind.

The emphasis is therefore on:

- robust descriptive statistics
- empirical distribution tools
- nonparametric tests
- transport/distance measures
- circular-data comparisons

---

# Statistical tests and distances in `utils`

This section explains the main comparison tools and when to use them.

---

## Kolmogorov–Smirnov two-sample test

### What it measures
The Kolmogorov–Smirnov (KS) test compares two empirical cumulative
distribution functions and finds the **maximum absolute difference** between
them.

### Use it when
- you want a general nonparametric comparison of two 1D distributions
- you want a simple, widely known test
- you care about the largest mismatch anywhere in the CDF

### Strengths
- simple interpretation
- nonparametric
- no normality assumption
- works directly on continuous 1D data

### Limitations
- less sensitive to tail differences than Anderson–Darling
- based on the largest single discrepancy, not overall mismatch
- not the right tool for circular data

### Good use cases
- comparing elevation distributions between two basins
- comparing slope distributions between two terrain classes
- checking whether two extracted raster subsets differ overall

---

## Kuiper test

### What it measures
The Kuiper test is a circular-data analogue of KS. It uses the sum of the
largest positive and negative differences between two circular empirical CDFs.

### Use it when
- your data are angular or circular
- values wrap around a period
- you want a circular test that is not biased by where the origin is chosen

### Strengths
- suitable for circular variables
- sensitive across the whole circle
- invariant to circular origin choice

### Limitations
- test statistic is unitless
- less directly interpretable than a transport distance
- for axial/orientation data you must use the axial version, not the ordinary
  directional version

### Good use cases
- comparing aspect distributions
- comparing flow-direction or bearing distributions
- comparing orientation distributions after converting them appropriately

---

## Directional versus axial circular data

This distinction is extremely important.

### Directional / bearing data
Use directional tests when **absolute direction matters**.

Examples:
- aspect as a bearing
- downslope azimuth
- wind direction
- channel azimuth

In directional data:

- `10°` and `190°` are different
- `0°` and `360°` are equivalent because of periodic wrap

### Axial / orientation data
Use axial tests when **facing does not matter**, only orientation.

Examples:
- line orientation
- ridge orientation
- fracture orientation
- any undirected axis

In axial data:

- `10°` and `190°` are the same orientation
- the effective period is `180°` rather than `360°`

### Rule
- for bearings: use **directional Kuiper** and **directional circular
  Wasserstein**
- for orientations: use **axial Kuiper** and **axial circular Wasserstein**

---

## Wasserstein distance

### What it measures
The Wasserstein distance measures how much “transport” is needed to move one
distribution into the other.

It is often called an earth-mover distance.

### Use it when
- you want to know **how different** two distributions are
- you want a result in the same units as the data
- you care about effect size, not only a hypothesis test

### Strengths
- directly interpretable
- reported in data units
- very useful as a “distance” rather than just a test statistic

### Limitations
- not itself a p-value-based test
- standard 1D Wasserstein does not work correctly for circular wrap-around

### Good use cases
- comparing two basin elevation distributions in metres
- comparing channel steepness distributions
- asking “how different are these datasets, in the units I care about?”

---

## Circular Wasserstein distance

### What it measures
Circular Wasserstein is the circle-aware counterpart of ordinary Wasserstein.
It measures transport on a circle, respecting the fact that angular data wrap.

### Use it when
- your data are periodic or circular
- you want an angular **distance**
- you want to quantify difference in angular units such as degrees or radians

### Strengths
- interpretable in angular units
- much better than linear Wasserstein for circular data
- complements Kuiper nicely

### Good use cases
- comparing aspect distributions
- comparing bearing distributions
- comparing orientation distributions (with the axial version)

---

## Rotation-invariant circular Wasserstein

### What it measures
This is a circular Wasserstein distance that allows one global rotation before
comparing the two distributions.

### Use it when
- you care about distribution **shape** on the circle
- absolute phase / absolute orientation is not important
- one distribution may simply be rotated relative to the other

### Strengths
- separates shape difference from simple angular offset
- useful when circular origin or absolute rotation is arbitrary

### Limitations
- answers a different question from ordinary circular Wasserstein
- should not be used if absolute bearing itself is scientifically important

### Good use cases
- comparing directional shapes up to rotation
- comparing orientation distributions where absolute reference direction is not
  the main interest

---

## Cramér–von Mises two-sample test

### What it measures
Cramér–von Mises compares two empirical distributions using the **overall
integrated difference** between their empirical CDFs.

### Use it when
- you care about overall distribution mismatch, not just the single largest gap
- you want a smoother whole-distribution comparison than KS

### Strengths
- sensitive to broad distributional differences
- often complements KS well

### Limitations
- less familiar to many users than KS
- interpretation of the statistic is less intuitive than Wasserstein

### Good use cases
- comparing overall basin elevation distributions
- comparing overall slope or curvature distributions

---

## Anderson–Darling two-sample test

### What it measures
Anderson–Darling is similar in spirit to Cramér–von Mises, but gives more
weight to the **tails** of the distribution.

### Use it when
- tail behaviour matters
- extremes are scientifically important
- you care about low-end or high-end differences more strongly

### Strengths
- more sensitive than KS to differences in extremes
- particularly useful when tails matter

### Limitations
- not as directly interpretable as a distance
- like all such tests, p-values need caution for spatially dependent data

### Good use cases
- comparing relief extremes
- comparing high-slope tails
- comparing extreme curvature or roughness behaviour

---

## Energy distance

### What it measures
Energy distance is a nonparametric distance between two distributions.

### Use it when
- you want a distance-style summary of whole-distribution difference
- you want something complementary to Wasserstein and KS

### Strengths
- nonparametric
- zero only when distributions match
- captures overall distributional difference

### Limitations
- less directly interpretable in domain units than Wasserstein
- best thought of as a distance measure rather than a simple hypothesis test

### Good use cases
- comparing two terrain-property distributions when you want a general
  distribution-distance measure
- complementing KS or Wasserstein in summary tables

---

# Which test or distance should I use?

## If I want a simple general comparison
Use:

- **Kolmogorov–Smirnov**

## If I want to compare overall distribution shape
Use:

- **Cramér–von Mises**
- **energy distance**

## If I care about tails or extremes
Use:

- **Anderson–Darling**

## If I want an interpretable distance in the same units as the data
Use:

- **Wasserstein distance**

## If my data are circular and absolute bearing matters
Use:

- **Kuiper**
- **circular Wasserstein**

## If my data are circular and only orientation matters
Use:

- **axial Kuiper**
- **axial circular Wasserstein**

## If I care about circular distribution shape but not absolute phase
Use:

- **rotation-invariant circular Wasserstein**
- or its axial counterpart for orientation data

---

# Important caution for raster and basin data

Many standard statistical tests assume observations are approximately
independent.

Raster-derived values, windowed values, and basin pixel values are often:

- spatially autocorrelated
- topographically correlated
- not independent in the classical statistical sense

Because of that:

- test statistics are often still useful
- distances are often very useful
- p-values should be interpreted cautiously unless dependence has been handled

In many terrain-analysis workflows, the most useful outputs are often:

- the test statistic
- the effect size or distance
- clear figures of the compared distributions

rather than a p-value alone.

---

# Recommended development rule for `utils`

When adding a new utility function, ask:

> Could this function operate on a vector of valid values from any source,
> without needing to know about raster masks, nodata values, CRS, or a backend?

If the answer is **yes**, it probably belongs in `utils`.

If the answer is **no**, it probably belongs in:

- `core`
- `cpu`
- or `gpu`

depending on what it depends on.

---

# Summary

`src/utils/` should contain:

- generic numerical functions
- generic statistics
- generic histogram/binning utilities
- nonparametric distribution comparison tools
- circular-data comparison tools
- lightweight reusable helper types

`src/utils/` should **not** contain:

- nodata filtering
- raster flattening
- CRS / GDAL logic
- backend-specific traversal
- CPU/GPU-specific extraction logic

The main rule is:

> code in `utils` assumes the data are already valid.

That rule keeps `utils` clean, reusable, and consistent with the overall
architecture of `lsdtt3`.
