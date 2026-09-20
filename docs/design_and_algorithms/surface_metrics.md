# Surface Metrics in lsdtt3

This page documents how `lsdtt3` computes surface metrics (i.e. gradient, aspect, curvature, etc), what sign
conventions it uses, and what benchmark evidence currently exists for method
cost.

---

## Scope

`lsdtt3` surface metrics are built around local derivatives of elevation and
local texture metrics. Unlike some GIS tools, we do not use a fixed 4 point or 8 point stencil, 
instead we use the elevation values in a window and fit a differentiable surface using either 
a polynomial fit or wavelets. Hurst at al. 2012 and the Roering et al. papers referenced therein give
the rationale for this approach. Also see comments in this document below.

Surface metrics include:
- Local derivatives/curvatures from a quadratic surface fit.
- Hillshade and moving-window statistics.
- Roughness/texture metrics (MADk2, residual IQR, NVD).

Primary implementation entry points:

- `src/cpu/lsdtt_raster_surface_metrics_cpu.hpp`
- `src/cpu/lsdtt_raster_surface_metrics_cpu.cpp`
- `src/cpu/lsdtt_surface_texture_cpu.hpp`

---

## Core derivative framework (Florinsky-style notation)

`lsdtt3` uses the standard local quadratic form

`z(x,y) = a x^2 + b y^2 + c x y + d x + e y + f`

with derivative notation `p=z_x`, `q=z_y`, `r=z_xx`, `s=z_xy`, `t=z_yy`.

This notation follows Florinsky (2017), and all polyfit-derived metrics are
computed from these derivatives.

Reference:

- Florinsky, I.V. (2017). *An illustrated introduction to general
  geomorphometry*. Progress in Physical Geography 41(6), 723-752.
  doi:10.1177/0309133317733667.

---

## Curvature sign conventions

`lsdtt3` uses the strict mathematical sign convention for curvature fields:

- Concave-up terrain (channels/hollows) -> positive curvature.
- Convex-up terrain (ridges/crests) -> negative curvature.

This is the convention that follows directly from second derivatives of the
fitted surface.

`lsdtt3` docs explicitly note this differs from legacy GIS usage that follows
the Zevenbergen and Thorne (1987) sign choice for some reported curvature
products. The authors of lsdtt3 (well, Simon) have remained flummoxed why Zevenbergen and Thorne (1987)
flipped the sign of a standard mathematical quantity.

References:

- Florinsky (2017) for derivative/curvature conventions.
- Zevenbergen, L.W., and Thorne, C.R. (1987). *Quantitative analysis of land
  surface topography*. Earth Surface Processes and Landforms 12(1), 47–56.
  doi:10.1002/esp.3290120107

---

## Surface shape classification (SMAP)

`classify_surface_shape()` reproduces the CurveCalc shape classification (SMAP)
of the TopoCurve MATLAB code, labelling each pixel by the **signs** of two
curvature invariants formed from the principal curvatures `k1` (= k_max) and
`k2` (= k_min):

- Gaussian curvature `KG = k1 · k2`
- Mean curvature `KM = 0.5 · (k1 + k2)`

Because the classification is driven entirely by these signs, it inherits the
sign convention above — concave-up positive. A reader importing TopoCurve
results, or comparing against a GIS that follows Zevenbergen & Thorne, must
account for that before the class labels will agree.

| Class | Code | `KG` | `KM` | Form |
|---|---|---|---|---|
| `perfect_saddle` | −4 | < 0 | = 0 | Saddle |
| `basin` | −3 | > 0 | > 0 | Concave bowl |
| `synformal_saddle` | −2 | < 0 | > 0 | Saddle, concave-dominated |
| `synform` | −1 | = 0 | > 0 | Cylindrical valley |
| `plane` | 0 | = 0 | = 0 | Planar |
| `antiform` | 1 | = 0 | < 0 | Cylindrical ridge |
| `antiformal_saddle` | 2 | < 0 | < 0 | Saddle, convex-dominated |
| `peak` | 3 | > 0 | < 0 | Convex dome |

The integer codes are the TopoCurve codes and are written as float32 in the
`shape_class` raster, so they can be read back and compared directly against
TopoCurve output.

### The `kt` threshold

Exact zeros do not occur in real curvature fields, so a threshold decides what
counts as zero. `kt` (CLI: `shape_curvature_threshold`, default 0.0, units
1/m) zeroes any principal curvature with `|k| <= kt` **before** `KG` and `KM`
are formed.

This is the parameter that controls what the classification is about. With
`kt = 0`, no curvature is ever exactly zero, so every pixel lands in one of the
four non-planar classes and the planar and cylindrical classes are empty. Raising
`kt` collapses weak or long-wavelength structure — a broad background basin, say
— toward the planar and cylindrical classes, leaving only curvature that exceeds
the threshold to drive the label. Choosing `kt` is therefore choosing the scale
of shape you want to see, and it should be set relative to the curvature
magnitudes your DEM and fitting window actually produce.

The returned `gaussian_curvature` and `mean_curvature` rasters carry the
**post-threshold** values, so they remain consistent with the labels rather than
with the raw principal curvatures. A pixel is classified only where both input
curvature rasters are valid.

Enabled with `write_shape_classification_raster`, writing
`<write_fname>_<method>_shape_class`.

---

## Why polynomial fit, not a 4-point/3x3 stencil alone

`lsdtt3` supports both local-stencil and polynomial-fit contexts, but the main
surface-metric path is polynomial fit because it gives a coherent derivative set
at a user-controlled smoothing scale.

Why this is preferred over a small fixed stencil (e.g., ArcGIS-style 3x3
derivative products):

- One fit yields many metrics (`slope`, `aspect`, multiple curvatures,
  `rotor`, `shape_index`, etc.) consistently from the same local model.
- Scale is explicit through `polyfit_lambda_metres`, rather than fixed to one
  cell-neighborhood.
- The same fitted coefficients can be reused for multiple outputs, so
  multi-metric workloads are efficient.

Historical context in the codebase:

- The polyfit lineage is explicitly tied to Hurst et al. (2012) in the
  `lsdtt_raster_surface_metrics_cpu.hpp` comments.
- Roering papers are core references for hillslope-form/curvature analysis in
  `lsdtt3` developer documentation (especially 1999/2007/2008), which is why
  curvature quality and convention are treated as first-order concerns.

References present in-repo:

- Hurst, M.D., Mudd, S.M., Walcott, R., Attal, M., and Yoo, K. (2012).
  *Using hilltop curvature to derive the spatial distribution of erosion rates*.
  Journal of Geophysical Research 117, F02017.
  doi:10.1029/2011JF002057.
- Roering, J.J., Kirchner, J.W., and Dietrich, W.E. (1999). WRR 35(3):853-870.
- Roering, J.J., Perron, J.T., and Kirchner, J.W. (2007). EPSL 264(1-2):245-258.
- Roering, J.J. (2008). GSA Bulletin 120(9-10):1248-1262.

---

## Wavelet vs direct polynomial fit

`lsdtt3` exposes two interchangeable coefficient generators feeding the same
metric engine:

- `wavelet_quadratic_fit(...)` in `src/cpu/lsdtt_raster_spectral_cpu.*`
- `direct_polyfit(...)` in `src/cpu/lsdtt_raster_surface_metrics_cpu.cpp`

Both return `poly_fit_coeffs {a,b,c,d,e,f}` and therefore support identical
downstream metrics via `compute_polyfit_metrics(...)`.

The key distinction is where the derivatives are taken:

- **Direct** computes local polynomial moments in physical space with rolling
  sums, then solves for coefficients.
- **Wavelet** computes weighted moment fields by FFT-based convolution, then
  solves the same local 6x6 coefficient system per valid pixel.

So both produce a polynomial fit; wavelet is not a separate curvature formula,
it is a different way of assembling the polynomial-fit normal-equation terms.

### Wavelet approach

- FFTW-based Gaussian wavelet workflow.
- Intended as a spectral route to local quadratic coefficients.
- Requires FFTW at build/runtime.

Wavelet context for hilltop curvature:

- Struble and Roering (2021) show a continuous-wavelet route to hilltop
  curvature (Laplacian) and report very large runtime advantages over 2D
  polynomial fitting in their MATLAB experiments, while also discussing
  systematic underestimation behavior at moderate/high erosion rates.

### Direct approach

- Separable rolling-moment two-pass solver.
- Explicitly documented/tested as O(n) in raster cell count for fixed
  implementation constants.
- Default method in CLI/docs due to speed and fine-scale behavior.

---

## Roughness / texture metrics

`lsdtt3` includes multiple roughness families:

- **MADk2** (`calculate_madk2`) from second-order directional differences;
  outputs OR, DMC, RA.
- **Residual IQR** (`moving_window_iqr`) on raw or detrended residual surfaces.
- **NVD** (`compute_nvd`) via orientation-tensor eigenvalues (`s1,s2,s3`).

References encoded in headers:

- Trevisani, S., Teza, G., and Guth, P.L. (2023). *A simplified geostatistical
  approach for characterizing key aspects of short-range roughness*. CATENA
  223, 106927. doi:10.1016/j.catena.2023.106927 (MADk2)
- Trevisani, S., and Guth, P.L. (2025). *Surface Roughness in Geomorphometry:
  From Basic Metrics Toward a Coherent Framework*. Remote Sensing 17(23), 3864.
  doi:10.3390/rs17233864.
- Woodcock, N.H. (1977). GSA Bulletin 88(9), 1231-1236. (NVD eigenvalue shape)

---

## Compute cost: what is currently measured

From `docs/developer_notes/testing_and_benchmarking.md` on `alpine_medium` (720x720, 30 m):

### Direct polyfit scaling properties

| Check | Result |
|---|---|
| Window-size scaling (W=1,5,15,40) | 39.8 / 43.3 / 42.5 / 47.5 ms (~1.19x spread across 729x window-cell change) |
| Multi-metric marginal cost | 4 metrics: 47.9 ms; 27 metrics: 41.7 ms (fit dominates, extra metrics mostly output cost) |
| Peak memory behavior | ~58 B/pixel measured working-set model in direct polyfit analysis |

### Speed and memory vs TopoToolbox context

| Comparison | lsdtt3 (1t) | TopoToolbox | Notes |
|---|---:|---:|---|
| Slope | 35.4 ms | 40.4 ms | near parity |
| Profile curvature | 31.5 ms | 27.0 ms | near parity |
| Mean curvature | 27.4 ms | 28.6 ms | near parity |
| Plan curvature | 28.1 ms | 27.8 ms | near parity |
| All four curvatures | 30.2 ms | 159.8 ms | one-fit reuse benefit |
| Peak RSS (surface workflow) | 83-86 MB | 207-233 MB | documented benchmark fixture |

!!! warning "These are FIT-ONLY timings, on a different fixture"
    The table above comes from `run_surface_benchmark.py`, whose timer reads
    the CLI's `Running direct fit... done (X s)` line only. That is the
    **shared quadratic fit**, and it excludes the per-metric derivation that
    follows — which on a 1000x1000 tile costs a further 20-44 ms per metric.
    So the lsdtt3 column understates lsdtt3's end-to-end cost, and the "all
    four curvatures" row in particular is not four metrics' worth of work.
    The fixture is also `alpine_medium` (720x720), not the 1000x1000 used
    elsewhere.

    For end-to-end numbers that include both phases, and for the marginal cost
    per metric and the crossover against TopoToolbox, see the surface-metrics
    section of [Testing and benchmarking](../developer_notes/testing_and_benchmarking.md). Summary
    of the measured position, curvatures only on 1000x1000:

    | | fixed | per metric |
    |---|---|---|
    | lsdtt3 | 50.6 ms (the fit) | 24.3 ms |
    | TopoToolbox | 8.0 ms | 50.2 ms |

    lsdtt3 is cheaper from **two** curvature metrics onward.

---

## 3-way algorithmic comparison (operation ordering)

Below is the code-derived ordering comparison you requested.

| Method | Implementation basis | Dominant work per output pixel | Ordering with window size |
|---|---|---|---|
| LSDTopoTools2 legacy polyfit | `LSDRaster::calculate_polyfit_surface_metrics` in `../LSDTT_Development/src/LSDRaster.cpp` | For each center pixel: scan local window, accumulate 6 RHS moments over support points, solve 6x6 LU system, then compute metrics | approximately O(window_area) per pixel (plus constant-size solve) |
| lsdtt3 wavelet polyfit | `wavelet_quadratic_fit` in `src/cpu/lsdtt_raster_spectral_cpu.cpp` | FFT convolutions over padded grid for 6 basis fields, then per-pixel 6x6 solve | approximately O(P log P) global FFT term + O(1) per-pixel local solve (P = padded pixel count) |
| lsdtt3 direct polyfit | `direct_polyfit` in `src/cpu/lsdtt_raster_surface_metrics_cpu.cpp` | Two separable rolling-moment passes + constant-cost per-pixel coefficient recovery | approximately O(1) per pixel for fixed implementation constants; whole raster O(N) |

### What this means practically

- **LSDTopoTools2 brute-force local fit** spends cost proportional to support
  points in the kernel for every pixel. If window radius grows, per-pixel work
  grows roughly with area.
- **lsdtt3 direct** amortizes window dependence into rolling updates; per-pixel
  updates stay constant-order as window size changes.
- **lsdtt3 wavelet** shifts most work into FFTs over the full padded grid, then
  solves a small local system per pixel.

### Code-level operation notes

- In legacy LSDTopoTools2, the central loop explicitly iterates over kernel
  cells for each center pixel and accumulates the 6-element RHS vector before
  `LU.solve(bb)`.
- In lsdtt3 wavelet, one DEM forward FFT is done, then 6 basis convolutions
  are computed (kernel FFT + pointwise complex multiply + inverse FFT), then
  each valid pixel solves a 6x6 system.
- In lsdtt3 direct, pass 1 computes horizontal rolling sums (`H0,H1,H2,Hn`),
  pass 2 uses vertical prefix sums plus a compact local solve; no per-pixel
  loop over all window cells is needed.

---

## 3-way memory and speed table (evidence status)

| Method | Published measured speed in repo docs | Published measured memory in repo docs |
|---|---|---|
| LSDTopoTools2 legacy polyfit | Not currently benchmarked in this repo | Not currently benchmarked in this repo |
| lsdtt3 wavelet polyfit | Not currently a maintained timing table in docs | Not currently a maintained memory table in docs |
| lsdtt3 direct polyfit | Yes (`docs/developer_notes/testing_and_benchmarking.md`) | Yes (`docs/developer_notes/testing_and_benchmarking.md`) |

So: we now have the operation-order comparison for all three methods from code,
and measured speed/memory numbers for direct polyfit from the current benchmark
docs. Fully numeric 3-way speed/memory requires running new benchmark rows.

---

## References

- Florinsky, I.V. (2017). *An illustrated introduction to general
  geomorphometry*. Progress in Physical Geography 41(6), 723-752.
  doi:10.1177/0309133317733667.
- Zevenbergen, L.W., and Thorne, C.R. (1987). *Quantitative analysis of land
  surface topography*. Earth Surface Processes and Landforms 12(1), 47–56.
  doi:10.1002/esp.3290120107
- Hurst, M.D., Mudd, S.M., Walcott, R., Attal, M., and Yoo, K. (2012).
  *Using hilltop curvature to derive the spatial distribution of erosion rates*.
  Journal of Geophysical Research 117, F02017.
  doi:10.1029/2011JF002057.
- Roering, J.J., Kirchner, J.W., and Dietrich, W.E. (1999).
  *Evidence for nonlinear, diffusive sediment transport on hillslopes and
  implications for landscape morphology*. Water Resources Research 35(3),
  853-870.
- Roering, J.J., Perron, J.T., and Kirchner, J.W. (2007).
  *Functional relationships between denudation and hillslope form and relief*.
  Earth and Planetary Science Letters 264(1-2), 245-258.
- Roering, J.J. (2008). *How well can hillslope evolution models explain
  topography?* GSA Bulletin 120(9-10), 1248-1262.
- Struble, W.T., and Roering, J.J. (2021). *Hilltop curvature as a proxy for
  erosion rate: wavelets enable rapid computation and reveal systematic
  underestimation*. Earth Surface Dynamics 9, 1279-1300.
  doi:10.5194/esurf-9-1279-2021.
- Trevisani, S., Teza, G., and Guth, P.L. (2023). *A simplified geostatistical
  approach for characterizing key aspects of short-range roughness*. CATENA
  223, 106927. doi:10.1016/j.catena.2023.106927
- Trevisani, S., and Guth, P.L. (2025). *Surface Roughness in Geomorphometry:
  From Basic Metrics Toward a Coherent Framework*. Remote Sensing 17(23), 3864.
  doi:10.3390/rs17233864.
- Woodcock, N.H. (1977). *Specification of fabric shapes using an eigenvalue
  method*. Geological Society of America Bulletin 88(9), 1231-1236.
