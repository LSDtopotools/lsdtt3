# Spectral analysis of gridded surfaces

This page describes the spectral method behind `lsdtt-spectral-analysis` and the
library code it calls (`src/cpu/lsdtt_raster_spectral_cpu.hpp` and `.cpp`). It
covers what the method measures, each processing step, how significance is
tested, how the orientation of a pattern is reported, how the Monte Carlo
background is made reproducible, how the test suite exercises the code, and what
every output file contains.

Parameter names, defaults and the generated parameter table are in the
[spectral analysis CLI reference](../cli_reference/cli_reference_lsdtt_spectral_analysis.md).
This page does not repeat that table.

---

## 1. Purpose

The method detects characteristic spatial scales in a gridded surface. A pattern
that repeats at a fixed spacing puts a concentration of power at one wavelength
in the two-dimensional spectrum. The method asks whether that concentration is
larger than would be expected from a power-law (fractal) background with the
same overall variance.

The classic application is ridge and valley spacing, where the spacing is a
hillslope-length or valley-spacing scale. Perron et al. (2008) describe the
spectral approach for landscapes, and Perron et al. (2009) describe the formation
of evenly spaced ridges and valleys that motivates looking for such a scale.
Full references are in section 9.

The same procedure applies to any quasi-periodic pattern on a raster:

- dunes and other bedforms;
- patterned vegetation such as banded or "tiger bush" vegetation;
- vegetation density or cover rasters derived from imagery;
- any other field with a repeating spacing and orientation.

### What the input can be

The input does not need to be elevation. The code accepts any scalar raster
with a projected coordinate system in metres. The metre requirement is enforced,
because wavelengths are reported in metres.

Units follow from the input. The periodogram carries the square of the input
units. An elevation raster gives m² per spectral cell. A dimensionless
vegetation-density raster gives dimensionless power. The normalised PSD
(section 2.11) is dimensionless in every case, so significance is always a
ratio to the local background, and the test does not depend on the input units.

For interpretation: the wavelength is the spacing of the pattern in the raster
values. It is not automatically a hillslope or valley length. Whether a detected
spacing has a landscape meaning depends on what the raster represents.

---

## 2. Method

The processing chain is:

1. Load the raster and check that it is in projected metres.
2. Remove a planar trend from the valid pixels.
3. Apply a two-dimensional window.
4. Zero-pad each dimension to the next power of two.
5. Take the two-dimensional FFT and form the periodogram.
6. Collapse the periodogram onto radial frequency and log-bin it.
7. Estimate the rollover frequency and the background exponent from the binned spectrum.
8. Build the Monte Carlo background spectrum with the same padding and window.
9. Normalise the observed spectrum by the background and apply the significance tests.

### 2.1 Detrending

A least-squares plane is fitted to all valid pixels:

```
z(r, c) ≈ a·c + b·r + d
```

where `r` is the row index and `c` is the column index, both counted from zero,
and `z` is the input value. The plane is subtracted from every valid pixel to give
the detrended field `z'`. Only a planar trend is removed. Larger-scale curvature
of the terrain is not removed and leaks into the lowest frequencies.

The detrended variance `σ²` is the population variance of `z'` over the valid
pixels. It is computed before windowing and is used later to scale the background.

Nodata pixels are excluded from the plane fit. In the FFT input they are set to
zero, which is the detrended mean. They still carry window weight, so a tile with
holes has less effective signal than its raster size suggests.

### 2.2 Window

The default window is a two-dimensional Hann window whose support is the ellipse
inscribed in the raster. For a raster with `n_rows` rows and `n_cols` columns:

```
a  = (n_cols - 1) / 2          (column semi-axis, pixels)
b  = (n_rows - 1) / 2          (row semi-axis, pixels)
ρ  = distance of the pixel from the centre
ρ' = a·b / sqrt(b²·cos²θ + a²·sin²θ)    (distance from the centre to the ellipse
                                          along the same direction θ)

w  = 0.5 · (1 + cos(π·ρ / ρ'))   for ρ < ρ'
w  = 0                            otherwise
```

For a square raster the window is radially symmetric. For a rectangular raster
it is elliptical. Pixels outside the inscribed ellipse, including the corners,
get zero weight.

The library also provides a unity window (`none`) and a Hamming window.
The command-line tool does not expose either option and always uses Hann.

The window sum of squares is

```
WSS = Σ w²
```

and it is used to normalise the periodogram.

### 2.3 Zero-padding

Each dimension is padded with zeros to the next power of two:

```
P_r = smallest power of two ≥ n_rows
P_c = smallest power of two ≥ n_cols
```

A 200-row raster is therefore transformed on a 256-row grid, with 56 rows of
zeros. The padding is independent in the two directions.

### 2.4 Periodogram

The windowed, padded field is transformed with a complex forward FFT (FFTW).
The periodogram is

```
X(k, l) = FFT₂[ w · z' ](k, l)

P(k, l) = |X(k, l)|² / (P_r · P_c · WSS)
```

Each value `P(k, l)` is the variance carried by one cell of the discrete
spectrum, in m² for an elevation input. This follows from Parseval's theorem:
the sum over all cells of the padded spectrum equals the window-weighted variance
of the detrended field,

```
Σ_{k,l} P(k, l) = Σ_{r,c} (w·z')² / WSS
```

The division by the frequency spacing that would turn this into a density per
(cycles/m)² is not applied, so the PSD is a variance per DFT cell (m² for an
elevation input), not a density per unit frequency.

### 2.5 Frequency axes and the fftshift layout

The frequency spacings and the Nyquist frequency are

```
dfx = 1 / (cell_size · P_c)      cycles/m, east
dfy = 1 / (cell_size · P_r)      cycles/m, north
f_N = 1 / (2 · cell_size)        cycles/m, Nyquist (both axes)
```

The 2D arrays are stored in fftshift order, so zero frequency (DC) sits at row
`P_r / 2` and column `P_c / 2` (integer division). For a cell at row `r` and
column `c` the frequency components are

```
f_east  = (c - P_c/2) · dfx
f_north = (P_r/2 - r) · dfy
```

Row 0 is the northern edge of the spectrum, because the raster is north-up and
northing increases as the row index decreases. The 2D rasters are written north-up
for the same reason (section 6).

### 2.6 Radial collapse

The periodogram is reduced to one radial frequency per cell:

```
f = sqrt(f_east² + f_north²)
```

The collapse uses the left half of the shifted spectrum (columns `0` to `P_c/2`
inclusive, all rows). It keeps only cells with `f ≤ f_N`. Cells in the corners of
the square spectrum with `f > f_N` are dropped here, although they are kept in the
2D output (section 3).

Each kept cell contributes a radial value equal to `2 · P(k, l)`. The factor of 2
is the legacy convention for a one-sided spectrum of a real field. It is applied to
every kept cell, including the DC column and the Nyquist column (column 0), which
have no separate mirror in the half spectrum. It changes the absolute radial
values, but the observed and background spectra get the same factor, so the
normalised values are unaffected.

### 2.7 Log binning

The radial values are binned in `log10(f)`. With bin width `w` (CLI parameter
`spectral_log_bin_width`), the bin of a point is `floor(log10(f) / w)`. Bin edges
therefore sit at fixed powers of `10^w`. Each occupied bin is reported at the mean
of the frequencies and the mean of the powers of its members. The centre of the
bin in log space is not used.

A bin of width `w` spans a factor of `10^w` in frequency. The CLI default of 0.05
gives bins that span a factor of about 1.12. The library default of 0.1 gives bins
that span a factor of about 1.26.

Only points with positive, finite frequency and power are binned. The DC point
is therefore excluded at this stage.

The number of occupied bins is roughly

```
log10( max(P_r, P_c) / 2 ) / w
```

because the lowest non-zero radial frequency is `min(dfx, dfy)` and the highest
is the Nyquist frequency `f_N`.

### 2.8 Rollover and the fractal fit

The rollover is the frequency at which the spectrum stops following a single
power law. The code finds it by repeated regression on suffixes of the binned
spectrum.

Let there be `n` binned points with frequency `f_j` and power `P_j`. The
regression variables are `x = log10 f` and `y = log10 P`.

1. For each start index `s` from `0` up to `n - floor(n/5) - 1`, fit an ordinary
   least-squares line to all bins from `s` to `n - 1`. The top fifth of the bins
   is never a start point, but it is always part of the suffix.
2. Define the exponent and its goodness of fit for each start:

   ```
   β_s = -slope_s        (dimensionless spectral exponent)
   R²_s                  (coefficient of determination of that fit)
   ```

3. The rollover index is the start with the largest `β_s`. If that index is the
   first or last bin, the start with the largest `R²_s` is used instead. If that
   is also the first bin, the last bin is used. The first bin is therefore never
   reported as the rollover when there is more than one bin.
4. The rollover frequency is `f` at the rollover index, and the rollover exponent
   is `β` at that index.

A sub-rollover exponent is then fitted to the bins below the rollover, from bin 2
up to and including the rollover index. Bins 0 and 1 are skipped when the rollover
index is above 2. This exponent is the one used for the background.

The code describes this search as legacy-inspired. The steps above are a
description of the code. They are not a transcription of the 2008 paper.

The sub-rollover fit is guarded to need at least two points. With at least three
bins and a rollover index of at least one, the fit always has at least two points,
so the guard does not trigger in normal use.

#### Fallback for the background exponent

The background exponent is chosen in this order:

1. the sub-rollover exponent, if its fit gives a finite positive value;
2. the rollover exponent, if it is finite and positive;
3. the constant 2.0.

The value written to the outputs as `sub_rollover_beta` is the value actually used,
so it can be a fallback. Section 7 covers small tiles, where the fallbacks apply.

### 2.9 Monte Carlo fractal background

The background is an average of synthetic fractal surfaces that share the observed
tile's dimensions, padding, window and variance. For iteration `i = 0 … N - 1`
(`N` is `n_background_iterations`), the seed is

```
s_i = base + i
base = background_seed        if background_seed ≥ 0
base = 1000                   if background_seed = -1 (default)
```

Each synthetic surface is generated as follows:

1. Fill the padded `P_r × P_c` grid with white noise `u(r, c)` uniform on `[-1, 1]`.
   Each value comes from a hash of `(s_i, r, c)`, so it does not depend on thread
   scheduling.
2. Take the complex forward FFT.
3. Multiply each coefficient by `f^(-β/2)`, where `f` is the radial frequency of
   that coefficient, `β` is the background exponent from section 2.8, and the DC
   coefficient is set to zero. The power spectrum of the synthetic field is then
   proportional to `f^(-β)`.
4. Take the inverse FFT, divide by `P_r · P_c`, and crop to `n_rows × n_cols`.
5. Normalise to zero mean and to standard deviation `σ`, the detrended standard
   deviation of the observed tile.
6. Compute the periodogram of the synthetic surface with the same padding, window
   and normalisation as the observed tile (section 2.4 to 2.6).

The background 2D periodogram `P̄(k, l)` is the mean of the `N` synthetic 2D
periodograms, and the background radial spectrum is the mean of their radial
values. Both are averages over the same cell or radial index.

#### Why average over iterations

Each periodogram cell of a Gaussian-like field is, to a good approximation,
exponentially distributed about its expected value. A single synthetic realisation
is therefore a noisy estimate of the expectation. If one realisation were used as
the denominator, the ratio would have the distribution of the ratio of two
independent exponential variables. That distribution has a much heavier upper tail
than the exponential, and it gives many false detections (section 2.12).

Averaging `N` realisations reduces the relative scatter of the background estimate
in each cell to about `1/√N`. That is roughly 10 % for the default of 100
iterations, and 22 % for 20 iterations (the CLI default before 2026-10-21). The scatter does not vanish, so the
significance threshold depends on `N` as well (section 2.12).

### 2.10 Normalised PSD

The normalised PSD is the ratio of the observed spectrum to the background
spectrum:

```
S(k, l) = P(k, l) / P̄(k, l)                       (2D cells)
S(f)    = radial_PSD(f) / background_PSD(f)       (radial points)
S_bin   = mean of S over the radial points in a log bin
```

For the binned product, a radial point is kept only if its frequency, observed
power, background power and ratio are all positive and finite. The three binned
series therefore share the same frequency grid.

### 2.11 Significance tests

Two tests are applied. They answer different questions and use different
thresholds.

#### (a) Binned radial spectrum

`_significant_frequencies.csv` (section 6) reports the log bins whose mean
normalised PSD exceeds the constant

```
S_bin > 2.9955
```

The constant is written in the CLI as `0.5 × 5.991`, which is the chi-squared value
for 95 % confidence with 2 degrees of freedom, halved. The exact value of
`-ln(1 - 0.95)` is 2.99573. The difference is in the fourth significant figure.

This test is applied to bin means of the per-point ratios, not to individual
cells, so the 95 % label does not give the false-alarm rate of this test. Under
the idealised null with independent cells, the mean of `n` unit-mean exponential
values has standard deviation `1/√n`. A bin with many points therefore exceeds
2.9955 by chance less often than 5 % of the time. Adjacent cells are correlated by
the window, so the effective count per bin is smaller than the number of cells.
The 95 % label records where the constant came from.

#### (b) Per-cell 2D test

!!! note "Where 9.21 comes from"

    9.21 is the value a single cell's normalised PSD exceeds by chance with
    probability 1 in 10,000 when the surface has no characteristic scale.
    For such a surface each periodogram cell is the squared magnitude of a
    complex Fourier coefficient whose real and imaginary parts are independent
    Gaussians of equal variance. A sum of two squared standard Gaussians is
    χ² with 2 degrees of freedom; dividing by its mean (2) gives an exponential
    distribution with mean 1. So the normalised PSD `S` of a cell satisfies

        P(S > t) = exp(-t)

    Setting this chance to `1 - p` gives the threshold `t = -ln(1 - p)`. At the
    default confidence `p = 0.9999`, `t = -ln(10⁻⁴) = 4 ln 10 = 9.2103`. The
    legacy LSDTopoTools code wrote the same number as half the χ²₂ value at
    99.99 %, `18.42185 / 2 = 9.21093`. The binned test's 2.9955 is the same
    construction at 95 %: `5.991 / 2`, or `-ln(0.05) = 2.9957`. The derivation
    and the correction for a background averaged from a finite number of
    synthetic surfaces follow below.

`_significant_frequencies_2d.csv` (section 6) reports each 2D cell whose ratio
exceeds

```
S(k, l) > -ln(1 - p)
```

where `p` is `significant_frequencies_2d_confidence` (default 0.9999). The
library function `spectral_normalised_psd_threshold` evaluates this threshold.
Values are 0.693147 at p = 0.5, 2.99573 at p = 0.95, and 9.21034 at p = 0.9999.
The legacy LSDTopoTools constant was 18.42185 / 2 = 9.21093 for the same 99.99 %
level.

The threshold follows from the distribution of a single cell. Under the null
hypothesis that the surface is a fractal background with expectation `B`, the
periodogram value is a sum of squares of two independent Gaussian components of
equal variance. So `P / B` has a chi-squared distribution with 2 degrees of
freedom, divided by 2:

```
P / B  ~  χ²₂ / 2  =  Exp(1)       (mean 1, CDF 1 - exp(-x))
```

The CDF of this distribution is `1 - exp(-x)`. Setting it equal to `p` gives the
threshold `x = -ln(1 - p)`. The code cites Perron et al. (2008), equation 9, for
this relation.

#### Finite background iterations

The derivation above treats the background `B` as known. In practice `B` is the
mean of `N` synthetic periodograms. Each is independent and exponential with mean
`B`, so `N·B̂/B` is Gamma-distributed with shape `N`. The ratio is then

```
S = P / B̂  ~  F(2, 2N)
```

an F distribution with 2 and `2N` degrees of freedom. The quantile at confidence
`p` is

```
threshold(p, N) = N · [ (1 - p)^(-1/N) - 1 ]
```

which tends to `-ln(1 - p)` as `N` grows. At p = 0.9999:

| `N` (iterations) | threshold |
|---|---|
| 1 | 9999 |
| 5 | 26.55 |
| 20 | 11.70 |
| 100 (default) | 9.65 |
| infinite | 9.21 |

These values come from the formula above. A Monte Carlo check with 2 × 10⁶ draws
reproduced them to within about 2 % (the tail at p = 0.9999 contains only about 200
draws, so the sampling error is of that order). In that idealised model, a threshold of
9.21 gives a per-cell false-alarm rate of `(1 + 9.21/N)^(-N)`: about 5 × 10⁻⁴ at
`N = 20`, five times the nominal 10⁻⁴, and about 1.5 × 10⁻⁴ at the default
`N = 100`. The default was raised from 20 to 100 for this reason. The reported
threshold stays at `-ln(1 - p)`, so it matches the legacy output. The model ignores the window correlations between cells, any
departure from Gaussian statistics, and the mismatch between the fitted power law
and the real spectrum. The actual false-alarm rate on a real tile may differ.

### 2.12 Multiple testing

The 2D test is applied to every cell, so many cells are tested. Under the idealised
null, the expected number of 2D entries is about `(P_r · P_c - 1) · (1 - p)`. This
count includes both members of each mirror pair (section 3), so the number of
independent chance passes is about half of it.

For the default `p = 0.9999` on a 128 × 128 padded grid, that expectation is about
1.6 entries per tile in the infinite-`N` limit, and about 2.4 with the default
`N = 100` (false-alarm rate 1.5 × 10⁻⁴, section 2.11). At `p = 0.95` it is about 800 entries, so most passes at that
level are expected by chance under the idealised null.

---

## 3. Bearings

Each significant 2D cell is reported with the direction of its wave vector.
The convention is clockwise from north:

```
θ = atan2(f_east, f_north)   reduced to [0, 2π)     (bearing_rad)
θ_deg = θ · 180/π                                    (bearing_deg, [0, 360))
f     = sqrt(f_east² + f_north²)                     (frequency_cycles_per_m)
λ     = 1 / f                                        (wavelength_m)
```

The bearing is the direction of the wave vector, the direction in which the phase
of the wave increases. The crests of the corresponding wave run perpendicular to
the wave vector, so a crest line lies at `θ + 90°` (equivalently `θ - 90°`). For
an orientation analysis that ignores direction, use `(θ + 90°) mod 180°` or
`θ mod 180°`, as appropriate.

The bearing uses the physical frequency components, `col_offset · dfx` for east
and `row_offset · dfy` for north. That is correct for a padded grid whose row and
column counts differ. The legacy code used index offsets, which agree only when
`dfx = dfy`.

### Point symmetry

The spectrum of a real surface is point-symmetric, so each significant frequency
normally appears twice. The two entries sit at `θ` and `θ + 180°`. Both have the
same frequency, and both have the same normalised value, because the observed and
background spectra are both symmetric. The DC cell is never reported.

Cells on the Nyquist edges do not follow this rule. The point-symmetric partner of
a cell is found by negating its offset from DC. On the Nyquist row (row 0, where
`f_north = +f_N`) and on the Nyquist column (column 0, where `f_east = -f_N`), the
negated offset wraps back onto the same edge, because `+f_N` and `-f_N` are the
same frequency in the discrete spectrum. The partner of a cell on row 0 is
therefore another cell on row 0, at the mirrored east component. Its bearing is
`360° - θ`, not `θ + 180°`. Similarly, a cell on column 0 pairs with another cell
on column 0 at the mirrored north component, at bearing `180° - θ`.

Three cells are their own partners: the corner at row 0, column 0 (bearing 315°),
the cell at row 0 in the DC column (bearing 0°), and the cell at column 0 in the DC
row (bearing 270°). The values at the edge pairs are equal, as they are for all
partners. Users who pair entries should allow for these edge cases.

### Orientation of the 2D rasters

The 2D PSD rasters are written north-up as they are stored (section 6). Row 0 is
the north edge, and the frequency axes run in the same directions as the bearings
above.

---

## 4. Reproducibility

The background is the only random part of the method. Reproducibility is handled
as follows.

- `background_seed = -1` (default) uses the built-in base seed of 1000. An
  explicit non-negative value is the base seed. Iteration `i` uses `base + i`.
- The effective base is written to the 2D CSV header (`# background_seed`). When
  the parameter is `-1`, the header shows 1000.
- The white-noise values are a hash of `(seed, row, col)`, so each value is
  independent of thread scheduling. The normalisation step is serial, and the loop
  over iterations is sequential, so the accumulation order is fixed.
- The observed periodogram does not depend on the seed or on the thread count.

The unit test `unit_spectral_analysis_significant_2d_seed` checks, on one build
and in one process:

- `background_psd_2d` is bit-identical at 1 and 4 threads for seed 42;
- the binned normalised PSD and the significant 2D cell list are identical at
  1 and 4 threads;
- seed 42 and seed 43 give different backgrounds;
- the default seed `-1` gives the same background as seed 1000;
- the observed `psd_2d` does not depend on the seed.

The test does not establish bit-identity across machines or FFTW builds. The
observed periodogram is planned with `FFTW_MEASURE`, which may choose different
algorithms on different machines, and different algorithms can round differently.

The legacy LSDTopoTools2 implementation did not have this property. Its background
iterations filled each random field through `rewrite_with_random_values`, which
seeds with `time(NULL)`. Legacy background spectra therefore changed from run to
run, and the significance result could change with them.

---

## 5. Test surfaces

The spectral code has three sources of test surfaces. Two are generated inside the
algorithm or by the project's synthetic-surface tools, and one is a set of
inline sinusoids in the unit test.

### 5.1 Backgrounds inside the algorithm

The Monte Carlo backgrounds in section 2.9 are the only synthetic surfaces that the
analysis itself generates. They are power-law fields with exponent `β`, normalised
to the detrended standard deviation of the observed tile. They are not
test fixtures in their own right. Their role is to define the expected spectrum.

### 5.2 Diamond-square fixture (`sb_diamond_square.tif`)

The spectral smoke tests use the DEM written by the stream-burning fixture. The
fixture parameter file is `tests/params/stream_burning_fixture_dem.param.in`. It
calls `lsdtt-synthetic-surfaces` with:

- `write_diamond_square_surface: true`;
- a 200 × 200 raster at 30 m cell size;
- `random_seed: 33`;
- `diamond_square_roughness: 0.6`;
- `diamond_square_stddev: 120.0`.

Diamond-square is a fractal generator. It builds the surface on the smallest
`2^k + 1` square that contains the requested raster, then crops it. For a 200 × 200
request that is a 257 × 257 grid. The surface has no characteristic wavelength,
so it is a useful negative control. The expected result is no significant cells at
the default 0.9999 level, and some significant cells at 0.5.

The fixture requires the Python test tier, because the stream-burning fixture is
written by a Python script.

### 5.3 Known-answer sinusoids

The `--bearings` mode of `lsdtt-test-spectral-significant-2d` builds its own
test surface. The function `make_wave_dem` writes a plane wave of unit amplitude
with a constant offset of 100, plus Gaussian noise with standard deviation 0.2
(bearings) or 0.5 (seed test) from a fixed generator seed. The cell size is 1 m,
with the synthetic georeference set by `set_default_synthetic_georeference`.
For a wave with `kx` cycles across the columns and `ky` cycles across the rows,
north positive:

```
z(r, c) = 100 + sin( 2π · (kx · c / n_cols  -  ky · r / n_rows) ) + noise
```

Four cases are checked:

| Case | Grid (rows × cols) | `kx`, `ky` | Expected bearing |
|---|---|---|---|
| north | 128 × 128 | 0, 8 | 0° (and 180°) |
| east | 128 × 128 | 8, 0 | 90° (and 270°) |
| north-east | 128 × 128 | 8, 8 | 45° (and 225°) |
| north-east, non-square | 64 × 128 | 8, 4 | 45° (and 225°) |

On the 64 × 128 grid, `dfy = 2 dfx`. The wave has equal physical components of
1/16 cycles/m, so the bearing is 45°. The index-offset bearing used by the legacy
code would give atan2(8, 4) = 63.43° on that grid, which is why the case is
included.

### 5.4 What each test asserts

| CTest name | Mode | What it checks |
|---|---|---|
| `unit_spectral_analysis_significant_2d_threshold` | `--threshold` | `spectral_normalised_psd_threshold(p)` equals `-ln(1 - p)` for p in {0.5, 0.95, 0.99, 0.9999} to 1e-12 relative; it is within 1e-3 of the legacy value 9.210925 at 0.9999; it throws `std::invalid_argument` for p = 0, 1, -0.1, 1.5 and NaN. Needs no FFTW. |
| `unit_spectral_analysis_significant_2d_bearings` | `--bearings` | For each of the four sinusoid cases: some cells are significant at 0.9999; the two strongest cells sit at the expected bearing and bearing + 180 (within 1e-6°); their frequency is the wave's frequency (within 1e-9 cycles/m); `wavelength × frequency = 1`; `bearing_rad` and `bearing_deg` agree and lie in [0, 2π); every significant interior cell has a partner at bearing + 180 with the same frequency. Nyquist-edge cells are excluded from the partner check, consistent with section 3. Needs FFTW. |
| `unit_spectral_analysis_significant_2d_seed` | `--seed` | The seed behaviour in section 4, on a 64 × 64 surface with 5 iterations. Needs FFTW. |
| `smoke_spectral_analysis_significant_2d` | CLI | On the diamond-square fixture at seed 17 with 5 iterations: the 2D CSV has the expected header and `# normalised_psd_threshold,9.21034` and `# background_seed,17`; `n_significant_frequencies` is 0; both the 2D and the 1D CSVs contain the "No significant frequencies" line. |
| `smoke_spectral_analysis_significant_2d_found` | CLI | The same DEM at `significant_frequencies_2d_confidence: 0.5`: the threshold is 0.693147; the 2D CSV has at least one cell; the "No significant frequencies" line is absent. |
| `smoke_spectral_analysis_valley_spacing` | CLI | The valley-spacing path runs on the stream-burning fixture, writes the two valley-spacing CSVs, and prints the valley-spacing banner and a hillslope-length line of 30 m. |
| `regression_spectral_analysis_valley_spacing_legacy_alias` | CLI | The deprecated key `valley_spacing_m` maps to `valley_spacing_area_exponent` and warns. |

The `--bearings` and `--seed` modes are registered only when FFTW is available.
The smoke tests that use the diamond-square fixture depend on the stream-burning
fixture, and therefore on the Python test tier. The legacy-alias test does not use
the fixture.

---

## 6. Outputs

All outputs are written to `write_path` with the prefix `write_fname`. Names below
use `<write_fname>` for that prefix. Raster outputs use the extension selected by
`write_output_format`: `.tif` for `cog`, `gtiff`, and `same_as_input` (which takes
the input extension), and `.png` for `png`.

| File | Written when (default) | Contents |
|---|---|---|
| `<write_fname>_psd_2d.<ext>` | `write_2d_spectrum_rasters` (true) | Observed 2D PSD, padded size, float, fftshift order, m² per cell for an elevation input. |
| `<write_fname>_background_psd_2d.<ext>` | `write_2d_spectrum_rasters` (true) | Background 2D PSD, same grid and units. |
| `<write_fname>_radial_spectrum.csv` | `write_radial_spectrum` (false) | Log-binned radial spectrum. Columns below. |
| `<write_fname>_significant_frequencies.csv` | `write_significant_frequencies` (true) | Rollover summary, the binned threshold and the significant bins. Columns below. |
| `<write_fname>_significant_frequencies_2d.csv` | `write_significant_frequencies_2d` (true) | Significant 2D cells with bearings. Columns below. |
| `<write_fname>_valley_spacing.csv` | `run_valley_spacing` (false) | Valley-spacing measures, with `quantity,value,units` columns. |
| `<write_fname>_valley_spacing_spectrum.csv` | `run_valley_spacing` (false) | The averaged valley-spacing spectrum. |
| `<write_fname>_ingested_params.param` and `<write_fname>_README` | Always | A record of the parameters for this run. |

### 2D PSD rasters

The rasters are stored in fftshift order at the padded size `P_r × P_c`. Their
georeferencing puts the DC pixel centre at coordinates (0, 0) in cycles/m. The
pixel centre at column `c` has east frequency `(c - P_c/2) · dfx`, and the pixel
centre at row `r` has north frequency `(P_r/2 - r) · dfy`. Row 0 is at the top, so
the rasters are north-up. Values are float32.

The pixel width is `dfx` and the pixel height `dfy`. They differ when the padded
row and column counts differ, which happens whenever the input's rows and
columns round up to different powers of two (for example 200 × 500 pads to
256 × 512, so `dfy = 2 dfx`), even though the input pixels are square. The
rasters are then written with non-square pixels, so both axes read in true
cycles/m in a GIS or in lsdviztools3. These are the only lsdtt3 outputs with
non-square pixels; lsdtt3 itself refuses non-square pixels on input, so they
cannot be fed back into another lsdtt3 tool.

### `_radial_spectrum.csv`

| Column | Units | Meaning |
|---|---|---|
| `frequency_cycles_per_m` | cycles/m | Mean frequency of the bin. |
| `wavelength_m` | m | 1 / frequency. |
| `psd` | m² | Mean of the binned radial PSD values (section 2.6). |
| `background_psd` | m² | Mean of the binned background values. |
| `normalised_psd` | - | Mean of the per-point ratios in the bin (section 2.10). |
| `beta_by_bin` | - | Suffix-regression exponent for that bin (section 2.8). |
| `r_sq_by_bin` | - | R² of that suffix regression. |

The `beta_by_bin` and `r_sq_by_bin` columns come from the rollover binning of the
observed spectrum. The rows come from the shared-valid binning (section 2.10),
which can drop a bin if none of its points has a positive background or ratio.
The two columns are matched to the rows by log bin (`spectral_align_bin_trace`),
so a dropped bin leaves the other rows' values in place. A row whose bin has no
rollover fit (the highest-frequency fifth of the bins, which the suffix
regression does not start from) holds `nan`.

### `_significant_frequencies.csv`

Header comment lines, in order:

```
# rollover_frequency_cycles_per_m,<value>
# rollover_wavelength_m,<value or inf>
# rollover_beta,<value>
# sub_rollover_beta,<value used for the background>
# normalised_psd_threshold,2.9955
# n_significant_frequencies,<count>
```

Then the column header and one row per significant bin:

| Column | Units | Meaning |
|---|---|---|
| `frequency_cycles_per_m` | cycles/m | Mean frequency of the bin. |
| `wavelength_m` | m | 1 / frequency. |
| `normalised_psd` | - | Mean per-point ratio in the bin. |

If no bin passes, the file has no data rows. A comment line follows the column
header:

```
# No significant frequencies were found in this raster: no bin of the radial spectrum exceeded the normalised PSD threshold of 2.9955 (95% confidence). This is a valid result, not an error.
```

### `_significant_frequencies_2d.csv`

Header comment lines, in order:

```
# confidence,<p>
# normalised_psd_threshold,<-ln(1-p)>
# background_seed,<effective base seed>
# n_background_iterations,<N>
# n_significant_frequencies,<count>
```

Then the column header and one row per significant cell:

```
bearing_rad,bearing_deg,frequency_cycles_per_m,wavelength_m,normalised_psd
```

| Column | Units | Meaning |
|---|---|---|
| `bearing_rad` | rad | Direction of the wave vector, clockwise from north, in [0, 2π). |
| `bearing_deg` | ° | The same direction, in [0, 360). |
| `frequency_cycles_per_m` | cycles/m | Radial frequency of the cell. |
| `wavelength_m` | m | 1 / frequency. |
| `normalised_psd` | - | Observed PSD divided by background PSD for that cell. |

If no cell passes, a comment line follows the column header:

```
# No significant frequencies were found in this raster: no cell of the 2D spectrum exceeded the normalised PSD threshold of <threshold> (confidence <p>). This is a valid result, not an error.
```

Read either CSV with `pandas.read_csv(path, comment="#")`.

### Valley-spacing outputs

These are written only when `run_valley_spacing` is true. The valley-spacing
measures are described in the CLI reference, not on this page.

`<write_fname>_valley_spacing.csv` has `quantity,value,units` columns, with these
rows: `hillslope_length`, `hillslope_length_mean`, `hillslope_length_median`,
`hillslope_length_p90` (all m); `n_hillslope_cells`, `n_convex_runs`,
`n_spectral_peaks`, `valleys_in_domain` (count); `convex_halfwidth`,
`convex_wavelength`, `spectral_wavelength_power_weighted`, `spectral_wavelength_peak`,
`domain_length` (m); `spectral_exponent_alpha`, `spectral_exponent_predicted`,
`channelization_index` (-); and `spectrum_under_sampled` (bool, 0 or 1).

`<write_fname>_valley_spacing_spectrum.csv` has `frequency_cycles_per_m,wavelength_m,power`
columns. It starts at the first non-zero frequency, so the DC term is excluded. The
source does not state the units of `power`.

### Standard output

The run prints the rollover frequency, the rollover wavelength, the rollover
exponent and the sub-rollover exponent. The last two are the values in the outputs
above.

### Plotting

The `lsdviztools3` package provides plotting for these outputs. Its loader is
`load_spectral_analysis`, and its plotting functions are
`plot_significant_frequencies_polar` and `plot_spectral_tile_summary`. This page
does not describe their interfaces.

---

## 7. Practical guidance and limitations

These points describe how the method behaves. They are not recommendations for
particular thresholds.

**The tile must contain several wavelengths of the pattern.** A wavelength `λ`
appears at about `(tile width) / λ` cycles across the tile. A pattern that repeats
only once or twice within the tile puts its power into the first few bins, where it
is hard to separate from the background and from the planar trend.

**Resolution and the Nyquist limit.** The Nyquist frequency is `1 / (2 · cell_size)`.
Wavelengths shorter than two cells cannot be represented. Along the diagonals of
the square spectrum, cells reach `√2` times the Nyquist frequency, that is a
wavelength of about `1.41` cells. Those diagonal cells appear in the 2D output but not
in the radial spectrum (section 2.6).

**Edges and windowing.** The Hann window tapers to zero at the inscribed ellipse,
so the corners of a square tile carry no weight and the tile's effective size is
smaller than its nominal size. The window reduces leakage from the edges, but it
also reduces the effective number of samples.

**Padding.** Each dimension is padded to a power of two, so a tile can be padded
by almost a factor of two in each direction. The padding is zero-filled. Because
the window tapers to zero at the tile edge, the periodic copies of the padded field
are separated by zero-valued gaps. The spectrum is still sampled on the padded grid,
so the cell spacing is set by the padded size, not the tile size.

**Nodata.** The command-line tool does not expose the valid-rectangle clipping that
the library provides. Nodata pixels are zero-filled in the observed periodogram.
The background surfaces have no holes. The observed and background spectra
therefore differ wherever the tile has holes. Supply a tile without holes, or crop
it to a rectangle of valid pixels before running the analysis.

**Small tiles and rollover at a low bin index.** The number of bins is roughly
`log10( max(P_r, P_c) / 2 ) / spectral_log_bin_width` (section 2.7). A small tile
or a wide bin gives few bins, and the rollover then rests on few points.

- With fewer than three binned points, the rollover search cannot run. The code
  throws `insufficient binned data for rollover fit`. The command-line tool then
  prints "Error running spectral analysis" and exits with a non-zero status.
- The rollover index is never the first bin when there is more than one bin
  (section 2.8).
- With at least three bins and a rollover index of at least one, the sub-rollover
  fit has at least two points, so the code never fits a single point. A rollover
  index of 1 or 2 gives a sub-rollover fit over two or three bins.
- A fit through exactly two points is exact, so its R² carries no information about
  scatter. A sub-rollover exponent from two bins is a rough value.
- The fallbacks in section 2.8 apply only when the fitted exponent is non-positive or
  not finite. The value in the outputs is the one used, so check `sub_rollover_beta`
  against the rollover exponent before interpreting a result from a small tile.

Earlier code crashed on small tiles when the sub-rollover fit had a single point.
lsdtt3 cannot form that fit, as described above.

**Header-only and no-significant results are genuine.** A result with no significant
cells is a valid outcome. The "No significant frequencies" lines are written so that
an empty table is not mistaken for a failed run. A tile with a featureless surface
should produce that outcome.

**Multiple testing.** The per-cell test applies a threshold at every cell. Section
2.12 gives the expected number of chance passes for the idealised null. The finite
number of background iterations raises the false-alarm rate above that idealised
value (section 2.11). For a sensitive search, use a threshold derived from the
actual `N` (the F-distribution quantile in section 2.11), or repeat the analysis
with a different background seed and check that the detections persist.

**Non-elevation inputs.** The planar detrend is applied to whatever the raster
holds. For a density raster, a large-scale gradient across the tile (for example,
a moisture gradient) is removed as a plane, and any remaining large-scale structure
is treated as background. Check that the input is in projected metres, and that the
pattern is visible in the raster before drawing conclusions about its spacing.

**Trend and curvature.** Only a plane is removed. Curvature in the tile leaks into
the low-frequency bins, and that can affect the rollover and the background
exponent.

---

## 8. Differences from the legacy code and known caveats

The items below are behaviour of the current code that differs from the legacy
LSDTopoTools implementation, or that a user could reasonably misread.

- **Background seed.** The legacy background was time-seeded. The current code is
  deterministic (section 4).
- **Bearing.** The legacy code used index offsets for the bearing. The current code
  uses physical frequency components, which differ on non-square padded grids.
- **Edge pairing.** Cells on the Nyquist row and column pair with cells on the same
  edge, not at `θ + 180°` (section 3). Three cells are their own partners.
- **Threshold constant.** The legacy 99.99 % constant was `18.42185 / 2 = 9.21093`.
  The current 2D test uses `-ln(1 - 0.9999) = 9.21034`. The difference is small.
- **Binned threshold.** The binned test uses the constant 2.9955 on bin means
  (section 2.11). It is not applied to individual cells.
- **Corner cells.** The 2D output includes corner cells outside the Nyquist circle,
  which the radial spectrum omits (section 2.6).
- **Default iterations and bin width.** `n_background_iterations` defaults to 100
  in both the CLI and the library (the CLI default was 20 before 2026-10-21). The
  library default bin width is 0.1, and the CLI default is 0.05. The CLI values are
  the ones used in a command-line run.
- **PSD units.** The PSD is the variance per DFT cell (section 2.4). No division by
  the frequency spacing is applied.
- **Finite-N threshold.** The per-cell threshold `-ln(1 - p)` is the limit for an
  infinite number of background iterations. For finite `N`, the F-distribution
  quantile in section 2.11 applies.
- **Non-square 2D rasters.** Before 2026-10-21 the 2D spectrum rasters used `dfx`
  for both axes, so when the padded grid was not square their north axis was
  scaled by `dfx / dfy`. They now carry a separate pixel height (section 6).
- **Radial bin alignment.** Before 2026-10-21 the `beta_by_bin` and `r_sq_by_bin`
  columns of `_radial_spectrum.csv` were copied to the rows by position, so a bin
  dropped by the shared-valid binning shifted every later value by one row. They
  are now matched by log bin (section 6).

---

## 9. Relation to the published method

The method is the spectral approach of Perron, Kirchner and Dietrich (2008): the
2D periodogram of a detrended, windowed DEM, normalised by the spectrum of
synthetic fractal surfaces, with significance judged from the chi-squared
distribution of periodogram values. The ridge-and-valley context is Perron et al.
(2009). The equations on this page describe the lsdtt3 implementation. Where the
code cites an equation of the 2008 paper (the per-cell threshold, cited as its
equation 9), the citation is kept as the code gives it; the page has not been
checked equation by equation against the paper.

**Perron, J. T., Kirchner, J. W., & Dietrich, W. E. (2008).** Spectral signatures of
characteristic spatial scales and nonfractal structure in landscapes. *Journal of
Geophysical Research: Earth Surface*, 113, F04003. doi:10.1029/2007JF000866

**Perron, J. T., Kirchner, J. W., & Dietrich, W. E. (2009).** Formation of evenly
spaced ridges and valleys. *Nature*, 460, 502–505. doi:10.1038/nature08174

The valley-spacing measures are described in the CLI reference, and they draw on
Bonetti et al. (2020) and Hooshyar et al. (2021), as cited in the valley-spacing
module.

---

## See also

- [Spectral analysis CLI reference](../cli_reference/cli_reference_lsdtt_spectral_analysis.md)
- [Synthetic surfaces: choosing a tool](synthetic_cli_cluster.md)
- [Surface metrics](surface_metrics.md), for the wavelet quadratic fit that shares the spectral module
