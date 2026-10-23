# Cosmogenic Nuclides API

This page is generated from the built Python runtime manifest. Do not edit signatures here.

## `CRNBurialResult`

Result of :func:`crn_burial_pair`: pre-burial erosion rate (g/cm^2/yr),
burial age (years) and ratio diagnostics; ``status`` is ``"Converged"``,
``"NoSolution"``, ``"NotUnique"`` or ``"NotConverged"``.

Returned by ``crn_burial_pair()``; not constructible directly.

Attributes
----------
apparent_erosion_rate_g_cm2_yr : float
burial_age_yr : float
burial_required : bool
measured_ratio : float
nuclide_fast : str
nuclide_slow : str
preburial_erosion_rate_g_cm2_yr : float
status : str
steady_state_ratio : float

Authors
-------
B. Gailleton

### Members

#### `apparent_erosion_rate_g_cm2_yr`

(self) -> float

#### `burial_age_yr`

(self) -> float

#### `burial_required`

(self) -> bool

#### `measured_ratio`

(self) -> float

#### `nuclide_fast`

(self) -> str

#### `nuclide_slow`

(self) -> str

#### `preburial_erosion_rate_g_cm2_yr`

(self) -> float

#### `status`

(self) -> str

#### `steady_state_ratio`

(self) -> float

---

## `CRNColumn`

One-dimensional cosmogenic-nuclide column model.

Set an erosion history, ordered from present backwards in time, then use
``transient`` to model depth profiles. Point steady-state and step-change
calculations require no history. Erosion rates for these point calculations
are effective rates in g/cm²/yr.

The history (``set_erosion_history`` or ``read_erosion_history``) has rows
``time`` (years before present at which the rate changes, increasing),
``rate`` (erosion rate in m/yr; 0 means an instantaneous removal event) and
``removal`` (metres removed, used only when rate is 0); the last row must be
``time = -9999`` with the positive long-term background rate.

Construct it as

    CRNColumn(rock_density_kg_m3: float = 2650.0)

Parameters
----------
rock_density_kg_m3 : float, default=2650.0
    Rock density (kg/m^3, positive) converting depths in m to effective
    depths in g/cm^2 (``depth_m * density * 0.1``).

Attributes
----------
depth_change_m : float
    Total depth (m) removed over the history: added to sampling depths to get their starting depths.
erosion_rate_m_yr : list[float]
    History erosion rates (m/yr).
removal_m : list[float]
    History removal depths (m).
rock_density_kg_m3 : float
    Rock density (kg/m^3).
time_yr_bp : list[float]
    History times (years before present).

Methods
-------
read_erosion_history(path: str) -> None
    Read the erosion history from a CSV file with header ``time,rate,removal``.
set_erosion_history(time_yr_bp: Sequence[float], erosion_rate_m_yr: Sequence[float], removal_metres: Sequence[float]) -> None
    Set the history. Last row must be ``time=-9999`` with a positive background rate.
steady_state(effective_depth_g_cm2: float = 0.0, effective_erosion_rate_g_cm2_yr: float = 0.01, production_scaling: float = 1.0, total_shielding: float = 1.0, muon_scheme: str = 'BraucherBorchers') -> CRNColumnSample
    Return steady-state concentrations at one effective depth.
step_change(effective_depth_g_cm2: float = 0.0, effective_erosion_rate_old_g_cm2_yr: float = 0.01, effective_erosion_rate_new_g_cm2_yr: float = 0.01, time_since_step_yr: float = 0.0, production_scaling: float = 1.0, total_shielding: float = 1.0, muon_scheme: str = 'BraucherBorchers') -> CRNColumnSample
    Return concentrations after a step change in erosion rate.
transient(sampling_depths_metres: Sequence[float], production_scaling: float = 1.0, total_shielding: float = 1.0, muon_scheme: str = 'BraucherBorchers', initial_10be: float = -1.0, initial_26al: float = -1.0, initial_14c: float = -1.0) -> list[CRNColumnSample]
    Run the erosion history and return concentrations at the given depths.

Authors
-------
B. Gailleton
S. M. Mudd

**Deprecated compatibility names:**
- `removal_m` -> `removal_metres`; earliest removal: following minor.
- `sampling_depths_m` -> `sampling_depths_metres`; earliest removal: following minor.

```python
def __init__(self, rock_density_kg_m3: float = 2650.0) -> None
```

### Members

#### `__init__`

__init__(self, rock_density_kg_m3: float = 2650.0) -> None

```python
def __init__(self, rock_density_kg_m3: float = 2650.0) -> None
```

#### `depth_change_m`

Total depth (m) removed over the history: added to sampling depths to get their starting depths.

#### `erosion_rate_m_yr`

History erosion rates (m/yr).

#### `read_erosion_history`

read_erosion_history(self, path: str) -> None

Read the erosion history from a CSV file with header ``time,rate,removal``.

Parameters
----------
path : str
    CSV file; columns as in ``set_erosion_history`` (years BP, m/yr, m).

Returns
-------
None

```python
def read_erosion_history(self, path: str) -> None
```

#### `removal_m`

History removal depths (m).

#### `rock_density_kg_m3`

Rock density (kg/m^3).

#### `set_erosion_history`

set_erosion_history(self, time_yr_bp: collections.abc.Sequence[float], erosion_rate_m_yr: collections.abc.Sequence[float], removal_metres: collections.abc.Sequence[float]) -> None

Set the history. Last row must be ``time=-9999`` with a positive background rate.

Parameters
----------
time_yr_bp : list of float
    Times (years before present) of each rate change, increasing; the last
    is -9999.
erosion_rate_m_yr : list of float
    Erosion rate (m/yr) from each time; 0 marks an instantaneous removal.
removal_metres : list of float
    Depth removed (m) at rows with rate 0, else 0.

Returns
-------
None

```python
def set_erosion_history(self, time_yr_bp: collections.abc.Sequence[float], erosion_rate_m_yr: collections.abc.Sequence[float], removal_metres: collections.abc.Sequence[float]) -> None
```

#### `steady_state`

steady_state(self, effective_depth_g_cm2: float = 0.0, effective_erosion_rate_g_cm2_yr: float = 0.01, production_scaling: float = 1.0, total_shielding: float = 1.0, muon_scheme: str = 'BraucherBorchers') -> lsdtt3.CRNColumnSample

Return steady-state concentrations at one effective depth.

Parameters
----------
effective_depth_g_cm2 : float, default=0.0
    Sample depth as effective depth (g/cm^2); 0 is the surface.
effective_erosion_rate_g_cm2_yr : float, default=0.01
    Steady erosion rate (g/cm^2/yr).
production_scaling : float, default=1.0
    Production-scaling factor (1 at sea level, high latitude).
total_shielding : float, default=1.0
    Combined shielding factor (0-1; 1 = unshielded).
muon_scheme : str, default="BraucherBorchers"
    Muon production parameterisation (``"Braucher"``, ``"BraucherBorchers"``,
    ``"Granger"``, ``"Schaller"`` or ``"newCRONUS"``).

Returns
-------
CRNColumnSample
    Be-10, Al-26 and C-14 concentrations (atoms/g).

```python
def steady_state(self, effective_depth_g_cm2: float = 0.0, effective_erosion_rate_g_cm2_yr: float = 0.01, production_scaling: float = 1.0, total_shielding: float = 1.0, muon_scheme: str = 'BraucherBorchers') -> lsdtt3.CRNColumnSample
```

#### `step_change`

step_change(self, effective_depth_g_cm2: float = 0.0, effective_erosion_rate_old_g_cm2_yr: float = 0.01, effective_erosion_rate_new_g_cm2_yr: float = 0.01, time_since_step_yr: float = 0.0, production_scaling: float = 1.0, total_shielding: float = 1.0, muon_scheme: str = 'BraucherBorchers') -> lsdtt3.CRNColumnSample

Return concentrations after a step change in erosion rate.

Parameters
----------
effective_depth_g_cm2 : float, default=0.0
    Sample depth as effective depth (g/cm^2).
effective_erosion_rate_old_g_cm2_yr, effective_erosion_rate_new_g_cm2_yr : float, default=0.01
    Erosion rate before and after the change (g/cm^2/yr).
time_since_step_yr : float, default=0.0
    Years since the change.
production_scaling : float, default=1.0
    Production-scaling factor (1 at sea level, high latitude).
total_shielding : float, default=1.0
    Combined shielding factor (0-1; 1 = unshielded).
muon_scheme : str, default="BraucherBorchers"
    Muon production parameterisation (``"Braucher"``, ``"BraucherBorchers"``,
    ``"Granger"``, ``"Schaller"`` or ``"newCRONUS"``).

Returns
-------
CRNColumnSample

```python
def step_change(self, effective_depth_g_cm2: float = 0.0, effective_erosion_rate_old_g_cm2_yr: float = 0.01, effective_erosion_rate_new_g_cm2_yr: float = 0.01, time_since_step_yr: float = 0.0, production_scaling: float = 1.0, total_shielding: float = 1.0, muon_scheme: str = 'BraucherBorchers') -> lsdtt3.CRNColumnSample
```

#### `time_yr_bp`

History times (years before present).

#### `transient`

transient(self, sampling_depths_metres: collections.abc.Sequence[float], production_scaling: float = 1.0, total_shielding: float = 1.0, muon_scheme: str = 'BraucherBorchers', initial_10be: float = -1.0, initial_26al: float = -1.0, initial_14c: float = -1.0) -> list[lsdtt3.CRNColumnSample]

Run the erosion history and return concentrations at the given depths.

Parameters
----------
sampling_depths_metres : list of float
    Depths (m) below the present surface.
production_scaling : float, default=1.0
    Production-scaling factor (1 at sea level, high latitude).
total_shielding : float, default=1.0
    Combined shielding factor (0-1; 1 = unshielded).
muon_scheme : str, default="BraucherBorchers"
    Muon production parameterisation (``"Braucher"``, ``"BraucherBorchers"``,
    ``"Granger"``, ``"Schaller"`` or ``"newCRONUS"``).
initial_10be, initial_26al, initial_14c : float, default=-1.0
    Uniform initial concentrations (atoms/g); negative starts at steady
    state with the background erosion rate.

Returns
-------
list of CRNColumnSample
    One sample per depth.

```python
def transient(self, sampling_depths_metres: collections.abc.Sequence[float], production_scaling: float = 1.0, total_shielding: float = 1.0, muon_scheme: str = 'BraucherBorchers', initial_10be: float = -1.0, initial_26al: float = -1.0, initial_14c: float = -1.0) -> list[lsdtt3.CRNColumnSample]
```

---

## `CRNColumnSample`

A modeled cosmogenic-nuclide sample within a column.

Concentrations are atoms/g and effective depth is g/cm².

Returned by ``CRNColumn.steady_state()``, ``CRNColumn.step_change()``, ``CRNColumn.transient()``; not constructible directly.

Attributes
----------
concentration_al26 : float
concentration_be10 : float
concentration_c14 : float
depth_m : float
effective_depth_g_cm2 : float

Authors
-------
B. Gailleton

### Members

#### `concentration_al26`

(self) -> float

#### `concentration_be10`

(self) -> float

#### `concentration_c14`

(self) -> float

#### `depth_m`

(self) -> float

#### `effective_depth_g_cm2`

(self) -> float

---

## `CRNStepPairResult`

Result of :func:`crn_step_change_pair`: new erosion rate (g/cm^2/yr), time
since the change (years) and the predicted concentrations (atoms/g);
``status`` is ``"Converged"``, ``"NoSolution"``, ``"NotUnique"`` or
``"NotConverged"``.

Returned by ``crn_step_change_pair()``; not constructible directly.

Attributes
----------
erosion_rate_g_cm2_yr : float
fast_concentration_predicted : float
n_brackets : int
slow_concentration_predicted : float
status : str
time_since_change_yr : float

Authors
-------
B. Gailleton

### Members

#### `erosion_rate_g_cm2_yr`

(self) -> float

#### `fast_concentration_predicted`

(self) -> float

#### `n_brackets`

(self) -> int

#### `slow_concentration_predicted`

(self) -> float

#### `status`

(self) -> str

#### `time_since_change_yr`

(self) -> float

---

## `crn_accumulated_concentration`

Route an upstream weight-averaged CRN concentration to every D8 pixel.

Each pixel gets the mean concentration of its upslope area weighted by
``weight`` (normally the erosion rate, so the result is the
sediment-flux-weighted concentration a river sample would carry).

Parameters
----------
concentration : Raster
    Concentration (atoms/g) on the flow grid.
weight : Raster
    Weights on the same grid, for example erosion rate (g/cm^2/yr).
flow : FlowInfo
    Flow topology on the same grid.

Returns
-------
Raster
    Weighted upslope-mean concentration (atoms/g).

Authors
-------
B. Gailleton

```python
def crn_accumulated_concentration(concentration: lsdtt3.Raster, weight: lsdtt3.Raster, flow: lsdtt3.FlowInfo) -> lsdtt3.Raster
```

---

## `crn_apparent_erosion_rate`

Invert a surface concentration for its apparent steady-state erosion rate.

Parameters
----------
nuclide : {"Be10", "Al26", "C14"}
    Measured nuclide.
concentration_atoms_g : float
    Measured concentration (atoms/g).
total_scaling : float, default=1.0
    Production scaling x shielding factor (dimensionless).
muon_scheme : str, default="BraucherBorchers"
    Muon production parameterisation (``"Braucher"``, ``"BraucherBorchers"``,
    ``"Granger"``, ``"Schaller"`` or ``"newCRONUS"``).

Returns
-------
dict
    ``"erosion_rate_g_cm2_yr"`` (float) and ``"status"`` (``"Converged"``,
    ``"NoSolution"``, ``"NotUnique"`` or ``"NotConverged"``).

Authors
-------
B. Gailleton

```python
def crn_apparent_erosion_rate(nuclide: str, concentration_atoms_g: float, total_scaling: float = 1.0, muon_scheme: str = 'BraucherBorchers') -> dict
```

---

## `crn_burial_pair`

Invert a two-nuclide burial pair. ``nuclide_fast`` must decay faster.

Parameters
----------
nuclide_fast, nuclide_slow : str
    Shorter- and longer-lived nuclides, for example ``"Al26"`` and ``"Be10"``.
concentration_fast_atoms_g, concentration_slow_atoms_g : float
    Measured concentrations (atoms/g).
total_scaling : float, default=1.0
    Production scaling x shielding factor (dimensionless).
muon_scheme : str, default="BraucherBorchers"
    Muon production parameterisation (``"Braucher"``, ``"BraucherBorchers"``,
    ``"Granger"``, ``"Schaller"`` or ``"newCRONUS"``).

Returns
-------
CRNBurialResult
    Pre-burial erosion rate, burial age and ratio diagnostics.

Authors
-------
B. Gailleton

```python
def crn_burial_pair(nuclide_fast: str, nuclide_slow: str, concentration_fast_atoms_g: float, concentration_slow_atoms_g: float, total_scaling: float = 1.0, muon_scheme: str = 'BraucherBorchers') -> lsdtt3.CRNBurialResult
```

---

## `crn_invert_basin_erosion`

Invert basin-mean erosion rates from a point GeoDataFrame.

Points are projected to the flow-grid CRS, snapped to the nearest D8 channel
meeting ``channel_threshold_pixels``, then inverted with the native CAIRN
solver over each snapped outlet's exact upslope pixel set. Failed rows remain
in the returned GeoDataFrame and are identified by ``crn_status``.

Parameters
----------
flow : FlowInfo
    Flow topology of the (conditioned) DEM.
samples : geopandas.GeoDataFrame
    Point geometries of the river samples (any declared CRS; reprojected to
    the flow grid's CRS) with a concentration column.
production_scaling : Raster
    Production-scaling factor on the flow grid (``crn_production_scaling``).
concentration_column : str, default="concentration"
    Column of measured concentrations in atoms/g (positive).
nuclide : {"Be10", "Al26", "C14"}, default="Be10"
    Cosmogenic nuclide.
muon_scheme : str, default="BraucherBorchers"
    Muon production parameterisation: ``"Braucher"``,
    ``"BraucherBorchers"``, ``"Granger"``, ``"Schaller"`` or
    ``"newCRONUS"``.
channel_threshold_pixels : int, default=1000
    Samples snap to the nearest pixel with at least this many contributing
    pixels (positive).
snap_radius_pixels : int, default=25
    Largest snapping search radius (pixels, non-negative).
topographic_shielding : Raster or None, default=None
    Topographic shielding factor (dimensionless, 0-1; 1 = unshielded) on the
    same grid; None means 1 everywhere.
snow_shielding : Raster or None, default=None
    Snow shielding as an effective depth (g/cm^2); None means 0 (no snow).
self_shielding : Raster or None, default=None
    Self (sample-thickness) shielding as an effective depth (g/cm^2); None
    means 0 (surface sample).
rock_density_kg_m3 : float, default=2650.0
    Rock density (kg/m^3) used to convert to mm/kyr.

Returns
-------
geopandas.GeoDataFrame
    Copy of ``samples`` with columns ``crn_status`` (``"ok"``,
    ``"invalid_geometry"``, ``"invalid_concentration"``,
    ``"outside_flow_grid"`` or ``"channel_not_found"``), ``input_x``,
    ``input_y``, ``input_row``, ``input_col``, ``snapped_row``,
    ``snapped_col``, ``snapped_x``, ``snapped_y``, ``outlet_flat``,
    ``channel_contributing_pixels``, ``snap_radius_pixels``,
    ``basin_pixels``, ``basin_area_m2``, ``erosion_rate_g_cm2_yr``,
    ``erosion_rate_mm_kyr``, ``predicted_concentration_atoms_g``.

Authors
-------
B. Gailleton

```python
def crn_invert_basin_erosion(flow: lsdtt3.FlowInfo, samples: object, production_scaling: lsdtt3.Raster, concentration_column: str = 'concentration', nuclide: str = 'Be10', muon_scheme: str = 'BraucherBorchers', channel_threshold_pixels: int = 1000, snap_radius_pixels: int = 25, topographic_shielding: lsdtt3.Raster | None = None, snow_shielding: lsdtt3.Raster | None = None, self_shielding: lsdtt3.Raster | None = None, rock_density_kg_m3: float = 2650.0) -> object
```

---

## `crn_isotope_ratio`

Return an isotope-ratio raster with nodata-accounting counts.

Parameters
----------
numerator, denominator : Raster
    Concentrations (atoms/g) on the same grid, for example Al26 and Be10.

Returns
-------
dict
    ``"raster"``: ratio Raster (nodata where either input is nodata, the
    denominator is 0, or the result is not finite); ``"n_valid"``,
    ``"n_zero_denominator"``, ``"n_one_input_masked"``, ``"n_non_finite"``:
    pixel counts.

Authors
-------
B. Gailleton

```python
def crn_isotope_ratio(numerator: lsdtt3.Raster, denominator: lsdtt3.Raster) -> object
```

---

## `crn_outlet_concentration`

Return erosion-weighted upstream CRN concentration at one outlet pixel.

Parameters
----------
concentration : Raster
    Concentration (atoms/g) on the flow grid.
erosion_rate : Raster
    Erosion rate (g/cm^2/yr) used as the weight.
flow : FlowInfo
    Flow topology on the same grid.
row, col : int
    Outlet pixel.

Returns
-------
float
    Basin-averaged concentration (atoms/g).

Authors
-------
B. Gailleton

```python
def crn_outlet_concentration(concentration: lsdtt3.Raster, erosion_rate: lsdtt3.Raster, flow: lsdtt3.FlowInfo, row: int, col: int) -> float
```

---

## `crn_pair_diagnostics`

Apply native two-nuclide burial or step-change diagnostics to GeoDataFrame rows.

``mode='burial'`` reports pre-burial erosion, burial age and ratio diagnostics.
``mode='step_change'`` reports the new erosion rate, time since change, and
forward-predicted pair concentrations. Invalid rows remain in place and carry
``crn_pair_status='invalid_concentration'``.

Parameters
----------
samples : geopandas.GeoDataFrame
    Rows with two concentration columns (atoms/g).
mode : {"burial", "step_change"}, default="burial"
    Diagnostic to compute.
nuclide_fast, nuclide_slow : str, default="Al26", "Be10"
    Shorter- and longer-lived nuclides of the pair.
concentration_fast_column, concentration_slow_column : str
    Columns holding the two concentrations (atoms/g).
total_scaling : float, default=1.0
    Combined production scaling x shielding factor (dimensionless,
    positive).
muon_scheme : str, default="BraucherBorchers"
    Muon production parameterisation: ``"Braucher"``,
    ``"BraucherBorchers"``, ``"Granger"``, ``"Schaller"`` or
    ``"newCRONUS"``.
effective_erosion_rate_old_g_cm2_yr : float, default=0.01
    ``"step_change"`` only: erosion rate before the change (g/cm^2/yr).

Returns
-------
geopandas.GeoDataFrame
    Copy of ``samples`` with ``crn_pair_status``,
    ``crn_fast_concentration_atoms_g``, ``crn_slow_concentration_atoms_g``,
    ``crn_pair_erosion_rate_g_cm2_yr`` (pre-burial or new rate) and
    ``crn_pair_time_yr`` (burial age or time since change), plus
    ``crn_pair_measured_ratio``, ``crn_pair_steady_state_ratio``,
    ``crn_pair_burial_required`` (burial) or
    ``crn_pair_fast_predicted_atoms_g``, ``crn_pair_slow_predicted_atoms_g``,
    ``crn_pair_n_brackets`` (step change).

Authors
-------
B. Gailleton

```python
def crn_pair_diagnostics(samples: object, mode: str = 'burial', nuclide_fast: str = 'Al26', nuclide_slow: str = 'Be10', concentration_fast_column: str = 'Al26_CONC', concentration_slow_column: str = 'Be10_CONC', total_scaling: float = 1.0, muon_scheme: str = 'BraucherBorchers', effective_erosion_rate_old_g_cm2_yr: float = 0.01) -> object
```

---

## `crn_pressure`

Return atmospheric pressure (hPa) for a DEM.

Parameters
----------
dem : Raster
    Elevation (m) in a projected metre CRS.
crn_parameters_path : str, default=""
    CRN parameter NetCDF file, as in ``crn_production_scaling``.
atmosphere : {"auto", "standard", "era40", "ncep"}, default="auto"
    Atmospheric model, as in ``crn_production_scaling``.

Returns
-------
Raster
    Pressure in hPa (mbar).

Authors
-------
B. Gailleton

```python
def crn_pressure(dem: lsdtt3.Raster, crn_parameters_path: str = '', atmosphere: str = 'auto') -> lsdtt3.Raster
```

---

## `crn_production_scaling`

Return a per-pixel CRN production-scaling raster. DEM coordinates and elevations must be metres.

The factor multiplies sea-level, high-latitude production; it depends on each
pixel's elevation and latitude (from the DEM's projected CRS).

Parameters
----------
dem : Raster
    Elevation (m) in a projected metre CRS.
crn_parameters_path : str, default=""
    Path to the lsdtt3 CRN parameter NetCDF file (atmosphere reanalysis,
    geomagnetic and LSDn grids). Empty uses the built-in standard
    atmosphere and modern field only.
sample_age_yr : float, default=0.0
    Exposure age (years) for the time-integrated palaeomagnetic correction
    (needs the parameter file); 0 uses the modern field.
atmosphere : {"auto", "standard", "era40", "ncep"}, default="auto"
    Atmospheric pressure model: ``"auto"`` uses ERA40 when the parameter file
    is loaded, else the standard atmosphere.
scaling : {"stone", "lsdn"}, default="stone"
    Stone (2000) scaling, or Lifton, Sato and Dunai (2014) LSDn scaling
    (needs the LSDn grids in the parameter file; Be10, Al26, C14 only).
nuclide : {"Be10", "Al26", "C14"}, default="Be10"
    Nuclide (affects LSDn scaling only).

Returns
-------
Raster
    Dimensionless production-scaling factor per pixel.

Authors
-------
B. Gailleton

```python
def crn_production_scaling(dem: lsdtt3.Raster, crn_parameters_path: str = '', sample_age_yr: float = 0.0, atmosphere: str = 'auto', scaling: str = 'stone', nuclide: str = 'Be10') -> lsdtt3.Raster
```

---

## `crn_steady_state`

Return whole-raster steady-state CRN concentration (atoms/g).

The surface concentration each pixel would have after eroding at a constant
rate for long enough to reach steady state.

Parameters
----------
erosion_rate : Raster
    Effective erosion rate in g/cm^2/yr (mass flux; divide by rock density
    in g/cm^3 for cm/yr, e.g. 0.01 g/cm^2/yr at 2.65 g/cm^3 is 37.7 mm/kyr).
production_scaling : Raster
    Production-scaling factor on the same grid (``crn_production_scaling``).
nuclide : {"Be10", "Al26", "C14"}, default="Be10"
    Cosmogenic nuclide.
muon_scheme : str, default="BraucherBorchers"
    Muon production parameterisation: ``"Braucher"``,
    ``"BraucherBorchers"``, ``"Granger"``, ``"Schaller"`` or
    ``"newCRONUS"``.
topographic_shielding : Raster or None, default=None
    Topographic shielding factor (dimensionless, 0-1; 1 = unshielded) on the
    same grid; None means 1 everywhere.
snow_shielding : Raster or None, default=None
    Snow shielding as an effective depth (g/cm^2); None means 0 (no snow).
self_shielding : Raster or None, default=None
    Self (sample-thickness) shielding as an effective depth (g/cm^2); None
    means 0 (surface sample).

Returns
-------
Raster
    Concentration in atoms/g.

Authors
-------
B. Gailleton

```python
def crn_steady_state(erosion_rate: lsdtt3.Raster, production_scaling: lsdtt3.Raster, nuclide: str = 'Be10', muon_scheme: str = 'BraucherBorchers', topographic_shielding: lsdtt3.Raster | None = None, snow_shielding: lsdtt3.Raster | None = None, self_shielding: lsdtt3.Raster | None = None) -> lsdtt3.Raster
```

---

## `crn_step_change`

Return a whole-raster surface concentration after a step change in erosion.

Each pixel was at steady state with ``erosion_rate_old_g_cm2_yr`` and has
eroded at ``erosion_rate_new_g_cm2_yr`` for ``time_since_change`` years.

Parameters
----------
time_since_change : Raster
    Years since the change reached each pixel; <= 0 where it has not
    arrived.
production_scaling : Raster
    Production-scaling factor on the same grid.
nuclide : {"Be10", "Al26", "C14"}, default="Be10"
    Cosmogenic nuclide.
muon_scheme : str, default="BraucherBorchers"
    Muon production parameterisation: ``"Braucher"``,
    ``"BraucherBorchers"``, ``"Granger"``, ``"Schaller"`` or
    ``"newCRONUS"``.
topographic_shielding : Raster or None, default=None
    Topographic shielding factor (0-1); None means 1.
erosion_rate_old_g_cm2_yr : float, default=0.01
    Erosion rate before the change (g/cm^2/yr).
erosion_rate_new_g_cm2_yr : float, default=0.01
    Erosion rate after the change (g/cm^2/yr).
pre_change_erosion_rate : Raster or None, default=None
    Optional per-pixel pre-change rate (g/cm^2/yr) used where
    ``time_since_change <= 0``; None uses ``erosion_rate_old_g_cm2_yr``.

Returns
-------
Raster
    Surface concentration in atoms/g.

Authors
-------
B. Gailleton

```python
def crn_step_change(time_since_change: lsdtt3.Raster, production_scaling: lsdtt3.Raster, nuclide: str = 'Be10', muon_scheme: str = 'BraucherBorchers', topographic_shielding: lsdtt3.Raster | None = None, erosion_rate_old_g_cm2_yr: float = 0.01, erosion_rate_new_g_cm2_yr: float = 0.01, pre_change_erosion_rate: lsdtt3.Raster | None = None) -> lsdtt3.Raster
```

---

## `crn_step_change_pair`

Invert a two-nuclide step-change pair for new erosion rate and timing.

Parameters
----------
nuclide_fast, nuclide_slow : str
    Shorter- and longer-lived nuclides.
concentration_fast_atoms_g, concentration_slow_atoms_g : float
    Measured concentrations (atoms/g).
effective_erosion_rate_old_g_cm2_yr : float
    Erosion rate before the change (g/cm^2/yr).
total_scaling : float, default=1.0
    Production scaling x shielding factor (dimensionless).
muon_scheme : str, default="BraucherBorchers"
    Muon production parameterisation (``"Braucher"``, ``"BraucherBorchers"``,
    ``"Granger"``, ``"Schaller"`` or ``"newCRONUS"``).

Returns
-------
CRNStepPairResult
    New erosion rate, time since the change and predicted concentrations.

Authors
-------
B. Gailleton

```python
def crn_step_change_pair(nuclide_fast: str, nuclide_slow: str, concentration_fast_atoms_g: float, concentration_slow_atoms_g: float, effective_erosion_rate_old_g_cm2_yr: float, total_scaling: float = 1.0, muon_scheme: str = 'BraucherBorchers') -> lsdtt3.CRNStepPairResult
```

---

## `crn_transient_history`

Predict Be-10, Al-26, and C-14 raster concentrations for an arbitrary erosion history.

Parameters
----------
column : CRNColumn
    Column with a validated present-to-past erosion/removal history.
production_scaling : Raster
    Per-pixel production-scaling factor.
topographic_shielding : Raster, optional
    Same-grid shielding factor. Unity is used when omitted.
sampling_depth_metres : float, optional
    Depth (m) below the present surface at which to predict concentration.
muon_scheme : str, default="BraucherBorchers"
    Muon production parameterisation (``"Braucher"``, ``"BraucherBorchers"``,
    ``"Granger"``, ``"Schaller"`` or ``"newCRONUS"``).

Returns
-------
dict[str, Raster]
    ``be10``, ``al26``, and ``c14`` concentration rasters in atoms g⁻¹.

Notes
-----
The calculation is O(number of valid pixels × history rows).

Authors
-------
B. Gailleton
S. M. Mudd

**Deprecated compatibility names:**
- `sampling_depth_m` -> `sampling_depth_metres`; earliest removal: following minor.

```python
def crn_transient_history(column: lsdtt3.CRNColumn, production_scaling: lsdtt3.Raster, topographic_shielding: lsdtt3.Raster | None = None, sampling_depth_metres: float = 0.0, muon_scheme: str = 'BraucherBorchers') -> object
```

---

## `crn_uniform_raster`

Create a masked constant raster for CRN erosion or shielding inputs.

Parameters
----------
reference : Raster
    Raster whose grid and nodata mask are copied.
value : float
    Constant value, for example an erosion rate in g/cm^2/yr.

Returns
-------
Raster
    ``value`` on every valid pixel of ``reference``.

Authors
-------
B. Gailleton

```python
def crn_uniform_raster(reference: lsdtt3.Raster, value: float) -> lsdtt3.Raster
```
