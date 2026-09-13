# Cosmogenic Nuclides API

This page is generated from the built Python runtime manifest. Do not edit signatures here.

## `CRNBurialResult`

No additional binding notes are available.

Parameters
----------
*args
    Constructor arguments described by the generated binding signature.
**kwargs
    Constructor keyword arguments described by the generated binding
    signature.

Attributes
----------
See the generated binding signature for available properties and methods.

Notes
-----
No additional binding notes are available.

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

    Parameters
    ----------
    *args
        Constructor arguments described by the generated binding signature.
    **kwargs
        Constructor keyword arguments described by the generated binding
        signature.

    Attributes
    ----------
    See the generated binding signature for available properties and methods.

    Notes
    -----
    One-dimensional cosmogenic-nuclide column model.

Set an erosion history, ordered from present backwards in time, then use
``transient`` to model depth profiles. Point steady-state and step-change
calculations require no history. Erosion rates for these point calculations
are effective rates in g/cm²/yr.

    Authors
    -------
    B. Gailleton


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

(self) -> float

#### `erosion_rate_m_yr`

(self) -> list[float]

#### `read_erosion_history`

read_erosion_history(self, path: str) -> None

```python
def read_erosion_history(self, path: str) -> None
```

#### `removal_m`

(self) -> list[float]

#### `rock_density_kg_m3`

(self) -> float

#### `set_erosion_history`

set_erosion_history(self, time_yr_bp: collections.abc.Sequence[float], erosion_rate_m_yr: collections.abc.Sequence[float], removal_m: collections.abc.Sequence[float]) -> None

Set the history. Last row must be ``time=-9999`` with a positive background rate.

```python
def set_erosion_history(self, time_yr_bp: collections.abc.Sequence[float], erosion_rate_m_yr: collections.abc.Sequence[float], removal_m: collections.abc.Sequence[float]) -> None
```

#### `steady_state`

steady_state(self, effective_depth_g_cm2: float = 0.0, effective_erosion_rate_g_cm2_yr: float = 0.01, production_scaling: float = 1.0, total_shielding: float = 1.0, muon_scheme: str = 'BraucherBorchers') -> lsdtt3.CRNColumnSample

```python
def steady_state(self, effective_depth_g_cm2: float = 0.0, effective_erosion_rate_g_cm2_yr: float = 0.01, production_scaling: float = 1.0, total_shielding: float = 1.0, muon_scheme: str = 'BraucherBorchers') -> lsdtt3.CRNColumnSample
```

#### `step_change`

step_change(self, effective_depth_g_cm2: float = 0.0, effective_erosion_rate_old_g_cm2_yr: float = 0.01, effective_erosion_rate_new_g_cm2_yr: float = 0.01, time_since_step_yr: float = 0.0, production_scaling: float = 1.0, total_shielding: float = 1.0, muon_scheme: str = 'BraucherBorchers') -> lsdtt3.CRNColumnSample

```python
def step_change(self, effective_depth_g_cm2: float = 0.0, effective_erosion_rate_old_g_cm2_yr: float = 0.01, effective_erosion_rate_new_g_cm2_yr: float = 0.01, time_since_step_yr: float = 0.0, production_scaling: float = 1.0, total_shielding: float = 1.0, muon_scheme: str = 'BraucherBorchers') -> lsdtt3.CRNColumnSample
```

#### `time_yr_bp`

(self) -> list[float]

#### `transient`

transient(self, sampling_depths_m: collections.abc.Sequence[float], production_scaling: float = 1.0, total_shielding: float = 1.0, muon_scheme: str = 'BraucherBorchers', initial_10be: float = -1.0, initial_26al: float = -1.0, initial_14c: float = -1.0) -> list[lsdtt3.CRNColumnSample]

```python
def transient(self, sampling_depths_m: collections.abc.Sequence[float], production_scaling: float = 1.0, total_shielding: float = 1.0, muon_scheme: str = 'BraucherBorchers', initial_10be: float = -1.0, initial_26al: float = -1.0, initial_14c: float = -1.0) -> list[lsdtt3.CRNColumnSample]
```

---

## `CRNColumnSample`

A modeled cosmogenic-nuclide sample within a column.

    Parameters
    ----------
    *args
        Constructor arguments described by the generated binding signature.
    **kwargs
        Constructor keyword arguments described by the generated binding
        signature.

    Attributes
    ----------
    See the generated binding signature for available properties and methods.

    Notes
    -----
    A modeled cosmogenic-nuclide sample within a column.

Concentrations are atoms/g and effective depth is g/cm².

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

No additional binding notes are available.

Parameters
----------
*args
    Constructor arguments described by the generated binding signature.
**kwargs
    Constructor keyword arguments described by the generated binding
    signature.

Attributes
----------
See the generated binding signature for available properties and methods.

Notes
-----
No additional binding notes are available.

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

crn_accumulated_concentration(concentration: lsdtt3.Raster, weight: lsdtt3.Raster, flow: lsdtt3.FlowInfo) -> lsdtt3.Raster

    Parameters
    ----------
    *args
        Positional arguments described by the generated binding signature.
    **kwargs
        Keyword arguments described by the generated binding signature.

    Returns
    -------
    object
        Result described by the generated binding signature.

    Notes
    -----
    crn_accumulated_concentration(concentration: lsdtt3.Raster, weight: lsdtt3.Raster, flow: lsdtt3.FlowInfo) -> lsdtt3.Raster

Route an upstream weight-averaged CRN concentration to every D8 pixel.

    Authors
    -------
    B. Gailleton


```python
def crn_accumulated_concentration(concentration: lsdtt3.Raster, weight: lsdtt3.Raster, flow: lsdtt3.FlowInfo) -> lsdtt3.Raster
```

---

## `crn_apparent_erosion_rate`

crn_apparent_erosion_rate(nuclide: str, concentration_atoms_g: float, total_scaling: float = 1.0, muon_scheme: str = 'BraucherBorchers') -> dict

    Parameters
    ----------
    *args
        Positional arguments described by the generated binding signature.
    **kwargs
        Keyword arguments described by the generated binding signature.

    Returns
    -------
    object
        Result described by the generated binding signature.

    Notes
    -----
    crn_apparent_erosion_rate(nuclide: str, concentration_atoms_g: float, total_scaling: float = 1.0, muon_scheme: str = 'BraucherBorchers') -> dict

Invert a surface concentration for its apparent steady-state erosion rate.

    Authors
    -------
    B. Gailleton


```python
def crn_apparent_erosion_rate(nuclide: str, concentration_atoms_g: float, total_scaling: float = 1.0, muon_scheme: str = 'BraucherBorchers') -> dict
```

---

## `crn_burial_pair`

crn_burial_pair(nuclide_fast: str, nuclide_slow: str, concentration_fast_atoms_g: float, concentration_slow_atoms_g: float, total_scaling: float = 1.0, muon_scheme: str = 'BraucherBorchers') -> lsdtt3.CRNBurialResult

    Parameters
    ----------
    *args
        Positional arguments described by the generated binding signature.
    **kwargs
        Keyword arguments described by the generated binding signature.

    Returns
    -------
    object
        Result described by the generated binding signature.

    Notes
    -----
    crn_burial_pair(nuclide_fast: str, nuclide_slow: str, concentration_fast_atoms_g: float, concentration_slow_atoms_g: float, total_scaling: float = 1.0, muon_scheme: str = 'BraucherBorchers') -> lsdtt3.CRNBurialResult

Invert a two-nuclide burial pair. ``nuclide_fast`` must decay faster.

    Authors
    -------
    B. Gailleton


```python
def crn_burial_pair(nuclide_fast: str, nuclide_slow: str, concentration_fast_atoms_g: float, concentration_slow_atoms_g: float, total_scaling: float = 1.0, muon_scheme: str = 'BraucherBorchers') -> lsdtt3.CRNBurialResult
```

---

## `crn_invert_basin_erosion`

crn_invert_basin_erosion(flow: lsdtt3.FlowInfo, samples: object, production_scaling: lsdtt3.Raster, concentration_column: str = 'concentration', nuclide: str = 'Be10', muon_scheme: str = 'BraucherBorchers', channel_threshold_pixels: int = 1000, snap_radius_pixels: int = 25, topographic_shielding: lsdtt3.Raster | None = None, snow_shielding: lsdtt3.Raster | None = None, self_shielding: lsdtt3.Raster | None = None, rock_density_kg_m3: float = 2650.0) -> object

    Parameters
    ----------
    *args
        Positional arguments described by the generated binding signature.
    **kwargs
        Keyword arguments described by the generated binding signature.

    Returns
    -------
    object
        Result described by the generated binding signature.

    Notes
    -----
    crn_invert_basin_erosion(flow: lsdtt3.FlowInfo, samples: object, production_scaling: lsdtt3.Raster, concentration_column: str = 'concentration', nuclide: str = 'Be10', muon_scheme: str = 'BraucherBorchers', channel_threshold_pixels: int = 1000, snap_radius_pixels: int = 25, topographic_shielding: lsdtt3.Raster | None = None, snow_shielding: lsdtt3.Raster | None = None, self_shielding: lsdtt3.Raster | None = None, rock_density_kg_m3: float = 2650.0) -> object

Invert basin-mean erosion rates from a point GeoDataFrame.

Points are projected to the flow-grid CRS, snapped to the nearest D8 channel
meeting ``channel_threshold_pixels``, then inverted with the native CAIRN
solver over each snapped outlet's exact upslope pixel set. Failed rows remain
in the returned GeoDataFrame and are identified by ``crn_status``.

    Authors
    -------
    B. Gailleton


```python
def crn_invert_basin_erosion(flow: lsdtt3.FlowInfo, samples: object, production_scaling: lsdtt3.Raster, concentration_column: str = 'concentration', nuclide: str = 'Be10', muon_scheme: str = 'BraucherBorchers', channel_threshold_pixels: int = 1000, snap_radius_pixels: int = 25, topographic_shielding: lsdtt3.Raster | None = None, snow_shielding: lsdtt3.Raster | None = None, self_shielding: lsdtt3.Raster | None = None, rock_density_kg_m3: float = 2650.0) -> object
```

---

## `crn_isotope_ratio`

crn_isotope_ratio(numerator: lsdtt3.Raster, denominator: lsdtt3.Raster) -> object

    Parameters
    ----------
    *args
        Positional arguments described by the generated binding signature.
    **kwargs
        Keyword arguments described by the generated binding signature.

    Returns
    -------
    object
        Result described by the generated binding signature.

    Notes
    -----
    crn_isotope_ratio(numerator: lsdtt3.Raster, denominator: lsdtt3.Raster) -> object

Return an isotope-ratio raster with nodata-accounting counts.

    Authors
    -------
    B. Gailleton


```python
def crn_isotope_ratio(numerator: lsdtt3.Raster, denominator: lsdtt3.Raster) -> object
```

---

## `crn_outlet_concentration`

crn_outlet_concentration(concentration: lsdtt3.Raster, erosion_rate: lsdtt3.Raster, flow: lsdtt3.FlowInfo, row: int, col: int) -> float

    Parameters
    ----------
    *args
        Positional arguments described by the generated binding signature.
    **kwargs
        Keyword arguments described by the generated binding signature.

    Returns
    -------
    object
        Result described by the generated binding signature.

    Notes
    -----
    crn_outlet_concentration(concentration: lsdtt3.Raster, erosion_rate: lsdtt3.Raster, flow: lsdtt3.FlowInfo, row: int, col: int) -> float

Return erosion-weighted upstream CRN concentration at one outlet pixel.

    Authors
    -------
    B. Gailleton


```python
def crn_outlet_concentration(concentration: lsdtt3.Raster, erosion_rate: lsdtt3.Raster, flow: lsdtt3.FlowInfo, row: int, col: int) -> float
```

---

## `crn_pair_diagnostics`

crn_pair_diagnostics(samples: object, mode: str = 'burial', nuclide_fast: str = 'Al26', nuclide_slow: str = 'Be10', concentration_fast_column: str = 'Al26_CONC', concentration_slow_column: str = 'Be10_CONC', total_scaling: float = 1.0, muon_scheme: str = 'BraucherBorchers', effective_erosion_rate_old_g_cm2_yr: float = 0.01) -> object

    Parameters
    ----------
    *args
        Positional arguments described by the generated binding signature.
    **kwargs
        Keyword arguments described by the generated binding signature.

    Returns
    -------
    object
        Result described by the generated binding signature.

    Notes
    -----
    crn_pair_diagnostics(samples: object, mode: str = 'burial', nuclide_fast: str = 'Al26', nuclide_slow: str = 'Be10', concentration_fast_column: str = 'Al26_CONC', concentration_slow_column: str = 'Be10_CONC', total_scaling: float = 1.0, muon_scheme: str = 'BraucherBorchers', effective_erosion_rate_old_g_cm2_yr: float = 0.01) -> object

Apply native two-nuclide burial or step-change diagnostics to GeoDataFrame rows.

``mode='burial'`` reports pre-burial erosion, burial age and ratio diagnostics.
``mode='step_change'`` reports the new erosion rate, time since change, and
forward-predicted pair concentrations. Invalid rows remain in place and carry
``crn_pair_status='invalid_concentration'``.

    Authors
    -------
    B. Gailleton


```python
def crn_pair_diagnostics(samples: object, mode: str = 'burial', nuclide_fast: str = 'Al26', nuclide_slow: str = 'Be10', concentration_fast_column: str = 'Al26_CONC', concentration_slow_column: str = 'Be10_CONC', total_scaling: float = 1.0, muon_scheme: str = 'BraucherBorchers', effective_erosion_rate_old_g_cm2_yr: float = 0.01) -> object
```

---

## `crn_pressure`

crn_pressure(dem: lsdtt3.Raster, crn_parameters_path: str = '', atmosphere: str = 'auto') -> lsdtt3.Raster

    Parameters
    ----------
    *args
        Positional arguments described by the generated binding signature.
    **kwargs
        Keyword arguments described by the generated binding signature.

    Returns
    -------
    object
        Result described by the generated binding signature.

    Notes
    -----
    crn_pressure(dem: lsdtt3.Raster, crn_parameters_path: str = '', atmosphere: str = 'auto') -> lsdtt3.Raster

Return atmospheric pressure (hPa) for a DEM.

    Authors
    -------
    B. Gailleton


```python
def crn_pressure(dem: lsdtt3.Raster, crn_parameters_path: str = '', atmosphere: str = 'auto') -> lsdtt3.Raster
```

---

## `crn_production_scaling`

crn_production_scaling(dem: lsdtt3.Raster, crn_parameters_path: str = '', sample_age_yr: float = 0.0, atmosphere: str = 'auto', scaling: str = 'stone', nuclide: str = 'Be10') -> lsdtt3.Raster

    Parameters
    ----------
    *args
        Positional arguments described by the generated binding signature.
    **kwargs
        Keyword arguments described by the generated binding signature.

    Returns
    -------
    object
        Result described by the generated binding signature.

    Notes
    -----
    crn_production_scaling(dem: lsdtt3.Raster, crn_parameters_path: str = '', sample_age_yr: float = 0.0, atmosphere: str = 'auto', scaling: str = 'stone', nuclide: str = 'Be10') -> lsdtt3.Raster

Return a per-pixel CRN production-scaling raster. DEM coordinates and elevations must be metres.

    Authors
    -------
    B. Gailleton


```python
def crn_production_scaling(dem: lsdtt3.Raster, crn_parameters_path: str = '', sample_age_yr: float = 0.0, atmosphere: str = 'auto', scaling: str = 'stone', nuclide: str = 'Be10') -> lsdtt3.Raster
```

---

## `crn_steady_state`

crn_steady_state(erosion_rate: lsdtt3.Raster, production_scaling: lsdtt3.Raster, nuclide: str = 'Be10', muon_scheme: str = 'BraucherBorchers', topographic_shielding: lsdtt3.Raster | None = None, snow_shielding: lsdtt3.Raster | None = None, self_shielding: lsdtt3.Raster | None = None) -> lsdtt3.Raster

    Parameters
    ----------
    *args
        Positional arguments described by the generated binding signature.
    **kwargs
        Keyword arguments described by the generated binding signature.

    Returns
    -------
    object
        Result described by the generated binding signature.

    Notes
    -----
    crn_steady_state(erosion_rate: lsdtt3.Raster, production_scaling: lsdtt3.Raster, nuclide: str = 'Be10', muon_scheme: str = 'BraucherBorchers', topographic_shielding: lsdtt3.Raster | None = None, snow_shielding: lsdtt3.Raster | None = None, self_shielding: lsdtt3.Raster | None = None) -> lsdtt3.Raster

Return whole-raster steady-state CRN concentration (atoms/g).

    Authors
    -------
    B. Gailleton


```python
def crn_steady_state(erosion_rate: lsdtt3.Raster, production_scaling: lsdtt3.Raster, nuclide: str = 'Be10', muon_scheme: str = 'BraucherBorchers', topographic_shielding: lsdtt3.Raster | None = None, snow_shielding: lsdtt3.Raster | None = None, self_shielding: lsdtt3.Raster | None = None) -> lsdtt3.Raster
```

---

## `crn_step_change`

crn_step_change(time_since_change: lsdtt3.Raster, production_scaling: lsdtt3.Raster, nuclide: str = 'Be10', muon_scheme: str = 'BraucherBorchers', topographic_shielding: lsdtt3.Raster | None = None, erosion_rate_old_g_cm2_yr: float = 0.01, erosion_rate_new_g_cm2_yr: float = 0.01, pre_change_erosion_rate: lsdtt3.Raster | None = None) -> lsdtt3.Raster

    Parameters
    ----------
    *args
        Positional arguments described by the generated binding signature.
    **kwargs
        Keyword arguments described by the generated binding signature.

    Returns
    -------
    object
        Result described by the generated binding signature.

    Notes
    -----
    crn_step_change(time_since_change: lsdtt3.Raster, production_scaling: lsdtt3.Raster, nuclide: str = 'Be10', muon_scheme: str = 'BraucherBorchers', topographic_shielding: lsdtt3.Raster | None = None, erosion_rate_old_g_cm2_yr: float = 0.01, erosion_rate_new_g_cm2_yr: float = 0.01, pre_change_erosion_rate: lsdtt3.Raster | None = None) -> lsdtt3.Raster

Return a whole-raster surface concentration after a step change in erosion.

    Authors
    -------
    B. Gailleton


```python
def crn_step_change(time_since_change: lsdtt3.Raster, production_scaling: lsdtt3.Raster, nuclide: str = 'Be10', muon_scheme: str = 'BraucherBorchers', topographic_shielding: lsdtt3.Raster | None = None, erosion_rate_old_g_cm2_yr: float = 0.01, erosion_rate_new_g_cm2_yr: float = 0.01, pre_change_erosion_rate: lsdtt3.Raster | None = None) -> lsdtt3.Raster
```

---

## `crn_step_change_pair`

crn_step_change_pair(nuclide_fast: str, nuclide_slow: str, concentration_fast_atoms_g: float, concentration_slow_atoms_g: float, effective_erosion_rate_old_g_cm2_yr: float, total_scaling: float = 1.0, muon_scheme: str = 'BraucherBorchers') -> lsdtt3.CRNStepPairResult

    Parameters
    ----------
    *args
        Positional arguments described by the generated binding signature.
    **kwargs
        Keyword arguments described by the generated binding signature.

    Returns
    -------
    object
        Result described by the generated binding signature.

    Notes
    -----
    crn_step_change_pair(nuclide_fast: str, nuclide_slow: str, concentration_fast_atoms_g: float, concentration_slow_atoms_g: float, effective_erosion_rate_old_g_cm2_yr: float, total_scaling: float = 1.0, muon_scheme: str = 'BraucherBorchers') -> lsdtt3.CRNStepPairResult

Invert a two-nuclide step-change pair for new erosion rate and timing.

    Authors
    -------
    B. Gailleton


```python
def crn_step_change_pair(nuclide_fast: str, nuclide_slow: str, concentration_fast_atoms_g: float, concentration_slow_atoms_g: float, effective_erosion_rate_old_g_cm2_yr: float, total_scaling: float = 1.0, muon_scheme: str = 'BraucherBorchers') -> lsdtt3.CRNStepPairResult
```

---

## `crn_transient_history`

crn_transient_history(column: lsdtt3.CRNColumn, production_scaling: lsdtt3.Raster, topographic_shielding: lsdtt3.Raster | None = None, sampling_depth_m: float = 0.0, muon_scheme: str = 'BraucherBorchers') -> object

    Parameters
    ----------
    *args
        Positional arguments described by the generated binding signature.
    **kwargs
        Keyword arguments described by the generated binding signature.

    Returns
    -------
    object
        Result described by the generated binding signature.

    Notes
    -----
    crn_transient_history(column: lsdtt3.CRNColumn, production_scaling: lsdtt3.Raster, topographic_shielding: lsdtt3.Raster | None = None, sampling_depth_m: float = 0.0, muon_scheme: str = 'BraucherBorchers') -> object

Predict Be-10, Al-26, and C-14 raster concentrations for an arbitrary erosion history.

Parameters
----------
column : CRNColumn
    Column with a validated present-to-past erosion/removal history.
production_scaling : Raster
    Per-pixel production-scaling factor.
topographic_shielding : Raster, optional
    Same-grid shielding factor. Unity is used when omitted.
sampling_depth_m : float, optional
    Depth below the present surface at which to predict concentration.

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

    Authors
    -------
    B. Gailleton


```python
def crn_transient_history(column: lsdtt3.CRNColumn, production_scaling: lsdtt3.Raster, topographic_shielding: lsdtt3.Raster | None = None, sampling_depth_m: float = 0.0, muon_scheme: str = 'BraucherBorchers') -> object
```

---

## `crn_uniform_raster`

crn_uniform_raster(reference: lsdtt3.Raster, value: float) -> lsdtt3.Raster

    Parameters
    ----------
    *args
        Positional arguments described by the generated binding signature.
    **kwargs
        Keyword arguments described by the generated binding signature.

    Returns
    -------
    object
        Result described by the generated binding signature.

    Notes
    -----
    crn_uniform_raster(reference: lsdtt3.Raster, value: float) -> lsdtt3.Raster

Create a masked constant raster for CRN erosion or shielding inputs.

    Authors
    -------
    B. Gailleton


```python
def crn_uniform_raster(reference: lsdtt3.Raster, value: float) -> lsdtt3.Raster
```
