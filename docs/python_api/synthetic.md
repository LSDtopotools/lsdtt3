# Synthetic Workflows API

This page is generated from the built Python runtime manifest. Do not edit signatures here.

## `SyntheticNetworkResult`

In-memory optimal-channel-network outputs.

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
In-memory optimal-channel-network outputs.

Authors
-------
B. Gailleton

### Members

#### `drainage_area_cells`

(self) -> lsdtt3.Raster

#### `elevation`

(self) -> lsdtt3.Raster

#### `final_energy`

(self) -> float

#### `initial_energy`

(self) -> float

#### `receivers`

(self) -> object

---

## `synthetic_hillslope_strip`

synthetic_hillslope_strip(metadata: lsdtt3.RasterMeta, channel_config: str = 'both_ends', init_profile: str = 'flat', cross_axis_is_x: bool = True, baselevel_elevation: float = 0.0, half_width: float = 0.0, uplift_rate: float = 0.001, diffusivity: float = 0.01, critical_slope: float = 0.8, nonlinear: bool = False, initial_amplitude: float = 50.0, seed: int = 1, n_threads: int = 0) -> dict

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
    synthetic_hillslope_strip(metadata: lsdtt3.RasterMeta, channel_config: str = 'both_ends', init_profile: str = 'flat', cross_axis_is_x: bool = True, baselevel_elevation: float = 0.0, half_width: float = 0.0, uplift_rate: float = 0.001, diffusivity: float = 0.01, critical_slope: float = 0.8, nonlinear: bool = False, initial_amplitude: float = 50.0, seed: int = 1, n_threads: int = 0) -> dict

Generate a one-dimensional hillslope test DEM and fixed channel mask.

    Authors
    -------
    B. Gailleton


```python
def synthetic_hillslope_strip(metadata: lsdtt3.RasterMeta, channel_config: str = 'both_ends', init_profile: str = 'flat', cross_axis_is_x: bool = True, baselevel_elevation: float = 0.0, half_width: float = 0.0, uplift_rate: float = 0.001, diffusivity: float = 0.01, critical_slope: float = 0.8, nonlinear: bool = False, initial_amplitude: float = 50.0, seed: int = 1, n_threads: int = 0) -> dict
```

---

## `synthetic_mountain_range`

synthetic_mountain_range(metadata: lsdtt3.RasterMeta, geometry: str = 'linear', edge: str = 'north', ridge_elevation: float = 2000.0, piedmont_elevation: float = 0.0, front_fraction: float = 0.2, noise_amplitude: float = 100.0, seed: int = 1, n_threads: int = 0) -> dict

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
    synthetic_mountain_range(metadata: lsdtt3.RasterMeta, geometry: str = 'linear', edge: str = 'north', ridge_elevation: float = 2000.0, piedmont_elevation: float = 0.0, front_fraction: float = 0.2, noise_amplitude: float = 100.0, seed: int = 1, n_threads: int = 0) -> dict

Generate a synthetic mountain range and its known mountain mask.

    Authors
    -------
    B. Gailleton


```python
def synthetic_mountain_range(metadata: lsdtt3.RasterMeta, geometry: str = 'linear', edge: str = 'north', ridge_elevation: float = 2000.0, piedmont_elevation: float = 0.0, front_fraction: float = 0.2, noise_amplitude: float = 100.0, seed: int = 1, n_threads: int = 0) -> dict
```

---

## `synthetic_network`

synthetic_network(metadata: lsdtt3.RasterMeta, outlet_mode: str = 'single_south', periodic_columns: bool = False, outlet_spacing: int = 1, gamma: float = 0.5, n_iterations: int = -1, temperature_scale: float = 1.0, m_exponent: float = 0.45, n_exponent: float = 1.0, relief_metres: float = 1000.0, seed: int = 1) -> lsdtt3.SyntheticNetworkResult

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
    synthetic_network(metadata: lsdtt3.RasterMeta, outlet_mode: str = 'single_south', periodic_columns: bool = False, outlet_spacing: int = 1, gamma: float = 0.5, n_iterations: int = -1, temperature_scale: float = 1.0, m_exponent: float = 0.45, n_exponent: float = 1.0, relief_metres: float = 1000.0, seed: int = 1) -> lsdtt3.SyntheticNetworkResult

Optimise a synthetic D8 optimal channel network and derive its chi steady-state DEM.

    Authors
    -------
    B. Gailleton


**Deprecated compatibility names:**
- `relief_m` -> `relief_metres`; earliest removal: following minor.

```python
def synthetic_network(metadata: lsdtt3.RasterMeta, outlet_mode: str = 'single_south', periodic_columns: bool = False, outlet_spacing: int = 1, gamma: float = 0.5, n_iterations: int = -1, temperature_scale: float = 1.0, m_exponent: float = 0.45, n_exponent: float = 1.0, relief_metres: float = 1000.0, seed: int = 1) -> lsdtt3.SyntheticNetworkResult
```

---

## `synthetic_noise`

synthetic_noise(metadata: lsdtt3.RasterMeta, kind: str = 'red', seed: int = 1, mean: float = 0.0, stddev: float = 1.0, blur_passes: int = 6, roughness: float = 0.55, beta: float = 2.0, n_threads: int = 0) -> lsdtt3.Raster

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
    synthetic_noise(metadata: lsdtt3.RasterMeta, kind: str = 'red', seed: int = 1, mean: float = 0.0, stddev: float = 1.0, blur_passes: int = 6, roughness: float = 0.55, beta: float = 2.0, n_threads: int = 0) -> lsdtt3.Raster

Generate a deterministic synthetic noise raster.

Parameters
----------
metadata : RasterMeta
    Output grid definition.
kind : {"uniform", "red", "blue", "diamond_square", "power_law"}, default="red"
    Surface generator.
seed : int, default=1
    Reproducible random seed.
mean, stddev : float, default=0, 1
    Target mean and spread; for uniform noise, values span mean ± stddev.
blur_passes, roughness, beta : int, float, float
    Generator-specific correlation or spectral controls.
n_threads : int, default=0
    OpenMP thread count.

Returns
-------
Raster
    Newly allocated synthetic raster.

Authors
-------
B. Gailleton

    Authors
    -------
    B. Gailleton


```python
def synthetic_noise(metadata: lsdtt3.RasterMeta, kind: str = 'red', seed: int = 1, mean: float = 0.0, stddev: float = 1.0, blur_passes: int = 6, roughness: float = 0.55, beta: float = 2.0, n_threads: int = 0) -> lsdtt3.Raster
```
