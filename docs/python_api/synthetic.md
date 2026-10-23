# Synthetic Workflows API

This page is generated from the built Python runtime manifest. Do not edit signatures here.

## `SyntheticNetworkResult`

In-memory optimal-channel-network outputs.

Returned by ``synthetic_network()``; not constructible directly.

Attributes
----------
drainage_area_cells : Raster
    Drainage area in cells (Raster).
elevation : Raster
    Steady-state elevation Raster (m) of the network.
final_energy : float
    Network energy after optimisation.
initial_energy : float
    Network energy sum(A^gamma L) before optimisation.
receivers : numpy.ndarray
    Int array (n_rows, n_cols) of each cell's receiver flat index.

Authors
-------
B. Gailleton

### Members

#### `drainage_area_cells`

Drainage area in cells (Raster).

#### `elevation`

Steady-state elevation Raster (m) of the network.

#### `final_energy`

Network energy after optimisation.

#### `initial_energy`

Network energy sum(A^gamma L) before optimisation.

#### `receivers`

Int array (n_rows, n_cols) of each cell's receiver flat index.

---

## `synthetic_hillslope_strip`

Generate a one-dimensional hillslope test DEM and fixed channel mask.

Elevation varies across one axis only, between channel lines at one or both
ends, for testing hillslope solvers against analytic profiles.

Parameters
----------
metadata : RasterMeta
    Output grid.
channel_config : {"both_ends", "single_end_low", "single_end_high"}, default="both_ends"
    Channels at both ends (divide in the middle) or at one end (divide at
    the other).
init_profile : {"flat", "analytic_parabola", "raised_plateau", "perturbed", "linear_ramp"}, default="flat"
    Initial hillslope surface: flat at base level, the analytic steady
    profile, a raised block, flat plus noise, or a triangular ridge.
cross_axis_is_x : bool, default=True
    True: elevation varies east-west and channels are north-south columns;
    False: it varies north-south and channels are rows.
baselevel_elevation : float, default=0.0
    Elevation (m) of the channel cells.
half_width : float, default=0.0
    Hillslope half-width L (m, divide to channel) of the analytic profiles;
    <= 0 derives it from the grid.
uplift_rate : float, default=0.001
    Uplift rate U (m/yr) of the analytic profile.
diffusivity : float, default=0.01
    Diffusivity D (m^2/yr) of the analytic profile.
critical_slope : float, default=0.8
    Critical gradient S_c (m/m) of the nonlinear analytic profile.
nonlinear : bool, default=False
    Use the nonlinear (Roering 2007) analytic profile instead of the linear
    parabola ``z = U / (2 D) (L^2 - x^2)``.
initial_amplitude : float, default=50.0
    Height (m) of the raised plateau, noise or ramp crest.
seed : int, default=1
    Random seed for ``"perturbed"``.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
dict
    ``"elevation"`` (Raster, m) and ``"channel_mask"`` (Raster, 0/1).

Authors
-------
B. Gailleton

```python
def synthetic_hillslope_strip(metadata: lsdtt3.RasterMeta, channel_config: str = 'both_ends', init_profile: str = 'flat', cross_axis_is_x: bool = True, baselevel_elevation: float = 0.0, half_width: float = 0.0, uplift_rate: float = 0.001, diffusivity: float = 0.01, critical_slope: float = 0.8, nonlinear: bool = False, initial_amplitude: float = 50.0, seed: int = 1, n_threads: int = 0) -> dict
```

---

## `synthetic_mountain_range`

Generate a synthetic mountain range and its known mountain mask.

Parameters
----------
metadata : RasterMeta
    Output grid (``RasterMeta.from_shape``).
geometry : {"linear", "halfspace", "circular", "arc"}, default="linear"
    ``linear``: a north-south Gaussian ridge through the centre with
    piedmonts on both sides; ``halfspace``: mountains against one edge;
    ``circular``: an elliptical dome in the middle; ``arc``: mountains
    against one edge with a front that curves into the piedmont.
edge : {"north", "south", "east", "west"}, default="north"
    Edge the mountains are against (``halfspace`` and ``arc`` only).
ridge_elevation : float, default=2000.0
    Crest elevation (m).
piedmont_elevation : float, default=0.0
    Base (piedmont) elevation (m).
front_fraction : float, default=0.2
    Size of the mountain zone as a fraction of the grid extent: the ridge
    half-width (``linear``), the ellipse semi-axes (``circular``) or the
    front's distance from the edge (``halfspace``, ``arc``).
noise_amplitude : float, default=100.0
    Amplitude (m) of the fractal noise added to the surface.
seed : int, default=1
    Random seed.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
dict
    ``"elevation"`` (Raster, m) and ``"mountain_mask"`` (Raster, 1 =
    mountain, 0 = piedmont).

Authors
-------
B. Gailleton

```python
def synthetic_mountain_range(metadata: lsdtt3.RasterMeta, geometry: str = 'linear', edge: str = 'north', ridge_elevation: float = 2000.0, piedmont_elevation: float = 0.0, front_fraction: float = 0.2, noise_amplitude: float = 100.0, seed: int = 1, n_threads: int = 0) -> dict
```

---

## `synthetic_network`

Optimise a synthetic D8 optimal channel network and derive its chi steady-state DEM.

An optimal channel network (OCN) is found by simulated annealing of the
energy ``sum(A**gamma * L)`` over D8 networks, and a steady-state
stream-power DEM is built on it from chi.

Parameters
----------
metadata : RasterMeta
    Output grid.
outlet_mode : {"single_south", "all_coast", "north_south"}, default="single_south"
    One outlet at the centre of the south edge; every edge pixel an outlet;
    or outlets along the north and south edges.
periodic_columns : bool, default=False
    Wrap the grid east-west (route the DEM with ``p`` on east and west).
outlet_spacing : int, default=1
    ``"north_south"`` only: column spacing of the outlets (larger values
    force wider basins).
gamma : float, default=0.5
    Energy exponent on drainage area.
n_iterations : int, default=-1
    Annealing iterations; <= 0 uses 40 x n_rows x n_cols.
temperature_scale : float, default=1.0
    Multiplier on the annealing temperature schedule.
m_exponent, n_exponent : float, default=0.45, 1.0
    Stream-power exponents used to build the DEM from chi.
relief_metres : float, default=1000.0
    Relief (m) of the output DEM.
seed : int, default=1
    Random seed.

Returns
-------
SyntheticNetworkResult

Authors
-------
B. Gailleton
S. M. Mudd

**Deprecated compatibility names:**
- `relief_m` -> `relief_metres`; earliest removal: following minor.

```python
def synthetic_network(metadata: lsdtt3.RasterMeta, outlet_mode: str = 'single_south', periodic_columns: bool = False, outlet_spacing: int = 1, gamma: float = 0.5, n_iterations: int = -1, temperature_scale: float = 1.0, m_exponent: float = 0.45, n_exponent: float = 1.0, relief_metres: float = 1000.0, seed: int = 1) -> lsdtt3.SyntheticNetworkResult
```

---

## `synthetic_noise`

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

```python
def synthetic_noise(metadata: lsdtt3.RasterMeta, kind: str = 'red', seed: int = 1, mean: float = 0.0, stddev: float = 1.0, blur_passes: int = 6, roughness: float = 0.55, beta: float = 2.0, n_threads: int = 0) -> lsdtt3.Raster
```
