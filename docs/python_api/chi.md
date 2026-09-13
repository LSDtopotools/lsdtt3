# Chi And Ksn API

This page is generated from the built Python runtime manifest. Do not edit signatures here.

## `ChiCollection`

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

#### `fit`

fit(self, dem: lsdtt3.Raster, critical_divergence_metres: float = 5.0, mainstem: str = 'max_flow_length') -> lsdtt3.ChiFit

Fit piecewise-linear channel profiles in chi-elevation space.

```python
def fit(self, dem: lsdtt3.Raster, critical_divergence_metres: float = 5.0, mainstem: str = 'max_flow_length') -> lsdtt3.ChiFit
```

#### `n_basins`

(self) -> int

#### `nodes`

nodes(self) -> object

Return extracted channel chi nodes as a pandas DataFrame.

```python
def nodes(self) -> object
```

#### `raster`

raster(self) -> lsdtt3.Raster

Return chi values as a masked Raster.

```python
def raster(self) -> lsdtt3.Raster
```

---

## `ChiFit`

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

#### `disorder`

disorder(self) -> object

Return full-network chi-disorder per basin.

```python
def disorder(self) -> object
```

#### `disorder_bootstrap`

disorder_bootstrap(self, replicates: int = 1000, seed: int = 5489, sample_full_tributary_count: bool = True, n_side_branches: int = 3, with_replacement: bool = True) -> object

Bootstrap disorder over side branches.

```python
def disorder_bootstrap(self, replicates: int = 1000, seed: int = 5489, sample_full_tributary_count: bool = True, n_side_branches: int = 3, with_replacement: bool = True) -> object
```

#### `disorder_jackknife`

disorder_jackknife(self) -> object

Leave-one-side-branch-out disorder samples.

```python
def disorder_jackknife(self) -> object
```

#### `disorder_subset_sensitivity`

disorder_subset_sensitivity(self, n_side_branches: int = 3) -> object

Evaluate disorder across mainstem-plus-tributary subsets.

```python
def disorder_subset_sensitivity(self, n_side_branches: int = 3) -> object
```

#### `invert_ksn`

invert_ksn(self, m_over_n_values: object, reference_area: float = 1.0, coarse_cell_size_metres: float = 2000.0, roughness_weight: float = 10000.0, diagonal_jitter: float = 1e-08, use_mae: bool = True) -> lsdtt3.ChiInversion

Run regularised spatial ksn/u-star inversion over m/n trials.

```python
def invert_ksn(self, m_over_n_values: object, reference_area: float = 1.0, coarse_cell_size_metres: float = 2000.0, roughness_weight: float = 10000.0, diagonal_jitter: float = 1e-08, use_mae: bool = True) -> lsdtt3.ChiInversion
```

#### `nodes`

(self) -> object

#### `segments`

(self) -> object

#### `summary`

(self) -> object

---

## `ChiInversion`

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

#### `nodes`

(self) -> object

#### `summary`

(self) -> object

#### `support_raster`

support_raster(self, basin_id: int = 0) -> lsdtt3.Raster

```python
def support_raster(self, basin_id: int = 0) -> lsdtt3.Raster
```

#### `trials`

(self) -> object

#### `u_star_raster`

u_star_raster(self, basin_id: int = 0) -> lsdtt3.Raster

```python
def u_star_raster(self, basin_id: int = 0) -> lsdtt3.Raster
```
