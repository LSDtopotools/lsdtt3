# Statistics And Utilities API

This page is generated from the built Python runtime manifest. Do not edit signatures here.

## `compare_bearings`

compare_bearings(first: object, second: object, axial: bool = False, rotation_invariant: bool = False) -> dict

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
compare_bearings(first: object, second: object, axial: bool = False, rotation_invariant: bool = False) -> dict

Compare two samples of bearings or axial orientations.

NaN and infinite values are ignored. The result combines the existing native
Kuiper two-sample test with circular Wasserstein distance.

Parameters
----------
first, second : array-like of float
    Bearing samples in degrees.
axial : bool, default=False
    Treat angles 180 degrees apart as equivalent orientations.
rotation_invariant : bool, default=False
    Allow a single global rotation of ``second`` when computing the
    Wasserstein distance. The Kuiper statistic always retains absolute phase.

Returns
-------
dict
    Sample sizes, period, circular Wasserstein distance, optional optimal
    rotation, and Kuiper statistic with its approximate p-value.

Authors
-------
B. Gailleton

Authors
-------
B. Gailleton

```python
def compare_bearings(first: object, second: object, axial: bool = False, rotation_invariant: bool = False) -> dict
```

---

## `helper_for_llm`

Return an ultra-compact, LLM-oriented API card.

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
Return an ultra-compact, LLM-oriented API card.

Parameters
----------
topic : str, default="index"
    Topic name; use ``"symbol:<name>"`` for one public symbol,
    ``"search:<words>"`` for a compact name search, or ``"all"``.

Returns
-------
str
    Machine-oriented API context with exact public identifiers.

Authors
-------
B. Gailleton

Authors
-------
B. Gailleton

```python
(topic='index')
```

---

## `linear_stats`

linear_stats(raster: lsdtt3.Raster, lines: object, spacing_metres: float = 10.0) -> object

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
linear_stats(raster: lsdtt3.Raster, lines: object, spacing_metres: float = 10.0) -> object

Sample a raster along line geometries and append descriptive statistics.

Parameters
----------
raster : Raster
    Raster to sample.
lines : geopandas.GeoDataFrame
    LineString or MultiLineString geometries in the raster CRS.
spacing_metres : float, default=10
    Spacing between samples along each line part.

Returns
-------
geopandas.GeoDataFrame
    Copy of ``lines`` with length, sample-count, and descriptive-statistic
    columns. Samples outside the grid or on nodata are ignored.

Authors
-------
B. Gailleton

Authors
-------
B. Gailleton

```python
def linear_stats(raster: lsdtt3.Raster, lines: object, spacing_metres: float = 10.0) -> object
```

---

## `relief_ratio`

relief_ratio(raster: lsdtt3.Raster, short_radius_metres: float = 1.8, long_radius_metres: float = 8.0) -> lsdtt3.Raster

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
relief_ratio(raster: lsdtt3.Raster, short_radius_metres: float = 1.8, long_radius_metres: float = 8.0) -> lsdtt3.Raster

Return short-scale local relief divided by broad-scale local relief.

Authors
-------
B. Gailleton

```python
def relief_ratio(raster: lsdtt3.Raster, short_radius_metres: float = 1.8, long_radius_metres: float = 8.0) -> lsdtt3.Raster
```
