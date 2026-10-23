# Statistics And Utilities API

This page is generated from the built Python runtime manifest. Do not edit signatures here.

## `compare_bearings`

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

```python
def compare_bearings(first: object, second: object, axial: bool = False, rotation_invariant: bool = False) -> dict
```

---

## `helper_for_llm`

Return a compact API guide for AI assistants (and for people in a hurry).

Parameters
----------
topic : str, default="index"
    ``"index"`` (the default) lists the topics and the notation;
    ``"workflow"`` is a runnable core analysis; ``"gotchas"`` lists common
    mistakes; ``"all"`` is every topic, sized for pasting into an
    assistant. ``"symbol:<name>"`` gives the documentation of one public
    symbol and ``"search:<words>"`` ranks public symbols by how many of
    the words their names and docstrings contain.

Returns
-------
str
    Plain text with exact public identifiers. It is a ``str`` subclass
    that notebooks display as text rather than as an escaped literal, so
    ``lsdtt3.helper_for_llm("all")`` as the last line of a cell reads
    well; ``print()`` works too.

Authors
-------
B. Gailleton
S. M. Mudd

```python
(topic='index')
```

---

## `linear_stats`

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

```python
def linear_stats(raster: lsdtt3.Raster, lines: object, spacing_metres: float = 10.0) -> object
```

---

## `relief_ratio`

Return short-scale local relief divided by broad-scale local relief.

Local relief is the elevation range in a circular window; narrow incised
features have high short-range relief compared with their surroundings.

Parameters
----------
raster : Raster
    Elevation raster (m).
short_radius_metres, long_radius_metres : float, default=1.8, 8.0
    Radii (m) of the short and broad windows.

Returns
-------
Raster
    Dimensionless relief ratio.

Authors
-------
B. Gailleton

```python
def relief_ratio(raster: lsdtt3.Raster, short_radius_metres: float = 1.8, long_radius_metres: float = 8.0) -> lsdtt3.Raster
```
