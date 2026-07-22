# CLI reference: lsdtt-bearing-comparison

Compares two angular (bearing / orientation) distributions read from CSV columns
and reports the **circular 1-Wasserstein distance** (optimal transport on the
circle) plus the **two-sample Kuiper test** (a rotation-invariant circular
analogue of the Kolmogorov–Smirnov test). It handles both directional data
(period 360°) and axial / orientation data (period 180°). The tests are provided
by the circular functions in `src/utils/lsdtt_statistical_tests.hpp`.

## Usage

```
lsdtt-bearing-comparison <path> <parameter_file>
lsdtt-bearing-comparison <parameter_file>        # path defaults to ./
lsdtt-bearing-comparison -h                      # print help and exit
lsdtt-bearing-comparison -v                      # print version and exit
lsdtt-bearing-comparison -c                      # print citation and exit
```

## Parameter file format

```
# Lines beginning with # are comments.
key: value
```

---

## Input parameters

Each bearing sample is read from a CSV column named by `bearing_column`. The two
files may be different outputs (e.g. two DEMs) but must share the column name.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `bearing_file_1` | string | `NULL` | **Required.** CSV file for the first bearing sample. Must contain the column named by `bearing_column`. |
| `bearing_file_2` | string | `NULL` | **Required.** CSV file for the second bearing sample. Must contain the column named by `bearing_column`. |
| `bearing_column` | string | `bearing_deg` | Name of the CSV column holding bearings in degrees. E.g. `bearing_deg` from a `_channel_link_bearings` or `_junction_com_bearings` output. |

---

## Comparison options

| Parameter | Type | Default | Description |
|---|---|---|---|
| `data_is_axial` | bool | `false` | Treat bearings as axial (orientation) data with period 180°, where θ and θ+180° are equivalent. Use for undirected line orientations. The default (directional) uses period 360°. |
| `rotation_invariant` | bool | `false` | Allow one global rotation of sample 2 before computing the Wasserstein distance, so the comparison reflects distribution **shape** rather than absolute orientation. The chosen rotation is reported as `optimal_rotation_deg`. |
| `stream_order_filter` | int | `-1` | If > 0, keep only rows whose `stream_order_column` equals this value. `-1` (default) uses all rows. |
| `stream_order_column` | string | `stream_order` | Column used by `stream_order_filter`. Consulted only when `stream_order_filter` > 0. |

---

## Read / write path parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | string | directory of parameter file | Directory prepended to relative `bearing_file_1` / `bearing_file_2` paths. |
| `write_path` | string | same as `read_path` | Directory for the output CSV. |
| `write_fname` | string | parameter file stem | Output filename stem. |

---

## Output

Output file: `<write_fname>_bearing_comparison.csv` — a single data row.

| Column | Description |
|---|---|
| `n_x` | Number of bearings read from `bearing_file_1` |
| `n_y` | Number of bearings read from `bearing_file_2` |
| `axial` | 1 if axial (period 180°), 0 if directional (period 360°) |
| `rotation_invariant` | 1 if a global rotation of sample 2 was allowed, else 0 |
| `period_deg` | Angular period used (180 or 360) |
| `circular_wasserstein_deg` | Circular 1-Wasserstein distance, in degrees |
| `optimal_rotation_deg` | Global rotation applied to sample 2 (degrees; meaningful only when `rotation_invariant: true`) |
| `kuiper_V` | Kuiper V statistic |
| `kuiper_p` | Approximate Kuiper p-value |

---

## Interpreting the results

- **Circular Wasserstein distance** gives the *magnitude* of the difference between the two distributions, expressed directly in **degrees**. It is a true metric (optimal transport on the circle) and carries no p-value.
- **Kuiper test** gives the *significance* of the difference. It is the rotation-(origin-)invariant circular analogue of the two-sample KS test, so it is unaffected by the arbitrary choice of zero bearing.
- Use `rotation_invariant: true` to compare distribution **shape** while ignoring absolute orientation (a single global rotation of sample 2 is optimised out).

> **Note.** The *linear* tests in `lsdtt_statistical_tests` (KS, Cramér–von Mises, Anderson–Darling, linear Wasserstein) are **not** appropriate for raw bearings: they ignore the wrap-around at 0°/360° (and 0°/180° for axial data). Always use the circular Wasserstein distance and the Kuiper test for angular data.

---

## Example: compare channel bearings between two DEMs

Uses the `_channel_link_bearings.csv` outputs of `lsdtt-channel-network-analysis`
(which carry a `bearing_deg` column).

```
# lsdtt-bearing-comparison parameter file
read_path:   /data/output/
write_path:  /data/output/
write_fname: dem_pair

bearing_file_1: dem_A_channel_link_bearings.csv
bearing_file_2: dem_B_channel_link_bearings.csv
bearing_column: bearing_deg

data_is_axial:      false
rotation_invariant: false
```

## Example: compare a single stream order, shape only

```
bearing_file_1:      dem_A_channel_link_bearings.csv
bearing_file_2:      dem_B_channel_link_bearings.csv
bearing_column:      bearing_deg
stream_order_filter: 2
stream_order_column: stream_order
rotation_invariant:  true
```
