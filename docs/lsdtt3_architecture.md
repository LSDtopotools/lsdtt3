# lsdtt3 Architecture Document

Paste this at the start of a new session to restore full context.

Last updated: session 8 (lsdtt_flow_engine→lsdtt_flow_base; lsdtt_flowinfo→lsdtt_flowinfo_cpu; lsdtt_flowgraph→lsdtt_flowgraph_cpu; make_flow_base factory moved to lsdtt_flow_base.cpp in lsdtt3_cpu to avoid circular dependency)

---

## Project identity

- **Package name:** lsdtt3
- **Repository:** https://github.com/LSDtopotools/lsdtt3
- **Copyright:** Simon M. Mudd
- **Licence:** GPL-3.0
- **C++ standard:** C++17
- **Build system:** CMake 3.15+

---

## Directory structure

```
lsdtt3/
├── src/
│   ├── core/
│   │   ├── lsdtt_raster_meta.hpp          — metadata struct + CRS validation
│   │   ├── lsdtt_raster_base.hpp/.cpp     — abstract raster dispatch interface
│   │   ├── lsdtt_data_loader.hpp/.cpp     — GDAL/OGR I/O, RAII, COG output
│   │   ├── lsdtt_backend_factory.hpp/.cpp — make_raster_engine(), gpu_available()
│   │   ├── lsdtt_param_parser.hpp/.cpp    — parameter file parser
│   │   ├── lsdtt_flow_types.hpp           — shared flow types (header-only)
│   │   ├── lsdtt_flow_base.hpp           — abstract D8 flow interface (header-only)
│   │   └── lsdtt_flowgraph_cpu.hpp                      — D-inf/MFD DAG routing (stub)
│   ├── cpu/
│   │   ├── lsdtt_raster_cpu.hpp/.cpp                    — concrete CPU raster class
│   │   ├── lsdtt_flowinfo_cpu.hpp/.cpp                      — D8 CPU flow routing
│   │   ├── lsdtt_raster_surface_metrics_cpu.hpp/.cpp    — hillshade, moving window
│   │   ├── lsdtt_raster_flow_preprocessing_cpu.hpp/.cpp — fill, breach, nodata
│   │   ├── lsdtt_raster_nodata_filling_cpu.hpp/.cpp     — interior void filling
│   │   ├── lsdtt_raster_spectral_cpu.hpp/.cpp           — wavelet quadratic fit
│   │   └── lsdtt_raster_synthetic_surfaces_cpu.hpp/.cpp — test surface generators
│   ├── gpu/
│   │   ├── lsdtt_raster_gpu.hpp/.cu       — GPU raster (LSDTT_GPU=ON)
│   │   └── lsdtt_flowinfo_gpu.hpp/.cpp    — GPU flow routing (future)
│   └── utils/
│       └── lsdtt_array2d.hpp              — header-only 2D array
├── cli/
│   ├── lsdtt-raster-preprocessing.cpp
│   ├── lsdtt-raster-synthetic-tests.cpp
│   └── lsdtt-test-flowinfo.cpp
├── bin/
├── CMakeLists.txt
└── README.md
```

---

## CMakeLists.txt summary

Two static libraries: `lsdtt3_core` and `lsdtt3_cpu`.

**`lsdtt3_core` sources** (`src/core/`):
- `lsdtt_raster_base.cpp`
- `lsdtt_data_loader.cpp`
- `lsdtt_backend_factory.cpp`
- `lsdtt_param_parser.cpp`
- *(removed — lsdtt_flow_base has no non-inline methods; factory moved to lsdtt3_cpu)

**`lsdtt3_cpu` sources** (`src/cpu/`):
- `lsdtt_raster_cpu.cpp`
- `lsdtt_raster_surface_metrics_cpu.cpp`
- `lsdtt_raster_flow_preprocessing_cpu.cpp`
- `lsdtt_raster_nodata_filling_cpu.cpp`
- `lsdtt_raster_spectral_cpu.cpp`
- `lsdtt_raster_synthetic_surfaces_cpu.cpp`
- `lsdtt_flowinfo_cpu.cpp`
- `lsdtt_flow_base.cpp`

`lsdtt3_cpu` depends on `lsdtt3_core`. All CLI executables link against
`lsdtt3_cpu` and `lsdtt3_core`. OpenMP and FFTW link to `lsdtt3_cpu`.

**Dependencies:**
- GDAL — required
- OpenMP — optional (warning if absent)
- FFTW3 (double, `fftw3`) — optional via pkg-config
- CUDA — optional (`-DLSDTT_GPU=ON`)

**FFTW:** Uses double API (`fftw_*`), NOT threaded FFTW, NOT float (`fftwf_*`).
Compile definition: `LSDTT_HAS_FFTW=1`. Stub in `lsdtt_raster_spectral_cpu.cpp`
throws at runtime if FFTW not available.

---

## Naming conventions

- All snake_case, `lsdtt_` prefix on all classes and files
- No camelCase anywhere
- Opening braces always on their own line
- All parameters use snake_case with no spaces
- Output raster parameters: `write_*_raster`
- Algorithm selection parameters: `*_algorithm_choice`
- Output format parameter: `write_output_format` (`cog` / `gtiff` / `same_as_input`)

**CLI executables use dashes, not underscores:**
- Object/library: `lsdtt_flowinfo_cpu`, `lsdtt_raster_cpu`, `lsdtt_param_parser`
- CLI executable: `lsdtt-test-flowinfo`, `lsdtt-raster-preprocessing`

**Output filename convention:**
`<write_fname>_<descriptor>.tif`. If `write_fname` not set, defaults to stem
of `read_fname`. E.g. input `my_DEM.tif` → `my_DEM_n_contributing_pixels.tif`.

**File suffix convention:**
CPU-specific analysis files use `_cpu` suffix: `lsdtt_raster_spectral_cpu.hpp`.
The CPU raster class (`lsdtt_raster_cpu`) and CPU flow routing (`lsdtt_flowinfo_cpu`)
do not use the suffix because they are the primary named objects, not variants.

---

## Architecture design principles

**`src/core/`** — abstract interfaces, factories, shared types, I/O.
Nothing backend-specific. No CPU arrays. No CUDA.

**`src/cpu/`** — everything that runs on host CPU: the CPU raster class,
all CPU analysis free functions, and the CPU flow routing object.

**`src/gpu/`** — GPU implementations. Currently `lsdtt_raster_gpu` placeholder.
Future: `lsdtt_flowinfo_gpu`.

**`src/utils/`** — header-only utilities with no backend dependency.

**Abstract interface pattern:**
```
src/core/lsdtt_raster_base    ← abstract interface for raster backends
src/cpu/lsdtt_raster_cpu      ← CPU implementation
src/gpu/lsdtt_raster_gpu      ← GPU implementation (placeholder)

src/core/lsdtt_flow_base      ← abstract interface for flow routing
src/cpu/lsdtt_flowinfo_cpu      ← CPU implementation
src/gpu/lsdtt_flowinfo_gpu    ← GPU implementation (future)
```

**Factory pattern:**
```
lsdtt_backend_factory  →  lsdtt_raster_cpu or lsdtt_raster_gpu
make_flow_base()       →  lsdtt_flowinfo_cpu or lsdtt_flowinfo_gpu (future)
```

**CLIs hold abstract interfaces via unique_ptr — never concrete types:**
```cpp
// Flow routing: include lsdtt_flow_base.hpp, never lsdtt_flowinfo_cpu.hpp
auto flow = lsdtt::make_flow_base(dem_data, dem_mask, meta, opts);
// returns unique_ptr<lsdtt_flow_base>
// currently wraps lsdtt_flowinfo_cpu; future: lsdtt_flowinfo_gpu
flow->n_valid_cells();
flow->upslope_nodes(flat);
```

Any CLI written against `lsdtt_flow_base` works without modification when the
GPU flow backend is added — only the factory changes.

---

## lsdtt_raster_meta

Plain struct in `src/core/lsdtt_raster_meta.hpp`.
Requires: `<cmath>`, `<stdexcept>`, `<string>`, `<gdal.h>`, `<ogr_spatialref.h>`

Key fields: `n_rows`, `n_cols`, `x_min`, `y_min`, `cell_size`, `no_data_value`,
`crs_wkt`, `epsg_code`, `data_type`.

**CRS validation methods:**
- `is_geographic()` — true if lat/lon CRS
- `linear_unit_name()` — returns e.g. "metre", "foot"
- `is_in_metres()` — true if linear unit is metres (within 1e-6)
- `check_projected_metres()` — throws `std::runtime_error` with messages:
  - *"You are using a geographic coordinate system but this program only works
    on projected data. Please convert your raster to a projected coordinate
    system like UTM."*
  - *"Your units are not in metres. Are you using feet? Please convert your
    raster to a unit that is not only used by lunatics."*

---

## lsdtt_raster_base

File: `src/core/lsdtt_raster_base.hpp/.cpp`

Pure dispatch interface. Does NOT expose `get_data()` or `get_mask()` —
those are CPU-specific types that have no meaning on a GPU backend.
Only `get_meta()` is on the base (metadata is always host-side).

CPU analysis free functions take `const lsdtt_raster_cpu&` directly — they
never go through the base class for data access.

**Pure virtual:** `get_meta()`, `write_to_loader()`, `read_from_loader()`,
`smooth_gaussian()`, `smooth_bilateral()`

**CPU fallback virtuals** (base implements, GPU overrides when ready):
`compute_flow_direction()`, `compute_flow_accumulation()`, `compute_fft()`,
`apply_wavelet_filter()`

---

## lsdtt_raster_cpu

File: `src/cpu/lsdtt_raster_cpu.hpp/.cpp`

**Storage:** `lsdtt_array2d_f` and `lsdtt_array2d_u8` (each wraps `std::vector`).

**Accessors:**
- `get_data_array()` → `const lsdtt_array2d_f&` — preferred for analysis, zero-copy
- `get_mask_array()` → `const lsdtt_array2d_u8&` — zero-copy
- `get_data()` → `const std::vector<float>&` — zero-copy into array's vector
- `get_mask()` → `const std::vector<uint8_t>&` — zero-copy
- `get_meta()` → `const lsdtt_raster_meta&` — implements base pure virtual

**Constructors:**
- `lsdtt_raster_cpu(const lsdtt_data_loader&)` — read full raster from GDAL
- `lsdtt_raster_cpu(const lsdtt_raster_tile&)` — from loader tile
- `lsdtt_raster_cpu(lsdtt_array2d_f&&, lsdtt_array2d_u8&&, meta)` — preferred for results
- `lsdtt_raster_cpu(std::vector<float>&&, std::vector<uint8_t>&&, meta)` — backwards compat

---

## lsdtt_param_parser

File: `src/core/lsdtt_param_parser.hpp/.cpp`

Parameter file format: `key: value` (snake_case keys, `#` for comments).

**File I/O parameters:**
- `read_fname` — full filename **including extension** (e.g. `my_dem.tif`)
- `write_fname` — stem only, no extension (defaults to read stem)
- `read_path` — directory (default: parameter file directory)
- `write_path` — directory (default: read_path)
- `boundary_conditions` — string of characters: `o`=open, `n`=no-flux, `p`=periodic

`get_full_read_path()` returns `read_path + read_fname + read_ext`.
`get_full_write_path()` returns `write_path + write_fname` (no extension).

**Static helpers:** `stem(filename)`, `extension(filename)` (returns dot + ext).

---

## lsdtt_data_loader

File: `src/core/lsdtt_data_loader.hpp/.cpp`

**`raster_output_format` enum:**
- `cog` — Cloud Optimised GeoTIFF (default; falls back to gtiff < 262144 px)
- `gtiff` — plain compressed GeoTIFF (tiled DEFLATE)
- `same_as_input` — uses the GDAL driver of the input raster

**Key methods:**
- `write_full(data, meta, path, format=cog)`
- `finalise_output(path)` — format-aware
- `create_output_raster(path, meta, format=cog, n_bands=1)`

---

## Flow routing — how to navigate the files

There are several flow-routing files. This table explains which one you need.

| You are... | Include this | Notes |
|---|---|---|
| Writing analysis code (chi, steepness, drainage area, etc.) | `lsdtt_flowinfo_cpu.hpp` | Only file you need. Construct `lsdtt_flowinfo`, call its methods directly. |
| Writing a CLI that builds a flow object | `lsdtt_flow_base.hpp` | Call `make_flow_base()`. Returns `unique_ptr<lsdtt_flow_base>`. GPU-compatible. |
| Writing D-inf / MFD distributed-flow analysis | `lsdtt_flowgraph_cpu.hpp` | Stub — implementation pending. |
| Adding a GPU flow routing backend | `lsdtt_flow_base.hpp` + new `lsdtt_flowinfo_gpu.hpp` | Implement `lsdtt_flow_base`, add branch in `src/cpu/lsdtt_flow_base.cpp`. |
| Changing `flat_index_type` (e.g. to int64_t) | `lsdtt_flow_types.hpp` | Change one typedef. All aliases follow. |

**Analysis code** never needs to know about `lsdtt_flow_base`, `lsdtt_flow_types`,
or `lsdtt_flow_base.cpp`. Those are infrastructure. The comment block at the top
of `lsdtt_flowinfo_cpu.hpp` repeats this guidance in-source.

---

## lsdtt_flow_types

File: `src/core/lsdtt_flow_types.hpp` (header-only, no `.cpp`)

Shared types included by both `lsdtt_flow_base.hpp` and `lsdtt_flowinfo_cpu.hpp`
to avoid circular includes. Neither of those headers includes the other — both
include this one.

**Include chain (acyclic):**
```
lsdtt_flow_types.hpp         (no lsdtt3 includes — just <cstdint> etc.)
       ↑
lsdtt_flow_base.hpp        (includes lsdtt_flow_types.hpp)
       ↑
lsdtt_flowinfo_cpu.hpp       (includes lsdtt_flow_types.hpp + lsdtt_flow_base.hpp)
       ↑
lsdtt_flow_base.cpp       (in lsdtt3_cpu — includes lsdtt_flowinfo_cpu.hpp, no cycle)
```

**Contents:**
- `flat_index_type = int32_t` — safe for grids up to ~46,340 × 46,340 (~2.15B cells).
  **To support larger grids: change this one typedef to `int64_t`** — all aliases
  and `k_max_flat_index` throughout the codebase track automatically.
- `donor_count_type = uint8_t` — max 8 donors/incoming edges on regular grid;
  correct for both D8 and MFD Kahn ordering
- `s_vector_size_type = uint32_t`
- `k_nodata_receiver = -1`, `k_outlet_receiver = -2`, `k_no_direction = 255`
- `k_max_flat_index` — max representable cell index for current `flat_index_type`
- `lsdtt_flat_index_slice` — zero-copy view into s_vector; returned by
  `upslope_nodes()`, valid for lifetime of the flow engine object
- `lsdtt_flowinfo_options` — construction options shared by flow routing objects

**`lsdtt_flowinfo_options` fields:**
```
boundary_conditions       = "oooo"   # north east south west
interior_nodata_is_outlet = true
require_acyclic_graph     = true
build_topological_order   = true
build_s_vector            = true
n_threads                 = 0
```
Boundary codes: `o`=open, `n`=no-flux, `p`=periodic (symmetric pairs),
`b`=forced outlet on that entire edge.

---

## lsdtt_flow_base

File: `src/core/lsdtt_flow_base.hpp` (header-only — factory impl in src/cpu/lsdtt_flow_base.cpp)

Abstract base class for D8 flow routing backends. Same role as
`lsdtt_raster_base` for raster backends.

**CLIs include this header, never `lsdtt_flowinfo_cpu.hpp` directly.**

**Type aliases** (from `lsdtt_flow_types.hpp`):
`flat_index_type`, `donor_count_type`, `s_vector_size_type`

**Sentinel values** (from `lsdtt_flow_types.hpp`):
`nodata_receiver = -1`, `outlet_receiver = -2`, `no_direction = 255`

**Pure virtual methods:** full set of metadata accessors, index conversion,
cell validity, routing queries (receiver, direction, donor count), s_vector
queries (start, size, upslope_nodes, is_upslope_of), and iteration helpers.

**`for_each_upslope_node`:** takes `const std::function<void(flat_index_type)>&`
not a template, to enable virtual dispatch. Small runtime cost at the
CLI/analysis boundary — acceptable.

**Factory function `make_flow_base()`**:
Declared in `lsdtt_flow_base.hpp`, implemented in `lsdtt_flow_base.cpp` (lsdtt3_cpu). Currently returns `std::make_unique<lsdtt_flowinfo_cpu>`.
When GPU backend is implemented: add `gpu_available()` branch here.
No other files change.

**Adding GPU flow routing in future:**
1. Write `src/gpu/lsdtt_flowinfo_gpu.hpp/.cpp` implementing `lsdtt_flow_base`
2. Add `gpu_available()` branch in `lsdtt_flow_base.cpp`
3. Add `lsdtt_flowinfo_gpu.cpp` to `lsdtt3_gpu` in CMakeLists
4. All existing CLIs work unchanged

---

## lsdtt_flowinfo_cpu

File: `src/cpu/lsdtt_flowinfo_cpu.hpp/.cpp`

CPU-only D8 flow routing. Inherits `lsdtt_flow_base`. All storage is
host-side (`lsdtt_array2d`, `std::vector`).
Requires a pit-free, flat-resolved DEM — use `lsdtt_raster_flow_preprocessing`
first.

**Direction code convention (receiver_dir_):**
```
7  0  1
6  x  2
5  4  3
```

**Build pipeline (in order):**
1. `validate_inputs()` — dimensions, metadata, overflow guard using `k_max_flat_index`
2. `parse_boundary_conditions()` — normalises, enforces periodic symmetry
3. `initialise_storage()` — allocates all arrays, counts valid cells
4. `build_receiver_graph_d8()` — steepest-descent D8; **safe to parallelise**
   with `#pragma omp parallel for collapse(2) schedule(static)` (no write conflicts)
5. `build_donor_count_d8()` — serial; parallel requires atomics on receiver
6. `build_outlets()` — row-major scan, serial
7. `build_topological_order()` — Kahn's BFS, serial
8. `build_s_vector()` — iterative DFS; serial. Parallel two-pass strategy
   (subtree sizes → prefix sum → fill) documented in source comments

**s_vector contract:**
- `s_vector_[s_vector_start_[i]] == i` — node is always first in its own slice
- `s_vector_[start : start+size)` — all upslope nodes of `i` including `i`,
  preorder DFS, donors sorted ascending by flat row-major index at each level
- `is_upslope_of(candidate, target)` — O(1) interval membership test
- `upslope_nodes(flat)` — returns zero-copy `lsdtt_flat_index_slice`

**Why s_vector is D8-only:**
D8 produces a forest of trees — one receiver per cell — so the upslope set of
any node is a strict subtree, representable as a contiguous interval. In
D-inf/MFD the graph is a DAG, upslope sets overlap and are fractional, and the
interval representation breaks. See `lsdtt_flowgraph_cpu` for the DAG approach.

---

## lsdtt_flowgraph_cpu

File: `src/cpu/lsdtt_flowgraph_cpu.hpp` (stub — implementation pending)

DAG-based flow routing for D-infinity and MFD. No s_vector.

**Data model (CSR):**
- `receiver_offsets_` — `n_cells+1` prefix array
- `receiver_indices_` — flat receiver indices
- `receiver_weights_` — `float` fractional weights summing to 1.0 per valid cell
- Cell `i`'s receivers: `[receiver_offsets_[i], receiver_offsets_[i+1])`

**Routing methods:** `d_infinity` (Tarboton 1997), `mfd` (Quinn 1991)

`flat_index_type` from `lsdtt_flow_types.hpp`. `donor_count_type = uint8_t`
correct for Kahn ordering. MFD weight accumulation uses a separate `float` array.

---

## CPU analysis free functions (src/cpu/)

All take `const lsdtt_raster_cpu&` directly — not `lsdtt_raster_base`.
All return `lsdtt_raster_cpu`. Input is never modified.

**`lsdtt_raster_surface_metrics_cpu`:**
- `compute_hillshade(raster, azimuth=315, altitude=45, z_factor=1.0)`
  — Zevenbergen & Thorne (1987), float32 in [0,255]
- `moving_window_mean(raster, radius=3, shape=square)`
- `window_shape` enum: `square`, `circular`

**`lsdtt_raster_flow_preprocessing_cpu`:**
- `mask_to_nodata_using_threshold(raster, threshold, is_below_nodata)`
- `fill_pits(raster, min_slope=0.0001f, n_threads=0)` — Barnes et al. 2014
- `fill_pits_with_boundary(...)` — serial only
- `breach_depressions(raster, n_threads=0)` — Lindsay 2016
- `breach_depressions_least_cost(raster, max_dist=1000, max_cost=max_float,
  min_dist=true, flat_increment=0, fill_remainder=true, n_threads=0)` — Lindsay 2020

**`lsdtt_raster_nodata_filling_cpu`:**
- `identify_interior_nodata(raster)` → `lsdtt_array2d<uint8_t>` (1=interior void)
- `fill_interior_nodata(raster, method=laplacian, max_iter=10000,
  tolerance=1e-6, n_threads=0)`
- `clip_to_maximum_rectangle(raster)` — O(n) largest-rectangle-in-histogram
- `nodata_fill_method` enum: `iterative_weighted_average`, `laplacian`

**`lsdtt_raster_spectral_cpu`:**
- `wavelet_quadratic_fit(raster, lambda=100.0, n_threads=0)`
  — FFTW double API, Gaussian kernel, conservative nodata masking
  — Returns `poly_fit_coeffs{a, b, c, d, e, f}` — six `lsdtt_raster_cpu` members
  — Fits: `z(x,y) = a*x² + b*y² + c*x*y + d*x + e*y + f`
  — Throws if FFTW not available at build time

**`lsdtt_raster_synthetic_surfaces_cpu`:**
- Takes `lsdtt_raster_meta` or raw dimensions — no raster input, no signature change needed
- Generators: uniform random, red noise, blue noise, diamond-square, periodic wave,
  global quadratic, random quadratic
- All stochastic generators deterministic for a given seed

---

## lsdtt-test-flowinfo CLI

File: `cli/lsdtt-test-flowinfo.cpp`

Builds D8 flow routing from a DEM and writes diagnostic rasters. Uses
`make_flow_base()` — compatible with
future GPU backend without modification.

**Includes** `lsdtt_flow_base.hpp`,
NOT `lsdtt_flowinfo_cpu.hpp`.

**Construction pattern:**
```cpp
std::unique_ptr<lsdtt::lsdtt_flow_base> flow =   
    lsdtt::make_flow_base(dem_data, dem_mask, meta, fi_opts);
flow->n_valid_cells();
flow->s_vector_size_unchecked(flat);
```

**Key parameters:**
```
read_fname: my_dem.tif
write_fname: my_dem
read_path: /data/
write_path: /output/
boundary_conditions: oooo
breach_before_fill: false
fill_remainder_after_breach: true
fill_before_routing: true
fill_algorithm_choice: barnes2014
breach_algorithm_choice: lindsay2016
fill_min_slope: 0.0001
minimum_elevation: 0.0
maximum_elevation: 30000.0
write_n_contributing_pixels_raster: false
write_flow_direction_raster: false
write_donor_count_raster: false
write_outlet_raster: false
write_fill_raster: false
write_breached_raster: false
write_output_format: cog
n_threads: 0
```

**Flow direction codes:**
`0`=N `1`=NE `2`=E `3`=SE `4`=S `5`=SW `6`=W `7`=NW; outlets = `-2`

---

## lsdtt-raster-preprocessing CLI

File: `cli/lsdtt-raster-preprocessing.cpp`

**Workflow order:** threshold mask → nodata fill → clip rectangle →
breach (optional) → fill → hillshade.

**Key parameters:**
```
read_fname: my_dem.tif
write_fname: my_dem
read_path: /data/
write_path: /output/
write_output_format: cog
write_hillshade_raster: false
write_fill_raster: false
fill_algorithm_choice: barnes2014
fill_min_slope: 0.0001
breach_before_fill: false
fill_remainder_after_breach: true
write_breached_raster: false
breach_algorithm_choice: lindsay2016
breach_lc_max_dist: 1000
breach_lc_max_cost: unlimited
write_filled_nodata_raster: false
nodata_fill_method_choice: laplacian
write_clipped_rectangle_raster: false
minimum_elevation: 0.0
maximum_elevation: 30000.0
hillshade_azimuth: 315.0
hillshade_altitude: 45.0
hillshade_z_factor: 1.0
n_threads: 0
```

---

## Known compile errors fixed (for reference)

| File | Issue | Fix |
|------|-------|-----|
| `lsdtt_raster_meta.hpp` | Missing `<stdexcept>` | Added include |
| `lsdtt_raster_meta.hpp` | Unused `to_metres` variable | Removed |
| `lsdtt_raster_flow_preprocessing_cpu.cpp` | Missing `<iostream>` | Added |
| `lsdtt_raster_flow_preprocessing_cpu.cpp` | Several unused params | Removed |
| `lsdtt_raster_flow_preprocessing_cpu.cpp` | `dest_backlink` set but unused | Removed |
| `lsdtt_raster_flow_preprocessing_cpu.cpp` | Operator precedence `&&` in `\|\|` | Added parens |
| `lsdtt_raster_nodata_filling_cpu.cpp` | `fill_laplacian` unused `meta` param | Removed |
| `lsdtt_raster_nodata_filling_cpu.cpp` | `fill_interior_nodata` unused `nodata` | Removed |
| `lsdtt_raster_nodata_filling_cpu.cpp` | No matching constructor for return | Fixed |
| `lsdtt-raster-preprocessing.cpp` | `applied_threshold` set but unused | Removed |
| `lsdtt-raster-preprocessing.cpp` | `const lsdtt_raster_base&` to cpu fn | Changed to `const lsdtt_raster_cpu&` |

---

## Key design principles

1. No `using namespace std` anywhere
2. All code in `namespace lsdtt`
3. `check_projected_metres()` called at entry of metric analysis functions
4. Mask array (`uint8_t`, 0=nodata, 1=valid) separate from data array
5. `meta.is_nodata(v)` handles both sentinel values and NaN
6. No warnings — `-Wall -Wextra -Wpedantic` enforced
7. FFTW double (`fftw_*`), not float; not threaded FFTW
8. Analysis free functions take `const lsdtt_raster_cpu&` — not the base class
9. CLIs hold flow routing via `unique_ptr<lsdtt_flow_base>` from factory —
   never construct `lsdtt_flowinfo_cpu` directly in a CLI
10. `src/core/` = abstract interfaces, factories, shared types, I/O only
11. `src/cpu/` = all CPU implementations (raster and flow)
12. `src/gpu/` = all GPU implementations
