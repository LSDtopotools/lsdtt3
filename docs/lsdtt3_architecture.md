# lsdtt3 Architecture Document

A guide to the architecture of lsdtt3. 
Use for context when interacting with AI  during code development. 
Can also be used as a developer guide.

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
│   │   ├── lsdtt_junction_network_base.hpp  - abstract channel network interface (header only)
│   │   ├── lsdtt_junction_network_types.hpp - shared junction network types (header-only)
│   │   ├── lsdtt_flow_base.hpp            — abstract D8 flow interface (header-only)
│   │   ├── lsdtt_flowgraph_cpu.hpp        — D-inf/MFD DAG routing (stub)
│   │   ├── lsdtt_crs_suggest.hpp/.cpp     — CRS inspection and projection recommendation
│   │   ├── lsdtt_raster_reproject.hpp/.cpp  — CRS reprojection
│   │   └── lsdtt_raster-_reproject.hpp/.cpp — GDAL-warp reprojection + grid snapping
│   ├── cpu/
│   │   ├── lsdtt_raster_cpu.hpp/.cpp                    — concrete CPU raster class
│   │   ├── lsdtt_flowinfo_cpu.hpp/.cpp                  — D8 CPU flow routing
│   │   ├── lsdtt_raster_algebra.hpp                     — element-wise raster algebra (header-only)
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
├── docs/
│   └── preprocessing.md                  — user guide for preprocessing CLI
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
- `lsdtt_crs_suggest.cpp`
- `lsdtt_raster_reproject.cpp`
- `lsdtt_vector_utils.cpp`

**`lsdtt3_cpu` sources** (`src/cpu/`):
- `lsdtt_raster_cpu.cpp`
- `lsdtt_raster_algebra.cpp`
- `lsdtt_raster_surface_metrics_cpu.cpp`
- `lsdtt_raster_flow_preprocessing_cpu.cpp`
- `lsdtt_raster_nodata_filling_cpu.cpp`
- `lsdtt_raster_spectral_cpu.cpp`
- `lsdtt_raster_synthetic_surfaces_cpu.cpp`
- `lsdtt_flowinfo_cpu.cpp`
- `lsdtt_flow_base.cpp`
- `lsdtt_junction_network_cpu.cpp`
- `lsdtt_raster_opencv_cpu.cpp` *(when OpenCV is available)*

`lsdtt3_cpu` depends on `lsdtt3_core`. All CLI executables link against
`lsdtt3_cpu` and `lsdtt3_core`. OpenMP, FFTW, and OpenCV link to `lsdtt3_cpu`.

**Dependencies:**
- GDAL — required
- OpenMP — optional
- FFTW3 — optional
- OpenCV — optional
- CUDA — optional


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

**Projection suffix convention** (applied by `lsdtt_raster_reproject`):
- `_utm`     — UTM zone of raster centre
- `_polar`   — NSIDC polar stereographic (|lat| > 75°)
- `_lambert` — custom Lambert Azimuthal Equal Area (wide-extent rasters, opt-in)

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
src/cpu/lsdtt_flowinfo_cpu    ← CPU implementation
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

## CLI executables

Current CLI programs are:

- **`lsdtt-raster-preprocessing.cpp`**  
  Raster preprocessing workflow: projection enforcement/reprojection, masking,
  nodata filling, breaching, pit filling, hillshade, and related outputs.

- **`lsdtt-raster-synthetic-tests.cpp`**  
  Synthetic raster generation and testing utilities.

- **`lsdtt-test-flowinfo.cpp`**  
  Builds D8 flow routing and writes diagnostic flow-routing rasters.

- **`lsdtt-surface-metrics.cpp`**  
  Surface-analysis workflow that does not require flow routing: hillshade,
  polynomial surface fit, Otsu segmentation, Hough line extraction.

- **`lsdtt-channel-extraction.cpp`**  
  Builds or reuses a pit-filled DEM, constructs D8 flow routing, extracts source
  pixels from a contributing-pixel threshold, builds the junction/channel
  network, and writes channel-network output to CSV, GeoPackage, or FlatGeobuf.

### CLI design rule

Every CLI uses `lsdtt_param_parser`.

CLI source files should follow a common layout:
1. usage/banner helpers
2. command-line flag parsing
3. parameter defaults grouped by theme
4. `help_map` entries immediately after each parameter declaration
5. parameter parsing
6. cross-parameter logic
7. input loading
8. processing pipeline
9. output writing
10. `print_parameters()` for reproducibility

### Parameter declaration rule

In CLI files:
- parameters should **not** be declared `const`
- parameters should be grouped by workflow theme
- every `help_map` entry should appear **immediately after** the corresponding
  parameter default declaration

This pattern is now the standard for all CLIs because it makes developer-facing
options easy to inspect and keeps the parameter contract readable in one place.

Try to avoid helper functions at the top of CLIs. See if helpers have a more logical position
within /cpu or /utils files. If necessary declare a new .hpp/.cpp pair instead of
helper functions since it is very likeley that a helper in one CLI will be recycled by another.


---

## lsdtt_crs_suggest

Files: `src/core/lsdtt_crs_suggest.hpp/.cpp`

Inspects a `lsdtt_raster_meta` and returns a `lsdtt_crs_suggestion` struct
describing the appropriate projected CRS for analysis. Does **not** reproject —
that is the job of `lsdtt_raster_reproject`.

**Decision logic:**
1. Already projected (non-geographic, non-Web Mercator) → return `already_projected`
2. All four corners are unprojected to WGS84 lon/lat via `OGRCoordinateTransformation`
3. `|centre_lat| > 75°` → NSIDC polar stereographic (EPSG:3413 north / EPSG:3976 south)
4. Otherwise → UTM zone of centre longitude + hemisphere
5. If raster spans ≥ 3 UTM zones: set `warn_multi_zone`, populate `laea_wkt`

**Key types:**
```cpp
enum class CrsRecommendationType {
    already_projected, polar_north, polar_south, utm, laea_custom
};

struct lsdtt_crs_suggestion {
    CrsRecommendationType type;
    std::string recommended_wkt;    // WKT for primary recommendation
    std::string recommended_name;   // human-readable label
    int         recommended_epsg;   // 0 for custom CRS
    double centre_lon, centre_lat;
    double min_lon, max_lon, min_lat, max_lat;
    int    utm_zone_min, utm_zone_max, utm_zone_count;
    bool   already_projected;
    bool   warn_multi_zone;
    std::string warning_message;
    std::string laea_wkt;           // populated when warn_multi_zone && offer_laea
    std::string laea_name;
};
```

**Public API:**
- `lsdtt_suggest_crs(meta, offer_laea=true)` → `lsdtt_crs_suggestion`
- `lsdtt_print_crs_suggestion(sug)` — pretty-print to stdout

Web Mercator (EPSG:3857) is treated as unsuitable alongside EPSG:4326 — both
trigger the projection recommendation path.

---

## lsdtt_raster_reproject

Files: `src/core/lsdtt_raster_reproject.hpp/.cpp`

Reprojects a raster file on disk to a suitable projected CRS using `GDALWarp()`.
Calls `lsdtt_suggest_crs()` internally. Returns a `lsdtt_reproject_result`
containing the output path, the CRS applied, and the grid spacing used.

**Key types:**
```cpp
struct lsdtt_reproject_options {
    bool        enforce_projection      = true;
    std::string resampling_method       = "cubic";
    double      grid_spacing            = -1.0;  // ≤ 0 = auto
    bool        allow_lambert_projection = false;
};

struct lsdtt_reproject_result {
    bool        reprojected;
    std::string out_path;        // empty if !reprojected
    lsdtt_crs_suggestion suggestion;
    double      grid_spacing_used;
    std::string target_wkt;
};
```

**Grid spacing snapping** (`snap_grid_spacing(pixel_m)`):

| Estimated pixel size | Snapping rule |
|---|---|
| ≥ 20 m | nearest 5 m |
| 5–20 m | nearest 2 m |
| < 5 m  | nearest 1 m |
| rounds to 0 | preserve raw pixel size (sub-metre guard) |
| ≤ 0 after guard | clamp to 0.01 m (malformed metadata safety net) |

**Output naming:** `write_path + write_stem + suffix + ".tif"`
where suffix is `_utm`, `_polar`, or `_lambert`.

**Resampling methods** (any GDAL-accepted string):
`near`, `bilinear`, `cubic` (default), `cubicspline`, `lanczos`,
`average`, `mode`, `max`, `min`, `med`, `q1`, `q3`, `sum`.

**LAMBERT / LAEA path:**
Only activated when `allow_lambert_projection == true` AND the raster spans
≥ 3 UTM zones. The custom LAEA is centred on (centre_lon, centre_lat) of the
raster. Output suffix: `_lambert`.

**Public API:**
- `lsdtt_reproject_if_needed(in_path, write_path, write_stem, meta, opts)` → `lsdtt_reproject_result`
- `snap_grid_spacing(pixel_m)` → `double`

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

`lsdtt_data_loader` handles both raster and vector I/O through GDAL/OGR.

### Raster output

**`raster_output_format` enum:**
- `cog` — Cloud Optimised GeoTIFF (default; falls back to gtiff for small rasters)
- `gtiff` — plain compressed GeoTIFF
- `same_as_input` — use the GDAL driver of the input raster

**Key raster methods:**
- `write_full(data, meta, path, format=cog)`
- `finalise_output(path)`
- `create_output_raster(path, meta, format=cog, n_bands=1)`

### Vector output

Various file formats are available for vector output but preference is for either csv
or FlatGeobuf. FlatGeobuf should be the default output since it is more storage efficient
than csv, includes georeferencing, and does not have the same file writing issues
as gpkg. GPKG struggles wien mapping across drives, which coould be a common problem if 
users are on a WSL2 implementation of linux writing to `/mnt/c`.  

Supported vector outputs are:
- GeoPackage (`.gpkg`)
- FlatGeobuf (`.fgb`)
- CSV

Vector features use:
- WKT geometry
- typed attribute fields
- GDAL/OGR layer creation and feature writing

Typical vector-writing methods include:
- `write_gpkg(...)`
- `write_flatgeobuf(...)`
- `write_csv(...)`
- `read_layer(...)`

GeoPackage and FlatGeobuf are written through GDAL/OGR; no direct SQLite coding
is required in `lsdtt3`.

Fields in FlatGeobuf should be typed and not strings.


---


## Vector I/O

Vector I/O is handled through `lsdtt_data_loader` and shared helpers in
`src/core/lsdtt_vector_utils.hpp/.cpp`.

### Preferred output formats

See above, FaltGeobuf is the preferred format. 

Supported vector formats are:

- **GeoPackage (`.gpkg`)**  
  Default binary single-file GIS output.

- **FlatGeobuf (`.fgb`)**  
  Lightweight binary single-file GIS output.

- **CSV**  
  Mainly for backward-compatible tabular export and plotting workflows.

### Attribute typing

Vector fields should preserve their natural types:
- integer attributes → integer fields
- floating-point attributes → real fields
- text attributes → string fields

Do not coerce all vector attributes to strings unless there is a compelling
backward-compatibility reason.

This is important for fields such as:
- stream order
- junction id
- receiver id
- row / col / flat index
- elevation
- flow distance

### CSV spatial convention

If point data are written to CSV for interoperability with older scripts, the
CSV should contain:
- `latitude`
- `longitude`

in **EPSG:4326**, using GDAL/OGR coordinate transformation from the source CRS.

### Shared helper split

- **`src/utils/lsdtt_raster_geometry.hpp`**  
  header-only grid geometry helpers such as row/col to cell-centre x/y

- **`src/core/lsdtt_vector_utils.hpp/.cpp`**  
  GDAL/OGR-backed vector/projection helpers such as reusable coordinate
  transformations and simple WKT generation

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
or `lsdtt_flow_base.cpp`. Those are infrastructure.

---

## lsdtt_flow_types

File: `src/core/lsdtt_flow_types.hpp` (header-only, no `.cpp`)

Shared types included by both `lsdtt_flow_base.hpp` and `lsdtt_flowinfo_cpu.hpp`
to avoid circular includes.

**Include chain (acyclic):**
```
lsdtt_flow_types.hpp         (no lsdtt3 includes — just <cstdint> etc.)
       ↑
lsdtt_flow_base.hpp          (includes lsdtt_flow_types.hpp)
       ↑
lsdtt_flowinfo_cpu.hpp       (includes lsdtt_flow_types.hpp + lsdtt_flow_base.hpp)
       ↑
lsdtt_flow_base.cpp          (in lsdtt3_cpu — includes lsdtt_flowinfo_cpu.hpp, no cycle)
```

**Contents:**

- `k_nodata_receiver = -1`, `k_outlet_receiver = -2`, `k_no_direction = 255` — all `std::int64_t` / `uint8_t`
- `lsdtt_flat_index_slice` — zero-copy view into s_vector; `value_type = std::int64_t`
- `lsdtt_flowinfo_options` — construction options shared by flow routing objects

All cell indices use `std::int64_t` directly throughout the flow routing layer.
There is no `flat_index_type` typedef — it was removed to simplify the code and
avoid explaining a flexible type that was rarely changed in practice.

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

Abstract base class for D8 flow routing backends.

**Factory function `make_flow_base()`**:
Declared in `lsdtt_flow_base.hpp`, implemented in `lsdtt_flow_base.cpp` (lsdtt3_cpu).
Currently returns `std::make_unique<lsdtt_flowinfo_cpu>`.

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
1. `validate_inputs()`
2. `parse_boundary_conditions()`
3. `initialise_storage()`
4. `build_receiver_graph_d8()` — parallelisable with OMP collapse(2)
5. `build_donor_count_d8()` — serial
6. `build_outlets()` — serial
7. `build_topological_order()` — Kahn's BFS, serial
8. `build_s_vector()` — iterative DFS, serial

**s_vector contract:**
- `s_vector_[s_vector_start_[i]] == i` — node is always first in its own slice
- `upslope_nodes(flat)` — returns zero-copy `lsdtt_flat_index_slice`
- `is_upslope_of(candidate, target)` — O(1) interval membership test

---

## lsdtt_flowgraph_cpu

File: `src/cpu/lsdtt_flowgraph_cpu.hpp` (stub — implementation pending)

DAG-based flow routing for D-infinity and MFD. No s_vector.

**Data model (CSR):**
- `receiver_offsets_`, `receiver_indices_`, `receiver_weights_`
- Cell `i`'s receivers: `[receiver_offsets_[i], receiver_offsets_[i+1])`

**Routing methods:** `d_infinity` (Tarboton 1997), `mfd` (Quinn 1991)

---


## Junction-network architecture

The junction-network layer is the compact channel-graph analogue of the pixel
flow-routing graph.

### Files

- `src/core/lsdtt_junction_network_types.hpp`
- `src/core/lsdtt_junction_network_base.hpp`
- `src/cpu/lsdtt_junction_network_cpu.hpp/.cpp`

### Type policy

- **pixel / raster flat indices** use `std::int64_t`
- **junction indices** use `std::int32_t`

This split is deliberate:
- full rasters may exceed the safe range of 32-bit indexing
- junction counts are much smaller than pixel counts

### Construction model

`lsdtt_junction_network_cpu` is built on top of an existing flow-routing object.

The flow object is not owned by the junction network and must outlive it.

### Persistent storage

The persistent raster-scale storage is intentionally small:
- channel mask
- junction-id raster

The compact graph stores:
- junction flat index
- receiver junction
- donor count
- stream order
- outlets
- topological order
- junction `s_vector`

Temporary per-pixel products, such as pixel stream order used during
construction or export, should remain local and not be stored permanently unless
there is a clear reason.

---


## CPU analysis free functions (src/cpu/)

All take `const lsdtt_raster_cpu&` directly — not `lsdtt_raster_base`.
All return `lsdtt_raster_cpu`. Input is never modified.

**`lsdtt_raster_algebra` (header-only, `src/cpu/lsdtt_raster_algebra.hpp`):**

Lives in `namespace lsdtt::algebra`. All operations are element-wise; output inherits
metadata from the left-hand operand. A cell is nodata if either input cell is nodata.
Compatibility (same n_rows, n_cols, cell_size, x_min, y_min, no_data_value) is
enforced by `assert_compatible()` — throws `std::invalid_argument` on mismatch.
Inner loops are OpenMP-parallelised.

- `assert_compatible(a, b)` — throws if rasters are not element-wise compatible
- `raster_add(a, b)` — result = a + b
- `raster_subtract(a, b)` — result = a − b
- `raster_multiply(a, b)` — result = a × b
- `raster_divide(a, b)` — result = a / b; additionally nodata where b == 0
- `raster_log(a)` — natural log; additionally nodata where a ≤ 0
- `raster_log10(a)` — base-10 log; additionally nodata where a ≤ 0

**`lsdtt_raster_surface_metrics_cpu`:**
- `compute_hillshade(raster, azimuth=315, altitude=45, z_factor=1.0)`
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
- `identify_interior_nodata(raster)` → `lsdtt_array2d<uint8_t>`
- `fill_interior_nodata(raster, method=laplacian, max_iter=10000, tolerance=1e-6, n_threads=0)`
- `clip_to_maximum_rectangle(raster)`
- `nodata_fill_method` enum: `iterative_weighted_average`, `laplacian`

**`lsdtt_raster_spectral_cpu`:**
- `wavelet_quadratic_fit(raster, lambda=100.0, n_threads=0)`
  — FFTW double API, Gaussian kernel
  — Returns `poly_fit_coeffs{a, b, c, d, e, f}` — six `lsdtt_raster_cpu` members

**`lsdtt_raster_synthetic_surfaces_cpu`:**
- Takes `lsdtt_raster_meta` or raw dimensions
- Generators: uniform random, red noise, blue noise, diamond-square, periodic wave,
  global quadratic, random quadratic

---

## lsdtt-raster-preprocessing CLI

File: `cli/lsdtt-raster-preprocessing.cpp`  |  Version: 0.2.0

**Preprocessing pipeline (in order):**
1. **Projection** — detect and reproject geographic / Web Mercator rasters
2. **Threshold masking** — remove pixels outside elevation range
3. **Interior nodata filling** — fill enclosed voids
4. **Clip to rectangle** — largest valid rectangle
5. **Depression breaching** — Lindsay 2016 or Lindsay 2020 least-cost
6. **Pit filling** — Barnes 2014 Priority-Flood
7. **Hillshade** — Zevenbergen & Thorne 1987

Steps 2–7 always operate on the (possibly reprojected) raster. If reprojection
occurs the new file is reopened before step 2.

**Threshold masking semantics:**
`minimum_elevation` masks pixels **at or below** the threshold (i.e. `value <= min_elev`).
Default is `0.0`, which catches ocean pixels that raster providers commonly set to exactly 0.
`maximum_elevation` masks pixels above the threshold (`value > max_elev`); default `30000.0`.

**Projection parameters:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enforce_projection` | bool | true | Reproject geographic/WebMercator input automatically. If false, exit with error. |
| `resampling_method` | string | cubic | Any GDAL warp resampling string. |
| `grid_spacing` | float | -1 | Target grid spacing (m). ≤ 0 = auto-snap from input pixel size. |
| `allow_lambert_projection` | bool | false | Use LAEA instead of UTM when raster spans ≥ 3 UTM zones. |

**All other parameters:**
```
read_fname, write_fname, read_path, write_path
write_output_format: cog / gtiff / same_as_input
write_hillshade_raster, write_fill_raster, write_fill_diff_raster
write_breached_raster, write_carve_diff_raster
write_filled_nodata_raster, write_clipped_rectangle_raster
breach_before_fill, fill_remainder_after_breach
fill_algorithm_choice: barnes2014 / wang_liu
breach_algorithm_choice: lindsay2016 / least_cost
nodata_fill_method_choice: laplacian / iterative_weighted_average
hillshade_azimuth, hillshade_altitude, hillshade_z_factor
fill_min_slope, minimum_elevation, maximum_elevation
breach_lc_max_dist, breach_lc_max_cost, breach_lc_flat_increment
n_threads, nodata_fill_max_iter
```

**Diff raster outputs** (use `lsdtt::algebra::raster_subtract` from `lsdtt_raster_algebra.hpp`):

| Parameter | Output suffix | Formula | Interpretation |
|---|---|---|---|
| `write_fill_diff_raster` | `_fill_diff.tif` | filled − original | Positive = surface raised by filling |
| `write_carve_diff_raster` | `_carve_diff.tif` | breached − original | Negative = surface lowered by breaching |

**Help / parameter printing:**
- `-h` flag prints usage to stdout and exits immediately (before any processing).
- `print_parameters()` is called at the end of every successful run to echo resolved parameters.
- `print_help()` (writes `_README` file) is **not** called during normal runs; it is only triggered by `-h`.

---

## lsdtt-test-flowinfo CLI

File: `cli/lsdtt-test-flowinfo.cpp`

Builds D8 flow routing from a DEM and writes diagnostic rasters. Uses
`make_flow_base()` — compatible with future GPU backend without modification.

**Includes** `lsdtt_flow_base.hpp`, NOT `lsdtt_flowinfo_cpu.hpp`.

**Construction pattern:**
```cpp
std::unique_ptr<lsdtt::lsdtt_flow_base> flow =
    lsdtt::make_flow_base(dem_data, dem_mask, meta, fi_opts);
flow->n_valid_cells();
flow->upslope_nodes(flat);
```

**Key parameters:**
```
read_fname, write_fname, read_path, write_path
boundary_conditions: oooo
breach_before_fill, fill_remainder_after_breach, fill_before_routing
fill_algorithm_choice, breach_algorithm_choice
fill_min_slope, minimum_elevation, maximum_elevation
write_n_contributing_pixels_raster, write_flow_direction_raster
write_donor_count_raster, write_outlet_raster
write_fill_raster, write_breached_raster
write_output_format: cog
n_threads: 0
```

**Flow direction codes:**
`0`=N `1`=NE `2`=E `3`=SE `4`=S `5`=SW `6`=W `7`=NW; outlets = `-2`

---

## Known compile errors fixed (for reference)

| File | Issue | Fix |
|------|-------|-----|
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
11. `src/cpu/` = all CPU implementations (raster, flow, and algebra)
12. `src/gpu/` = all GPU implementations
13. `src/core/lsdtt_crs_suggest` inspects CRS; `src/core/lsdtt_raster_reproject`
    acts on it — these are deliberately separate so suggestion logic can be used
    without triggering I/O
14. Reprojected rasters always receive a projection suffix (_utm / _polar / _lambert)
    so the provenance of every output file is unambiguous from its filename
15. `lsdtt_raster_algebra.hpp` is header-only (no `.cpp`); include from `src/cpu/`.
    Use `lsdtt::algebra::raster_subtract(a, b)` etc. for diff rasters in CLIs.
16. Every CLI uses `lsdtt_param_parser`
17. CLI parameters are grouped by theme, not declared in one unstructured block
18. In CLIs, `help_map` entries come immediately after the corresponding
    parameter default declaration
19. CLI parameter variables are intentionally not `const` because
    cross-parameter logic may adjust them after parsing
20. Pure grid-geometry helpers belong in `src/utils/`; GDAL/OGR-backed vector
    and projection helpers belong in `src/core/`
21. Vector outputs should prefer typed attributes and single-file binary formats
    such as FlatGeobuf

