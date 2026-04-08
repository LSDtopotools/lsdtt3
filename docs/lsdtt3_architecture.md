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
├── apptainer_def/
│   └── apptainer_setup.md
├── cli/
│   ├── lsdtt-channel-extraction.cpp
│   ├── lsdtt-feature-finder.cpp
│   ├── lsdtt-raster-basin-stats.cpp
│   ├── lsdtt-raster-preprocessing.cpp
│   ├── lsdtt-raster-synthetic-tests.cpp
│   ├── lsdtt-surface-metrics.cpp
│   └── lsdtt-test-flowinfo.cpp
├── docs/
│   ├── basin_extraction_guide.md
│   ├── debian_packaging.md
│   ├── developer_onboarding.md
│   ├── lsdtt3_architecture.md
│   ├── lsdtt_utils_guide.md
│   └── using_lsdtt_flowinfo.md
├── python_helpers/
│   ├── dem_downloader_instructions.md
│   └── download_dem_tiles.py
├── src/
│   ├── core/
│   │   ├── lsdtt_types.hpp                  — all lightweight shared primitive types (see below)
│   │   ├── lsdtt_raster_meta.hpp            — metadata struct + CRS validation
│   │   ├── lsdtt_raster_base.hpp/.cpp       — abstract raster dispatch interface
│   │   ├── lsdtt_data_loader.hpp/.cpp       — GDAL/OGR I/O, RAII, COG output
│   │   ├── lsdtt_backend_factory.hpp/.cpp   — make_raster_engine(), gpu_available()
│   │   ├── lsdtt_param_parser.hpp/.cpp      — parameter file parser
│   │   ├── lsdtt_basin_base.hpp             — abstract basin interface
│   │   ├── lsdtt_junction_network_base.hpp  — abstract channel network interface + lsdtt_junction_network_options
│   │   ├── lsdtt_flow_base.hpp              — abstract D8 flow interface + lsdtt_flowinfo_options
│   │   ├── lsdtt_crs_suggest.hpp/.cpp       — CRS inspection and projection recommendation
│   │   ├── lsdtt_raster_reproject.hpp/.cpp  — CRS reprojection
│   │   ├── lsdtt_raster_warp.hpp/.cpp       — GDAL warp helpers
│   │   ├── lsdtt_point_input.hpp/.cpp       — point/coordinate input handling
│   │   └── lsdtt_vector_utils.hpp/.cpp      — GDAL/OGR-backed vector helpers
│   ├── cpu/
│   │   ├── lsdtt_raster_cpu.hpp/.cpp                       — concrete CPU raster class
│   │   ├── lsdtt_flow_base.cpp                             — factory impl for make_flow_base()
│   │   ├── lsdtt_flowinfo_cpu.hpp/.cpp                     — D8 CPU flow routing
│   │   ├── lsdtt_flowgraph_cpu.hpp                         — D-inf/MFD DAG routing (stub)
│   │   ├── lsdtt_depression_hierarchy_cpu.hpp/.cpp         — depression hierarchy construction
│   │   ├── lsdtt_fill_spill_merge_cpu.hpp/.cpp             — fill-spill-merge depression routing
│   │   ├── lsdtt_flow_snapping_cpu.hpp/.cpp                — snapping points to channel network
│   │   ├── lsdtt_junction_network_cpu.hpp/.cpp             — CPU channel junction network
│   │   ├── lsdtt_basin_cpu.hpp/.cpp                        — CPU basin objects
│   │   ├── lsdtt_junction_basin_cpu.hpp/.cpp               — junction-based basin extraction
│   │   ├── lsdtt_junction_basin_selection_cpu.hpp/.cpp     — basin selection helpers
│   │   ├── lsdtt_raster_algebra.hpp/.cpp                   — element-wise raster algebra
│   │   ├── lsdtt_raster_surface_metrics_cpu.hpp/.cpp       — hillshade, moving window, polyfit metrics
│   │   ├── lsdtt_raster_surface_fitting_cpu.hpp/.cpp       — plane and polynomial surface fitting to masks
│   │   ├── lsdtt_raster_flow_preprocessing_cpu.hpp/.cpp    — fill, breach, nodata
│   │   ├── lsdtt_raster_nodata_filling_cpu.hpp/.cpp        — interior void filling
│   │   ├── lsdtt_raster_spectral_cpu.hpp/.cpp              — wavelet quadratic fit; home for future spectral surface fitting
│   │   ├── lsdtt_raster_opencv_cpu.hpp/.cpp                — OpenCV-based raster processing
│   │   ├── lsdtt_raster_feature_finder_cpu.hpp/.cpp        — feature-specific raster workflows
│   │   └── lsdtt_raster_synthetic_surfaces_cpu.hpp/.cpp    — test surface generators including mountain range geometries
│   └── utils/
│       ├── lsdtt_array2d.hpp              — header-only 2D array
│       ├── lsdtt_raster_geometry.hpp      — header-only raster geometry helpers
│       ├── lsdtt_statistics.hpp/.cpp      — descriptive statistics (median, MAD, stddev, etc.)
│       ├── lsdtt_histogram.hpp/.cpp       — histogram, PDF, CDF
│       ├── lsdtt_statistical_tests.hpp/.cpp — nonparametric tests, KS, circular statistics
│       ├── lsdtt_splines.hpp/.cpp         — 1D spline interpolation
│       └── lsdtt_regression.hpp/.cpp      — linear and polynomial regression helpers
├── CMakeLists.txt
├── diagnose_fftw.sh
└── README.md
```

---

## CMakeLists.txt summary

Three main static libraries:

- **`lsdtt3_core`** — abstract interfaces, factories, GDAL/OGR I/O, shared types
- **`lsdtt3_utils`** — generic mathematical / statistical utilities on valid flat data
- **`lsdtt3_cpu`** — concrete CPU implementations, raster operations, flow routing

`lsdtt3_cpu` depends on `lsdtt3_core` and `lsdtt3_utils`.
All CLI executables link against `lsdtt3_cpu` and therefore also gain access to
the utility-layer numerical functions.



**`lsdtt3_core` sources** (`src/core/`):
- `lsdtt_raster_base.cpp`
- `lsdtt_data_loader.cpp`
- `lsdtt_backend_factory.cpp`
- `lsdtt_param_parser.cpp`
- `lsdtt_crs_suggest.cpp`
- `lsdtt_raster_reproject.cpp`
- `lsdtt_raster_warp.cpp`
- `lsdtt_vector_utils.cpp`
- `lsdtt_point_input.cpp`

**`lsdtt3_cpu` sources** (`src/cpu/`):
- `lsdtt_raster_cpu.cpp`
- `lsdtt_raster_algebra.cpp`
- `lsdtt_raster_surface_metrics_cpu.cpp`
- `lsdtt_raster_surface_fitting_cpu.cpp`
- `lsdtt_raster_flow_preprocessing_cpu.cpp`
- `lsdtt_raster_nodata_filling_cpu.cpp`
- `lsdtt_raster_spectral_cpu.cpp`
- `lsdtt_raster_synthetic_surfaces_cpu.cpp`
- `lsdtt_raster_feature_finder_cpu.cpp`
- `lsdtt_flowinfo_cpu.cpp`
- `lsdtt_flow_base.cpp`
- `lsdtt_junction_network_cpu.cpp`
- `lsdtt_depression_hierarchy_cpu.cpp`
- `lsdtt_fill_spill_merge_cpu.cpp`
- `lsdtt_flow_snapping_cpu.cpp`
- `lsdtt_basin_cpu.cpp`
- `lsdtt_junction_basin_cpu.cpp`
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

**`docs/lsdtt_utils_guide.md`**  
Developer guide to the `src/utils/` layer. Explains what belongs in utils,
what should remain in `core` / `cpu` / `gpu`, the valid-data contract for
utility functions, and when to use the main statistical tests and distance
measures.


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

- **`lsdtt-raster-basin-stats.cpp`**  
  Extracts basins from a junction network and computes raster statistics per basin.

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

lsdtt3 holds vector layers in memory using an `lsdtt_vector_layer` structure, which contains:
- an `lsdtt_vector_schema` structure to hold attribute field names and types
- a number of features, each stored as an `lsdtt_vector_record`. Records contain an `lsdtt_geom` geometry (which can be point, line or polygon) and a vector of attribute table values for that feature.

The `loader` object in `lsdtt_data_loader.hpp/cpp` manages vector writes, using:
- `begin_vector_write(...)` to start a writing session, takes parameters for paths, layer names, the schema, the crs, the format (gpkg/flatgeobuf) and whether or not to write a spatial index. Including the spatial index makes loading data and spatial queries faster in a GIS, but to calculate the (optional) spatial index, the entire vector layer needs to sit in memory (no feature-by-feature writing) which results in large overheads when writing large vector datasets.
- `write_point/line/polygon(..)` write individual features, which are streamed if no spatial index required and buffered in memory if spatial index required.
- `finish_vector_write()` closes the writing session, and calculates the spatial index if required.

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

- **`src/core/lsdtt_types.hpp`**  
  all lightweight shared primitive types (spatial, raster geometry, flow, junction)

- **`src/utils/lsdtt_raster_geometry.hpp`**  
  header-only grid geometry helpers such as row/col to cell-centre x/y;
  primitive types (lsdtt_xy, lsdtt_extent, lsdtt_rowcol_window) come from lsdtt_types.hpp

- **`src/core/lsdtt_vector_utils.hpp/.cpp`**  
  GDAL/OGR-backed vector/projection helpers such as reusable coordinate
  transformations

---


## Flow routing — how to navigate the files

There are several flow-routing files. This table explains which one you need.

| You are... | Include this | Notes |
|---|---|---|
| Writing analysis code (chi, steepness, drainage area, etc.) | `lsdtt_flowinfo_cpu.hpp` | Only file you need. Construct `lsdtt_flowinfo`, call its methods directly. |
| Writing a CLI that builds a flow object | `lsdtt_flow_base.hpp` | Call `make_flow_base()`. Returns `unique_ptr<lsdtt_flow_base>`. GPU-compatible. |
| Writing D-inf / MFD distributed-flow analysis | `lsdtt_flowgraph_cpu.hpp` | Stub — implementation pending. |
| Adding a GPU flow routing backend | `lsdtt_flow_base.hpp` + new `lsdtt_flowinfo_gpu.hpp` | Implement `lsdtt_flow_base`, add branch in `src/cpu/lsdtt_flow_base.cpp`. |
| Understanding cell-index types | `lsdtt_types.hpp` | Flow-routing cell indices use `std::int64_t` directly throughout; there is no `flat_index_type` typedef. |

**Analysis code** never needs to know about `lsdtt_flow_base`, `lsdtt_flow_types`,
or `lsdtt_flow_base.cpp`. Those are infrastructure.

---

## lsdtt_types

File: `src/core/lsdtt_types.hpp` (header-only, no `.cpp`)

**Single home for all lightweight primitive types shared across lsdtt3.**

Dependencies: only `<cstdint>`, `<cstddef>`. No GDAL, no OGR, no other lsdtt3 header.

**Rule:** Do not create a new `*_types.hpp` file. If you need a new lightweight
shared type, add it here. A separate types file is only justified when including
the type everywhere would create an unwanted compile-time dependency — that
situation should be extremely rare in this codebase.

**Include chain (acyclic):**
```
lsdtt_types.hpp                          (no lsdtt3 includes — just <cstdint> etc.)
       ↑
lsdtt_flow_base.hpp                      (includes lsdtt_types.hpp; owns lsdtt_flowinfo_options)
       ↑
lsdtt_flowinfo_cpu.hpp                   (includes lsdtt_flow_base.hpp)
       ↑
lsdtt_flow_base.cpp                      (in lsdtt3_cpu — includes lsdtt_flowinfo_cpu.hpp, no cycle)

lsdtt_types.hpp
       ↑
lsdtt_junction_network_base.hpp          (includes lsdtt_types.hpp; owns lsdtt_junction_network_options)
       ↑
lsdtt_junction_network_cpu.hpp           (includes lsdtt_junction_network_base.hpp)

lsdtt_types.hpp
       ↑
lsdtt_raster_geometry.hpp               (includes lsdtt_types.hpp + lsdtt_raster_meta.hpp)
```

**Contents:**

*Spatial primitives*
- `lsdtt_pixel_rc` — integer row/col location; negative = unset/invalid
- `lsdtt_xy` — double x/y projected coordinate pair

*Raster geometry*
- `lsdtt_extent` — axis-aligned bounding box (lower-left origin)
- `lsdtt_rowcol_window` — half-open `[row_begin, row_end)` / `[col_begin, col_end)` window

*Flow-routing types*
- `k_nodata_receiver = -1`, `k_outlet_receiver = -2` — `std::int64_t` sentinels
- `k_no_direction = 255` — `std::uint8_t` sentinel
- `lsdtt_flat_index_slice` — zero-copy view into s_vector; `value_type = std::int64_t`

*Junction-network types*
- `k_nodata_junction = -1`, `k_outlet_junction = -2` — `std::int32_t` sentinels
- `lsdtt_junction_index_slice` — zero-copy view into junction s_vector; `value_type = std::int32_t`

**`lsdtt_flowinfo_options`** lives in `src/core/lsdtt_flow_base.hpp` — it is
domain-specific to flow routing and would only clutter the shared types file.

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

**`lsdtt_junction_network_options`** lives in `src/core/lsdtt_junction_network_base.hpp`.

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

- `src/core/lsdtt_types.hpp` — shared primitive types (sentinels, index slices)
- `src/core/lsdtt_junction_network_base.hpp` — abstract interface + `lsdtt_junction_network_options`
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

**`lsdtt_raster_algebra` (`src/cpu/lsdtt_raster_algebra.hpp/.cpp`):**

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
- `moving_window_mean(raster, radius_metres, shape=circular)`
- `moving_window_local_relief(raster, radius_metres, shape=circular)` — max minus min in window
- `moving_window_stddev(raster, radius_metres, shape=circular)` — population stddev
- `moving_window_median(raster, radius_metres, shape=circular)`
- `moving_window_mad(raster, radius_metres, shape=circular)` — median absolute deviation
- `moving_window_tpi(raster, radius_metres, shape=circular)` — Topographic Position Index (Weiss 2001)
- `calculate_madk2(raster, window_radius_metres, lag_pixels=1, n_threads=0)` → `madk2_result` — Trevisani et al. (2023) directional roughness
- `compute_polyfit_metrics(polycoeff, request)` → `polyfit_metric_result` — batch computation of curvature and derivative rasters from wavelet quadratic fit coefficients
- `window_shape` enum: `square`, `circular`
- `polyfit_metric` enum: slope, aspect, curvatures (mean, Gaussian, max, min, profile, plan, tangential, etc.), rotor

**`lsdtt_raster_surface_fitting_cpu`:**

Plane and polynomial surface fitting to masked or labelled raster regions.
Use for piedmont / fan fitting, terrace treads, floodplain surfaces.
For spectral (long-wavelength) surface fitting see `lsdtt_raster_spectral_cpu`.

- `fit_plane_to_mask(dem, mask)` → `lsdtt_plane_fit_result` — OLS plane, z = a*x + b*y + c
- `fit_plane_to_mask_ransac(dem, mask, inlier_tolerance, max_iter=1000, seed=1)` → `lsdtt_plane_fit_result` — robust RANSAC fit; tolerates incised channels and outliers
- `evaluate_fitted_plane(meta, coeffs, n_threads=0)` → `lsdtt_raster_cpu`
- `residual_from_fitted_plane(dem, coeffs, n_threads=0)` → `lsdtt_raster_cpu`
- `fit_planes_to_labels(dem, label_raster)` → `std::map<int, lsdtt_plane_fit_result>` — one plane per label
- `fit_poly2_to_mask(dem, mask)` → `lsdtt_poly2_fit_result` — OLS degree-2 polynomial, z = a*x² + b*y² + c*x*y + d*x + e*y + f
- `evaluate_fitted_poly2(meta, coeffs, n_threads=0)` → `lsdtt_raster_cpu`
- `residual_from_fitted_poly2(dem, coeffs, n_threads=0)` → `lsdtt_raster_cpu`
- `poly2_max_principal_curvature(coeffs)` → `double` — max Hessian eigenvalue for curvature screening

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
- Takes `lsdtt_raster_meta` or raw dimensions (n_rows, n_cols, cell_size)
- Basic generators: uniform random, red noise, blue noise, diamond-square, periodic wave,
  global quadratic, random quadratic
- Mountain range generators (return `lsdtt_mountain_range_surface{dem, front_mask}`):
  - `generate_linear_range_surface` — bilateral piedmont, ridge in DEM interior
  - `generate_halfspace_range_surface` — ridge at DEM edge (Himalaya half-space analogy)
  - `generate_circular_range_surface` — elliptical dome in DEM interior
  - `generate_arc_range_surface` — parabolic curved front at one edge
- `lsdtt_range_edge` enum: `north`, `south`, `east`, `west`

---


## Utility-layer design (`src/utils/`)

The `src/utils/` layer is for generic mathematical, numerical, and lightweight
helper code that is reusable across the project and is not tied to GDAL/OGR,
file I/O, or a specific backend implementation.

Typical contents include:
- lightweight containers
- metadata-only raster geometry helpers
- descriptive statistics
- histogram / PDF / CDF tools
- binning utilities
- nonparametric statistical tests
- circular-data comparison tools

### Valid-data contract

Utility functions in `src/utils/` assume that the values passed to them are
already valid.

This means `utils` code should **not** perform:
- nodata filtering
- raster flattening
- moving-window extraction
- basin extraction
- backend-specific mask handling

Those tasks belong in:
- `src/core/` when they are shared infrastructure or GDAL/OGR-aware helpers
- `src/cpu/` for concrete CPU-side extraction and traversal
- `src/gpu/` for concrete GPU-side extraction and traversal

In practice, raster, moving-window, and basin code should extract valid values
first and then pass the resulting flat sample into `src/utils/` routines.




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
- `print_help()` on `lsdtt_param_parser` writes CSV and HTML help files, but this CLI does not call it during `-h` handling or normal runs.

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
13. All lightweight shared primitive types live in `src/core/lsdtt_types.hpp`.
    Do not create a new `*_types.hpp` file — add to lsdtt_types.hpp instead.
    A separate types file is only justified when inclusion would create an
    unwanted compile-time dependency; that situation should be extremely rare.
14. `src/core/lsdtt_crs_suggest` inspects CRS; `src/core/lsdtt_raster_reproject`
    acts on it — these are deliberately separate so suggestion logic can be used
    without triggering I/O
15. Reprojected rasters always receive a projection suffix (_utm / _polar / _lambert)
    so the provenance of every output file is unambiguous from its filename
16. `lsdtt_raster_algebra` lives in `src/cpu/lsdtt_raster_algebra.hpp/.cpp`; include from `src/cpu/`.
    Use `lsdtt::algebra::raster_subtract(a, b)` etc. for diff rasters in CLIs.
17. Every CLI uses `lsdtt_param_parser`
18. CLI parameters are grouped by theme, not declared in one unstructured block
19. In CLIs, `help_map` entries come immediately after the corresponding
    parameter default declaration
20. CLI parameter variables are intentionally not `const` because
    cross-parameter logic may adjust them after parsing
21. Pure grid-geometry helpers belong in `src/utils/`; GDAL/OGR-backed vector
    and projection helpers belong in `src/core/`
22. Vector outputs should prefer typed attributes and single-file binary formats
    such as FlatGeobuf

