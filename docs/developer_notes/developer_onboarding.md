# lsdtt3 Developer Onboarding

This document is for developers porting analysis code from earlier versions of
LSDTopoTools (LSDRaster, LSDFlowInfo, etc.) or writing new analyses against the
lsdtt3 architecture.

!!! note "Pull requests are from collaborators only, for now"

    Until the code base is stable, as decided by Simon M. Mudd (SMM), only
    repository collaborators can open pull requests. Anyone can open an
    [issue](https://github.com/LSDtopotools/lsdtt3/issues) to report a bug or
    ask for a feature, or fork the code for their own use. See `CONTRIBUTING.md`
    in the repository.

IMPORTANT NOTE: The codebase has become large, with interdependencies. 
We have aimed to cover a reasonable amount of this code with ctest, 
but if you try to modify the code without the assistance of a (good) LLM there is a very high chance your will break something.
If you want to modify things you should feed your LLM both `lsdtt3_architecture.md` and `cli_parameter_conventions.md` at a minimum.  

---

## Code style

All names use `snake_case` throughout — classes, functions, variables, parameters,
and filenames. No camelCase anywhere. All code lives in `namespace lsdtt`. No
`using namespace std` anywhere.

Opening braces always go on their own line:

```cpp
if (condition)
{
    do_something();
}
```

And we use a 4 space tab for indentation. 

Code text should be ASCII-only in source files, build scripts, and machine-read
config files. Replace Unicode punctuation/symbols with ASCII equivalents (for
example `-`, `->`, `<=`, `>=`, `x`).

Output raster parameters follow the convention `write_*_raster`. Algorithm
selection parameters follow `*_algorithm_choice`. These conventions exist so
parameter files are self-documenting and consistent across CLIs.

---
## CLI design rule

Every CLI uses `lsdtt_param_parser`. The full required layout is the
"CLI design rule" in [`lsdtt3_architecture.md`](lsdtt3_architecture.md#cli-design-rule).
In outline:

1. file header, includes, and `k_tool_name` / `k_status` / `k_citation`
2. `print_banner()`, skipped when `lsdtt_param_parser::help_json_requested()`
   is true so `-h json` prints only JSON
3. `lsdtt_param_parser::parse_args_deferred()`, which handles `-v`, `-c`, no
   arguments and bad arguments, and records `-h`
4. parameter defaults grouped by theme, each followed immediately by its
   `help_map` entry
5. if help was requested, `print_terminal_help()` and return. Help is rendered
   from the declarations, so no parameter file is needed. Then create the run
   record, `lsdtt::lsdtt_run_record run_record(k_tool_name, argc, argv);`
6. parameter parsing and cross-parameter logic
7. input loading, processing pipeline, output writing
8. `print_help()` and `print_parameters()` for reproducibility, then
   `return run_record.succeed();` on every successful path

Do not add command-line flags. Every option goes in the parameter file. The
`-h` output must be enough for someone with only the binary to write a
parameter file.

### Adding a new CLI

A new `cli/lsdtt-foo.cpp` must be registered in up to three places:

1. `lsdtt_cli(lsdtt-foo cli/lsdtt-foo.cpp)` in the top-level `CMakeLists.txt`.
   This is required.
2. The pip wheel, which declares its CLIs in four places that must agree:
   `lsdtt_python_cli(lsdtt-foo "cli/lsdtt-foo.cpp")` in
   `python_bindings/CMakeLists.txt`; `lsdtt-foo = "lsdtt3._cli:main"` under
   `[project.scripts]` in `python_bindings/pyproject.toml`; the `_DRIVERS` set
   in `python_bindings/src/lsdtt3/_cli.py`; and `console_scripts` in
   `python_bindings/api/public_api.json`, which you regenerate with
   `python python_bindings/tools/generate_public_api_manifest.py --update`.
   The CTest case `unit_wheel_cli_inventory` fails if any of the four differs
   from the root `lsdtt_cli()` targets. So does the pytest file
   `python_bindings/tests/test_cli_packaging.py`, which runs only in the bindings
   suite.
3. A profile in `tools/cli_profiles.txt`, so that it ships in the release
   bundles. The `core` profile requires at least one CTest test that runs the
   CLI. `qgis` is for CLIs wired into the QGIS plugin.

Then re-run `cmake`, because `unit_cli_parameter_conventions` refuses to run if
any `cli/` source has no binary. Also add its reference page under
`docs/cli_reference/`, with the generated-parameters markers, and run
`tools/generate_cli_reference.py --update`.

### Parameter declaration rule

In CLI files:
- parameters should **not** be declared `const`
- parameters should be grouped by workflow theme (where possible)
- every `help_map` entry should appear **immediately after** the corresponding
  parameter default declaration

This pattern is now the standard for all CLIs because it makes developer-facing
options easy to inspect and keeps the parameters readable in one place.

### Parameter naming conventions

Every parameter you add must follow the naming standard in
`docs/developer_notes/cli_parameter_conventions.md`. The standard ensures that parameter names are
self-documenting across all CLIs.

**Checklist when adding a parameter:**

1. **Look for a shared name first.** Is this a threshold, an output format, a
   raster file, or a window size? Check the conventions table to see if it already
   has a canonical name. If so, use it exactly. Do not invent a variant.

2. **Use the shared declaration layer if available.** If your parameter is part of
   a semantic group (flow routing, output formats, basin outlets, fitting parameters,
   synthetic domains), check whether your CLI can call the corresponding
   `declare_*_parameters()` function from `src/core/lsdtt_common_parameters.hpp`.
   If yes, use it — do not hand-declare the parameter.

3. **Follow the unit-suffix rule.** Always spell out units: `_metres`, `_pixels`,
   `_degrees`, never `_m`, `_px`, `_deg`. A threshold on area is
   `<what>_area_threshold_m2`, never `_area_threshold` or `_area` alone.

4. **Use the right prefix.** Output switches start with `write_` (not `print_` for
   files), stage switches start with `run_`, boolean rasters end in `_raster` or
   `_rasters`. Bounds use `min_`/`max_` unless the full word `minimum` reads more
   clearly (e.g. `basin_minimum_contributing_pixels`).

5. **Distinguish the two fill concepts.** If your CLI does filling, know which one
   you need: `raster_is_filled` (input assertion, compute-saving) or
   `fill_before_routing` (analytical choice). Write help text that says which
   question the parameter answers.

6. **Document the parameter's role in the help entry.** Every parameter must have
   a short description (one line, including units) and a longer explanation
   (covering valid range, cross-references to related parameters, caveats). Do not
   skip the units — a developer will find this in the CSV reference, and they may
   not have the source code to hand.

7. **Plan for migration.** Even if you are not renaming anything today, know that
   if your parameter is later renamed, the old name will be added to an alias map
   (see `src/core/lsdtt_param_parser.hpp`). Existing `.param` files will keep
   working with a deprecation warning. This means naming is not a permanent choice —
   you can fix a bad name later without breaking users. A CLI that has never been
   in a release tag can simply be renamed.

8. **Build and check.** Rebuild, then run
   `ctest -R "unit_cli_parameter_conventions|unit_cli_reference_is_current"`.
   The first test fails on any naming violation; there is no allow-list to add
   to. The second fails until you regenerate the reference pages with
   `python3 tools/generate_cli_reference.py --update --bin-dir build`. Also update
   any hand-written examples in `docs/cli_reference/` that use the parameter.

### File I/O parameter note

For CLIs using `lsdtt_param_parser`:

- `read_fname` may include an extension (for example `my_dem.tif`)
- `write_fname` should remain a stem with no extension

If `read_fname` includes an extension, the parser splits it internally into
stem + extension before reconstructing the full input path.


---

## Where to put new code

### The short version

If you are porting or writing analysis code, **work entirely in `src/cpu/`**.
You do not need to touch `src/core/` at all.

`src/core/` contains the abstract interfaces (`lsdtt_raster_base`,
`lsdtt_flow_base`), I/O (`lsdtt_data_loader`), metadata (`lsdtt_raster_meta`),
the parameter parser, and CRS/reprojection utilities. These are infrastructure.
New analysis code does not belong there and should never need to modify them.

`src/cpu/` is where all analysis lives. It contains the concrete CPU raster
class (`lsdtt_raster_cpu`), the D8 flow routing object (`lsdtt_flowinfo`), and
all analysis grouped into thematic free-function files (this is a partial list):

```
src/cpu/lsdtt_raster_surface_metrics_cpu.hpp/.cpp   — hillshade, moving window
src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp/.cpp — fill, breach, masking
src/cpu/lsdtt_raster_nodata_filling_cpu.hpp/.cpp     — void filling, clipping
src/cpu/lsdtt_raster_algebra.hpp/.cpp               — element-wise arithmetic
src/cpu/lsdtt_raster_spectral_cpu.hpp/.cpp          — spectral / wavelet
src/cpu/lsdtt_raster_synthetic_surfaces_cpu.hpp/.cpp — test surface generators
```

To add a new group of analysis functions — say, glacial erosion — you create
`src/cpu/lsdtt_glacial_erosion_cpu.hpp/.cpp`, write your free functions there, and
add the `.cpp` to `LSDTT3_CPU_SOURCES` in `CMakeLists.txt`. That is the entire
change required outside your new files.

### Free functions vs member functions

In LSDTopoTools2 the objects (especially LSDRaster) contained vast numbers of member function which made finding analyses a challenge.

Here default choice for new analysis is a **free function** in a thematic
`.hpp/.cpp` pair, taking `const lsdtt_raster_cpu&` or `const lsdtt_flowinfo&`, or any new class
as arguments.

Add a **member function** to `lsdtt_raster_cpu` or `lsdtt_flowinfo` only when
the function needs direct access to private data members — i.e. data that is not
reachable through the existing public accessors. For almost all analysis code
this will not be the case.

The practical test: if you can write the function using only `get_data()`,
`get_mask()`, `get_meta()`, `upslope_nodes()`, `receiver_of_unchecked()` and
similar public methods, write it as a free function. If you find yourself needing
`receiver_[flat]` or `s_vector_start_[flat]` directly, it can go on the class.

---

## Example: adding a member function to lsdtt_raster_cpu

Note: you really only want to do this in extreme cases, preference is for free functions!

### 1. Declare it in `src/cpu/lsdtt_raster_cpu.hpp`

Add the declaration in the public section:

```cpp
// Returns a new raster with val added to every valid pixel.
// Nodata cells remain nodata.
[[nodiscard]] lsdtt_raster_cpu add_constant(float val) const;
```

### 2. Implement it in `src/cpu/lsdtt_raster_cpu.cpp`

```cpp
lsdtt_raster_cpu lsdtt_raster_cpu::add_constant(float val) const
{
    const int   n      = meta_.n_rows * meta_.n_cols;
    const float nodata = static_cast<float>(meta_.no_data_value);

    std::vector<float>   out_data(n);
    std::vector<uint8_t> out_mask(n);

    const float* src_data = data_.data();
    const uint8_t* src_mask = mask_.data();

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        if (src_mask[i])
        {
            out_data[i] = src_data[i] + val;
            out_mask[i] = 1u;
        }
        else
        {
            out_data[i] = nodata;
            out_mask[i] = 0u;
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta_);
}
```

No changes to `CMakeLists.txt` are needed because `lsdtt_raster_cpu.cpp` is
already in `LSDTT3_CPU_SOURCES`.

### 3. Call it from a CLI

In a CLI source file (e.g. `cli/lsdtt-raster-preprocessing.cpp`):

```cpp
#include "../src/cpu/lsdtt_raster_cpu.hpp"

// dem is an lsdtt_raster_cpu already loaded from disk
lsdtt_raster_cpu raised = dem.add_constant(1.0f);

std::string o = out_path(pp, "_raised.tif");
loader.write_full(raised.get_data(), raised.get_meta(), o, out_fmt);
```

If instead this were a thematic free function (e.g. a hillshade variant living
in `lsdtt_raster_surface_metrics_cpu`), the call looks identical from the CLI's
perspective — the only difference is the include path and that the raster is
passed as an argument:

```cpp
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"

lsdtt_raster_cpu hs = lsdtt::compute_hillshade(dem, 315.0f, 45.0f, 1.0f);
```

---

## Notes on types

### `[[nodiscard]]`

You will see `[[nodiscard]]` on many functions throughout the headers, especially
query functions on `lsdtt_flowinfo`:

```cpp
[[nodiscard]] std::int64_t receiver_of_unchecked(std::int64_t flat) const noexcept;
[[nodiscard]] lsdtt_flat_index_slice upslope_nodes(std::int64_t flat) const;
[[nodiscard]] lsdtt_raster_cpu add_constant(float val) const;
```

This attribute tells the compiler to emit a warning if the return value is
discarded by the caller. It is used on query and analysis functions, that is functions
with no side effects, because discarding their return value is almost certainly
a mistake. For example:

```cpp
fi.upslope_nodes(flat);      // compiler warns — result thrown away, probably a bug
dem.add_constant(1.0f);      // compiler warns — new raster thrown away, definitely a bug
```

Without `[[nodiscard]]` these would silently compile and do nothing useful (and potentially something bad).
It is not used on functions that act through side effects (e.g. `clear()`,
`print_raster_info()`), where a void return or ignored return is normal.

### Cell indices and sentinel values

For a handful of data members in `lsdtt_flowinfo` we need to be able to hold very large integer numbers. 
For example a 50000x500000 DEM will have more pixels than there are integers available using a 32 bit integer. So we use `std::int64_t` for these. 

All flow routing code uses `std::int64_t` for cell indices throughout. There is
no `flat_index_type` typedef — using the concrete type directly makes the code
easier to read and avoids explaining a flexible typedef that was rarely changed
in practice. `int64_t` is safe for any raster you can load into memory.

`src/core/lsdtt_types.hpp` defines the sentinel values you will encounter when
iterating the flow graph. That file is the single home for shared types — see
design principle 13 in the architecture note; do not add a new `*_types.hpp`:

```cpp
k_nodata_receiver = -1   // cell is nodata / outside the valid mask
k_outlet_receiver = -2   // cell drains off-grid or to a boundary outlet
k_no_direction    = 255  // direction code for nodata cells
```

Always test against these constants rather than the raw values:

```cpp
if (fi.receiver_of_unchecked(flat) == lsdtt::k_outlet_receiver)
{
    // this cell is a basin outlet
}
```

Donor counts are stored as `std::uint8_t`: a cell on a regular grid can have
at most 8 direct donors under D8 routing, so a byte is sufficient and keeps
the donor-count array compact.

---

## Adding tests

Every test must say how it knows its answer. Give each new CTest test a row in
`tests/cmake/evidence_labels.cmake`, unless the CLI harness or the name rules
there already label it. Give each new Python test pytest markers. Either way, a
test needs one or more evidence labels (`exact_solution`, `property`,
`metamorphic`, `differential`, `characterisation`, `smoke`, `negative`,
`conformance`, `performance`) and exactly one level (`unit` or `system`).
Configure and pytest collection both stop on a test without them. The
categories are defined in the "Testing categories" section of
[`lsdtt3_architecture.md`](lsdtt3_architecture.md). Prefer the strongest
evidence the test can give: an exact solution, a property or a metamorphic
relation, rather than a stored snapshot.

---

## Summary checklist for new analysis code

- Work in `src/cpu/` 
- Add tests labelled by evidence and level (see "Adding tests" above)
- Write free functions unless you need private data members
- Add new `.cpp` files to `LSDTT3_CPU_SOURCES` in `CMakeLists.txt`
- Use `std::int64_t` for all cell indices, and the sentinel constants for
  boundary/nodata tests
- Mark query functions `[[nodiscard]]`
- No `using namespace std`, no camelCase, braces on their own line
