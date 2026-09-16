# Python Bindings Architecture And API Contracts

The `lsdtt3` Python package exposes native terrain-analysis operations as an
in-memory API. It is intended for Python workflows, notebooks, tests, and
higher-level applications that need to compose several operations without
writing a parameter file and intermediate raster after every step.

This guide explains how the Python backend is implemented, what the public API
means, and how to change either safely. For callable details, use the
[Python API reference](../python_api/index.md). For task-oriented examples, use
the [Python guides](../python/index.md).

## When To Use Python Or A CLI

Use the Python API when:

- rasters and result objects should remain in memory;
- several native operations need to be composed programmatically;
- an analysis belongs in a notebook, Python application, or runtime test;
- NumPy, pandas, or GeoPandas integration is useful.

Use the vanilla `lsdtt-*` CLIs when:

- the workflow is file-oriented and driven by parameter files;
- it runs from a shell, scheduler, container, or Snakemake rule;
- reproducible command-line execution is more important than in-memory
  composition;
- the user should not need a Python geospatial environment.

A binding should expose reusable C++ behavior with a coherent Python contract.
Do not bind a CLI's `main` function or duplicate orchestration merely to make
the command appear in Python. Wheels separately bundle vanilla CLI executables
and dispatch their console scripts through `_cli.py`.

## Layers And Responsibilities

| Layer | Responsibility |
|---|---|
| `src/core`, `src/cpu`, `src/utils` | Native implementation and scientific behavior |
| `python_bindings/bindings/` | C++ nanobind adapters that expose native behavior through `_core` |
| `python_bindings/src/lsdtt3/` | Public package exports, Python wrappers, accelerators, runtime setup, and CLI dispatch |
| Public Python API | Supported names, signatures, defaults, objects, and behavior reached through `import lsdtt3` |
| `python_bindings/api/` | Machine-readable API snapshots, baselines, schemas, and reviewed decisions |
| `python_bindings/tools/` | Developer programs that generate or check those contracts |
| Generated `docs/python_api/*.md` | Complete runtime-manifest reference; signatures are not handwritten |
| Bundled `lsdtt-*` drivers | File-oriented executables inside wheels; not in-memory bindings |

The backend is how behavior is implemented and exposed. The API is the public
contract callers depend on. `python_bindings/api/` records that contract; it is
not a second backend implementation.

## The Native Backend

`python_bindings/CMakeLists.txt` builds a shared native runtime from the same C++
sources used by the main project. `python_bindings/bindings/CMakeLists.txt` then
builds the nanobind `_core` extension.

Feature-specific files under `bindings/` perform narrow adaptation:

- convert `Raster` and metadata objects;
- validate grids, masks, arrays, strings, and options at the Python boundary;
- map Python objects to native requests/options;
- return Python-facing result classes without reimplementing algorithms.

Nanobind declarations are authoritative for native runtime signatures. Public
native functions expose structured overload data through `__nb_signature__`.

## The Public Package

`python_bindings/src/lsdtt3/__init__.py` imports supported native objects, adds
Python accelerators, and defines `__all__`. `__all__` is the supported top-level
namespace and is checked against the committed API manifest.

Python modules in the package have separate roles:

- `_graphflood_accelerators.py` composes native GraphFlood backends;
- `_compat.py` contains private keyword-migration plumbing;
- `_cli.py` forwards installed console scripts to bundled executables;
- `_llm.py` provides compact, checked API cards;
- `generate_api_reference.py` renders API pages from the committed runtime manifest.

Runtime nanobind signatures and Python functions are authoritative. MkDocs reads
their committed generated representation and does not import `_core`.

## API Contract Files

`python_bindings/api/` contains reviewed build artifacts:

- `public_api.json`: all exported classes/functions, public class members,
  overload signatures, defaults, documentation coverage, console scripts, and
  declared drivers;
- `decisions.json`: release policy, approved CLI/Python differences, completed
  migrations, intentional default differences, known failures, and optional
  features;
- `graphflood_backend_options.json`: accepted built-in backend options and the
  accelerator migration policy;
- `naming_baseline.json`: current naming violations; it may only shrink;
- generated domain pages under `docs/python_api/`: every export exactly once.

These files make public changes visible in review. An implementation-only change
should leave them unchanged. A deliberate API change updates code, tests,
decision records, generated contracts, documentation, and release notes
together.

## Building Locally

The binding must use matching GDAL headers and libraries. Mixing system GDAL
headers with a Pixi/Conda GDAL library can link successfully and still produce
an extension that fails at import with unresolved symbols.

One verified development configuration is:

```bash
pixi run -e docs cmake -S python_bindings -B build/python-api \
  -DCMAKE_BUILD_TYPE=Release \
  -DLSDTT_PYTHON_CLI=OFF \
  -DGDAL_INCLUDE_DIR=/usr/include/gdal \
  -DGDAL_LIBRARY=/usr/lib/libgdal.so
pixi run -e docs cmake --build build/python-api -j2
pixi run -e docs cmake --install build/python-api \
  --prefix /tmp/lsdtt3-python-api
```

Run against that direct install by activating its import path (and, when using
Pixi's Python with system GDAL, the matching runtime libraries):

```bash
export PYTHONPATH=/tmp/lsdtt3-python-api
export LD_LIBRARY_PATH=/path/to/the/python/environment/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
python -c "import lsdtt3; print(lsdtt3.__version__)"
```

Use paths appropriate to the machine, but keep header and library provenance
consistent. Release wheels use the platform build in cibuildwheel instead.

`LSDTT_PYTHON_COMPATIBILITY_PROBE=ON` enables private Phase 0 test hooks for
keyword migration and supplied-valley-mask characterization. It is off by
default and does not enter release wheels.

## Optional Dependencies

The standalone binding build discovers optional native capabilities at CMake
configure time:

- **FFTW** enables spectral routines and the bundled
  `lsdtt-spectral-analysis` driver. The Python console-script declaration is
  stable, but installed-wheel checks allow the native driver to be absent only
  when `spectral_available()` reports false.
- **Eigen** enables chi/ksn inversion functionality.
- **netCDF** enables atmospheric/palaeomagnetic cosmogenic support where the
  corresponding native data are available.
- **OpenMP** enables native parallel execution; it may be disabled on platforms
  where the wheel toolchain cannot provide a compatible runtime.

Optional-feature decisions live in `api/decisions.json`. Public symbols must not
silently disappear without a manifest/decision update. Runtime availability
functions or direct actionable errors are preferable to an import that fails
halfway through package initialization.

## Using An Existing Binding

After installing the package:

```python
import lsdtt3

dem = lsdtt3.Raster.read("dem.tif")
filled = lsdtt3.fill_pits(dem)
flow = lsdtt3.FlowInfo(filled, boundary_conditions="oooo")
area = flow.drainage_area()
area.write("drainage_area.tif")
```

Prefer result objects and native composition over writing and re-reading
temporary files.

## Adding A Binding

1. Confirm the underlying behavior belongs in reusable C++, not in a CLI.
2. Add the smallest adapter in the appropriate `bindings/*.cpp` file.
3. Validate Python-facing arrays, units, grids, enums, and optional values at
   the binding seam.
4. Export the object deliberately through `__init__.py` and `__all__`.
5. Add runtime tests that exercise behavior, errors, signatures, and defaults.
6. Add authored workflow guidance where useful and regenerate the API pages.
7. Run the manifest and naming checks. Update generated files only for an
   intentional public change.
8. Add a decision record and release note for compatibility-sensitive changes.

Do not add a second Python implementation of a native algorithm to avoid writing
a binding. Python composition is appropriate; duplicated scientific numerics
are not.

## Changing A Public Signature Or Default

1. Add a sensitive pre-change behavioral oracle.
2. Record the current signature/default in `public_api.json` and the decision in
   `decisions.json`.
3. Follow the shipped-API compatibility policy: canonical and legacy names,
   warning, both-key conflict, and an earliest removal release.
4. Preserve positional behavior and native class identity.
5. Keep renames separate from default/numerical changes.
6. Update runtime docs, facade, API pages, examples, LLM cards, and tests.
7. Regenerate/check all contracts and include a release note.

Functions and native constructors use different compatibility seams. Python
wrappers can normalize kwargs before forwarding to an unwrapped native
function. Native classes require nanobind constructor overloads/placement
construction so `type` and `isinstance` remain correct.

## Verification

Source-only checks do not need `_core`:

```bash
python3 -m pytest -q \
  python_bindings/tests/test_binding_source_contract.py \
  python_bindings/tests/test_cli_packaging.py
```

With a current installed extension:

```bash
python3 python_bindings/tools/generate_public_api_manifest.py --check
python3 python_bindings/tools/check_api_naming.py
python3 python_bindings/tools/generate_api_reference.py --check
python3 python_bindings/tools/check_api_examples.py
python3 python_bindings/tools/check_llm_cards.py
python3 -m pytest -q python_bindings/tests
```

Documentation:

```bash
pixi run docs-build
```

GitHub builds run only through the manually triggered Python version-increment
workflow. Before requesting one, run these checks locally. Cibuildwheel then
repeats the runtime contracts against installed platform wheels and additionally
checks bundled executables and console scripts.
