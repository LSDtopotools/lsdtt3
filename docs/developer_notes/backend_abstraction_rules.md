# Backend Abstraction Rules for lsdtt3

This document defines the architectural rules for backend abstraction in `lsdtt3`.

The project does **not** aim to make every analysis uniformly CPU/GPU-agnostic.
Instead, it uses **two different patterns** for two different kinds of components.

These rules are intended to keep the codebase consistent while allowing GPU
support only where it is worthwhile.

---

## Why these rules exist

`lsdtt3` contains several very different kinds of components:

1. **Broad compute substrates**, especially rasters, where many CPU analyses are
   implemented as free functions taking `const lsdtt_raster_cpu&` and returning
   new `lsdtt_raster_cpu` objects.

2. **Persistent domain engines**, such as flow routing and channel/junction
   networks, which are built once and then queried repeatedly through a fairly
   small public API.

These two kinds of components should **not** be forced into the same abstraction
pattern.

The result is two rules:

- **Rule 1 — selective dispatch**
- **Rule 2 — narrow abstract bases**

---

# Rule 1 — Selective dispatch for broad compute substrates

## Summary

Use **concrete CPU objects plus selective dispatch** for components that expose
many analyses or many free functions, where only a **small subset** of those
analyses is likely to get GPU implementations.

The main current example is:

- `lsdtt_raster_cpu`

This rule exists because `lsdtt_raster_cpu` is the substrate for a large number
of CPU-side raster analyses. It does **not** fit naturally into the same narrow,
query-style abstract base pattern used by `lsdtt_flow_base`.

## What this means in practice

For broad compute substrates:

- keep the concrete CPU object
- keep most CPU free functions as ordinary CPU free functions
- add dispatch wrappers **only** around the expensive kernels that are realistic
  GPU targets
- do **not** declare virtual methods for every raster analysis
- do **not** try to force every CLI to become backend-neutral at the level of
  individual raster operations

## Typical examples

Good candidates for selective dispatch include:

- FFT / spectral routines
- wavelet fitting
- large moving-window kernels
- image-style filters such as Frangi-style routines
- any other computationally dominant raster kernel

Poor candidates for selective dispatch include:

- small helper transforms
- low-cost arithmetic or bookkeeping
- niche CPU-only analyses
- orchestration logic in CLIs
- algorithms with no realistic GPU plan

---

## Design pattern for Rule 1

### Public-facing pattern

The public object remains concrete, for example:

```cpp
lsdtt::lsdtt_raster_cpu dem(loader);
```

Heavy kernels are then called through optional dispatch wrappers, for example:

```cpp
auto coeffs = lsdtt::wavelet_quadratic_fit_dispatch(
    dem,
    request.polyfit_lambda_small,
    request.n_threads,
    dispatch_options
);
```

The CPU upload/download and GPU fallback logic remain **inside** the wrapper.

### File placement

For a selectively dispatched raster kernel:

- **Declaration** goes in `src/core/`
- **Dispatch implementation** goes in `src/cpu/`
- **CPU implementation** stays in `src/cpu/`
- **GPU implementation** goes in `src/gpu/` when present

### Example layout

```text
src/core/lsdtt_raster_spectral_dispatch.hpp
src/cpu/lsdtt_raster_spectral_dispatch.cpp
src/cpu/lsdtt_raster_spectral_cpu.hpp
src/cpu/lsdtt_raster_spectral_cpu.cpp
src/gpu/lsdtt_raster_spectral_gpu.hpp
src/gpu/lsdtt_raster_spectral_gpu.cu
```

### Why dispatch implementation belongs in `src/cpu/`

The dispatch implementation usually needs to know about:

- the CPU implementation
- the optional GPU implementation

Putting that implementation in `src/core/` would create an unwanted dependency
from `core` to backend-specific code.

This is the same reason `make_flow_base()` is declared in `src/core/` but
implemented in `src/cpu/`.

---

## Suggested dispatch options type

A small shared dispatch-options struct is useful when an algorithm may choose
between CPU, GPU, or automatic fallback.

For example:

```cpp
#pragma once

#include <cstddef>

#include "lsdtt_backend_factory.hpp"

namespace lsdtt
{

struct lsdtt_dispatch_options
{
    backend_type backend = backend_type::auto_detect;
    std::size_t  min_pixels_for_gpu = 1024ULL * 1024ULL;
    bool         verbose = false;
};

} // namespace lsdtt
```

This can live in `src/core/lsdtt_dispatch_options.hpp`.

### Semantics

- `backend_type::cpu`
  - always use CPU
- `backend_type::gpu`
  - require GPU or throw
- `backend_type::auto_detect`
  - use GPU only if it exists, the algorithm actually has a GPU backend, and the
    problem size is large enough to justify transfer overhead

---

## Example: Rule 1 wrapper declaration

Example declaration in `src/core/lsdtt_raster_spectral_dispatch.hpp`:

```cpp
#pragma once

#include "lsdtt_dispatch_options.hpp"
#include "../cpu/lsdtt_raster_cpu.hpp"
#include "../cpu/lsdtt_raster_spectral_cpu.hpp"

namespace lsdtt
{

poly_fit_coeffs wavelet_quadratic_fit_dispatch(
    const lsdtt_raster_cpu&       raster,
    double                        lambda = 100.0,
    int                           n_threads = 0,
    const lsdtt_dispatch_options& dispatch = lsdtt_dispatch_options{}
);

} // namespace lsdtt
```

---

## Example: Rule 1 wrapper implementation

Example implementation in `src/cpu/lsdtt_raster_spectral_dispatch.cpp`:

```cpp
#include "../core/lsdtt_raster_spectral_dispatch.hpp"
#include "lsdtt_raster_spectral_cpu.hpp"

#ifdef LSDTT_GPU_KERNELS
#include "../gpu/lsdtt_raster_spectral_gpu.hpp"
#endif

namespace lsdtt
{

poly_fit_coeffs wavelet_quadratic_fit_dispatch(
    const lsdtt_raster_cpu&       raster,
    double                        lambda,
    int                           n_threads,
    const lsdtt_dispatch_options& dispatch
)
{
    const std::size_t n_pixels = raster.full_grid_size();

    if (dispatch.backend == backend_type::cpu)
    {
        return wavelet_quadratic_fit(raster, lambda, n_threads);
    }

#ifdef LSDTT_GPU_KERNELS
    const bool gpu_ok =
        gpu_available() &&
        n_pixels >= dispatch.min_pixels_for_gpu;

    if (dispatch.backend == backend_type::gpu)
    {
        return wavelet_quadratic_fit_gpu(raster, lambda, n_threads);
    }

    if (dispatch.backend == backend_type::auto_detect && gpu_ok)
    {
        return wavelet_quadratic_fit_gpu(raster, lambda, n_threads);
    }
#else
    if (dispatch.backend == backend_type::gpu)
    {
        throw std::runtime_error(
            "wavelet_quadratic_fit_dispatch: GPU backend requested but GPU support is not compiled in."
        );
    }
#endif

    return wavelet_quadratic_fit(raster, lambda, n_threads);
}

} // namespace lsdtt
```

---

## Example: Rule 1 integration inside a workflow

A workflow such as channel-bottom detection should usually remain
**CPU-native orchestration**, while selectively dispatching the heavy kernels.

For example:

```cpp
lsdtt::channel_bottom_result find_channel_bottoms(
    const lsdtt_raster_cpu&       dem,
    const channel_bottom_request& request,
    const lsdtt_dispatch_options& dispatch
)
{
    poly_fit_coeffs small_fit =
        wavelet_quadratic_fit_dispatch(
            dem,
            request.polyfit_lambda_small,
            request.n_threads,
            dispatch
        );

    poly_fit_coeffs large_fit =
        wavelet_quadratic_fit_dispatch(
            dem,
            request.polyfit_lambda_large,
            request.n_threads,
            dispatch
        );

    lsdtt_raster_cpu local_mean =
        moving_window_mean_dispatch(
            dem,
            static_cast<double>(request.depth_window_radius) * dem.cell_size(),
            request.depth_window_shape,
            request.n_threads,
            dispatch
        );

    // Remaining scoring logic can stay ordinary CPU code.
    // Thresholding, diagnostics, and small helper steps do not need to be virtual.

    ...
}
```

### Important point

This pattern does **not** mean the whole workflow has become a GPU engine.
It only means the expensive internal steps may be accelerated.

---

## Rule 1 checklist

Use selective dispatch when all of the following are true:

- the component is a broad compute substrate, not a narrow engine
- the CPU API is already large or analysis-heavy
- only a subset of operations is likely to get GPU support
- the GPU target is a kernel or a small family of kernels
- the result can sensibly come back to ordinary host-side objects

If these are not true, consider Rule 2 instead.

---

# Rule 2 — Narrow abstract base for persistent stateful engines

## Summary

Use a **narrow abstract base plus factory** for long-lived domain engines that:

- own significant internal state
- are built or initialised once
- are queried or run through a small stable API
- should hide backend-specific storage from CLIs and high-level workflows

This is the correct pattern for:

- `lsdtt_flow_base`
- `lsdtt_junction_network_base`
- a future basin-manager base
- a future particle-tracer base

## Existing example

`lsdtt_flow_base` already follows this model:

```cpp
std::unique_ptr<lsdtt::lsdtt_flow_base> flow =
    lsdtt::make_flow_base(
        routed_dem.get_data_array(),
        routed_dem.get_mask_array(),
        routed_dem.get_meta(),
        fi_opts
    );
```

The CLI holds the abstract base and does not need to know whether the object is
CPU or GPU-backed.

This is the pattern that should also be applied more consistently to
junction-network, basin-manager, and particle-tracer layers.

---

## What belongs in Rule 2

Use Rule 2 for objects that behave like engines or graphs, for example:

- flow-routing engines
- junction-network engines
- basin-selection / basin-management engines
- particle tracing engines

These objects usually have a lifecycle such as:

1. construct or create
2. build internal state
3. run or select
4. query
5. write outputs

That is a strong sign that they want a base-class interface.

---

## What does **not** belong in Rule 2

Do **not** use a narrow abstract base for:

- every raster analysis
- giant families of unrelated CPU free functions
- simple value containers
- small result structs
- trivial host-side geometry records
- one-off helper algorithms

---

## Design pattern for Rule 2

### Public-facing pattern

The CLI and high-level workflows should hold only the abstract base:

```cpp
std::unique_ptr<lsdtt::lsdtt_junction_network_base> network;
```

or:

```cpp
std::unique_ptr<lsdtt::lsdtt_particle_tracer_base> tracer;
```

Concrete CPU classes should **not** appear in CLI code when the object belongs to
Rule 2.

### Base interface location

The abstract base class lives in `src/core/`.

### Concrete backends

- CPU implementation lives in `src/cpu/`
- GPU implementation lives in `src/gpu/` when present

### Factory

- factory declaration lives in `src/core/`
- factory implementation lives in `src/cpu/`

Again, this avoids a `core -> cpu` dependency.

---

## Rule 2: interface design principles

The public base interface should be:

- narrow
- stable
- behaviour-oriented
- not storage-oriented

Good public methods are things like:

- create
- run
- size
- empty
- query receiver / donor / outlet relationships
- select basins
- write outputs
- export summaries

Bad public methods are things like:

- exposing backend-specific internal arrays
- exposing device buffers
- hundreds of tiny per-element accessors just because CPU has them
- forcing the CLI to know storage details

---

## Example: junction network should follow Rule 2

### Desired CLI usage

Long term, this should look more like:

```cpp
std::unique_ptr<lsdtt::lsdtt_junction_network_base> channel_network =
    lsdtt::make_junction_network_base(
        *flow,
        source_flats,
        jn_opts
    );
```

rather than:

```cpp
lsdtt::lsdtt_junction_network_cpu channel_network;
channel_network.create(*flow, source_flats, jn_opts);
```

### Suggested file layout

```text
src/core/lsdtt_junction_network_base.hpp
src/cpu/lsdtt_junction_network_cpu.hpp
src/cpu/lsdtt_junction_network_cpu.cpp
src/cpu/lsdtt_junction_network_base.cpp
src/gpu/lsdtt_junction_network_gpu.hpp
src/gpu/lsdtt_junction_network_gpu.cpp
```

### Factory declaration example

In `src/core/lsdtt_junction_network_base.hpp`:

```cpp
std::unique_ptr<lsdtt_junction_network_base> make_junction_network_base(
    const lsdtt_flow_base&                    flow,
    const std::vector<std::int64_t>&          source_flats,
    const lsdtt_junction_network_options&     options = lsdtt_junction_network_options{}
);
```

### Factory implementation example

In `src/cpu/lsdtt_junction_network_base.cpp`:

```cpp
#include "../core/lsdtt_junction_network_base.hpp"
#include "lsdtt_junction_network_cpu.hpp"

namespace lsdtt
{

std::unique_ptr<lsdtt_junction_network_base> make_junction_network_base(
    const lsdtt_flow_base&                    flow,
    const std::vector<std::int64_t>&          source_flats,
    const lsdtt_junction_network_options&     options
)
{
    auto network = std::make_unique<lsdtt_junction_network_cpu>();
    network->create(flow, source_flats, options);
    return network;
}

} // namespace lsdtt
```

A future GPU implementation can be added in the same way `make_flow_base()` can
later choose a GPU flow engine.

---

## Example: basin layer should follow Rule 2 at the manager level

For basin work, the public engine is usually the **manager**, not necessarily the
single basin object.

That means the better abstraction target is usually a future:

- `lsdtt_junction_basin_manager_base`

rather than forcing the CLI to work directly with a concrete
`lsdtt_junction_basin_manager_cpu`.

### Why manager-level abstraction is preferable

The manager is the thing that typically handles:

- basin selection
- nested basin expansion
- outlet-based selection
- output writing
- collection-level query methods

That makes it the natural engine object.

A single basin object can remain concrete if it is mostly a data holder.

### Desired usage

```cpp
std::unique_ptr<lsdtt::lsdtt_junction_basin_manager_base> basin_manager =
    lsdtt::make_junction_basin_manager_base(
        *flow,
        *channel_network
    );

basin_manager->select_all_junction_basins();
```

---

## Example: particle tracing should start with Rule 2 from day one

Particle tracing is expected to be one of the most compute-heavy parts of the
codebase, and is a strong future GPU target.

That makes it a poor candidate for ad hoc raster-style selective dispatch at the
top level. Instead, the **particle tracer itself** should be designed as a
backend-hidden engine from the beginning.

### Why particle tracing belongs in Rule 2

A particle tracer is typically:

- stateful
- long-lived during a run
- backend-sensitive in storage layout
- likely to manage large internal arrays or device buffers
- naturally driven by a small number of engine methods

### Suggested base class name

- `lsdtt_particle_tracer_base`

### Suggested file layout

```text
src/core/lsdtt_particle_tracer_base.hpp
src/cpu/lsdtt_particle_tracer_cpu.hpp
src/cpu/lsdtt_particle_tracer_cpu.cpp
src/cpu/lsdtt_particle_tracer_base.cpp
src/gpu/lsdtt_particle_tracer_gpu.hpp
src/gpu/lsdtt_particle_tracer_gpu.cu
```

### Suggested usage

```cpp
std::unique_ptr<lsdtt::lsdtt_particle_tracer_base> tracer =
    lsdtt::make_particle_tracer_base(
        *flow,
        tracer_options
    );

tracer->seed_particles(seed_points);
tracer->run();
tracer->write_particle_points(out_path, format);
```

### Important interface advice

Do **not** make the particle-tracer base enormous. It should expose a narrow
engine API, for example:

- `seed_particles(...)`
- `run()`
- `empty()`
- `n_particles()`
- `n_active_particles()`
- `write_particle_points(...)`
- `write_particle_lines(...)`
- `particle_density_raster()`
- `travel_time_raster()`

Avoid exposing large numbers of fine-grained per-particle getters unless they
are truly necessary for all backends.

---

## Rule 2 checklist

Use a narrow abstract base when all of the following are true:

- the object is long-lived and stateful
- the object is built or initialised once, then queried or run repeatedly
- CPU and GPU storage layouts may differ significantly
- the public API can be kept reasonably small and stable
- the CLI does not need to know the concrete backend type

If these are not true, consider Rule 1 instead.

---

# Choosing between Rule 1 and Rule 2

Use this decision guide.

## Choose Rule 1 if:

- the main public object is a broad compute substrate
- many analyses already exist as CPU free functions
- only a small subset of those analyses may ever use GPU
- forcing all analyses behind virtual methods would create an oversized API

**Example:** `lsdtt_raster_cpu`

## Choose Rule 2 if:

- the object is a graph, engine, manager, or simulator
- the object owns important hidden state
- the public API is stable and query-oriented or run-oriented
- the object should be usable in CLIs through `std::unique_ptr<..._base>`

**Examples:** flow, junction network, basin manager, particle tracer

---

# Common mistakes to avoid

## Mistake 1 — trying to make `lsdtt_raster_cpu` look like `lsdtt_flow_base`

Do **not** try to put virtual methods for every raster analysis onto a raster
base or analysis base. The API surface is too wide and the value is too low.

## Mistake 2 — exposing concrete CPU engines in CLIs for Rule 2 objects

If the object is a flow-like engine, the CLI should not construct the concrete
CPU class directly.

Bad:

```cpp
lsdtt::lsdtt_junction_network_cpu network;
```

Better:

```cpp
std::unique_ptr<lsdtt::lsdtt_junction_network_base> network =
    lsdtt::make_junction_network_base(...);
```

## Mistake 3 — inventing GPU parity for everything

GPU support is expected for a **small subset** of the codebase. Do not create
abstraction layers for routines that have no credible GPU path.

## Mistake 4 — making base interfaces too wide

Rule 2 works only when the base interface stays narrow and stable.

## Mistake 5 — mixing storage access into public base interfaces

Public bases should expose behaviour and queries, not backend-specific buffers.

---

# CLI guidance

## For Rule 1 components

CLIs may still use concrete CPU objects such as:

```cpp
lsdtt::lsdtt_raster_cpu dem(loader);
```

and then pass those into selectively dispatched kernels or workflows.

This is acceptable and expected.

## For Rule 2 components

CLIs should hold abstract pointers:

```cpp
std::unique_ptr<lsdtt::lsdtt_flow_base> flow;
std::unique_ptr<lsdtt::lsdtt_junction_network_base> network;
std::unique_ptr<lsdtt::lsdtt_junction_basin_manager_base> basin_manager;
std::unique_ptr<lsdtt::lsdtt_particle_tracer_base> tracer;
```

This keeps CLIs backend-agnostic where it actually matters.

---

# CMake guidance

## Rule 1 dispatch wrappers

Add dispatch implementation sources to `lsdtt3_cpu`, for example:

```text
src/cpu/lsdtt_raster_spectral_dispatch.cpp
src/cpu/lsdtt_raster_surface_metrics_dispatch.cpp
```

Add optional GPU implementations to `lsdtt3_gpu`, for example:

```text
src/gpu/lsdtt_raster_spectral_gpu.cu
src/gpu/lsdtt_raster_surface_metrics_gpu.cu
```

## Rule 2 factories

Factory implementation sources should also live in `lsdtt3_cpu`, because they
instantiate CPU implementations and may later dispatch to GPU implementations.

---

# Final rule statement

`lsdtt3` uses **two backend-abstraction patterns**.

## Rule 1
**Broad compute substrates use concrete CPU objects plus selective dispatch of a
small number of heavy kernels.**

Current main example:
- `lsdtt_raster_cpu`

## Rule 2
**Persistent stateful engines use narrow abstract bases in `src/core/` plus
backend factories implemented in `src/cpu/`.**

Current and intended examples:
- `lsdtt_flow_base`
- `lsdtt_junction_network_base`
- future `lsdtt_junction_basin_manager_base`
- future `lsdtt_particle_tracer_base`

These two rules are deliberate. The codebase should not try to force both kinds
of component into one abstraction style.

---

## Relationship to `gpu_compute.md`

This note is **design-ahead prescription**: it states how backend dispatch
should be structured as GPU work lands. [`gpu_compute.md`](gpu_compute.md)
**describes what is actually built** — the process-global
`set_compute_backend()` / `get_compute_backend()` preference resolved inside
`direct_polyfit`, and the `LSDTT_GPU_KERNELS` compile definition that is the
only GPU macro CMake ever defines.

Where the two differ, `gpu_compute.md` is authoritative for present
behaviour. In particular the `lsdtt_dispatch_options` struct and the
`wavelet_quadratic_fit_dispatch` entry point used as examples above do not
exist in the tree: they illustrate the shape of the pattern, not code you can
call today.
