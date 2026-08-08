# GPU compute backend

lsdtt3 can offload a small, growing set of heavy per-pixel compute kernels to a
CUDA GPU while everything else — file I/O, orchestration, and any unported
analysis — continues to run on the CPU. This document describes how that works
and how to use it.

> **Status (2026-06):** the first kernel ported is `direct_polyfit` (local
> quadratic surface fitting). The GPU code is *design-ahead*: it is wired into
> the build and the dispatch path, but the CUDA kernel has not yet been compiled
> or parity-tested on real hardware. The default build (`-DLSDTT_GPU=OFF`) does
> not compile any CUDA and is unaffected.

---

## Two independent GPU axes

There are two separate, unrelated notions of "GPU" in the codebase. Do not
confuse them:

| Axis | Macro | What it is | Status |
|---|---|---|---|
| GPU **compute kernels** | `LSDTT_GPU_KERNELS` | Selective-dispatch kernels for heavy free functions (`direct_polyfit`, future roughness metrics) | being built |
| GPU **raster-engine backend** | `LSDTT_GPU_ENABLED` | A full `lsdtt_raster_base` GPU implementation returned by `make_raster_engine()` | future, not built |

This document is about the **compute kernel** axis. Both are enabled by the same
CMake option (`-DLSDTT_GPU=ON`) when the raster backend eventually lands, but for
now only the kernel library (`lsdtt3_gpu`) is built and only `LSDTT_GPU_KERNELS`
is defined.

`src/core/` stays free of any CUDA dependency: the dispatch decision is made in
`lsdtt3_cpu`, which links the kernel library `lsdtt3_gpu`.

---

## The global compute-backend toggle

Backend selection is a process-wide **preference**, not a hard mode switch:

```cpp
#include "lsdtt_backend_factory.hpp"

lsdtt::set_compute_backend(lsdtt::backend_type::gpu);   // prefer GPU
lsdtt::set_compute_backend(lsdtt::backend_type::cpu);   // force CPU
lsdtt::set_compute_backend(lsdtt::backend_type::auto_detect);  // default
```

Each dispatched function resolves the preference like this:

```
if (preference != cpu) AND (a GPU kernel exists for THIS function)
   AND (a device is present)   ->  run the GPU kernel
else                           ->  run the CPU implementation
```

Consequences:

- A function with **no GPU port ignores the toggle** and runs on the CPU.
- **File I/O never consults the toggle** — `lsdtt_data_loader` is always CPU/GDAL.
- A failed GPU launch (no device, allocation error, CUDA fault) **falls back to
  the CPU** automatically; correctness never depends on the GPU path.
- A typical "GPU run" is a **mix**: CPU does I/O and orchestration, GPU runs the
  ported kernels, unported work stays on CPU.

The preference is global mutable state. Set it once at start-up (e.g. from a CLI
parameter); do not flip it from multiple threads during a dispatched call.

### Host ↔ device transfers

`lsdtt_raster_cpu` still lives in host memory, but the direct-polyfit path now
uses a staged GPU bundle to cut redundant transfers:

- `direct_polyfit` can keep its coefficient rasters (`a..f` + mask) resident on
  device in a `polyfit_device_bundle`.
- `compute_polyfit_metrics` then computes requested metrics on that bundle and
  downloads **one metric at a time** for writing.
- The bundle owns reusable per-metric device output buffers, so metric loops do
  not `cudaMalloc`/`cudaFree` temporary outputs every iteration.
- If the bundle path is unavailable, the code falls back to the stateless
  per-call launcher path and then to CPU if needed.

This keeps memory bounded for `polyfit_metrics: all` on large rasters while
avoiding repeated host->device uploads of the same coefficient rasters.

---

## Using it from a CLI

`lsdtt-surface-metrics` exposes the toggle via the `compute_backend` parameter:

```
compute_backend: auto    # auto (default) | cpu | gpu
```

`auto`/`gpu` route `direct_polyfit` (used by the `polyfit_method: direct` path)
to the GPU when available. Polyfit-derived metrics (`compute_polyfit_metrics`)
are also GPU-dispatched when possible (metric-by-metric download), with CPU
fallback on failure. In a `-DLSDTT_GPU=OFF` build the parameter is accepted but
has no effect.

---

## Building with GPU support

```bash
cmake -S . -B build -DLSDTT_GPU=ON
cmake --build build -j
```

Requires a CUDA toolkit (`nvcc`) and `CUDAToolkit` discoverable by CMake. The
default (`-DLSDTT_GPU=OFF`) needs no CUDA and compiles no `.cu` files.

---

## Architecture of the dispatch seam

```
src/core/lsdtt_backend_factory.{hpp,cpp}
    backend_type, set_compute_backend(), get_compute_backend()   (pure host)

src/gpu/lsdtt_surface_gpu.hpp        POD interface (no CUDA/lsdtt types)
    lsdtt::gpu::gpu_available()
    lsdtt::gpu::direct_polyfit_gpu(raw pointers + dims) -> bool
    lsdtt::gpu::direct_polyfit_gpu_bundle_create(...) -> bundle handle
    lsdtt::gpu::compute_polyfit_metric_from_bundle_gpu(...) -> bool
    lsdtt::gpu::destroy_polyfit_device_bundle(...)

src/gpu/lsdtt_surface_gpu.cu         CUDA kernels + launchers   (lsdtt3_gpu)

src/cpu/lsdtt_raster_surface_metrics_cpu.cpp
    direct_polyfit(): GPU dispatch + optional bundle cache
    compute_polyfit_metrics(): GPU metric dispatch (bundle first, fallbacks)
```

The kernel interface takes only raw row-major pointers and dimensions, so the
GPU library has no dependency on `lsdtt_raster_cpu` or GDAL — the one-way
dependency is `lsdtt3_cpu -> lsdtt3_gpu`.

### Numerical parity

The GPU `direct_polyfit` kernel accumulates its per-window OLS moments in
**double precision** and uses the same decoupled normal-equation algebra as the
CPU path (symmetric centred window). Results should match the CPU output to
floating-point tolerance. A tolerance-based CPU-vs-GPU comparison test must be
added and run on hardware before the GPU path is trusted in production.

---

## Porting another kernel (checklist)

1. Add a POD launcher to `src/gpu/lsdtt_surface_gpu.hpp` + kernel in the `.cu`.
2. In the CPU free function, add an `#ifdef LSDTT_GPU_KERNELS` dispatch block at
   the top that checks `get_compute_backend() != cpu && gpu::gpu_available()`,
   calls the launcher, wraps the result, and returns; otherwise falls through to
   the existing CPU body.
3. No CMake change needed (the kernel goes in the existing `lsdtt3_gpu` sources).
4. Add a parity test against the CPU reference.
