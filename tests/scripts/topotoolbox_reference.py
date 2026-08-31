#!/usr/bin/env python3
"""TopoToolbox reference implementations for the flow-routing and
surface-metrics benchmarks.

Sibling of `richdem_reference.py`, and for flow routing the more trustworthy of
the two.

WHY BOTH REFERENCES

RichDEM's D8 selects the neighbour with the LOWEST ELEVATION, not the steepest
DESCENT (drop / distance). Wherever the lowest neighbour is diagonal but a
cardinal neighbour is steeper, the two rules disagree, so its accumulation is
kept as a *report-only* reference.

(Corrected 2026-08-21: this previously said RichDEM has a diagonal-preferring
tie-break and fails a planar ramp. Both are wrong -- its tie-break prefers
cardinals, and on a 50x50 north-south ramp it returns cardinal south for all
2304 interior cells. The divergence comes from the selection rule and needs no
tie at all. See docs/testing_and_benchmarking.md 5.1.)

TopoToolbox computes true steepest descent -- `libtopotoolbox/src/flow_routing.c`
divides the drop by a chamfer distance (1 for cardinals, sqrt(2) for diagonals)
before comparing -- which is the same rule lsdtt3 uses. It reproduces the planar
ramp exactly, so it *is* a valid correctness reference for D8 and its
comparisons are asserted.

See docs/testing_and_benchmarking.md for the measurements behind that.

INSTALL

topotoolbox is not on conda-forge; it comes from PyPI via the `reference`
feature in pixi.toml, which also pins Python to >=3.11,<3.14 because the package
ships wheels for cp311-cp313 only. Run under:

    pixi run -e benchmark python tests/scripts/topotoolbox_reference.py ...

Routine -> TopoToolbox mapping:

  Flow routing (docs/testing_and_benchmarking.md, flow plan):
    fill               GridObject.fillsinks()   plain fill, no epsilon
    accum_d8           FlowObject(method='d8').flow_accumulation()

  Surface metrics (surface-metrics plan Tier B -- SPEED ONLY):
    evansslope         GridObject.evansslope(modified=False)
    evansslope_partials  ...(partial_derivatives=True) -> (fx, fy)
    curvature_profc    GridObject.curvature(ctype='profc')
    curvature_planc    GridObject.curvature(ctype='planc')
    curvature_tangc    GridObject.curvature(ctype='tangc')
    curvature_meanc    GridObject.curvature(ctype='meanc')
    curvature_all4     the four above, back to back, one shared timing
    aspect             GridObject.aspect()                 REPORT-ONLY
    hillshade_fused    GridObject.hillshade(fused=True)    libtopotoolbox
    hillshade_unfused  GridObject.hillshade(fused=False)   numpy

SURFACE-METRIC ROUTINES ARE A **TIMING** REFERENCE, NOT A CORRECTNESS ONE.
Correctness for the polyfit chain comes from the analytic Taylor oracle in
`lsdtt_raster_synthetic_surfaces_cpu.hpp` (Tier A), which is exact and needs no
network. Nothing here should ever be turned into a value assertion: sign
conventions and aspect handedness deliberately differ between the two tools and
matching them is explicitly out of scope.

What these routines actually run, verified against topotoolbox 0.0.12 sources
(grid_object.py), because it decides what a timing ratio MEANS:

  * `evansslope` and `curvature` are `scipy.ndimage.convolve` calls with the
    Evans 3x3 kernels, in Python, SINGLE-THREADED. The mathematics is the same
    computation lsdtt3's `direct` polyfit does at half-width W = 1, so this is
    a genuine equivalence of work -- but it is scipy-serial vs OpenMP-C++, so
    the single-thread lsdtt3 number is the honest headline and the threaded one
    is the practical one. See run_surface_benchmark.py.
  * `curvature` runs FIVE convolutions (fx, fy, fxx, fyy, fxy) per call and
    throws four of the five products away when you ask for the next ctype;
    `evansslope` runs two (fx, fy). `n_convolutions` is reported so the caller
    can normalise. lsdtt3 solves for all six quadratic coefficients once.
  * `evansslope` additionally runs `distance_transform_edt` to in-fill NaNs
    before convolving -- real work lsdtt3 does not do. It is timed inside the
    call because it is not separable from it; `seconds_nan_fill` reports the
    same operation measured on its own so it can be discounted.
  * `aspect` is `np.gradient(edge_order=2)` -- plain central differences, less
    work than a quadratic fit. REPORT-ONLY for timing; never quote a ratio.
  * `hillshade(fused=True)` goes through libtopotoolbox (compiled). It is the
    only compiled-vs-compiled like-for-like comparison in the whole API.
  * `gradient8` is deliberately absent: it is the D8 steepest-neighbour drop,
    a different quantity, so a timing ratio against a fitted-surface gradient
    would compare two unlike computations.

Usage:
    python topotoolbox_reference.py --routine accum_d8 --in dem.tif --out ref.tif
    python topotoolbox_reference.py --routine curvature_profc --in dem.tif \\
        --out ref.tif --repeats 7
"""

import argparse
import contextlib
import json
import os
import resource
import statistics
import sys
import time

FLOW_ROUTINES = ("fill", "accum_d8")

# Surface-metric routines. Each is timed best-of-N after a discarded warm-up
# (see --repeats); the flow routines are single-shot because they are the
# expensive ones and the driver caches their output.
SURFACE_ROUTINES = (
    "evansslope",
    "evansslope_partials",
    "curvature_profc",
    "curvature_planc",
    "curvature_tangc",
    "curvature_meanc",
    "curvature_all4",
    "aspect",
    "hillshade_fused",
    "hillshade_unfused",
)

ROUTINES = FLOW_ROUTINES + SURFACE_ROUTINES

# Routines whose timing must never be quoted as a ratio against lsdtt3: they do
# demonstrably less work than a quadratic surface fit.
REPORT_ONLY_ROUTINES = frozenset({"aspect"})

# Number of scipy convolutions each surface routine performs, so a caller can
# normalise "seconds per convolution" rather than pretending one call of
# `curvature` and one call of `evansslope` are the same amount of work.
N_CONVOLUTIONS = {
    "evansslope": 2,
    "evansslope_partials": 2,
    "curvature_profc": 5,
    "curvature_planc": 5,
    "curvature_tangc": 5,
    "curvature_meanc": 5,
    "curvature_all4": 20,
    "aspect": 0,          # np.gradient, not a convolution
    "hillshade_fused": 0,
    "hillshade_unfused": 0,
}


def peak_rss_mb() -> float:
    """Peak resident set size of this process, in MiB.

    ru_maxrss is in kibibytes on Linux and bytes on macOS. The Linux
    interpretation is used because that is where the benchmark is run; on macOS
    the number will read ~1024x high, which is obvious rather than subtle.
    """
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


@contextlib.contextmanager
def suppressed_output():
    """Silence any progress output from the underlying C library.

    Redirects the file descriptors, not just sys.stdout, because the work
    happens in libtopotoolbox and would otherwise write past a Python-level
    capture and corrupt the JSON line this script emits.
    """
    saved = os.dup(1), os.dup(2)
    devnull = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull, 1)
        os.dup2(devnull, 2)
        yield
    finally:
        os.dup2(saved[0], 1)
        os.dup2(saved[1], 2)
        os.close(devnull)
        os.close(saved[0])
        os.close(saved[1])


def _surface_call(routine: str, dem, azimuth: float, altitude: float):
    """Return a zero-argument callable performing one surface-metric routine.

    Building the callable outside the timing loop keeps attribute lookup and
    argument marshalling out of the measurement, which matters here because the
    routines themselves are milliseconds on a small fixture.
    """
    if routine == "evansslope":
        return lambda: dem.evansslope(modified=False)
    if routine == "evansslope_partials":
        # modified=False for the same reason as above: True adds a Shary
        # smoothing convolution that lsdtt3 does not perform, which would show
        # up as an unearned lsdtt3 win.
        return lambda: dem.evansslope(partial_derivatives=True, modified=False)
    if routine.startswith("curvature_") and routine != "curvature_all4":
        ctype = routine.split("_", 1)[1]
        return lambda: dem.curvature(ctype=ctype, meanfilt=False)
    if routine == "curvature_all4":
        # The point of this routine: TopoToolbox recomputes all five Evans
        # convolutions for EVERY curvature type, where one lsdtt3 polyfit
        # yields all of them. Four calls is what a user actually pays to get
        # profile + plan + tangential + mean curvature.
        return lambda: [dem.curvature(ctype=c, meanfilt=False)
                        for c in ("profc", "planc", "tangc", "meanc")]
    if routine == "aspect":
        return lambda: dem.aspect()
    if routine == "hillshade_fused":
        return lambda: dem.hillshade(azimuth=azimuth, altitude=altitude,
                                     fused=True)
    if routine == "hillshade_unfused":
        return lambda: dem.hillshade(azimuth=azimuth, altitude=altitude,
                                     fused=False)
    raise ValueError(f"unknown surface routine {routine!r}")


def run_surface(routine: str, dem, repeats: int,
                azimuth: float = 315.0, altitude: float = 60.0):
    """Best-of-N timing for one surface routine, after a discarded warm-up.

    Best-of-N, not the mean: the fastest run is the one least polluted by
    scheduler noise and page faults, and a mean over a cold first call measures
    the allocator. The median is reported alongside so a single lucky run
    cannot be quoted on its own.
    """
    import numpy as np

    call = _surface_call(routine, dem, azimuth, altitude)

    with suppressed_output():
        result = call()                      # warm-up, discarded

    timings = []
    for _ in range(max(1, repeats)):
        t0 = time.perf_counter()
        with suppressed_output():
            result = call()
        timings.append(time.perf_counter() - t0)

    info = {
        "routine": routine,
        "seconds": min(timings),
        "seconds_median": statistics.median(timings),
        "seconds_all": timings,
        "repeats": len(timings),
        "n_convolutions": N_CONVOLUTIONS.get(routine),
        "report_only": routine in REPORT_ONLY_ROUTINES,
        "threads": 1,               # scipy convolve and np.gradient are serial
        "peak_rss_mb": peak_rss_mb(),
    }

    # evansslope in-fills NaNs with distance_transform_edt before convolving.
    # That is real work lsdtt3 does not do, and on a fixture with a large
    # nodata margin it can dominate. Time it separately so the caller can say
    # how much of the evansslope number is not the convolution.
    if routine.startswith("evansslope"):
        from scipy.ndimage import distance_transform_edt
        z = np.asarray(dem.z)
        # Best-of-N, matching how `seconds` itself is measured. Timing this
        # once while the total is a best-of-5 would mix a slow sample with a
        # fast one and the difference -- which is the quantity the caller
        # actually uses -- could come out near zero or negative purely from
        # that mismatch.
        fills = []
        for _ in range(max(1, repeats)):
            t0 = time.perf_counter()
            distance_transform_edt(np.isnan(z), return_distances=False,
                                   return_indices=True)
            fills.append(time.perf_counter() - t0)
        info["seconds_nan_fill"] = min(fills)
        info["seconds_nan_fill_median"] = statistics.median(fills)

    return info, result


def _repeat_time(call, repeats: int) -> dict:
    """Best-of-N timing for a flow routine, with the result carried along.

    Returns the timing fields plus `_result`, which the caller pops. Bundling
    the result with the timing keeps the callable's construction and its
    invocation in one place, so a future routine cannot accidentally time the
    construction too.

    repeats=1 does exactly one call with no warm-up, preserving the historical
    single-shot numbers that run_flow_benchmark.py's tolerances were pinned
    against.
    """
    import statistics as _st

    if repeats > 1:
        with suppressed_output():
            call()                       # warm-up, discarded

    timings = []
    for _ in range(max(1, repeats)):
        t0 = time.perf_counter()
        with suppressed_output():
            result = call()
        timings.append(time.perf_counter() - t0)

    return {
        "seconds": min(timings),
        "seconds_median": _st.median(timings),
        "seconds_all": timings,
        "seconds_stdev": _st.stdev(timings) if len(timings) > 1 else 0.0,
        "repeats": len(timings),
        "_result": result,
    }


def run(routine: str, in_path: str, out_path: str, repeats: int = 1,
        azimuth: float = 315.0, altitude: float = 60.0):
    import numpy as np
    import topotoolbox as ttb

    t_read = time.perf_counter()
    with suppressed_output():
        dem = ttb.read_tif(in_path)
    seconds_read = time.perf_counter() - t_read

    info = {"reference": "topotoolbox", "routine": routine,
            "seconds_read": seconds_read}

    if routine in SURFACE_ROUTINES:
        # update, not rebind: run_surface returns a fresh dict, and rebinding
        # would drop the `reference` and `seconds_read` fields set above --
        # which run_benchmark_suite.py reads to keep I/O out of compute.
        surface_info, result = run_surface(routine, dem, repeats, azimuth,
                                           altitude)
        info.update(surface_info)
        info["azimuth"] = azimuth
        info["altitude"] = altitude
        # Multi-output routines: score the correlation guard on one
        # representative field. fx for the partials, mean curvature for the
        # four-curvature sweep.
        if isinstance(result, (tuple, list)):
            result = result[0] if routine == "evansslope_partials" else result[-1]
        info["dem_dtype"] = str(np.asarray(dem.z).dtype)
        info["cellsize"] = float(dem.cellsize)

    elif routine == "fill":
        # repeats=1 reproduces the historical single-shot behaviour exactly, so
        # run_flow_benchmark.py's pinned tolerances and cached numbers are
        # untouched. repeats>1 adds a discarded warm-up and reports best-of-N,
        # which is what run_benchmark_suite.py asks for so that this reference
        # is timed on the same basis as the other two rather than being
        # single-shot against their best-of-N.
        info.update(_repeat_time(lambda: dem.fillsinks(), repeats))
        result = info.pop("_result")

    elif routine == "accum_d8":
        # Timed in two phases because they are separately meaningful and
        # because RichDEM's single call does directions+accumulation without
        # sink resolution -- so `seconds` (build + accumulate) is the
        # like-for-like number, while seconds_accumulate_only is comparable
        # with lsdtt3's in-process accumulation timing.
        #
        # sink_resolution='carve' is a no-op on an already-filled DEM, which is
        # what the benchmark feeds it, so it does not silently modify the
        # surface being compared.
        build = _repeat_time(
            lambda: ttb.FlowObject(dem, method="d8", sink_resolution="carve"),
            repeats)
        flow = build.pop("_result")
        accum = _repeat_time(lambda: flow.flow_accumulation(), repeats)
        result = accum.pop("_result")

        # The phase totals are summed on the SAME statistic (both medians, or
        # both bests) rather than mixing one phase's best with the other's
        # median, which would produce a total no single run ever achieved.
        info["seconds"] = build["seconds"] + accum["seconds"]
        info["seconds_median"] = build["seconds_median"] + accum["seconds_median"]
        info["seconds_all"] = [b + a for b, a
                               in zip(build["seconds_all"], accum["seconds_all"])]
        info["repeats"] = accum["repeats"]
        info["seconds_build_only"] = build["seconds"]
        info["seconds_accumulate_only"] = accum["seconds"]

    else:
        raise ValueError(f"unknown routine {routine!r}")

    if out_path:
        with suppressed_output():
            ttb.write_tif(result, out_path)

    arr = np.asarray(result.z)
    finite = arr[np.isfinite(arr)]
    info.update({
        # Serial for these routines (scipy convolve, np.gradient and
        # libtopotoolbox's single-threaded paths). Stated rather than omitted
        # so a campaign cannot silently attribute a thread count to it.
        "threads": 1,
        "peak_rss_mb": peak_rss_mb(),
        "shape": list(arr.shape),
        "n_cells": int(arr.size),
        # Equivalence-of-work check (plan Sec2.4): a tool that produced far
        # fewer valid cells did less work, and its timing must not be averaged
        # in with the other's as though the comparison were fair.
        "n_valid": int(finite.size),
        "out_dtype": str(arr.dtype),
        "min": float(finite.min()) if finite.size else None,
        "max": float(finite.max()) if finite.size else None,
    })
    return info


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--routine", required=True, choices=ROUTINES)
    ap.add_argument("--in", dest="in_path", required=True)
    ap.add_argument("--out", dest="out_path", default=None,
                    help="write the result raster here. Optional for the "
                         "surface routines, which are timing-only unless the "
                         "caller wants the correlation guard.")
    ap.add_argument("--azimuth", type=float, default=315.0,
                    help="hillshade illumination azimuth, degrees. NOTE: "
                         "TopoToolbox and lsdtt3 do not use the same azimuth "
                         "origin -- measured on a real fixture, "
                         "ttb_azimuth = lsdtt3_azimuth + 90 (mod 360) gives "
                         "Spearman rho = 1.0000. See run_surface_benchmark.py.")
    ap.add_argument("--altitude", type=float, default=60.0,
                    help="hillshade illumination altitude, degrees")
    ap.add_argument("--repeats", type=int, default=1,
                    help="timed repetitions after a discarded warm-up "
                         "(surface routines only). The reported `seconds` is "
                         "the best of these; `seconds_median` the median.")
    args = ap.parse_args()

    try:
        import topotoolbox  # noqa: F401
    except ImportError:
        # Exit 77 is CTest's skip code; the driver propagates it. An optional
        # reference implementation must never turn CI red.
        print("SKIP: topotoolbox not importable. Install via "
              "`pixi run -e benchmark ...` (PyPI, pinned in the reference "
              "feature of pixi.toml).", file=sys.stderr)
        return 77

    info = run(args.routine, args.in_path, args.out_path, args.repeats,
               args.azimuth, args.altitude)
    print(json.dumps(info))
    return 0


if __name__ == "__main__":
    sys.exit(main())
