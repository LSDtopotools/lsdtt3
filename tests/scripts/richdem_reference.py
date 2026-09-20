#!/usr/bin/env python3
"""RichDEM reference implementations for the flow-routing benchmark.

docs/developer_notes/testing_and_benchmarking.md.

Runs one RichDEM routine over a GeoTIFF and writes the result as a GeoTIFF.
Kept as a standalone script, invoked as a subprocess by run_flow_benchmark.py,
for two reasons:

  1. RichDEM is an optional dependency.  Isolating it here means the driver can
     detect its absence and skip cleanly rather than failing to import.
  2. It lets the reference run under a different interpreter than the driver
     (`--python`), which matters because RichDEM from PyPI pins numpy<2.  Use
     the conda-forge build via pixi instead: `pixi run -e benchmark ...`.

Results are cached by the caller; RichDEM is deterministic, and the 2400^2
fixture is slow.

Routine -> RichDEM mapping:
    fill              rd.FillDepressions(epsilon=False)   Barnes 2014
    fill_epsilon      rd.FillDepressions(epsilon=True)    Barnes 2014 + epsilon
    breach            rd.BreachDepressions()              Lindsay 2016
    accum_d8          rd.FlowAccumulation(method='D8')
    accum_dinf        rd.FlowAccumulation(method='Dinf')  Tarboton 1997
    accum_mfd         rd.FlowAccumulation(method='Freeman', exponent=...)

Usage:
    python richdem_reference.py --routine fill --in dem.tif --out ref.tif
"""

import argparse
import contextlib
import io
import json
import os
import resource
import statistics
import sys
import time

FLOW_ROUTINES = ("fill", "fill_epsilon", "breach",
                 "accum_d8", "accum_dinf", "accum_mfd")

# Surface metrics -- TIMING ONLY.
#
# These went missing from the benchmark for a bad reason. The case registry in
# run_benchmark_suite.py carried a comment asserting "RichDEM has no slope
# tool", and on the strength of that comment RichDEM was left out of the
# surface suite entirely. It is wrong: rd.TerrainAttribute in 2.3.0 offers
# slope_riserun, slope_percentage, slope_degrees, slope_radians, aspect,
# curvature, planform_curvature and profile_curvature.
#
# WHAT THESE ACTUALLY COMPUTE, because the four engines are NOT computing the
# same quantity and the rows must not be read as if they were:
#
#   RichDEM      Zevenbergen & Thorne (1987): a 9-term quadratic fitted
#                analytically to the 3x3 neighbourhood, curvatures from its
#                partial derivatives. Fixed 3x3 -- there is no window size.
#   lsdtt3       a 6-term quadratic least-squares fit over a moving window of
#                polyfit_lambda_metres, every metric from the ONE fit.
#   TopoToolbox  scipy convolutions with the Evans 3x3 kernels, in Python.
#   Whitebox     3x3 finite-difference operators in Rust.
#
# RichDEM is the closest of the three references to lsdtt3 in FORM -- both fit
# a quadratic and differentiate it -- and still not comparable in VALUE,
# because RichDEM's is an exact 3x3 interpolation where lsdtt3's is a
# least-squares fit over a window the user chooses. Sign and handedness
# conventions for the curvatures were not checked and are not assumed to
# agree. Every surface case is report_only; no value is compared.
SURFACE_ROUTINES = (
    "slope",
    "aspect",
    "curvature_profc",
    "curvature_planc",
    "curvature_all2",   # profile + planform back to back, one shared timing
)

ROUTINES = FLOW_ROUTINES + SURFACE_ROUTINES


@contextlib.contextmanager
def suppressed_output():
    """Silence RichDEM's progress bars.

    RichDEM writes multi-kilobyte ANSI progress bars straight to the C-level
    stdout/stderr, which swamps test logs.  Redirect the file descriptors, not
    just sys.stdout, or the C++ layer writes past the Python-level capture.
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


def peak_rss_mb() -> float:
    """Peak resident set size of this process, in MiB.

    ru_maxrss is kibibytes on Linux and bytes on macOS; the Linux reading is
    used because that is where the benchmark runs.  Matches the identical
    helper in topotoolbox_reference.py and whitebox_reference.py so the three
    references report memory on the same basis.
    """
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def _call_for(routine: str, dem, mfd_exponent: float):
    """Zero-argument callable for one routine, built outside the timing loop."""
    import richdem as rd

    if routine == "fill":
        return lambda: rd.FillDepressions(dem, epsilon=False, in_place=False)
    if routine == "fill_epsilon":
        return lambda: rd.FillDepressions(dem, epsilon=True, in_place=False)
    if routine == "breach":
        return lambda: rd.BreachDepressions(dem, in_place=False)
    if routine == "accum_d8":
        return lambda: rd.FlowAccumulation(dem, method="D8")
    if routine == "accum_dinf":
        return lambda: rd.FlowAccumulation(dem, method="Dinf")
    if routine == "accum_mfd":
        return lambda: rd.FlowAccumulation(dem, method="Freeman",
                                           exponent=mfd_exponent)

    # -- surface -----------------------------------------------------------
    # Degrees for slope, matching whitebox_reference.py's units="degrees", so
    # the two references at least agree with each other on what they return.
    # The unit does not affect cost.
    if routine == "slope":
        return lambda: rd.TerrainAttribute(dem, attrib="slope_degrees")
    if routine == "aspect":
        return lambda: rd.TerrainAttribute(dem, attrib="aspect")
    if routine == "curvature_profc":
        return lambda: rd.TerrainAttribute(dem, attrib="profile_curvature")
    if routine == "curvature_planc":
        return lambda: rd.TerrainAttribute(dem, attrib="planform_curvature")
    if routine == "curvature_all2":
        # What a user pays to get both curvatures. RichDEM, like Whitebox and
        # TopoToolbox, refits per curvature type; one lsdtt3 polyfit yields
        # every coefficient at once. Two calls rather than four, mirroring
        # whitebox_reference.py's curvature_all2: RichDEM's `curvature` is
        # total curvature, not a tangential/mean pair, so there is no
        # four-metric equivalent to compare against TopoToolbox's
        # curvature_all4.
        return lambda: [rd.TerrainAttribute(dem, attrib="profile_curvature"),
                        rd.TerrainAttribute(dem, attrib="planform_curvature")]

    raise ValueError(f"unknown routine {routine!r}")


def run(routine: str, in_path: str, out_path, mfd_exponent: float,
        repeats: int = 1):
    """Time one RichDEM routine.

    repeats=1 is the historical behaviour EXACTLY: one timed call, no warm-up,
    so run_flow_benchmark.py's numbers and cached JSON are unchanged.  Only
    when repeats > 1 is a warm-up discarded and a best-of-N taken -- which is
    what run_benchmark_suite.py asks for, so that RichDEM rows are produced on
    the same statistical basis as the Whitebox and TopoToolbox rows rather than
    being single-shot against their best-of-N.

    I/O is deliberately outside the timed region: the load happens before the
    loop and the save after it, so `seconds` is the algorithm.
    """
    import numpy as np

    t0 = time.perf_counter()
    with suppressed_output():
        dem = rd_load(in_path)
    seconds_read = time.perf_counter() - t0

    call = _call_for(routine, dem, mfd_exponent)

    if repeats > 1:
        with suppressed_output():
            call()                       # warm-up, discarded

    timings = []
    for _ in range(max(1, repeats)):
        t0 = time.perf_counter()
        with suppressed_output():
            result = call()
        timings.append(time.perf_counter() - t0)

    # A list result (curvature_all2) has no single raster to stand for it;
    # report the first and write the first, as whitebox_reference.py does.
    # Stacking them instead would put a leading axis into `shape` and make the
    # row look like a different-sized problem.
    primary = result[0] if isinstance(result, list) else result

    seconds_write = None
    if out_path:
        import richdem as rd
        t0 = time.perf_counter()
        with suppressed_output():
            rd.SaveGDAL(out_path, primary)
        seconds_write = time.perf_counter() - t0

    arr = np.asarray(primary)
    finite = arr[np.isfinite(arr)]
    return {
        "reference": "richdem",
        "routine": routine,
        "seconds": min(timings),
        "seconds_median": statistics.median(timings),
        "seconds_all": timings,
        "seconds_stdev": statistics.stdev(timings) if len(timings) > 1 else 0.0,
        "repeats": len(timings),
        # RichDEM's FillDepressions/BreachDepressions/FlowAccumulation and
        # TerrainAttribute are all serial, so the thread count is 1 whatever
        # the caller asked for. Stated rather than omitted so a campaign cannot
        # silently attribute a thread count to it.
        "threads": 1,
        "seconds_read": seconds_read,
        "seconds_write": seconds_write,
        "peak_rss_mb": peak_rss_mb(),
        "shape": list(arr.shape),
        "min": float(finite.min()) if finite.size else None,
        "max": float(finite.max()) if finite.size else None,
    }


def rd_load(in_path: str):
    import richdem as rd
    return rd.LoadGDAL(in_path)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--routine", required=True, choices=ROUTINES)
    ap.add_argument("--in", dest="in_path", required=True)
    ap.add_argument("--out", dest="out_path", default=None,
                    help="optional output raster. Omit for a timing-only run: "
                         "run_benchmark_suite.py never compares reference "
                         "values, so writing them would be pure I/O cost.")
    ap.add_argument("--mfd-exponent", type=float, default=1.0)
    ap.add_argument("--repeats", type=int, default=1,
                    help="timed repeats. 1 (the default) reproduces the "
                         "historical single-shot behaviour exactly; >1 adds a "
                         "discarded warm-up and reports best-of-N.")
    args = ap.parse_args()

    try:
        import richdem  # noqa: F401
    except ImportError:
        # Exit 77 is CTest's skip code; the driver propagates it.
        print("SKIP: richdem not importable. Install via "
              "`pixi run -e benchmark ...` (conda-forge build).",
              file=sys.stderr)
        return 77

    info = run(args.routine, args.in_path, args.out_path,
               args.mfd_exponent, args.repeats)
    # Single JSON line on stdout: the driver parses this, and suppressed_output()
    # guarantees RichDEM has not polluted the stream.
    print(json.dumps(info))
    return 0


if __name__ == "__main__":
    sys.exit(main())
