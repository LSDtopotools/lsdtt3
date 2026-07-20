#!/usr/bin/env python3
"""RichDEM reference implementations for the flow-routing benchmark.

docs/testing_and_benchmarking.md.

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
import sys
import time

ROUTINES = ("fill", "fill_epsilon", "breach",
            "accum_d8", "accum_dinf", "accum_mfd")


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


def run(routine: str, in_path: str, out_path: str, mfd_exponent: float):
    import numpy as np
    import richdem as rd

    with suppressed_output():
        dem = rd.LoadGDAL(in_path)

    t0 = time.perf_counter()
    with suppressed_output():
        if routine == "fill":
            result = rd.FillDepressions(dem, epsilon=False, in_place=False)
        elif routine == "fill_epsilon":
            result = rd.FillDepressions(dem, epsilon=True, in_place=False)
        elif routine == "breach":
            result = rd.BreachDepressions(dem, in_place=False)
        elif routine == "accum_d8":
            result = rd.FlowAccumulation(dem, method="D8")
        elif routine == "accum_dinf":
            result = rd.FlowAccumulation(dem, method="Dinf")
        elif routine == "accum_mfd":
            result = rd.FlowAccumulation(dem, method="Freeman",
                                         exponent=mfd_exponent)
        else:
            raise ValueError(f"unknown routine {routine!r}")
    elapsed = time.perf_counter() - t0

    with suppressed_output():
        rd.SaveGDAL(out_path, result)

    arr = np.asarray(result)
    finite = arr[np.isfinite(arr)]
    return {
        "routine": routine,
        "seconds": elapsed,
        "shape": list(arr.shape),
        "min": float(finite.min()) if finite.size else None,
        "max": float(finite.max()) if finite.size else None,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--routine", required=True, choices=ROUTINES)
    ap.add_argument("--in", dest="in_path", required=True)
    ap.add_argument("--out", dest="out_path", required=True)
    ap.add_argument("--mfd-exponent", type=float, default=1.0)
    args = ap.parse_args()

    try:
        import richdem  # noqa: F401
    except ImportError:
        # Exit 77 is CTest's skip code; the driver propagates it.
        print("SKIP: richdem not importable. Install via "
              "`pixi run -e benchmark ...` (conda-forge build).",
              file=sys.stderr)
        return 77

    info = run(args.routine, args.in_path, args.out_path, args.mfd_exponent)
    # Single JSON line on stdout: the driver parses this, and suppressed_output()
    # guarantees RichDEM has not polluted the stream.
    print(json.dumps(info))
    return 0


if __name__ == "__main__":
    sys.exit(main())
