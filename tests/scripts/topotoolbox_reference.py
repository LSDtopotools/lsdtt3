#!/usr/bin/env python3
"""TopoToolbox reference implementations for the flow-routing benchmark.

Sibling of `richdem_reference.py`, and for flow routing the more trustworthy of
the two.

WHY BOTH REFERENCES

RichDEM's D8 selects the neighbour with the LOWEST ELEVATION with a
diagonal-preferring tie-break, not the steepest DESCENT (drop / distance). On a
planar ramp with a known analytic answer it routes diagonally and gets the
accumulation wrong. It is therefore kept as a *report-only* reference for
accumulation.

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
    fill      GridObject.fillsinks()            plain fill, no epsilon
    accum_d8  FlowObject(method='d8').flow_accumulation()

Usage:
    python topotoolbox_reference.py --routine accum_d8 --in dem.tif --out ref.tif
"""

import argparse
import contextlib
import json
import os
import sys
import time

ROUTINES = ("fill", "accum_d8")


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


def run(routine: str, in_path: str, out_path: str):
    import numpy as np
    import topotoolbox as ttb

    with suppressed_output():
        dem = ttb.read_tif(in_path)

    info = {"routine": routine}

    if routine == "fill":
        t0 = time.perf_counter()
        with suppressed_output():
            result = dem.fillsinks()
        info["seconds"] = time.perf_counter() - t0

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
        t0 = time.perf_counter()
        with suppressed_output():
            flow = ttb.FlowObject(dem, method="d8", sink_resolution="carve")
        build_s = time.perf_counter() - t0

        t1 = time.perf_counter()
        with suppressed_output():
            result = flow.flow_accumulation()
        accum_s = time.perf_counter() - t1

        info["seconds"] = build_s + accum_s
        info["seconds_build_only"] = build_s
        info["seconds_accumulate_only"] = accum_s

    else:
        raise ValueError(f"unknown routine {routine!r}")

    with suppressed_output():
        ttb.write_tif(result, out_path)

    arr = np.asarray(result.z)
    finite = arr[np.isfinite(arr)]
    info.update({
        "shape": list(arr.shape),
        "min": float(finite.min()) if finite.size else None,
        "max": float(finite.max()) if finite.size else None,
    })
    return info


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--routine", required=True, choices=ROUTINES)
    ap.add_argument("--in", dest="in_path", required=True)
    ap.add_argument("--out", dest="out_path", required=True)
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

    info = run(args.routine, args.in_path, args.out_path)
    print(json.dumps(info))
    return 0


if __name__ == "__main__":
    sys.exit(main())
