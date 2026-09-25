#!/usr/bin/env python3
"""Check that lsdtt3 binaries launched by the test harness use their OWN GDAL.

Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3

Copyright (C) 2026 Simon M. Mudd
SPDX-License-Identifier: AGPL-3.0-or-later

The test scripts start lsdtt3 CLIs through lsdtt_test_env. When the tests run
inside Pixi/conda (`pixi run -e test ctest ...`), that environment's GDAL_DATA,
PROJ_DATA and GDAL_DRIVER_PATH belong to ITS GDAL, which is not the one a
system-built CLI links. lsdtt_test_env strips them for lsdtt3 children in that
case; this test proves it, using an existing CLI's `--version`, which reports
the GDAL/PROJ versions loaded at run time and where their data came from.

Checks:
  1. the CLI's run-time GDAL version is the one the build compiled against
     (--expect-gdal) -- in any environment;
  2. inside a Pixi/conda environment whose GDAL is not the build's: the CLI's
     PROJ and GDAL data do not come from that environment;
  3. and, as a control, the same CLI launched WITHOUT the launcher's knowledge
     of the build's GDAL does pick the environment's data up -- so check 2 is
     not passing vacuously;
  4. the parent Python keeps its own GDAL/PROJ working (when osgeo is
     importable): the environment is stripped for the child only.

Outside any such environment, only check 1 applies, and the test says so.

Usage: check_test_env_isolation.py --exe <lsdtt CLI> --expect-gdal <x.y.z>
"""

from __future__ import annotations

import argparse
import os
import re
import sys

import lsdtt_test_env


def version_report(exe: str, env: dict[str, str] | None = None) -> dict[str, str]:
    kwargs = {"capture_output": True, "text": True, "check": False}
    if env is not None:
        kwargs["env"] = env
    result = lsdtt_test_env.run([exe, "--version"], **kwargs)
    text = result.stdout + result.stderr
    report = {"text": text}
    m = re.search(r"GDAL (\S+), PROJ (\S+)", text)
    if m:
        report["gdal"], report["proj"] = m.group(1), m.group(2)
    for key, label in (("proj_data", "proj data"), ("gdal_data", "GDAL data")):
        m = re.search(rf"{label}:\s*(.+)", text)
        if m:
            report[key] = m.group(1).strip()
    return report


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", required=True, help="an lsdtt3 CLI")
    ap.add_argument("--expect-gdal", required=True,
                    help="the GDAL version the build compiled against")
    args = ap.parse_args()

    failures: list[str] = []

    def check(name: str, ok: bool, detail: str = "") -> None:
        print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"   [{detail}]" if detail else ""))
        if not ok:
            failures.append(name)

    child = version_report(args.exe)
    check("the CLI reports its run-time GDAL version", "gdal" in child,
          child["text"].strip().splitlines()[-1] if child["text"].strip() else "no output")
    if "gdal" in child:
        check("run-time GDAL is the version the build compiled against",
              child["gdal"] == args.expect_gdal,
              f"run time {child['gdal']}, compiled {args.expect_gdal}")

    env = dict(os.environ)
    if not lsdtt_test_env.build_gdal_is_foreign_to_env(env):
        why = ("no Pixi/conda environment is active" if not env.get("CONDA_PREFIX")
               else "the build's GDAL is the active environment's own")
        print(f"NOTE  {why}; the isolation checks do not apply")
    else:
        prefix = os.path.realpath(env["CONDA_PREFIX"])
        print(f"NOTE  inside {prefix}; build GDAL {env[lsdtt_test_env.GDAL_LIBRARY_ENV]}")

        for key in ("proj_data", "gdal_data"):
            value = child.get(key, "")
            check(f"the CLI's {key.replace('_', ' ')} is not the environment's",
                  prefix not in os.path.realpath(value) if value.startswith("/")
                  else prefix not in value,
                  value or "not reported")

        # Control: without LSDTT_TEST_GDAL_LIBRARY the launcher cannot know the
        # build's GDAL is foreign, so it leaves the environment alone.
        control_env = dict(env)
        control_env.pop(lsdtt_test_env.GDAL_LIBRARY_ENV, None)
        control = version_report(args.exe, control_env)
        leaked = any(prefix in control.get(k, "") for k in ("proj_data", "gdal_data"))
        check("control: unstripped, the same CLI does see the environment's data",
              leaked, control.get("proj_data", "not reported"))

        try:
            from osgeo import osr
        except ImportError:
            print("NOTE  osgeo not importable in the parent; parent check skipped")
        else:
            srs = osr.SpatialReference()
            check("the parent's own GDAL/PROJ still work (EPSG:4326)",
                  srs.ImportFromEPSG(4326) == 0)

    print(f"RESULT: {'FAIL' if failures else 'PASS'}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
