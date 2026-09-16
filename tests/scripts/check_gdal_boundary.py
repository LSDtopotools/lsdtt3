#!/usr/bin/env python3
# check_gdal_boundary.py
#
# Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
#
# Copyright (C) 2026 Simon M. Mudd
# SPDX-License-Identifier: AGPL-3.0-or-later
# Released under the GNU Affero General Public Licence v3 or later
# (AGPL-3.0-or-later). See LICENSE.
#
# -----------------------------------------------------------------------------
# Enforces the GDAL facade established in Phase 1 of the GDAL/PROJ friction work.
#
# Two rules:
#
#   1. No lsdtt3 HEADER may include a GDAL, OGR or CPL header.
#      Headers are included by hundreds of translation units. Once <gdal.h> is
#      reachable from one of them it is reachable from most of the tree, which
#      is the state Phase 1 undid: compile times, and any hope of swapping the
#      backend later, both depend on this staying true.
#
#      The one exception is the bridge header lsdtt_pixel_type_gdal.hpp, whose
#      whole purpose is to name GDALDataType, and which only .cpp files at the
#      boundary may include.
#
#   2. No file in cli/ may include a GDAL, OGR or CPL header.
#      A CLI's job is argument parsing and orchestration. When a CLI needs GDAL
#      it means a library function is missing; add it to src/core/ instead.
#
# A .cpp under src/ may include GDAL freely -- that is where the boundary lives.
# -----------------------------------------------------------------------------

import pathlib
import re
import sys

GDAL_INCLUDE = re.compile(r'^\s*#\s*include\s*[<"](gdal|ogr|cpl)[A-Za-z0-9_/.]*[>"]',
                          re.MULTILINE)

# The bridge header exists precisely to name GDAL types; see its own comment.
HEADER_ALLOWLIST = {"lsdtt_pixel_type_gdal.hpp"}

HEADER_SUFFIXES = {".hpp", ".h", ".cuh"}


def main() -> int:
    repo = pathlib.Path(__file__).resolve().parents[2]

    violations = []

    # Rule 1: headers anywhere under src/
    for path in sorted((repo / "src").rglob("*")):
        if path.suffix not in HEADER_SUFFIXES or not path.is_file():
            continue
        if path.name in HEADER_ALLOWLIST:
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        for m in GDAL_INCLUDE.finditer(text):
            line = text[: m.start()].count("\n") + 1
            violations.append(
                f"{path.relative_to(repo)}:{line}: header includes "
                f"{m.group(0).strip()}"
            )

    # Rule 2: anything under cli/
    for path in sorted((repo / "cli").rglob("*.cpp")):
        text = path.read_text(encoding="utf-8", errors="replace")
        for m in GDAL_INCLUDE.finditer(text):
            line = text[: m.start()].count("\n") + 1
            violations.append(
                f"{path.relative_to(repo)}:{line}: CLI includes "
                f"{m.group(0).strip()}"
            )

    if violations:
        print("GDAL boundary violations:\n", file=sys.stderr)
        for v in violations:
            print(f"    {v}", file=sys.stderr)
        print(
            "\nlsdtt3 keeps GDAL/OGR/CPL confined to .cpp files under src/.\n"
            "  - A header that needs a GDAL type should use lsdtt_pixel_type,\n"
            "    an opaque void*, or a forward declaration (see\n"
            "    src/core/lsdtt_data_loader.hpp for all three).\n"
            "  - A CLI that needs GDAL is missing a library function; add it\n"
            "    to src/core/ and call that instead.\n",
            file=sys.stderr,
        )
        return 1

    print("GDAL boundary clean: no GDAL/OGR/CPL includes in headers or cli/.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
