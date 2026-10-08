#!/usr/bin/env python3
# make_marsh_topography_fixture.py
#
# Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
#
# Copyright (C) 2026 Simon M. Mudd
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# -----------------------------------------------------------------------------
# Writes the synthetic salt-marsh DEM used by the lsdtt-marsh-topography tests.
#
# lsdtt-synthetic-surfaces has no marsh surface, and a marsh test needs the three
# things the workflow looks for: a high, nearly flat platform; a scarp down to a
# low tidal flat; and channels cut into the platform that drain across the scarp.
# Standard library only, fixed seed, so the file is the same on every machine and
# every Python 3 (random.Random(seed).random() is stable across versions).
#
#   rows   0-59   platform, 1.60 m falling 0.002 m/m to the south, +-0.03 m noise
#   rows  60-63   scarp, linear drop to the flat
#   rows  64-99   tidal flat, 0.30 m falling to 0.10 m at the south edge,
#                 +-0.02 m noise
#   channels      trunk down column 60 from row 5 across the scarp, 3 px wide and
#                 0.8 m deep; two tributaries, 2 px wide and 0.6 m deep, joining
#                 the trunk at rows 40 and 35; parabolic cross profile
#
# Output: <out_dir>/marsh_fixture_dem.asc and .prj (WGS 84 / UTM 30N, 1 m cells).
#
#   python3 tests/scripts/make_marsh_topography_fixture.py <out_dir>
# -----------------------------------------------------------------------------

import math
import random
import sys
from pathlib import Path

N_ROWS = 100
N_COLS = 120
CELL = 1.0
XLL = 400000.0
YLL = 6200000.0
SEED = 7

PLATFORM_Z = 1.60
PLATFORM_FALL = 0.002
SCARP_TOP_ROW = 60
FLAT_TOP_ROW = 64
FLAT_Z_TOP = 0.30
FLAT_Z_BOTTOM = 0.10

# (polyline in (row, col), half-width in pixels, depth in metres)
CHANNELS = [
    ([(5, 60), (70, 60)], 1.5, 0.8),
    ([(10, 25), (25, 45), (40, 60)], 1.0, 0.6),
    ([(15, 100), (35, 60)], 1.0, 0.6),
]

UTM30N_WKT = (
    'PROJCS["WGS 84 / UTM zone 30N",GEOGCS["WGS 84",DATUM["WGS_1984",'
    'SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],'
    'UNIT["degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],'
    'PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-3],'
    'PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],'
    'PARAMETER["false_northing",0],UNIT["metre",1]]'
)


def distance_to_segment(r, c, a, b):
    (r0, c0), (r1, c1) = a, b
    dr, dc = r1 - r0, c1 - c0
    L2 = dr * dr + dc * dc
    t = 0.0 if L2 == 0 else max(0.0, min(1.0, ((r - r0) * dr + (c - c0) * dc) / L2))
    return math.hypot(r - (r0 + t * dr), c - (c0 + t * dc))


def surface(r, c, rng):
    if r < SCARP_TOP_ROW:
        z = PLATFORM_Z - PLATFORM_FALL * r * CELL + rng.uniform(-0.03, 0.03)
    elif r < FLAT_TOP_ROW:
        top = PLATFORM_Z - PLATFORM_FALL * SCARP_TOP_ROW * CELL
        f = (r - SCARP_TOP_ROW + 1) / (FLAT_TOP_ROW - SCARP_TOP_ROW + 1)
        z = top + f * (FLAT_Z_TOP - top)
    else:
        f = (r - FLAT_TOP_ROW) / (N_ROWS - 1 - FLAT_TOP_ROW)
        z = FLAT_Z_TOP + f * (FLAT_Z_BOTTOM - FLAT_Z_TOP) + rng.uniform(-0.02, 0.02)
    cut = 0.0
    for line, half_width, depth in CHANNELS:
        d = min(distance_to_segment(r, c, line[i], line[i + 1])
                for i in range(len(line) - 1))
        if d < half_width:
            cut = max(cut, depth * (1.0 - (d / half_width) ** 2))
    # A channel never cuts below the flat it drains to.
    return max(z - cut, FLAT_Z_BOTTOM - 0.05) if cut > 0 else z


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: make_marsh_topography_fixture.py <out_dir>")
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED)
    lines = [
        f"ncols {N_COLS}", f"nrows {N_ROWS}",
        f"xllcorner {XLL:.1f}", f"yllcorner {YLL:.1f}",
        f"cellsize {CELL:g}", "NODATA_value -9999",
    ]
    for r in range(N_ROWS):
        lines.append(" ".join(f"{surface(r, c, rng):.4f}" for c in range(N_COLS)))
    (out / "marsh_fixture_dem.asc").write_text("\n".join(lines) + "\n")
    (out / "marsh_fixture_dem.prj").write_text(UTM30N_WKT + "\n")
    print(f"wrote {out / 'marsh_fixture_dem.asc'} ({N_ROWS} x {N_COLS})")


if __name__ == "__main__":
    main()
