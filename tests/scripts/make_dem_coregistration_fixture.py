#!/usr/bin/env python3
# make_dem_coregistration_fixture.py
#
# Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
#
# Copyright (C) 2026 Simon M. Mudd
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# -----------------------------------------------------------------------------
# Writes the DEM pair used by the lsdtt-raster-calculator coregistration test.
#
# Nuth and Kaab coregistration needs two overlapping DEMs that differ by a
# horizontal shift, on terrain whose aspects cover the compass. Both grids here
# sample one analytic surface on the same 120 x 120 grid of 10 m cells:
#
#   coreg_reference.asc   z = f(x, y)
#   coreg_shifted.asc     z = f(x - DX, y - DY) + DZ
#
# with f a smooth sum of a long-wavelength sinusoid, a Gaussian hill and a
# diagonal ripple, so every aspect is present and slopes are moderate. The
# known shift (DX, DY, DZ) = (4, -3, 0.5) m is what the coregistration should
# recover. Standard library only, no randomness.
#
# Output: <out_dir>/coreg_reference.asc|.prj and coreg_shifted.asc|.prj
# (WGS 84 / UTM 30N).
#
#   python3 tests/scripts/make_dem_coregistration_fixture.py <out_dir>
# -----------------------------------------------------------------------------

import math
import sys
from pathlib import Path

N_ROWS = 120
N_COLS = 120
CELL = 10.0
XLL = 400000.0
YLL = 6200000.0
DX, DY, DZ = 4.0, -3.0, 0.5

UTM30N_WKT = (
    'PROJCS["WGS 84 / UTM zone 30N",GEOGCS["WGS 84",DATUM["WGS_1984",'
    'SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],'
    'UNIT["degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],'
    'PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-3],'
    'PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],'
    'PARAMETER["false_northing",0],UNIT["metre",1]]'
)


def f(x, y):
    """Surface height at grid-relative coordinates (metres from the SW corner)."""
    return (200.0
            + 40.0 * math.sin(x / 170.0) * math.cos(y / 230.0)
            + 25.0 * math.exp(-((x - 600.0) ** 2 + (y - 500.0) ** 2) / (2.0 * 150.0 ** 2))
            + 10.0 * math.sin((x + y) / 90.0))


def write_grid(path: Path, dx: float, dy: float, dz: float):
    lines = [
        f"ncols {N_COLS}", f"nrows {N_ROWS}",
        f"xllcorner {XLL:.1f}", f"yllcorner {YLL:.1f}",
        f"cellsize {CELL:g}", "NODATA_value -9999",
    ]
    for r in range(N_ROWS):
        y = (N_ROWS - 1 - r + 0.5) * CELL   # cell-centre height above the south edge
        lines.append(" ".join(
            f"{f((c + 0.5) * CELL - dx, y - dy) + dz:.4f}" for c in range(N_COLS)))
    path.write_text("\n".join(lines) + "\n")
    path.with_suffix(".prj").write_text(UTM30N_WKT + "\n")


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: make_dem_coregistration_fixture.py <out_dir>")
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    write_grid(out / "coreg_reference.asc", 0.0, 0.0, 0.0)
    write_grid(out / "coreg_shifted.asc", DX, DY, DZ)
    print(f"wrote coreg_reference.asc and coreg_shifted.asc in {out}")


if __name__ == "__main__":
    main()
