#!/usr/bin/env python3
"""Append the NCEP2 atmosphere climatology to an existing crn_parameters.nc.

The CAIRN / OCTOPUS cosmogenic calculators use the NCEP2 reanalysis atmosphere
(Balco's NCEPatm_2). lsdtt3's default is ERA40; adding the NCEP grids lets the
cosmo tools run with `atmosphere: NCEP` for a like-for-like comparison.

Source: the legacy binary `NCEP2.bin` shipped with LSDTopoTools_CRNBasinwide
(little-endian float64, column-major lon-outer/lat-inner):
    slp[145x73]  +  t1000[145x73]  +  lat[73]  +  lon[145]   = 21388 doubles.
`slp` is mean sea-level pressure (hPa); `t1000` is mean 1000 hPa temperature
(degC); lat runs 90 -> -90, lon 0 -> 360 (2.5 deg grid, with a wrap column).

The companion `NCEP_hgt.bin` (geopotential heights) is DELIBERATELY NOT used:
the legacy NCEPatm_2 pressure calculation reads only slp + t1000, and the loader
overwrites the geometric-height array with the geopotential one, so the entire
1.35 MB height file is dead weight. We store only what the calculation needs.

Encoding: float32 + zlib-9 (shuffle). Lossless to ~6e-8 relative vs the float64
source (NCEP data is ~5 significant figures), and ~30x smaller than the raw .bin.

Adds variables: ncep_lat (73), ncep_lon (145), ncep_meanP (73x145, hPa),
ncep_meanT (73x145, degC). It does NOT touch any existing variable, so the ERA40
/ Muscheler data (and the numerics they drive) are unchanged.

Usage:
    python3 python_helpers/add_ncep_to_crn_nc.py \
        --ncep2 /path/to/NCEP2.bin \
        --nc data/crn_parameters.nc
"""
import argparse
import sys
import numpy as np
import netCDF4 as nc

NR, NC = 73, 145  # lat, lon


def decode_ncep2(path):
    a = np.fromfile(path, dtype="<f8")
    expect = NC * NR * 2 + NR + NC
    if a.size != expect:
        raise ValueError(f"{path}: {a.size} doubles, expected {expect} "
                         f"(is this the right NCEP2.bin?)")
    slp = a[0:NC * NR].reshape(NC, NR).T            # -> (lat, lon)
    t1000 = a[NC * NR:2 * NC * NR].reshape(NC, NR).T
    lat = a[2 * NC * NR:2 * NC * NR + NR]
    lon = a[2 * NC * NR + NR:2 * NC * NR + NR + NC]
    # sanity
    if not (900 < slp.min() < slp.max() < 1100):
        raise ValueError(f"slp out of range ({slp.min()}..{slp.max()}) - wrong file/layout?")
    return lat, lon, slp, t1000


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ncep2", required=True, help="path to legacy NCEP2.bin")
    ap.add_argument("--nc", default="data/crn_parameters.nc",
                    help="crn_parameters.nc to append to (in place)")
    args = ap.parse_args()

    lat, lon, slp, t1000 = decode_ncep2(args.ncep2)
    print(f"decoded NCEP2: slp {slp.min():.1f}-{slp.max():.1f} hPa, "
          f"t1000 {t1000.min():.1f}-{t1000.max():.1f} degC, grid {NR}x{NC}")

    ds = nc.Dataset(args.nc, "a")
    try:
        if "ncep_meanP" in ds.variables:
            print("ncep_meanP already present; overwriting values in place.")
        else:
            ds.createDimension("ncep_lat", NR)
            ds.createDimension("ncep_lon", NC)
        kw = dict(zlib=True, complevel=9, shuffle=True)

        def put(name, dims, arr, units, longname):
            if name not in ds.variables:
                v = ds.createVariable(name, "f4", dims, **kw)
            else:
                v = ds.variables[name]
            v[:] = arr.astype("f4")
            v.units = units
            v.long_name = longname
            v.source = "NCEP2.bin (LSDTopoTools_CRNBasinwide); Balco NCEPatm_2 atmosphere"

        put("ncep_lat", ("ncep_lat",), lat, "degrees_north", "NCEP2 latitude (90 to -90)")
        put("ncep_lon", ("ncep_lon",), lon, "degrees_east", "NCEP2 longitude (0 to 360)")
        put("ncep_meanP", ("ncep_lat", "ncep_lon"), slp, "hPa",
            "NCEP2 mean sea-level pressure")
        put("ncep_meanT", ("ncep_lat", "ncep_lon"), t1000, "degC",
            "NCEP2 mean 1000 hPa temperature")
    finally:
        ds.close()
    print(f"appended ncep_lat/ncep_lon/ncep_meanP/ncep_meanT to {args.nc}")


if __name__ == "__main__":
    sys.exit(main())
