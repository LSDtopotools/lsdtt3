#!/usr/bin/env python3
"""
build_crn_parameters_nc.py

Assemble `crn_parameters.nc` for LSDTopoTools3's cosmogenic-nuclide tools
(`lsdtt-cosmo-tool`, `lsdtt-cosmo-prediction`). This file supplies the atmosphere
grids, the Lifton lapse-rate polynomial, and the geomagnetic virtual dipole moment
(VDM) record needed for the ERA40-style atmosphere and the time-integrated
palaeomagnetic (time-varying production) correction. Any user who wants the
`crn_nc_fname` + `sample_age_yr` options needs this file; it is not test-only.

It writes to the EXACT schema read by
`lsdtt::lsdtt_crn_parameters::load_nc_data` (src/cpu/lsdtt_crn_parameters.cpp).
It is orchestration-only (house rule: Python for data prep / diff, never the numeric
core). The three inputs can be supplied independently, so you can gather the
atmosphere now and add the VDM record / lapse-rate coefficients later.

------------------------------------------------------------------------------
Schema required by load_nc_data (variable NAMES and DTYPES are load-bearing; dim
lengths are read from the file, not hardcoded, canonical grid 73 x 145):

  era40_lat       double  (nlat,)          deg, DESCENDING  (+90 -> -90)
  era40_lon       double  (nlon,)          deg, ASCENDING   (0 -> 360)
  era40_meanP     float   (nlat, nlon)     hPa, SEA-LEVEL mean pressure, row-major [ilat][ilon]
  era40_meanT     float   (nlat, nlon)     K,   1000 hPa (near-sea-level) mean temperature, row-major
  lifton_lr       double  (7,)             Lifton CIRA-86 lapse-rate polynomial coeffs, K m^-1
  muscheler_time  double  (ntime,)         yr BP, ASCENDING
  muscheler_vdm   float   (ntime,)         VDM, units 1e22 A m^2

`loaded` becomes true iff era40_lat, era40_meanP, lifton_lr, muscheler_vdm are all
present and non-empty.

WHAT era40_meanP / era40_meanT ARE (confirmed against the code, which ports
Balco/Lifton's ERA40atm from lsd.py): meanP is SEA-LEVEL pressure; meanT is the
1000 hPa temperature used as the sea-level reference in the barometric formula
P = P_sl*exp((gmr/dtdz)*(ln T_sl - ln T_site)), T_site = T_sl - dtdz*elev. meanT is
NOT 2 m surface temperature.

ATMOSPHERE SOURCE (default `--fetch-era40`): NCEP-DOE Reanalysis 2 (Kanamitsu et al.
2002, doi:10.1175/BAMS-83-11-1631) from NOAA PSL -- public domain, no API key, and
already on the native 2.5 deg 73x144 grid (this script adds the 360 wrap column to
reach 145). This is the faithful lineage of the legacy LSDTopoTools NCEP2 atmosphere.
`mslp.mon.mean.nc` -> meanP (Pascals -> hPa); `air.mon.mean.nc` @ 1000 hPa -> meanT.
Alternatively `--era40 <file>` accepts a pre-built NetCDF (see load_era40_grids).

VDM RECORD (chosen): Knudsen et al. (2008), EPSL 272:319,
doi:10.1016/j.epsl.2008.04.048 (Muscheler is a coauthor; the NetCDF variables keep
the names muscheler_time/_vdm because that is what load_nc_data reads). It is a 1-D
time series (time_yrBP, VDM[1e22 A m^2]); there is no clean DOI'd data file, so it is
supplied via `--vdm <2-col table>` (derive from GEOMAGIA50 / the paper supplement).

LIFTON coefficients: Lifton, Sato & Dunai (2014), EPSL 386:149,
doi:10.1016/j.epsl.2013.10.052 (supplementary MATLAB 'lr' field); supply via
`--lifton <7-number text file>`.

------------------------------------------------------------------------------
REBUILDING THE COMMITTED data/crn_parameters.nc -- THIS SCRIPT IS ONLY STEP 1.

The shipped file also carries the NCEP atmosphere and the LSDn scaling data, both
appended in place afterwards. If you regenerate from scratch you must re-run BOTH
appends or the committed file will lose `atmosphere: NCEP` and `scaling_scheme:
lsdn`:

  1. python_helpers/build_crn_parameters_nc.py  --fetch-era40 --out crn_parameters.nc
  2. python_helpers/add_ncep_to_crn_nc.py       --ncep2 NCEP2.bin --nc crn_parameters.nc
  3. python_helpers/extract_lsdn_reference.py   --data-dir <Balco v3 data/> --out-dir <scratch>
     python_helpers/add_lsdn_to_crn_nc.py       --npz <scratch>/lsdn_reference.npz \
                                                --nc crn_parameters.nc --verify

Steps 2 and 3 are non-destructive, idempotent appends and never modify what this
script wrote. Step 3 needs the Balco v3 reference `.mat` files, which are not in
this repository; see data/README.md for what the LSDn variables are and where they
come from.
------------------------------------------------------------------------------
"""
from __future__ import annotations

import argparse
import os
import sys
import tempfile
import urllib.request

import numpy as np

try:
    from netCDF4 import Dataset, num2date
except ImportError:
    sys.exit("netCDF4 not available -- install into your env (conda-forge netcdf4).")

CANONICAL_NLAT = 73
CANONICAL_NLON = 145
N_LIFTON = 7

# Lifton lapse-rate fit to COSPAR CIRA-86 (<10 km), used by Balco/Lifton ERA40atm:
#   dtdz(lat) = -(Sum_k LIFTON_LR[k] * lat^k),  lat SIGNED, gmr = -0.03417.
# Verbatim from expage ERA40atm.m (Balco/Lifton, Mar 2011), line 72; matches the
# k=0..6 loop in lsdtt_crn_parameters.cpp::atmospheric_pressure_era40. Hardcoded
# because it is a fixed constant of the scaling framework; override with --lifton.
LIFTON_LR_CIRA86 = np.array([
    -6.1517e-03, -3.1831e-06, -1.5014e-07, 1.8097e-09,
    1.1791e-10, -6.5359e-14, -9.5209e-15,
], dtype=np.float64)

# Knudsen et al. (2008) mean VADM record, EPSL 272:319, doi:10.1016/j.epsl.2008.04.048,
# supplementary table mmc1 (columns "Year BP" [before 2000 AD] and "VADM" [1e22 A m^2]).
# 107 points, 250-49500 yr BP ascending. stone2000sp_paleomag reads muscheler_vdm as
# absolute VDM in 1e22 A m^2 and takes index 0 as the modern reference, so this is the
# correct form. (The ~50 yr offset between "BP before 2000 AD" and BP-1950 is negligible.)
# Hardcoded so the file is self-contained; override with --vdm. Cite Knudsen et al. 2008.
KNUDSEN2008_VDM = np.array([
    (250, 8.42), (350, 8.63), (450, 9.01), (550, 9.26), (650, 9.58), (750, 9.85),
    (850, 10.1), (950, 10.28), (1050, 10.51), (1150, 10.68), (1250, 10.85), (1350, 10.88),
    (1450, 10.85), (1550, 10.75), (1650, 10.69), (1750, 10.49), (1850, 10.37), (1950, 10.3),
    (2050, 10.31), (2150, 10.47), (2250, 10.77), (2350, 11.06), (2450, 11.26), (2550, 11.4),
    (2650, 11.47), (2750, 11.24), (2850, 10.98), (2950, 10.9), (3050, 10.7), (3150, 10.26),
    (3250, 10.12), (3350, 9.91), (3450, 9.77), (3550, 9.47), (3650, 9.05), (3750, 8.84),
    (3850, 8.9), (3950, 8.83), (4050, 8.89), (4250, 8.49), (4500, 8.29), (4750, 7.81),
    (5000, 7.61), (5250, 7.32), (5500, 7.19), (5750, 7.1), (6000, 7.01), (6250, 7.01),
    (6500, 7.04), (6750, 7.0), (7000, 7.04), (7250, 7.13), (7500, 7.31), (7750, 7.5),
    (8000, 7.89), (8250, 8.0), (8500, 8.16), (8750, 8.58), (9000, 8.7), (9250, 8.61),
    (9500, 8.9), (9750, 8.19), (10000, 8.2), (10250, 8.44), (10500, 7.85), (10750, 8.17),
    (11000, 8.32), (11250, 7.86), (11500, 7.76), (12500, 6.06), (13500, 5.97), (14500, 6.48),
    (15500, 6.77), (16500, 6.65), (17500, 6.33), (18500, 5.77), (19500, 4.07), (20500, 4.09),
    (21500, 4.65), (22500, 5.22), (23500, 6.8), (24500, 5.69), (25500, 5.31), (26500, 4.78),
    (27500, 5.22), (28500, 5.8), (29500, 5.93), (30500, 6.41), (31500, 6.43), (32500, 5.6),
    (33500, 5.45), (34500, 5.51), (35500, 5.26), (36500, 5.37), (37500, 5.75), (38500, 6.12),
    (39500, 5.3), (40500, 5.16), (41500, 4.65), (42500, 3.45), (43500, 3.58), (44500, 3.5),
    (45500, 3.92), (46500, 3.99), (47500, 4.48), (48500, 4.41), (49500, 4.02),
], dtype=np.float64)

# True present-day (t = 0 yr BP) reference prepended to the Knudsen record so index 0
# is the actual modern field, not the earliest sample (250 yr BP). Knudsen's "yr BP"
# is before 2000 AD, so the matching value is the IGRF-2000 geomagnetic dipole moment,
# 7.79e22 A m^2 (computed from the IGRF-2000 degree-1 Gauss coefficients). This is the
# modern reference stone2000sp_paleomag divides by (sf_modern at index 0).
PRESENT_DAY_VDM_TIME = 0.0
PRESENT_DAY_VDM = 7.79

# NCEP-DOE Reanalysis 2 monthly means at NOAA PSL (public domain, no auth).
NCEP_R2_MSLP_URL = ("https://downloads.psl.noaa.gov/Datasets/"
                    "ncep.reanalysis2/Monthlies/surface/mslp.mon.mean.nc")
NCEP_R2_AIR_URL = ("https://downloads.psl.noaa.gov/Datasets/"
                   "ncep.reanalysis2/Monthlies/pressure/air.mon.mean.nc")


# ---------------------------------------------------------------------------
# Atmosphere: fetch NCEP-DOE R2, compute climatology, reduce to the 73x145 grid.
# ---------------------------------------------------------------------------

def _download(url: str, dest: str) -> bool:
    """Download url to dest unless already present. Returns True if downloaded."""
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        print(f"  using cached {dest} ({os.path.getsize(dest)/1e6:.1f} MB)")
        return False
    print(f"  downloading {url}\n           -> {dest} ...")
    urllib.request.urlretrieve(url, dest)
    print(f"     got {os.path.getsize(dest)/1e6:.1f} MB")
    return True


def _climatology_mean(nc_path, varname, y0, y1, level_hpa=None):
    """Annual climatological mean of varname over calendar years [y0, y1].

    Returns (lat, lon, mean2d, units). If level_hpa is given, first selects that
    vertical level. netCDF4 auto-applies scale_factor/add_offset, so values are
    physical.
    """
    ds = Dataset(nc_path)
    try:
        t = ds.variables["time"]
        dates = num2date(t[:], t.units, only_use_cftime_datetimes=False)
        years = np.array([d.year for d in dates])
        sel = (years >= y0) & (years <= y1)
        if not sel.any():
            raise ValueError(f"{nc_path}: no months in {y0}-{y1} "
                             f"(available {years.min()}-{years.max()})")
        var = ds.variables[varname]
        units = getattr(var, "units", "")
        if level_hpa is not None:
            lev = np.asarray(ds.variables["level"][:], dtype=np.float64)
            li = int(np.argmin(np.abs(lev - level_hpa)))
            if abs(lev[li] - level_hpa) > 1e-6:
                raise ValueError(f"{nc_path}: level {level_hpa} hPa not found; "
                                 f"have {lev}")
            data = np.asarray(var[sel, li, :, :], dtype=np.float64)
        else:
            data = np.asarray(var[sel, :, :], dtype=np.float64)
        mean2d = data.mean(axis=0)
        lat = np.asarray(ds.variables["lat"][:], dtype=np.float64)
        lon = np.asarray(ds.variables["lon"][:], dtype=np.float64)
        print(f"    {varname}: {int(sel.sum())} months "
              f"({years[sel].min()}-{years[sel].max()}), "
              f"grid {lat.size}x{lon.size}, units '{units}'")
        return lat, lon, mean2d, units
    finally:
        ds.close()


def _add_wrap_column(lon, grid):
    """Append the 360-degree wrap column (== 0 deg) so lon spans [0, 360]."""
    lon2 = np.concatenate([lon, [360.0]])
    grid2 = np.concatenate([grid, grid[:, :1]], axis=1)
    return lon2, grid2


def _to_hPa(field, units):
    u = (units or "").lower()
    if "pa" in u and "hpa" not in u and "mb" not in u and "mbar" not in u:
        return field / 100.0     # Pascals -> hPa
    return field                  # already hPa / mbar


def _to_K(field, units):
    u = (units or "").lower()
    if "degc" in u or "celsius" in u:
        return field + 273.15
    return field                  # degK / Kelvin


def fetch_ncep_era40_grids(workdir, y0, y1, keep):
    """Download NCEP-DOE R2, build (lat73, lon145, meanP hPa, meanT K)."""
    os.makedirs(workdir, exist_ok=True)
    mslp_f = os.path.join(workdir, "mslp.mon.mean.nc")
    air_f = os.path.join(workdir, "air.mon.mean.nc")
    got_p = _download(NCEP_R2_MSLP_URL, mslp_f)
    got_t = _download(NCEP_R2_AIR_URL, air_f)

    lat, lon, mslp_mean, punits = _climatology_mean(mslp_f, "mslp", y0, y1)
    latT, lonT, air_mean, tunits = _climatology_mean(air_f, "air", y0, y1, level_hpa=1000.0)
    if not (np.array_equal(lat, latT) and np.array_equal(lon, lonT)):
        raise ValueError("mslp and air grids differ; cannot combine.")

    meanP = _to_hPa(mslp_mean, punits).astype(np.float32)
    meanT = _to_K(air_mean, tunits).astype(np.float32)

    # NCEP native: lat descending (90->-90), lon ascending (0..357.5). Wrap to 145.
    lon145, meanP = _add_wrap_column(lon, meanP)
    _, meanT = _add_wrap_column(lon, meanT)

    if not keep:
        for f, downloaded in ((mslp_f, got_p), (air_f, got_t)):
            if downloaded:
                try:
                    os.remove(f)
                except OSError:
                    pass
        print("  removed large downloads (use --keep-downloads to retain)")

    prov = (f"NCEP-DOE Reanalysis 2 (Kanamitsu et al. 2002), NOAA PSL, "
            f"annual climatology {y0}-{y1}; mslp->meanP, air@1000hPa->meanT")
    return (lat.astype(np.float64), lon145.astype(np.float64), meanP, meanT, prov)


# ---------------------------------------------------------------------------
# Other source loaders.
# ---------------------------------------------------------------------------

def load_lifton_coeffs(path: str) -> np.ndarray:
    """7 CIRA-86 lapse-rate polynomial coefficients (K m^-1), whitespace/comma text."""
    coeffs = np.loadtxt(path).astype(np.float64).ravel()
    if coeffs.size != N_LIFTON:
        raise ValueError(f"lifton_lr: expected {N_LIFTON} coeffs, got {coeffs.size}")
    return coeffs


def load_vdm_record(path: str) -> tuple[np.ndarray, np.ndarray]:
    """VDM record: 2 columns (time_yrBP, vdm[1e22 A m^2]); coerced to ascending time."""
    data = np.loadtxt(path, delimiter=None).astype(np.float64)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError("vdm record: expected >=2 columns (time_yrBP, vdm)")
    order = np.argsort(data[:, 0])
    return data[order, 0].astype(np.float64), data[order, 1].astype(np.float32)


def load_era40_grids(path: str):
    """Pre-built ERA40-style NetCDF -> (lat, lon, meanP, meanT).

    Adapter for a file you already have. Edit the variable names to match yours;
    axis order is fixed automatically (lat descending, lon ascending).
    """
    ds = Dataset(path, "r")
    try:
        lat = np.asarray(ds.variables["lat"][:], dtype=np.float64)
        lon = np.asarray(ds.variables["lon"][:], dtype=np.float64)
        meanP = np.asarray(ds.variables["meanP"][:], dtype=np.float32)
        meanT = np.asarray(ds.variables["meanT"][:], dtype=np.float32)
    finally:
        ds.close()
    if lat[0] < lat[-1]:
        lat, meanP, meanT = lat[::-1], meanP[::-1, :], meanT[::-1, :]
    if lon[0] > lon[-1]:
        lon, meanP, meanT = lon[::-1], meanP[:, ::-1], meanT[:, ::-1]
    return lat, lon, meanP, meanT


def synthetic_era40(nlat=CANONICAL_NLAT, nlon=CANONICAL_NLON):
    """Uniform standard atmosphere (P=1013.25 hPa, T=288.15 K). Schema check only."""
    lat = np.linspace(90.0, -90.0, nlat).astype(np.float64)
    lon = np.linspace(0.0, 360.0, nlon).astype(np.float64)
    meanP = np.full((nlat, nlon), 1013.25, dtype=np.float32)
    meanT = np.full((nlat, nlon), 288.15, dtype=np.float32)
    return lat, lon, meanP, meanT


def synthetic_vdm(n=6000, span_yr=60000.0):
    """Constant modern VDM (8.0e22 A m^2), ascending 0..span_yr BP. Limit-check input."""
    t = np.linspace(0.0, span_yr, n).astype(np.float64)
    v = np.full(n, 8.0, dtype=np.float32)
    return t, v


# ---------------------------------------------------------------------------
# Validation + writer (components are optional; validate/write what is present).
# ---------------------------------------------------------------------------

def validate(era40, lifton, vdm) -> list[str]:
    errs: list[str] = []
    if era40 is not None:
        lat, lon, meanP, meanT = era40
        if lat.size < 2 or lon.size < 2:
            errs.append("era40 grid needs nlat>=2 and nlon>=2 (bilinear_era40).")
        if not np.all(np.diff(lat) < 0):
            errs.append("era40_lat must be strictly DESCENDING (+90 -> -90).")
        if not np.all(np.diff(lon) > 0):
            errs.append("era40_lon must be strictly ASCENDING (0 -> 360).")
        if meanP.shape != (lat.size, lon.size):
            errs.append(f"era40_meanP shape {meanP.shape} != ({lat.size},{lon.size}).")
        if meanT.shape != (lat.size, lon.size):
            errs.append(f"era40_meanT shape {meanT.shape} != ({lat.size},{lon.size}).")
        for name, arr in (("meanP", meanP), ("meanT", meanT)):
            if not np.all(np.isfinite(arr)):
                errs.append(f"era40_{name} has non-finite values.")
        if lat.size != CANONICAL_NLAT or lon.size != CANONICAL_NLON:
            print(f"  note: grid {lat.size}x{lon.size}, canonical is "
                  f"{CANONICAL_NLAT}x{CANONICAL_NLON} (allowed).")
    if lifton is not None and lifton.size != N_LIFTON:
        errs.append(f"lifton_lr must have {N_LIFTON} coeffs, got {lifton.size}.")
    if vdm is not None:
        mtime, mvdm = vdm
        if mtime.size != mvdm.size:
            errs.append(f"muscheler_time ({mtime.size}) != _vdm ({mvdm.size}).")
        if mtime.size < 2:
            errs.append("VDM record needs >=2 points.")
        elif not np.all(np.diff(mtime) > 0):
            errs.append("muscheler_time must be strictly ASCENDING (yr BP).")
        if not np.all(np.isfinite(mvdm)):
            errs.append("muscheler_vdm has non-finite values.")
    return errs


def write_nc(out_path, era40, lifton, vdm, provenance: str):
    ds = Dataset(out_path, "w", format="NETCDF4")
    try:
        if era40 is not None:
            lat, lon, meanP, meanT = era40
            ds.createDimension("nlat", lat.size)
            ds.createDimension("nlon", lon.size)
            ds.createVariable("era40_lat", "f8", ("nlat",))[:] = lat
            ds.createVariable("era40_lon", "f8", ("nlon",))[:] = lon
            ds.createVariable("era40_meanP", "f4", ("nlat", "nlon"))[:] = meanP
            ds.createVariable("era40_meanT", "f4", ("nlat", "nlon"))[:] = meanT
            ds.variables["era40_lat"].units = "degrees_north (descending)"
            ds.variables["era40_lon"].units = "degrees_east (ascending)"
            ds.variables["era40_meanP"].units = "hPa"
            ds.variables["era40_meanT"].units = "K"
        if lifton is not None:
            ds.createDimension("nlifton", lifton.size)
            ds.createVariable("lifton_lr", "f8", ("nlifton",))[:] = lifton
            ds.variables["lifton_lr"].units = "K m-1"
        if vdm is not None:
            mtime, mvdm = vdm
            ds.createDimension("ntime", mtime.size)
            ds.createVariable("muscheler_time", "f8", ("ntime",))[:] = mtime
            ds.createVariable("muscheler_vdm", "f4", ("ntime",))[:] = mvdm
            ds.variables["muscheler_time"].units = "yr BP (ascending)"
            ds.variables["muscheler_vdm"].units = "1e22 A m2"

        ds.title = "crn_parameters for LSDTopoTools3 (atmosphere + Lifton lapse rate + VDM)"
        ds.provenance = provenance
        ds.vdm_record = "Knudsen et al. (2008) EPSL 272:319, doi:10.1016/j.epsl.2008.04.048"
        ds.schema = "matches lsdtt_crn_parameters::load_nc_data (2026-07-23)"
    finally:
        ds.close()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, help="output crn_parameters.nc")
    src = ap.add_argument_group("atmosphere source (choose one)")
    src.add_argument("--fetch-era40", action="store_true",
                     help="download NCEP-DOE R2 and build the meanP/meanT grids")
    src.add_argument("--era40", help="pre-built ERA40-style NetCDF (load_era40_grids)")
    ap.add_argument("--era40-years", nargs=2, type=int, default=(1981, 2010),
                    metavar=("Y0", "Y1"), help="climatology window (default 1981 2010)")
    ap.add_argument("--era40-workdir", default=None,
                    help="where to stage downloads (default: a temp dir)")
    ap.add_argument("--keep-downloads", action="store_true",
                    help="do not delete the large NCEP files after building")
    ap.add_argument("--lifton", help="override the hardcoded CIRA-86 lapse-rate "
                    "coeffs with a text file of 7 numbers (default: hardcoded)")
    ap.add_argument("--vdm", help="override the hardcoded Knudsen 2008 VDM record with "
                    "a 2-col table: time_yrBP, vdm[1e22 A m2] (default: hardcoded)")
    ap.add_argument("--synthetic", action="store_true",
                    help="fill ALL components with schema-valid placeholders (smoke test)")
    args = ap.parse_args()

    era40 = lifton = vdm = None
    prov_bits = []

    if args.synthetic:
        era40 = synthetic_era40()
        lifton = np.zeros(N_LIFTON, dtype=np.float64)
        vdm = synthetic_vdm()
        prov_bits.append("SYNTHETIC placeholders (smoke test only)")
    else:
        if args.fetch_era40 and args.era40:
            ap.error("choose either --fetch-era40 or --era40, not both.")
        if args.fetch_era40:
            workdir = args.era40_workdir or tempfile.mkdtemp(prefix="ncep_r2_")
            lat, lon, meanP, meanT, prov = fetch_ncep_era40_grids(
                workdir, args.era40_years[0], args.era40_years[1], args.keep_downloads)
            era40 = (lat, lon, meanP, meanT)
            prov_bits.append(prov)
        elif args.era40:
            era40 = load_era40_grids(args.era40)
            prov_bits.append(f"era40={args.era40}")
        if args.lifton:
            lifton = load_lifton_coeffs(args.lifton)
            prov_bits.append(f"lifton={args.lifton}")
        else:
            lifton = LIFTON_LR_CIRA86.copy()
            prov_bits.append("lifton=CIRA-86 hardcoded (Balco/Lifton ERA40atm)")
        if args.vdm:
            vdm = load_vdm_record(args.vdm)
            prov_bits.append(f"vdm={args.vdm}")
        else:
            # Prepend the true present-day point so index 0 is the modern reference.
            mtime = np.concatenate([[PRESENT_DAY_VDM_TIME], KNUDSEN2008_VDM[:, 0]])
            mvdm = np.concatenate([[PRESENT_DAY_VDM],
                                   KNUDSEN2008_VDM[:, 1]]).astype(np.float32)
            vdm = (mtime, mvdm)
            prov_bits.append("vdm=Knudsen et al. 2008 (mmc1) + IGRF-2000 present-day point")

    if era40 is None and lifton is None and vdm is None:
        ap.error("nothing to write: give --fetch-era40/--era40, --lifton, --vdm, "
                 "or --synthetic")

    errs = validate(era40, lifton, vdm)
    if errs:
        print("VALIDATION FAILED:")
        for e in errs:
            print("  -", e)
        return 1

    write_nc(args.out, era40, lifton, vdm, "; ".join(prov_bits))

    # Report completeness against the loader's "loaded" requirement.
    have = {
        "era40_lat/meanP": era40 is not None,
        "lifton_lr": lifton is not None,
        "muscheler_vdm": vdm is not None,
    }
    complete = all(have.values())
    print(f"wrote {args.out}")
    if era40 is not None:
        lat, lon, meanP, meanT = era40
        print(f"  atmosphere {lat.size}x{lon.size}: "
              f"meanP {meanP.min():.0f}-{meanP.max():.0f} hPa, "
              f"meanT {meanT.min():.1f}-{meanT.max():.1f} K")
    if vdm is not None:
        print(f"  VDM {vdm[0].size} steps [{vdm[0][0]:.0f}..{vdm[0][-1]:.0f} yr BP]")
    print(f"  components present: " +
          ", ".join(k for k, v in have.items() if v) or "(none)")
    print(f"  file will load as {'COMPLETE (loaded=true)' if complete else 'PARTIAL'}"
          + ("" if complete else " -- add: "
             + ", ".join(k for k, v in have.items() if not v)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
