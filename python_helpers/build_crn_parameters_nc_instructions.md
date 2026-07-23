# Building `crn_parameters.nc`

`build_crn_parameters_nc.py` assembles the `crn_parameters.nc` file that the
cosmogenic tools (`lsdtt-cosmo-tool`, `lsdtt-cosmo-prediction`) read via the
`crn_nc_file` option. It is required for the ERA40-style atmosphere and the
time-integrated palaeomagnetic (time-varying production) correction
(`sample_age_yr > 0`). Without it, the tools use the standard (SLHL) atmosphere
with no palaeomagnetic correction.

The script writes to the exact schema read by
`lsdtt::lsdtt_crn_parameters::load_nc_data` (variable names and dtypes are
load-bearing). **Two of the three inputs are built in** (Lifton coefficients and the
Knudsen VDM record are hardcoded), so a complete file needs only the atmosphere:

```bash
python build_crn_parameters_nc.py --fetch-era40 --out crn_parameters.nc
```

## Requirements

- Python with `netCDF4` and `numpy` (conda-forge).

## The three inputs

| Input | Flag | Default |
|-------|------|---------|
| Atmosphere (mean P + T grids) | `--fetch-era40` (download) or `--era40 <file>` | none — you must choose one |
| Lifton lapse-rate coeffs | `--lifton <file>` (override) | **hardcoded** CIRA-86 coeffs (Balco/Lifton ERA40atm) |
| VDM record | `--vdm <file>` (override) | **hardcoded** Knudsen et al. (2008) record |

The two hardcoded constants:
- **Lifton lapse-rate polynomial** — the 7 CIRA-86 coefficients verbatim from
  Balco/Lifton's `ERA40atm` (Lifton, Sato & Dunai 2014, doi:10.1016/j.epsl.2013.10.052).
- **Knudsen VDM record** — the 107-point mean VADM series (250–49500 yr BP, 10²² A m²)
  from the Knudsen et al. (2008) supplementary table (EPSL 272:319,
  doi:10.1016/j.epsl.2008.04.048; Muscheler is a coauthor, hence the `muscheler_*`
  variable names the loader reads), with a **true present-day point prepended at
  t = 0 yr BP** (IGRF-2000 dipole moment, 7.79×10²² A m²) → 108 points. `stone2000sp_paleomag`
  reads it as absolute VDM and divides by index 0 as the modern reference, so index 0
  is now the actual present-day field rather than the 250 yr BP sample.

Override either with a text file if you want a different record (e.g. Muscheler et al.
2005, doi:10.1016/j.quascirev.2005.01.012, via `--vdm`).

### Atmosphere: `--fetch-era40`

Downloads **NCEP-DOE Reanalysis 2** monthly means from NOAA PSL (public domain,
**no API key**), computes an annual climatological mean, and reduces them to the
canonical 2.5° **73×145** grid:

- `mslp.mon.mean.nc` (~9 MB) → `era40_meanP` (Pascals → hPa), **sea-level pressure**.
- `air.mon.mean.nc` (~190 MB) at the **1000 hPa** level → `era40_meanT` (K).

**Why 1000 hPa temperature, not 2 m temperature:** the C++ atmosphere routine
(`atmospheric_pressure_era40`) is a port of Balco/Lifton's `ERA40atm()`; `meanT` is
the near-sea-level (1000 hPa) reference temperature in the barometric formula, not a
surface field. NCEP-DOE R2 is the faithful lineage of the legacy LSDTopoTools NCEP2
atmosphere and is already on the 2.5° grid (the script only adds the 360° wrap
column). Citation: Kanamitsu et al. (2002), doi:10.1175/BAMS-83-11-1631.

The large `air` file is deleted after the grids are built (the "download big, save
small" pattern). Options:
- `--era40-years Y0 Y1` — climatology window (default `1981 2010`, a WMO normal).
- `--era40-workdir DIR` — where to stage downloads (default: a temp dir).
- `--keep-downloads` — keep the raw NCEP files (e.g. to build several climatologies).

If you already have a suitable NetCDF, use `--era40 <file>` instead and edit the
`load_era40_grids` adapter's variable names to match it.

### VDM record

The Knudsen (2008) record is hardcoded (from the paper's `mmc1` supplementary table),
so no download or manual table is needed. There is no clean DOI'd data file for it
(the raw paleointensities derive from GEOMAGIA50, http://geomagia.gfz-potsdam.de/), so
it is embedded rather than fetched. To use a different record, pass `--vdm <2-col
table>` (`time_yrBP  vdm`, VDM in 10²² A m²); the NetCDF variables keep the names
`muscheler_time` / `muscheler_vdm` because that is what the loader reads.

## Usage

Full canonical build (Lifton + Knudsen come from the built-in defaults):

```bash
python build_crn_parameters_nc.py --fetch-era40 --out crn_parameters.nc
```

Override a default and/or set the climatology window:

```bash
python build_crn_parameters_nc.py \
    --fetch-era40 --era40-years 1981 2010 \
    --vdm my_other_vdm.txt \
    --out crn_parameters.nc
```

Smoke test (schema-valid placeholders, no downloads):

```bash
python build_crn_parameters_nc.py --synthetic --out crn_synth.nc
```

The script validates inputs before writing (axis ordering, grid shape, 7 Lifton
coeffs, matched/ascending VDM record, finite values), writes only the components
supplied, and prints whether the file will load as COMPLETE (`loaded=true`) or
PARTIAL.

## Schema written

| Variable | dtype | dims | units |
|----------|-------|------|-------|
| `era40_lat` | double | (nlat,) | deg, descending (+90→−90) |
| `era40_lon` | double | (nlon,) | deg, ascending (0→360) |
| `era40_meanP` | float | (nlat, nlon) | hPa (sea-level pressure) |
| `era40_meanT` | float | (nlat, nlon) | K (1000 hPa temperature) |
| `lifton_lr` | double | (7,) | K m⁻¹ |
| `muscheler_time` | double | (ntime,) | yr BP, ascending |
| `muscheler_vdm` | float | (ntime,) | 10²² A m² |

Canonical grid is 73×145; other grid sizes are accepted (a note is printed).
`loaded` requires `era40_lat`, `era40_meanP`, `lifton_lr`, `muscheler_vdm`.
