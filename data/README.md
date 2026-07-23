# lsdtt3 bundled data

Read-only data files shipped with lsdtt3 so the tools work out of the box without
any downloads. Committed to the repository and installed alongside the binaries.

## `crn_parameters.nc`

Parameter file for the cosmogenic-nuclide tools (`lsdtt-cosmo-tool`,
`lsdtt-cosmo-prediction`, `lsdtt-cosmo-column-tool`), read via the `crn_nc_file`
option. It supplies the ERA40-style atmosphere and the time-integrated
palaeomagnetic (time-varying production) correction (`sample_age_yr > 0`). Without
it the tools fall back to the standard (SLHL) atmosphere with no palaeomagnetic
correction.

Built by `python_helpers/build_crn_parameters_nc.py` (see that script and its
`_instructions.md` for the full schema and rebuild procedure). Regenerating it
requires a ~200 MB reanalysis download, so the ~100 kB product is committed here.

### Contents and provenance

| Variable | Meaning | Source |
|----------|---------|--------|
| `era40_lat`, `era40_lon` | 2.5° grid axes, 73×145 (lat 90→−90, lon 0→360) | — |
| `era40_meanP` | sea-level mean pressure (hPa) | NCEP-DOE Reanalysis 2, 1981–2010 annual climatology |
| `era40_meanT` | 1000 hPa mean temperature (K) — the sea-level reference in the barometric formula | NCEP-DOE Reanalysis 2, 1981–2010 annual climatology |
| `lifton_lr` | Lifton CIRA-86 lapse-rate polynomial (7 coeffs) | Balco/Lifton `ERA40atm` (Lifton, Sato & Dunai 2014) |
| `muscheler_time`, `muscheler_vdm` | virtual (axial) dipole moment vs age; absolute VDM in 10²² A m², 108 points, 0–49500 yr BP | Knudsen et al. (2008) supplement + a present-day (t=0) point (IGRF-2000 dipole moment, 7.79×10²² A m²) |

The atmosphere routine (`atmospheric_pressure_era40`) ports Balco/Lifton's
`ERA40atm`; `muscheler_vdm` is read as absolute VDM by `stone2000sp_paleomag`, which
takes index 0 as the modern reference (hence the prepended present-day point).

### Citations

- Kanamitsu, M., et al. (2002). NCEP-DOE AMIP-II Reanalysis (R-2). *Bull. Amer. Meteor. Soc.* 83, 1631–1643. doi:10.1175/BAMS-83-11-1631. (NCEP-DOE R2 is public domain; NOAA PSL.)
- Lifton, N., Sato, T., Dunai, T. (2014). Scaling in situ cosmogenic nuclide production rates using analytical approximations to atmospheric cosmic-ray fluxes. *Earth Planet. Sci. Lett.* 386, 149–160. doi:10.1016/j.epsl.2013.10.052.
- Knudsen, M.F., et al. (2008). Variations in the geomagnetic dipole moment during the Holocene and the past 50 kyr. *Earth Planet. Sci. Lett.* 272, 319–329. doi:10.1016/j.epsl.2008.04.048.
