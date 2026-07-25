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
`_instructions.md` for the full schema and rebuild procedure), then extended in
place by `add_ncep_to_crn_nc.py` and `add_lsdn_to_crn_nc.py`. Regenerating it
requires a ~200 MB reanalysis download plus the Balco v3 reference `.mat` files, so
the ~1.25 MB product is committed here.

### Contents and provenance

| Variable | Meaning | Source |
|----------|---------|--------|
| `era40_lat`, `era40_lon` | 2.5° grid axes, 73×145 (lat 90→−90, lon 0→360) | — |
| `era40_meanP` | sea-level mean pressure (hPa) | NCEP-DOE Reanalysis 2, 1981–2010 annual climatology |
| `era40_meanT` | 1000 hPa mean temperature (K) — the sea-level reference in the barometric formula | NCEP-DOE Reanalysis 2, 1981–2010 annual climatology |
| `lifton_lr` | Lifton CIRA-86 lapse-rate polynomial (7 coeffs) | Balco/Lifton `ERA40atm` (Lifton, Sato & Dunai 2014) |
| `muscheler_time`, `muscheler_vdm` | virtual (axial) dipole moment vs age; absolute VDM in 10²² A m², 108 points, 0–49500 yr BP | Knudsen et al. (2008) supplement + a present-day (t=0) point (IGRF-2000 dipole moment, 7.79×10²² A m²) |
| `ncep_lat`, `ncep_lon` | 2.5° grid axes, 73×145 (lat 90→−90, lon 0→360) | — |
| `ncep_meanP` | NCEP2 mean sea-level pressure (hPa) | legacy `NCEP2.bin` (LSDTopoTools_CRNBasinwide); the atmosphere CAIRN/OCTOPUS use |
| `ncep_meanT` | NCEP2 mean 1000 hPa temperature (**°C**) | legacy `NCEP2.bin` |
| `lsdn_rc`, `lsdn_p` | LSDn grid axes: cutoff rigidity 0–20 GV step 0.1 (201), pressure 300–1058 hPa step 2 (380) | Balco v3 `LSDn2015_grid.mat` (`sfgrid.Rc`, `sfgrid.p`) |
| `lsdn_b_Be10`, `lsdn_m_Be10` | LSDn scaling grid for Be-10 in quartz, 380×201 as (p, Rc); raw scaling = `b + S·m` for solar modulation `S` (MV) | Balco v3 `LSDn2015_N10quartz.mat` (`sfN10quartz.b`, `.m`); Lifton et al. (2014), spectra from Sato et al. (2008) |
| `lsdn_b_Al26`, `lsdn_m_Al26` | same, Al-26 in quartz | `LSDn2015_N26quartz.mat` |
| `lsdn_b_C14`, `lsdn_m_C14` | same, C-14 in quartz | `LSDn2015_N14quartz.mat` |
| `lsdn_sphi`, `lsdn_sphi_inf` | solar modulation potential per Holocene step (91, MV) and its long-term mean (462.54 MV) | `LiftonS2015.mat` (`SPhi`, `SPhiInf`); Usoskin et al. (2011) as used by Lifton et al. (2014) |
| `pavon_mm01`, `pavon_t1min`, `pavon_t1max` | Holocene dipole-moment ratio `M_t/M_0` and its time-step bounds (91 steps, 0–13960 yr) | `PavonDip.mat` (`m.MM01`, `m.t1min`, `m.t1max`) |
| `pavon_pole_lat`, `pavon_pole_lon` | paleopole position per Holocene step (rad) | `PavonDip.mat` (`m.latr_pp`, `m.lonr_pp`) |
| `pavon_mm02`, `pavon_t2min`, `pavon_t2max` | long-record (GAD) dipole-moment ratio and time-step bounds (1710 steps, 13960 yr – 2.02044×10⁷ yr) | `PavonDip.mat` (`m.MM02`, `m.t2min`, `m.t2max`) |
| `pavon_t2` | long-record step **centres** — the `get_DipRc.m` branch threshold only, never an averaging weight (see below) | `PavonDip.mat` (`m.t2`) |
| `pavon_tzero` | reference year of the paleomagnetic records (2010 CE) | `PavonDip.mat` (`m.tzero`) |

The default atmosphere routine (`atmospheric_pressure_era40`) ports Balco/Lifton's
`ERA40atm`; `muscheler_vdm` is read as absolute VDM by `stone2000sp_paleomag`, which
takes index 0 as the modern reference (hence the prepended present-day point).

The `ncep_*` grids (added via `python_helpers/add_ncep_to_crn_nc.py`) supply the
alternative NCEP atmosphere (`atmospheric_pressure_ncep`, a port of Balco's
`NCEPatm_2`) selectable with `atmosphere: NCEP` in `lsdtt-cosmo-tool` — the
atmosphere CAIRN/OCTOPUS use, for like-for-like comparison. Only the SLP + 1000 hPa
temperature grids are stored (float32 + zlib): the legacy `NCEP_hgt.bin`
geopotential heights are never used by `NCEPatm_2`. **Note the units differ from
ERA40**: `ncep_meanT` is in °C, `era40_meanT` in K.

The `lsdn_*` and `pavon_*` variables (added via
`python_helpers/add_lsdn_to_crn_nc.py`) supply the **LSDn** (Lifton–Sato–Dunai 2014)
production scaling scheme. lsdtt3 evaluates LSDn by bilinearly interpolating the
same precomputed `b`/`m` grids that Balco's v3 online calculator interpolates, rather
than by porting the underlying Sato/PARMA spectral integration — so the grids, not a
quadrature, are the thing that has to be right. Storage order is `(lsdn_p, lsdn_rc)`
(pressure slow, rigidity fast), matching the source 380×201 shape with no transpose
at read time. The grids are float32 + zlib-9 + shuffle (round-trip ≈ 6×10⁻⁸
relative, far below the bilinear-interpolation error); the small records are float64
and uncompressed. The grids inherit the fixed ground water content used when they
were generated — lsdtt3 exposes no water-content knob.

`pavon_t2` deserves a specific warning. `get_DipRc.m` uses a step **centre** as the
branch threshold for "exposure longer than the long record" (`in.t <= m.t2(end)`,
= 2.00074×10⁷ yr), whereas `pavon_t2max[-1]` = 2.02044×10⁷ yr — about 1 % apart. Use
`pavon_t2` for **that branch condition only**; the duration-weighted time averaging
uses the bounds (`pavon_t1min/t1max`, `pavon_t2min/t2max`). Never use `pavon_t2` as
an averaging weight, and never use `pavon_t2max` as the branch threshold. The short
and long records abut exactly (`pavon_t1max[-1] == pavon_t2min[0]` = 13960 yr), which
is a property `get_DipRc.m` relies on.

The Balco v3 reference implementation these data come from is GPL-v2. What is
redistributed here is **data** — grid axes, coefficient tables, and published
paleomagnetic and solar records — cited as factual reference data on the same footing
as the Lifton lapse-rate polynomial and the Knudsen VDM record above. No MATLAB
source is copied into this repository. Each variable carries a `source` attribute
naming its `.mat` file and citation, and the file's global `lsdn_provenance`
attribute records the sha256 of every source `.mat`.

### Citations

- Kanamitsu, M., et al. (2002). NCEP-DOE AMIP-II Reanalysis (R-2). *Bull. Amer. Meteor. Soc.* 83, 1631–1643. doi:10.1175/BAMS-83-11-1631. (NCEP-DOE R2 is public domain; NOAA PSL.)
- Lifton, N., Sato, T., Dunai, T. (2014). Scaling in situ cosmogenic nuclide production rates using analytical approximations to atmospheric cosmic-ray fluxes. *Earth Planet. Sci. Lett.* 386, 149–160. doi:10.1016/j.epsl.2013.10.052.
- Knudsen, M.F., et al. (2008). Variations in the geomagnetic dipole moment during the Holocene and the past 50 kyr. *Earth Planet. Sci. Lett.* 272, 319–329. doi:10.1016/j.epsl.2008.04.048.
- Sato, T., Yasuda, H., Niita, K., Endo, A., Sihver, L. (2008). Development of PARMA: PHITS-based analytical radiation model in the atmosphere. *J. Nucl. Sci. Technol.* 45, 392–400. doi:10.1080/18811248.2008.9711443. (The neutron/proton spectra underlying the `lsdn_b_*` / `lsdn_m_*` grids.)
- Pavón-Carrasco, F.J., Osete, M.L., Torta, J.M., De Santis, A. (2014). A geomagnetic field model for the Holocene based on archaeomagnetic and lava flow data. *Earth Planet. Sci. Lett.* 388, 98–109. doi:10.1016/j.epsl.2013.11.046. (Source of the `pavon_*` Holocene record; compiled into the v3 form by N. Lifton.)
- Usoskin, I.G., Bazilevskaya, G.A., Kovaltsov, G.A. (2011). Solar modulation parameter for cosmic rays since 1936 reconstructed from ground-based neutron monitors and ionization chambers. *J. Geophys. Res.* 116, A02104. doi:10.1029/2010JA016105. (Basis of the `lsdn_sphi` solar modulation record as used by Lifton et al. 2014.)
