# HYPE oracle files

Reference output from the HYPE hydrological model, used by the lsdtt3 HYPE-column
tests. CTest reads these files only; it never downloads or builds HYPE.

## Provenance

| | |
|---|---|
| HYPE version | 5.36.0 |
| Source archive | `hype_5_36_0_src.tgz` |
| Download URL | https://sourceforge.net/projects/hype/files/release_hype_5_36_0/hype_5_36_0_src.tgz/download |
| Archive SHA-256 | `d90bc3670f341bb80105de94a23531775497ab472dbbc3560cd25b2afb03f43a` |
| SourceForge git master at generation | `fa9def0d359685c9cde0763f7040e253c345bb65` |
| Compiler | GNU Fortran (Ubuntu 11.4.0-1ubuntu1~22.04.3) 11.4.0 |
| Oracle build | shipped makefile, `make hype` (default REAL = single precision) |
| Precision-floor build | `make hype switch="-O2 -fbacktrace -ffree-line-length-0 -fno-automatic -std=gnu -cpp -fdefault-real-8 -fdefault-double-8"` |
| Output precision | 9 significant figures (`signfigures`); 9 is the most HYPE's output format accepts |

HYPE is developed by SMHI and released under the GNU Lesser General Public
Licence v3. These files are HYPE *outputs* for input sets written for lsdtt3;
no HYPE source is included.

## Regenerating

    python3 tests/scripts/generate_hype_oracle.py            # regenerate
    python3 tests/scripts/generate_hype_oracle.py --check    # verify committed files

The script uses the latest HYPE release that ships source. If that differs from
the version above, it stops until rerun with `--accept-new-version`.
Regenerating is a reviewable event: commit the changed files and state why.

## Layout

    cases/<case>.json             case definitions (the single source for both views)
    <case>/hype_input/            HYPE input files
    <case>/lsdtt_input/           the same set-up as lsdtt3 tables
    <case>/expected/class_<v>.csv per-class output; columns <subid>_<classid>
    <case>/expected/basin_<v>.csv per-sub-basin output; columns <subid>
    manifest.json                 provenance and precision floors

All runs use daily steps, one precipitation and one temperature series for the
whole domain, and HYPE's default model options unless a case lists others.
Initial state is HYPE's default: soil at field capacity, no snow, rivers at
their dead volume.

## Precision floor (tolerance basis)

The oracle is HYPE's released, single-precision build. The same cases were run
with a double-precision build; the largest difference per variable is the
floor below which a double-precision reimplementation cannot be judged against
the single-precision oracle.

| Case | Variable | Max abs diff single vs double | Max abs value |
|---|---|---|---|
| default_all | basin_cloc | 3.534e-06 | 9.02802 |
| default_all | basin_cout | 2.990e-06 | 13.752 |
| default_all | basin_crun | 1.210e-05 | 26.0972 |
| default_all | class_crun | 6.770e-05 | 31.5173 |
| default_all | class_evap | 2.310e-06 | 3.77008 |
| default_all | class_infi | 9.270e-05 | 43.1361 |
| default_all | class_melt | 1.025e-04 | 24.8043 |
| default_all | class_ros1 | 1.178e-05 | 4.15246 |
| default_all | class_ros2 | 1.023e-05 | 6.38587 |
| default_all | class_sml1 | 6.710e-05 | 144 |
| default_all | class_sml2 | 1.000e-04 | 198 |
| default_all | class_sml3 | 1.230e-04 | 440 |
| default_all | class_snow | 1.110e-04 | 506.702 |
| diversion | basin_prec | 0.000e+00 | 40 |
| diversion | basin_temp | 9.000e-07 | 20 |
| diversion | class_evap | 6.350e-06 | 4.257 |
| diversion | class_infi | 4.000e-07 | 30.4 |
| diversion | class_ros2 | 3.800e-07 | 9.6 |
| diversion | class_sml1 | 3.980e-05 | 50 |
| diversion | class_sml2 | 1.330e-05 | 102.556 |
| diversion | class_sml3 | 3.100e-05 | 280 |
| et | basin_prec | 0.000e+00 | 10 |
| et | basin_temp | 9.000e-07 | 17 |
| et | class_epot | 8.800e-07 | 8.31707 |
| et | class_evap | 8.590e-07 | 2.55174 |
| et | class_sml1 | 1.270e-05 | 49.8103 |
| et | class_sml2 | 2.350e-05 | 104.139 |
| et | class_sml3 | 1.500e-05 | 210 |
| percolation | basin_prec | 0.000e+00 | 30 |
| percolation | basin_temp | 5.000e-07 | 16 |
| percolation | class_evap | 4.000e-07 | 3.9 |
| percolation | class_infi | 0.000e+00 | 30 |
| percolation | class_sml1 | 1.480e-05 | 100 |
| percolation | class_sml2 | 1.700e-05 | 150 |
| percolation | class_sml3 | 3.100e-05 | 200 |
| routing | basin_cloc | 7.100e-07 | 3.29977 |
| routing | basin_cout | 8.600e-07 | 2.69779 |
| routing | basin_crun | 1.650e-06 | 10.209 |
| routing | class_crun | 1.650e-06 | 10.209 |
| snow | basin_prec | 0.000e+00 | 20 |
| snow | basin_temp | 4.600e-07 | 9 |
| snow | class_cprc | 8.000e-07 | 18 |
| snow | class_ctmp | 8.400e-07 | 9.9 |
| snow | class_infi | 7.770e-05 | 32.9217 |
| snow | class_melt | 7.770e-05 | 16.2065 |
| snow | class_snow | 9.400e-05 | 397.455 |
| snow | class_snwc | 1.070e-05 | 31.5201 |
| soil_runoff | basin_crun | 9.120e-06 | 14.6762 |
| soil_runoff | class_cro1 | 1.038e-05 | 12.0153 |
| soil_runoff | class_cro2 | 1.506e-05 | 7.19536 |
| soil_runoff | class_cro3 | 3.450e-06 | 3.23989 |
| soil_runoff | class_crun | 1.814e-05 | 18.2907 |
| soil_runoff | class_evap | 1.710e-06 | 3.9 |
| soil_runoff | class_infi | 1.100e-06 | 41.5 |
| soil_runoff | class_ros1 | 4.125e-06 | 3.7115 |
| soil_runoff | class_ros2 | 7.200e-08 | 3.5 |
| soil_runoff | class_sml1 | 2.370e-05 | 45 |
| soil_runoff | class_sml2 | 1.010e-04 | 135 |
| soil_runoff | class_sml3 | 1.260e-04 | 495 |

## Cases

- **default_all** (365 daily steps): All kept processes together: two sub-basins in series, three classes (3-, 2- and 1-layer soils, two land uses, two soils), lapse rate, interception, snow, PET/AET, diversion with macropores, percolation, soil runoff with slope term, saturation excess, local and main river routing. One hydrological year.
- **diversion** (365 daily steps): Threshold diversion of rain to infiltration-excess surface runoff (srrate, mactrinf, mactrsm), with ET on so layer 1 crosses the soil-moisture threshold. Macropore flow, percolation and soil runoff off; no snow.
- **et** (365 daily steps): Potential evaporation (option 0, seasonal coefficient) and actual evapotranspiration from the top two layers, with the PET layer split. No snow (temperature always above ttmp); diversion, percolation and soil runoff off.
- **percolation** (365 daily steps): Percolation between three layers with mperc1/mperc2 binding and non-binding, and the lowest layer filling so the receiving-layer limit binds. ET on; diversion, macropore flow and soil runoff off; no snow.
- **routing** (365 daily steps): Local and main river translation plus damping in two sub-basins in series (1 -> 2), with fractional translation steps (rivvel 0.5 m/s, damp 0.4). Column as in soil_runoff without snow.
- **snow** (181 daily steps): Class forcing (lapse-rate temperature correction, pcluse interception), rain/snow ramp and snowpack with liquid-water retention. ET, diversion, percolation and soil runoff are off.
- **soil_runoff** (365 daily steps): Soil-layer runoff with depth-decaying recession coefficients, stream depth in layer 1, 2 and 3 (one class each), two sub-basin slopes (2 % and 25 %) for the slope term; saturation-excess surface runoff; macropore placement and the second percolation call. ET on; no snow. Sub-basins are not connected.
