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
    <case>/lsdtt_input/           the set-up as lsdtt3 tables (the HYPE input files are
                                  regenerated from cases/<case>.json; --keep-work keeps them)
    <case>/expected/class_<v>.csv per-class output; columns <subid>_<classid>
    <case>/expected/basin_<v>.csv per-sub-basin output; columns <subid>
    manifest.json                 provenance and precision floors

Values are HYPE's printed 9 significant figures, written in shortest form
('0', '12.3456789'). Forcing uses a compressed seasonal cycle
(temp_period_steps) so each case is short.

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
| default_all | basin_cout | 2.831e-06 | 5.69663 |
| default_all | basin_crun | 1.088e-05 | 11.156 |
| default_all | class_crun | 5.434e-05 | 14.3524 |
| default_all | class_evap | 7.800e-07 | 2.88719 |
| default_all | class_infi | 4.900e-06 | 35.9755 |
| default_all | class_melt | 4.930e-06 | 15.75 |
| default_all | class_ros1 | 1.760e-06 | 2.09236 |
| default_all | class_ros2 | 9.130e-07 | 5.56703 |
| default_all | class_sml1 | 2.900e-05 | 144 |
| default_all | class_sml2 | 7.000e-05 | 198 |
| default_all | class_sml3 | 1.080e-04 | 440 |
| default_all | class_snow | 1.010e-05 | 99.3692 |
| diversion | class_evap | 9.500e-07 | 4.23 |
| diversion | class_infi | 4.000e-07 | 40 |
| diversion | class_ros2 | 3.800e-07 | 9.6 |
| diversion | class_sml1 | 5.600e-06 | 50 |
| diversion | class_sml2 | 1.390e-05 | 102.589 |
| diversion | class_sml3 | 3.100e-05 | 280 |
| et | class_epot | 8.500e-07 | 8.22266 |
| et | class_evap | 1.420e-06 | 4.00584 |
| et | class_sml1 | 9.400e-06 | 49.6348 |
| et | class_sml2 | 1.200e-05 | 104.16 |
| et | class_sml3 | 1.500e-05 | 210 |
| percolation | class_evap | 3.300e-07 | 3.9 |
| percolation | class_infi | 0.000e+00 | 30 |
| percolation | class_sml1 | 4.310e-05 | 100 |
| percolation | class_sml2 | 2.300e-05 | 150 |
| percolation | class_sml3 | 3.100e-05 | 200 |
| routing | basin_cloc | 6.700e-07 | 3.27818 |
| routing | basin_cout | 5.500e-07 | 2.67926 |
| routing | basin_crun | 1.400e-06 | 10.1827 |
| snow | class_cprc | 8.000e-07 | 18 |
| snow | class_ctmp | 8.400e-07 | 9.9 |
| snow | class_infi | 2.200e-06 | 24.887 |
| snow | class_melt | 1.810e-06 | 14.9022 |
| snow | class_snow | 1.200e-05 | 114.057 |
| snow | class_snwc | 1.740e-06 | 9.11974 |
| soil_runoff | class_cro1 | 9.190e-06 | 11.942 |
| soil_runoff | class_cro2 | 6.715e-06 | 7.05802 |
| soil_runoff | class_cro3 | 2.320e-06 | 2.56553 |
| soil_runoff | class_crun | 1.190e-05 | 18.1885 |
| soil_runoff | class_ros1 | 3.120e-06 | 3.66238 |
| soil_runoff | class_ros2 | 7.200e-08 | 3.5 |
| soil_runoff | class_sml1 | 2.350e-05 | 45 |
| soil_runoff | class_sml2 | 4.900e-05 | 135 |
| soil_runoff | class_sml3 | 7.300e-05 | 495 |

## Cases

- **default_all** (150 daily steps): All kept processes together: two sub-basins in series, three classes (3-, 2- and 1-layer soils, two land uses, two soils), lapse rate, interception, snow, PET/AET, diversion with macropores, percolation, soil runoff with slope term, saturation excess, local and main river routing. 150 daily steps spanning one compressed seasonal cycle (snow accumulation and melt included).
- **diversion** (120 daily steps): Threshold diversion of rain to infiltration-excess surface runoff (srrate, mactrinf, mactrsm), with ET on so layer 1 crosses the soil-moisture threshold. Macropore flow, percolation and soil runoff off; no snow.
- **et** (120 daily steps): Potential evaporation (option 0, seasonal coefficient) and actual evapotranspiration from the top two layers, with the PET layer split. No snow (temperature always above ttmp); diversion, percolation and soil runoff off.
- **percolation** (120 daily steps): Percolation between three layers with mperc1/mperc2 binding and non-binding, and the lowest layer filling so the receiving-layer limit binds. ET on; diversion, macropore flow and soil runoff off; no snow.
- **routing** (120 daily steps): Local and main river translation plus damping in two sub-basins in series (1 -> 2), with fractional translation steps (rivvel 0.5 m/s, damp 0.4). Column as in soil_runoff without snow.
- **snow** (120 daily steps): Class forcing (lapse-rate temperature correction, pcluse interception), rain/snow ramp and snowpack with liquid-water retention. ET, diversion, percolation and soil runoff are off.
- **soil_runoff** (120 daily steps): Soil-layer runoff with depth-decaying recession coefficients, stream depth in layer 1, 2 and 3 (one class each), two sub-basin slopes (2 % and 25 %) for the slope term; saturation-excess surface runoff; macropore placement and the second percolation call. ET on; no snow. Sub-basins are not connected.
