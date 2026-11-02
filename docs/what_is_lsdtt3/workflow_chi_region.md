# Workflow: chi analysis of a region

The `chi_region` workflow goes from a box on the map to a chi analysis of the
basins in it:

1. [`lsdtt-fetch-raster`](../cli_reference/cli_reference_lsdtt_fetch_raster.md)
   fetches a DEM for the box (or the workflow reads your own DEM);
2. [`lsdtt-basin-extraction`](../cli_reference/cli_reference_lsdtt_basin_extraction.md)
   picks the basins and writes their outlets. By default these are every
   complete basin (not cut by the DEM edge), not nested in another, and above
   a minimum area;
3. [`lsdtt-chi-analysis`](../cli_reference/cli_reference_lsdtt_chi_analysis.md)
   runs on exactly those basins: the chi coordinate, k_sn, the m/n disorder
   sweep, and piecewise fits of the chi-elevation profiles;
4. if [lsdviztools3](workflows.md#figures) is installed, draws the chi
   coordinate over the shaded DEM. Figures need `m_over_n` set in the
   `chi_analysis` block, because the chi raster is named after it.

Picking the basins in a separate step matters on a fetched box: rivers run
out of the box, and a chi analysis is only meaningful on basins whose whole
drainage area is inside the DEM.

Set up the Snakemake environment first: see [Workflows](workflows.md#setting-up).

## Quick start

```bash
pixi run -e workflows snakemake --cores 1 --snakefile workflows/Snakefile \
    --directory workflows --configfile workflows/config/chi_region.yaml
```

The example fetches a 30 m Copernicus DEM of part of the Eastern Carpathians.
Results are in `workflows/results/chi-region-example/carpathians/`.

## Your own region

Copy `workflows/config/chi_region.yaml` and, under `regions:`, give each
region a name and either a `fetch:` block or a `file:`:

```yaml
chi_region:
  regions:
    my_range:
      fetch:                      # any lsdtt-fetch-raster parameters
        west: 86.5
        south: 27.8
        east: 87.0
        north: 28.1
        dem_source: "cop30_aws"   # or a national source, e.g. us_3dep
        grid_spacing: 30
    my_lidar:
      file: "/data/my_lidar_dem.tif"
      shared_parameters:
        threshold_contributing_pixels: 20000
```

Several regions run in one go; each gets its own folder of results. The
workflow sets `basin_outlet_fname` and `basin_outlet_snapping_method` in the
chi analysis, so it analyses the basins from step 2 exactly.

**Basins.** Change the `basin_extraction` block to choose others: raise
`basins_min_area_pixels` for fewer, larger basins, or replace the selection
with `select_basins_from_outlets: true` and `basin_outlet_fname: <your
outlets file>` to analyse basins you have chosen.

**Chi.** `m_over_n` fixes the concavity used for k_sn; the disorder sweep
(`fit_chi_networks: true`) tests the range `mn_min` to `mn_max` and writes
`<region>_chi_disorder_sweep.csv`. Every other parameter of
`lsdtt-chi-analysis` can be added to the `chi_analysis` block. See
[Changing the parameters](workflows.md#changing-the-parameters).

## Outputs

In `workflows/results/<run_id>/<region>/`:

| File | Contents |
|---|---|
| `dem/<region>_DEM.tif` | The fetched DEM (with its provenance in `<region>_README.html`) |
| `basins/<region>_selected_basin_outlets.csv` | The outlets analysed, with basin areas |
| `basins/<region>_junction_basin_polygons.fgb` | Their outlines |
| `chi/<region>_chi_ksn_basin_summary_mn_<m/n>.csv` | k_sn per basin (mean, median, range, main stem) |
| `chi/<region>_chi_fit_segments_mn_<m/n>.csv`, `_fit_segment_lines_mn_<m/n>.fgb` | Piecewise chi fits, as a table and as lines |
| `chi/<region>_chi_fit_nodes_mn_<m/n>.csv` | Every channel node with chi, elevation and fitted k_sn |
| `chi/<region>_chi_coordinate_mn_<m/n>.tif` | Chi for every pixel of the analysed basins |
| `chi/<region>_chi_disorder_sweep.csv` | Disorder against m/n |
| `figures/<region>_chi_coordinate.png` | Chi over the shaded DEM (when lsdviztools3 is found) |

`<m/n>` is written with `p` for the decimal point (`0p450`).
