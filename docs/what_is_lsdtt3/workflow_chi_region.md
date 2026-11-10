# Workflow: chi analysis of a region

The `chi_region` workflow goes from a box on the map to a chi analysis of the
basins in it:

1. [`lsdtt-fetch-raster`](../cli_reference/cli_reference_lsdtt_fetch_raster.md)
   fetches a DEM for the box (or the workflow reads your own DEM);
2. [`lsdtt-basin-extraction`](../cli_reference/cli_reference_lsdtt_basin_extraction.md)
   makes the basins: one above each of your own outlet points, or, if you give
   none, every complete basin (not cut by the DEM edge), not nested in another,
   and above a minimum area;
3. [`lsdtt-chi-analysis`](../cli_reference/cli_reference_lsdtt_chi_analysis.md)
   runs on exactly those basins: the chi coordinate, k_sn, the m/n disorder
   sweep, and piecewise fits of the chi-elevation profiles;
4. if [lsdviztools3](workflows.md#figures) is installed, draws a map of k_sn
   along the fitted channel segments, a stacked chi-elevation plot of all the
   basins, and the chi coordinate over the shaded DEM.

Picking the basins in a separate step matters on a fetched box: rivers run
out of the box, and a chi analysis is only meaningful on basins whose whole
drainage area is inside the DEM.

Set up the Snakemake environment first: see [Workflows](workflows.md#setting-up).

!!! warning "Figures need lsdviztools3"
    The figures listed on this page are drawn only if
    [lsdviztools3](https://lsdtopotools.github.io/lsdviztools3/getting_started/) is installed
    and the workflow can find it ([how to set it up](workflows.md#setting-up),
    step 4). Without it the run says the figures were skipped and writes
    everything else.

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

## Your own outlets

To analyse basins you have chosen, give their outlet points as
`basin_outlet_fname`: for every region, at the top of the `chi_region`
section, or under one region for that region only.

```yaml
chi_region:
  basin_outlet_fname: "/data/my_outlets.gpkg"   # every region
  regions:
    my_range:
      fetch: {west: 86.5, south: 27.8, east: 87.0, north: 28.1}
      basin_outlet_fname: "/data/my_range_outlets.csv"   # this region only
```

The file can be:

| Format | What it needs |
|---|---|
| CSV (`.csv`) | `latitude` and `longitude` columns in WGS84 decimal degrees, one row per outlet |
| FlatGeobuf (`.fgb`) | point features, in any coordinate system |
| GeoPackage (`.gpkg`) | point features, in any coordinate system |

Give an absolute path. One basin is made above each point. A point is rarely
exactly on the D8 channel, so it is first moved onto it; the
`basin_extraction_from_outlets` block sets how. The default
(`basin_outlet_snapping_method: downstream`) follows the flow downstream from
the point until it reaches a pixel with at least
`basin_outlet_contributing_pixels_threshold` (1000) contributing pixels; the
other methods are `expanding_square`, `largest_area` and `none` (see the
[lsdtt-basin-extraction reference](../cli_reference/cli_reference_lsdtt_basin_extraction.md)).
Points outside the DEM are skipped, with a warning in
`logs/basin_extraction.log`; if none fall on the DEM the run stops. The outlets actually used are written to
`basins/<region>_selected_basin_outlets.csv`. Basins are numbered
(`basin_index` in the chi outputs) in the order they appear in that file.

**Basins picked automatically.** Without outlets, the `basin_extraction` block
decides: raise `basins_min_area_pixels` for fewer, larger basins.

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
| `chi/<region>_chi_tributaries_mn_<m/n>.csv` | One row per tributary (a branch of the fit, from its head to where it joins its parent): its `source_key`, its drainage area just before the junction, that area as a fraction of the basin, its rank by area, and its k_sn range |
| `chi/<region>_chi_knickpoints_mn_<m/n>.csv` | Every node where two fitted segments meet: k_sn either side, `delta_ksn` (upstream minus downstream), and `at_junction` (1 at a tributary mouth) |
| `chi/<region>_chi_knickpoint_thresholds_mn_<m/n>.csv` | Suggested thresholds for a major knickpoint (90th and 95th percentiles, median + 3 MAD, Otsu) and an outlier k_sn for segments |
| `chi/<region>_chi_coordinate_mn_<m/n>.tif` | Chi for every pixel of the analysed basins |
| `chi/<region>_chi_disorder_sweep.csv` | Disorder against m/n |
| `figures/<region>_ksn_map.png` | k_sn (log scale) along the fitted channel segments, over a hillshade |
| `figures/<region>_chi_stacked.png` | Chi-elevation profiles of every basin side by side, coloured by k_sn |
| `figures/<region>_chi_profile_basin<N>.png` | One chi-elevation profile per basin, with `figures_per_basin_profiles: true` |
| `figures/<region>_chi_coordinate.png` | Chi over the shaded DEM |

The figures are drawn when lsdviztools3 is found, and need `m_over_n` set in
the `chi_analysis` block. They also switch on the chi outputs they draw
(`fit_chi_networks`, `write_fit_tables`, `write_chi_fit_segment_lines`,
`write_chi_coordinate_raster`).

`<m/n>` is written with `p` for the decimal point (`0p450`).

Every fit table carries `source_key`, the channel key of the chi data maps,
so a tributary can be found by the same key in every output. The tributary
and knickpoint tables are switched by `write_chi_tributary_table` and
`write_chi_knickpoint_table` (both on by default). lsdviztools3 draws them:
- `lsdviz3-chi` with `--remove-ksn-above`, `--min-tributary-area-m2`,
  `--min-tributary-area-fraction` or `--top-tributaries` leaves tributaries
  out of the plots (for example links routed over a dam onto a steep face),
  and maps what was removed;
- `lsdviz3-knickpoints` shows the distribution of k_sn changes with the
  suggested thresholds, and the major knickpoints on the profiles, along
  chi and on a map.
