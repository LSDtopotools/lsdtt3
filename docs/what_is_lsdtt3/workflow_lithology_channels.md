# Workflow: lithology and channels

The `lithology_channels` workflow asks whether channels behave differently on
different rocks. It fetches a DEM and the [Macrostrat](https://macrostrat.org)
geologic map over the same area, then compares channel steepness and the
shape of the channel network between lithologic units:

1. [`lsdtt-fetch-raster`](../cli_reference/cli_reference_lsdtt_fetch_raster.md)
   fetches the DEM (or the workflow reads your own);
2. [`lsdtt-fetch-lithology`](../cli_reference/cli_reference_lsdtt_fetch_lithology.md)
   fetches the Macrostrat map units over the same box, in the DEM's
   coordinate system;
3. [`lsdtt-basin-extraction`](../cli_reference/cli_reference_lsdtt_basin_extraction.md)
   picks the complete basins;
4. [`lsdtt-chi-analysis`](../cli_reference/cli_reference_lsdtt_chi_analysis.md)
   computes k_sn in those basins and summarises it by lithologic unit: per-unit
   medians with bootstrap intervals, histograms, and unit-by-basin tables;
5. [`lsdtt-channel-network-analysis`](../cli_reference/cli_reference_lsdtt_channel_network_analysis.md)
   groups channel-link bearings, junction angles and Strahler statistics by
   lithologic unit, and compares each unit with the rest;
6. if [lsdviztools3](workflows.md#figures) is installed, draws the lithology
   map, k_sn by unit (box plots and histograms), and bearing roses and
   junction-angle summaries by unit.

How the statistics are built, and their caveats (a knickzone crossing a unit,
strike-parallel drainage, mixed reaches), are in
[Lithology and channel steepness statistics](../design_and_algorithms/lithology_channel_statistics.md)
and [Channel-link bearing statistics](../design_and_algorithms/channel_bearing_statistics.md).

Set up the Snakemake environment first: see [Workflows](workflows.md#setting-up).

!!! warning "Figures need lsdviztools3"
    The figures listed on this page are drawn only if
    [lsdviztools3](https://github.com/LSDtopotools/lsdviztools3) is installed
    and the workflow can find it ([how to set it up](workflows.md#setting-up),
    step 4). Without it the run says the figures were skipped and writes
    everything else.

## Quick start

```bash
pixi run -e workflows snakemake --cores 1 --snakefile workflows/Snakefile \
    --directory workflows --configfile workflows/config/lithology_channels.yaml
```

The example covers part of the Blue Ridge in North Carolina at 30 m. Results
are in `workflows/results/lithology-channels-example/blue_ridge/`.
Macrostrat coverage is most complete in North America; elsewhere check
`lithology/<region>_macrostrat_units.fgb` in QGIS before reading the statistics.

## Your own region

Regions, and your own outlet points (`basin_outlet_fname`: CSV, FlatGeobuf
or GeoPackage), are given as in the
[chi workflow](workflow_chi_region.md#your-own-region); see
[Your own outlets](workflow_chi_region.md#your-own-outlets).
lsdtt lays the lithology polygons on the DEM grid without reprojecting them,
so the map has to come back in the DEM's coordinate system:

- a fetched region passes its box and its `target_epsg` to the lithology
  fetch (when `target_epsg` is not given, both tools choose the UTM zone of the
  box centre, so they agree);
- a region read from a file must give its EPSG code:

    ```yaml
    lithology_channels:
      regions:
        my_area:
          file: "/data/my_dem.tif"
          epsg: 32617
    ```

**Which units.** `lithology_field_name` picks the map attribute that names a
unit. The default, `lith_primary`, is the first lithology Macrostrat lists for
the unit. `<region>_lithology_vocabulary.csv` in `lithology/` lists every value
found. To merge values into fewer units, write a crosswalk CSV with columns
`value,code` and optionally `name`, one row per map value, giving the values
of one unit the same code. `lithology/<region>_crosswalk_template.csv` lists
every value with its area. Give the crosswalk as `lithology_crosswalk_fname`
in both the `chi_analysis` and `channel_network_analysis` blocks, as an
absolute path: a bare file name is looked for in the DEM's folder.

## Outputs

In `workflows/results/<run_id>/<region>/`:

| File | Contents |
|---|---|
| `lithology/<region>_macrostrat_units.fgb` | The geologic map, with `lith_primary` and the other Macrostrat fields |
| `lithology/<region>_ATTRIBUTION.txt` | The map sources to cite |
| `chi_lithology/<region>_chi_ksn_by_lithology_mn_<m/n>.csv` | k_sn per lithologic unit: node and reach medians, quartiles, bootstrap interval, sample counts |
| `chi_lithology/<region>_chi_ksn_histogram_by_lithology_mn_<m/n>.csv` | k_sn histograms per unit, on shared bins |
| `chi_lithology/<region>_chi_ksn_by_lithology_and_basin_mn_<m/n>.csv` | The same per unit and basin |
| `chi_lithology/<region>_chi_ksn_lithology_nodes_mn_<m/n>.csv` | Every channel node with its unit and k_sn |
| `chi_lithology/<region>_lithology_units.csv` | Unit codes and names |
| `network/<region>_channel_bearing_group_summary.csv`, `_comparisons.csv` | Channel-link bearing statistics per unit and Strahler order, and each unit against the rest |
| `network/<region>_junction_angle_group_summary.csv`, `_comparisons.csv` | The same for junction angles |
| `network/<region>_strahler_channels.fgb`, `_channel_link_bearings.fgb` | The network, for maps |
| `*/<region>_lithology_code.tif` | The unit of every pixel, as each tool read it |
| `figures/<region>_lithology_map.png` | Units over the hillshade (when lsdviztools3 is found) |
| `figures/<region>_ksn_*.png` | k_sn by unit: box plots, histograms, unit by basin |
| `figures/<region>_lithology_rose_*.png`, `_junction_alpha_*.png`, `_bearing_comparisons_*.png` | Bearing roses per unit and Strahler order; junction angles; unit-against-rest comparisons |

Units with few channel nodes are flagged `low_sample` rather than dropped.
Map polygons with an empty lithology field form a unit with no name in the
tables. The lithology map labels it "no lithology in map".
