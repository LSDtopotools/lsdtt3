# Workflow: valley-floor widths

The `valley_widths` workflow takes one DEM or a folder of them, and for each one:

1. finds the floodplain with
   [`lsdtt-fluvial-surfaces`](../cli_reference/cli_reference_lsdtt_fluvial_surfaces.md);
2. traces a valley centreline down every channel segment, and casts cross
   sections across it, with
   [`lsdtt-valley-centreline`](../cli_reference/cli_reference_lsdtt_valley_centreline.md)
   (`run_valley_width_network`). Each section's width is the span of the
   floodplain it crosses;
3. draws maps and profiles with lsdviztools3, if it is installed;
4. gathers every DEM's widths into one table and a per-DEM summary.

Set up the Snakemake environment first: see [Workflows](workflows.md#setting-up).
The methods themselves are described in
[Floodplains, terraces and valleys](../design_and_algorithms/fluvial_surfaces_and_valleys.md).

---

## Quick start

From the repository root, the shipped example fetches a 30 m Copernicus DEM
of a stretch of the River Spey in Scotland and runs the whole chain:

```bash
pixi run -e workflows snakemake --cores 1 --snakefile workflows/Snakefile \
    --directory workflows --configfile workflows/config/valley_widths.yaml
```

Results are in `workflows/results/valley-widths-example/`.

## Running it on your own DEMs

1. Put the DEMs in one folder. They must be in a projected coordinate
   system in metres (UTM, a national grid), with nodata set. The file stem
   names the DEM in every output (`glen_a.tif` -> `glen_a`), so use names
   without spaces.
2. Copy the config:

    ```bash
    cp workflows/config/valley_widths.yaml workflows/config/my_valleys.yaml
    ```

3. In the copy, set `run_id`, set `dem_dir` to your folder (an absolute path is
   safest), and **delete the `dems:` block**, which holds the example. With no
   `dems`, every file in `dem_dir` matching `dem_glob` is processed.
4. Set the channel threshold for your resolution (next section).
5. Dry run with `-n`, then run:

    ```bash
    pixi run -e workflows snakemake -n --cores 1 --snakefile workflows/Snakefile \
        --directory workflows --configfile workflows/config/my_valleys.yaml
    pixi run -e workflows snakemake --cores 1 --snakefile workflows/Snakefile \
        --directory workflows --configfile workflows/config/my_valleys.yaml
    ```

Add DEMs to the folder later and run again: only the new ones are processed,
and then the combined table is rebuilt.

---

## The parameters that matter most

Every block in the config is copied into the tool's parameter file as written;
[Changing the parameters](workflows.md#changing-the-parameters) explains how,
and how to give one DEM its own values.

**Channel threshold** (`shared_parameters.threshold_contributing_pixels`).
It sets where channels start, so it decides which valleys are measured. It is
given to both tools so they use the same network. It counts pixels, so the
same drainage area needs a different number at each resolution:

| Channel heads at | 30 m DEM | 10 m DEM | 5 m DEM | 1 m DEM |
|---|---|---|---|---|
| 0.1 km² | 111 | 1000 | 4000 | 100000 |
| 1 km² | 1111 | 10000 | 40000 | 1000000 |

If the DEMs have different resolutions, give each its own value under `dems:`.

**Floodplain method** (`fluvial_surfaces.fluvial_surface_algorithm_choice`).
The default config uses `clubb2017`: floodplain = pixels low above the nearest
channel (HAND) and gentle, with both bounds read off quantile-quantile plots
(Clubb et al., 2017), in patches that touch a channel.
`floodplain_second_pass_algorithm_choice: qq_upper_departure` then trims fans
and low terraces from that first pass. Set it to `none` to keep the first-pass
floodplain. The bounds used are written to
`<dem>_floodplain_thresholds.csv`. `clubb2017` defines channels by Strahler
order (`min_channel_stream_order`, 3), and stops with a message if the network
never reaches that order. Lower the order, or the channel threshold, if it does.

**Which mask the width is measured against** (`width_mask`). `floodplain`
(the default) measures the valley floor. `valley` measures the wider valley
mask, which takes in terraces.

**Cross sections** (`valley_centreline`):

| Parameter | Default here | Effect |
|---|---|---|
| `cross_section_spacing_metres` | 100 | Distance between sections along each centreline |
| `write_valley_cross_sections` | `true` | Section lines and the elevation profile along each one |
| `write_cross_section_features` | `false` | Valley walls, channel crossings and slope breaks along each section |
| `write_cross_section_hydraulics` | `false` | Stage-by-stage flow width, area and hydraulic depth of each section |
| `network_centreline_algorithm_choice` | `medial_path` | `medial_path` follows the middle of the valley floor; `smoothed_channel` smooths the channel itself, kept inside the mask |
| `write_relative_elevation_raster` | `false` | Elevation above the centreline's baseline, within the mask |

**Figures.** When [lsdviztools3](https://github.com/LSDtopotools/lsdviztools3)
is found, each DEM also gets a map of the sections coloured by
`excluded_reason`, a width profile down the largest trunk, a map of the
channels and centrelines, the floodplain over a hillshade, and the HAND QQ
plots its bounds were read from. Without lsdviztools3 the widths are still
computed, and the run says the figures were skipped. See
[Figures](workflows.md#figures) for how lsdviztools3 is found.

---

## Outputs

For each DEM, in `workflows/results/<run_id>/<dem>/`:

| File | Contents |
|---|---|
| `floodplain/<dem>_floodplain_mask.tif` | The floodplain the widths are measured against |
| `floodplain/<dem>_floodplain_polygons.fgb` | The same as polygons |
| `floodplain/<dem>_floodplain_thresholds.csv` | HAND and gradient bounds used |
| `floodplain/<dem>_reference_field.tif` | HAND (m) |
| `floodplain/<dem>_surface_class.tif` | 0 none, 1 floodplain, 2 terrace |
| `widths/<dem>_valley_widths.csv` | **One row per cross section**: segment, distance along the centreline, left and right half-widths, `valley_width_m`, channel elevation and drainage area at the station, and `excluded_reason` |
| `widths/<dem>_valley_segments.csv` | One row per channel segment: parent segment, nesting depth, head and end, a width estimate |
| `widths/<dem>_valley_centreline.fgb` | One centreline per segment |
| `widths/<dem>_valley_segment_channels.fgb` | Each segment's channel |
| `widths/<dem>_valley_cross_sections.fgb` | Section lines, with the same fields as the width table |
| `widths/<dem>_valley_cross_section_samples.csv` | Elevation samples along every section: offset from the centreline, elevation, whether in the valley. These are the cross-section profiles |
| `figures/<dem>_valley_widths_map.png`, `_profile.png`, `_channels.png`, `<dem>_floodplain_map.png`, `<dem>_floodplain_hand_qq.png` | When lsdviztools3 is found |
| `logs/` | Each tool's log |

And for the whole run, in `workflows/results/<run_id>/reports/`:

| File | Contents |
|---|---|
| `all_valley_widths.csv` | Every DEM's width table, with a `dem` column first |
| `valley_width_summary.csv` | Per DEM: number of sections, how many were kept and why the others were excluded, kept sections with zero width, and the quartiles of the kept widths |

### Which sections to use

Filter the width table on `excluded_reason == "none"`. The other values mark
sections that do not measure a valley floor:

| `excluded_reason` | Meaning |
|---|---|
| `not_in_valley` | The channel runs outside the floodplain mask at this station. Usually the most common value: steep reaches have no floodplain |
| `near_parent_channel` | A tributary's section inside its parent's valley; the parent measures that valley |
| `tributary_crossing` | A parent's section that runs up a tributary valley |
| `side_not_found` | No valley edge found within the section's length |

Some kept sections can have `valley_width_m` 0. The summary counts them
(`n_kept_zero_width`), so check that number before using the quartiles.

---

## Your own DEMs at different resolutions

```yaml
valley_widths:
  dem_dir: "/data/valleys"
  dems:
    upper_reach:
      file: "upper_reach_5m.tif"
      shared_parameters:
        threshold_contributing_pixels: 40000
    lower_reach:
      file: "lower_reach_10m.tif"
      shared_parameters:
        threshold_contributing_pixels: 10000
      valley_centreline:
        cross_section_spacing_metres: 200
```

Once `dems:` lists any DEM, only the listed DEMs are processed.
