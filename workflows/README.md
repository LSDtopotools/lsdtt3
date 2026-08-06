# LSDTT workflows

This directory is the workflow layer for multi-stage LSDTT pipelines.

Planned stage model:

1. fetch,
2. preprocess,
3. extract,
4. analyse,
5. report.

Current implemented modules:

1. `check_basins` (bbox fetch -> preprocess -> outlet-driven basin extraction),
2. `cosmo_landslides`,
3. `cosmo_step_change_random_network`,
4. `cosmo_landslides_synthetic`.

Currently configured default module in `config/config.yaml`:

- `analysis.module: check_basins`
- chain: bounding-box DEM fetch -> preprocess -> extract basins from supplied
  outlets CSV/FGB -> write basin and channel FGB outputs.

## Quick start (check basins)

1. Set your run ID and outlet file in `workflows/config/config.yaml`.
2. Set your bounding box (`fetch.west/south/east/north`).
3. Run:

   `pixi run -e workflows snakemake --cores 1 --snakefile workflows/Snakefile --directory workflows`

4. Open these outputs in QGIS:

   - `workflows/results/<run_id>/extract/check_junction_basin_outlines.fgb`
   - `workflows/results/<run_id>/extract/check_channel_network_lines.fgb`

This module is designed to verify that supplied outlet points select the
expected basins.

## Synthetic landslide workflow quick start

Use the ready config at `workflows/config/cosmo_landslides_synthetic.yaml`:

`pixi run -e workflows snakemake --cores 1 --snakefile workflows/Snakefile --directory workflows --configfile workflows/config/cosmo_landslides_synthetic.yaml`

This run writes outputs under:

- `data/test_landslide_cosmo/workflow_results/test-landslide-cosmo/`

## Environment

Workflows use the dedicated Pixi environment defined in the repo root
`pixi.toml`.

- install: `pixi install -e workflows`
- run Snakemake: `pixi run -e workflows snakemake --help`
- dry-run check: `pixi run -e workflows snakemake -n --snakefile workflows/Snakefile --directory workflows`
- execute: `pixi run -e workflows snakemake --cores 1 --snakefile workflows/Snakefile --directory workflows`

## Layout

- `Snakefile` - top-level entry point.
- `config/` - workflow config and sweep tables.
- `rules/` - stage and domain rule modules.
- `scripts/` - manifest/replay/report helpers.
- `params/` - parameter templates by stage.
- `profiles/` - local/HPC execution profiles.
- `manifests/` - generated run manifests (not tracked by default).
- `results/` - generated workflow outputs (not tracked by default).

## Configure modules

Edit `workflows/config/config.yaml`:

- `analysis.module` selects the chain to run.
- `fetch.west/south/east/north` define the download area for fetch-based modules.
- `check_basins.outlets_file` accepts either CSV (`latitude`,`longitude`) or FGB/GPKG point vectors.
- `check_basins.*` controls outlet snapping and channel-source threshold.
- `cosmo_landslides.*` configures landslide + CRN analysis.
- `cosmo_landslides_synthetic.*`, `cosmo_landslides_single.*`,
  `cosmo_landslides_series.*`, and `cosmo_landslides_diagnostics.*` configure
  the synthetic landslide workflow (single event + series diagnostics).
- `synthetic_network.*` and `cosmo_step_change.*` configure random-network
  step-change analysis.

For `check_basins`, outputs include:

- `workflows/results/<run_id>/extract/check_junction_basin_outlines.fgb`
- `workflows/results/<run_id>/extract/check_channel_network_lines.fgb`

Optional debug output from basin-extraction can be enabled in
`rules/domain_check_basins.smk` by turning on `write_basin_outlet_points`.

For `cosmo_step_change_random_network`, outputs include:

- `workflows/results/<run_id>/generate/rn_ocn_dem.tif`
- `workflows/results/<run_id>/analysis/cosmo_step_change/sc_Be10_step_crn_conc.tif`
- `workflows/results/<run_id>/analysis/cosmo_step_change/sc_Al26_step_crn_conc.tif`
- `workflows/results/<run_id>/analysis/cosmo_step_change/sc_C14_step_crn_conc.tif`
- `workflows/results/<run_id>/analysis/cosmo_step_change/sc_R2610_step.tif`
- `workflows/results/<run_id>/analysis/cosmo_step_change/sc_R1410_step.tif`

For `cosmo_landslides`, outputs include:

- `workflows/results/<run_id>/analysis/cosmo_landslides/ls_LS_INVENTORY.csv`
- `workflows/results/<run_id>/analysis/cosmo_landslides/ls_LS_TFAIL.tif`
- `workflows/results/<run_id>/analysis/cosmo_landslides/ls_LS_EXCAVATION.tif`
- `workflows/results/<run_id>/analysis/cosmo_landslides/ls_LS_CONC.tif`

For `cosmo_landslides_synthetic`, outputs include:

- `workflows/results/<run_id>/preprocess/ls_synth_hillshade.tif`
- `workflows/results/<run_id>/analysis/cosmo_landslides_single/ls_single_LS_INVENTORY.csv`
- `workflows/results/<run_id>/analysis/cosmo_landslides_series/ls_series_LS_SAMPLE_EFFECT.csv`
- `workflows/results/<run_id>/reports/cosmo_landslides/ls_single_landslide_outlines.fgb`
- `workflows/results/<run_id>/reports/cosmo_landslides/ls_single_larsen_size_comparison.png`
- `workflows/results/<run_id>/reports/cosmo_landslides/ls_series_apparent_vs_true_erosion.png`
- `workflows/results/<run_id>/reports/cosmo_landslides/ls_series_sample_bias_ratio.png`
