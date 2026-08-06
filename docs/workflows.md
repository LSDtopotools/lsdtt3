# Workflows

`lsdtt3` includes a repository-level workflow layer under `workflows/` for
multi-stage, reproducible analysis chains.

The workflow runner is Snakemake, installed in a dedicated Pixi environment
(`workflows`) so workflow tooling is not mixed into the default Python
environment.

---

## Install and run

From the repository root:

```bash
pixi install -e workflows
pixi run -e workflows snakemake --help
pixi run -e workflows snakemake -n --snakefile workflows/Snakefile --directory workflows
```

Run a workflow:

```bash
pixi run -e workflows snakemake --cores 1 --snakefile workflows/Snakefile --directory workflows
```

---

## Current modules

Select a module with `analysis.module` in `workflows/config/config.yaml`.

### `check_basins`

Purpose: verify outlet points select the intended basins.

Chain:

1. fetch DEM from bounding box,
2. preprocess DEM,
3. extract basins from supplied outlets,
4. write basin and channel vectors as FlatGeobuf.

Inputs:

- Bounding box (`fetch.west/south/east/north`),
- outlet points file (`check_basins.outlets_file`) as:
  - CSV with `latitude` and `longitude` columns, or
  - FGB/GPKG point vector.

Outputs:

- `workflows/results/<run_id>/extract/check_junction_basin_outlines.fgb`
- `workflows/results/<run_id>/extract/check_channel_network_lines.fgb`

### `cosmo_landslides`

Purpose: run `lsdtt-cosmo-landslides` on a preprocessed DEM and write
landslide/CRN products.

### `cosmo_step_change_random_network`

Purpose: generate a seeded synthetic network DEM and run multi-isotope
step-change prediction.

### `cosmo_landslides_synthetic`

Purpose: generate a seeded synthetic network DEM, run a single-event
`lsdtt-cosmo-landslides` simulation, run a series+sample-ensemble diagnostic,
and produce landslide outlines and summary charts.

Key outputs:

- `workflows/results/<run_id>/preprocess/ls_synth_hillshade.tif`
- `workflows/results/<run_id>/analysis/cosmo_landslides_single/ls_single_LS_INVENTORY.csv`
- `workflows/results/<run_id>/analysis/cosmo_landslides_series/ls_series_LS_SAMPLE_EFFECT.csv`
- `workflows/results/<run_id>/reports/cosmo_landslides/ls_single_landslide_outlines.fgb`
- `workflows/results/<run_id>/reports/cosmo_landslides/ls_single_larsen_size_comparison.png`
- `workflows/results/<run_id>/reports/cosmo_landslides/ls_series_apparent_vs_true_erosion.png`

---

## Reproducibility

Each run writes a manifest JSON under `workflows/results/<run_id>/reports/`.
The manifest records key inputs and output hashes for replay/audit.

---

## Directory layout

- `workflows/Snakefile` - entry point
- `workflows/config/` - run configuration
- `workflows/rules/` - stage and module rules
- `workflows/scripts/` - helper scripts
- `workflows/results/` - generated outputs

For active details and module-specific parameters, see `workflows/README.md`.
