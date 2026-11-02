# LSDTT workflows

Snakemake pipelines that chain several lsdtt3 command-line tools. The user
guide is `docs/what_is_lsdtt3/workflows.md` (setting up the `workflows` pixi
environment, running, and changing parameters), with one page per module:

| Module | Config | Guide |
|---|---|---|
| `valley_widths` | `config/valley_widths.yaml` | `docs/what_is_lsdtt3/workflow_valley_widths.md` |
| `chi_region` | `config/chi_region.yaml` | `docs/what_is_lsdtt3/workflow_chi_region.md` |
| `lithology_channels` | `config/lithology_channels.yaml` | `docs/what_is_lsdtt3/workflow_lithology_channels.md` |
| `check_basins` | `config/config.yaml` (the default module) | `docs/what_is_lsdtt3/workflows.md` |
| `cosmo_landslides`, `cosmo_step_change_random_network`, `cosmo_landslides_synthetic`, `cosmo_transient_landslides` | `config/config.yaml`, `config/cosmo_landslides_synthetic.yaml` | below |

Quick start, from the repository root:

```bash
pixi install -e workflows
pixi run -e workflows snakemake -n --cores 1 --snakefile workflows/Snakefile \
    --directory workflows --configfile workflows/config/valley_widths.yaml
pixi run -e workflows snakemake --cores 1 --snakefile workflows/Snakefile \
    --directory workflows --configfile workflows/config/valley_widths.yaml
```

`config/config.yaml` is always loaded first; a `--configfile` is merged over
it. Relative paths in a config are relative to `workflows/`.

## Layout of the rules

- `rules/common_params.smk`: helpers for the per-item modules. A config block
  of CLI parameters is copied verbatim into the CLI's parameter file; the
  workflow sets only file I/O and the keys that chain one CLI to the next
  (`render_param_file` refuses a config that sets one of those).
- `rules/common_regions.smk`: DEM fetch, basin selection and Macrostrat fetch
  per region, shared by `chi_region` and `lithology_channels`.
- `rules/domain_*.smk`: one file per module.
- `rules/common_fetch.smk`, `common_preprocess.smk`, `common_extract.smk`:
  the single-DEM stages used by `check_basins` and `cosmo_landslides`.
- The Snakefile's `shell.prefix` unsets the pixi environment's GDAL/PROJ
  variables before every command (they make the CLIs load the wrong
  plugins), and sources `lsdtt_setenv_script` when a config names one.

- Figures: `figures_enabled()` in `common_params.smk` probes for lsdviztools3
  (`lsdviz3-drape --help` behind `figures.lsdviz3_command_prefix`) once per run
  and caches the answer in the environment for child processes; figure rules
  use `run_lsdviz3`, which runs without the GDAL/PROJ variables and without the
  setenv script. `figures.enabled`: auto (default) | true | false.

Every CLI's log is kept under `results/<run_id>/<item>/logs/`; on failure
its last lines are printed.

## Synthetic landslide workflow quick start

Use the ready config at `workflows/config/cosmo_landslides_synthetic.yaml`:

`pixi run -e workflows snakemake --cores 1 --snakefile workflows/Snakefile --directory workflows --configfile workflows/config/cosmo_landslides_synthetic.yaml`

This run writes outputs under:

- `data/test_landslide_cosmo/workflow_results/test-landslide-cosmo/`

## Cosmogenic transient + landslide ensemble

`analysis.module: cosmo_transient_landslides` is the Snakemake translation of
the hand-run cosmogenic transient-landslide experiment. It runs

1. a step-change landscape per magnitude (`lsdtt-cosmo-prediction`), giving
   `SC_DEM`, `SC_ERATE`, `SC_STEPTIME`, `SC_TAG`, `SC_CHI`;
2. the channel-network sample points on that surface
   (`lsdtt-channel-extraction`);
3. a landslide field per (magnitude, seed) confined to `SC_TAG == 4` -- the
   post-knickpoint hillslopes -- with a Larsen area-volume scaling and a
   Malamud inverse-gamma size distribution (`lsdtt-cosmo-landslides`);
4. a forward multi-isotope prediction per (magnitude, seed, coupling) and three
   inversions per run, one each for `Be10`, `Al26` and `C14`
   (`lsdtt-cosmo-prediction`, `lsdtt-cosmo-tool`);
5. the paired join and the seed-ensemble statistics.

The ensemble is declarative: seed, step magnitude and coupling are wildcards,
so `cosmo_transient_landslides.seeds`, `.magnitudes` and `.scenarios` in
`config/config.yaml` are the only place the matrix size is stated. The shipped
configuration is 10 seeds x 2 magnitudes x 4 couplings plus one control per
magnitude = 82 forward runs and 246 inversions.

Run it with:

`pixi run -e workflows snakemake --cores 1 --snakefile workflows/Snakefile --directory workflows --config analysis='{"module":"cosmo_transient_landslides"}'`

These runs take minutes each; check the DAG with `-n` first.

Points that are easy to get wrong, and are pinned in the config:

- `step_change.route_on_generated_dem` is **false**. The generated step-change
  DEM has a ~56 m elevation cliff along drainage divides, so rebuilding D8 on
  it sends 5.4 % of tag-4 pixels draining into the unadjusted catchment.
- `crn_nc_fname` resolves against `read_path`, not the working directory, so the
  CRN netCDF is staged next to every DEM a run reads. The same is true of the
  landslide `*_raster_prefix` parameters.
- `crn.rock_density_kg_m3` (2650) and `landslides.scar_density_kg_m3` (1500)
  are different quantities and are not interchangeable.
- the comparison is paired per sample point against the matching control at the
  same magnitude; the per-seed statistic is the median of those paired deltas
  and the ensemble statistic is the median over seeds with an IQR.

Outputs:

- `workflows/results/<run_id>/reports/cosmo_transient_landslides/ensemble_per_seed.csv`
- `workflows/results/<run_id>/reports/cosmo_transient_landslides/ensemble_summary.csv`
- `workflows/results/<run_id>/reports/cosmo_transient_landslides/ensemble_anchor_check.csv`
- `workflows/results/<run_id>/reports/cosmo_transient_landslides/ensemble_report.txt`

`ensemble_anchor_check.csv` is the validation output: the 32 samples that lie
entirely upstream of the knickpoint cannot have a landslide above them, so
their paired delta must be zero to floating point. The rule reports
max|delta| per combination and fails if any exceeds
`analysis.anchor_alarm`. The same rule asserts that S1 (`both`) and S1f
(`flush`) are identical row for row, and that every scenario sample pairs 1:1
with the control.

## Environment

Workflows use the dedicated Pixi environment defined in the repo root
`pixi.toml`.

- install: `pixi install -e workflows`
- run Snakemake: `pixi run -e workflows snakemake --help`
- dry-run check: `pixi run -e workflows snakemake -n --snakefile workflows/Snakefile --directory workflows`
- execute: `pixi run -e workflows snakemake --cores 1 --snakefile workflows/Snakefile --directory workflows`

CLIs whose `n_threads` is `0` use all available hardware threads. On shared
hosts, set an explicit `n_threads` no greater than the physical cores allocated
to each job; Snakemake's `--cores` limit does not by itself rewrite CLI thread
parameters.

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
- `cosmo_transient_landslides.*` configures the transient step-change +
  landslide seed ensemble.

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
