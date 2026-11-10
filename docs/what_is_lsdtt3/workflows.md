# Workflows

Some analyses take more than one command-line tool: one tool's output is the
next tool's input. The `workflows/` directory of the repository chains these
tools with [Snakemake](https://snakemake.readthedocs.io). You describe the run
in one configuration file, and Snakemake writes each tool's parameter file,
runs the tools in order, and reruns only what a change affects.

| Module | Chain | Guide |
|---|---|---|
| `valley_widths` | DEM(s) -> floodplain -> valley centrelines -> widths and cross sections | [Valley-floor widths](workflow_valley_widths.md) |
| `chi_region` | fetch DEM -> basins -> chi analysis | [Chi analysis of a region](workflow_chi_region.md) |
| `lithology_channels` | fetch DEM + Macrostrat geology -> basins -> k_sn and network statistics by rock unit | [Lithology and channels](workflow_lithology_channels.md) |
| `check_basins` | fetch DEM -> preprocess -> basins from outlet points | [below](#check_basins) |
| `cosmo_*` | cosmogenic-nuclide scenario ensembles | `workflows/README.md` |

---

## Setting up

You need three things to run the analyses: the lsdtt3 repository, the
command-line tools, and the Snakemake environment. **For figures you also need
lsdviztools3** (step 4), a separate install. Without it the workflows still run
and write every table, raster and vector file, but draw no figures.

**1. The repository.** The workflows are files in the repository, not part of
the installed package:

```bash
git clone https://github.com/LSDtopotools/lsdtt3.git
cd lsdtt3
```

**2. The command-line tools.** Either build them from this checkout
(see [Building from source](../developer_notes/building_from_source.md)); the workflows look
for them in `build/` by default. Or install a
[package](install_lsdtt3.md) and tell the workflow where its `setenv.sh` is
(see [Using an installed package](#using-an-installed-package)).

**3. The Snakemake environment.** The workflows run in a
[pixi](https://pixi.sh) environment called `workflows`, kept apart from the
default Python environment. Install pixi once:

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

Open a new terminal so `pixi` is on your `PATH`, then, from the repository root:

```bash
pixi install -e workflows
pixi run -e workflows snakemake --version
```

The second command prints the Snakemake version, which means the environment is ready.
`pixi install` is needed only once; `pixi run -e workflows` uses the
environment without activating anything.

**4. lsdviztools3, for figures (optional).** The Snakemake environment holds no
plotting software. Every figure is drawn by a command from
[lsdviztools3](https://lsdtopotools.github.io/lsdviztools3/cli_reference/)
(`lsdviz3-valley-widths`, `lsdviz3-drape`, `lsdviz3-lithology-map`, ...), which
lives in its own pixi environment with GMT, PyGMT and geopandas, roughly 2 GB.
Install it next to lsdtt3:

```bash
cd ..        # the folder that holds lsdtt3
git clone https://github.com/LSDtopotools/lsdviztools3.git
cd lsdviztools3
pixi install
pixi run install
pixi run lsdviz3-drape --help      # prints the help: lsdviztools3 is ready
```

Then tell the workflow how to run it, in the `figures` block of your config:

```yaml
figures:
  lsdviz3_command_prefix: "pixi run --manifest-path /path/to/lsdviztools3/pixi.toml"
```

If lsdviztools3's commands are already on your `PATH`, leave the prefix empty.
See [Figures](#figures) for what happens when lsdviztools3 is not found.

---

## Running a workflow

Every command below is run from the repository root. Each module has a
configuration file in `workflows/config/`. Copy it, edit the copy, and give
the copy to Snakemake with `--configfile`:

```bash
cp workflows/config/valley_widths.yaml workflows/config/my_valleys.yaml
# edit workflows/config/my_valleys.yaml

# Dry run: list the jobs without running them
pixi run -e workflows snakemake -n --cores 1 --snakefile workflows/Snakefile \
    --directory workflows --configfile workflows/config/my_valleys.yaml

# Run
pixi run -e workflows snakemake --cores 4 --snakefile workflows/Snakefile \
    --directory workflows --configfile workflows/config/my_valleys.yaml
```

- `--cores` is how many jobs run at once. Each tool also uses several threads
  of its own, so with several large DEMs start with `--cores 1` or `2`.
- Outputs go to `workflows/results/<run_id>/`. Change `run_id` in the config
  to keep runs apart.
- Run the same command again after editing the config, and only the steps
  whose parameters changed run again, with everything downstream of them.
  `--forceall` reruns everything.
- **Relative paths in a config are relative to `workflows/`**, not to where
  you typed the command. Absolute paths are safest for your own data.

When a tool fails, Snakemake prints the end of that tool's log. The full logs
are in `workflows/results/<run_id>/<item>/logs/`.

---

## Changing the parameters

The configuration files of `valley_widths`, `chi_region` and
`lithology_channels` hold one block per tool. **Each block is copied into
that tool's parameter file as written**, so any parameter from the tool's
reference page (or from `lsdtt-<tool> -h`) can be added, changed, or removed
there. Nothing has to change in the workflow itself.

```yaml
valley_widths:
  fluvial_surfaces:                     # -> lsdtt-fluvial-surfaces
    fluvial_surface_algorithm_choice: "clubb2017"
    min_channel_stream_order: 3
  valley_centreline:                    # -> lsdtt-valley-centreline
    cross_section_spacing_metres: 100
```

Values are written as `key: value`. A file parameter you add yourself (an
outlets file, a crosswalk) is read relative to the DEM's folder, as the tool
does, so give it as an absolute path. `true`/`false` become `true`/`false`, and
a YAML list becomes a comma-separated list (`[11, 31]` -> `11,31`). A
misspelt parameter stops the tool with a message naming it, because every
tool rejects parameters it does not know.

The parameters are merged in this order, later ones winning:

1. `shared_parameters`, given to every tool of the module (for example the
   channel threshold, which must be the same in each tool);
2. the tool's own block (`fluvial_surfaces`, `chi_analysis`, ...);
3. the item's own `shared_parameters` and tool blocks, under that DEM or
   region. Use these when one DEM needs different settings, such as a
   different resolution:

    ```yaml
    dems:
      fine_lidar:
        file: "lidar_2m.tif"
        shared_parameters:
          threshold_contributing_pixels: 250000
    ```

A parameter set to `null` is dropped, so an item can return one parameter
to the tool's default.

The workflow sets the file names and paths (`read_path`, `read_fname`,
`write_path`, `write_fname`), and the parameters that pass one tool's output
to the next. Each module's config file lists them. Setting one of them in
the config stops the run with an error saying so.

The parameter file actually used for each step is written next to its
outputs (for example `results/<run_id>/<dem>/floodplain/fluvial_surfaces.param`),
and each tool also writes `<name>_README.html`, listing every parameter with
the value it ran with.

### Using an installed package

Snakemake runs inside the pixi environment, which brings its own GDAL and
PROJ settings. The workflow clears those before every tool, which is right
for tools built from this checkout. A package needs its own settings back:
point `lsdtt_setenv_script` at the package's `setenv.sh`, and give the tools
by name, since `setenv.sh` puts them on the `PATH`:

```yaml
lsdtt_setenv_script: "/home/me/lsdtt3-backend-linux-x86_64-core-v0.7.1/setenv.sh"
executables:
  fluvial_surfaces: "lsdtt-fluvial-surfaces"
  valley_centreline: "lsdtt-valley-centreline"
```

The package must include every tool the module uses.

---

## Figures

`valley_widths`, `chi_region` and `lithology_channels` can also draw figures,
using the `lsdviz3-*` commands of
[lsdviztools3](https://lsdtopotools.github.io/lsdviztools3/). **The figures
need lsdviztools3 installed** (see [Setting up](#setting-up), step 4); it is a
separate install with a much larger environment (GMT, PyGMT, geopandas), so
the workflows do not depend on it for anything else. They look for it once per run, and the
`figures` block of the config decides what happens:

```yaml
figures:
  enabled: "auto"                 # auto | true | false
  lsdviz3_command_prefix: ""
```

| `enabled` | lsdviztools3 found | lsdviztools3 not found |
|---|---|---|
| `auto` (default) | figures are drawn | figures are skipped; the analysis runs, with a message saying why |
| `true` | figures are drawn | the run stops before starting, saying how to fix it |
| `false` | not looked for; no figures | not looked for; no figures |

"Found" means `lsdviz3-drape --help` runs. If lsdviztools3 is installed in an
environment whose commands are on your `PATH`, nothing needs setting. If it
lives in its own pixi project (a checkout of lsdviztools3 with its
`pixi.toml`), give the command that runs it:

```yaml
figures:
  lsdviz3_command_prefix: "pixi run --manifest-path /path/to/lsdviztools3/pixi.toml"
```

The figure commands run without the workflows environment's GDAL/PROJ
settings and without `lsdtt_setenv_script`, so the lsdviztools3 environment
uses its own. Figures land in `results/<run_id>/<item>/figures/`, and their
logs in `logs/figures_*.log`. Turning figures on later and running again
draws them without rerunning the analysis, except in `valley_widths` and
`chi_region`: there the figures need a hillshade or chi raster that is only
written when figures are on, so the last analysis step runs again once.

---

## Reproducibility

With `stages: report: true` (the default in the shipped configs) each run
writes `workflows/results/<run_id>/reports/<run_id>_RUN_MANIFEST.json`. It
lists the parameter files, inputs and main outputs, each with its SHA-256 hash.

---

## check_basins

Checks that a set of outlet points selects the basins you meant: it fetches a
DEM for a box, fills it, extracts the basin above each outlet, and writes the
basin outlines and channel network as FlatGeobuf for QGIS. Its settings live in
`workflows/config/config.yaml` (the box under `fetch`, the outlets file under
`check_basins.outlets_file`, a CSV with `latitude` and `longitude` or a point
vector file). It is the default module, so it runs without `--configfile`:

```bash
pixi run -e workflows snakemake --cores 1 --snakefile workflows/Snakefile --directory workflows
```

Outputs: `workflows/results/<run_id>/extract/check_junction_basin_outlines.fgb`
and `check_channel_network_lines.fgb`.

---

## Directory layout

| Path | Contents |
|---|---|
| `workflows/Snakefile` | Entry point |
| `workflows/config/` | One configuration file per module; `config.yaml` holds the defaults every run starts from |
| `workflows/rules/` | The rules of each module and the shared stages |
| `workflows/scripts/` | Helper scripts (table gathering, ensemble statistics) |
| `workflows/results/` | Outputs, one folder per `run_id` (not tracked by git) |
