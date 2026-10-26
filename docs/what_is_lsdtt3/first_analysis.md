# Your first analysis

This page takes you from a DEM sitting on disk to a hillshade and a channel
network. It assumes you have lsdtt3 installed and the tools on your `PATH` —
if not, start with [Installing lsdtt3](install_lsdtt3.md).

The one thing to understand before anything else: **lsdtt3 tools are driven by
a parameter file, not by command-line flags.** There is no `--dem` or
`--output-directory`. You write a short text file that says what to read, where
to write it, and what to compute; then you hand that file to a tool.

---

## The parameter file

A parameter file is plain text, one setting per line, in the form `key: value`:

```
# hillshade.param
# Lines beginning with # are ignored.

read_path: /home/alice/data/
read_fname: my_dem.tif
write_path: /home/alice/output/
write_fname: my_dem

write_hillshade_raster: true
```

These rules cover everything:

- **`key: value`, one per line.** Keys are `snake_case` and never contain
  spaces.
- **`#` starts a comment** anywhere on a line: on its own line, after a value,
  or even in the middle of one (`run#2` reads as `run`).
- **A value must not contain spaces.** There is no quoting mechanism, so a path
  containing a space will not parse; lists are comma-separated with no spaces
  (`25,49,99`). Put your data somewhere without spaces in the path.
- **Values are never quoted.** lsdtt3 keeps quote characters as part of the
  value, so `"my_dem.tif"` would look for a file with quotes in its name; the
  tool stops and tells you to remove them.
- **Booleans are `true` or `false`.** (`True`, `TRUE`, `t`, `T`, `1` and
  `False`, `FALSE`, `f`, `F`, `0` also work.) Anything else, including `yes`
  and `no`, stops the run.
- **Numbers are plain numbers, with no units.** A whole-number parameter takes
  a whole number: `1000`, not `1e3` or `1000.0`. Other numbers can be written
  `0.5`, `2` or `1e-4`.
- **Blank lines are fine.** So is any order — the file is a set of settings,
  not a script. If a key appears twice, the last value is used.

A value the tool cannot read stops the run before anything is computed, with a
message naming the parameter, the value and the fix. Every problem in the file
is reported at once:

```
FATAL: parameter 'n_iterations' value '1.5' is not a whole number. Write a whole number, e.g. 1000.
FATAL: parameter 'methd' is not a parameter of this tool. Did you mean 'method'?
  2 problem(s) in run.param. Fix them and run again; "-h <parameter>" explains any parameter.
```

The file extension does not matter. `.param` and `.driver` are both used in
lsdtt3's own examples and neither is special.

### The four I/O parameters

Every tool reads and writes through the same four parameters:

| Parameter | Meaning |
|---|---|
| `read_path` | Directory the input lives in |
| `read_fname` | Input filename |
| `write_path` | Directory outputs go to |
| `write_fname` | **Stem** for output filenames — no extension |

`read_fname` may include the extension (`my_dem.tif`) and lsdtt3 will split it
off. `write_fname` is a stem only: the tool appends its own descriptor and
extension, so `write_fname: my_dem` produces
`my_dem_hillshade.tif`, `my_dem_fill.tif` and so on. If you leave `write_fname` out it defaults to the
`read_fname` stem, and if you leave `write_path` out it defaults to
`read_path`.

Outputs are Cloud Optimised GeoTIFF by default. `write_output_format` changes
that to `gtiff`, `png` or `same_as_input`.

---

## Running a tool

Two forms work:

```bash
# Give the tool one argument: the parameter file
lsdtt-raster-preprocessing hillshade.param

# The one argument can be a path, relative or absolute
lsdtt-raster-preprocessing /home/alice/params/hillshade.param

# Or two: a directory, then the filename within it
lsdtt-raster-preprocessing /home/alice/params/ hillshade.param
```

If you omit `read_path` from the file itself, it defaults to the directory the
parameter file is in, however you named the file — which is why keeping the
parameter file next to the data is a common and convenient habit.

The only bare flags any lsdtt3 tool accepts are `-h` / `--help`, `-v` /
`--version` and `-c` (citation). Everything else goes in the file.

### How a run ends

A run that finishes ends with one line saying so, how many files it wrote, and
where the run summary is:

```
lsdtt3: lsdtt-raster-preprocessing finished OK; 4 output files; summary: /home/alice/output/my_dem_run_summary.json
```

A run that fails ends with a box on the screen, holding the error:

```
=========================================================
|| RUN FAILED: lsdtt-raster-preprocessing
|| FATAL: parameter 'write_hillshade_rastr' is not a parameter of this tool. Did you mean 'write_hillshade_raster'?
||   1 problem(s) in hillshade.param. Fix them and run again; "-h <parameter>" explains any parameter.
|| Summary: /home/alice/output/my_dem_run_summary.json
=========================================================
```

and with a non-zero exit code, so a script can tell. If `write_path` itself is
the problem, the box says that no summary file could be written.

Every run that knows its `write_path` leaves two files beside its outputs:

- **`<write_fname>_run_summary.json`** — what happened: the tool and lsdtt3
  version, the command line, the resolved paths, start and end times,
  `status` (`ok` or `failed`), the error and warning messages, and every
  output file the run wrote, with its size.
- **`<write_fname>_ingested_params.param`** — every value the run used, the
  ones you set and the defaults you did not. Run the tool on this file to
  repeat the run exactly.

The summary lists files named `<write_fname>_...` that the run wrote, so give
runs that share a `write_path` at the same time different `write_fname`s.

---

## Worked example 1: a hillshade

Make a file called `hillshade.param`:

```
read_path: /home/alice/data/
read_fname: my_dem.tif
write_path: /home/alice/output/
write_fname: my_dem

write_hillshade_raster: true
```

Run it:

```bash
lsdtt-raster-preprocessing hillshade.param
```

You get `/home/alice/output/my_dem_hillshade.tif`. Open it in QGIS, or any image
viewer if you wrote PNG. A hillshade is the fastest way to confirm that your
DEM loaded with the geometry you expected — if it looks inverted or striped,
something is wrong with the data before any analysis starts.

!!! warning "Check your vertical units" (applies to people working in the USA, which inexplicably still uses imperial units)

    A DEM in feet on a metre grid is undetectable from the file alone and
    produces plausible-looking but wrong slopes everywhere downstream. If your
    data might be in feet, set `vertical_units` explicitly rather than hoping.

---

## Worked example 2: a channel network

Channel extraction needs a DEM without pits, so the pipeline is: fill the
pits, then extract. `lsdtt-channel-extraction` does both itself: it fills the
pits and routes the flow, so you just point it at the DEM. If your DEM is
already filled (for example by `lsdtt-raster-preprocessing`), set
`raster_is_filled: true` to skip the fill.

```
# channels.param
read_path: /home/alice/data/
read_fname: my_dem.tif
write_path: /home/alice/output/
write_fname: my_dem

# Where a channel begins, as a count of upslope pixels
threshold_contributing_pixels: 1000

write_channel_network_lines: true
write_drainage_area_raster: true
```

```bash
lsdtt-channel-extraction channels.param
```

`threshold_contributing_pixels` is the parameter you will actually tune. It is
a **pixel count, not an area**, so it depends on your grid resolution: the
drainage area it corresponds to is `threshold_contributing_pixels × cell area`.
On a 10 m DEM the default of 1000 is 100,000 m²; on a 1 m DEM it is 1,000 m².
Lower it for a denser network, raise it for a sparser one.

The channel network comes out as a vector file (FlatGeobuf by default,
`write_vector_format: gpkg` for GeoPackage) that you can load straight into
QGIS on top of your hillshade.

---

## Finding out what a tool can do

Every tool documents its own parameters:

```bash
lsdtt-channel-extraction --help
```

This prints every parameter the tool declares, with its type, default and
description, grouped by theme. It is generated from the same declarations the
tool parses, so it cannot drift out of date.

The same content is on this site under
[Command line reference](../cli_reference/cli_reference_lsdtt_channel_extraction.md),
which is easier to search.

To see what a tool actually used at runtime — including defaults you never set
— look at the `_ingested_params.param` file it writes next to your outputs.
That file is a complete, re-runnable record of the run, which makes it the
right thing to keep for reproducibility. (A default that means "not set", such
as `nan`, is written as a comment: leaving the key out gives the same default.)

`lsdtt-channel-extraction -h json` prints the same parameter list as a single
JSON document and nothing else, for programs that read it.

---

## When something goes wrong

**"No parameter file specified."** You ran the tool with no arguments.

**The run stopped with "Did you mean ...?".** A key in your file looks like a
misspelling of one of the tool's parameters (the same apart from capitals, or a
letter or two out). Correct it; running on would have silently used the
default for the parameter you meant.

**A parameter you set had no effect.** Most likely it is not a parameter of
that tool. A key that is not close to any of the tool's parameters is reported
as a warning, not stopped on — read the warnings the tool prints on startup,
or the `warnings` list in the run summary, and check the key against `--help`.

**The tool cannot find your DEM.** `read_path` and `read_fname` are
concatenated, so a missing trailing slash on the path is the usual cause,
though lsdtt3 tries to repair that. Check the `lsdtt_param_parser: read :`
line the tool prints — it shows the full path it resolved.

**Outputs are not where you expected.** `write_path` defaults to `read_path`,
not to the current directory.

**"has a value containing spaces".** Values cannot contain spaces. Remove the
spaces from a comma-separated list (`a,b,c`), or move the data to a path
without spaces.

---

## Writing parameter files from a script or a coding agent

Everything above applies; a few points matter more when a program writes the
file and reads the result:

- Write one `key: value` line per setting, with booleans as `true`/`false`
  and numbers as plain numbers. In Python, `f"{key}: {value}"` is enough;
  write `true`/`false` rather than Python's `True`/`False` for clarity (both
  are accepted).
- The format looks like YAML but is not YAML: there is no quoting, nesting or
  list syntax, and a YAML library will read some values differently (`yes`,
  dates, `1e-4`, empty values). If you read a parameter file with PyYAML, use
  `yaml.load(f, Loader=yaml.BaseLoader)`, which returns every value as a
  string, exactly as lsdtt3 sees it.
- Run the tool and check its exit code; then read
  `<write_fname>_run_summary.json` for the status, the messages and the list of
  output files, rather than parsing the screen output.
- `<tool> -h json` gives every parameter's name, type, default and description
  as JSON on stdout, with nothing else printed.

---

## Where to go next

- [What lsdtt3 does](overview.md) — the capabilities, and what it is not
- [Raster data](raster_data.md) — supported formats, and the conventions on
  nodata and projections
- [Vector data](vector_data.md) — the vector output formats and the CSV path
- [Command line reference](../cli_reference/cli_reference_lsdtt_raster_preprocessing.md)
  — every parameter of every tool
- [Design and algorithms](../design_and_algorithms/surface_metrics.md) — what
  the numbers actually mean, once you are producing them
