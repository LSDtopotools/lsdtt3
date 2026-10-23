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

Four rules cover almost everything:

- **`key: value`, one per line.** Keys are `snake_case` and never contain
  spaces.
- **`#` starts a comment**, either on its own line or after a value.
- **A value ends at the first space.** There is no quoting mechanism, so a path
  containing a space will not parse. Put your data somewhere without spaces in
  the path.
- **Blank lines are fine.** So is any order — the file is a set of settings,
  not a script.

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

# Or two: a directory, then the filename within it
lsdtt-raster-preprocessing /home/alice/params/ hillshade.param
```

If you omit `read_path` from the file itself, it defaults to the directory the
parameter file is in — which is why keeping the parameter file next to the data
is a common and convenient habit.

The only bare flags any lsdtt3 tool accepts are `-h` / `--help`, `-v` /
`--version` and `-c` (citation). Everything else goes in the file.

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
— look at the parameter echo it prints on exit, or the `_ingested_params.param`
file it writes next to your outputs. That file is a complete, re-runnable
record of the run, which makes it the right thing to keep for reproducibility.

---

## When something goes wrong

**"No parameter file specified."** You ran the tool with no arguments.

**A parameter you set had no effect.** Most likely it is not a parameter of
that tool. By default unknown keys are reported rather than silently ignored —
read the warnings the tool prints on startup. Check the spelling against
`--help`; a key that looks right but is not declared does nothing.

**The tool cannot find your DEM.** `read_path` and `read_fname` are
concatenated, so a missing trailing slash on the path is the usual cause,
though lsdtt3 tries to repair that. Check the `lsdtt_param_parser: read :`
line the tool prints — it shows the full path it resolved.

**Outputs are not where you expected.** `write_path` defaults to `read_path`,
not to the current directory.

**A value with a space in it was truncated.** Expected: values end at the first
whitespace. Move the data to a path without spaces.

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
