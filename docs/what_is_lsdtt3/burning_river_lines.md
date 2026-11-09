# Burning mapped rivers into a DEM

A DEM and a mapped river network often disagree. A bridge or culvert puts a
dam across the channel. A road embankment blocks a valley. A river runs
through a tunnel. On a coarse DEM, the lowest cells may follow a different
route from the river.

Routing flow on such a DEM sends water the wrong way, and every channel,
basin and drainage-area result built on it inherits the error.
[`lsdtt-stream-burning`](../cli_reference/cli_reference_lsdtt_stream_burning.md)
conditions the DEM so that flow follows the mapped network.

This page covers `stream_burn_mode: burn_to_corridor`, the mode for when
**the vector network is the truth**. Use it when you trust where the mapped
river goes more than you trust the DEM's channel, for example a national river
network on a 5–30 m DEM.

---

## Two situations, one setting

How much to trust the DEM is a continuum:

- **The DEM is the truth** (for example a lidar DTM with a detailed water
  network). The DEM's channels are in the right place. Only short
  obstructions, such as bridges and culverts, need cutting, and the vector
  supplies the network's topology. The other modes of `lsdtt-stream-burning`
  (`snap_to_accumulation`, `snap_to_network`) are built for this.
- **The vector is the truth** (a coarser DEM with a national river network).
  Long stretches of channel must follow the mapped line, wherever the DEM's
  lowest cells go. Use `burn_to_corridor`.

Within `burn_to_corridor`, one parameter sets where on that continuum you are:
`vector_corridor_half_width_metres` (X).

- **X = 0** (the default) burns the line's own cells. The channel is exactly
  where the map says.
- **X > 0** lets the channel move up to X metres from the line, to the
  lowest route through the valley floor within that corridor. That suits a
  map line that is a few cells off the valley floor.

## Preparing the river lines

The lines must be connected and oriented downstream.
[`lsdtt-fetch-river-lines`](../cli_reference/cli_reference_lsdtt_fetch_river_lines.md)
does this for a file you already have (`river_line_source: local_file`, with
your field names), or by downloading OS Open Rivers for Great Britain
(`gb_os_open_rivers`, no key needed). It clips the lines to your area, joins
their ends, votes each network's flow direction, and writes
`<write_fname>_river_lines.fgb`. That file is ready to use as
`stream_vector_path`. Use the DEM as the area (`bounds_fname`), so the lines
come out in the DEM's coordinate system.

## A parameter file

```text
# burn.param
read_path: /data/my_site/
read_fname: my_site_dem.tif
write_path: /data/my_site/results/
write_fname: my_site

stream_vector_path: /data/my_site/results/my_site_river_lines.fgb
stream_burn_mode: burn_to_corridor
run_stream_burn: true

# how far (m) the channel may move from the mapped line
vector_corridor_half_width_metres: 0
# the river names, for the outputs and the profile plots
vector_name_field: name

write_burn_depth_raster: true
```

Run it with `lsdtt-stream-burning burn.param`.

## What the burn does

For every path through the network (main stems first, longest first):

1. **Small pits are filled for sampling.** DEM noise makes small pits that
   are not real. Every depression of fewer than `burn_profile_max_pit_pixels`
   cells (default 16) is filled to its spill level. The filled DEM is used
   only to choose and sample the channel. A large real depression, such as
   the valley behind an embankment, is left alone, while noise pits inside
   it are filled.

   Because the count is in cells, scale it with resolution: 16 cells is
   14,400 m² at 30 m but 400 m² at 5 m. 0 gives the pure minima profile.

2. **The channel is chosen.** At X = 0 this is the line's own cells. At
   X > 0 it is the lowest-cost route through the corridor on the
   pit-filled DEM.

3. **The profile joins the minima.** Walking downstream, every cell lower
   than everything upstream of it on the path is a minimum. The profile
   joins successive minima, so it is never above the DEM under the line. It
   only ever cuts down, through embankments and other obstructions, and it
   falls downstream by at least `fill_min_slope` per metre. The minima are
   joined:
   - in a straight line against distance (`burn_profile_fit_algorithm_choice:
     linear`, the default); or
   - in χ (`chi`), with concavity `m_over_n` (default 0.45) and reference
     area `chi_A0`. This gives the concave-up shape of a graded river between
     the same minima. The drainage area for χ is taken from the burned
     network, not from the unburned DEM.

4. **Deep cuts are capped.** A minimum that would make the profile cut more
   than `burn_profile_max_cut_metres` (default 2000 m, so effectively off)
   below the DEM is not used.

   The cap guards against a DEM artefact, such as a deep hole, dragging a
   long stretch of river down. It cannot tell an obstruction (an embankment
   the river does pass through) from a hill the mapped line wrongly crosses.
   Check large cuts in the burn depth before trusting them.

5. **Junctions.** A main stem is never changed by a tributary. A tributary is
   joined to the main stem's elevation where they meet, and raised near the
   junction if its own floor is lower.

6. **The rest of the DEM drains to the channels.** The channel cells are
   fixed. Everything else is breached and then filled, with the channels and
   the DEM edge as outlets, so each hillslope drains to the network.

   **Banks:** a neighbour lower than a channel cell, which would take the
   channel's flow (an offset line at X = 0), is raised into a bank.

   **Ends inside the DEM:** a channel that ends inside the DEM, neither at
   the edge nor on another channel, is given an exit route. This is the
   cheapest route, lowered to drain to the DEM edge, to NoData, or to a lower
   channel.

7. **Final conditioning.** The whole DEM, channels included, is conditioned
   once more as the other modes do: breached if `carve_before_fill` is true
   (`breach_algorithm_choice`), then filled with `fill_min_slope`. So flow
   routing, drainage area and channel extraction on the output all follow the
   burned network. This pass should change nothing, and the log says how
   much it changed.

**Lines crossing NoData.** At X = 0 a channel cannot cross a gap of NoData,
for example a strip of removed sea. Such a line is cut into its pieces over
valid cells:

- each piece is burned on its own;
- the piece before a gap ends at the gap (it drains into the NoData);
- the piece after starts as a new source.

At X > 0 the corridor usually routes round a small gap, and the line is not
cut.

**The sea.** Many DEMs store the sea as valid cells near 0 m, so a river
reaching the coast ends "inside" the DEM and is given an exit route across the
sea.

`min_elevation` turns every cell at or below a value into NoData before the
burn. NoData is an outlet, so the river then ends at the coast. It is off by
default (-9999) in this tool.

Do not set it to 0 on a DEM with real land below sea level (reclaimed land,
fens, a DTM of a coast with a large tidal range). That land would be removed
and its rivers lost. Either set it just below the lowest land, or leave it off
and remove the sea from a vector instead, before burning, with
[`lsdtt-raster-preprocessing`](../cli_reference/cli_reference_lsdtt_raster_preprocessing.md):

```text
min_elevation: -9999
sea_mask_fname: sea_below_low_water.fgb
sea_mask_side_choice: auto
write_sea_mask_raster: true
```

`sea_mask_fname` takes sea (or land) polygons, or a coastline as lines that
need not be closed. `sea_mask_side_choice` says which side is sea:

- `sea_inside` or `land_inside` for polygons;
- `sea_left` or `sea_right` for lines, looking along the direction they were
  digitised;
- `auto` decides from the DEM itself. The side already NoData much more
  often, or else the side with the lower median elevation, is sea, and the
  log gives the evidence. It refuses when the two sides look alike.

For lines, the sea is what a flood from the sea side reaches without crossing
them. If the flood leaks to the land side, because the line has a gap or does
not cross the DEM, only the sea side within `sea_mask_max_distance_metres` is
removed, and the log says so. The removed sea is always an outlet: an interior
NoData fill never fills it back in.

In Great Britain, the sea below the low water mark can be built from the OS
OpenMap Local or VectorMap District tiles (OS OpenData, keyless) as their
`TidalWater` polygons minus their `Foreshore` polygons. The low water mark
there is mean low water in England and Wales, and mean low water springs in
Scotland.

## Checking a run

The log reports what was done. Two of its lines should always read 0:

```text
Final conditioning (carve_before_fill false) changed 0 line cell(s) (largest change 0 m) and 0 other cell(s)
Downhill check, final DEM: 0 of 111 path(s) rise anywhere
```

Other lines report:
- how many reaches were fitted or skipped, and why;
- pits filled, minima used and minima rejected by the cut cap;
- cells cut and raised on the line;
- channel ends drained;
- bank cells;
- cells whose flow still leaves the line.

## Outputs

| File | What it holds |
|---|---|
| `<write_fname>_burned.tif` | The conditioned DEM: route flow on this. |
| `<write_fname>_burned_profiles.csv` | One row per channel cell, downstream along each path: `path_id`, `path_rank` (0 = the longest), `name`, distances, `z_original`, `z_cleaned` (pits filled), `z_profile` (the fit), `z_final`, `pit_filled`, the minima flags, `drainage_area_m2` and `chi` (χ fit), and how the path joins others. |
| `<write_fname>_corridor_channels.fgb` | One line per reach: the channel used (or the mapped line if skipped), `status`, `skip_reason`, the largest cut and raise, and `line_index` / `piece` / `n_pieces` for lines cut at NoData. |
| `<write_fname>_exit_routes.fgb` | One line per channel end inside the DEM: the route lowered to drain it. `exit_type` is `edge`, `line` (to another channel) or `exit_route` (to another end's route), and `max_cut_m` is its deepest cut. A long or deep route is worth a look: it may mean the line stops short, or the DEM holds the sea. |
| `<write_fname>_cross_path_receivers.fgb` | Points: channel cells whose steepest descent goes to another channel's cell, which could not be raised without making that channel rise. The flow leaves the line there. This typically happens where a tributary runs beside the river it joins, one cell apart. |
| `<write_fname>_burn_depth.tif` | (`write_burn_depth_raster`) The input minus the output DEM (m): positive where lowered, negative where raised (filled pits, banks), NoData where unchanged. |

All vector outputs use `write_vector_format`.

## Plotting the profiles

[lsdviztools3](https://github.com/LSDtopotools/lsdviztools3) draws the
burned long profiles from the CSV. It shows:

- the input DEM and the burned channel;
- the minima used and rejected;
- the filled pits;
- the cut and raised stretches;
- where tributaries join.

One run per panel, which makes it easy to compare runs at different X or
`burn_profile_max_pit_pixels`:

```bash
lsdviz3-burned-profiles --csv results/my_site_burned_profiles.csv \
    --name "River Name" --out figures/my_site_profiles.png
```

## From Python

The same routines, with the same parameter names:

```python
import geopandas as gpd
import lsdtt3

lines = lsdtt3.fetch_river_lines(
    river_line_source="gb_os_open_rivers", bounds_fname="my_site_dem.tif").lines
dem = lsdtt3.Raster.read("my_site_dem.tif")
result = lsdtt3.burn_to_corridor(dem, lines, vector_name_field="name",
                                 vector_corridor_half_width_metres=0.0)
burned = result.dem                 # route flow on this
profiles = result.profiles          # the CSV's columns
print(result.report["n_paths_rising_final"], result.report["final_other_cells_changed"])
```

`result` also holds `burned` (before the final conditioning), `burn_depth`,
`channels`, `exit_routes` and `cross_path_receivers`.
