# lsdtt-fluvial-surfaces

Sorts flat ground into fluvial levels — floodplain, terrace treads, fan lobes —
by banding a height reference and labelling connected components within each
band.

---

## Where this sits

Three tools, each doing one job:

| tool | question | needs |
|---|---|---|
| `lsdtt-flat-surface-finder` | which ground is flat? | a DEM |
| **`lsdtt-fluvial-surfaces`** | **which fluvial level is each flat surface on?** | a channel, or a supplied height field |
| `lsdtt-valley-centreline` | trace a centreline through this valley | a valley mask |

The split is at the point where a channel network becomes necessary.
`lsdtt-flat-surface-finder` is pure morphometry: no flow routing, no channel, no
height reference. This tool does everything that needs one, and writes the
valley mask `lsdtt-valley-centreline` consumes.

The normal two-step workflow:

```
lsdtt-flat-surface-finder  →  <stem>_combined_flat_mask.tif
                                        ↓  flat_mask_fname
lsdtt-fluvial-surfaces     →  <stem>_valley_mask.tif, _band_index.tif, …
                                        ↓  valley_mask_fname
lsdtt-valley-centreline    →  <stem>_valley_centreline.fgb
```

You can run this tool on its own — leave `flat_mask_fname` empty and it falls
back to a built-in relief or slope admission. That is for a first look, not for
a result you intend to publish; the texture metrics in
`lsdtt-flat-surface-finder` are much better at deciding what is flat.

---

## The idea

### Band first, then label

The obvious way to segment flat surfaces is to label connected components and
then work out how high each one is. That forces the riser network between
terraces to be **continuous**: a road, a gully or a stock track cutting through
a riser lets two treads at different heights merge into a single component, and
no amount of repair afterwards can separate them again.

This tool inverts the order. It bands the height reference first, then labels
components **within each band**. Two pixels can only share a label if they share
a band, and two treads at different heights are by construction in different
bands — so the merge cannot happen. The riser never has to be found, and never
has to be continuous.

That is the whole reason the tool exists. Everything else is in service of it.

### There is no valley mode and no terrace mode

Finding a valley floor and segmenting a terrace flight are the same computation
at different band counts. The valley floor is the lowest band of a two-band cut;
a terrace flight is an n-band ladder on the same population. So there is one
parameter, `ladder_level_count`:

- `2` — the valley cut. Band 0 is the valley floor.
- `n` — an n-level ladder.
- `0` (default) — automatic: the largest level count the DEM can actually
  resolve, see below.

### Which reference field

**`hand`** (default) — height above the in-scope channel network. Terraces along
a reach sit at a roughly constant height above the channel, so HAND makes a
terrace flight a step function. This is the field for fluvial terraces and
floodplains.

**`fitted_surface`** — the signed residual to a surface fitted to the flat
ground. A fan is *not* a step function in HAND: it climbs continuously from toe
to apex, so a HAND ladder cuts one surface into arcs banded by distance from the
apex rather than by age. What *is* roughly constant on a fan is the depth of an
inset lobe below the projection of the older surface — which is what the
residual measures.

Choose the fitted surface with `fitted_surface_algorithm_choice`:

| value | when |
|---|---|
| `cone` (default) | **a fan or a debris cone.** z falls linearly with distance from an apex. |
| `plane_ransac` | a piedmont, a terrace tread, a broadly planar surface. RANSAC because incised channels and inset lobes are large-amplitude outliers that would tilt an OLS fit. |
| `poly2` | a broadly curved surface with no single centre. |

**Use `cone` for a fan.** A fan is radial, and a cone is not a low-order
polynomial in x and y. Measured on a synthetic fan (2 km radius, 0.03 gradient,
lobes inset 4 m), the spread of the older surface's own residual — which should
be flat, so its spread is pure misfit — was about 34 m under a RANSAC plane and
10 m under a poly2, against a 4 m signal. Under a cone it was 0.3 m.

The cone's apex is found automatically as the centroid of the highest 0.1 % of
fit-mask cells, or you can state it with `fitted_surface_apex_x` /
`_apex_y`. It is supplied rather than solved for because it is a strongly
non-linear parameter and you can usually point at a fan head more reliably than
a solver can find it.

### How many levels, and when to believe them

A level count is **admissible** when the tightest gap between consecutive band
centres clears `ladder_resolvability_kappa * vertical_accuracy_metres`. The
answer is the largest admissible count.

`vertical_accuracy_metres` is metadata you must supply — it cannot be derived
from the DEM. Leaving it at 0 disables the check, and then **every window
reports `resolved` and the ladder carries no instrument guarantee at all**. The
tool warns when you do. (Note that a DEM whose vertical units are feet on a
metre grid has an effective accuracy 3.28 times larger than it looks.)

Where no count qualifies, the window reports `unresolved` rather than returning
a number. **Read `<stem>_ladder_windows.csv` before describing your levels**:
whether they may be called resolved levels at all depends on which branch the
ladder took, and that is not visible in any raster.

An explicit `ladder_level_count` is honoured and reported honestly — if that
partition is not admissible the window still records the count you asked for,
with status `unresolved`, rather than silently substituting a different one.

---

## Composability: `hand_fname`

HAND is the expensive stage, and the level count is what you iterate on. So:

```
# once
write_hand_raster: true

# then, as often as you like
hand_fname: <stem>_hand.tif
ladder_level_count: 5
```

The supplied-HAND path produces **bit-identical** bands and labels to the
computed path, which is pinned by a test.

One limit worth knowing: reach-relative windowing needs each pixel's HAND
*donor node*, and a raster cannot carry those. A supplied HAND is therefore
always banded globally. Iterate on the level count with a supplied HAND, then do
the final run computing HAND if you want the windowed ladder.

---

## Outputs

| file | what |
|---|---|
| `_band_index.tif` | Int16, nodata −1. Which level each pixel is on. Authoritative. |
| `_banded_labels.tif` | UInt32. Patch label, a global running counter 1..N; 0 is background. |
| `_valley_mask.tif` | Band 0. The input to `lsdtt-valley-centreline`. |
| `_surface_class.tif` | floodplain / terrace class. `hand` field only. |
| `_surface_patches.csv` | One row per patch: label, band, class, source_id, n_pixels, area, median field value, aspect. **The join key for everything downstream.** |
| `_ladder_windows.csv` | The resolvability report. Read it first. |
| `_reference_field.tif` | The field that was actually banded, whichever it was. |
| `_hand.tif` | The HAND raster, for feeding back as `hand_fname`. |

The band is **not** packed into the label. The label is a plain running counter
and the band travels in `_band_index.tif` and in the patch table, so nothing has
to divide a label to learn its level.

`_surface_class.tif` and `_surface_class.csv` are written only for the `hand`
field. Floodplain and terrace are statements about height above a channel; they
mean nothing on a fitted-surface residual, so no class is invented for one.

---

## Things that will bite you

**`antimode` cannot band a fitted surface.** The antimode break finder builds its
histogram on `log10` of the field, which suits HAND (non-negative by
construction) and destroys a signed one: every negative residual is floored into
the lowest bin, so the entire below-surface population — the inset lobes, the
thing you are looking for — collapses into a single spike. Selecting
`fitted_surface` switches the break finder to `jenks` automatically if you did
not name one, and rejects an explicit `antimode` with that explanation.

**A terrace flight often wants `jenks`, not the default `antimode`.** A flight is
spaced roughly linearly in height above the channel, and the log10 binning
compresses the upper treads until they merge.

**`min_component_pixels` is per band here**, unlike in
`lsdtt-flat-surface-finder` where the same name filters components of the whole
mask. Here it is a floor on the size of one surface on one level.

**Cleaning the valley mask can make it bigger.** `clean_valley_mask` erodes,
drops patches that do not reach the channel, then regrows — and the regrow fills
enclosed holes. It is not a filter.

**The tool needs flow routing unless you supply `hand_fname`**, so the
single-command case is not cheap. That is what the two-step workflow above is
for.

---

## Full parameter list

Run `lsdtt-fluvial-surfaces -h` for every parameter with its type, default and
description, grouped in pipeline order. `-h html` writes a browsable reference.
Every run also writes `<write_fname>_README.csv` and `.html` alongside its
outputs, recording exactly what it was given.
