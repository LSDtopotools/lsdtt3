# Hillslope-channel tagging

[`lsdtt-hillslope-channel-tagging`](../cli_reference/cli_reference_lsdtt_hillslope_channel_tagging.md)
moves channel information onto hillslopes. This page describes the method and
how it is tested. The [Python guide](../python/hillslope_channel_tagging.md)
shows the same steps in Python.

## The channel link

Every hillslope pixel is linked to the first channel pixel on its D8 path. This
link does not depend on the ridge method, so HAND and every tag are the same
under all three methods.

$$
\text{relief}(p) = z_{\text{ridge}(p)} - z_{\text{channel}(p)},
\qquad
\text{relative height}(p) = \frac{\text{HAND}(p)}{\text{relief}(p)}
$$

Heights come from the routing surface. Undefined values are nodata.

## Why D8 ridges are shared within a group

D8 flow converges, so about half of the hillslope pixels of a typical DEM have
no donors. Taking each pixel's own farthest upslope pixel would make each of
those its own ridge, with a relative height of exactly 1. `farthest` and
`highest` therefore choose one ridge per group:

- `channel_pixel`: all pixels draining to one channel pixel;
- `entry_pixel`: all pixels entering the channel through the same last
  hillslope pixel, which separates the slopes joining a channel pixel;
- `pixel`: the pixel and those draining through it.

Ties are resolved exactly, with path lengths held as cardinal and diagonal step
counts. `farthest` then prefers the higher pixel, and `highest` the longer path.

`highest` gives relative heights between 0 and 1 on a surface where elevation
does not rise downstream. With `farthest` a pixel on a shorter, steeper branch
can stand above the group's farthest pixel, so its relative height exceeds 1.
This is reported and counted rather than clipped.

## The divide network and aspect rays

Each channel pixel belongs to a junction-to-junction link, with a junction
belonging to the link below it. Each hillslope pixel takes the link its D8 path
reaches. A divide edge is a pixel edge between two valid non-channel pixels with
different labels. No divide is filtered out.

An `aspect_ray` runs straight upslope from the pixel centre, opposite the
downslope aspect relative to grid north. It visits every cell it touches
(Amanatides and Woo 1987), so it cannot pass through a diagonal channel or
divide. Link catchments also meet along lines running upslope from junctions,
where the ground keeps rising past the divide, so a ray stops only at the first
divide edge that satisfies both:

1. its higher pixel is not below the starting pixel;
2. the next cell beyond the crossing is not higher than that pixel. Nodata,
   off-grid and channel cells count as not higher.

The ridge is the higher pixel of that edge. If the ray reaches a channel, nodata
or the grid edge first, the ridge is the highest cell on the ray. A pixel with
no aspect has no ridge.

The crest test looks only one cell ahead, so a small bump on a rising slope can
stop a ray early. The default aspect comes from a polynomial surface on the
input DEM. D-infinity is available and gives noisier ray directions.

## Chi and k_sn

Chi and k_sn have no meaning without a true baselevel. A sea mask and the collar
of a clipped DEM are both nodata, so geometry cannot separate them, and the user
declares the valid basins by outlet file or maximum outlet elevation. Basins
with a pixel on the grid edge are always excluded. Of nested outlets only the
outermost is kept, so each channel pixel has one chi value.

A channel node takes the k_sn of the fitted segment whose chi interval contains
it, within its own branch.

## Imported attributes

Records join by exact pixel when they carry `row`, `col` and `flat`, and
otherwise snap to the nearest channel pixel centre. Lines are sampled every half
cell. Several records on one pixel give their median, since smoothed planform
geometry rarely maps one to one onto channel pixels, and a line has one vote per
pixel.

A channel pixel with no record takes the nearest record on its own link, by D8
channel distance and up to a limit. Equidistant records give their median.
Filled values are never used as sources, and filling never crosses a junction.

## How it is tested

The driver `lsdtt-test-hillslope-channel-tagging` checks:

- exact HAND, relief and relative height on planar valleys;
- ties, ratios above 1 and channel fragments on hand-built D8 paths;
- each ridge method against a brute-force search of its group;
- exact divide edges on hand-built networks;
- rays in every octant against a fine-step reference ray;
- tags, basin selection, attribute joins, reprojection, medians, gap filling and
  every refusal.

CLI tests carry the file handling, and Python tests the bindings.

## References

Amanatides, J. and Woo, A. (1987). A fast voxel traversal algorithm for ray
tracing. *Eurographics* 87, 3-10.

Nobre, A. D., et al. (2011). Height Above the Nearest Drainage: a hydrologically
relevant new terrain model. *Journal of Hydrology* 404, 13-29.
