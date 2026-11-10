# Oblique analysis frames

This page describes how lsdtt3 builds the oblique frames used to analyse a DEM
along a feature, how rasters are warped into a frame and back, what happens to
bearings, and how far results agree between frames and UTM. Usage is in
[Analysing along a feature](../what_is_lsdtt3/oblique_frames.md) and
[Oblique frames in Python](../python/oblique_frames.md).

The code is in `src/core/lsdtt_oblique_frame.{hpp,cpp}` (frame geometry and the
feature-file reader), `src/core/lsdtt_rotation_frame_file.{hpp,cpp}` (the frame
file), `src/core/lsdtt_raster_warp.{hpp,cpp}` (warping into a frame and back),
`src/core/lsdtt_raster_fetch.cpp` (fetching into a frame) and
`src/core/lsdtt_true_north.{hpp,cpp}` (grid convergence and true-north
bearings).

## The projection

A frame is a Hotine Oblique Mercator projection, variant B (PROJ `omerc`), whose
central line runs along the feature:

```
+proj=omerc +lat_0=<centre lat> +lonc=<centre lon> +alpha=<azimuth> +gamma=90
            +k=1 +x_0=0 +y_0=0 +<datum of the source> +units=m +no_defs
```

- **Feature along +x.** The axis azimuth is folded into (−90°, 90°] (a line has
  no direction), and `gamma=90` turns the rectified grid so the central line
  is the x axis. Without the fold, an azimuth of 135° put the line on the y
  axis. Azimuths within 1e-4° above −90° are taken as +90°: PROJ places the
  frame centre half a meridian away for alpha within about 1e-6° of −90.
- **Centre at the origin.** The frame centre (the geodesic midpoint of the axis
  end points) maps to (0, 0), so the extent is symmetric:
  x in [−L/2 − b, L/2 + b], y in [−W/2, W/2] for axis length L, end buffer b
  and width W (b also on the long sides for a polygon or a data footprint),
  snapped outward to whole cells.
- **Scale.** The scale factor is 1 on the central line and grows as about
  d²/(2R²) at distance d from it: 3 × 10⁻⁵ at 50 km, against 4 × 10⁻⁴ on a UTM
  zone's central meridian.
- **Datum.** The frame is built on the geographic CRS of its source: WGS84 for
  a fetch (Copernicus GLO-30 is WGS84), the source DEM's own datum when an
  existing DEM is rotated, so no datum transformation enters the warp.
- **Self-check.** Every frame is checked when built and when read: its centre
  must map to (0, 0) within 1 mm, and a point placed along the azimuth
  (`geod_direct` through an azimuthal equidistant projection) must land on +x
  within 0.01 %.

A raster in a frame is an ordinary north-up raster in that CRS, so the DEM
loader, the flow routing and every analysis read it unchanged.

## Building a frame from a feature

| Feature | Axis | Width |
|---|---|---|
| Line or multi-part line | the two furthest vertices of all parts | given (`*_width`) |
| Two points | the two points | given |
| Polygon | longest side of the minimum-area enclosing rectangle | the rectangle's short side |
| Valid data of a raster | as a polygon, from the outer corners of the first and last valid cell of every row | the rectangle's short side |

The end buffer defaults to ten cells. A polygon whose rectangle has
long/short below 1.1 is warned as having an ill-defined axis. Geodesics
(lengths, azimuths, midpoints) are ellipsoidal on the frame's datum.

## Warping into a frame and back

**Into a frame.** Fetching warps the source tiles straight into the frame; a
DEM on disk is warped from its own CRS. Either way the DEM is resampled once.
The cell size is `grid_spacing` if given; otherwise a DEM projected in metres
keeps its own, and any other DEM gets the size reprojection would choose.

**Rotate-back target.** The frame file records the grid that products go back
onto. After rotating a DEM projected in metres, it is that DEM's own CRS, cell
size and lower-left corner, so products land on its pixels (a whole number of
cells away). After a fetch, or for a geographic source, it is the UTM (or
NSIDC polar stereographic beyond 75°) CRS of the frame centre on the 0, 0
lattice that plain fetches use. A raster at the frame's cell size goes to the
target cell size; any other keeps its own.

**Fixed resampling kernel.** GDAL sizes its resampling kernel per warp chunk
from the ratio of the chunk's source bounding box to its destination size. A
rotated grid's source box is larger than the chunk, so an equal-cell rotation
is treated as downsampling and the kernel is widened by an amount that depends
on the chunk size, the rotation and the output driver. Frame warps therefore
fix the kernel scale (GDAL's `XSCALE`/`YSCALE`: 1, or source/target pixel size
when the frame cell is coarser) and use the exact transform (`-et 0`). For a
30 m Copernicus DEM rotated by 36° (5.6 million cells), the warp then agrees
with a fixed-scale exact `gdalwarp` to 0 m whatever the warp memory or output
format, where the default kernel had differed from it by 0.87 m median and up
to 55 m on cliffs.

**Round trips.** Warping a DEM into a frame and back (cubic both ways) on that
DEM changed elevations by a mean of 0.0002 m, 0.36 m RMS, 1.3 m at the 99th
percentile and 18 m at most; on a smooth synthetic DEM, by 0.0007 m median and
0.22 m at most. Categorical rasters rotated with `near` keep their values: 184
of 184 IDs survived a round trip with no new values.

## The rotation frame file

`<name>_rotation_frame.txt` holds one `key: value` per line, written so the
file is also a flat YAML mapping (`yaml.safe_load` reads it): text values
single-quoted, WKT on one line, numbers with 17 significant digits and always
a decimal point, booleans as `true`/`false`. It records provenance
(`frame_file_version`, `lsdtt3_version`, `created_utc`, `created_by`,
`rotated_dem_fname`), the frame (`centre_lon`, `centre_lat`, `azimuth_deg`,
`cell_size`, `frame_x_min` … `frame_y_max`, `base_crs_wkt`, `frame_crs_wkt`),
how it was defined (`frame_source`, `feature_fname`, `feature_layer`, the axis
end points, `axis_length_m`, `width_m`, `end_buffer_m`, `buffer_all_sides`) and
the rotate-back target (`target_crs_epsg` or `target_crs_wkt`,
`target_cell_size`, `target_origin_x`, `target_origin_y`, `source_dem_fname`).

It has its own reader rather than the parameter-file parser, which ends a value
at the first space and caps a line at 1023 characters. Reading refuses a
missing or newer version, duplicate or missing keys, a `frame_crs_wkt` that
does not match the stored centre and azimuth, and a frame that fails the axis
self-check; it warns on unknown keys and when the frame rebuilt from the stored
definition differs from the stored one.

Without the file, a frame is rebuilt from a raster's own CRS and extent, and
the target falls back to the analysis CRS of the frame centre on the 0, 0
lattice, with a warning: an unsnapped source origin cannot be recovered.

## Bearings and grid convergence

Every bearing lsdtt3 writes is clockwise from grid up. In a frame with the
feature along +x, 90° is along the feature. The grid convergence c(x, y) is the
angle of true north measured clockwise from grid up, computed by stepping
1e-4° north through the CRS's own geographic CRS; a grid bearing g is the true
bearing g − c. At a frame's centre c = 90° − azimuth; on UTM, |c| is at most
about 3°.

With `bearings_relative_to_true_north`, each output is converted at its own
location: a raster cell by cell (c on a lattice every 32 cells, interpolated
bilinearly: within 2e-6° of a per-cell conversion on a 45 km × 12 km frame), a
point output at its point, a bearing between two points at their midpoint, a
basin or patch axis at its centroid, and a whole-raster 2D spectrum with c at
the grid centre (written in the CSV header). Angles between two directions are
not converted. The D8 diagnostic `fraction_on_d8_multiples` in the
channel-bearing statistics is always taken from grid bearings. Hillshade
illumination is never converted.

Measured against pyproj: c to within 1e-3° in UTM north and south, NSIDC polar
stereographic and frames; a plane dipping to true azimuth 200° gives a
converted polyfit aspect of 200.000° in UTM and in frames at azimuths −60°, 0°,
30° and 89°.

## Agreement between frames and UTM

Analytic surfaces were defined in UTM 12N metres at 38.4° N, −110° (600 × 600
cells of 10 m, centred off every lattice node) and analysed on the UTM grid and
on frames at azimuths 0°, 15°, 30°, 45°, 60°, 75° and 89.9°:

| Quantity | UTM | Frames, all azimuths |
|---|---|---|
| Plane slope (true 0.2), median | 0.2 | 0.199939 |
| Paraboloid Laplacian (true −1e-3 /m), median | −1.000e-3 | −0.99939e-3 |
| Funnel D8 outlet area (km²) | 19.6357 | 19.6464 to 19.6472 |
| Funnel D-infinity and MFD outlet area (km²) | 19.4361 | 19.4468 to 19.4476 |
| Funnel D8 path length from the rim / radius | 1.0548 | 1.0550 to 1.0558 |
| V-valley D8 channel length / straight length | 1.0008 | 1.005 (0°, 45°, 89.9°); 1.073 to 1.076 (15°, 30°, 60°, 75°) |

Frames agree with each other to better than 4 × 10⁻⁵ in slope, curvature and
area. The offset between frames and UTM (−3.1 × 10⁻⁴ in slope, +6.1 × 10⁻⁴ in
area) is the UTM scale factor 1° from the zone's central meridian
(k = 0.99969), not the rotation. D8 channel lengths are not invariant: they
depend on the channel's angle to the grid (+0.5 % along an axis or diagonal,
+7.3 to 7.6 % at 15° off), the known D8 path bias, which is the same on UTM for
a channel at that angle.

## Storage

For a 100 km × 10 km strip at 45° on the Colorado Plateau (Copernicus 30 m),
the frame holds 3354 × 334 = 1,120,236 cells (GeoTIFF 3.7 MB); the same strip
on a UTM grid needs 2629 × 2604 = 6,845,916 cells (22.0 MB), 6.1 times as many.
Channel extraction on the frame peaked at 143 MB resident memory, against
618 MB on the UTM grid.

## Figures

[lsdviztools3](https://lsdtopotools.github.io/lsdviztools3/) draws a frame raster in the Cartesian `X` projection in frame km,
so the strip is a rectangle. GMT's own oblique Mercator (`-JO`) is not used: it
is spherical. `mapproject -Joa` gives the same coordinates whatever ellipsoid
GMT is told to use, and matches PROJ's `omerc` on a sphere to the centimetre;
against the ellipsoidal frames lsdtt3 writes, a point 50 km along a 30° axis
at 38.4° N is 33 m long and 89 m off the axis (177 m at 100 km). The north
arrow is turned by the grid convergence at its own position, and the graticule
is the reference CRS's easting/northing grid transformed into the frame.
