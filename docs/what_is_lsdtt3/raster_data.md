# Raster Data

Rasters are the main data model in lsdtt3. A raster combines a regular grid of
values with its dimensions, cell size, origin, coordinate reference system
(CRS), data type, and nodata value.

## Input Formats

lsdtt3 reads raster formats supported by GDAL. GeoTIFF is the usual input
format, including Cloud Optimized GeoTIFF (COG). Input filenames are supplied
with `read_path` and `read_fname`; individual tools document any additional
raster inputs they require.

## Coordinate Reference Systems

Most geomorphic analyses measure lengths, areas, gradients, or search radii and
therefore require a projected CRS with horizontal units in metres. A geographic
latitude/longitude raster should normally be reprojected before analysis. The
`lsdtt-raster-preprocessing` tool can inspect, reproject, crop, and prepare a DEM.

Synthetic rasters receive a projected default georeference unless the generating
tool provides another CRS option.

## Nodata and Masks

Nodata cells are excluded from analysis. Internally, lsdtt3 keeps raster values
and a validity mask together; output nodata cells carry the output raster's
declared nodata value. Be careful around raster edges and interior voids because
many flow-routing and neighbourhood operations depend on which cells are valid.

`lsdtt-raster-preprocessing` provides tools for filling interior nodata holes,
masking values, and preparing hydrologically suitable surfaces.

## Output Formats

Most raster-producing CLIs use `write_output_format`:

| Value | Result |
|---|---|
| `cog` | Cloud Optimized GeoTIFF; the standard default |
| `gtiff` | Conventional GeoTIFF |
| `same_as_input` | Use the input raster format where supported |
| `png` | PNG plus georeferencing sidecars where supported |

Output filenames use `write_path` and `write_fname`, followed by a descriptor
specific to the result. The generated CLI reference gives the exact outputs and
format support for each executable.
