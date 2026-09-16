# Raster Data

Rasters are the main data model in lsdtt3. A raster combines a regular grid of
values with its dimensions, cell size, origin, coordinate reference system
(CRS), data type, and nodata value. I'm not sure why we are telling you this because it is very difficult to imagine you looking for this software if you don't know what a raster is. But better safe than sorry, we guess. 

## Input Formats

lsdtt3 reads raster formats supported by [GDAL](https://gdal.org/en/stable/). GeoTIFF is the usual input
format, including Cloud Optimized GeoTIFF (COG). Input filenames are supplied
with the parameters `read_path` and `read_fname`; individual tools document any additional
raster inputs they require.

## Coordinate Reference Systems

Most geomorphic analyses measure lengths, areas, gradients, or search radii and
therefore require a projected CRS with horizontal units in metres. A geographic
latitude/longitude raster should normally be reprojected before analysis, and in fact our `lsdtt-raster-preprocessing` to enforces this. 
You should use `lsdtt-raster-preprocessing` before you do any serious analyis, it insures your raster is in the right format and it can additionaly can inspect, reproject, crop, and prepare a DEM.

We also have tools to make synthetic rasters, which receive a projected default georeference unless the generating
tool (you have three choices, see command line reference) provides another CRS option.

## Nodata and Masks

Nodata cells are excluded from analysis. Internally, lsdtt3 keeps raster values
and a validity mask together; output nodata cells carry the output raster's
declared nodata value. Be careful around raster edges and interior voids because
many flow-routing and neighbourhood operations depend on which cells are valid. Not literally around. This is not Tron and you don't like inside a computer. Just be careful if your raster has some holes, okay?

`lsdtt-raster-preprocessing` provides tools for filling interior nodata holes,
masking values, and preparing hydrologically suitable surfaces. We already told you to use `lsdtt-raster-preprocessing` before doing any serious analysis. Don't make us tell you that a third time. 

## Output Formats

Most raster-producing CLIs use `write_output_format`:

| Value | Result |
|---|---|
| `cog` | Cloud Optimized GeoTIFF; the standard default |
| `gtiff` | Conventional GeoTIFF |
| `same_as_input` | Use the input raster format where supported |
| `png` | PNG plus georeferencing where supported. This is nice if you want to look at the raster you made with no fancy software. |

Output filenames use `write_path` and `write_fname`, followed by a descriptor
specific to the result. The generated CLI reference gives the exact outputs and
format support for each executable.
