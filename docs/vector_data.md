# Vector data in lsdtt3

LSDTopoTools3 reads and writes vector data through GDAL/OGR. GeoPackage,
FlatGeobuf, GeoJSON, and other GDAL-supported vector formats are suitable for
most workflows.

## Coordinate reference systems

Keep input vectors in a defined CRS. Tools that combine vectors with a DEM
expect compatible projected coordinates where distances, areas, or raster cells
are involved. Reproject deliberately before analysis rather than relying on an
implicit geographic-coordinate calculation.

## Typical uses

- define a study-area polygon for DEM fetching or clipping;
- provide outlet/source points for flow and basin workflows;
- burn or repair a stream/vector network;
- export derived channel, basin, contour, or feature layers.

For Python workflows, GeoPandas inputs are accepted where a vector argument is
documented. Several result objects return GeoDataFrames, allowing direct
plotting or `GeoDataFrame.to_file` export.

See also [Fetch a DEM](python/fetching.md) and
[Fetching and vector API](python_api/fetch.md).
