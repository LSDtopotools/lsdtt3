# Vector Data

`lsdtt3` has much more support of vector data format than `LSDTopoTools2`, and you can make nice river network data or polygons to your heart's content. 

`lsdtt3` reads vector data through [GDAL/OGR](https://gdal.org/en/stable/). You can read a Shapefile, which is a bizarre data format that people came up with because they were delirious from living in [Redlands](https://en.wikipedia.org/wiki/Esri). Modern formats, like [GeoPackage](https://en.wikipedia.org/wiki/GeoPackage) and [FlatGeobuf](https://flatgeobuf.org/) are the
preferred formats for `lsdtt3`; Shapefile remains supported where its field-name and data-
type limitations are acceptable.

## Multipart Geometry

CLI parameters that accept lines or points also accept their multipart
counterparts, such as `MultiLineString` and `MultiPoint`. Each feature is split
into individual parts on read, its attributes are copied to every part, and the
split is reported in the log. This avoids the silent loss of secondary parts
that occurred in versions before August 2026.

## Coordinate Systems

Vector inputs that are documented as accepting any CRS are reprojected to the
raster CRS internally. Other vector inputs must already use the raster's
projected coordinate system. Check the individual CLI parameter guidance for
the applicable contract.
