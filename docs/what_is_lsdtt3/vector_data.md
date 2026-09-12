# Vector Data

lsdtt3 reads vector data through GDAL/OGR. GeoPackage and FlatGeobuf are the
preferred formats; Shapefile remains supported where its field-name and data-
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
