# Vector Data

`lsdtt3` has much more support of vector data format than `LSDTopoTools2`, and you can make nice river network data or polygons to your heart's content. 

`lsdtt3` reads vector data through [GDAL/OGR](https://gdal.org/en/stable/). You can read a Shapefile, which is a bizarre data format that people came up with because they were delirious from living in [Redlands](https://en.wikipedia.org/wiki/Esri). Modern formats, like [GeoPackage](https://en.wikipedia.org/wiki/GeoPackage) and [FlatGeobuf](https://flatgeobuf.org/) are the
preferred formats for `lsdtt3`; Shapefile remains supported where its field-name and data-
type limitations are acceptable.

## Which formats, exactly

Read and write are not the same list, and it is worth knowing which is which
before you write a parameter file expecting an output you cannot have.

| Format | Extension | Read | Write |
|---|---|:---:|:---:|
| [FlatGeobuf](https://flatgeobuf.org/) | `.fgb` | yes | yes |
| [GeoPackage](https://en.wikipedia.org/wiki/GeoPackage) | `.gpkg` | yes | yes |
| GeoJSON | `.geojson` | yes | **no** |
| Shapefile | `.shp` | yes | **no** |
| CSV | `.csv` | yes | yes |

Two things that surprise people:

**GeoJSON and Shapefile are read-only.** lsdtt3 will happily take your basin
outlines or channel heads in either, but it writes FlatGeobuf or GeoPackage. If
you need a Shapefile at the end, `ogr2ogr` it — we are not going to add an
output format whose field names are capped at ten characters.

**CSV is deliberately plain.** It carries no CRS and gets no `.prj` sidecar. It
is the simple output: coordinates and attributes, nothing else. If you want
georeferencing that survives the trip into a GIS, write FlatGeobuf or
GeoPackage, both of which carry the CRS properly. The point-writing path emits
latitude and longitude in WGS84 so a bare CSV is at least unambiguous.

One practical note on GeoPackage: it is SQLite underneath, and SQLite does not
always behave on network or Windows-mounted filesystems. If you are working
under WSL2 on a `/mnt/c` path and a GeoPackage write fails, that is why —
lsdtt3 warns you when it spots this. Use FlatGeobuf instead.

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
