# Synthetic Macrostrat tile fixture

One Mapbox Vector Tile (zoom 10, x 540, y 364; about 9.84-10.20 E,
45.83-46.10 N) in the layout `lsdtt-fetch-lithology` reads with
`macrostrat_tile_dir`: `<dir>/{z}/{x}/{y}.pbf`. It holds a layer `units` with
the fields of Macrostrat's `carto` tiles that the fetch reads. The data are
invented, not Macrostrat's.

The tile is generated from `units_source.geojson`:

```
ogr2ogr -f MVT out units_source.geojson -nln units \
    -dsco MINZOOM=10 -dsco MAXZOOM=10 -dsco COMPRESS=NO
cp out/10/540/364.pbf 10/540/364.pbf
```

| map_id | lith | Inside the test box (10.00-10.08 E, 46.00-46.05 N)? | Expected lith_primary (source) |
|---|---|---|---|
| 1 | `Major:{Granite,diorite}, Minor:{aplite}` | yes | `granite` (major) |
| 2 | `  Sandstone   and SHALE ` | yes | `sandstone and shale` (text) |
| 3 | `limestone` | no, dropped | - |
| 4 | empty | yes | empty (none) |

Used by `unit_fetch_lithology` and `smoke_fetch_lithology_offline`.
