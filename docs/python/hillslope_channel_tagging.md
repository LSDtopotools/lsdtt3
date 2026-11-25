# Hillslope-channel tagging

These functions run the routines of
[`lsdtt-hillslope-channel-tagging`](../cli_reference/cli_reference_lsdtt_hillslope_channel_tagging.md)
on rasters in memory. The method, and what each ridge option means, are in
[Hillslope-channel tagging](../design_and_algorithms/hillslope_channel_tagging.md).

The options have the CLI's names without the `hillslope_` and
`channel_attribute_` prefixes. Python and the CLI call the same C++ routines
and the same file reader. The CLI holds floats in single precision, so to
reproduce a CLI run bit for bit pass `m_over_n=float(numpy.float32(0.45))`.

!!! note "Experimental"
    These bindings are new. `helper_for_llm` marks them `[experimental]` until
    they are promoted.

| Python | Returns | CLI equivalent |
|---|---|---|
| `FlowInfo.hillslope_relief(elevation, channels, ...)` | dict of link arrays, Rasters and counts | `hillslope_ridge_method` and the relief outputs |
| `FlowInfo.hillslope_chi_tags(routing_surface, channel_flat, ...)` | dict with `chi`, `k_sn` and `basins` | `tag_chi`, `tag_k_sn` |
| `FlowInfo.channel_attributes_to_raster(source, channels, columns, ...)` | dict with `values` on channel pixels and counts | `channel_attributes_fname` |
| `lsdtt3.tag_hillslopes_from_channel(channel_flat, channel_values, dem)` | dict with `tag` and counts | `_tag_<column>` |
| `lsdtt3.summarise_hillslopes_by_channel_pixel(channel_flat, values)` | DataFrame | `write_channel_pixel_summary` |

## Relief and relative height

The CLI carves, then fills, then routes, and measures heights on that routing
surface. With `aspect_ray`, pass the input DEM as `aspect_surface`, which is the
CLI's default.

```python
import numpy as np
import lsdtt3

dem = lsdtt3.Raster.read("my_dem.tif")
routing = lsdtt3.fill_pits(lsdtt3.breach_least_cost(dem))
flow = lsdtt3.FlowInfo(routing, boundary_conditions="oooo")

area = flow.contributing_pixels()
is_channel = area.mask_numpy().astype(bool) & (area.to_numpy() >= 500)
channels = lsdtt3.Raster.from_numpy(
    is_channel.astype(np.float32), routing.mask_numpy().astype(bool), routing.metadata
)

relief = flow.hillslope_relief(
    routing, channels, ridge_method="aspect_ray", aspect_surface=dem
)
relative_height = relief["relative_height"]
```

`relief["channel_flat"]` is each pixel's D8 channel pixel, and every tag uses it.

## Chi and k_sn tags

Give exactly one of `outlet_flats` and `chi_max_outlet_elevation_metres`.

```python
tags = flow.hillslope_chi_tags(
    routing,
    relief["channel_flat"],
    chi_max_outlet_elevation_metres=5.0,
    threshold_contributing_pixels=500,
    m_over_n=0.45,
)
k_sn = tags["k_sn"]["tag"]
print(tags["k_sn"]["n_hillslope_tagged"], tags["k_sn"]["n_excluded"])
```

Check `n_excluded`: hillslopes of basins left out are nodata.

## Imported channel attributes

`source` is a CSV or vector file, a DataFrame or a GeoDataFrame. A table is
written to a temporary file and read by the CLI's reader.

```python
joined = flow.channel_attributes_to_raster(
    "my_channel_planform.fgb",
    channels,
    ["sinuosity"],
    max_snap_distance_metres=60.0,
    max_gap_fill_distance_metres=150.0,
)
print(joined["counts"]["sinuosity"])

sinuosity = lsdtt3.tag_hillslopes_from_channel(
    relief["channel_flat"], joined["values"]["sinuosity"], routing
)["tag"]
```
