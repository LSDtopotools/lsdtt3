# Channel-head extraction

The simplest channel definition uses a drainage-area threshold. For more
specialised source detection, `extract_channel_heads` returns a `ChannelHeads`
object that can be converted back into the shared flow/channel workflow.

```python
heads = lsdtt3.extract_channel_heads(
    dem,
    flow,
    method="wiener",
    pruning_drainage_area_m2=2_000.0,
    connected_components_threshold_pixels=10,
)

scope = flow.scope_from_heads(heads)
network = flow.channel_network(scope)
source_points = heads.table
channel_points = network.points_geodataframe()
```

The `threshold` method is a transparent accumulation threshold. The `wiener`
method supplies a terrain-informed alternative with diagnostic products on the
result object. Compare several settings before treating detected heads as a
geomorphic result; `python_bindings/examples/channel_extraction_comparison.py`
does exactly this.

## Channel sources where a valley floor crosses the DEM edge

A channel that enters the DEM across its edge has no head inside the DEM, so
the threshold and the head detectors above cannot start it.
`find_perimeter_sources` traces the DEM perimeter, finds the gentle, wide, low
runs where a valley floor crosses it, and keeps each run whose path inward is
gentle and drains to one of the largest outlets:

```python
filled = lsdtt3.fill_pits(dem)
flow = lsdtt3.FlowInfo(filled, boundary_conditions="oooo")
perimeter = lsdtt3.find_perimeter_sources(dem, filled, flow, n_outlets=1)

scope = flow.channel_scope(threshold_contributing_pixels=500,
                           sources=perimeter.pixels)
```

`dem` is the DEM after nodata handling, and `filled` is the DEM the flow was
routed on. The flow needs open boundaries on every edge, `boundary_conditions="oooo"`.

The result holds one table per stage:

- `sources`: one row per source, with its cell, coordinates, elevation and
  the crossing it came from.
- `crossings`: every gentle run on the perimeter, with its `verdict` (for
  example `source`, `outlet_crossing` or `failed_inward`).
- `outlets`: the outer outlet runs, in rank order by summed drainage area,
  with a `selected` flag.
- `profile`: the perimeter trace, one row per visit.
- `report`: the resolved automatic parameters and the counts.

`min_valley_width_metres` (default 300 m) is the narrowest run that counts as a
valley-floor crossing, so lower it for narrower floors. The other keywords
match the `perimeter_*` parameters of `lsdtt-channel-extraction` without the
prefix: `max_gradient` is `perimeter_max_gradient`, and `n_outlets` is
`n_perimeter_outlets`. `help(lsdtt3.find_perimeter_sources)` gives each
default.
