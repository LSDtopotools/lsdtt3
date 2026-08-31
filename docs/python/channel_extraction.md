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
