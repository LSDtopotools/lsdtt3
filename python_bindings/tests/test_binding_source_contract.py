from pathlib import Path


def test_road_flow_binding_uses_only_canonical_keywords():
    source = (
        Path(__file__).parents[1] / "bindings" / "raster_feature_finder.cpp"
    ).read_text(encoding="utf-8")
    binding = source.split(
        'm.def("find_road_flow_discrimination"', 1
    )[1].split('m.def("extract_roads_smooth"', 1)[0]

    assert 'nb::arg("threshold_contributing_pixels")' in binding
    assert 'nb::arg("raster_is_filled")' in binding
    assert 'nb::arg("flow_accumulation_threshold_pixels")' not in binding
    assert 'nb::arg("dem_is_filled")' not in binding
