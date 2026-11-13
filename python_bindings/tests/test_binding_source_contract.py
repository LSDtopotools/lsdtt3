from pathlib import Path

import pytest

pytestmark = [pytest.mark.conformance, pytest.mark.unit]


def test_road_flow_binding_uses_only_canonical_keywords():
    """The C++ binding of find_road_flow_discrimination declares
    threshold_contributing_pixels and raster_is_filled, not the old
    flow_accumulation_threshold_pixels and dem_is_filled.
    """
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


def test_mountain_front_outlet_keeps_approved_python_keyword():
    """The C++ binding of mountain_front_outlets declares min_contributing_pixels, not
    basin_minimum_contributing_pixels.
    """
    source = (
        Path(__file__).parents[1] / "bindings" / "raster_mountain_front.cpp"
    ).read_text(encoding="utf-8")
    binding = source.split(
        'm.def("mountain_front_outlets"', 1
    )[1]

    assert 'nb::arg("min_contributing_pixels")' in binding
    assert 'nb::arg("basin_minimum_contributing_pixels")' not in binding
