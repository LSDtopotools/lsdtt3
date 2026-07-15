import os

from lsdtt_qgis.naming import (
    expected_output_path,
    extension_for_format,
    split_source,
    strip_layer_suffix,
)


def test_strip_layer_suffix_removes_pipe_options():
    assert strip_layer_suffix("/data/basins.gpkg|layername=outlets") == "/data/basins.gpkg"
    assert strip_layer_suffix("/data/dem.tif") == "/data/dem.tif"


def test_split_source_returns_directory_and_filename():
    assert split_source("/data/dem.tif") == ("/data", "dem.tif")
    assert split_source("/data/basins.gpkg|layername=outlets") == ("/data", "basins.gpkg")


def test_split_source_rejects_bare_filename():
    assert split_source("dem.tif") is None
    assert split_source("") is None


def test_extension_for_format_maps_flatgeobuf_to_fgb():
    assert extension_for_format("flatgeobuf") == "fgb"
    assert extension_for_format("gpkg") == "gpkg"
    assert extension_for_format("csv") == "csv"


def test_expected_output_path_concatenates_prefix_and_suffix():
    assert expected_output_path("/out", "demo", "_fill.tif") == os.path.join(
        "/out", "demo_fill.tif"
    )
    assert expected_output_path("/out", "demo", "_depressions.fgb") == os.path.join(
        "/out", "demo_depressions.fgb"
    )
