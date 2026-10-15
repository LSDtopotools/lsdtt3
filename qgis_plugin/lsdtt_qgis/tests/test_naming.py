import os

from lsdtt_qgis.naming import (
    build_backend_help_index,
    expected_output_path,
    extension_for_format,
    lookup_param_help,
    next_qgis_driver_name,
    split_source,
    strip_layer_suffix,
    source_layer_name,
    write_line_geojson,
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


def test_next_qgis_driver_name_starts_at_001(tmp_path):
    assert next_qgis_driver_name(str(tmp_path), "my_dem") == (
        "my_dem_qgis_analysis_001.driver"
    )


def test_next_qgis_driver_name_increments_past_existing(tmp_path):
    (tmp_path / "my_dem_qgis_analysis_001.driver").write_text("x")
    (tmp_path / "my_dem_qgis_analysis_007.driver").write_text("x")
    assert next_qgis_driver_name(str(tmp_path), "my_dem") == (
        "my_dem_qgis_analysis_008.driver"
    )


def test_next_qgis_driver_name_ignores_other_stems_and_junk(tmp_path):
    (tmp_path / "other_dem_qgis_analysis_004.driver").write_text("x")
    (tmp_path / "my_dem_qgis_analysis_abc.driver").write_text("x")
    (tmp_path / "my_dem_qgis_analysis_002.txt").write_text("x")
    assert next_qgis_driver_name(str(tmp_path), "my_dem") == (
        "my_dem_qgis_analysis_001.driver"
    )


def test_next_qgis_driver_name_survives_missing_directory():
    assert next_qgis_driver_name("/nonexistent/dir", "dem").endswith("_001.driver")


def _sample_help_payload():
    return {
        "tool": "lsdtt-raster-preprocessing",
        "version": "0.2.0",
        "parameters": [
            {
                "name": "minimum_elevation",
                "description": "Mask pixels at or below this elevation [m].",
                "guidance": "Default 0.0 catches ocean pixels.",
            },
            {
                "name": "write_fill_raster",
                "description": "Write the pit-filled DEM.",
                "guidance": "",
            },
            {
                "name": "fill_algorithm_choice",
                "description": "Pit-filling algorithm.",
                "guidance": "barnes2014 or wang_liu.",
            },
            {"name": "n_threads", "description": "Thread count.", "guidance": ""},
            {"name": "no_text", "description": "", "guidance": ""},
        ],
    }


def test_build_backend_help_index_handles_bad_payloads():
    assert build_backend_help_index(None) == {}
    assert build_backend_help_index({"parameters": "nope"}) == {}
    index = build_backend_help_index(_sample_help_payload())
    assert "minimum_elevation" in index


def test_lookup_param_help_exact_and_suffix_matches():
    index = build_backend_help_index(_sample_help_payload())
    text = lookup_param_help(index, "MINIMUM_ELEVATION")
    assert text.startswith("Mask pixels")
    assert "ocean pixels" in text
    assert lookup_param_help(index, "WRITE_FILL") == "Write the pit-filled DEM."
    assert lookup_param_help(index, "FILL_ALGORITHM").startswith("Pit-filling")


def test_lookup_param_help_alias_override_and_misses():
    index = build_backend_help_index(_sample_help_payload())
    assert lookup_param_help(index, "THREADS") == "Thread count."
    assert lookup_param_help(index, "SOMETHING_ELSE") is None
    assert lookup_param_help(index, "SOMETHING_ELSE", override="n_threads") == (
        "Thread count."
    )
    assert lookup_param_help(index, "NO_TEXT") is None


def test_write_line_geojson_writes_one_linestring(tmp_path):
    import json

    path = write_line_geojson(
        str(tmp_path / "baseline.geojson"), [(500000, 6000250.5), (500500, 6000250.5)],
        "EPSG:32633",
    )
    data = json.loads(open(path, encoding="utf-8").read())
    assert data["crs"]["properties"]["name"] == "EPSG:32633"
    (feature,) = data["features"]
    assert feature["geometry"] == {
        "type": "LineString",
        "coordinates": [[500000.0, 6000250.5], [500500.0, 6000250.5]],
    }


def test_write_line_geojson_rejects_a_single_point(tmp_path):
    import pytest

    with pytest.raises(ValueError, match="at least two points"):
        write_line_geojson(str(tmp_path / "b.geojson"), [(0, 0)])


def test_source_layer_name_prefers_layername_option():
    assert source_layer_name("/d/soils.gpkg|layername=hutton") == "hutton"
    assert source_layer_name("/d/soils.gpkg|subset=x|layername=lu") == "lu"


def test_source_layer_name_falls_back_to_file_stem():
    assert source_layer_name("/d/hutton_series_code_min.fgb") == "hutton_series_code_min"
    assert source_layer_name("/d/landuse.shp|layerid=0") == "landuse"
