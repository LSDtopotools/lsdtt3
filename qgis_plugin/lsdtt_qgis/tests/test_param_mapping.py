from lsdtt_qgis.templates.lsdtt_params import write_parameter_file


def test_parameter_file_is_sorted_and_boolean_lowercase(tmp_path):
    output = tmp_path / "params.driver"
    write_parameter_file(
        str(output),
        {
            "write_fill_raster": True,
            "read_fname": "dem.tif",
            "n_threads": 0,
        },
    )
    lines = output.read_text(encoding="utf-8").splitlines()
    assert lines == [
        "n_threads: 0",
        "read_fname: dem.tif",
        "write_fill_raster: true",
    ]
