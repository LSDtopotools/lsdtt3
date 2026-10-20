"""Raster Preprocessing and Surface Metrics options added to close CLI gaps.

Each test runs processAlgorithm with QGIS stubbed out and checks the keys and
values written to the backend parameter file, plus the outputs collected.
"""

import os

import pytest

from test_phase4_algorithm_contract import (  # noqa: F401  (installs the qgis stubs)
    _Feedback,
    _Layer,
    _configure_parameters,
    _load,
    _touch_outputs,
)


def _run(module_name, class_name, tmp_path, monkeypatch, extra):
    module, algorithm_class = _load(module_name, class_name)
    monkeypatch.setattr(module, "attach_backend_help", lambda *args, **kwargs: None)
    captured = {}
    monkeypatch.setattr(
        module,
        "run_lsdtt_tool",
        lambda tool, values, folder, prefix, feedback: captured.update(values),
    )
    algorithm = algorithm_class()
    algorithm.initAlgorithm()
    _configure_parameters(algorithm)
    defaults = {
        p.name(): p.defaultValue()
        for p in algorithm.parameterDefinitions()
        # The stub reads a layer parameter's type list as its default; only
        # scalar defaults are real ones.
        if isinstance(p.defaultValue(), (bool, int, float, str))
    }
    params = dict(defaults)
    params.update({
        algorithm.INPUT_DEM: _Layer(str(tmp_path / "dem.tif")),
        algorithm.OUTPUT_FOLDER: str(tmp_path),
        algorithm.OUTPUT_PREFIX: "run",
    })
    params.update(extra(algorithm))
    return algorithm, params, captured


def test_preprocessing_defaults_match_the_cli(tmp_path, monkeypatch):
    algorithm, params, captured = _run(
        "raster_preprocessing", "RasterPreprocessingAlgorithm", tmp_path, monkeypatch,
        lambda a: {a.WRITE_FILL: False},
    )
    algorithm.processAlgorithm(params, None, _Feedback())
    # The CLI's own defaults (lsdtt-raster-preprocessing -h json).
    expected = {
        "smooth_algorithm_choice": "none",
        "perona_malik_iterations": 10,
        "perona_malik_lambda": 0.2,
        "perona_malik_conductance": "exponential",
        "perona_malik_k_quantile": 0.9,
        "perona_malik_k_metres_per_cell": 0.0,
        "write_stencil_difference_raster": False,
        "coherence_iterations": 10,
        "coherence_time_step": 0.15,
        "coherence_derivative_radius_metres": 0.0,
        "coherence_integration_radius_metres": 0.0,
        "coherence_contrast": 0.2,
        "coherence_min_diffusivity": 0.001,
        "vertical_units": "metres",
        "check_vertical_units_against_usgs": True,
        "priority_flood_queue": "binary",
        "fill_flat_resolution_choice": "inline_gradient",
        "fill_spill_merge_rainfall_rate": 0.1,
        "fsm_boundary_conditions": "oooo",
        "write_sink_candidates": False,
        "sink_candidate_min_pixels": 1,
    }
    for key, value in expected.items():
        if isinstance(value, float):
            assert captured[key] == pytest.approx(value), key
        else:
            assert captured[key] == value, key
    # Pixel replacement is only written when a point layer is given.
    assert "replace_pixels_fname" not in captured


def test_preprocessing_smoothing_replacement_and_sinks(tmp_path, monkeypatch):
    points = tmp_path / "pixels.fgb"
    points.write_text("", encoding="utf-8")
    algorithm, params, captured = _run(
        "raster_preprocessing", "RasterPreprocessingAlgorithm", tmp_path, monkeypatch,
        lambda a: {
            a.WRITE_FILL: False,
            a.SMOOTH_ALGORITHM: a.SMOOTH_ALGORITHM_OPTIONS.index("perona_malik"),
            a.WRITE_STENCIL_DIFF: True,
            a.FLAT_RESOLUTION: a.FLAT_RESOLUTION_OPTIONS.index("barnes2014b_gradient"),
            a.VERTICAL_UNITS: a.VERTICAL_UNITS_OPTIONS.index("us_survey_feet"),
            a.FSM_BOUNDARY_CONDITIONS: "NnoO",
            a.REPLACE_PIXELS: _Layer(str(points) + "|layername=pixels"),
            a.REPLACE_PIXELS_VALUE_FIELD: "z",
            a.WRITE_REPLACED_PIXELS: True,
            a.WRITE_SINK_CANDIDATES: True,
        },
    )
    _touch_outputs(str(tmp_path), "run", [
        "_smoothed.tif", "_stencil_difference.tif", "_pixels_replaced.tif",
        "_sink_candidates.fgb",
    ])
    outputs = algorithm.processAlgorithm(params, None, _Feedback())

    assert captured["smooth_algorithm_choice"] == "perona_malik"
    assert captured["write_stencil_difference_raster"] is True
    assert captured["fill_flat_resolution_choice"] == "barnes2014b_gradient"
    assert captured["vertical_units"] == "us_survey_feet"
    assert captured["fsm_boundary_conditions"] == "nnoo"
    assert captured["replace_pixels_fname"] == str(points)  # layer suffix stripped
    assert captured["replace_pixels_value_field"] == "z"
    assert captured["write_replaced_pixels_raster"] is True
    assert captured["write_sink_candidates"] is True
    for key, suffix in [
        (algorithm.OUTPUT_SMOOTHED, "_smoothed.tif"),
        (algorithm.OUTPUT_STENCIL_DIFF, "_stencil_difference.tif"),
        (algorithm.OUTPUT_REPLACED_PIXELS, "_pixels_replaced.tif"),
        (algorithm.OUTPUT_SINK_CANDIDATES, "_sink_candidates.fgb"),
    ]:
        assert os.path.basename(outputs[key]) == "run" + suffix


@pytest.mark.parametrize("extra, message", [
    (lambda a: {a.WRITE_STENCIL_DIFF: True}, "Perona-Malik"),
    (lambda a: {a.FSM_BOUNDARY_CONDITIONS: "oox"}, "four letters"),
    (lambda a: {a.WRITE_REPLACED_PIXELS: True}, "point layer"),
])
def test_preprocessing_rejects_inconsistent_options(tmp_path, monkeypatch, extra, message):
    algorithm, params, _ = _run(
        "raster_preprocessing", "RasterPreprocessingAlgorithm", tmp_path, monkeypatch, extra,
    )
    with pytest.raises(RuntimeError, match=message):
        algorithm.processAlgorithm(params, None, _Feedback())


def test_surface_metrics_image_outputs(tmp_path, monkeypatch):
    algorithm, params, captured = _run(
        "surface_metrics", "SurfaceMetricsAlgorithm", tmp_path, monkeypatch,
        lambda a: {
            a.WRITE_HILLSHADE: False,
            a.POLYFIT_METRICS: [],
            a.WRITE_FRANGI: True,
            a.FRANGI_TARGET: a.FRANGI_TARGET_OPTIONS.index("valleys"),
            a.WRITE_OTSU: True,
            a.WRITE_HOUGH_LINES: True,
            a.WRITE_CURVATURE_SKELETONS: True,
            a.WRITE_VEC_ADAPTIVE: True,
            a.WRITE_VEC_ADAPTIVE_SKELETONS: True,
            a.COMPUTE_BACKEND: a.COMPUTE_BACKEND_OPTIONS.index("cpu"),
        },
    )
    algorithm.parameterAsEnums = lambda p, key, context: p.get(key) or []
    suffixes = [suffix for _, _, suffix, _ in algorithm.IMAGE_OUTPUTS
                if suffix != "_frangi_scale.tif"]
    _touch_outputs(str(tmp_path), "run", suffixes)
    outputs = algorithm.processAlgorithm(params, None, _Feedback())

    assert captured["write_frangi_raster"] is True
    assert captured["write_frangi_scale_raster"] is False
    assert captured["frangi_target"] == "valleys"
    assert captured["compute_backend"] == "cpu"
    assert captured["vec_adaptive_block_size_pixels"] == 51
    # Defaults are the CLI's (lsdtt-surface-metrics -h json).
    assert captured["hough_threshold"] == 50
    assert captured["canny_threshold"] == pytest.approx(50.0)
    assert captured["curvature_skeleton_lambda_metres"] == pytest.approx(100.0)
    assert captured["frangi_n_scales"] == 4
    assert algorithm.OUTPUT_FRANGI_SCALE not in outputs
    for suffix in suffixes:
        assert any(os.path.basename(path) == "run" + suffix for path in outputs.values()), suffix


@pytest.mark.parametrize("extra, message", [
    (lambda a: {a.WRITE_VEC_ADAPTIVE_SKELETONS: True}, "adaptive VEC threshold masks"),
    (lambda a: {a.WRITE_VEC_ADAPTIVE: True, a.VEC_ADAPTIVE_BLOCK_SIZE: 50}, "odd"),
])
def test_surface_metrics_rejects_inconsistent_options(tmp_path, monkeypatch, extra, message):
    algorithm, params, _ = _run(
        "surface_metrics", "SurfaceMetricsAlgorithm", tmp_path, monkeypatch, extra,
    )
    algorithm.parameterAsEnums = lambda p, key, context: p.get(key) or []
    with pytest.raises(RuntimeError, match=message):
        algorithm.processAlgorithm(params, None, _Feedback())
