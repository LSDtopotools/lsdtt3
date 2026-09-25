import importlib
import os
import sys
import types

import pytest


class _Parameter:
    FlagAdvanced = 1
    Integer = 0
    Double = 1

    def __init__(self, name=None, *args, **kwargs):
        self._name = name
        self._flags = 0
        self._description = args[0] if args else None
        self._default = kwargs.get("defaultValue", args[1] if len(args) > 1 else None)

    def name(self):
        return self._name

    def flags(self):
        return self._flags

    def setFlags(self, flags):
        self._flags = flags

    def description(self):
        return self._description

    def defaultValue(self):
        return self._default


class _Algorithm:
    def __init__(self, *args, **kwargs):
        self._parameters = []
        self._outputs = []

    def addParameter(self, parameter):
        self._parameters.append(parameter)

    def addOutput(self, output):
        self._outputs.append(output)

    def parameterDefinitions(self):
        return self._parameters

    def outputDefinitions(self):
        return self._outputs


def _install_qgis_stubs():
    qgis = types.ModuleType("qgis")
    core = types.ModuleType("qgis.core")
    core.QgsProcessingAlgorithm = _Algorithm
    core.QgsProcessingException = RuntimeError
    core.QgsProcessingParameterDefinition = _Parameter
    core.QgsProcessingParameterNumber = _Parameter
    core.QgsProcessing = types.SimpleNamespace(
        TypeRaster=1,
        TypeVectorLine=2,
        TypeVectorPoint=3,
        TypeVectorPolygon=4,
        TypeVectorAnyGeometry=5,
    )
    core.__getattr__ = lambda name: _Parameter
    qgis.core = core
    sys.modules.setdefault("qgis", qgis)
    sys.modules.setdefault("qgis.core", core)


_install_qgis_stubs()

ALGORITHMS = [
    ("basin_extraction", "BasinExtractionAlgorithm", "basin_extraction"),
    ("channel_extraction", "ChannelExtractionAlgorithm", "channel_extraction"),
    (
        "channel_network_analysis",
        "ChannelNetworkAnalysisAlgorithm",
        "channel_network_analysis",
    ),
    ("chi_analysis", "ChiAnalysisAlgorithm", "chi_analysis"),
    (
        "side_tributary_basins",
        "SideTributaryBasinsAlgorithm",
        "side_tributary_basins",
    ),
]


def _load(module_name, class_name):
    module = importlib.import_module("lsdtt_qgis.algorithms." + module_name)
    return module, getattr(module, class_name)


class _Layer:
    def __init__(self, source):
        self._source = source

    def source(self):
        return self._source


class _Feedback:
    def __init__(self):
        self.messages = []

    def pushInfo(self, message):
        self.messages.append(message)


def _configure_parameters(algorithm):
    algorithm.parameterAsRasterLayer = lambda params, key, context: params.get(key)
    algorithm.parameterAsVectorLayer = lambda params, key, context: params.get(key)
    algorithm.parameterAsFile = lambda params, key, context: params.get(key, "")
    algorithm.parameterAsString = lambda params, key, context: params.get(key, "")
    algorithm.parameterAsBool = lambda params, key, context: bool(params.get(key, False))
    algorithm.parameterAsInt = lambda params, key, context: int(params.get(key, 0))
    algorithm.parameterAsDouble = lambda params, key, context: float(params.get(key, 0.0))
    algorithm.parameterAsEnum = lambda params, key, context: int(params.get(key, 0))


def _base_parameters(algorithm, tmp_path):
    return {
        algorithm.INPUT_DEM: _Layer(str(tmp_path / "dem.tif")),
        algorithm.OUTPUT_FOLDER: str(tmp_path),
        algorithm.OUTPUT_PREFIX: "run",
        algorithm.BOUNDARY_CONDITIONS: "oooo",
    }


def _touch_outputs(folder, prefix, suffixes):
    for suffix in suffixes:
        path = os.path.join(folder, prefix + suffix)
        with open(path, "w", encoding="utf-8"):
            pass


def test_phase4_public_ids_and_parameter_id_are_pinned():
    for module_name, class_name, public_id in ALGORITHMS:
        _, algorithm_class = _load(module_name, class_name)
        algorithm = algorithm_class()
        assert algorithm.name() == public_id
        assert algorithm.THRESHOLD_CONTRIBUTING_PIXELS == (
            "THRESHOLD_CONTRIBUTING_PIXELS"
        )
        assert not hasattr(algorithm, "SOURCE_THRESHOLD")


def test_phase4_backend_help_maps_threshold_to_canonical_key(monkeypatch):
    for module_name, class_name, _ in ALGORITHMS:
        module, algorithm_class = _load(module_name, class_name)
        captured = {}

        def capture(algorithm, tool_name, overrides=None):
            captured.update(overrides or {})

        monkeypatch.setattr(module, "attach_backend_help", capture)
        algorithm = algorithm_class()
        algorithm.initAlgorithm()
        assert captured[algorithm.THRESHOLD_CONTRIBUTING_PIXELS] == (
            "threshold_contributing_pixels"
        )
        assert algorithm.THRESHOLD_CONTRIBUTING_PIXELS in {
            parameter.name() for parameter in algorithm.parameterDefinitions()
        }


def test_phase4_emitted_key_construction_is_canonical():
    for module_name, class_name, _ in ALGORITHMS:
        _, algorithm_class = _load(module_name, class_name)
        constants = set()
        pending = list(algorithm_class.processAlgorithm.__code__.co_consts)
        while pending:
            value = pending.pop()
            if isinstance(value, tuple):
                pending.extend(value)
            elif isinstance(value, str):
                constants.add(value)
        assert "threshold_contributing_pixels" in constants
        assert "source_pixels_threshold" not in constants


def test_phase5_window_and_outlet_snap_mappings_are_canonical(monkeypatch):
    outlet_algorithms = [
        ("basin_extraction", "BasinExtractionAlgorithm"),
        ("basin_stats", "BasinStatsAlgorithm"),
        ("channel_network_analysis", "ChannelNetworkAnalysisAlgorithm"),
        ("chi_analysis", "ChiAnalysisAlgorithm"),
    ]
    for module_name, class_name in outlet_algorithms:
        module, algorithm_class = _load(module_name, class_name)
        captured = {}
        monkeypatch.setattr(
            module,
            "attach_backend_help",
            lambda algorithm, tool_name, overrides=None: captured.update(overrides or {}),
        )
        algorithm = algorithm_class()
        algorithm.initAlgorithm()
        assert captured[algorithm.SNAPPING_MIN_RADIUS] == (
            "basin_outlet_snap_min_search_radius_pixels"
        )
        assert captured[algorithm.SNAPPING_MAX_RADIUS] == (
            "basin_outlet_snap_max_search_radius_pixels"
        )

    module, algorithm_class = _load(
        "channel_network_analysis", "ChannelNetworkAnalysisAlgorithm"
    )
    captured = {}
    monkeypatch.setattr(
        module,
        "attach_backend_help",
        lambda algorithm, tool_name, overrides=None: captured.update(overrides or {}),
    )
    algorithm = algorithm_class()
    algorithm.initAlgorithm()
    assert captured[algorithm.VECTOR_FORMAT] == "write_vector_format"

    module, algorithm_class = _load("channel_extraction", "ChannelExtractionAlgorithm")
    captured = {}
    monkeypatch.setattr(
        module,
        "attach_backend_help",
        lambda algorithm, tool_name, overrides=None: captured.update(overrides or {}),
    )
    algorithm = algorithm_class()
    algorithm.initAlgorithm()
    assert captured[algorithm.WIENER_FITTING_RADIUS] == (
        "wiener_surface_fitting_radius_metres"
    )


def test_phase4_backend_version_gate_runs_before_execution(monkeypatch):
    support = importlib.import_module("lsdtt_qgis.algorithms.support")

    class Backend:
        def backend_path(self):
            return "/old-backend"

    monkeypatch.setattr(support, "BackendManager", Backend)
    monkeypatch.setattr(
        support, "load_manifest", lambda path: {"backend_version": "0.4.0"}
    )
    with pytest.raises(RuntimeError, match="requires version 0.5.0 or newer"):
        support.run_lsdtt_tool("lsdtt-basin-extraction", {}, "/out", "run", None)


def test_backend_plugin_version_gate_runs_before_execution(monkeypatch):
    support = importlib.import_module("lsdtt_qgis.algorithms.support")

    class Backend:
        def backend_path(self):
            return "/incompatible-backend"

    monkeypatch.setattr(support, "BackendManager", Backend)
    monkeypatch.setattr(
        support,
        "load_manifest",
        lambda path: {
            "backend_version": "0.5.0",
            "plugin_version_range": ">=0.1.0 <0.2.0",
        },
    )
    with pytest.raises(RuntimeError, match="current plugin version is 0.3.0"):
        support.run_lsdtt_tool("lsdtt-basin-extraction", {}, "/out", "run", None)


def test_surface_metrics_uses_canonical_slope_threshold_key(monkeypatch):
    module, algorithm_class = _load("surface_metrics", "SurfaceMetricsAlgorithm")
    captured = {}
    monkeypatch.setattr(
        module,
        "attach_backend_help",
        lambda algorithm, tool_name, overrides=None: captured.update(overrides or {}),
    )
    algorithm = algorithm_class()
    algorithm.initAlgorithm()
    assert captured[algorithm.SLOPE_THRESHOLD_LAYERS] == "slope_thresholds"

    constants = set(algorithm.processAlgorithm.__code__.co_consts)
    assert "slope_thresholds" in constants
    assert "slope_threshold_layers" not in constants


def test_basin_stats_threshold_default_and_label(monkeypatch):
    module, algorithm_class = _load("basin_stats", "BasinStatsAlgorithm")
    monkeypatch.setattr(module, "attach_backend_help", lambda *args, **kwargs: None)
    algorithm = algorithm_class()
    algorithm.initAlgorithm()
    threshold = next(
        parameter
        for parameter in algorithm.parameterDefinitions()
        if parameter.name() == algorithm.THRESHOLD_CONTRIBUTING_PIXELS
    )
    assert threshold.defaultValue() == 1000
    assert threshold.description() == "Threshold contributing pixels"


def test_basin_stats_area_ratios_do_not_require_channel_stats(monkeypatch, tmp_path):
    module, algorithm_class = _load("basin_stats", "BasinStatsAlgorithm")
    captured = {}
    monkeypatch.setattr(module, "export_points_as_latlon_csv", lambda *args: None)
    monkeypatch.setattr(module, "require_output", lambda *args: str(tmp_path / "output"))
    monkeypatch.setattr(
        module,
        "run_lsdtt_tool",
        lambda tool, values, folder, prefix, feedback: captured.update(values),
    )

    class Layer:
        def __init__(self, source):
            self._source = source

        def source(self):
            return self._source

    algorithm = algorithm_class()
    values = {
        algorithm.INPUT_DEM: Layer(str(tmp_path / "dem.tif")),
        algorithm.INPUT_OUTLETS: Layer(str(tmp_path / "outlets.gpkg")),
        algorithm.OUTPUT_FOLDER: str(tmp_path),
        algorithm.OUTPUT_PREFIX: "basins",
        algorithm.CALCULATE_CHANNEL_STATS: False,
        algorithm.CALCULATE_JUNCTION_AREA_RATIOS: True,
        algorithm.WRITE_BASIN_POLYGONS: False,
    }
    algorithm.parameterAsRasterLayer = lambda params, key, context: params.get(key)
    algorithm.parameterAsVectorLayer = lambda params, key, context: params.get(key)
    algorithm.parameterAsString = lambda params, key, context: params.get(key, "")
    algorithm.parameterAsBool = lambda params, key, context: bool(params.get(key, False))
    algorithm.parameterAsInt = lambda params, key, context: int(params.get(key, 0))
    algorithm.parameterAsDouble = lambda params, key, context: float(params.get(key, 0.0))
    algorithm.parameterAsEnum = lambda params, key, context: int(params.get(key, 0))

    outputs = algorithm.processAlgorithm(values, None, None)
    assert captured["calculate_channel_stats"] is False
    assert captured["calculate_junction_area_ratios"] is True
    assert algorithm.OUTPUT_AREA_RATIOS in outputs


@pytest.mark.parametrize(
    "module_name,class_name,parameter_name,backend_key",
    [
        (
            "clip_raster_to_polygon",
            "ClipRasterToPolygonAlgorithm",
            "INPUT_POLYGON",
            "polygon_clip_fname",
        ),
        (
            "clip_vector_to_raster",
            "ClipVectorToRasterAlgorithm",
            "INPUT_VECTOR",
            "vector_data_fname",
        ),
        (
            "mosaic_rasters",
            "MosaicRastersAlgorithm",
            "INPUT_RASTERS",
            "mosaic_file_list",
        ),
        (
            "sample_rasters_to_points",
            "SampleRastersToPointsAlgorithm",
            "RASTER_ATTRIBUTE_CSV",
            "raster_attribute_fname",
        ),
        (
            "sample_rasters_to_points",
            "SampleRastersToPointsAlgorithm",
            "INPUT_POINTS",
            "point_data_fname",
        ),
        (
            "sample_rasters_to_points",
            "SampleRastersToPointsAlgorithm",
            "LAYER_NAME",
            "point_data_layer_name",
        ),
    ],
)
def test_utility_backend_help_overrides(
    monkeypatch, module_name, class_name, parameter_name, backend_key
):
    module, algorithm_class = _load(module_name, class_name)
    captured = {}
    monkeypatch.setattr(
        module,
        "attach_backend_help",
        lambda algorithm, tool_name, overrides=None: captured.update(overrides or {}),
    )
    algorithm = algorithm_class()
    algorithm.initAlgorithm()
    assert captured[getattr(algorithm, parameter_name)] == backend_key


def test_channel_network_all_exposed_outputs_are_declared_and_returned(monkeypatch, tmp_path):
    module, algorithm_class = _load(
        "channel_network_analysis", "ChannelNetworkAnalysisAlgorithm"
    )
    monkeypatch.setattr(module, "attach_backend_help", lambda *args, **kwargs: None)
    algorithm = algorithm_class()
    algorithm.initAlgorithm()
    declared = {output.name() for output in algorithm.outputDefinitions()}
    expected_keys = {
        algorithm.OUTPUT_LINK_BEARINGS,
        algorithm.OUTPUT_JUNCTION_ANGLES,
        algorithm.OUTPUT_STRAHLER_CHANNELS,
        algorithm.OUTPUT_TOKUNAGA_CHANNELS,
        algorithm.OUTPUT_TOKUNAGA_CHANNEL_DATA,
        algorithm.OUTPUT_TOKUNAGA_SUMMARY,
        algorithm.OUTPUT_TOKUNAGA_PAIR_COUNTS,
        algorithm.OUTPUT_TOKUNAGA_OFFSET_STATS,
        algorithm.OUTPUT_AREA_RATIOS,
        algorithm.OUTPUT_AREA_RATIO_SUMMARY,
        algorithm.OUTPUT_AREA_RATIO_DECADES,
        algorithm.OUTPUT_COM_BEARINGS,
        algorithm.OUTPUT_LINK_BEARING_VECTORS,
        algorithm.OUTPUT_JUNCTION_ANGLE_VECTORS,
    }
    assert expected_keys <= declared

    suffixes = [
        "_channel_link_bearings.csv",
        "_junction_angles.csv",
        "_strahler_channels.csv",
        "_tokunaga_channels.csv",
        "_tokunaga_channel_data.csv",
        "_tokunaga_summary_pair_counts.csv",
        "_tokunaga_summary_offset_stats.csv",
        "_tokunaga_summary_fit.csv",
        "_junction_area_ratios.csv",
        "_junction_area_ratio_summary.csv",
        "_junction_area_ratio_decades.csv",
        "_junction_com_bearings.csv",
        "_channel_link_bearing_vectors.fgb",
        "_junction_angle_bearing_vectors.fgb",
    ]
    captured = {}

    def fake_run(tool, values, folder, prefix, feedback):
        captured.update(values)
        _touch_outputs(folder, prefix, suffixes)

    monkeypatch.setattr(module, "run_lsdtt_tool", fake_run)
    _configure_parameters(algorithm)
    parameters = _base_parameters(algorithm, tmp_path)
    parameters.update({
        algorithm.VECTOR_FORMAT: 2,
        algorithm.WRITE_LINK_BEARINGS: True,
        algorithm.WRITE_JUNCTION_ANGLES: True,
        algorithm.WRITE_STRAHLER_CHANNELS: True,
        algorithm.WRITE_TOKUNAGA_CHANNELS: True,
        algorithm.WRITE_TOKUNAGA_CHANNEL_DATA: True,
        algorithm.WRITE_TOKUNAGA_SUMMARY: True,
        algorithm.WRITE_AREA_RATIO_DATA: True,
        algorithm.WRITE_AREA_RATIO_SUMMARY: True,
        algorithm.WRITE_COM_BEARINGS: True,
        algorithm.WRITE_LINK_BEARING_VECTORS: True,
        algorithm.WRITE_JUNCTION_ANGLE_VECTORS: True,
    })
    outputs = algorithm.processAlgorithm(parameters, None, _Feedback())
    assert set(outputs) == expected_keys
    assert outputs[algorithm.OUTPUT_TOKUNAGA_SUMMARY].endswith(
        "run_tokunaga_summary_fit.csv"
    )
    assert outputs[algorithm.OUTPUT_AREA_RATIO_DECADES].endswith(
        "run_junction_area_ratio_decades.csv"
    )
    assert outputs[algorithm.OUTPUT_LINK_BEARING_VECTORS].endswith(".fgb")
    assert captured["write_channel_link_bearing_vectors"] is True


def test_channel_network_qa_only_request_passes_output_gate(monkeypatch, tmp_path):
    module, algorithm_class = _load(
        "channel_network_analysis", "ChannelNetworkAnalysisAlgorithm"
    )
    algorithm = algorithm_class()
    _configure_parameters(algorithm)
    parameters = _base_parameters(algorithm, tmp_path)
    parameters[algorithm.WRITE_LINK_BEARINGS] = False
    parameters[algorithm.WRITE_LINK_BEARING_VECTORS] = True
    monkeypatch.setattr(
        module,
        "run_lsdtt_tool",
        lambda tool, values, folder, prefix, feedback: _touch_outputs(
            folder, prefix, ["_channel_link_bearing_vectors.fgb"]
        ),
    )
    outputs = algorithm.processAlgorithm(parameters, None, _Feedback())
    assert set(outputs) == {algorithm.OUTPUT_LINK_BEARING_VECTORS}


def test_basin_grouped_flatgeobuf_outputs_are_numerically_ordered(monkeypatch, tmp_path):
    module, algorithm_class = _load("basin_extraction", "BasinExtractionAlgorithm")
    algorithm = algorithm_class()
    _configure_parameters(algorithm)
    parameters = _base_parameters(algorithm, tmp_path)
    parameters.update({
        algorithm.WRITE_BASIN_OUTLINES: True,
        algorithm.OUTLINES_AS_POLYGON: True,
        algorithm.OUTLINES_PER_ORDER: True,
    })

    def fake_run(tool, values, folder, prefix, feedback):
        _touch_outputs(folder, prefix, [
            "_junction_basin_polygons_order_10.fgb",
            "_junction_basin_polygons_order_2.fgb",
        ])

    monkeypatch.setattr(module, "run_lsdtt_tool", fake_run)
    feedback = _Feedback()
    outputs = algorithm.processAlgorithm(parameters, None, feedback)
    assert [os.path.basename(path) for path in outputs[algorithm.OUTPUT_BASIN_OUTLINES]] == [
        "run_junction_basin_polygons_order_2.fgb",
        "run_junction_basin_polygons_order_10.fgb",
    ]
    assert "2 output file(s)" in feedback.messages[-1]


def test_basin_grouped_gpkg_is_one_multilayer_result(monkeypatch, tmp_path):
    module, algorithm_class = _load("basin_extraction", "BasinExtractionAlgorithm")
    algorithm = algorithm_class()
    _configure_parameters(algorithm)
    parameters = _base_parameters(algorithm, tmp_path)
    parameters.update({
        algorithm.WRITE_BASIN_OUTLINES: True,
        algorithm.OUTLINES_AS_POLYGON: False,
        algorithm.OUTLINES_PER_ORDER: True,
        algorithm.BASIN_OUTLINES_FORMAT: 1,
    })
    monkeypatch.setattr(
        module,
        "run_lsdtt_tool",
        lambda tool, values, folder, prefix, feedback: _touch_outputs(
            folder, prefix, ["_junction_basin_outlines.gpkg"]
        ),
    )
    outputs = algorithm.processAlgorithm(parameters, None, _Feedback())
    assert outputs[algorithm.OUTPUT_BASIN_OUTLINES] == [
        str(tmp_path / "run_junction_basin_outlines.gpkg")
    ]


def test_disabled_basin_outlines_force_dependent_flags_false(monkeypatch, tmp_path):
    module, algorithm_class = _load("basin_extraction", "BasinExtractionAlgorithm")
    algorithm = algorithm_class()
    _configure_parameters(algorithm)
    parameters = _base_parameters(algorithm, tmp_path)
    parameters.update({
        algorithm.WRITE_BASIN_OUTLINES: False,
        algorithm.OUTLINES_AS_POLYGON: True,
        algorithm.OUTLINES_PER_ORDER: True,
        algorithm.WRITE_CHANNEL_LINES: True,
    })
    captured = {}

    def fake_run(tool, values, folder, prefix, feedback):
        captured.update(values)
        _touch_outputs(folder, prefix, ["_channel_network_lines.fgb"])

    monkeypatch.setattr(module, "run_lsdtt_tool", fake_run)
    algorithm.processAlgorithm(parameters, None, _Feedback())
    assert captured["write_basin_outlines"] is False
    assert captured["write_basin_outlines_as_polygon"] is False
    assert captured["write_basin_outlines_per_order"] is False


def test_channel_extraction_burned_profiles_only_passes_gate(monkeypatch, tmp_path):
    module, algorithm_class = _load("channel_extraction", "ChannelExtractionAlgorithm")
    algorithm = algorithm_class()
    _configure_parameters(algorithm)
    parameters = _base_parameters(algorithm, tmp_path)
    parameters.update({
        algorithm.WRITE_CHANNEL_LINES: False,
        algorithm.WRITE_BURNED_PROFILES: True,
        algorithm.STREAM_VECTOR: _Layer(str(tmp_path / "streams.gpkg")),
    })
    monkeypatch.setattr(
        module,
        "run_lsdtt_tool",
        lambda tool, values, folder, prefix, feedback: _touch_outputs(
            folder, prefix, ["_burned_stream_profiles.csv"]
        ),
    )
    outputs = algorithm.processAlgorithm(parameters, None, _Feedback())
    assert set(outputs) == {algorithm.OUTPUT_BURNED_PROFILES}


def test_channel_extraction_wiener_diagnostics_only_pass_gate(monkeypatch, tmp_path):
    module, algorithm_class = _load("channel_extraction", "ChannelExtractionAlgorithm")
    capabilities = importlib.import_module("lsdtt_qgis.backend.capabilities")
    monkeypatch.setattr(capabilities, "require_feature", lambda *args, **kwargs: None)
    algorithm = algorithm_class()
    _configure_parameters(algorithm)
    parameters = _base_parameters(algorithm, tmp_path)
    parameters.update({
        algorithm.SOURCE_ALGORITHM: 1,
        algorithm.WRITE_CHANNEL_LINES: False,
        algorithm.WIENER_WRITE_FILTERED_DEM: True,
        algorithm.WIENER_WRITE_CURVATURE: True,
        algorithm.WIENER_WRITE_QQ: True,
    })
    monkeypatch.setattr(
        module,
        "run_lsdtt_tool",
        lambda tool, values, folder, prefix, feedback: _touch_outputs(folder, prefix, [
            "_wiener_filtered_dem.tif",
            "_wiener_tangential_curvature.tif",
            "_wiener_qq.csv",
        ]),
    )
    outputs = algorithm.processAlgorithm(parameters, None, _Feedback())
    assert set(outputs) == {
        algorithm.OUTPUT_WIENER_FILTERED_DEM,
        algorithm.OUTPUT_WIENER_CURVATURE,
        algorithm.OUTPUT_WIENER_QQ,
    }


def test_channel_extraction_rejects_wiener_diagnostic_in_threshold_mode(tmp_path):
    _, algorithm_class = _load("channel_extraction", "ChannelExtractionAlgorithm")
    algorithm = algorithm_class()
    _configure_parameters(algorithm)
    parameters = _base_parameters(algorithm, tmp_path)
    parameters.update({
        algorithm.WRITE_CHANNEL_LINES: False,
        algorithm.WIENER_WRITE_QQ: True,
    })
    with pytest.raises(RuntimeError, match="Wiener diagnostics require"):
        algorithm.processAlgorithm(parameters, None, _Feedback())


def test_chi_segment_lines_auto_enable_fitting_and_are_required(monkeypatch, tmp_path):
    module, algorithm_class = _load("chi_analysis", "ChiAnalysisAlgorithm")
    algorithm = algorithm_class()
    _configure_parameters(algorithm)
    parameters = _base_parameters(algorithm, tmp_path)
    parameters.update({
        algorithm.FIT_NETWORKS: False,
        algorithm.WRITE_SEGMENT_LINES: True,
    })
    captured = {}

    def fake_run(tool, values, folder, prefix, feedback):
        captured.update(values)
        _touch_outputs(folder, prefix, [
            "_chi_disorder_sweep.csv",
            "_chi_fit_segment_lines_mn_0p45.fgb",
        ])

    monkeypatch.setattr(module, "run_lsdtt_tool", fake_run)
    feedback = _Feedback()
    outputs = algorithm.processAlgorithm(parameters, None, feedback)
    assert captured["fit_chi_networks"] is True
    assert any("enabling fit_chi_networks" in message for message in feedback.messages)
    assert outputs[algorithm.OUTPUT_SEGMENT_LINES].endswith(
        "run_chi_fit_segment_lines_mn_0p45.fgb"
    )

    def fake_run_without_segments(tool, values, folder, prefix, feedback):
        _touch_outputs(folder, prefix, ["_chi_disorder_sweep.csv"])

    os.remove(outputs[algorithm.OUTPUT_SEGMENT_LINES])
    monkeypatch.setattr(module, "run_lsdtt_tool", fake_run_without_segments)
    with pytest.raises(RuntimeError, match="Expected chi fit segment lines output not found"):
        algorithm.processAlgorithm(parameters, None, _Feedback())
