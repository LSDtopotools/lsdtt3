"""Re-running an algorithm into the same output folder and prefix.

QGIS keeps loaded layers' files open, and on Windows an open file cannot be
replaced, so a re-run has to close those layers first (prepareAlgorithm), and
an output the run could not replace must not be passed off as new.
"""

import os
import time

import pytest

from test_phase4_algorithm_contract import (  # noqa: F401  (installs the qgis stubs)
    _Feedback,
    _Layer,
    _configure_parameters,
    _load,
)
from lsdtt_qgis.naming import existing_prefixed_outputs, normalised_path, written_since


def _touch(path, mtime=None):
    with open(path, "w", encoding="utf-8"):
        pass
    if mtime is not None:
        os.utime(path, (mtime, mtime))
    return str(path)


def test_existing_prefixed_outputs_ignores_other_prefixes_and_bookkeeping(tmp_path):
    for name in ["run_fill.tif", "run_hillshade.tif", "run_raster_preprocessing.log",
                 "dem_qgis_analysis_001.driver", "run.tif", "runner_fill.tif", "other_fill.tif"]:
        _touch(tmp_path / name)
    found = [os.path.basename(p) for p in existing_prefixed_outputs(str(tmp_path), "run")]
    assert found == ["run_fill.tif", "run_hillshade.tif"]
    assert existing_prefixed_outputs(str(tmp_path / "missing"), "run") == []


def test_written_since_detects_a_leftover_file(tmp_path):
    started = time.time()
    old = _touch(tmp_path / "old.tif", mtime=started - 3600)
    new = _touch(tmp_path / "new.tif")
    assert not written_since(old, started)
    assert written_since(new, started)


class _MapLayer:
    def __init__(self, layer_id, source):
        self._id, self._source = layer_id, source

    def id(self):
        return self._id

    def name(self):
        return self._id

    def source(self):
        return self._source


class _Project:
    def __init__(self, layers):
        self.layers = {layer.id(): layer for layer in layers}
        self.removed = []

    def mapLayers(self):
        return dict(self.layers)

    def removeMapLayers(self, ids):
        self.removed.extend(ids)
        for layer_id in ids:
            self.layers.pop(layer_id)


class _Context:
    def __init__(self, project):
        self._project = project

    def project(self):
        return self._project


def _algorithm(monkeypatch):
    module, algorithm_class = _load("raster_preprocessing", "RasterPreprocessingAlgorithm")
    monkeypatch.setattr(module, "attach_backend_help", lambda *args, **kwargs: None)
    algorithm = algorithm_class()
    algorithm.initAlgorithm()
    _configure_parameters(algorithm)
    algorithm.parameterAsLayer = lambda params, key, context: params.get(key)
    algorithm.parameterAsLayerList = lambda params, key, context: params.get(key) or []
    for pdef in algorithm.parameterDefinitions():
        kind = "raster" if pdef.name() in (algorithm.INPUT_DEM, algorithm.BINARY_MASK) else (
            "vector" if pdef.name() == algorithm.REPLACE_PIXELS else "other")
        pdef.type = (lambda k: (lambda: k))(kind)
    return algorithm


def test_every_algorithm_offers_overwrite_on_by_default():
    import glob
    from lsdtt_qgis.algorithms.support import OVERWRITE_OUTPUTS
    folder = os.path.join(os.path.dirname(__file__), "..", "algorithms")
    modules = [os.path.splitext(os.path.basename(p))[0] for p in glob.glob(os.path.join(folder, "*.py"))
               if os.path.basename(p) not in ("__init__.py", "support.py")]
    import importlib
    checked = 0
    for name in modules:
        module = importlib.import_module("lsdtt_qgis.algorithms." + name)
        module.attach_backend_help = lambda *args, **kwargs: None
        classes = [v for v in vars(module).values() if isinstance(v, type)
                   and v.__module__ == module.__name__ and hasattr(v, "createInstance")]
        for cls in classes:
            algorithm = cls()
            algorithm.initAlgorithm()
            overwrite = [p for p in algorithm.parameterDefinitions() if p.name() == OVERWRITE_OUTPUTS]
            assert len(overwrite) == 1, cls.__name__
            assert overwrite[0].defaultValue() is True, cls.__name__
            checked += 1
    # 14 algorithm files; drainage_analysis holds five algorithms and its base.
    assert checked >= 18, checked


def test_rerun_closes_the_layers_it_will_overwrite(tmp_path, monkeypatch):
    algorithm = _algorithm(monkeypatch)
    fill = _touch(tmp_path / "run_fill.tif")
    _touch(tmp_path / "dem.tif")
    project = _Project([
        _MapLayer("previous fill", fill),
        _MapLayer("the input", str(tmp_path / "dem.tif")),
        _MapLayer("elsewhere", str(tmp_path / "other_fill.tif")),
    ])
    params = {
        algorithm.INPUT_DEM: _Layer(str(tmp_path / "dem.tif")),
        algorithm.OUTPUT_FOLDER: str(tmp_path),
        algorithm.OUTPUT_PREFIX: "run",
        "OVERWRITE_OUTPUTS": True,
    }
    feedback = _Feedback()
    assert algorithm.prepareAlgorithm(params, _Context(project), feedback) is True
    assert project.removed == ["previous fill"]
    assert any("previous fill" in message for message in feedback.messages)


def test_overwrite_off_refuses_to_clobber(tmp_path, monkeypatch):
    algorithm = _algorithm(monkeypatch)
    _touch(tmp_path / "run_fill.tif")
    params = {
        algorithm.INPUT_DEM: _Layer(str(tmp_path / "dem.tif")),
        algorithm.OUTPUT_FOLDER: str(tmp_path),
        algorithm.OUTPUT_PREFIX: "run",
        "OVERWRITE_OUTPUTS": False,
    }
    with pytest.raises(RuntimeError, match="tick 'Overwrite existing outputs'"):
        algorithm.prepareAlgorithm(params, _Context(_Project([])), _Feedback())
    # A fresh prefix is fine with the option off.
    params[algorithm.OUTPUT_PREFIX] = "run2"
    assert algorithm.prepareAlgorithm(params, _Context(_Project([])), _Feedback()) is True


def test_an_input_that_is_also_an_output_is_refused(tmp_path, monkeypatch):
    algorithm = _algorithm(monkeypatch)
    previous = _touch(tmp_path / "run_fill.tif")
    params = {
        algorithm.INPUT_DEM: _Layer(previous),
        algorithm.OUTPUT_FOLDER: str(tmp_path),
        algorithm.OUTPUT_PREFIX: "run",
        "OVERWRITE_OUTPUTS": True,
    }
    project = _Project([_MapLayer("input", previous)])
    with pytest.raises(RuntimeError, match="Choose another output prefix"):
        algorithm.prepareAlgorithm(params, _Context(project), _Feedback())
    assert project.removed == []


def test_require_output_rejects_a_file_the_run_did_not_rewrite(tmp_path, monkeypatch):
    from lsdtt_qgis.algorithms import support
    monkeypatch.setitem(support._RUN_STARTED, support._run_key(str(tmp_path), "run"), time.time())
    _touch(tmp_path / "run_fill.tif", mtime=time.time() - 3600)
    with pytest.raises(RuntimeError, match="not rewritten by this run"):
        support.require_output(str(tmp_path), "run", "_fill.tif", "fill")
    _touch(tmp_path / "run_fill.tif")
    assert normalised_path(support.require_output(str(tmp_path), "run", "_fill.tif", "fill")) == \
        normalised_path(str(tmp_path / "run_fill.tif"))
