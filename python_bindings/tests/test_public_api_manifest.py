import importlib.util
import importlib.metadata
import json
from pathlib import Path
import re


def load_generator():
    path = Path(__file__).parents[1] / "tools" / "generate_public_api_manifest.py"
    spec = importlib.util.spec_from_file_location("public_api_manifest", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_public_api_manifest_matches_runtime():
    root = Path(__file__).parents[1]
    expected = json.loads((root / "api" / "public_api.json").read_text(encoding="utf-8"))
    actual = load_generator().build_manifest()
    if actual != expected:
        raise AssertionError(
            "public API manifest differs from the installed runtime; "
            "run generate_public_api_manifest.py --update after reviewing the change"
        )


def test_api_decision_manifest_has_release_gate():
    root = Path(__file__).parents[1]
    decisions = json.loads((root / "api" / "decisions.json").read_text(encoding="utf-8"))
    assert decisions["schema_version"] == 1
    load_generator().validate_decisions(
        load_generator().build_manifest(), decisions
    )


def test_runtime_exposes_package_version():
    import lsdtt3

    try:
        installed = importlib.metadata.distribution("lsdtt3")
        expected_path = Path(installed.locate_file("lsdtt3")).resolve()
        expected = installed.version if expected_path == Path(lsdtt3.__file__).resolve().parent else "0.0.0"
    except importlib.metadata.PackageNotFoundError:
        expected = "0.0.0"
    assert lsdtt3.__version__ == expected


def test_graphflood_backend_option_schema_matches_runtime_manifest():
    root = Path(__file__).parents[1]
    manifest = json.loads((root / "api" / "public_api.json").read_text(encoding="utf-8"))
    schema = json.loads(
        (root / "api" / "graphflood_backend_options.json").read_text(encoding="utf-8")
    )
    assert schema["schema_version"] == 1
    assert set(schema["accelerators"]) == {
        "run_graphflood_multigrid", "run_graphflood_induced"
    }
    assert schema["migration"]["null_means"]
    assert schema["migration"]["legacy_conflicts_are_checked_before_backend_invocation"] is True
    assert schema["custom_callable_policy"]
    exports = {entry["name"]: entry for entry in manifest["exports"]}
    required_inputs = {"bed", "node_state", "rainfall_rate", "friction_coefficient"}
    assert set(schema["accelerators"]).issubset(exports)
    assert set(schema["migration"]["preferred_public_tolerance_parameters"]) == {
        "median_dhw_convergence_tolerance_metres",
        "median_dhw_relative_convergence_tolerance",
    }
    for backend in schema["built_in_backends"].values():
        assert len(backend["options"]) == len(set(backend["options"]))
        definition = exports[backend["callable"]]["signatures"][0]["definition"]
        runtime_names = set(re.findall(r"(?:\(|, )([A-Za-z_]\w*):", definition))
        assert required_inputs.issubset(runtime_names)
        assert runtime_names - required_inputs == set(backend["options"])


def test_manifest_distinguishes_static_methods_and_absent_constructors():
    manifest = load_generator().build_manifest()
    exports = {entry["name"]: entry for entry in manifest["exports"]}
    raster_members = {member["name"]: member for member in exports["Raster"]["members"]}
    assert raster_members["read"]["kind"] == "static_method"
    assert raster_members["from_numpy"]["kind"] == "static_method"
    assert raster_members["to_numpy"]["kind"] == "method"
    assert exports["BasinCollection"]["signatures"] == []
    assert all(
        member["name"] != "__init__"
        for member in exports["BasinCollection"]["members"]
    )
