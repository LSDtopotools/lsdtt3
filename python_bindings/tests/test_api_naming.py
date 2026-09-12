import importlib.util
import json
from pathlib import Path


def load_checker():
    path = Path(__file__).parents[1] / "tools" / "check_api_naming.py"
    spec = importlib.util.spec_from_file_location("check_api_naming", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_api_naming_matches_reviewed_baseline():
    root = Path(__file__).parents[1]
    manifest = json.loads((root / "api" / "public_api.json").read_text(encoding="utf-8"))
    baseline = json.loads((root / "api" / "naming_baseline.json").read_text(encoding="utf-8"))
    assert load_checker().collect_violations(manifest) == baseline["violations"]
