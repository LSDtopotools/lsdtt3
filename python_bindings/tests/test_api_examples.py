import importlib.util
from pathlib import Path


def test_examples_compile_and_use_public_exports():
    path = Path(__file__).parents[1] / "tools" / "check_api_examples.py"
    spec = importlib.util.spec_from_file_location("check_api_examples", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.main([]) == 0
