import ast
import json
from pathlib import Path
import re
import importlib.util


def load_signature_generator():
    path = Path(__file__).parents[1] / "tools" / "generate_documentation_signature_contract.py"
    spec = importlib.util.spec_from_file_location("documentation_signatures", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_documentation_exclusions_exactly_cover_missing_facade_symbols():
    root = Path(__file__).parents[1]
    manifest = json.loads((root / "api" / "public_api.json").read_text(encoding="utf-8"))
    exclusions = json.loads(
        (root / "api" / "documentation_exclusions.json").read_text(encoding="utf-8")
    )
    assert exclusions["schema_version"] == 1
    assert exclusions["symbols"] == manifest["documentation_facade"]["missing"]


def test_python_api_pages_reference_real_facade_symbols():
    repository = Path(__file__).parents[2]
    facade_path = repository / "python_bindings" / "src" / "lsdtt3" / "_api_docs.py"
    tree = ast.parse(facade_path.read_text(encoding="utf-8"))
    facade_symbols = {
        node.name for node in tree.body
        if isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef))
    }
    referenced = set()
    for page in (repository / "docs" / "python_api").glob("*.md"):
        referenced.update(re.findall(r"::: lsdtt3\._api_docs\.([A-Za-z_]\w*)", page.read_text(encoding="utf-8")))
    assert referenced
    assert referenced.issubset(facade_symbols)


def test_documentation_signature_contract_is_current():
    root = Path(__file__).parents[1]
    expected = json.loads(
        (root / "api" / "documentation_signatures.json").read_text(encoding="utf-8")
    )
    assert load_signature_generator().build_contract(root) == expected
    exceptions = json.loads(
        (root / "api" / "documentation_signature_exceptions.json").read_text(encoding="utf-8")
    )
    mismatches = [
        record["symbol"] for record in expected["symbols"]
        if not record["matches_runtime_signature"]
    ]
    assert exceptions["symbols"] == mismatches
