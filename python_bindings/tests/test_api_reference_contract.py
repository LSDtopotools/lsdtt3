import json
import importlib.util
from pathlib import Path
import re

import pytest

pytestmark = [pytest.mark.conformance, pytest.mark.unit]


def test_generated_reference_assigns_every_export_exactly_once():
    """The generated docs/python_api pages are exactly the generator's domains, and every
    export in the manifest appears on exactly one page, the one its domain says.
    """
    repository = Path(__file__).parents[2]
    manifest = json.loads(
        (repository / "python_bindings/api/public_api.json").read_text(encoding="utf-8")
    )
    generator_path = repository / "python_bindings/tools/generate_api_reference.py"
    spec = importlib.util.spec_from_file_location("api_reference_generator", generator_path)
    assert spec is not None and spec.loader is not None
    generator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(generator)
    pages = {
        page.stem: page for page in (repository / "docs/python_api").glob("*.md")
        if page.name != "index.md"
    }
    assert set(pages) == set(generator.DOMAINS)
    names = []
    assignments = {}
    for page_name, page in pages.items():
        page_names = re.findall(r"^## `([^`]+)`$", page.read_text(encoding="utf-8"), re.MULTILINE)
        names.extend(page_names)
        assignments.update((name, page_name) for name in page_names)
    assert sorted(names) == sorted(export["name"] for export in manifest["exports"])
    assert len(names) == len(set(names))
    assert all(assignments[name] == generator.domain(name) for name in names)
