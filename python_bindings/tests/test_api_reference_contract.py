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


def test_package_and_reference_share_one_topic_rule():
    """lsdtt3._topics is the one domain rule: the API reference generator uses it, and
    the installed package's copy gives every export the same domain (plan_llm_guide_restructure.md,
    step 1).
    """
    import lsdtt3
    from lsdtt3 import _topics

    repository = Path(__file__).parents[2]
    generator_path = repository / "python_bindings/tools/generate_api_reference.py"
    spec = importlib.util.spec_from_file_location("api_reference_generator", generator_path)
    generator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(generator)
    assert generator.DOMAINS == _topics.DOMAINS
    for name in lsdtt3.__all__:
        assert generator.domain(name) == _topics.domain(name), name
        assert _topics.domain(name) in _topics.DOMAINS, name
