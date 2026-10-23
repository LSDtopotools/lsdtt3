"""Tests for compact LLM API context retrieval.

Authors
-------
B. Gailleton
S. M. Mudd
"""

import ast
import importlib.util
from pathlib import Path
import re
import textwrap

import pytest

import lsdtt3
from lsdtt3._llm import _TOPICS, WORKFLOW_CODE


def test_index_is_small_and_exposes_retrieval_protocol():
    """The starter card stays small enough for routine prompt injection."""
    card = lsdtt3.helper_for_llm()
    assert card.startswith("LSDTT3/API2")
    assert "'symbol:<exact name>'" in card and "'search:<words>'" in card and "'all'" in card
    topics = re.search(r"^Topics: (.*)$", card, re.MULTILINE).group(1).split()
    assert sorted(topics) == sorted(_TOPICS)
    assert topics[:2] == ["workflow", "gotchas"]
    assert len(card.split()) < 160


def test_helper_returns_text_that_displays_as_text():
    """In a notebook the last expression is shown with repr(): keep it readable."""
    card = lsdtt3.helper_for_llm()
    assert isinstance(card, str)
    assert repr(card) == str(card)
    assert "\\n" not in repr(card)
    assert repr(lsdtt3.helper_for_llm("all")).count("\n") > 50


def test_topic_card_keeps_exact_code_identifiers():
    """A topic card has the API names required for generated Python code."""
    card = lsdtt3.helper_for_llm("flow")
    assert "FlowInfo" in card
    assert "N=F.channel_network(S)" in card
    assert "F.basins(outlets?,min_upslope_pixels=0)->B" in card
    assert "F.receiver_of(row,col)->(row,col) or None at an outlet" in card


def test_symbol_cards_give_the_full_documentation():
    """symbol: returns the whole docstring, signature first for functions."""
    symbol = lsdtt3.helper_for_llm("symbol:fill_pits")
    assert symbol.startswith("@fill_pits\nfill_pits(raster: Raster, min_slope: float = 0.0001")
    assert "Parameters\n----------\nraster : Raster" in symbol
    assert "Authors" not in symbol
    flow = lsdtt3.helper_for_llm("symbol:FlowInfo")
    assert "FlowInfo(dem: Raster, boundary_conditions: str = 'oooo'" in flow
    chi_fit = lsdtt3.helper_for_llm("symbol:ChiFit")
    for member in ("nodes", "segments", "summary", "invert_ksn("):
        assert re.search(rf"^{re.escape(member)}", chi_fit, re.MULTILINE), member
    assert "ChiCollection.fit()" in chi_fit
    assert len(chi_fit) > 400


def test_search_ranks_names_and_docstrings():
    """search: matches any query word in names or docstrings, best first."""
    search = lsdtt3.helper_for_llm("search:flow")
    assert "FlowGraph" in search
    assert "FlowInfo" in search
    chi_ksn = lsdtt3.helper_for_llm("search:chi ksn")
    names = re.findall(r"^(\w+): ", chi_ksn, re.MULTILINE)
    assert names, chi_ksn
    assert "ChiFit" in names[:5]
    assert all("chi" in name.lower() or "ksn" in name.lower() for name in names[:3])
    assert "none" in lsdtt3.helper_for_llm("search:zzzqqq")


def test_workflow_and_gotchas_topics():
    """The workflow is runnable Python and the gotchas cover the known traps."""
    workflow = lsdtt3.helper_for_llm("workflow")
    code = re.search(r"```python\n(.*?)```", workflow, re.DOTALL).group(1)
    assert code == WORKFLOW_CODE
    compile(ast.parse(code), "workflow", "exec")
    for step in ("lsdtt3.fetch_dem(", "lsdtt3.mask_to_nodata(dem, 0.0)", "lsdtt3.fill_pits(",
                 "lsdtt3.FlowInfo(", "F.channel_scope(threshold_contributing_pixels=",
                 "F.snap_points(", 'method="downstream"',
                 "threshold_contributing_pixels=threshold,",
                 "F.basins(outlet_pixels)", "B.polygons",
                 "critical_divergence_metres=10.0", 'fit.summary[["basin_id", "median_ksn"]]',
                 '.write("'):
        assert step in code, step
    gotchas = lsdtt3.helper_for_llm("gotchas")
    for text in ("mask_to_nodata(R,0.0)", "cp313", "R.write(path)", "R.metadata.crs_wkt",
                 "basin_id+1", "silently drops duplicate", "F.row_col(flat_index)",
                 "m^2", "pixel count", "NORTH edge", "basin_index", "k_sn",
                 "PROJ_DATA", "lsdtt-*", "F.snap_points(", "M.xy_to_row_col(x,y)",
                 "M.row_col_to_xy(row,col)", "threshold_drainage_area_m2"):
        assert text in gotchas, text
    # The workflow no longer hand-rolls the coordinate arithmetic.
    for text in ("Transformer", "// meta.cell_size", "receiver_of"):
        assert text not in code, text
    everything = lsdtt3.helper_for_llm("all")
    for topic in _TOPICS.values():
        assert topic in everything
    assert everything.index("@workflow") < everything.index("@terrain")
    # Kept small enough to paste into an assistant along with a question.
    assert len(everything) < 25_000


def test_module_docstring_is_a_task_oriented_entry_point():
    """help(lsdtt3) explains the domains, conventions and the core workflow."""
    doc = lsdtt3.__doc__
    for text in ("Conventions", "Core workflow", "row 0 is the NORTH edge",
                 "x = meta.x_min + (col + 0.5) * meta.cell_size",
                 "y = meta.y_max - (row + 0.5) * meta.cell_size",
                 "R.write(path)", "zero-based", 'helper_for_llm("all")',
                 "FlowInfo.snap_points", "meta.xy_to_row_col(x, y)",
                 "Authors\n-------\nB. Gailleton"):
        assert text in doc, text
    for topic in _TOPICS:
        assert f'"{topic}"' in doc, topic
    code = doc.split("ksn::\n\n", 1)[1].split("\nHelp for AI assistants", 1)[0]
    assert textwrap.dedent(code).strip() == WORKFLOW_CODE.strip()
    # Every function, class and member the domain list names is real.
    domains = doc.split("What it does", 1)[1].split("Conventions", 1)[0]
    for class_name, member in re.findall(r"\b([A-Z]\w+)\.([a-z_]\w*)", domains):
        assert hasattr(getattr(lsdtt3, class_name), member), f"{class_name}.{member}"
    for name in re.findall(r"(?<![.\w])([a-z]+_[a-z0-9_]+)\b", domains):
        assert name in lsdtt3.__all__, name
    for name in re.findall(r"\b([A-Z][a-z]+(?:[A-Z][a-z]*)+)\b", domains):
        assert name in lsdtt3.__all__ or name == "GraphFlood", name


def test_unknown_topic_or_symbol_is_actionable():
    """Typos tell an agent how to recover rather than returning empty context."""
    with pytest.raises(ValueError, match="Unknown lsdtt3 LLM topic"):
        lsdtt3.helper_for_llm("volcanoes")
    with pytest.raises(ValueError, match="Unknown lsdtt3 public symbol"):
        lsdtt3.helper_for_llm("symbol:not_a_binding")


def test_high_risk_cards_match_runtime_contracts():
    path = Path(__file__).parents[1] / "tools" / "check_llm_cards.py"
    spec = importlib.util.spec_from_file_location("check_llm_cards", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.main() == 0
