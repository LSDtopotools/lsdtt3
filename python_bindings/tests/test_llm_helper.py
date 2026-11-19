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
from lsdtt3._llm import _GUIDANCE, WORKFLOW_CODE

pytestmark = [pytest.mark.conformance, pytest.mark.unit]


def test_index_is_small_and_exposes_retrieval_protocol():
    """The starter card stays small, says how to work, and lists every topic once."""
    card = lsdtt3.helper_for_llm()
    assert card.startswith("LSDTT3/API3")
    for text in ("How to work:", "helper_for_llm('gotchas')", "'symbol:<name>'",
                 "'symbol:<Class>.<method>'", "'search:<words>'", "'all'", "'catalogue'",
                 "[experimental]"):
        assert text in card, text
    block = card.split("\nTopics:\n", 1)[1].split("\nAlso:", 1)[0]
    topics = [re.match(r"  (\w+): ", line).group(1) for line in block.splitlines()]
    assert sorted(topics) == sorted(_GUIDANCE)
    assert topics[:2] == ["workflow", "gotchas"]
    assert len(card) < 2500


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
    for text in ("mask_to_nodata(R,0.0)", "cp312-abi3", "R.write(path)", "R.metadata.crs_wkt",
                 "basin_id+1", "silently drops duplicate", "F.row_col(flat_index)",
                 "m^2", "pixel count", "NORTH edge", "basin_index", "k_sn",
                 "PROJ_DATA", "lsdtt-*", "F.snap_points(", "M.xy_to_row_col(x,y)",
                 "M.row_col_to_xy(row,col)", "threshold_drainage_area_m2"):
        assert text in gotchas, text
    # The workflow no longer hand-rolls the coordinate arithmetic.
    for text in ("Transformer", "// meta.cell_size", "receiver_of"):
        assert text not in code, text
    everything = lsdtt3.helper_for_llm("all")
    for topic in _GUIDANCE.values():
        assert topic in everything
    assert everything.index("@workflow") < everything.index("@terrain")
    # Kept small enough to paste into an assistant along with a question
    # (limit 40000, Simon 2026-11-19, plan_llm_guide_restructure.md).
    assert len(everything) < 40_000


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
    for topic in _GUIDANCE:
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


@pytest.mark.negative
def test_unknown_topic_or_symbol_is_actionable():
    """Typos tell an agent how to recover rather than returning empty context."""
    with pytest.raises(ValueError, match="Unknown lsdtt3 LLM topic"):
        lsdtt3.helper_for_llm("volcanoes")
    with pytest.raises(ValueError, match="Unknown lsdtt3 public symbol"):
        lsdtt3.helper_for_llm("symbol:not_a_binding")


def test_high_risk_cards_match_runtime_contracts():
    """python_bindings/tools/check_llm_cards.py passes: the LLM helper's cards match the
    runtime contracts.
    """
    path = Path(__file__).parents[1] / "tools" / "check_llm_cards.py"
    spec = importlib.util.spec_from_file_location("check_llm_cards", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.main() == 0


def test_guide_status_lists_name_real_exports():
    """READY and EXPERIMENTAL_OVERRIDES name only public exports and their members."""
    from lsdtt3._llm_status import EXPERIMENTAL_OVERRIDES, READY

    public = set(lsdtt3.__all__)
    assert READY <= public, sorted(READY - public)
    for name in EXPERIMENTAL_OVERRIDES:
        owner, _, member = name.partition(".")
        assert owner in public, name
        assert not member or hasattr(getattr(lsdtt3, owner), member), name


def test_guide_status_resolution(monkeypatch):
    """Ready if listed, experimental otherwise; methods inherit; overrides demote."""
    from lsdtt3 import _llm
    from lsdtt3._llm import _guide_status

    namespace = vars(lsdtt3)
    assert _guide_status("fill_pits", namespace) == "ready"
    assert _guide_status("FlowInfo", namespace) == "ready"
    assert _guide_status("FlowInfo.low_concavity", namespace) == "ready"
    assert _guide_status("slope_area_bins", namespace) == "ready"
    # mentioned only as a polyfit metric string in the cards, so not frozen as ready
    assert _guide_status("plan_curvature", namespace) == "experimental"
    statuses = {name: _guide_status(name, namespace) for name in lsdtt3.__all__}
    assert set(statuses.values()) == {"ready", "experimental"}
    assert sum(value == "ready" for value in statuses.values()) == len(_llm.READY)
    with pytest.raises(ValueError, match="Unknown lsdtt3 public symbol"):
        _guide_status("no_such_export", namespace)
    with pytest.raises(ValueError, match="has no member"):
        _guide_status("FlowInfo.no_such_method", namespace)
    monkeypatch.setattr(_llm, "EXPERIMENTAL_OVERRIDES", frozenset({"FlowInfo.low_concavity", "fill_pits"}))
    assert _guide_status("FlowInfo.low_concavity", namespace) == "experimental"
    assert _guide_status("FlowInfo.d8_flow_path_gradient", namespace) == "ready"
    assert _guide_status("fill_pits", namespace) == "experimental"


def _manifest_exports():
    path = Path(__file__).parents[1] / "api" / "public_api.json"
    import json
    return {entry["name"]: entry for entry in json.loads(path.read_text(encoding="utf-8"))["exports"]}


def test_every_export_has_one_guide_topic():
    """guide_topic gives every export one catalogue topic, and the placements are not stale."""
    from lsdtt3._topics import GUIDE_TOPIC_PLACEMENTS, GUIDE_TOPICS, guide_topic

    assert set(GUIDE_TOPICS) == set(_GUIDANCE) - {"workflow", "gotchas"}
    assert set(GUIDE_TOPIC_PLACEMENTS) <= set(lsdtt3.__all__), sorted(
        set(GUIDE_TOPIC_PLACEMENTS) - set(lsdtt3.__all__))
    assert set(GUIDE_TOPIC_PLACEMENTS.values()) <= set(GUIDE_TOPICS)
    for name in lsdtt3.__all__:
        assert guide_topic(name) in GUIDE_TOPICS, name


def _catalogue_heads(topic):
    """Export names catalogued in a topic, and the member names under each class."""
    text = lsdtt3.helper_for_llm(topic).split(f"@catalogue {topic}:", 1)[1].splitlines()[1:]
    heads, members, current = [], {}, None
    for line in text:
        if line.startswith("  ."):
            members[current].append(re.match(r"  \.(\w+)", line).group(1))
        elif line.strip():
            current = re.match(r"(\w+)", line).group(1)
            heads.append(current)
            members[current] = []
    return heads, members


def test_catalogues_list_every_ready_export_once_and_no_experimental_one():
    """Coverage (plan 5.1): each ready export is in its topic's catalogue, once;
    experimental exports are in none."""
    from lsdtt3._llm import _guide_status
    from lsdtt3._topics import GUIDE_TOPICS, guide_topic

    namespace = vars(lsdtt3)
    seen = {}
    for topic in GUIDE_TOPICS:
        heads, _members = _catalogue_heads(topic)
        for name in heads:
            assert name not in seen, f"{name} in {seen.get(name)} and {topic}"
            seen[name] = topic
    ready = {name for name in lsdtt3.__all__ if _guide_status(name, namespace) == "ready"}
    assert set(seen) == ready, (sorted(ready - set(seen)), sorted(set(seen) - ready))
    assert all(seen[name] == guide_topic(name) for name in seen)


def test_catalogues_list_every_public_member_of_a_ready_class():
    """Methods (plan 5.2): a ready class lists the members the API manifest has
    (except __init__), unless one is demoted on its own."""
    from lsdtt3._topics import GUIDE_TOPICS

    exports = _manifest_exports()
    for topic in GUIDE_TOPICS:
        _heads, members = _catalogue_heads(topic)
        for name, listed in members.items():
            if exports[name]["kind"] != "class":
                assert not listed, name
                continue
            expected = sorted(member["name"] for member in exports[name]["members"]
                              if member["name"] != "__init__")
            assert sorted(listed) == expected, name


def test_catalogue_lines_carry_argument_names_and_a_summary():
    """A catalogue line is name(args): summary; ? marks a parameter with a default."""
    channel = lsdtt3.helper_for_llm("channel")
    assert "slope_area_bins(drainage_area,slope,*,n_log_bins?,exclude?): Bin slope by drainage area" in channel
    flow = lsdtt3.helper_for_llm("flow")
    assert "  .low_concavity(drainage_area,gradient,*,threshold_contributing_pixels?," in flow
    assert "\nChannelScope: " in channel and "ChannelScope(*args" not in channel
    assert "FlowInfo(dem,boundary_conditions?," in flow
    # experimental exports are left out (aspect is not READY)
    assert "\naspect(" not in lsdtt3.helper_for_llm("terrain")


def test_symbol_works_for_every_export_and_public_member():
    """Lookup (plan 5.5): symbol: gives a card for every export and every member
    the manifest lists; an experimental one says so; unknown names raise."""
    from lsdtt3._llm import _EXPERIMENTAL_NOTE, _guide_status

    namespace = vars(lsdtt3)
    for name, entry in _manifest_exports().items():
        names = [name] + [f"{name}.{member['name']}" for member in entry["members"]
                          if member["name"] != "__init__"]
        for symbol in names:
            card = lsdtt3.helper_for_llm(f"symbol:{symbol}")
            assert card.startswith(f"@{symbol}\n"), symbol
            experimental = _guide_status(symbol, namespace) == "experimental"
            assert (_EXPERIMENTAL_NOTE in card) == experimental, symbol
    method = lsdtt3.helper_for_llm("symbol:FlowInfo.low_concavity")
    assert method.startswith("@FlowInfo.low_concavity\nlow_concavity(self, drainage_area: Raster,")
    assert "run_low_concavity_detection" in method
    assert "Authors" not in method
    assert lsdtt3.helper_for_llm("symbol:BasinCollection.table").startswith(
        "@BasinCollection.table\nOne row per basin")
    with pytest.raises(ValueError, match="has no public member"):
        lsdtt3.helper_for_llm("symbol:FlowInfo.no_such_method")
    with pytest.raises(ValueError, match="Unknown lsdtt3 public symbol"):
        lsdtt3.helper_for_llm("symbol:no_such_export.table")


@pytest.mark.parametrize("query, first", [
    ("low concavity", "FlowInfo.low_concavity"),
    ("run_low_concavity_detection", "FlowInfo.low_concavity"),
    ("fill pits", "fill_pits"),
    ("basins", "FlowInfo.basins"),
    ("hillshade", "hillshade"),
    ("burial", "[experimental] crn_burial_pair"),
])
def test_search_puts_the_specific_routine_first(query, first):
    """search: indexes methods, matches CLI parameter names, and ranks the
    specific routine above the class holding it; experimental entries are marked."""
    lines = lsdtt3.helper_for_llm(f"search:{query}").splitlines()
    assert lines[1].split(": ", 1)[0] == first, lines[:4]


def test_search_finds_a_function_behind_its_result_class():
    lines = lsdtt3.helper_for_llm("search:perimeter sources").splitlines()[1:3]
    assert "[experimental] find_perimeter_sources" in [line.split(": ", 1)[0] for line in lines]


def test_search_finds_the_routine_documenting_a_cli_parameter():
    lines = lsdtt3.helper_for_llm("search:basin_outlet_snapping_method").splitlines()[1:3]
    assert "FlowInfo.snap_points" in [line.split(": ", 1)[0] for line in lines]


def test_sizes_stay_within_the_limits():
    """Size (plan 5.3): index < 2500, each topic < 8000, all < 40000 characters.
    A topic over the limit is split (or its guidance shortened), never its catalogue."""
    assert len(lsdtt3.helper_for_llm("index")) < 2500
    for topic in _GUIDANCE:
        size = len(lsdtt3.helper_for_llm(topic))
        assert size < 8000, f"topic {topic} is {size} characters: split it in lsdtt3/_topics.py"
    size = len(lsdtt3.helper_for_llm("all"))
    assert size < 40_000, f"all is {size} characters: shorten the hand-written guidance"


def test_all_names_every_ready_export_and_the_key_class_methods():
    """all = index + guidance + every ready export's name by topic + key class methods."""
    from lsdtt3._llm import _ALL_KEY_CLASSES, _guide_status
    from lsdtt3._topics import GUIDE_TOPICS, guide_topic

    everything = lsdtt3.helper_for_llm("all")
    assert everything.startswith(lsdtt3.helper_for_llm("index"))
    names = everything.split("\n@names ", 1)[1].split("\n@methods ", 1)[0].splitlines()[1:]
    listed = {line.split(": ", 1)[0]: line.split(": ", 1)[1].split() for line in names}
    assert list(listed) == list(GUIDE_TOPICS)
    namespace = vars(lsdtt3)
    for name in lsdtt3.__all__:
        ready = _guide_status(name, namespace) == "ready"
        assert (name in listed[guide_topic(name)]) == ready, name
    methods = everything.split("\n@methods ", 1)[1].split("\n@end ", 1)[0].splitlines()[1:]
    assert [line.split(": ", 1)[0] for line in methods] == list(_ALL_KEY_CLASSES)
    assert "low_concavity" in methods[0].split(": ", 1)[1].split()
    assert everything.rstrip().endswith("paste the result.")
    assert "@catalogue" not in everything


def test_catalogue_mode_joins_every_topic_catalogue():
    from lsdtt3._topics import GUIDE_TOPICS

    catalogue = lsdtt3.helper_for_llm("catalogue")
    for topic in GUIDE_TOPICS:
        assert f"@catalogue {topic}:" in catalogue, topic
    assert "slope_area_bins(drainage_area,slope,*,n_log_bins?,exclude?)" in catalogue
    assert lsdtt3.helper_for_llm("catalog") == catalogue
