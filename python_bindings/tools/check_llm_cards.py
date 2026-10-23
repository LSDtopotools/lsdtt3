#!/usr/bin/env python3

"""Validate LLM API cards against current public contracts.

Every function, constructor and method a card mentions must exist in the
public API manifest, and every argument name a card gives it (``name=...``,
``name:Type``, ``name?`` or a bare lower-case positional name) must be one of
its real parameters. The runnable ``workflow`` code is parsed and checked the
same way, following the types its calls return.
"""

import ast
import json
from pathlib import Path
import re

import lsdtt3
from lsdtt3._llm import _INDEX, _TOPICS, WORKFLOW_CODE

# Card aliases for public classes (see the notation line of the index card).
ALIASES = {
    "R": "Raster", "M": "RasterMeta", "F": "FlowInfo", "G": "FlowGraph",
    "H": "RasterHierarchy", "S": "ChannelScope", "N": "ChannelNetwork",
    "B": "BasinCollection", "C": "ChiCollection", "L": "LandscapeModel",
    "Q": "CRNColumn", "fit": "ChiFit", "inv": "ChiInversion",
    "b": "VectorBounds",
}

FUNCTION_CALL = re.compile(
    r"(?<![^\s;=,(])([a-z][a-z0-9_]*(?:\|[a-z][a-z0-9_]*)*)\(([^()]*)\)", re.MULTILINE
)
CONSTRUCTOR_CALL = re.compile(r"(?<![\w.])([A-Z][A-Za-z0-9]+)\(([^()]*)\)")
MEMBER_CALL = re.compile(r"(?<![\w.])([A-Za-z_]\w*)\.([a-z_]\w*)\(([^()]*)\)")
MEMBER_ACCESS = re.compile(r"(?<![\w.])([A-Za-z_]\w*)\.([a-z_]\w*)\b(?!\()")


def split_arguments(text):
    """Split card arguments at commas outside brackets and quotes."""
    parts, depth, start, quote = [], 0, 0, None
    for index, character in enumerate(text):
        if quote:
            if character == quote:
                quote = None
        elif character in "'\"":
            quote = character
        elif character in "([{":
            depth += 1
        elif character in ")]}":
            depth -= 1
        elif character == "," and depth == 0:
            parts.append(text[start:index].strip())
            start = index + 1
    parts.append(text[start:].strip())
    return [part for part in parts if part]


def parameter_names(signatures):
    names = set()
    for signature in signatures:
        names.update(re.findall(r"\b([A-Za-z_]\w*)(?=:)", signature["definition"]))
    return names


def is_variadic(signatures):
    return any("**" in signature["definition"] for signature in signatures)


def check_arguments(label, arguments, signatures):
    """Every lower-case argument name in a card call must be a real parameter."""
    allowed = parameter_names(signatures)
    if is_variadic(signatures):
        return
    for token in split_arguments(arguments):
        if token.startswith("...") or not re.match(r"[a-z_]", token):
            continue  # an alias value (R, F, V), a literal, or *unpacking
        match = re.match(r"[a-z_][a-z0-9_]*", token)
        candidate, rest = match.group(0), token[match.end():]
        if rest and rest[0] not in ":=?|":
            continue  # an expression such as b.bbox
        alternatives = [candidate]
        if rest.startswith("|"):
            alternatives += [
                item for item in rest[1:].split("|") if re.fullmatch(r"[a-z_][a-z0-9_]*", item)
            ]
        for name in alternatives:
            assert name in allowed, f"LLM card uses unknown argument: {label}.{name}"


def member(exports, class_name, name):
    members = {item["name"]: item for item in exports[class_name]["members"]}
    assert name in members, f"unknown LLM-card member: {class_name}.{name}"
    return members[name]


def returned_class(definitions_or_doc, classes):
    """The public class named after the last ``->`` of a signature, if any."""
    match = re.search(r"->\s*(?:lsdtt3\.)?([A-Za-z_]\w*)\s*$", definitions_or_doc.strip())
    if match and match.group(1) in classes:
        return match.group(1)
    return None


def check_cards(cards, exports):
    classes = {name for name, entry in exports.items() if entry["kind"] == "class"}
    for names, arguments in FUNCTION_CALL.findall(cards):
        for name in names.split("|"):
            assert name in exports, f"unknown LLM-card function: {name}"
            assert exports[name]["kind"] == "function", f"LLM card calls {name} as a function"
            check_arguments(name, arguments, exports[name]["signatures"])
    for name, arguments in CONSTRUCTOR_CALL.findall(cards):
        if name not in classes:
            continue
        init = member(exports, name, "__init__")
        assert init["signatures"], f"LLM card constructs {name}, which has no constructor"
        check_arguments(name, arguments, init["signatures"])
    for alias, name, arguments in MEMBER_CALL.findall(cards):
        class_name = ALIASES.get(alias, alias if alias in classes else None)
        if class_name is None:
            continue
        item = member(exports, class_name, name)
        assert item["kind"] != "property", f"LLM card calls property as method: {alias}.{name}"
        check_arguments(f"{alias}.{name}", arguments, item["signatures"])
    for alias, name in MEMBER_ACCESS.findall(cards):
        class_name = ALIASES.get(alias, alias if alias in classes else None)
        if class_name is not None:
            member(exports, class_name, name)


class WorkflowChecker(ast.NodeVisitor):
    """Check lsdtt3 calls in example code, following the types they return."""

    def __init__(self, exports):
        self.exports = exports
        self.classes = {name for name, entry in exports.items() if entry["kind"] == "class"}
        self.types = {}

    def check_keywords(self, label, node, signatures):
        if is_variadic(signatures):
            return
        allowed = parameter_names(signatures)
        for keyword in node.keywords:
            if keyword.arg is not None:
                assert keyword.arg in allowed, f"workflow uses unknown keyword: {label}.{keyword.arg}"

    def type_of(self, node):
        if isinstance(node, ast.Name):
            return self.types.get(node.id)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
            owner = node.func.value
            if isinstance(owner, ast.Name) and owner.id == "lsdtt3":
                entry = self.exports.get(node.func.attr)
                assert entry is not None, f"workflow calls unknown export lsdtt3.{node.func.attr}"
                if entry["kind"] == "class":
                    self.check_keywords(node.func.attr, node, member(self.exports, node.func.attr, "__init__")["signatures"])
                    return node.func.attr
                self.check_keywords(node.func.attr, node, entry["signatures"])
                return self.returned(entry["signatures"])
            class_name = self.type_of(owner)
            if class_name is None:
                return None
            item = member(self.exports, class_name, node.func.attr)
            assert item["kind"] != "property", f"workflow calls property {class_name}.{node.func.attr}"
            self.check_keywords(f"{class_name}.{node.func.attr}", node, item["signatures"])
            return self.returned(item["signatures"])
        if isinstance(node, ast.Attribute):
            class_name = self.type_of(node.value)
            if class_name is None:
                return None
            item = member(self.exports, class_name, node.attr)
            if item["kind"] == "property":
                return returned_class(item["doc"].splitlines()[0] if item["doc"] else "", self.classes)
        return None

    def returned(self, signatures):
        for signature in signatures:
            found = returned_class(signature["definition"], self.classes)
            if found:
                return found
        return None

    def visit_Assign(self, node):
        self.visit(node.value)
        if len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            self.types[node.targets[0].id] = self.type_of(node.value)

    def visit_Call(self, node):
        self.type_of(node)
        self.generic_visit(node)

    def visit_Attribute(self, node):
        self.type_of(node)
        self.generic_visit(node)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads((root / "api" / "public_api.json").read_text(encoding="utf-8"))
    exports = {entry["name"]: entry for entry in manifest["exports"]}

    # Card notation everywhere except the code block, which is Python.
    cards = "\n".join([_INDEX] + list(_TOPICS.values()))
    cards = re.sub(r"```python\n.*?```", "", cards, flags=re.DOTALL)
    check_cards(cards, exports)

    tree = ast.parse(WORKFLOW_CODE)
    compile(tree, "workflow", "exec")
    WorkflowChecker(exports).visit(tree)
    workflow = lsdtt3.helper_for_llm("workflow")
    assert WORKFLOW_CODE in workflow, "the workflow card must carry WORKFLOW_CODE verbatim"

    hydraulic = lsdtt3.helper_for_llm("hydraulic")
    assert "fill_hydraulic_surface(surface:R,node_state:R" in hydraulic
    assert "boundary_conditions='nnnn'" in hydraulic
    assert "fill_hydraulic_surface(z:R,h:R" not in hydraulic
    for function in (
        "run_graphflood", "run_graphflood_optimised",
        "run_graphflood_optimised_v2", "run_graphflood_multigrid",
        "run_graphflood_induced",
    ):
        assert f"{function}(bed:R,node_state:R" in hydraulic

    landscape = lsdtt3.helper_for_llm("landscape")
    assert "coupling_regime='automatic'|'coupled'|'channels_as_boundaries'" in landscape
    assert "hillslope_model='none'|'linear'|'roering'|'roering_implicit'|'cidre'|'adi'|'cellular'" in landscape
    assert "L.step(compute_erosion_rate=False)" in landscape
    assert "L.run(n_steps,compute_erosion_rate=False)" in landscape
    assert "dt_years" not in landscape

    cosmo = lsdtt3.helper_for_llm("cosmo")
    assert "set_erosion_history(time_yr_bp,erosion_rate_m_yr,removal_metres)" in cosmo
    assert "add_erosion_rate" not in cosmo
    assert "add_burial" not in cosmo
    flow = lsdtt3.helper_for_llm("flow")
    assert "method='downstream'|'expanding_square'|'largest_area'" in flow
    assert "nearest_channel" not in flow
    assert "threshold_contributing_pixels?,threshold_drainage_area_m2?" in flow
    assert "F.snap_points(x?,y?," in flow
    assert 'method="downstream"' in WORKFLOW_CODE
    assert "threshold_contributing_pixels=threshold" in WORKFLOW_CODE
    assert "F.snap_points(" in WORKFLOW_CODE and "Transformer" not in WORKFLOW_CODE
    terrain = lsdtt3.helper_for_llm("terrain")
    for text in ("M.xy_to_row_col(x,y)", "M.row_col_to_xy(row,col)",
                 "M.lonlat_to_row_col(longitude,latitude,crs='EPSG:4326')"):
        assert text in terrain, text
    assert "topographic_wetness_index(drainage_area:R,slope:R)" in flow
    assert "F.height_above_channel(dem:R,channels:R)" in flow
    basin = lsdtt3.helper_for_llm("basin")
    assert "B.labels->R" in basin and "B.polygons->V" in basin
    assert "B.labels()" not in basin and "B.polygons()" not in basin
    assert "outline_ridgelines(dem:R,split_method='longest_channel')->V" in basin
    assert "list[Ridgeline]" not in basin
    channel = lsdtt3.helper_for_llm("channel")
    assert "extract_channel_heads(dem:R,flow:F," in channel
    vector = lsdtt3.helper_for_llm("vector")
    assert "extract_valley_centrelines(dem:R,scope:S|channels," in vector
    chi = lsdtt3.helper_for_llm("chi")
    assert "N.chi(outlets?,m_over_n=0.45,reference_area=1.0,denominator?)" in chi
    assert "fit.invert_ksn(" in chi and "chi_ksn_inversion" not in chi
    chi_symbol = lsdtt3.helper_for_llm("symbol:ChiCollection")
    assert "ChannelNetwork.chi()" in chi_symbol
    assert "divide_migration(flow:F,dem:R,outlets," in basin
    assert "route_sediment(network:N,dem:R," in landscape
    assert "crn_uniform_raster(R,value)->R" in cosmo
    assert "crn_uniform_raster(R,value=1)" not in cosmo
    fetch = lsdtt3.helper_for_llm("fetch")
    assert "fetch_dem(west,south,east,north,source='cop30_aws',target_epsg=0,grid_spacing=30.0" in fetch
    assert "bounds_from_vector(path,layer='')" in fetch and "fetch_dem(*b.bbox)" in fetch
    assert "utm_epsg_for_point(longitude,latitude)" in fetch
    assert "dataset='cop30'" not in fetch and "bounds|V" not in fetch
    stats = lsdtt3.helper_for_llm("stats")
    assert "add|subtract|multiply|divide(a:R,b:R)->R" in stats and "scale(raster:R,factor)" in stats
    assert "R|scalar" not in stats
    print("LLM API cards match current high-risk contracts")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
