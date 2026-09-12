#!/usr/bin/env python3

"""Validate LLM API cards against current public contracts."""

import json
from pathlib import Path
import re

import lsdtt3
from lsdtt3._llm import _TOPICS


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads((root / "api" / "public_api.json").read_text(encoding="utf-8"))
    exports = {entry["name"]: entry for entry in manifest["exports"]}
    cards = "\n".join(_TOPICS.values())
    mentioned_functions = set(
        re.findall(r"(?:^|[; ])([a-z][a-z0-9_]*)\(", cards, re.MULTILINE)
    )
    unknown_functions = mentioned_functions - set(exports)
    assert not unknown_functions, f"unknown LLM-card functions: {sorted(unknown_functions)}"

    aliases = {
        "R": "Raster", "M": "RasterMeta", "F": "FlowInfo",
        "G": "FlowGraph", "S": "ChannelScope", "N": "ChannelNetwork",
        "B": "BasinCollection", "C": "ChiCollection", "L": "LandscapeModel",
        "Q": "CRNColumn",
    }
    for alias, method, arguments in re.findall(
        r"\b([A-Z])\.([a-z_][a-z0-9_]*)\(([^)]*)\)", cards
    ):
        if alias not in aliases:
            continue
        members = {item["name"]: item for item in exports[aliases[alias]]["members"]}
        assert method in members, f"unknown LLM-card method: {alias}.{method}"
        assert members[method]["kind"] != "property", (
            f"LLM card calls property as method: {alias}.{method}"
        )
        allowed = set()
        for signature in members[method]["signatures"]:
            allowed.update(re.findall(r"\b([A-Za-z_]\w*)(?=:)", signature["definition"]))
        for token in arguments.split(","):
            token = token.strip()
            if not token or token == "..." or token[:1].isupper():
                continue
            if "=" not in token and not token.endswith("?"):
                continue
            candidate = token.split("=", 1)[0].rstrip("?")
            if re.fullmatch(r"[a-z_][a-z0-9_]*", candidate):
                assert candidate in allowed, (
                    f"LLM card uses unknown argument: {alias}.{method}.{candidate}"
                )

    hydraulic = lsdtt3.helper_for_llm("hydraulic")
    assert "fill_hydraulic_surface(surface:R,node_state:R" in hydraulic
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
    assert "set_erosion_history(time_yr_bp,erosion_rate_m_yr,removal_m)" in cosmo
    assert "add_erosion_rate" not in cosmo
    assert "add_burial" not in cosmo
    flow = lsdtt3.helper_for_llm("flow")
    assert "topographic_wetness_index(drainage_area:R,slope:R)" in flow
    basin = lsdtt3.helper_for_llm("basin")
    assert "B.labels->R" in basin and "B.polygons->V" in basin
    assert "B.labels()" not in basin and "B.polygons()" not in basin
    channel = lsdtt3.helper_for_llm("channel")
    assert "extract_channel_heads(R,F," in channel
    vector = lsdtt3.helper_for_llm("vector")
    assert "extract_valley_centrelines(R,S|ordered_channels," in vector
    chi = lsdtt3.helper_for_llm("chi")
    assert "N.chi(outlets?,m_over_n=.45,reference_area=1,denominator?)" in chi
    chi_symbol = lsdtt3.helper_for_llm("symbol:ChiCollection")
    assert "N.chi(outlets?,m_over_n=.45,reference_area=1,denominator?)" in chi_symbol
    assert "divide_migration(F,R,outlets," in basin
    assert "route_sediment(N,R," in landscape
    assert "crn_uniform_raster(R,value)->R" in cosmo
    assert "crn_uniform_raster(R,value=1)" not in cosmo
    print("LLM API cards match current high-risk contracts")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
