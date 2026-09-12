"""Tests for compact LLM API context retrieval.

Authors
-------
B. Gailleton
"""

import pytest
import importlib.util
from pathlib import Path

import lsdtt3


def test_index_is_small_and_exposes_retrieval_protocol():
    """The starter card stays small enough for routine prompt injection."""
    card = lsdtt3.helper_for_llm()
    assert card.startswith("LSDTT3/API1")
    assert "?: helper_for_llm(T); symbol:<exact>; search:<words>; all" in card
    assert len(card.split()) < 100


def test_topic_card_keeps_exact_code_identifiers():
    """A topic card has the API names required for generated Python code."""
    card = lsdtt3.helper_for_llm("flow")
    assert "FlowInfo" in card
    assert "F.channel_network(S)->N" in card
    assert "F.basins(outlets?)->B" in card


def test_symbol_and_search_retrieval_use_public_api_only():
    """Direct lookup avoids injecting an unrelated full documentation page."""
    symbol = lsdtt3.helper_for_llm("symbol:fill_pits")
    assert symbol.startswith("@fill_pits\nfill_pits(")
    assert "FlowInfo(dem,boundary_conditions='oooo'" in lsdtt3.helper_for_llm("symbol:FlowInfo")
    search = lsdtt3.helper_for_llm("search:flow")
    assert "FlowGraph" in search
    assert "FlowInfo" in search


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
