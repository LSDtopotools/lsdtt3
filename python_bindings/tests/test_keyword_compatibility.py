import inspect
import warnings

import pytest

from lsdtt3._compat import MISSING, resolve_deprecated_keyword
import lsdtt3._core as native


def native_function(prefix, *, value=5):
    return prefix, value


def compatibility_function(prefix, *, canonical=MISSING, legacy=MISSING):
    value = resolve_deprecated_keyword(
        "canonical", canonical, "legacy", legacy, stacklevel=3
    )
    if value is MISSING:
        return native_function(prefix)
    return native_function(prefix, value=value)


compatibility_function.__wrapped__ = native_function
compatibility_function.__signature__ = inspect.signature(
    compatibility_function, follow_wrapped=False
)


def test_function_keyword_compatibility_prototype_preserves_signature_and_value():
    parameters = inspect.signature(compatibility_function).parameters
    assert list(parameters) == ["prefix", "canonical", "legacy"]
    assert parameters["prefix"].kind is inspect.Parameter.POSITIONAL_OR_KEYWORD
    assert parameters["canonical"].kind is inspect.Parameter.KEYWORD_ONLY
    assert parameters["legacy"].kind is inspect.Parameter.KEYWORD_ONLY
    assert parameters["canonical"].default is MISSING
    assert parameters["legacy"].default is MISSING
    assert inspect.unwrap(compatibility_function) is native_function
    assert compatibility_function("p") == ("p", 5)
    assert compatibility_function("p", canonical=7) == ("p", 7)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        warning_line = inspect.currentframe().f_lineno + 1
        assert compatibility_function("p", legacy=9) == ("p", 9)
    assert len(caught) == 1
    assert caught[0].category is DeprecationWarning
    assert caught[0].filename == __file__
    assert caught[0].lineno == warning_line
    assert "use 'canonical'" in str(caught[0].message)


def test_function_keyword_compatibility_prototype_rejects_conflict():
    with pytest.raises(TypeError, match="cannot both be supplied"):
        compatibility_function("p", canonical=7, legacy=9)


@pytest.mark.skipif(
    not hasattr(native, "_KeywordCompatibilityProbe")
    and not bool(__import__("os").environ.get("LSDTT_REQUIRE_COMPATIBILITY_PROBE")),
    reason="private nanobind compatibility probe is disabled",
)
def test_native_constructor_keyword_compatibility_prototype():
    probe_type = native._KeywordCompatibilityProbe
    signatures = probe_type.__init__.__nb_signature__
    assert len(signatures) == 2
    assert "canonical" in signatures[0][0]
    assert "legacy" in signatures[1][0]

    canonical = probe_type(canonical=7)
    assert type(canonical) is probe_type
    assert isinstance(canonical, probe_type)
    assert canonical.value == 7
    positional = probe_type(8)
    assert positional.value == 8

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        warning_line = inspect.currentframe().f_lineno + 1
        legacy = probe_type(legacy=9)
    assert type(legacy) is probe_type
    assert legacy.value == 9
    assert len(caught) == 1
    assert caught[0].category is DeprecationWarning
    assert caught[0].filename == __file__
    assert caught[0].lineno == warning_line

    with pytest.raises(TypeError):
        probe_type(canonical=7, legacy=9)

    class Child(probe_type):
        pass

    child = Child(canonical=11)
    assert isinstance(child, probe_type)
    assert child.value == 11
