"""Regression checks for public Python API documentation.

Authors
-------
B. Gailleton
"""

import lsdtt3


def test_exported_functions_have_numpy_style_docstrings_and_author():
    """Require NumPy headings and an Authors section on the public API.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Authors
    -------
    B. Gailleton
    """
    missing = []
    function_sections = (
        "Parameters\n    ----------",
        "Returns\n    -------",
        "Authors\n    -------\n    ",
        "B. Gailleton",
    )
    class_sections = (
        "Parameters\n    ----------",
        "Attributes\n    ----------",
        "Authors\n    -------\n    ",
        "B. Gailleton",
    )

    for name in lsdtt3.__all__:
        value = getattr(lsdtt3, name)
        if not callable(value):
            continue
        docstring = value.__doc__ or ""
        required_sections = class_sections if isinstance(value, type) else function_sections
        if not all(section in docstring for section in required_sections):
            missing.append(name)

    assert not missing, f"Missing NumPy-style public docstrings: {missing}"


def test_additional_authors_follow_the_author_list():
    """Each export lists the authors :mod:`lsdtt3._authors` gives it, in order.

    Every callable export must have exactly one ``Authors`` section (a
    binding's own section is folded into the generated one), and its authors
    must be the ones :mod:`lsdtt3._authors` gives for that name.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Authors
    -------
    S. M. Mudd
    """
    from lsdtt3._authors import AUTHORS, CHANGED_BY_MUDD, WRITTEN_BY_MUDD, authors_for

    assert set(AUTHORS) <= set(lsdtt3.__all__), (
        f"names in lsdtt3._authors that are not exports: {set(AUTHORS) - set(lsdtt3.__all__)}")
    assert not set(WRITTEN_BY_MUDD) & set(CHANGED_BY_MUDD), "an export is in both lists"

    wrong = []
    for name in lsdtt3.__all__:
        value = getattr(lsdtt3, name)
        if not callable(value):
            continue
        docstring = value.__doc__ or ""
        expected = "Authors\n    -------\n    " + "\n    ".join(authors_for(name))
        if docstring.count("Authors\n") != 1 or expected not in docstring:
            wrong.append(name)
    assert not wrong, f"Exports whose Authors section is missing, doubled or wrong: {wrong}"

    assert "Authors\n    -------\n    S. M. Mudd\n    B. Gailleton" in lsdtt3.coregister_dem.__doc__
    assert "Authors\n    -------\n    B. Gailleton\n    S. M. Mudd" in lsdtt3.LandscapeModel.__doc__
    assert "S. M. Mudd" not in lsdtt3.FlowInfo.__doc__
