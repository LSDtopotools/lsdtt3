"""Regression checks for public Python API documentation.

Authors
-------
B. Gailleton
"""

import lsdtt3


def test_exported_functions_have_numpy_style_docstrings_and_author():
    """Require NumPy headings and B.G. attribution on the public API.

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
        "Authors\n    -------\n    B. Gailleton",
    )
    class_sections = (
        "Parameters\n    ----------",
        "Attributes\n    ----------",
        "Authors\n    -------\n    B. Gailleton",
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
