"""Regression checks for public Python API documentation.

Authors
-------
B. Gailleton
"""

import lsdtt3


def test_exported_functions_have_numpy_style_docstrings_and_author():
    """Require NumPy headings and B.G. attribution on public functions.

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
    required_sections = ("Parameters\n    ----------", "Returns\n    -------", "Authors\n    -------\n    B. Gailleton")

    for name in lsdtt3.__all__:
        value = getattr(lsdtt3, name)
        if not callable(value) or isinstance(value, type):
            continue
        docstring = value.__doc__ or ""
        if not all(section in docstring for section in required_sections):
            missing.append(name)

    assert not missing, f"Missing NumPy-style public docstrings: {missing}"
