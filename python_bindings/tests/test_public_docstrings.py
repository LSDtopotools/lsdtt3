"""Regression checks for public Python API documentation.

These tests check content, not only headings: the generated ``Parameters``
section of every function must name its real parameters,
``inspect.signature`` must show them, and every class must say how to get an
instance and list its members.

Authors
-------
B. Gailleton
S. M. Mudd
"""

import inspect
import re

import pytest

import lsdtt3

pytestmark = [pytest.mark.conformance, pytest.mark.unit]


# Core-workflow classes, documented in their bindings: none may carry the
# placeholder.
CORE_WORKFLOW_CLASSES = (
    "Raster", "RasterMeta", "FlowInfo", "ChannelScope", "ChannelNetwork",
    "BasinCollection", "ChiCollection", "ChiFit",
)

# Functions of the core (Colab) workflow: every parameter is described.
CORE_WORKFLOW_FUNCTIONS = (
    "fetch_dem", "fetch_dem_for_outlets", "bounds_from_vector", "utm_epsg_for_point",
    "mask_to_nodata", "fill_pits", "hillshade", "multidirectional_hillshade",
    "slope", "polyfit_metrics", "swath_profile", "junction_basins",
    "extract_channel_heads", "sample_points", "run_graphflood",
    "crn_invert_basin_erosion", "reproject", "reproject_to_grid", "mosaic",
    "clip_to_polygons", "rasterize_polygons", "add", "subtract", "multiply",
    "divide", "scale", "threshold_above", "threshold_below", "mask_by_binary",
)

# Methods of the core-workflow classes: each binding docstring has a
# Parameters section naming every argument, and a Returns section.
CORE_WORKFLOW_METHODS = {
    "Raster": ("read", "from_numpy", "write", "to_numpy", "mask_numpy"),
    "RasterMeta": ("from_shape", "xy_to_row_col", "row_col_to_xy", "lonlat_to_row_col"),
    "FlowInfo": (
        "channel_scope", "channel_network", "basins", "snap_points", "receiver_of",
        "contributing_pixels", "drainage_area", "flow_distance_to_outlet",
        "row_col", "flat_index", "upslope_pixels", "is_outlet", "sources",
        "channel_mask", "height_above_channel", "accumulate",
    ),
    "ChannelScope": ("nodes", "table"),
    "ChannelNetwork": (
        "lines_geodataframe", "points_geodataframe", "chi", "junctions", "links",
    ),
    "BasinCollection": ("mask", "pixels", "membership", "statistics", "outline_ridgelines"),
    "ChiCollection": ("nodes", "raster", "fit"),
    "ChiFit": ("disorder", "invert_ksn"),
    "LandscapeModel": ("step", "run"),
}

# Properties of the core-workflow classes that must carry a description.
CORE_WORKFLOW_PROPERTIES = {
    "Raster": ("metadata", "shape", "valid_pixel_count"),
    "FlowInfo": ("outlets", "metadata", "shape"),
    "ChannelNetwork": ("channel_mask", "n_channel_cells", "n_junctions"),
    "BasinCollection": ("table", "polygons", "labels"),
    "ChiFit": ("nodes", "segments", "summary"),
}

UNDESCRIBED = "Not yet described."

PLACEHOLDER = "No additional binding notes are available."


def _exports(kind):
    """Callable exports that are classes (``"class"``) or not (``"function"``)."""
    for name in lsdtt3.__all__:
        value = getattr(lsdtt3, name)
        if callable(value) and isinstance(value, type) == (kind == "class"):
            yield name, value


def _section(doc, heading):
    """The body of one NumPy section of a cleaned docstring, or None."""
    match = re.search(
        rf"^{re.escape(heading)}\n-{{{len(heading)}}}\n(.*?)(?=^\S[^\n]*\n-{{3,}}\n|\Z)",
        doc, re.MULTILINE | re.DOTALL,
    )
    return match.group(1) if match else None


def _native_parameter_names(value):
    """Parameter names of every call form, read from the binding itself."""
    native = inspect.unwrap(value)
    signatures = getattr(native, "__nb_signature__", None)
    names = []
    if signatures is not None:
        for definition, _doc, _defaults in signatures:
            for name in re.findall(r"(?:\(|, )\**([A-Za-z_]\w*)(?=[:,)])", definition):
                if name != "self" and name not in names:
                    names.append(name)
        return names
    return [
        name for name in inspect.signature(native).parameters if name != "self"
    ]


def test_exported_functions_have_numpy_style_docstrings_and_author():
    """Require NumPy sections and an Authors section on the public API.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Authors
    -------
    B. Gailleton
    S. M. Mudd
    """
    missing = []
    for name, value in _exports("function"):
        docstring = value.__doc__ or ""
        required = ("Parameters\n    ----------", "Returns\n    -------", "Authors\n    -------\n    ")
        if not all(section in docstring for section in required):
            missing.append(name)
    for name, value in _exports("class"):
        if "Authors\n    -------\n    " not in (value.__doc__ or ""):
            missing.append(name)
    assert not missing, f"Missing NumPy-style public docstrings: {missing}"


def test_generated_docstrings_have_no_generic_envelope():
    """The old envelope's boilerplate is gone from every export."""
    boilerplate = (
        "described by the generated binding signature",
        "Constructor arguments",
        "See the generated binding signature",
    )
    stale = [
        name for name in lsdtt3.__all__
        if any(text in (getattr(lsdtt3, name).__doc__ or "") for text in boilerplate)
    ]
    assert not stale, f"Exports with generic generated sections: {stale}"
    for name, value in _exports("function"):
        doc = inspect.getdoc(value)
        assert doc.count("\nParameters\n----------\n") <= 1, name
        assert doc.count("\nReturns\n-------\n") <= 1, name
        assert not doc.startswith(f"{name}("), f"{name}: summary is the raw signature"


def test_every_parameter_is_named_in_the_parameters_section():
    """The Parameters section names every real parameter of every function."""
    undocumented = []
    for name, value in _exports("function"):
        section = _section(inspect.getdoc(value), "Parameters") or ""
        for parameter in _native_parameter_names(value):
            if not re.search(rf"(?<![\w]){re.escape(parameter)}\b", section):
                undocumented.append(f"{name}.{parameter}")
    assert not undocumented, f"Parameters missing from generated sections: {undocumented}"


def test_generated_parameters_carry_types_and_defaults():
    """Generated docstrings give each parameter its type and default (fill_pits), keep and
    complete a binding's own Parameters section without duplicating it (sample_points),
    and mark keyword-only parameters.
    """
    doc = inspect.getdoc(lsdtt3.fill_pits)
    parameters = _section(doc, "Parameters")
    assert "raster : Raster\n    Elevation raster (m)." in parameters
    assert "min_slope : float, default=0.0001" in parameters
    assert 'priority_flood_backend : {"binary", "radix"}, default="binary"' in parameters
    assert _section(doc, "Returns").startswith("Raster\n")
    assert doc.splitlines()[0] == "Fill depressions while imposing an optional minimum drainage slope."
    # A binding's own Parameters section is kept and completed, not duplicated.
    sample = inspect.getdoc(lsdtt3.sample_points)
    assert sample.count("Parameters\n----------") == 1
    assert "Point geometries with a CRS exactly matching raster." in sample
    assert re.search(r"^raster : Raster$", _section(sample, "Parameters"), re.MULTILINE)
    # Keyword-only parameters are marked.
    assert "valley_mask : Raster, keyword-only" in inspect.getdoc(lsdtt3.extract_valley_centrelines)


def test_overloaded_functions_list_every_call_form():
    """The docstring of an overloaded function (classify_surface_shape) lists both of its
    call forms.
    """
    doc = inspect.getdoc(lsdtt3.classify_surface_shape)
    assert "This function has 2 call forms:" in doc
    assert "1. ``classify_surface_shape(raster: Raster" in doc
    assert "2. ``classify_surface_shape(maximum_curvature: Raster, minimum_curvature: Raster" in doc


def test_inspect_signature_shows_the_real_parameters():
    """inspect.signature (and so notebook tooltips) shows real parameters."""
    variadic = []
    for name, value in _exports("function"):
        signature = inspect.signature(value)
        if [*signature.parameters] == ["args", "kwargs"]:
            variadic.append(name)
            continue
        native = [parameter for parameter in _native_parameter_names(value)]
        assert [*signature.parameters][: len(native)] == native[: len(signature.parameters)], name
    assert not variadic, f"Exports whose signature is (*args, **kwargs): {variadic}"

    fill_pits = inspect.signature(lsdtt3.fill_pits)
    assert str(fill_pits) == (
        "(raster: Raster, min_slope: float = 0.0001, n_threads: int = 0, "
        "priority_flood_backend: str = 'binary') -> Raster"
    )
    # The displayed default is shortened; the value is unchanged.
    assert fill_pits.parameters["min_slope"].default == 9.999999747378752e-05
    assert inspect.signature(lsdtt3.extract_valley_centrelines).parameters[
        "valley_mask"].kind is inspect.Parameter.KEYWORD_ONLY
    assert str(inspect.signature(lsdtt3.FlowInfo)).startswith("(dem: Raster, boundary_conditions: str = 'oooo'")


def test_classes_say_how_to_get_an_instance_and_list_members():
    """Each class gives its constructor or what returns it, and its members."""
    problems = []
    for name, cls in _exports("class"):
        doc = inspect.getdoc(cls)
        init = inspect.getattr_static(cls, "__init__")
        constructible = hasattr(init, "__nb_signature__") or inspect.isfunction(init)
        if constructible:
            if f"Construct it as\n\n    {name}(" not in doc:
                problems.append(f"{name}: no constructor signature")
        elif "not constructible directly" not in doc.lower():
            problems.append(f"{name}: does not say it is not constructible")
        for member in dir(cls):
            if member.startswith("_"):
                continue
            raw = inspect.getattr_static(cls, member)
            if isinstance(raw, property) or callable(getattr(cls, member)):
                if not re.search(rf"^(?:{re.escape(name)}\.)?{re.escape(member)}\b", doc, re.MULTILINE):
                    problems.append(f"{name}.{member}: not listed")
    assert not problems, problems

    chi_fit = inspect.getdoc(lsdtt3.ChiFit)
    assert "Returned by ``ChiCollection.fit()``; not constructible directly." in chi_fit
    for member in ("nodes", "segments", "summary"):
        assert re.search(rf"^{member} : ", _section(chi_fit, "Attributes"), re.MULTILINE)
    assert "invert_ksn(m_over_n_values: object" in _section(chi_fit, "Methods")
    raster = inspect.getdoc(lsdtt3.Raster)
    assert "create one with ``Raster.from_numpy()`` or ``Raster.read()``" in raster
    assert "Raster.read(path: str) -> Raster" in raster


def test_core_workflow_classes_have_binding_documentation():
    """The classes of the core workflow do not carry the placeholder."""
    placeholders = [
        name for name in CORE_WORKFLOW_CLASSES
        if PLACEHOLDER in (getattr(lsdtt3, name).__doc__ or "")
        or UNDESCRIBED in (getattr(lsdtt3, name).__doc__ or "")
    ]
    assert not placeholders, f"Core classes still undocumented: {placeholders}"


def test_core_workflow_functions_describe_every_parameter():
    """No parameter or return value of a core-workflow function is a placeholder."""
    undescribed = []
    for name in CORE_WORKFLOW_FUNCTIONS:
        doc = inspect.getdoc(getattr(lsdtt3, name))
        for heading in ("Parameters", "Returns"):
            if UNDESCRIBED in (_section(doc, heading) or ""):
                undescribed.append(f"{name}: {heading}")
    assert not undescribed, f"Core functions with undescribed entries: {undescribed}"


def test_core_workflow_methods_document_parameters_and_returns():
    """Each core-workflow method names every argument and its return value."""
    problems = []
    for class_name, methods in CORE_WORKFLOW_METHODS.items():
        cls = getattr(lsdtt3, class_name)
        for method in methods:
            value = getattr(cls, method)
            doc = inspect.cleandoc(inspect.unwrap(value).__doc__ or "")
            parameters = _section(doc, "Parameters")
            names = [
                parameter for parameter in _native_parameter_names(value)
                if parameter != "self"
            ]
            if names and parameters is None:
                problems.append(f"{class_name}.{method}: no Parameters section")
            for parameter in names:
                if parameters is not None and not re.search(
                    rf"(?<![\w]){re.escape(parameter)}\b", parameters
                ):
                    problems.append(f"{class_name}.{method}: {parameter} not described")
            if _section(doc, "Returns") is None:
                problems.append(f"{class_name}.{method}: no Returns section")
    assert not problems, problems


def test_core_workflow_properties_are_described():
    """Core-workflow properties have a docstring and DataFrame tables say so."""
    missing = []
    for class_name, members in CORE_WORKFLOW_PROPERTIES.items():
        attributes = _section(inspect.getdoc(getattr(lsdtt3, class_name)), "Attributes") or ""
        for member in members:
            entry = re.search(rf"^{member}\b[^\n]*\n    \S", attributes, re.MULTILINE)
            if entry is None:
                missing.append(f"{class_name}.{member}")
    assert not missing, f"Properties without a description: {missing}"
    chi_fit = _section(inspect.getdoc(lsdtt3.ChiFit), "Attributes")
    for member in ("nodes", "segments", "summary"):
        assert re.search(rf"^{member} : pandas\.DataFrame$", chi_fit, re.MULTILINE), member
    assert re.search(
        r"^table : pandas\.DataFrame$",
        _section(inspect.getdoc(lsdtt3.BasinCollection), "Attributes"), re.MULTILINE,
    )
    # The documented column lists match what the tables really contain.
    raster_meta = inspect.getdoc(lsdtt3.RasterMeta)
    assert "row 0 is the TOP (north) row" in raster_meta
    assert "mask_to_nodata(dem, 0.0)" in inspect.getdoc(lsdtt3.fetch_dem)


@pytest.mark.metamorphic
def test_generated_text_does_not_depend_on_docstring_indentation():
    """Python 3.12 keeps docstring indentation and 3.13 strips it; same result.

    A Python-defined export (``run_graphflood_multigrid``, ``helper_for_llm``)
    has its docstring indented on 3.12 and dedented on 3.13+. Both forms must
    give the same generated docstring.
    """
    def indented(raster, threshold: float = 0.5):
        return raster

    def dedented(raster, threshold: float = 0.5):
        return raster

    indented.__doc__ = (
        "Summary line.\n\n    More prose.\n\n    Parameters\n    ----------\n"
        "    raster : Raster\n        Input.\n\n    Authors\n    -------\n    X. Y.\n    "
    )
    dedented.__doc__ = (
        "Summary line.\n\nMore prose.\n\nParameters\n----------\n"
        "raster : Raster\n    Input.\n\nAuthors\n-------\nX. Y.\n"
    )
    first = lsdtt3._numpy_style_docstring(indented, "fill_pits").__doc__
    second = lsdtt3._numpy_style_docstring(dedented, "fill_pits").__doc__
    assert first == second
    doc = inspect.cleandoc(first)
    assert "threshold : float, default=0.5\n    Not yet described." in doc
    assert doc.count("Authors\n") == 1 and "X. Y." not in doc


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
    assert "Authors\n-------\nB. Gailleton\nS. M. Mudd" in inspect.cleandoc(lsdtt3.FlowInfo.__doc__)
    assert "S. M. Mudd" not in lsdtt3.Raster.__doc__
