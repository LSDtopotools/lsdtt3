"""Python bindings for LSDTopoTools3."""

import os
from pathlib import Path


def _configure_bundled_geospatial_data() -> None:
    """Point GDAL and PROJ at data shipped in a standalone wheel."""
    package_dir = Path(__file__).resolve().parent
    gdal_data = package_dir / "gdal_data"
    proj_data = package_dir / "proj_data"

    if gdal_data.is_dir():
        os.environ.setdefault("GDAL_DATA", str(gdal_data))
    if proj_data.is_dir():
        os.environ.setdefault("PROJ_DATA", str(proj_data))
        # Older GDAL/PROJ combinations use the legacy spelling.
        os.environ.setdefault("PROJ_LIB", str(proj_data))


_configure_bundled_geospatial_data()

from ._core import (
    Raster,
    RasterMeta,
    add,
    add_uniform_noise,
    aspect,
    breach_depressions,
    breach_least_cost,
    carve_depth,
    clip_to_maximum_rectangle,
    clip_to_polygons,
    divide,
    fill_depth,
    fill_pits,
    fill_pits_boundary_aware,
    fill_pits_with_outlets,
    fill_pits_with_sinks,
    fill_interior_nodata,
    hillshade,
    invert_binary_mask,
    least_cost_carve_depth,
    make_boundary_conditions,
    make_outlet_mask,
    mask_by_binary,
    mask_from_points,
    mask_to_nodata,
    multidirectional_hillshade,
    mosaic,
    multiply,
    normalize_unit_interval,
    plan_curvature,
    rasterize_polygons,
    reproject,
    reproject_to_grid,
    scale,
    set_nodata_rectangle,
    sample_points,
    slope,
    subtract,
    threshold_above,
    threshold_below,
    total_curvature,
)

__all__ = [
    "Raster",
    "RasterMeta",
    "add",
    "add_uniform_noise",
    "aspect",
    "breach_depressions",
    "breach_least_cost",
    "carve_depth",
    "clip_to_maximum_rectangle",
    "clip_to_polygons",
    "divide",
    "fill_depth",
    "fill_pits",
    "fill_pits_boundary_aware",
    "fill_pits_with_outlets",
    "fill_pits_with_sinks",
    "fill_interior_nodata",
    "hillshade",
    "invert_binary_mask",
    "least_cost_carve_depth",
    "make_boundary_conditions",
    "make_outlet_mask",
    "mask_by_binary",
    "mask_from_points",
    "mask_to_nodata",
    "multidirectional_hillshade",
    "mosaic",
    "multiply",
    "normalize_unit_interval",
    "plan_curvature",
    "rasterize_polygons",
    "reproject",
    "reproject_to_grid",
    "scale",
    "set_nodata_rectangle",
    "sample_points",
    "slope",
    "subtract",
    "threshold_above",
    "threshold_below",
    "total_curvature",
]
