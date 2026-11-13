"""FlowInfo.channel_point_attributes: lsdtt-channel-extraction run_channel_point_attributes.

The DEM is the CLI smoke test's tiny_dem.asc (5 x 5, 1 m, z = 100 - 2*row - 5*col):
each row drains east, then down the east edge, and a cell has col + 1
contributing pixels. With threshold 3 the channel cells are columns 2-4. The
same point as the CLI test, at (0.5, 2.4) in column 0 of row 2, must get the
same values: flat 12, sqrt(4.01) m, 3 pixels, flow distance 4 m, elevation 86,
slope 5, smoothed elevation 86, window 2 m. A point off the grid keeps its
position with empty values.
"""

import math
import warnings

import geopandas as gpd
import numpy as np
import pytest
from pyproj import CRS
from shapely.geometry import LineString, Point

import lsdtt3

pytestmark = [pytest.mark.unit]

ADDED = ["original_x", "original_y", "snap_distance_m", "snap_status", "snapped_flat",
         "contributing_pixels", "drainage_area_m2", "flow_distance_m", "elevation_m",
         "smoothed_elevation_m", "slope_m_per_m", "slope_valid", "slope_window_length_m"]
WKT = CRS.from_epsg(32630).to_wkt()


@pytest.fixture(scope="module")
def grid():
    rows, cols = np.indices((5, 5))
    z = (100.0 - 2.0 * rows - 5.0 * cols).astype(np.float32)
    meta = lsdtt3.RasterMeta.from_shape(5, 5, x_min=0.0, y_min=0.0, cell_size=1.0, crs_wkt=WKT,
                                        epsg_code=32630)
    dem = lsdtt3.Raster.from_numpy(z, np.ones((5, 5), dtype=bool), meta)
    filled = lsdtt3.fill_pits(dem)
    return filled, lsdtt3.FlowInfo(filled, boundary_conditions="oooo")


def _sites():
    return gpd.GeoDataFrame({"site_id": ["01", "02"], "d50_mm": [39.6, np.nan]},
                            geometry=[Point(0.5, 2.4), Point(-50.0, 2.4)], crs="EPSG:32630")


@pytest.mark.exact_solution
@pytest.mark.conformance
def test_same_values_as_the_cli(grid):
    """channel_point_attributes gives the same columns and values as the CLI: the first
    site snaps to channel cell 12 with its exact snap distance, area, flow distance,
    elevation and slope, and keeps its own attributes; the second, off the grid, is
    input_out_of_bounds, keeps its position and gets empty values.
    """
    filled, flow = grid
    out = flow.channel_point_attributes(filled, _sites(), threshold_contributing_pixels=3,
                                        slope_window_half_length_metres=2.0)
    assert list(out.columns) == ["site_id", "d50_mm", "geometry"] + ADDED
    a = out.iloc[0]
    assert a["site_id"] == "01" and a["d50_mm"] == 39.6
    assert a["snap_status"] == "success" and a["snapped_flat"] == 12
    assert a["snap_distance_m"] == pytest.approx(math.sqrt(4.01), abs=1e-12)
    assert a["contributing_pixels"] == 3 and a["drainage_area_m2"] == 3.0
    assert a["flow_distance_m"] == 4.0 and a["elevation_m"] == 86.0
    assert a["slope_m_per_m"] == pytest.approx(5.0) and a["smoothed_elevation_m"] == pytest.approx(86.0)
    assert a["slope_valid"] == 1 and a["slope_window_length_m"] == 2.0
    assert (a["original_x"], a["original_y"]) == (0.5, 2.4)
    assert (a.geometry.x, a.geometry.y) == (2.5, 2.5)          # the snapped cell centre
    b = out.iloc[1]
    assert b["snap_status"] == "input_out_of_bounds" and b["slope_valid"] == 0
    assert b["snapped_flat"] is pd_na() and math.isnan(b["slope_m_per_m"])
    assert (b.geometry.x, b.geometry.y) == (-50.0, 2.4)        # keeps its position
    assert out.crs.equals(CRS.from_wkt(WKT))


def pd_na():
    import pandas as pd
    return pd.NA


@pytest.mark.metamorphic
def test_points_in_another_crs_are_reprojected(grid):
    """Sites given in WGS84 are reprojected and snap to the same cell, with their original
    position reported in the DEM's CRS.
    """
    filled, flow = grid
    ll = _sites().to_crs("EPSG:4326")
    out = flow.channel_point_attributes(filled, ll, threshold_contributing_pixels=3,
                                        slope_window_half_length_metres=2.0)
    assert out.iloc[0]["snapped_flat"] == 12
    assert out.iloc[0]["original_x"] == pytest.approx(0.5, abs=1e-6)


@pytest.mark.exact_solution
def test_clashing_column_renamed_with_warning(grid):
    """A site column named like an added field (slope_m_per_m) is kept unchanged; the added
    field becomes slope_m_per_m_lsdtt, with a warning.
    """
    filled, flow = grid
    sites = _sites().assign(slope_m_per_m=[0.1, 0.2])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        out = flow.channel_point_attributes(filled, sites, threshold_contributing_pixels=3,
                                            slope_window_half_length_metres=2.0)
    assert any("slope_m_per_m_lsdtt" in str(w.message) for w in caught)
    assert out.iloc[0]["slope_m_per_m"] == 0.1                  # the user's column, untouched
    assert out.iloc[0]["slope_m_per_m_lsdtt"] == pytest.approx(5.0)


@pytest.mark.negative
def test_refusals(grid):
    """channel_point_attributes refuses sites with no CRS, non-point geometry, a slope
    window shorter than two cells, and input that is not a GeoDataFrame.
    """
    filled, flow = grid
    with pytest.raises(ValueError, match="CRS"):
        flow.channel_point_attributes(filled, _sites().set_crs(None, allow_override=True),
                                      threshold_contributing_pixels=3, slope_window_half_length_metres=2.0)
    lines = gpd.GeoDataFrame(geometry=[LineString([(0, 0), (1, 1)])], crs="EPSG:32630")
    with pytest.raises(ValueError, match="Point"):
        flow.channel_point_attributes(filled, lines, threshold_contributing_pixels=3,
                                      slope_window_half_length_metres=2.0)
    with pytest.raises(ValueError, match="two cell sizes"):
        flow.channel_point_attributes(filled, _sites(), threshold_contributing_pixels=3,
                                      slope_window_half_length_metres=1.0)
    with pytest.raises(ValueError, match="GeoDataFrame"):
        flow.channel_point_attributes(filled, [(0.5, 2.4)], threshold_contributing_pixels=3)
