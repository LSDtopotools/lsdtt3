"""Focused regression tests for tidal raster workflows.

Authors
-------
B. Gailleton
"""

import numpy as np

import lsdtt3


def trough_dem():
    """Create a shallow channel that reaches both grid edges.

    Parameters
    ----------
    None

    Returns
    -------
    lsdtt3.Raster
        Small projected raster with an east-west trough.

    Authors
    -------
    B. Gailleton
    """
    meta = lsdtt3.RasterMeta.from_shape(31, 41, cell_size=2.0, epsg_code=32631)
    values = np.full((31, 41), 2.0, dtype=np.float32)
    values[14:17, :] = 0.0
    return lsdtt3.Raster.from_numpy(values, np.ones_like(values, dtype=bool), meta)


def test_flood_mask_extracts_channel_skeleton_and_upl():
    """Check the simple channel mode creates all expected raster products.

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
    result = lsdtt3.find_tidal_channels(
        trough_dem(), use_flood_mask=True, flood_elevation_metres=0.5,
        extract_skeleton=True, compute_upl=True,
    )
    assert result.channel_mask.shape == (31, 41)
    assert result.skeleton is not None
    assert result.upl is not None
    assert np.count_nonzero(result.channel_mask.to_numpy() == 1.0) == 3 * 41


def test_skeleton_metrics_return_station_table_and_orders():
    """Check channel skeletons can be converted into ordered stations.

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
    dem = trough_dem()
    channels = lsdtt3.find_tidal_channels(
        dem, use_flood_mask=True, flood_elevation_metres=0.5, extract_skeleton=True,
    )
    metrics = lsdtt3.analyse_tidal_skeleton(
        dem, channels.channel_mask, channels.skeleton,
        outlet_mode="edge_nodes", station_spacing_metres=4.0,
    )
    assert metrics.n_sections >= 1
    assert metrics.n_stations >= 2
    assert {"x", "y", "strahler_order", "width_metres"}.issubset(metrics.stations.columns)
    assert metrics.section_ids.shape == dem.shape
