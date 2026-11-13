import hashlib
import os

import numpy as np
import pytest
from pyproj import CRS

import lsdtt3
import lsdtt3._core as native

pytestmark = [pytest.mark.exact_solution, pytest.mark.characterisation, pytest.mark.unit]


PROBE_AVAILABLE = hasattr(native, "_extract_valley_centrelines_with_mask_probe")


@pytest.mark.skipif(
    not PROBE_AVAILABLE and not bool(os.environ.get("LSDTT_REQUIRE_COMPATIBILITY_PROBE")),
    reason="private supplied-mask probe is disabled",
)
def test_ordered_channels_explicit_mask_oracle():
    """Valley centrelines from an ordered channel and an explicit valley mask: 429 valley
    pixels and a 33-node centreline, routing succeeds, and the mask and centreline match
    their frozen sha256; a narrower mask changes the result.
    """
    metadata = lsdtt3.RasterMeta.from_shape(41, 41, cell_size=10.0)
    rows, cols = np.indices((41, 41))
    elevation = (np.abs(cols - 20) * 2.0 + rows * 0.05).astype(np.float32)
    dem = lsdtt3.Raster.from_numpy(
        elevation, np.ones(elevation.shape, dtype=bool), metadata
    )
    channel = np.column_stack((np.arange(4, 37), np.full(33, 20, dtype=int)))
    mask_values = np.zeros((41, 41), dtype=np.float32)
    mask_values[4:37, 14:27] = 1.0
    valley_mask = lsdtt3.Raster.from_numpy(
        mask_values, np.ones(mask_values.shape, dtype=bool), metadata
    )

    result = native._extract_valley_centrelines_with_mask_probe(
        dem, [channel], valley_mask, corridor_width_metres=120.0,
        valley_dilation_pixels=1, centreline_loops=1, n_threads=1,
    )
    output_mask = result.valley_mask.to_numpy()
    pixels = result.centreline_pixels(0)
    row = result.table.iloc[0]

    assert np.count_nonzero(output_mask) == 429
    assert len(pixels) == 33
    assert bool(row["routing_succeeded"])
    assert int(row["n_valley_pixels"]) == 429
    assert np.isnan(row["valley_hand_break_metres"])
    assert hashlib.sha256(output_mask.tobytes()).hexdigest() == (
        "42592fe5ed1b199b98d8f892ab846b310d40a62cd96423bddd8948cb80795ee0"
    )
    assert hashlib.sha256(pixels.tobytes()).hexdigest() == (
        "3a31ba0512de39c881c7de0a6ce076b85ec133c3eb0dd1af608844679b5992b5"
    )

    narrower_values = mask_values.copy()
    narrower_values[:, 14] = 0.0
    narrower_mask = lsdtt3.Raster.from_numpy(
        narrower_values, np.ones(narrower_values.shape, dtype=bool), metadata
    )
    perturbed = native._extract_valley_centrelines_with_mask_probe(
        dem, [channel], narrower_mask, corridor_width_metres=120.0,
        valley_dilation_pixels=1, centreline_loops=1, n_threads=1,
    )
    assert hashlib.sha256(perturbed.valley_mask.to_numpy().tobytes()).hexdigest() != (
        "42592fe5ed1b199b98d8f892ab846b310d40a62cd96423bddd8948cb80795ee0"
    )


@pytest.mark.skipif(
    not PROBE_AVAILABLE and not bool(os.environ.get("LSDTT_REQUIRE_COMPATIBILITY_PROBE")),
    reason="private supplied-mask probe is disabled",
)
def test_channel_scope_explicit_mask_oracle():
    """Valley centrelines from a channel scope with an all-valley mask: 36 valley pixels
    for each of the 6 channels, 4 centreline nodes each, and routing not succeeding for
    any.
    """
    crs = CRS.from_epsg(32630)
    metadata = lsdtt3.RasterMeta.from_shape(
        6, 6, cell_size=10.0, crs_wkt=crs.to_wkt(), epsg_code=32630
    )
    rows, cols = np.indices((6, 6), dtype=np.float32)
    dem = lsdtt3.Raster.from_numpy(
        -rows - 0.01 * cols, np.ones((6, 6), dtype=bool), metadata
    )
    flow = lsdtt3.FlowInfo(dem, boundary_conditions="nnbo", n_threads=1)
    scope = flow.channel_scope(3)
    valley_mask = lsdtt3.Raster.from_numpy(
        np.ones((6, 6), dtype=np.float32),
        np.ones((6, 6), dtype=bool), metadata,
    )

    result = native._extract_valley_centrelines_with_mask_probe(
        dem, scope, valley_mask, corridor_width_metres=100.0,
        centreline_loops=1, n_threads=1,
    )
    table = result.table

    assert np.count_nonzero(result.valley_mask.to_numpy()) == 36
    assert len(table) == 6
    assert set(table["n_valley_pixels"]) == {36}
    assert set(table["n_centreline_nodes"]) == {4}
    assert not table["routing_succeeded"].any()
    assert all(len(result.centreline_pixels(index)) == 4 for index in range(6))
