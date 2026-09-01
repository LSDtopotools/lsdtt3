"""Focused tests for in-memory synthetic workflow bindings.

Authors
-------
B. Gailleton
"""

import numpy as np

import lsdtt3


def meta():
    """Build a small projected test grid.

    Parameters
    ----------
    None

    Returns
    -------
    lsdtt3.RasterMeta
        Synthetic grid metadata.

    Authors
    -------
    B. Gailleton
    """
    return lsdtt3.RasterMeta.from_shape(24, 32, cell_size=10.0, epsg_code=32631)


def test_noise_is_reproducible_and_has_requested_shape():
    """Check seeded stochastic surfaces are deterministic.

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
    first = lsdtt3.synthetic_noise(meta(), kind="red", seed=9, stddev=4.0)
    second = lsdtt3.synthetic_noise(meta(), kind="red", seed=9, stddev=4.0)
    assert first.shape == (24, 32)
    np.testing.assert_array_equal(first.to_numpy(), second.to_numpy())


def test_synthetic_surfaces_supply_matching_truth_masks():
    """Check range and hillslope products remain co-registered.

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
    mountain = lsdtt3.synthetic_mountain_range(meta(), geometry="arc", seed=4)
    strip = lsdtt3.synthetic_hillslope_strip(meta(), init_profile="linear_ramp")
    assert mountain["elevation"].shape == mountain["mountain_mask"].shape
    assert strip["elevation"].shape == strip["channel_mask"].shape
    assert np.any(strip["channel_mask"].to_numpy() == 1.0)


def test_ocn_returns_dem_area_and_receiver_grid():
    """Check a short optimal-channel-network run has coherent outputs.

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
    result = lsdtt3.synthetic_network(meta(), n_iterations=200, seed=5)
    assert result.elevation.shape == (24, 32)
    assert result.drainage_area_cells.shape == (24, 32)
    assert result.receivers.shape == (24, 32)
    assert result.final_energy > 0.0
