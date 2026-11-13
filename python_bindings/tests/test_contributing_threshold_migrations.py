import inspect
import warnings

import numpy as np
import pytest

import lsdtt3

pytestmark = [pytest.mark.metamorphic, pytest.mark.conformance, pytest.mark.negative, pytest.mark.unit]


def flow_fixture():
    meta = lsdtt3.RasterMeta.from_shape(6, 6, cell_size=10.0)
    rows, cols = np.indices((6, 6), dtype=np.float32)
    dem = lsdtt3.Raster.from_numpy(
        -rows - 0.01 * cols, np.ones((6, 6), dtype=bool), meta
    )
    return dem, lsdtt3.FlowInfo(dem, boundary_conditions="nnbo", n_threads=1)


def test_flow_threshold_methods_accept_legacy_keyword_without_signature_drift():
    """FlowInfo.channel_mask, channel_scope, distance_to_channel and sources expose
    threshold_contributing_pixels; the old threshold_pixels works with a
    DeprecationWarning each and gives the same results; giving both is refused.
    """
    _, flow = flow_fixture()
    for method_name in ("channel_mask", "channel_scope", "distance_to_channel", "sources"):
        parameters = inspect.signature(getattr(flow, method_name)).parameters
        assert "threshold_contributing_pixels" in parameters
        assert "threshold_pixels" not in parameters

    canonical_mask = flow.channel_mask(threshold_contributing_pixels=3)
    canonical_scope = flow.channel_scope(threshold_contributing_pixels=3)
    canonical_distance = flow.distance_to_channel(threshold_contributing_pixels=3)
    canonical_sources = flow.sources(threshold_contributing_pixels=3)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        legacy_mask = flow.channel_mask(threshold_pixels=3)
        legacy_scope = flow.channel_scope(threshold_pixels=3)
        legacy_distance = flow.distance_to_channel(threshold_pixels=3)
        legacy_sources = flow.sources(threshold_pixels=3)
    assert [warning.category for warning in caught] == [DeprecationWarning] * 4
    assert all(warning.filename == __file__ for warning in caught)
    np.testing.assert_array_equal(canonical_mask.to_numpy(), legacy_mask.to_numpy())
    assert canonical_scope.table().equals(legacy_scope.table())
    np.testing.assert_array_equal(canonical_distance.to_numpy(), legacy_distance.to_numpy())
    np.testing.assert_array_equal(canonical_sources, legacy_sources)
    with pytest.raises(TypeError, match="cannot both be supplied"):
        flow.channel_mask(threshold_contributing_pixels=3, threshold_pixels=3)


def test_channel_head_threshold_legacy_keyword_preserves_results():
    """extract_channel_heads exposes threshold_contributing_pixels; the old
    threshold_pixels works with a DeprecationWarning and finds the same heads; giving
    both is refused.
    """
    dem, flow = flow_fixture()
    parameters = inspect.signature(lsdtt3.extract_channel_heads).parameters
    assert "threshold_contributing_pixels" in parameters
    assert "threshold_pixels" not in parameters
    canonical = lsdtt3.extract_channel_heads(
        dem, flow, threshold_contributing_pixels=3
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        legacy = lsdtt3.extract_channel_heads(dem, flow, threshold_pixels=3)
    assert [warning.category for warning in caught] == [DeprecationWarning]
    assert caught[0].filename == __file__
    np.testing.assert_array_equal(canonical.pixels, legacy.pixels)
    with pytest.raises(TypeError, match="cannot both be supplied"):
        lsdtt3.extract_channel_heads(
            dem, flow, threshold_contributing_pixels=3, threshold_pixels=3
        )
