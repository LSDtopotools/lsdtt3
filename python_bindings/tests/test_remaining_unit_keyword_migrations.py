import inspect
import warnings

import numpy as np
import pytest

import lsdtt3


def test_crn_depth_keyword_compatibility_preserves_column_results():
    column_parameters = inspect.signature(lsdtt3.CRNColumn.set_erosion_history).parameters
    transient_parameters = inspect.signature(lsdtt3.CRNColumn.transient).parameters
    assert "removal_metres" in column_parameters
    assert "removal_m" not in column_parameters
    assert "sampling_depths_metres" in transient_parameters
    assert "sampling_depths_m" not in transient_parameters

    canonical = lsdtt3.CRNColumn()
    canonical.set_erosion_history(
        [0.0, 1_000.0, -9_999.0], [1.0e-4] * 3, removal_metres=[0.0] * 3
    )
    legacy = lsdtt3.CRNColumn()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        legacy.set_erosion_history(
            [0.0, 1_000.0, -9_999.0], [1.0e-4] * 3, removal_m=[0.0] * 3
        )
        legacy_profile = legacy.transient(sampling_depths_m=[0.0, 1.0])
    canonical_profile = canonical.transient(sampling_depths_metres=[0.0, 1.0])
    assert [warning.category for warning in caught] == [DeprecationWarning] * 2
    assert all(warning.filename == __file__ for warning in caught)
    np.testing.assert_allclose(
        [sample.concentration_be10 for sample in canonical_profile],
        [sample.concentration_be10 for sample in legacy_profile],
    )
    with pytest.raises(TypeError, match="cannot both be supplied"):
        canonical.transient(sampling_depths_metres=[0.0], sampling_depths_m=[0.0])


def test_crn_raster_history_depth_keyword_compatibility():
    meta = lsdtt3.RasterMeta.from_shape(2, 2, cell_size=30.0)
    production = lsdtt3.Raster.from_numpy(
        np.ones((2, 2), dtype=np.float32), np.ones((2, 2), dtype=bool), meta
    )
    column = lsdtt3.CRNColumn()
    column.set_erosion_history(
        [1_000.0, -9_999.0], [1.0e-3, 1.0e-3], removal_metres=[0.0, 0.0]
    )
    parameters = inspect.signature(lsdtt3.crn_transient_history).parameters
    assert "sampling_depth_metres" in parameters
    assert "sampling_depth_m" not in parameters
    canonical = lsdtt3.crn_transient_history(
        column, production, sampling_depth_metres=1.0
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        legacy = lsdtt3.crn_transient_history(column, production, sampling_depth_m=1.0)
    assert [warning.category for warning in caught] == [DeprecationWarning]
    assert caught[0].filename == __file__
    np.testing.assert_allclose(canonical["be10"].to_numpy(), legacy["be10"].to_numpy())


def test_fetch_size_and_synthetic_relief_keyword_compatibility():
    fetch_parameters = inspect.signature(lsdtt3.fetch_dem).parameters
    assert "size_pixels" in fetch_parameters
    assert "size_px" not in fetch_parameters
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        with pytest.raises(RuntimeError, match="west must be less"):
            lsdtt3.fetch_dem(2.0, 48.0, 1.0, 48.1, size_px=10)
    assert [warning.category for warning in caught] == [DeprecationWarning]
    assert caught[0].filename == __file__
    with pytest.raises(TypeError, match="cannot both be supplied"):
        lsdtt3.fetch_dem(2.0, 48.0, 1.0, 48.1, size_pixels=10, size_px=10)

    metadata = lsdtt3.RasterMeta.from_shape(12, 12, cell_size=10.0)
    network_parameters = inspect.signature(lsdtt3.synthetic_network).parameters
    assert "relief_metres" in network_parameters
    assert "relief_m" not in network_parameters
    canonical = lsdtt3.synthetic_network(
        metadata, n_iterations=50, relief_metres=800.0, seed=7
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        legacy = lsdtt3.synthetic_network(
            metadata, n_iterations=50, relief_m=800.0, seed=7
        )
    assert [warning.category for warning in caught] == [DeprecationWarning]
    assert caught[0].filename == __file__
    np.testing.assert_allclose(canonical.elevation.to_numpy(), legacy.elevation.to_numpy())
