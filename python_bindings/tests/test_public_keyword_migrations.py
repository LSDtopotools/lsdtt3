import warnings

import numpy as np
import pytest

import lsdtt3

pytestmark = [pytest.mark.unit]


def hydraulic_inputs():
    meta = lsdtt3.RasterMeta.from_shape(3, 3, cell_size=1.0)
    valid = np.ones((3, 3), dtype=bool)
    bed = lsdtt3.Raster.from_numpy(
        np.broadcast_to(np.arange(3, dtype=np.float32), (3, 3)).copy(), valid, meta
    )
    states = np.ones((3, 3), dtype=np.float32)
    states[1, 0] = 3.0
    node_state = lsdtt3.Raster.from_numpy(states, valid, meta)
    rainfall = lsdtt3.Raster.from_numpy(np.zeros((3, 3), dtype=np.float32), valid, meta)
    friction = lsdtt3.Raster.from_numpy(np.full((3, 3), 30.0, dtype=np.float32), valid, meta)
    return bed, node_state, rainfall, friction


@pytest.mark.parametrize(
    "solver", [
        lsdtt3.run_graphflood,
        lsdtt3.run_graphflood_optimised,
        lsdtt3.run_graphflood_optimised_v2,
    ],
)
@pytest.mark.metamorphic
@pytest.mark.conformance
@pytest.mark.negative
def test_graphflood_legacy_tolerances_warn_and_match(solver):
    inputs = hydraulic_inputs()
    common = dict(max_steps=2, run_until_converged=False)
    canonical = solver(
        *inputs, median_dhw_convergence_tolerance_metres=2.0e-7,
        median_dhw_relative_convergence_tolerance=3.0e-5, **common,
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        legacy = solver(
            *inputs, tol_median_abs_dhw=2.0e-7,
            tol_median_rel_dhw=3.0e-5, **common,
        )
    assert [warning.category for warning in caught] == [DeprecationWarning] * 2
    assert all(warning.filename == __file__ for warning in caught)
    np.testing.assert_array_equal(
        canonical.water_depth.to_numpy(), legacy.water_depth.to_numpy()
    )
    with pytest.raises(TypeError, match="cannot both be supplied"):
        solver(
            *inputs, median_dhw_convergence_tolerance_metres=2.0e-7,
            tol_median_abs_dhw=2.0e-7, **common,
        )


@pytest.mark.parametrize("accelerator", [
    lsdtt3.run_graphflood_multigrid,
    lsdtt3.run_graphflood_induced,
])
@pytest.mark.conformance
@pytest.mark.negative
def test_accelerator_legacy_tolerances_normalize_before_dispatch(accelerator):
    calls = []

    def backend(*args, **kwargs):
        calls.append(kwargs.copy())
        return lsdtt3.run_graphflood_optimised_v2(*args, **kwargs)

    bed, states, rainfall, friction = hydraulic_inputs()
    kwargs = dict(
        backend=backend, point_inflows=None, max_steps=1,
        run_until_converged=False,
    )
    if accelerator is lsdtt3.run_graphflood_multigrid:
        kwargs.update(factors=(3, 1), steps_per_level=1)
    else:
        kwargs.update(qin_threshold=0.0, halo_cells=0, max_expansions=0)
    with pytest.warns(DeprecationWarning) as caught:
        accelerator(
            bed, states, rainfall, friction,
            tol_median_abs_dhw=2.0e-7,
            tol_median_rel_dhw=3.0e-5, **kwargs,
        )
    assert len(caught) == 2
    assert calls
    assert all("tol_median_abs_dhw" not in call for call in calls)
    assert all(call["median_dhw_convergence_tolerance_metres"] == 2.0e-7 for call in calls)
    with pytest.raises(TypeError, match="cannot both be supplied"):
        accelerator(
            bed, states, rainfall, friction,
            median_dhw_convergence_tolerance_metres=2.0e-7,
            tol_median_abs_dhw=2.0e-7, **kwargs,
        )
@pytest.mark.metamorphic
@pytest.mark.conformance
@pytest.mark.negative
def test_landscape_threshold_legacy_keyword_warns_and_conflicts():
    meta = lsdtt3.RasterMeta.from_shape(5, 5, cell_size=10.0)
    rows, _ = np.indices((5, 5), dtype=np.float32)
    dem = lsdtt3.Raster.from_numpy(
        -rows, np.ones((5, 5), dtype=bool), meta
    )
    canonical = lsdtt3.LandscapeModel(
        dem, incision_threshold_specific_area=20.0, hillslope_model="none"
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        legacy = lsdtt3.LandscapeModel(
            dem, incision_threshold_area=20.0, hillslope_model="none"
        )
    assert len(caught) == 1
    assert caught[0].category is DeprecationWarning
    assert caught[0].filename == __file__
    np.testing.assert_array_equal(
        canonical.channel_mask().to_numpy(), legacy.channel_mask().to_numpy()
    )
    positional = lsdtt3.LandscapeModel(
        dem, 100.0, 1.0e-3, 0.45, 1.0, 1.0e-4, 1.0e-3,
        "none", False, 0.6, 20.0,
    )
    assert type(positional) is lsdtt3.LandscapeModel

    class Child(lsdtt3.LandscapeModel):
        pass

    child = Child(dem, incision_threshold_specific_area=20.0, hillslope_model="none")
    assert isinstance(child, lsdtt3.LandscapeModel)
    with pytest.raises(TypeError):
        lsdtt3.LandscapeModel(
            dem, incision_threshold_specific_area=20.0,
            incision_threshold_area=20.0,
        )


@pytest.mark.metamorphic
@pytest.mark.conformance
@pytest.mark.negative
def test_flat_relief_window_legacy_keyword_warns_and_conflicts():
    meta = lsdtt3.RasterMeta.from_shape(9, 9, cell_size=1.0)
    values = np.zeros((9, 9), dtype=np.float32)
    values[4, 4] = 10.0
    dem = lsdtt3.Raster.from_numpy(values, np.ones((9, 9), dtype=bool), meta)
    common = dict(
        use_madk2=True, madk2_window_radius_metres=1.0,
        use_residual_iqr=False, use_nvd=False, min_component_pixels=1,
        max_relief_metres=20.0, n_threads=1,
    )
    canonical = lsdtt3.find_flat_surfaces(
        dem, relief_window_radius_metres=2.0, **common
    )
    with pytest.warns(DeprecationWarning, match="relief_window_side_pixels"):
        side = lsdtt3.find_flat_surfaces(
            dem, relief_window_side_pixels=5, **common
        )
    np.testing.assert_array_equal(
        canonical.relief.to_numpy(), side.relief.to_numpy()
    )
    with pytest.warns(DeprecationWarning, match="relief_threshold_metres"):
        lsdtt3.find_flat_surfaces(
            dem, relief_window_radius_metres=2.0, relief_threshold_metres=20.0,
            **{k: v for k, v in common.items() if k != "max_relief_metres"}
        )
    with pytest.warns(DeprecationWarning, match="relief_window_side_pixels"):
        legacy = lsdtt3.find_flat_surfaces(
            dem, relief_window_pixels=5, **common
        )
    np.testing.assert_array_equal(
        canonical.relief.to_numpy(), legacy.relief.to_numpy()
    )
    with pytest.raises(TypeError, match="cannot both be supplied"):
        lsdtt3.find_flat_surfaces(
            dem, relief_window_side_pixels=5, relief_window_pixels=5, **common
        )


@pytest.mark.exact_solution
@pytest.mark.negative
def test_valley_mask_is_required_validated_and_thresholded():
    meta = lsdtt3.RasterMeta.from_shape(9, 9, cell_size=1.0)
    rows, cols = np.indices((9, 9), dtype=np.float32)
    dem = lsdtt3.Raster.from_numpy(
        rows + (cols - 4.0) ** 2, np.ones((9, 9), dtype=bool), meta
    )
    channel = np.column_stack((np.arange(1, 8), np.full(7, 4)))
    with pytest.raises(TypeError):
        lsdtt3.extract_valley_centrelines(dem, [channel])

    wrong_meta = lsdtt3.RasterMeta.from_shape(8, 8, cell_size=1.0)
    wrong = lsdtt3.Raster.from_numpy(
        np.ones((8, 8), dtype=np.float32), np.ones((8, 8), dtype=bool), wrong_meta
    )
    with pytest.raises(ValueError, match="valley_mask must share the DEM grid"):
        lsdtt3.extract_valley_centrelines(dem, [channel], valley_mask=wrong)

    values = np.full((9, 9), 0.4, dtype=np.float32)
    values[:, 2:7] = 0.6
    mask = lsdtt3.Raster.from_numpy(values, np.ones((9, 9), dtype=bool), meta)
    result = lsdtt3.extract_valley_centrelines(
        dem, [channel], valley_mask=mask, centreline_loops=1, n_threads=1
    )
    assert np.count_nonzero(result.valley_mask.to_numpy()) == 45

    empty = lsdtt3.scale(mask, 0.0)
    empty_result = lsdtt3.extract_valley_centrelines(
        dem, [channel], valley_mask=empty, centreline_loops=1, n_threads=1
    )
    assert np.count_nonzero(empty_result.valley_mask.to_numpy()) == 0


@pytest.mark.smoke
def test_classification_floodplain_mask_composes_with_centreline():
    meta = lsdtt3.RasterMeta.from_shape(9, 9, cell_size=1.0)
    rows, cols = np.indices((9, 9), dtype=np.float32)
    dem = lsdtt3.Raster.from_numpy(
        rows + (cols - 4.0) ** 2, np.ones((9, 9), dtype=bool), meta
    )
    labels = np.zeros((9, 9), dtype=np.float32)
    labels[:, 2:7] = 1.0
    hand = np.full((9, 9), 5.0, dtype=np.float32)
    hand[:, 2:7] = 0.5
    sources = np.full((9, 9), -1.0, dtype=np.float32)
    sources[1:8, 4] = 0.0
    valid = np.ones((9, 9), dtype=bool)
    classified = lsdtt3.classify_flat_surfaces(
        lsdtt3.Raster.from_numpy(labels, valid, meta),
        lsdtt3.Raster.from_numpy(hand, valid, meta),
        lsdtt3.Raster.from_numpy(sources, valid, meta),
        max_floodplain_hand_metres=1.0, n_threads=1,
    )
    channel = np.column_stack((np.arange(1, 8), np.full(7, 4)))
    result = lsdtt3.extract_valley_centrelines(
        dem, [channel], valley_mask=classified.floodplain_mask,
        centreline_loops=1, n_threads=1,
    )
    assert result.valley_mask.shape == dem.shape
    assert len(result.centreline_pixels(0)) > 0
