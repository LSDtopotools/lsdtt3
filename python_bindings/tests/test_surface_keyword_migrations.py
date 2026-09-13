import inspect
import warnings

import numpy as np
import pytest
from pyproj import CRS

import lsdtt3


def raster_and_flow():
    meta = lsdtt3.RasterMeta.from_shape(
        64, 64, cell_size=10.0, crs_wkt=CRS.from_epsg(32630).to_wkt(), epsg_code=32630
    )
    rows, cols = np.indices((64, 64), dtype=np.float32)
    dem = lsdtt3.Raster.from_numpy(
        -rows + 2.0 * np.sin(cols / 4.0), np.ones((64, 64), dtype=bool), meta
    )
    return dem, lsdtt3.FlowInfo(dem, boundary_conditions="nnbo", n_threads=1)


def test_surface_keyword_migrations_preserve_native_signatures():
    channel_head_parameters = inspect.signature(lsdtt3.extract_channel_heads).parameters
    assert "surface_fitting_radius_metres" in channel_head_parameters
    assert "surface_fitting_radius_m" not in channel_head_parameters

    wiener_parameters = inspect.signature(lsdtt3.wiener_filter).parameters
    assert "signal_wavelength_long_metres" in wiener_parameters
    assert "signal_wavelength_short_metres" in wiener_parameters
    assert "signal_wavelength_long_m" not in wiener_parameters
    assert "signal_wavelength_short_m" not in wiener_parameters


def test_channel_head_legacy_radius_warns_and_matches_canonical():
    dem, flow = raster_and_flow()
    canonical = lsdtt3.extract_channel_heads(
        dem, flow, threshold_pixels=3, surface_fitting_radius_metres=8.0
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        legacy = lsdtt3.extract_channel_heads(
            dem, flow, threshold_pixels=3, surface_fitting_radius_m=8.0
        )
    assert len(caught) == 1
    assert caught[0].category is DeprecationWarning
    assert caught[0].filename == __file__
    np.testing.assert_array_equal(canonical.pixels, legacy.pixels)
    with pytest.raises(TypeError, match="cannot both be supplied"):
        lsdtt3.extract_channel_heads(
            dem, flow, surface_fitting_radius_metres=8.0,
            surface_fitting_radius_m=8.0,
        )


@pytest.mark.skipif(not lsdtt3.spectral_available(), reason="requires an FFTW-enabled build")
def test_channel_head_wiener_radius_legacy_keyword_preserves_results():
    dem, flow = raster_and_flow()
    common = dict(
        method="wiener", pruning_drainage_area_m2=100.0,
        connected_components_threshold_pixels=1, diagnostics=True,
    )
    canonical = lsdtt3.extract_channel_heads(
        dem, flow, surface_fitting_radius_metres=20.0, **common
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        legacy = lsdtt3.extract_channel_heads(
            dem, flow, surface_fitting_radius_m=20.0, **common
        )
    assert [warning.category for warning in caught] == [DeprecationWarning]
    np.testing.assert_array_equal(canonical.pixels, legacy.pixels)
    np.testing.assert_allclose(
        canonical.tangential_curvature.to_numpy(), legacy.tangential_curvature.to_numpy()
    )

    perturbed = lsdtt3.extract_channel_heads(
        dem, flow, surface_fitting_radius_metres=40.0, **common
    )
    assert not np.allclose(
        canonical.tangential_curvature.to_numpy(), perturbed.tangential_curvature.to_numpy()
    )


@pytest.mark.skipif(not lsdtt3.spectral_available(), reason="requires an FFTW-enabled build")
def test_wiener_legacy_wavelengths_warn_and_match_canonical():
    dem, _ = raster_and_flow()
    canonical = lsdtt3.wiener_filter(
        dem, signal_wavelength_long_metres=300.0,
        signal_wavelength_short_metres=30.0, n_threads=1,
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        legacy = lsdtt3.wiener_filter(
            dem, signal_wavelength_long_m=300.0,
            signal_wavelength_short_m=30.0, n_threads=1,
        )
    assert [warning.category for warning in caught] == [DeprecationWarning] * 2
    assert all(warning.filename == __file__ for warning in caught)
    np.testing.assert_allclose(canonical.to_numpy(), legacy.to_numpy())


def test_wiener_canonical_and_legacy_wavelengths_conflict():
    with pytest.raises(TypeError, match="cannot both be supplied"):
        lsdtt3.wiener_filter(
            None, signal_wavelength_long_metres=300.0,
            signal_wavelength_long_m=300.0,
        )
