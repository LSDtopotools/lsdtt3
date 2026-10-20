"""DEM conditioning: smoothing, unit conversion, coregistration, raster algebra.

Authors
-------
B. Gailleton
"""

import unittest

import numpy as np

import lsdtt3


def _raster(values, cell_size=10.0, mask=None):
    values = np.asarray(values, dtype=np.float32)
    rows, cols = values.shape
    metadata = lsdtt3.RasterMeta.from_shape(rows, cols, x_min=500000.0,
                                            y_min=4000000.0, cell_size=cell_size)
    if mask is None:
        mask = np.ones(values.shape, dtype=bool)
    return lsdtt3.Raster.from_numpy(values, mask, metadata)


def _step_with_noise(seed=1):
    rng = np.random.default_rng(seed)
    values = np.where(np.arange(40)[None, :] < 20, 100.0, 200.0) * np.ones((40, 1))
    return values + rng.normal(0.0, 0.5, values.shape)


class RasterAlgebraTest(unittest.TestCase):
    def test_divide_is_nodata_where_the_divisor_is_zero(self):
        a = _raster([[6.0, 8.0], [1.0, 9.0]])
        b = _raster([[3.0, 0.0], [2.0, -3.0]])
        result = lsdtt3.divide(a, b)
        np.testing.assert_array_equal(result.mask_numpy(), [[True, False], [True, True]])
        values = result.to_numpy()
        self.assertAlmostEqual(float(values[0, 0]), 2.0)
        self.assertAlmostEqual(float(values[1, 0]), 0.5)
        self.assertAlmostEqual(float(values[1, 1]), -3.0)

    def test_absolute_gives_the_absolute_difference(self):
        a = _raster([[1.0, 5.0]])
        b = _raster([[4.0, 2.0]])
        np.testing.assert_allclose(lsdtt3.absolute(lsdtt3.subtract(a, b)).to_numpy(), [[3.0, 3.0]])


class SmoothingTest(unittest.TestCase):
    def test_perona_malik_smooths_noise_but_keeps_the_step(self):
        dem = _raster(_step_with_noise())
        result = lsdtt3.perona_malik_smoothing(dem, n_iterations=20)
        smoothed = result["raster"].to_numpy()
        self.assertEqual(result["n_iterations"], 20)
        self.assertGreater(result["k_metres_per_cell"], 0.0)
        # Noise on each plateau shrinks.
        self.assertLess(np.std(smoothed[:, 2:18]), 0.5 * np.std(dem.to_numpy()[:, 2:18]))
        # The 100 m step is an edge: it survives.
        self.assertGreater(float(np.mean(smoothed[:, 21] - smoothed[:, 18])), 95.0)

    def test_perona_malik_k_can_be_stated(self):
        dem = _raster(_step_with_noise())
        result = lsdtt3.perona_malik_smoothing(dem, k_quantile=0.0, k_metres_per_cell=2.0)
        self.assertAlmostEqual(result["k_metres_per_cell"], 2.0)

    def test_perona_malik_rejects_bad_options(self):
        dem = _raster(_step_with_noise())
        with self.assertRaises(ValueError):
            lsdtt3.perona_malik_smoothing(dem, conductance="linear")
        with self.assertRaises(RuntimeError):
            lsdtt3.perona_malik_smoothing(dem, time_step=0.3)

    def test_stencil_difference_is_the_largest_neighbour_step(self):
        ramp = _raster(np.tile(np.arange(6, dtype=np.float32) * 2.0, (6, 1)))
        np.testing.assert_allclose(lsdtt3.stencil_difference_magnitude(ramp).to_numpy()[2:4, 2:4], 2.0)

    def test_coherence_diffusion_zero_radii_mean_one_and_three_cells(self):
        dem = _raster(_step_with_noise(seed=2))
        by_default = lsdtt3.coherence_diffusion(dem)
        explicit = lsdtt3.coherence_diffusion(dem, derivative_radius_metres=10.0,
                                              integration_radius_metres=30.0)
        np.testing.assert_array_equal(by_default.to_numpy(), explicit.to_numpy())
        self.assertLess(np.std(by_default.to_numpy()[:, 2:18]), np.std(dem.to_numpy()[:, 2:18]))


class VerticalUnitsTest(unittest.TestCase):
    def test_feet_are_converted_to_metres(self):
        dem = _raster([[100.0, 1000.0]])
        np.testing.assert_allclose(lsdtt3.vertical_units_to_metres(dem, "feet").to_numpy(),
                                   [[30.48, 304.8]], rtol=1e-6)
        np.testing.assert_allclose(lsdtt3.vertical_units_to_metres(dem, "metres").to_numpy(),
                                   dem.to_numpy())
        with self.assertRaises(RuntimeError):
            lsdtt3.vertical_units_to_metres(dem, "furlongs")


class CoregistrationTest(unittest.TestCase):
    @staticmethod
    def _hills(x_shift=0.0, y_shift=0.0, z_shift=0.0, n=120, cell=10.0):
        x = (np.arange(n) + 0.5) * cell + x_shift
        y = (np.arange(n) + 0.5) * cell + y_shift
        xx, yy = np.meshgrid(x, y)
        surface = (200.0 * np.sin(xx / 170.0) * np.cos(yy / 230.0)
                   + 80.0 * np.sin((xx + yy) / 90.0) + z_shift)
        return _raster(surface, cell_size=cell)

    def test_a_known_shift_is_recovered_and_removed(self):
        reference = self._hills()
        moving = self._hills(x_shift=15.0, y_shift=-10.0, z_shift=3.0)
        before = np.median(np.abs(moving.to_numpy() - reference.to_numpy()))
        result = lsdtt3.coregister_dem(reference, moving, max_iterations=15)
        after_values = result["raster"].to_numpy()
        valid = result["raster"].mask_numpy() & reference.mask_numpy()
        after = np.median(np.abs(after_values[valid] - reference.to_numpy()[valid]))
        self.assertLess(after, 0.1 * before)
        self.assertAlmostEqual(abs(result["dx_metres"]), 15.0, delta=2.0)
        self.assertAlmostEqual(abs(result["dy_metres"]), 10.0, delta=2.0)
        self.assertGreaterEqual(len(result["diagnostics"]), 1)
        self.assertIn("nmad_residual_metres", result["diagnostics"].columns)

    def test_rasters_must_share_a_grid(self):
        reference = self._hills()
        other = _raster(np.zeros((10, 10)), cell_size=20.0)
        with self.assertRaisesRegex(ValueError, "reproject_to_grid"):
            lsdtt3.coregister_dem(reference, other)
