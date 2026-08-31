import unittest

import numpy as np

import lsdtt3


class HillslopeTransportTest(unittest.TestCase):
    def setUp(self):
        self.meta = lsdtt3.RasterMeta.from_shape(21, 21, cell_size=5.0)
        rows, cols = np.indices((21, 21), dtype=np.float32)
        elevation = 80.0 - 0.1 * rows + 0.02 * (cols - 10.0) ** 2
        self.dem = lsdtt3.Raster.from_numpy(
            elevation, np.ones((21, 21), dtype=bool), self.meta
        )
        channel_data = np.zeros((21, 21), dtype=np.float32)
        channel_data[[0, -1], :] = 1.0
        self.channels = lsdtt3.Raster.from_numpy(
            channel_data, np.ones((21, 21), dtype=bool), self.meta
        )

    def test_steady_state_keeps_channel_elevations_and_reports_convergence(self):
        result = lsdtt3.hillslope_steady_state(
            self.dem, self.channels, diffusivity=0.05, uplift_rate=4.0e-4,
            model="roering", boundary_conditions="nnnn", tolerance=1.0e-4,
            max_iterations=50000, n_threads=1,
        )
        initial = self.dem.to_numpy()
        solved = result.elevation.to_numpy()
        self.assertTrue(result.converged)
        self.assertLessEqual(result.final_residual, 1.0e-4)
        np.testing.assert_allclose(solved[[0, -1]], initial[[0, -1]])

    def test_explicit_and_implicit_transport_keep_channels_fixed(self):
        explicit = lsdtt3.hillslope_transport_explicit(
            self.dem, self.channels, 0.05, 10.0, boundary_conditions="nnnn", n_threads=1
        )
        implicit = lsdtt3.hillslope_transport_implicit(
            self.dem, self.channels, 10.0, diffusivity=0.05,
            boundary_conditions="nnnn", n_threads=1,
        )
        initial = self.dem.to_numpy()
        np.testing.assert_allclose(explicit.to_numpy()[[0, -1]], initial[[0, -1]])
        np.testing.assert_allclose(implicit.elevation.to_numpy()[[0, -1]], initial[[0, -1]])
        self.assertGreaterEqual(implicit.n_substeps, 1)

    def test_rejects_nonbinary_channels(self):
        bad = lsdtt3.Raster.from_numpy(
            np.full((21, 21), 2.0, dtype=np.float32),
            np.ones((21, 21), dtype=bool), self.meta,
        )
        with self.assertRaisesRegex(ValueError, "exactly 0 or 1"):
            lsdtt3.hillslope_transport_explicit(self.dem, bad, 0.05, 1.0)


if __name__ == "__main__":
    unittest.main()
