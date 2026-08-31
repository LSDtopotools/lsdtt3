import unittest

import numpy as np

import lsdtt3


class LandscapeModelTest(unittest.TestCase):
    def setUp(self):
        self.meta = lsdtt3.RasterMeta.from_shape(15, 15, cell_size=10.0)
        rows, cols = np.indices((15, 15), dtype=np.float32)
        data = 200.0 - 0.8 * rows - 0.05 * (cols - 7.0) ** 2
        self.dem = lsdtt3.Raster.from_numpy(
            data, np.ones(data.shape, dtype=bool), self.meta
        )

    def test_fluvial_model_advances_and_returns_erosion_diagnostics(self):
        model = lsdtt3.LandscapeModel(
            self.dem, timestep_years=10.0, erodibility=1.0e-4,
            uplift_rate=0.0, hillslope_model="none", boundary_conditions="nnon",
        )
        result = model.step(compute_erosion_rate=True)
        self.assertEqual(result.iteration, 1)
        self.assertEqual(model.iteration, 1)
        self.assertEqual(result.time_years, 10.0)
        self.assertIsInstance(result.erosion_rate, lsdtt3.Raster)
        self.assertGreaterEqual(result.max_local_elevation_change, 0.0)
        np.testing.assert_allclose(
            model.elevation.to_numpy()[-1], self.dem.to_numpy()[-1]
        )

    def test_hillslope_only_run_reuses_model_state(self):
        channels = np.zeros((15, 15), dtype=np.float32)
        channels[-1] = 1.0
        channel_raster = lsdtt3.Raster.from_numpy(
            channels, np.ones(channels.shape, dtype=bool), self.meta
        )
        model = lsdtt3.LandscapeModel(
            self.dem, timestep_years=2.0, fluvial_method="none",
            uplift_rate=0.0, hillslope_model="linear", hillslope_diffusivity=0.05,
            channel_mask=channel_raster, hold_channel_mask_fixed=True,
            boundary_conditions="nnnn",
        )
        result = model.run(2)
        self.assertEqual(result.iteration, 2)
        self.assertEqual(model.time_years, 4.0)
        np.testing.assert_allclose(
            result.elevation.to_numpy()[-1], self.dem.to_numpy()[-1]
        )


if __name__ == "__main__":
    unittest.main()
