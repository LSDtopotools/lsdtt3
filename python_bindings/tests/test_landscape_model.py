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

    def test_cellular_model_reports_resolved_coupling(self):
        model = lsdtt3.LandscapeModel(
            self.dem, timestep_years=1.0, fluvial_method="none", uplift_rate=0.0,
            hillslope_model="cellular", boundary_conditions="nnnn",
        )
        self.assertEqual(model.hillslope_model, "cellular")
        self.assertEqual(model.coupling_regime, "coupled")
        model.set_forcing(1.0e-4, baselevel_fall_rate=2.0e-4)

    def test_threshold_controls_automatic_coupling_and_channel_mask(self):
        common = dict(
            timestep_years=1.0, fluvial_method="none", hillslope_model="none",
            uplift_rate=0.0, boundary_conditions="nnon", n_threads=1,
        )
        zero = lsdtt3.LandscapeModel(
            self.dem, incision_threshold_area=0.0, **common
        )
        positive = lsdtt3.LandscapeModel(
            self.dem, incision_threshold_area=20.0, **common
        )
        zero_mask = zero.channel_mask().to_numpy()
        positive_mask = positive.channel_mask().to_numpy()
        self.assertEqual(zero.coupling_regime, "coupled")
        self.assertEqual(positive.coupling_regime, "channels_as_boundaries")
        self.assertEqual(np.count_nonzero(zero_mask), 0)
        # 20 is compared with A/dx in metres. On this 10 m grid that is a
        # two-contributing-pixel threshold; pixel-count or m2 interpretations
        # produce different exact masks.
        self.assertEqual(np.count_nonzero(positive_mask), 197)

    def test_default_critical_slope_is_point_six_not_point_eight(self):
        meta = lsdtt3.RasterMeta.from_shape(7, 7, cell_size=1.0)
        rows, cols = np.indices((7, 7), dtype=np.float32)
        values = 20.0 - 0.7 * rows - 0.01 * (cols - 3.0) ** 2
        dem = lsdtt3.Raster.from_numpy(
            values, np.ones(values.shape, dtype=bool), meta
        )
        common = dict(
            timestep_years=0.1, erodibility=0.0, fluvial_method="none",
            uplift_rate=0.0, hillslope_model="cidre",
            hillslope_diffusivity=0.05, nonlinear_hillslopes=True,
            incision_threshold_area=0.0, boundary_conditions="nnnn",
            fixed_pixels=[], n_threads=1,
        )
        default = lsdtt3.LandscapeModel(dem, **common).step().elevation.to_numpy()
        point_six = lsdtt3.LandscapeModel(
            dem, critical_hillslope_slope=0.6, **common
        ).step().elevation.to_numpy()
        point_eight = lsdtt3.LandscapeModel(
            dem, critical_hillslope_slope=0.8, **common
        ).step().elevation.to_numpy()
        np.testing.assert_array_equal(default, point_six)
        self.assertGreater(np.max(np.abs(point_six - point_eight)), 1.0e-6)


if __name__ == "__main__":
    unittest.main()
