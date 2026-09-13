import inspect
import unittest
import warnings

import numpy as np

import lsdtt3


class SedimentRoutingTest(unittest.TestCase):
    def setUp(self):
        meta = lsdtt3.RasterMeta.from_shape(20, 20, cell_size=10.0)
        rows, cols = np.indices((20, 20), dtype=np.float32)
        data = 100.0 - rows - 0.02 * (cols - 10.0) ** 2
        self.dem = lsdtt3.Raster.from_numpy(
            data, np.ones(data.shape, dtype=bool), meta
        )
        flow = lsdtt3.FlowInfo(
            self.dem, boundary_conditions="nnon", allow_pits_as_outlets=True
        )
        self.network = flow.channel_network(flow.channel_scope(5))

    def test_routes_and_returns_reach_and_raster_diagnostics(self):
        result = lsdtt3.route_sediment(
            self.network, self.dem, n_timesteps=2, timestep_seconds=3600.0
        )
        self.assertEqual(result.n_timesteps, 2)
        self.assertEqual(result.n_reaches, self.network.n_junctions)
        self.assertGreaterEqual(result.total_flux_out_m3, 0.0)
        self.assertEqual(result.reaches.shape[0], self.network.n_junctions)
        self.assertEqual(result.deposition.shape, self.dem.shape)
        self.assertEqual(result.active_layer_d50.shape, self.dem.shape)

    def test_dimension_keyword_compatibility_preserves_routing_results(self):
        parameters = inspect.signature(lsdtt3.route_sediment).parameters
        for canonical, legacy in (
            ("active_layer_depth_metres", "active_layer_depth_m"),
            ("d50_metres", "d50_m"),
            ("d84_metres", "d84_m"),
            ("d90_metres", "d90_m"),
        ):
            self.assertIn(canonical, parameters)
            self.assertNotIn(legacy, parameters)

        values = dict(
            active_layer_depth_metres=0.08,
            d50_metres=0.02,
            d84_metres=0.05,
            d90_metres=0.06,
        )
        common = dict(n_timesteps=2, timestep_seconds=3600.0)
        canonical = lsdtt3.route_sediment(self.network, self.dem, **values, **common)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            legacy = lsdtt3.route_sediment(
                self.network, self.dem,
                active_layer_depth_m=0.08, d50_m=0.02, d84_m=0.05, d90_m=0.06,
                **common,
            )
        self.assertEqual([warning.category for warning in caught], [DeprecationWarning] * 4)
        self.assertTrue(all(warning.filename == __file__ for warning in caught))
        self.assertEqual(canonical.total_flux_out_m3, legacy.total_flux_out_m3)
        np.testing.assert_allclose(canonical.deposition.to_numpy(), legacy.deposition.to_numpy())
        np.testing.assert_allclose(
            canonical.active_layer_d50.to_numpy(), legacy.active_layer_d50.to_numpy()
        )
        with self.assertRaisesRegex(TypeError, "cannot both be supplied"):
            lsdtt3.route_sediment(
                self.network, self.dem, active_layer_depth_metres=0.08,
                active_layer_depth_m=0.08, **common,
            )

    def test_rejects_invalid_discharge_shape(self):
        with self.assertRaisesRegex(ValueError, "n_timesteps, n_junctions"):
            lsdtt3.route_sediment(
                self.network, self.dem, n_timesteps=2,
                discharge=np.ones((1, self.network.n_junctions)),
            )


if __name__ == "__main__":
    unittest.main()
