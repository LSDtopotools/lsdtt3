import unittest

import numpy as np

import lsdtt3


class HydraulicRoutingTest(unittest.TestCase):
    def setUp(self):
        self.meta = lsdtt3.RasterMeta.from_shape(5, 5, cell_size=10.0)
        _, cols = np.indices((5, 5), dtype=np.float32)
        self.bed = lsdtt3.Raster.from_numpy(
            cols, np.ones((5, 5), dtype=bool), self.meta
        )
        states = np.ones((5, 5), dtype=np.float32)
        states[2, 0] = 3  # west edge; the only permitted hydraulic outlet
        self.states = lsdtt3.Raster.from_numpy(
            states, np.ones((5, 5), dtype=bool), self.meta
        )

    def test_fill_and_stationary_route(self):
        filled = lsdtt3.fill_hydraulic_surface(
            self.bed, self.states, boundary_conditions="nnno"
        )
        self.assertEqual(filled.filled_surface.shape, self.bed.shape)
        self.assertEqual(filled.added_depth.shape, self.bed.shape)
        self.assertGreaterEqual(filled.added_volume, 0.0)

        result = lsdtt3.run_graphflood(
            self.bed,
            self.states,
            boundary_conditions="nnno",
            rainfall_rate=1.0e-6,
            point_inflows=[[0, 4, 0.01]],
            max_steps=10,
            run_until_converged=False,
            n_threads=1,
        )
        self.assertEqual(result.step_count, 10)
        self.assertEqual(result.water_depth.shape, self.bed.shape)
        self.assertEqual(result.total_discharge.shape, self.bed.shape)
        self.assertGreaterEqual(result.total_water_volume, 0.0)

    def test_outlet_must_be_on_an_open_edge(self):
        with self.assertRaisesRegex(ValueError, "open"):
            lsdtt3.run_graphflood(
                self.bed, self.states, boundary_conditions="nnnn", max_steps=1
            )


if __name__ == "__main__":
    unittest.main()
