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
        valid = np.ones((5, 5), dtype=bool)
        self.rainfall_rate = lsdtt3.Raster.from_numpy(
            np.full((5, 5), 1.0e-6, dtype=np.float32), valid, self.meta
        )
        self.friction_coefficient = lsdtt3.Raster.from_numpy(
            np.full((5, 5), 1.0 / 0.033, dtype=np.float32), valid, self.meta
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
            self.rainfall_rate,
            self.friction_coefficient,
            boundary_conditions="nnno",
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
                self.bed, self.states, self.rainfall_rate,
                self.friction_coefficient, boundary_conditions="nnnn", max_steps=1
            )

    def test_hydraulic_fill_keeps_a_representable_slope(self):
        """Tiny epsilon must not round flats back to zero hydraulic slope."""
        meta = lsdtt3.RasterMeta.from_shape(3, 3, cell_size=1.0)
        bed = lsdtt3.Raster.from_numpy(
            np.full((3, 3), 1000.0, dtype=np.float32),
            np.ones((3, 3), dtype=bool),
            meta,
        )
        states = np.ones((3, 3), dtype=np.float32)
        states[1, 0] = 3.0
        node_state = lsdtt3.Raster.from_numpy(
            states, np.ones((3, 3), dtype=bool), meta
        )

        result = lsdtt3.run_graphflood(
            bed,
            node_state,
            lsdtt3.Raster.from_numpy(
                np.ones((3, 3), dtype=np.float32), np.ones((3, 3), dtype=bool), meta
            ),
            lsdtt3.Raster.from_numpy(
                np.full((3, 3), 1.0 / 0.033, dtype=np.float32),
                np.ones((3, 3), dtype=bool), meta
            ),
            boundary_conditions="nnno",
            max_steps=1,
            run_until_converged=False,
            n_threads=1,
        )

        # Without nextafter, 1000 + 1e-8 rounds to 1000 in float32 and all
        # non-outlet cells become sinks. The outlet must receive upstream flow.
        self.assertGreater(result.qin.to_numpy()[1, 0], 1.0)


if __name__ == "__main__":
    unittest.main()
