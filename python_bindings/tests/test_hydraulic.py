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

    def test_each_active_island_requires_a_state_three_sink(self):
        no_outlet = lsdtt3.Raster.from_numpy(
            np.ones((5, 5), dtype=np.float32), np.ones((5, 5), dtype=bool), self.meta
        )
        with self.assertRaisesRegex(ValueError, "outlet"):
            lsdtt3.run_graphflood(
                self.bed, no_outlet, self.rainfall_rate,
                self.friction_coefficient, boundary_conditions="nnnn", max_steps=1
            )

    def test_internal_outlet_is_an_explicit_sink(self):
        states = np.ones((5, 5), dtype=np.float32)
        states[2, 2] = 3.0
        internal_state = lsdtt3.Raster.from_numpy(
            states, np.ones((5, 5), dtype=bool), self.meta
        )
        result = lsdtt3.run_graphflood_optimised_v2(
            self.bed, internal_state, self.rainfall_rate, self.friction_coefficient,
            boundary_conditions="nnnn",
            max_steps=1, run_until_converged=False,
        )
        self.assertGreater(result.qin.to_numpy()[2, 2], 0.0)

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

    def test_separate_optimised_graphflood(self):
        result = lsdtt3.run_graphflood_optimised(
            self.bed,
            self.states,
            self.rainfall_rate,
            self.friction_coefficient,
            boundary_conditions="nnno",
            point_inflows=[[0, 4, 0.01]],
            max_steps=1,
            run_until_converged=False,
        )

        expected = 25 * 1.0e-6 * 100.0 + 0.01
        self.assertAlmostEqual(result.qin.to_numpy()[2, 0], expected, places=7)
        self.assertEqual(result.step_count, 1)
        self.assertGreater(result.performance.routing_seconds, 0.0)
        self.assertGreater(result.performance.persistent_working_bytes, 0)
        self.assertGreaterEqual(result.performance.topology_rebuilds, 1)

    def test_v2_multigrid_and_induced_drivers(self):
        options = dict(
            boundary_conditions="nnno",
            point_inflows=[[0, 4, 0.01]],
            max_steps=2,
            run_until_converged=False,
        )
        v2 = lsdtt3.run_graphflood_optimised_v2(
            self.bed, self.states, self.rainfall_rate,
            self.friction_coefficient, dirty_tile_size=2, **options,
        )
        expected = 25 * 1.0e-6 * 100.0 + 0.01
        self.assertAlmostEqual(v2.qin.to_numpy()[2, 0], expected, places=7)
        self.assertGreater(v2.performance.radix_redistributions, 0)
        self.assertGreater(v2.performance.interior_routed_cells, 0)

        multigrid = lsdtt3.run_graphflood_multigrid(
            self.bed, self.states, self.rainfall_rate,
            self.friction_coefficient, factors=(5, 1), steps_per_level=1,
            run_until_converged=False, boundary_conditions="nnno",
        )
        self.assertEqual(multigrid.factors, (5, 1))
        self.assertEqual(len(multigrid.level_results), 2)
        self.assertEqual(multigrid.water_depth.shape, self.bed.shape)

        induced = lsdtt3.run_graphflood_induced(
            self.bed, self.states, self.rainfall_rate,
            self.friction_coefficient, qin_threshold=1.0e-3,
            halo_cells=1, max_expansions=0, **options,
        )
        self.assertLess(induced.active_fraction, 1.0)
        self.assertLess(induced.relative_mass_error, 1.0e-5)
        self.assertAlmostEqual(
            induced.qin.to_numpy()[2, 0], expected, places=6
        )

    def test_multigrid_rejects_unregistered_fields_and_invalid_points(self):
        shifted = lsdtt3.RasterMeta.from_shape(5, 5, x_min=1.0, cell_size=10.0)
        shifted_rain = lsdtt3.Raster.from_numpy(
            np.ones((5, 5), dtype=np.float32), np.ones((5, 5), dtype=bool), shifted
        )
        with self.assertRaisesRegex(ValueError, "rainfall_rate.*bed grid"):
            lsdtt3.run_graphflood_multigrid(
                self.bed, self.states, shifted_rain, self.friction_coefficient,
                factors=(5, 1), steps_per_level=1,
            )

        invalid_friction = lsdtt3.Raster.from_numpy(
            np.ones((5, 5), dtype=np.float32),
            np.pad(np.ones((4, 5), dtype=bool), ((0, 1), (0, 0)), constant_values=False),
            self.meta,
        )
        with self.assertRaisesRegex(ValueError, "friction_coefficient.*active"):
            lsdtt3.run_graphflood_multigrid(
                self.bed, self.states, self.rainfall_rate, invalid_friction,
                factors=(5, 1), steps_per_level=1,
            )

        with self.assertRaisesRegex(ValueError, "integral coordinates"):
            lsdtt3.run_graphflood_multigrid(
                self.bed, self.states, self.rainfall_rate, self.friction_coefficient,
                factors=(5, 1), steps_per_level=1, point_inflows=[[0.5, 1, 0.1]],
            )


if __name__ == "__main__":
    unittest.main()
