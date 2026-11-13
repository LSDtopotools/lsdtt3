import inspect
import unittest

import numpy as np
import pytest

import lsdtt3

pytestmark = [pytest.mark.unit]


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

    @pytest.mark.property
    def test_fill_and_stationary_route(self):
        """fill_hydraulic_surface returns the filled surface and added depth at the bed's
        shape with a non-negative added volume; run_graphflood with a point inflow runs
        exactly the 10 requested steps and returns depth and discharge at the bed's
        shape.
        """
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

    @pytest.mark.negative
    def test_each_active_island_requires_a_state_three_sink(self):
        """run_graphflood refuses a domain with closed edges and no outlet cell (no state-3
        sink), with a message naming the outlet.
        """
        no_outlet = lsdtt3.Raster.from_numpy(
            np.ones((5, 5), dtype=np.float32), np.ones((5, 5), dtype=bool), self.meta
        )
        with self.assertRaisesRegex(ValueError, "outlet"):
            lsdtt3.run_graphflood(
                self.bed, no_outlet, self.rainfall_rate,
                self.friction_coefficient, boundary_conditions="nnnn", max_steps=1
            )

    @pytest.mark.property
    def test_internal_outlet_is_an_explicit_sink(self):
        """An interior cell given state 3 is an outlet: after one step it receives inflow.
        """
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

    @pytest.mark.exact_solution
    @pytest.mark.property
    @pytest.mark.metamorphic
    @pytest.mark.negative
    def test_graphflood_v2_ghost_state(self):
        """Flux-only ghosts replace incoming flow but evolve their depth."""
        meta = lsdtt3.RasterMeta.from_shape(3, 5, cell_size=1.0)
        valid = np.ones((3, 5), dtype=bool)
        bed_values = np.broadcast_to(
            np.arange(5, dtype=np.float32), (3, 5)
        ).copy()
        bed = lsdtt3.Raster.from_numpy(bed_values, valid, meta)
        states = np.zeros((3, 5), dtype=np.float32)
        states[1, :] = 1.0
        states[1, 0] = 3.0
        states[1, 2] = 2.0
        node_state = lsdtt3.Raster.from_numpy(states, valid, meta)
        zero = lsdtt3.Raster.from_numpy(
            np.zeros((3, 5), dtype=np.float32), valid, meta
        )
        friction = lsdtt3.Raster.from_numpy(
            np.full((3, 5), 1.0 / 0.033, dtype=np.float32), valid, meta
        )
        initial_values = np.zeros((3, 5), dtype=np.float32)
        initial_values[1, 2] = 0.125
        initial = lsdtt3.Raster.from_numpy(initial_values, valid, meta)
        prescribed_values = np.zeros((3, 5), dtype=np.float32)
        prescribed_values[1, 2] = 0.25
        prescribed = lsdtt3.Raster.from_numpy(prescribed_values, valid, meta)

        result = lsdtt3.run_graphflood_optimised_v2(
            bed, node_state, zero, friction, boundary_conditions="nnnn",
            initial_water_depth=initial, prescribed_qin=prescribed,
            point_inflows=[[1, 3, 10.0]], max_steps=3,
            run_until_converged=False,
        )
        qin = result.qin.to_numpy()
        qout = result.qout.to_numpy()
        depth = result.water_depth.to_numpy()

        # The 10 m3/s donor reaches the ghost, then is discarded.
        self.assertEqual(qin[1, 2], np.float32(0.25))
        self.assertGreater(qout[1, 2], qin[1, 2])
        self.assertEqual(qin[1, 1], np.float32(0.25))
        self.assertNotEqual(depth[1, 2], initial_values[1, 2])

        # Existing v2 calls have the same output with None or explicit zeros.
        kwargs = dict(boundary_conditions="nnno", max_steps=2,
                      run_until_converged=False)
        baseline = lsdtt3.run_graphflood_optimised_v2(
            self.bed, self.states, self.rainfall_rate,
            self.friction_coefficient, **kwargs
        )
        explicit_zero = lsdtt3.Raster.from_numpy(
            np.zeros((5, 5), dtype=np.float32), np.ones((5, 5), dtype=bool),
            self.meta,
        )
        with_zero_field = lsdtt3.run_graphflood_optimised_v2(
            self.bed, self.states, self.rainfall_rate,
            self.friction_coefficient, prescribed_qin=explicit_zero, **kwargs
        )
        for name in ("water_depth", "qin", "qout"):
            self.assertEqual(
                getattr(baseline, name).to_numpy().tobytes(),
                getattr(with_zero_field, name).to_numpy().tobytes(),
            )

        invalid_values = prescribed_values.copy()
        invalid_values[1, 2] = -1.0
        invalid = lsdtt3.Raster.from_numpy(invalid_values, valid, meta)
        with self.assertRaisesRegex(ValueError, "prescribed_qin"):
            lsdtt3.run_graphflood_optimised_v2(
                bed, node_state, zero, friction, prescribed_qin=invalid,
                max_steps=1,
            )
        invalid_values[1, 2] = np.inf
        invalid = lsdtt3.Raster.from_numpy(invalid_values, valid, meta)
        with self.assertRaisesRegex(ValueError, "prescribed_qin"):
            lsdtt3.run_graphflood_optimised_v2(
                bed, node_state, zero, friction, prescribed_qin=invalid,
                max_steps=1,
            )

    @pytest.mark.property
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

    @pytest.mark.exact_solution
    @pytest.mark.property
    def test_separate_optimised_graphflood(self):
        """run_graphflood_optimised, one step: the inflow at the outlet cell equals the
        rain on all 25 cells plus the point inflow (exact), and the performance counters
        are filled.
        """
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

    @pytest.mark.metamorphic
    def test_direct_graphflood_current_defaults(self):
        """For each GraphFlood solver, a run with the defaults left out equals a run with
        the current defaults written out, in water depth, inflow, outflow and depth
        change.
        """
        controls = dict(max_steps=2, run_until_converged=False)
        cases = (
            (lsdtt3.run_graphflood, dict(
                method="mfd", boundary_conditions="oooo", mfd_exponent=1.0,
                condition_bed="none", fill_slope_epsilon=1.0e-8,
                check_valid_outlet=True, friction_depth_exponent=5.0 / 3.0,
                friction_slope_exponent=0.5, dt=1.0e-3, min_slope=1.0e-4,
                point_inflows=None, convergence_check_every=10,
                median_dhw_convergence_tolerance_metres=1.0e-6, median_dhw_relative_convergence_tolerance=1.0e-4,
                n_threads=0, initial_water_depth=None,
            )),
            (lsdtt3.run_graphflood_optimised, dict(
                boundary_conditions="oooo", mfd_exponent=1.0,
                condition_bed="none", fill_slope_epsilon=1.0e-8,
                check_valid_outlet=True, friction_depth_exponent=5.0 / 3.0,
                friction_slope_exponent=0.5, dt=1.0e-3, min_slope=1.0e-4,
                point_inflows=None, convergence_check_every=10,
                median_dhw_convergence_tolerance_metres=1.0e-6, median_dhw_relative_convergence_tolerance=1.0e-4,
                reuse_topology=True, initial_water_depth=None,
            )),
            (lsdtt3.run_graphflood_optimised_v2, dict(
                boundary_conditions="oooo", mfd_exponent=1.0,
                condition_bed="none", fill_slope_epsilon=1.0e-8,
                check_valid_outlet=True, friction_depth_exponent=5.0 / 3.0,
                friction_slope_exponent=0.5, dt=1.0e-3, min_slope=1.0e-4,
                point_inflows=None, convergence_check_every=10,
                median_dhw_convergence_tolerance_metres=1.0e-6, median_dhw_relative_convergence_tolerance=1.0e-4,
                reuse_topology=True, dirty_tile_size=64,
                full_validation_every=128, initial_water_depth=None,
                prescribed_qin=None,
            )),
        )
        for solver, explicit_defaults in cases:
            with self.subTest(solver=solver.__name__):
                implicit = solver(
                    self.bed, self.states, self.rainfall_rate,
                    self.friction_coefficient, **controls,
                )
                explicit = solver(
                    self.bed, self.states, self.rainfall_rate,
                    self.friction_coefficient, **controls, **explicit_defaults,
                )
                self.assertEqual(implicit.step_count, 2)
                for field in ("water_depth", "qin", "qout", "dhw"):
                    np.testing.assert_array_equal(
                        getattr(implicit, field).to_numpy(),
                        getattr(explicit, field).to_numpy(),
                    )

    @pytest.mark.conformance
    def test_direct_graphflood_defaults_match_cli(self):
        """The three GraphFlood solvers' defaults for max_steps, convergence_check_every
        and the two convergence tolerances are the CLI's values.
        """
        for solver in (
            lsdtt3.run_graphflood,
            lsdtt3.run_graphflood_optimised,
            lsdtt3.run_graphflood_optimised_v2,
        ):
            parameters = inspect.signature(solver).parameters
            self.assertEqual(parameters["max_steps"].default, 100)
            self.assertEqual(parameters["convergence_check_every"].default, 10)
            self.assertEqual(
                parameters["median_dhw_convergence_tolerance_metres"].default, 1.0e-6
            )
            self.assertEqual(
                parameters["median_dhw_relative_convergence_tolerance"].default, 1.0e-4
            )

    @pytest.mark.exact_solution
    @pytest.mark.metamorphic
    @pytest.mark.characterisation
    def test_graphflood_convergence_defaults_are_sensitive(self):
        """With zero rain, each GraphFlood solver stops on convergence at the first check
        (10 steps) by default and runs the full 100 steps without the convergence test;
        from a uniform initial depth the default tolerances stop after one step, loose
        ones likewise, and strict ones run all 1000 steps.
        """
        zero_rain = lsdtt3.Raster.from_numpy(
            np.zeros((5, 5), dtype=np.float32),
            np.ones((5, 5), dtype=bool), self.meta,
        )
        for solver in (
            lsdtt3.run_graphflood,
            lsdtt3.run_graphflood_optimised,
            lsdtt3.run_graphflood_optimised_v2,
        ):
            with self.subTest(solver=solver.__name__):
                result = solver(
                    self.bed, self.states, zero_rain,
                    self.friction_coefficient,
                )
                self.assertEqual(result.step_count, 10)

        for solver in (
            lsdtt3.run_graphflood,
            lsdtt3.run_graphflood_optimised,
            lsdtt3.run_graphflood_optimised_v2,
        ):
            with self.subTest(max_steps_solver=solver.__name__):
                result = solver(
                    self.bed, self.states, zero_rain,
                    self.friction_coefficient, run_until_converged=False,
                )
                self.assertEqual(result.step_count, 100)

        initial_depth = lsdtt3.Raster.from_numpy(
            np.full((5, 5), 0.01, dtype=np.float32),
            np.ones((5, 5), dtype=bool), self.meta,
        )
        common = dict(
            max_steps=1000, run_until_converged=True,
            convergence_check_every=1, initial_water_depth=initial_depth,
        )
        for solver in (
            lsdtt3.run_graphflood,
            lsdtt3.run_graphflood_optimised,
            lsdtt3.run_graphflood_optimised_v2,
        ):
            with self.subTest(tolerance_solver=solver.__name__):
                default_tolerance = solver(
                    self.bed, self.states, zero_rain,
                    self.friction_coefficient, **common,
                )
                loose = solver(
                    self.bed, self.states, zero_rain,
                    self.friction_coefficient,
                    median_dhw_convergence_tolerance_metres=1.0,
                    median_dhw_relative_convergence_tolerance=1.0, **common,
                )
                strict = solver(
                    self.bed, self.states, zero_rain,
                    self.friction_coefficient,
                    median_dhw_convergence_tolerance_metres=1.0e-7,
                    median_dhw_relative_convergence_tolerance=1.0e-5, **common,
                )
                self.assertEqual(default_tolerance.step_count, 1)
                self.assertEqual(loose.step_count, default_tolerance.step_count)
                self.assertEqual(strict.step_count, 1000)

    @pytest.mark.exact_solution
    def test_accelerators_forward_backend_options(self):
        """run_graphflood_multigrid and run_graphflood_induced pass every solver option
        through to the backend unchanged, with the steps per level (multigrid) and a
        first short run without convergence (induced).
        """
        calls = []

        def recording_backend(*args, **kwargs):
            calls.append(kwargs.copy())
            return lsdtt3.run_graphflood_optimised_v2(*args, **kwargs)

        forwarded = dict(
            boundary_conditions="nnno", mfd_exponent=1.25,
            condition_bed="none", fill_slope_epsilon=2.0e-7,
            check_valid_outlet=True, friction_depth_exponent=1.5,
            friction_slope_exponent=0.4, dt=2.0e-3, min_slope=3.0e-4,
            run_until_converged=False, convergence_check_every=7,
            median_dhw_convergence_tolerance_metres=2.0e-7,
            median_dhw_relative_convergence_tolerance=3.0e-5,
            reuse_topology=False, dirty_tile_size=2,
            full_validation_every=0,
        )
        lsdtt3.run_graphflood_multigrid(
            self.bed, self.states, self.rainfall_rate,
            self.friction_coefficient, backend=recording_backend,
            factors=(5, 1), steps_per_level=(1, 2),
            point_inflows=[[0, 4, 0.01]], **forwarded,
        )
        self.assertEqual([call["max_steps"] for call in calls], [1, 2])
        for call in calls:
            for key, value in forwarded.items():
                self.assertEqual(call[key], value)

        calls.clear()
        lsdtt3.run_graphflood_induced(
            self.bed, self.states, self.rainfall_rate,
            self.friction_coefficient, qin_threshold=1.0e-3,
            backend=recording_backend, halo_cells=1, max_expansions=0,
            point_inflows=[[0, 4, 0.01]], max_steps=2, **forwarded,
        )
        self.assertEqual(len(calls), 2)
        self.assertEqual(calls[0]["max_steps"], 1)
        self.assertFalse(calls[0]["run_until_converged"])
        self.assertEqual(calls[1]["max_steps"], 2)
        for call in calls:
            for key, value in forwarded.items():
                if key == "run_until_converged" and call is calls[0]:
                    continue
                self.assertEqual(call[key], value)

    @pytest.mark.exact_solution
    @pytest.mark.property
    def test_v2_multigrid_and_induced_drivers(self):
        """run_graphflood_optimised_v2 gives the exact outlet inflow (rain plus point
        inflow) and fills its counters; multigrid returns one result per level at the
        bed's shape; induced activates only part of the domain, conserves mass to 1e-5
        and gives the same outlet inflow.
        """
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

    @pytest.mark.negative
    def test_multigrid_rejects_unregistered_fields_and_invalid_points(self):
        """run_graphflood_multigrid refuses rain on a shifted grid, a friction raster with
        nodata in the active domain, and point inflows at non-integer cell coordinates,
        each with its message.
        """
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
