import inspect
import unittest
import warnings

import numpy as np
import pytest

import lsdtt3

pytestmark = [pytest.mark.unit]


class LandslideTest(unittest.TestCase):
    def setUp(self):
        meta = lsdtt3.RasterMeta.from_shape(40, 40, cell_size=10.0)
        row, col = np.indices((40, 40), dtype=np.float32)
        self.dem = lsdtt3.Raster.from_numpy(
            200.0 - row - 0.01 * (col - 20.0) ** 2,
            np.ones((40, 40), dtype=bool), meta,
        )
        self.flow = lsdtt3.FlowInfo(
            self.dem, boundary_conditions="nnon", allow_pits_as_outlets=True
        )

    @pytest.mark.property
    def test_event_and_cosmogenic_response(self):
        event = lsdtt3.trigger_landslides(
            self.dem, self.flow, 1_000.0, channel_threshold_pixels=100, seed=4
        )
        self.assertGreater(event.n_landslides, 0)
        self.assertGreater(event.achieved_volume_m3, 0.0)
        self.assertEqual(event.elevation.shape, self.dem.shape)
        concentration = lsdtt3.landslide_crn_concentration(
            event.excavation_depth, event.failure_time
        )
        self.assertEqual(concentration.shape, self.dem.shape)
        self.assertEqual(
            set(lsdtt3.landslide_crn_concentrations(
                event.excavation_depth, event.failure_time
            )), {"Be10", "Al26", "C14"}
        )

    @pytest.mark.exact_solution
    @pytest.mark.property
    def test_annealing_aware_series(self):
        series = lsdtt3.run_landslide_series(
            self.dem, self.flow,
            end_time_yr=1_000.0,
            mean_recurrence_yr=250.0,
            poisson_timing=False,
            target_volume_per_event_m3=1_000.0,
            channel_threshold_pixels=100,
            anneal_time_yr=0.0,
            seed=4,
        )
        self.assertEqual(series.n_events, 4)
        self.assertGreater(series.n_landslides, 0)
        self.assertGreater(series.total_volume_m3, 0.0)
        self.assertEqual(series.elevation.shape, self.dem.shape)
        self.assertEqual(len(series.event_times_yr), series.n_events)
        be10 = lsdtt3.landslide_crn_concentration(
            series.excavation_depth, series.failure_time,
            report_time_yr=1_000.0,
        )
        self.assertEqual(be10.shape, self.dem.shape)

    @pytest.mark.metamorphic
    @pytest.mark.negative
    @pytest.mark.conformance
    def test_event_budget_tolerance_keyword_compatibility(self):
        parameters = inspect.signature(lsdtt3.run_landslide_series).parameters
        self.assertIn("event_budget_relative_tolerance", parameters)
        self.assertNotIn("budget_tol", parameters)

        common = dict(
            end_time_yr=1_000.0, mean_recurrence_yr=250.0, poisson_timing=False,
            target_volume_per_event_m3=1_000.0, channel_threshold_pixels=100,
            seed=4,
        )
        canonical = lsdtt3.run_landslide_series(
            self.dem, self.flow, event_budget_relative_tolerance=0.01, **common
        )
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            legacy = lsdtt3.run_landslide_series(
                self.dem, self.flow, budget_tol=0.01, **common
            )
        self.assertEqual([warning.category for warning in caught], [DeprecationWarning])
        self.assertEqual(caught[0].filename, __file__)
        self.assertEqual(canonical.total_volume_m3, legacy.total_volume_m3)
        np.testing.assert_allclose(canonical.elevation.to_numpy(), legacy.elevation.to_numpy())
        with self.assertRaisesRegex(TypeError, "cannot both be supplied"):
            lsdtt3.run_landslide_series(
                self.dem, self.flow, event_budget_relative_tolerance=0.01,
                budget_tol=0.01, **common,
            )


if __name__ == "__main__":
    unittest.main()
