import unittest

import numpy as np

import lsdtt3


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


if __name__ == "__main__":
    unittest.main()
