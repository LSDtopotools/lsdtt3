import unittest

import numpy as np
from pyproj import CRS

import lsdtt3


class InternalDrainageTest(unittest.TestCase):
    def setUp(self):
        metadata = lsdtt3.RasterMeta.from_shape(
            9, 9, cell_size=10.0, crs_wkt=CRS.from_epsg(32630).to_wkt(), epsg_code=32630
        )
        surface = np.full((9, 9), 10.0, dtype=np.float32)
        surface[1:8, 1:8] = 5.0
        surface[4, 4] = 0.0
        self.dem = lsdtt3.Raster.from_numpy(surface, np.ones_like(surface, dtype=bool), metadata)

    def test_hierarchy_candidates_and_polygons(self):
        result = lsdtt3.find_internal_drainage(self.dem)

        self.assertEqual(result.n_candidates, 1)
        self.assertEqual(len(result.candidates), 1)
        self.assertEqual(len(result.polygons), 1)
        candidate = result.candidates.iloc[0]
        self.assertEqual((candidate.pit_row, candidate.pit_col), (4, 4))
        self.assertGreater(candidate.depth_m, 0.0)
        self.assertGreater(candidate.volume_m3, 0.0)
        self.assertEqual(result.candidates.crs.to_epsg(), 32630)

    def test_detected_pit_can_be_retained_for_routing(self):
        result = lsdtt3.find_internal_drainage(self.dem, include_polygons=False)
        pits = result.candidates[["pit_row", "pit_col"]].to_numpy(dtype=np.int64)
        retained = lsdtt3.fill_pits_with_sinks(self.dem, pits, min_slope=0.001)

        self.assertEqual(retained.to_numpy()[4, 4], self.dem.to_numpy()[4, 4])
        self.assertTrue(retained.mask_numpy()[4, 4])

    def test_empty_result_keeps_candidate_schema(self):
        result = lsdtt3.find_internal_drainage(self.dem, min_depth_metres=100.0)

        self.assertEqual(result.n_candidates, 0)
        self.assertTrue(result.candidates.empty)
        self.assertIn("sink_id", result.candidates.columns)
        self.assertIn("volume_m3", result.candidates.columns)


if __name__ == "__main__":
    unittest.main()
