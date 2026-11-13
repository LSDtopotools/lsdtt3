import inspect
import unittest
import warnings

import numpy as np
import pytest
from pyproj import CRS

import lsdtt3

pytestmark = [pytest.mark.unit]


class InternalDrainageTest(unittest.TestCase):
    def setUp(self):
        metadata = lsdtt3.RasterMeta.from_shape(
            9, 9, cell_size=10.0, crs_wkt=CRS.from_epsg(32630).to_wkt(), epsg_code=32630
        )
        surface = np.full((9, 9), 10.0, dtype=np.float32)
        surface[1:8, 1:8] = 5.0
        surface[4, 4] = 0.0
        self.dem = lsdtt3.Raster.from_numpy(surface, np.ones_like(surface, dtype=bool), metadata)

    @pytest.mark.exact_solution
    @pytest.mark.property
    def test_hierarchy_candidates_and_polygons(self):
        """find_internal_drainage finds the one pit, at (4, 4), with positive depth and
        volume, one candidate and one polygon, in the DEM's CRS.
        """
        result = lsdtt3.find_internal_drainage(self.dem)

        self.assertEqual(result.n_candidates, 1)
        self.assertEqual(len(result.candidates), 1)
        self.assertEqual(len(result.polygons), 1)
        candidate = result.candidates.iloc[0]
        self.assertEqual((candidate.pit_row, candidate.pit_col), (4, 4))
        self.assertGreater(candidate.depth_m, 0.0)
        self.assertGreater(candidate.volume_m3, 0.0)
        self.assertEqual(result.candidates.crs.to_epsg(), 32630)

    @pytest.mark.exact_solution
    def test_detected_pit_can_be_retained_for_routing(self):
        """The detected pit, given to fill_pits_with_sinks, is kept at its own elevation
        and stays valid.
        """
        result = lsdtt3.find_internal_drainage(self.dem, include_polygons=False)
        pits = result.candidates[["pit_row", "pit_col"]].to_numpy(dtype=np.int64)
        retained = lsdtt3.fill_pits_with_sinks(self.dem, pits, min_slope=0.001)

        self.assertEqual(retained.to_numpy()[4, 4], self.dem.to_numpy()[4, 4])
        self.assertTrue(retained.mask_numpy()[4, 4])

    @pytest.mark.exact_solution
    def test_empty_result_keeps_candidate_schema(self):
        """With a minimum depth no pit reaches, the candidate table is empty but keeps its
        columns.
        """
        result = lsdtt3.find_internal_drainage(self.dem, min_depth_metres=100.0)

        self.assertEqual(result.n_candidates, 0)
        self.assertTrue(result.candidates.empty)
        self.assertIn("sink_id", result.candidates.columns)
        self.assertIn("volume_m3", result.candidates.columns)

    @pytest.mark.metamorphic
    @pytest.mark.conformance
    @pytest.mark.negative
    def test_area_and_volume_keyword_compatibility_preserves_candidates(self):
        """find_internal_drainage's signature has min_area_m2 and min_volume_m3; the old
        names still work with DeprecationWarnings and give the same candidates; giving
        both forms is refused.
        """
        parameters = inspect.signature(lsdtt3.find_internal_drainage).parameters
        self.assertIn("min_area_m2", parameters)
        self.assertIn("min_volume_m3", parameters)
        self.assertNotIn("min_area_metres2", parameters)
        self.assertNotIn("min_volume_metres3", parameters)

        canonical = lsdtt3.find_internal_drainage(
            self.dem, min_area_m2=100.0, min_volume_m3=100.0
        )
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            legacy = lsdtt3.find_internal_drainage(
                self.dem, min_area_metres2=100.0, min_volume_metres3=100.0
            )
        self.assertEqual([warning.category for warning in caught], [DeprecationWarning] * 2)
        self.assertTrue(all(warning.filename == __file__ for warning in caught))
        self.assertEqual(canonical.n_candidates, legacy.n_candidates)
        self.assertEqual(
            canonical.candidates[["pit_row", "pit_col"]].values.tolist(),
            legacy.candidates[["pit_row", "pit_col"]].values.tolist(),
        )
        with self.assertRaisesRegex(TypeError, "cannot both be supplied"):
            lsdtt3.find_internal_drainage(
                self.dem, min_area_m2=100.0, min_area_metres2=100.0
            )

    @pytest.mark.property
    @pytest.mark.smoke
    def test_complete_hierarchy_exposes_labels_and_parent_tree(self):
        """depression_hierarchy has at least the ocean and one leaf, a label raster at the
        DEM's shape and a node table with parent ids.
        """
        hierarchy = lsdtt3.depression_hierarchy(self.dem)

        self.assertGreaterEqual(hierarchy.n_nodes, 2)  # ocean plus one leaf
        self.assertEqual(hierarchy.labels.shape, self.dem.shape)
        self.assertIn("parent_id", hierarchy.nodes.columns)
        self.assertTrue(hierarchy.nodes.is_leaf.any())


if __name__ == "__main__":
    unittest.main()
