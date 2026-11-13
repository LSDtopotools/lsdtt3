import unittest

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import LineString

import lsdtt3

pytestmark = [pytest.mark.unit]


class SwathTest(unittest.TestCase):
    def setUp(self):
        self.crs = gpd.GeoSeries([], crs=32612).crs
        metadata = lsdtt3.RasterMeta.from_shape(
            5, 5, x_min=500000.0, y_min=4400000.0, cell_size=10.0,
            crs_wkt=self.crs.to_wkt(),
        )
        self.raster = lsdtt3.Raster.from_numpy(
            np.arange(25, dtype=np.float32).reshape(5, 5),
            np.ones((5, 5), dtype=bool), metadata,
        )
        self.baseline = gpd.GeoDataFrame(
            geometry=[LineString([(500005.0, 4400025.0), (500045.0, 4400025.0)])],
            crs=self.crs,
        )

    @pytest.mark.exact_solution
    def test_swath_distances(self):
        distances = lsdtt3.swath_distances(self.raster, self.baseline, n_threads=1)

        self.assertEqual(
            set(distances), {"along_axis", "perpendicular", "signed_perpendicular"}
        )
        np.testing.assert_allclose(
            distances["along_axis"].to_numpy()[2], [0.0, 10.0, 20.0, 30.0, 40.0]
        )
        np.testing.assert_allclose(
            distances["perpendicular"].to_numpy()[2], 0.0
        )
        self.assertTrue(np.all(distances["signed_perpendicular"].to_numpy()[1] > 0.0))
        self.assertTrue(np.all(distances["signed_perpendicular"].to_numpy()[3] < 0.0))

    @pytest.mark.exact_solution
    def test_swath_profile_table_and_width(self):
        profile = lsdtt3.swath_profile(
            self.raster, self.baseline, half_width_metres=5.0,
            bin_width_metres=10.0, n_threads=1,
        )

        self.assertEqual(profile.along_axis.shape, self.raster.shape)
        self.assertEqual(profile.quantile_levels.shape, (21,))
        self.assertEqual(profile.table.columns.tolist()[:3], [
            "along_axis_centre", "n_pixels", "mean",
        ])
        self.assertEqual(profile.table.columns.tolist()[-1], "p100")
        self.assertEqual(profile.table["n_pixels"].tolist(), [1, 1, 1, 2])
        np.testing.assert_allclose(profile.table["mean"], [10.0, 11.0, 12.0, 13.5])
        np.testing.assert_allclose(profile.table["p50"], profile.table["mean"])

    @pytest.mark.exact_solution
    def test_swath_profile_leaves_out_pixels_beyond_the_ends(self):
        # The baseline runs from the centre of column 1 to the centre of
        # column 3 of the middle row. Columns 0 and 4 lie behind its start and
        # past its end: they are within half_width of an end point but not
        # alongside the line, so they are left out (before 2026-11-13 they
        # were clamped into the first and last bins, giving [2, 3] pixels).
        short = gpd.GeoDataFrame(
            geometry=[LineString([(500015.0, 4400025.0), (500035.0, 4400025.0)])],
            crs=self.crs,
        )
        profile = lsdtt3.swath_profile(
            self.raster, short, half_width_metres=5.0,
            bin_width_metres=10.0, n_threads=1,
        )

        self.assertEqual(profile.table["n_pixels"].tolist(), [1, 2])
        np.testing.assert_allclose(profile.table["mean"], [11.0, 12.5])

    @pytest.mark.negative
    def test_swath_rejects_mismatched_value_grid(self):
        other = lsdtt3.Raster.from_numpy(
            np.ones((5, 5), dtype=np.float32), np.ones((5, 5), dtype=bool),
            lsdtt3.RasterMeta.from_shape(
                5, 5, x_min=500010.0, y_min=4400000.0, cell_size=10.0,
                crs_wkt=self.crs.to_wkt(),
            ),
        )
        with self.assertRaisesRegex(ValueError, "identical grids"):
            lsdtt3.swath_profile(self.raster, self.baseline, values=other)


if __name__ == "__main__":
    unittest.main()
