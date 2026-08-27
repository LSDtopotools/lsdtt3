import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import lsdtt3


class RasterFetchTest(unittest.TestCase):
    def test_vector_bounds_and_local_helpers(self):
        self.assertEqual(lsdtt3.utm_epsg_for_point(-3.0, 50.0), 32630)
        self.assertEqual(
            lsdtt3.cop30_tiles_for_bbox(-3.1, 50.1, -2.9, 50.2),
            ["N50_00_W004_00", "N50_00_W003_00"],
        )
        self.assertTrue(lsdtt3.opentopography_datasets())

        with TemporaryDirectory() as directory:
            vector = Path(directory) / "bounds.geojson"
            vector.write_text(
                """{
                    "type": "FeatureCollection",
                    "name": "bounds",
                    "crs": {"type": "name", "properties": {"name": "EPSG:4326"}},
                    "features": [{"type": "Feature", "properties": {}, "geometry": {
                        "type": "Polygon",
                        "coordinates": [[[2.0, 48.0], [2.2, 48.0], [2.2, 48.1], [2.0, 48.1], [2.0, 48.0]]]
                    }}]
                }"""
            )
            bounds = lsdtt3.bounds_from_vector(str(vector))

        self.assertTrue(bounds.ok, bounds.message)
        self.assertEqual(bounds.geometry_kind, "polygon")
        self.assertEqual(bounds.n_features, 1)
        self.assertEqual(bounds.bbox, (2.0, 48.0, 2.2, 48.1))

    def test_invalid_fetch_raises(self):
        with self.assertRaisesRegex(RuntimeError, "west must be less"):
            lsdtt3.fetch_dem(2.0, 48.0, 1.0, 48.1)


if __name__ == "__main__":
    unittest.main()
