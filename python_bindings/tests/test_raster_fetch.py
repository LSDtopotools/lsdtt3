import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import geopandas as gpd
import pytest
from shapely.geometry import Point

import lsdtt3

pytestmark = [pytest.mark.unit]


class RasterFetchTest(unittest.TestCase):
    @pytest.mark.exact_solution
    def test_vector_bounds_and_local_helpers(self):
        """Offline helpers: the UTM EPSG code for a point, the Copernicus tiles covering a
        box, a non-empty OpenTopography dataset list, and the box, geometry type and
        feature count read from a polygon GeoJSON.
        """
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

    @pytest.mark.negative
    def test_invalid_fetch_raises(self):
        """fetch_dem refuses a box whose west edge is east of its east edge."""
        with self.assertRaisesRegex(RuntimeError, "west must be less"):
            lsdtt3.fetch_dem(2.0, 48.0, 1.0, 48.1)

    @pytest.mark.negative
    def test_outlet_fetch_validates_inputs_before_network_access(self):
        """fetch_dem_for_outlets refuses a negative initial radius and an unknown source,
        before any download.
        """
        outlets = gpd.GeoDataFrame(geometry=[Point(25.5, 47.1)], crs="EPSG:4326")
        with self.assertRaisesRegex(RuntimeError, "initial outlet radius must be positive"):
            lsdtt3.fetch_dem_for_outlets(outlets, initial_radius_km=-1.0)
        with self.assertRaisesRegex(ValueError, "source must be"):
            lsdtt3.fetch_dem_for_outlets(outlets, source="not-a-source")

    @pytest.mark.exact_solution
    def test_dem_source_catalogue(self):
        """dem_sources lists the national sources with their surface models (AHN: DSM and
        DTM), collections (Scottish phases; none for 3DEP), the LERC requirement (New
        Zealand) and attribution (Japan).
        """
        sources = {row["dem_source"]: row for row in lsdtt3.dem_sources()}
        for name in ("cop30_aws", "us_3dep", "arctic_arcticdem", "nl_ahn",
                     "gb_sct_lidar", "nz_lidar", "jp_gsi"):
            self.assertIn(name, sources)
        self.assertEqual(sorted(sources["nl_ahn"]["surface_models"]), ["dsm", "dtm"])
        self.assertIn("phase-1", sources["gb_sct_lidar"]["collections"])
        self.assertEqual(sources["us_3dep"]["collections"], [])
        self.assertTrue(sources["nz_lidar"]["requires_lerc"])
        self.assertTrue(sources["jp_gsi"]["attribution"])

    @pytest.mark.negative
    def test_unknown_source_lists_the_catalogue(self):
        """fetch_dem with an unknown source is refused with a message listing the known
        ones.
        """
        with self.assertRaisesRegex(ValueError, "gb_sct_lidar"):
            lsdtt3.fetch_dem(2.0, 48.0, 2.1, 48.1, source="not-a-source")

    @pytest.mark.negative
    def test_source_options_are_checked_before_any_download(self):
        """fetch_dem refuses, before any download, a collection for a source without
        collections, an unknown Scottish phase, and a DTM choice for a source that has
        only a DSM.
        """
        with self.assertRaisesRegex(RuntimeError, "source_collection_choice"):
            lsdtt3.fetch_dem(8.0, 61.0, 8.01, 61.01, source="no_dtm",
                             source_collection_choice="phase-1")
        with self.assertRaisesRegex(RuntimeError, "phase-9"):
            lsdtt3.fetch_dem(-4.9, 55.7, -4.89, 55.71, source="gb_sct_lidar",
                             source_collection_choice="phase-9")
        with self.assertRaisesRegex(RuntimeError, "DSM"):
            lsdtt3.fetch_dem(2.0, 48.0, 2.1, 48.1, surface_model_choice="dtm")

    @pytest.mark.negative
    def test_lerc_source_warns_before_fetching(self):
        """Asking for New Zealand lidar warns about LERC before the fetch fails (outside
        New Zealand).
        """
        # Outside New Zealand, so it fails in pre-flight whether or not this
        # GDAL can read LERC; the warning comes first either way.
        with self.assertWarnsRegex(UserWarning, "LERC"):
            with self.assertRaises(RuntimeError):
                lsdtt3.fetch_dem(2.0, 48.0, 2.1, 48.1, source="nz_lidar")

    @pytest.mark.negative
    def test_outlet_fetch_refuses_high_resolution_sources(self):
        """fetch_dem_for_outlets refuses a high-resolution source such as us_3dep."""
        outlets = gpd.GeoDataFrame(geometry=[Point(-105.5, 40.0)], crs="EPSG:4326")
        with self.assertRaisesRegex(ValueError, "us_3dep"):
            lsdtt3.fetch_dem_for_outlets(outlets, source="us_3dep")


if __name__ == "__main__":
    unittest.main()
