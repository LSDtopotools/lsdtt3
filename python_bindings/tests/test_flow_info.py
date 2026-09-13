import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import geopandas as gpd
import numpy as np
from pyproj import CRS
from shapely.geometry import box

import lsdtt3


class FlowInfoTest(unittest.TestCase):
    def setUp(self):
        crs = CRS.from_epsg(32630)
        self.meta = lsdtt3.RasterMeta.from_shape(
            6, 6, cell_size=10.0, crs_wkt=crs.to_wkt(), epsg_code=32630
        )
        rows, cols = np.indices((6, 6), dtype=np.float32)
        self.dem = lsdtt3.Raster.from_numpy(
            -rows - 0.01 * cols, np.ones((6, 6), dtype=bool), self.meta
        )
        self.flow = lsdtt3.FlowInfo(self.dem, boundary_conditions="nnbo", n_threads=1)

    def test_topology_and_core_products(self):
        self.assertEqual(self.flow.shape, self.dem.shape)
        self.assertEqual(self.flow.n_outlets, 6)
        self.assertEqual(self.flow.outlets.shape, (6, 2))
        self.assertEqual(self.flow.flat_index(2, 3), 15)
        self.assertEqual(self.flow.row_col(15), (2, 3))
        self.assertEqual(self.flow.receiver_of(0, 0), (1, 0))
        self.assertTrue(self.flow.is_upslope(0, 0, 5, 0))

        area = self.flow.drainage_area().to_numpy()
        np.testing.assert_allclose(area[-1], 600.0)
        self.assertEqual(self.flow.contributing_pixels().shape, self.dem.shape)
        self.assertEqual(self.flow.flow_distance_to_outlet().shape, self.dem.shape)
        self.assertEqual(self.flow.accumulate().shape, self.dem.shape)

    def test_channel_hand_and_reductions(self):
        channels = self.flow.channel_mask(3)
        self.assertEqual(self.flow.sources(3).shape, (6, 2))
        hand = self.flow.height_above_channel(self.dem, channels)
        handed = self.flow.height_above_channel_with_donor(self.dem, channels)
        self.assertEqual(hand.shape, self.dem.shape)
        self.assertEqual(handed["channel_donor"].shape, self.dem.shape)
        self.assertEqual(self.flow.distance_to_channel(3).shape, self.dem.shape)
        self.assertEqual(self.flow.upslope_sum(self.dem).shape, self.dem.shape)
        self.assertEqual(self.flow.upslope_mean(self.dem).shape, self.dem.shape)
        self.assertEqual(self.flow.downstream_sum(self.dem, channels).shape, self.dem.shape)

    def test_forced_outlet_coordinates(self):
        flow = lsdtt3.FlowInfo(
            self.dem, boundary_conditions="nnbo", forced_outlets=[[3, 2]], n_threads=1
        )
        self.assertTrue(flow.is_outlet(3, 2))

    def test_channel_scope_and_junction_network(self):
        scope = self.flow.channel_scope(3)
        self.assertEqual(scope.mode, "all_channels_threshold")
        self.assertEqual(scope.mask.shape, self.dem.shape)
        self.assertEqual(scope.source_id.shape, self.dem.shape)
        self.assertEqual(scope.nodes(0).shape[1], 2)
        self.assertGreater(len(scope.table()), 0)

        network = self.flow.channel_network(scope)
        self.assertEqual(network.channel_mask.shape, self.dem.shape)
        self.assertEqual(network.junction_id.shape, self.dem.shape)
        self.assertGreater(network.n_channel_cells, 0)
        self.assertGreater(network.n_junctions, 0)
        self.assertIn("stream_order", network.junctions().columns)
        self.assertIn("length_metres", network.links().columns)
        self.assertIn("value", network.link_profiles(self.dem).columns)
        self.assertGreater(len(network.link_bearings()), 0)
        self.assertIn("donor_opening_angle_deg", network.junction_angles().columns)
        self.assertEqual(len(network.basin_com_bearings()), network.n_junctions)
        chi = network.chi()
        self.assertGreater(chi.n_basins, 0)
        self.assertIn("chi", chi.nodes().columns)
        self.assertEqual(chi.raster().shape, self.dem.shape)
        fit = chi.fit(self.dem)
        self.assertIn("ksn", fit.segments.columns)
        self.assertIn("disorder", fit.disorder().columns)
        self.assertEqual(self.flow.chi_from_outlets().shape, self.dem.shape)

    def test_fractional_flow_graphs(self):
        for method in ("mfd", "dinf"):
            graph = lsdtt3.FlowGraph(self.dem, method=method, boundary_conditions="nnbo")
            self.assertEqual(graph.shape, self.dem.shape)
            self.assertEqual(graph.drainage_area().shape, self.dem.shape)
            routing = graph.receivers_of(0, 0)
            self.assertGreater(len(routing["pixels"]), 0)
            self.assertAlmostEqual(float(routing["weights"].sum()), 1.0, places=6)
            if method == "dinf":
                self.assertEqual(graph.slope().shape, self.dem.shape)
                self.assertEqual(graph.direction().shape, self.dem.shape)

    def test_wetness_connectivity_and_sti_tables(self):
        area = self.flow.drainage_area()
        gradient = lsdtt3.slope(self.dem)
        twi = lsdtt3.topographic_wetness_index(area, gradient)
        self.assertEqual(twi.shape, self.dem.shape)
        self.assertTrue(np.all(np.isfinite(twi.to_numpy()[twi.mask_numpy()])))

        connectivity = lsdtt3.sediment_connectivity_index(
            self.flow, gradient, gradient, area, area_threshold_m2=300.0
        )
        self.assertEqual(connectivity.index.shape, self.dem.shape)
        self.assertEqual(connectivity.weight.shape, self.dem.shape)
        self.assertEqual(connectivity.upstream.shape, self.dem.shape)
        self.assertEqual(connectivity.downstream.shape, self.dem.shape)

        with TemporaryDirectory() as directory:
            root = Path(directory)
            ksd = root / "ksd.csv"
            ksd.write_text(
                "SERIES_CODE,SSKIB_PROFILE_ID,LUSE,HORZ_TOP,HORZ_BOTTOM,COND_0_HYPRES\n"
                "AB1,p1,CULT,0,100,2\n"
            )
            spr = root / "spr.csv"
            spr.write_text("SERCDE,SPR_PRC\nAB1,20\n")
            lnkd_table = lsdtt3.load_sti_lnkd_table(str(ksd))
            spr_table = lsdtt3.load_sti_spr_table(str(spr))
            self.assertEqual(lnkd_table.summary["pairs_built"], 1)
            self.assertEqual(spr_table.summary["rows_used"], 1)

            vector = root / "categories.geojson"
            gpd.GeoDataFrame(
                {"series": ["AB1"], "land_use": ["CULT"]},
                geometry=[box(0.0, 0.0, 60.0, 60.0)],
                crs=self.meta.crs_wkt,
            ).to_file(vector, driver="GeoJSON")
            series = lsdtt3.rasterize_sti_categories(
                self.dem, str(vector), "series", "categories"
            )
            land_use = lsdtt3.rasterize_sti_categories(
                self.dem, str(vector), "land_use", "categories"
            )
            sti = lsdtt3.soil_topographic_index(
                twi, series, land_use, lnkd_table, spr_table
            )
            self.assertEqual(sti.sti.shape, self.dem.shape)
            self.assertGreater(sti.stats["valid_output_pixels"], 0)

    def test_scope_workflows_and_basin_collection(self):
        scope = self.flow.channel_scope(3)
        valley_mask = lsdtt3.Raster.from_numpy(
            np.ones(self.dem.shape, dtype=np.float32),
            np.ones(self.dem.shape, dtype=bool), self.meta,
        )
        valleys = lsdtt3.extract_valley_centrelines(
            self.dem, scope, corridor_width_metres=100.0, n_threads=1,
            valley_mask=valley_mask,
        )
        self.assertEqual(valleys.valley_mask.shape, self.dem.shape)

        labels = lsdtt3.Raster.from_numpy(
            np.ones(self.dem.shape, dtype=np.float32),
            np.ones(self.dem.shape, dtype=bool), self.meta,
        )
        flat_classes = lsdtt3.classify_flat_surfaces_with_scope(
            labels, self.dem, scope, n_threads=1
        )
        self.assertEqual(flat_classes.surface_class.shape, self.dem.shape)

        basins = self.flow.basins([[5, 0], [3, 0]])
        self.assertEqual(basins.size, 2)
        self.assertEqual(basins.mask(0).shape, self.dem.shape)
        self.assertEqual(basins.pixels(0).shape[1], 2)
        self.assertEqual(int(basins.table.loc[1, "parent_id"]), 0)
        np.testing.assert_array_equal(basins.membership(0, 0), [0, 1])


if __name__ == "__main__":
    unittest.main()
