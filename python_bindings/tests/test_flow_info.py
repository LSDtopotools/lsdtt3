import unittest

import numpy as np

import lsdtt3


class FlowInfoTest(unittest.TestCase):
    def setUp(self):
        self.meta = lsdtt3.RasterMeta.from_shape(6, 6, cell_size=10.0)
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


if __name__ == "__main__":
    unittest.main()
