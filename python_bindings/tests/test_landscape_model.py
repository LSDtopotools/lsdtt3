import unittest

import numpy as np
import pytest

import lsdtt3

pytestmark = [pytest.mark.unit]


class LandscapeModelTest(unittest.TestCase):
    def setUp(self):
        self.meta = lsdtt3.RasterMeta.from_shape(15, 15, cell_size=10.0)
        rows, cols = np.indices((15, 15), dtype=np.float32)
        data = 200.0 - 0.8 * rows - 0.05 * (cols - 7.0) ** 2
        self.dem = lsdtt3.Raster.from_numpy(
            data, np.ones(data.shape, dtype=bool), self.meta
        )

    @pytest.mark.property
    def test_fluvial_model_advances_and_returns_erosion_diagnostics(self):
        model = lsdtt3.LandscapeModel(
            self.dem, timestep_years=10.0, erodibility=1.0e-4,
            uplift_rate=0.0, hillslope_model="none", boundary_conditions="nnon",
        )
        result = model.step(compute_erosion_rate=True)
        self.assertEqual(result.iteration, 1)
        self.assertEqual(model.iteration, 1)
        self.assertEqual(result.time_years, 10.0)
        self.assertIsInstance(result.erosion_rate, lsdtt3.Raster)
        self.assertGreaterEqual(result.max_local_elevation_change, 0.0)
        np.testing.assert_allclose(
            model.elevation.to_numpy()[-1], self.dem.to_numpy()[-1]
        )

    @pytest.mark.property
    def test_hillslope_only_run_reuses_model_state(self):
        channels = np.zeros((15, 15), dtype=np.float32)
        channels[-1] = 1.0
        channel_raster = lsdtt3.Raster.from_numpy(
            channels, np.ones(channels.shape, dtype=bool), self.meta
        )
        model = lsdtt3.LandscapeModel(
            self.dem, timestep_years=2.0, fluvial_method="none",
            uplift_rate=0.0, hillslope_model="linear", hillslope_diffusivity=0.05,
            channel_mask=channel_raster, hold_channel_mask_fixed=True,
            boundary_conditions="nnnn",
        )
        result = model.run(2)
        self.assertEqual(result.iteration, 2)
        self.assertEqual(model.time_years, 4.0)
        np.testing.assert_allclose(
            result.elevation.to_numpy()[-1], self.dem.to_numpy()[-1]
        )

    @pytest.mark.exact_solution
    def test_cellular_model_reports_resolved_coupling(self):
        model = lsdtt3.LandscapeModel(
            self.dem, timestep_years=1.0, fluvial_method="none", uplift_rate=0.0,
            hillslope_model="cellular", boundary_conditions="nnnn",
        )
        self.assertEqual(model.hillslope_model, "cellular")
        self.assertEqual(model.coupling_regime, "coupled")
        model.set_forcing(1.0e-4, baselevel_fall_rate=2.0e-4)

    @pytest.mark.exact_solution
    @pytest.mark.characterisation
    def test_threshold_controls_automatic_coupling_and_channel_mask(self):
        common = dict(
            timestep_years=1.0, fluvial_method="none", hillslope_model="none",
            uplift_rate=0.0, boundary_conditions="nnon", n_threads=1,
        )
        zero = lsdtt3.LandscapeModel(
            self.dem, incision_threshold_specific_area=0.0, **common
        )
        positive = lsdtt3.LandscapeModel(
            self.dem, incision_threshold_specific_area=20.0, **common
        )
        zero_mask = zero.channel_mask().to_numpy()
        positive_mask = positive.channel_mask().to_numpy()
        self.assertEqual(zero.coupling_regime, "coupled")
        self.assertEqual(positive.coupling_regime, "channels_as_boundaries")
        self.assertEqual(np.count_nonzero(zero_mask), 0)
        # 20 is compared with A/dx in metres. On this 10 m grid that is a
        # two-contributing-pixel threshold; pixel-count or m2 interpretations
        # produce different exact masks.
        self.assertEqual(np.count_nonzero(positive_mask), 197)

    @pytest.mark.metamorphic
    def test_default_critical_slope_is_point_eight(self):
        meta = lsdtt3.RasterMeta.from_shape(7, 7, cell_size=1.0)
        rows, cols = np.indices((7, 7), dtype=np.float32)
        values = 20.0 - 0.7 * rows - 0.01 * (cols - 3.0) ** 2
        dem = lsdtt3.Raster.from_numpy(
            values, np.ones(values.shape, dtype=bool), meta
        )
        common = dict(
            timestep_years=0.1, fluvial_method="none",
            uplift_rate=0.0, hillslope_model="cidre",
            hillslope_diffusivity=0.05, nonlinear_hillslopes=True,
            incision_threshold_specific_area=0.0, boundary_conditions="nnnn",
            fixed_pixels=[], n_threads=1,
        )
        default = lsdtt3.LandscapeModel(dem, **common).step().elevation.to_numpy()
        point_six = lsdtt3.LandscapeModel(
            dem, critical_hillslope_slope=0.6, **common
        ).step().elevation.to_numpy()
        point_eight = lsdtt3.LandscapeModel(
            dem, critical_hillslope_slope=0.8, **common
        ).step().elevation.to_numpy()
        np.testing.assert_array_equal(default, point_eight)
        self.assertGreater(np.max(np.abs(point_six - point_eight)), 1.0e-6)




class LandscapeModelErodibilityTest(unittest.TestCase):
    """Per-pixel erodibility (K) fields on LandscapeModel."""

    def setUp(self):
        self.meta = lsdtt3.RasterMeta.from_shape(15, 15, cell_size=10.0)
        rows, cols = np.indices((15, 15), dtype=np.float32)
        data = 200.0 - 0.8 * rows - 0.05 * (cols - 7.0) ** 2
        self.dem = lsdtt3.Raster.from_numpy(
            data, np.ones(data.shape, dtype=bool), self.meta
        )
        self.common = dict(
            timestep_years=100.0, erodibility=2.0e-4, uplift_rate=1.0e-4,
            hillslope_model="none", boundary_conditions="nnon", n_threads=1,
        )

    def _run(self, steps=3, **kwargs):
        model = lsdtt3.LandscapeModel(self.dem, **{**self.common, **kwargs})
        return model.run(steps).elevation.to_numpy()

    @pytest.mark.metamorphic
    def test_erodibility_uniform_field_equals_scalar_bitwise(self):
        for method in ("d8", "dinf"):
            scalar = self._run(fluvial_method=method)
            uniform = self._run(
                fluvial_method=method,
                erodibility_field=np.full((15, 15), 2.0e-4, dtype=np.float32),
            )
            self.assertEqual(scalar.tobytes(), uniform.tobytes(), method)
            as_raster = self._run(
                fluvial_method=method,
                erodibility_field=lsdtt3.Raster.from_numpy(
                    np.full((15, 15), 2.0e-4, dtype=np.float32),
                    np.ones((15, 15), dtype=bool), self.meta,
                ),
            )
            self.assertEqual(scalar.tobytes(), as_raster.tobytes(), method)

    @pytest.mark.metamorphic
    def test_set_erodibility_field_applies_to_next_step_and_none_restores(self):
        field = np.full((15, 15), 2.0e-4, dtype=np.float32)
        field[:7] = 8.0e-4
        model = lsdtt3.LandscapeModel(self.dem, **self.common)
        reference = lsdtt3.LandscapeModel(self.dem, **self.common)
        first = model.step().elevation.to_numpy()
        self.assertEqual(first.tobytes(), reference.step().elevation.to_numpy().tobytes())
        model.set_erodibility_field(field)
        changed = model.step().elevation.to_numpy()
        self.assertNotEqual(changed.tobytes(), reference.step().elevation.to_numpy().tobytes())
        model.set_erodibility_field(None)
        self.assertIsNone(model.erodibility_field)
        scalar_again = lsdtt3.LandscapeModel(
            lsdtt3.Raster.from_numpy(changed, np.ones((15, 15), dtype=bool), self.meta),
            **self.common,
        )
        self.assertEqual(
            model.step().elevation.to_numpy().tobytes(),
            scalar_again.step().elevation.to_numpy().tobytes(),
        )

    @pytest.mark.property
    def test_erodibility_field_property_round_trips(self):
        field = np.linspace(1.0e-4, 5.0e-4, 225, dtype=np.float32).reshape(15, 15)
        model = lsdtt3.LandscapeModel(self.dem, erodibility_field=field, **self.common)
        self.assertIsInstance(model.erodibility_field, lsdtt3.Raster)
        self.assertEqual(model.erodibility_field.to_numpy().tobytes(), field.tobytes())

    @pytest.mark.negative
    def test_bad_erodibility_fields_raise_value_error(self):
        bad = {}
        bad["shape"] = np.full((14, 15), 2.0e-4, dtype=np.float32)
        for name, value in (("zero", 0.0), ("negative", -1.0e-4), ("nan", np.nan)):
            f = np.full((15, 15), 2.0e-4, dtype=np.float32)
            f[3, 4] = value
            bad[name] = f
        other = lsdtt3.RasterMeta.from_shape(15, 15, cell_size=12.0)
        bad["other_grid"] = lsdtt3.Raster.from_numpy(
            np.full((15, 15), 2.0e-4, dtype=np.float32), np.ones((15, 15), dtype=bool), other
        )
        for name, field in bad.items():
            with self.subTest(name=name), self.assertRaises(ValueError):
                lsdtt3.LandscapeModel(self.dem, erodibility_field=field, **self.common)
        with self.assertRaises(ValueError):
            lsdtt3.LandscapeModel(
                self.dem, erodibility_field=np.full((15, 15), 2.0e-4, dtype=np.float32),
                erodibility_m_exponent=0.5, **self.common,
            )

    @pytest.mark.negative
    def test_scalar_erodibility_must_be_positive(self):
        for value in (0.0, -1.0e-3):
            for key in ("incision_threshold_specific_area", "incision_threshold_area"):
                kwargs = {**self.common, "erodibility": value, key: 0.0}
                with self.subTest(value=value, overload=key), \
                        self.assertRaisesRegex(ValueError, "erodibility must be > 0"):
                    lsdtt3.LandscapeModel(self.dem, **kwargs)

    @pytest.mark.exact_solution
    def test_erodibility_strip_steady_slopes_follow_k(self):
        # A1 at 40x3: m = 0, so at steady state S_i = U / K_i either side of the
        # contact (n = 1); K1 upstream, K2 = K1/4 downstream.
        nr, nc, dx = 40, 3, 10.0
        meta = lsdtt3.RasterMeta.from_shape(nr, nc, cell_size=dx)
        rows, cols = np.indices((nr, nc), dtype=np.float32)
        z = (0.02 * dx * (nr - 1 - rows) + 0.01 * np.abs(cols - 1)).astype(np.float32)
        dem = lsdtt3.Raster.from_numpy(z, np.ones(z.shape, dtype=bool), meta)
        U, K1, dt = 1.0e-3, 1.0e-2, 1000.0
        field = np.where(rows < nr // 2, K1, K1 / 4).astype(np.float32)
        model = lsdtt3.LandscapeModel(
            dem, timestep_years=dt, erodibility=K1, m_exponent=0.0, uplift_rate=U,
            hillslope_model="none", boundary_conditions="nnon", n_threads=1,
            erodibility_field=field,
        )
        for _ in range(20000):
            if model.step().max_local_elevation_change < 1.0e-3 * U * dt:
                break
        else:
            self.fail("did not reach steady state")
        e = model.elevation.to_numpy()[:, 1].astype(np.float64)
        c = nr // 2
        s1 = (e[c - 3] - e[c - 2]) / dx
        s2 = (e[c + 2] - e[c + 3]) / dx
        self.assertAlmostEqual(s1 / (U / K1), 1.0, delta=2.0e-3)
        self.assertAlmostEqual(s2 / (U / (K1 / 4)), 1.0, delta=2.0e-3)

if __name__ == "__main__":
    unittest.main()
