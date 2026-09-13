import unittest

import numpy as np

import lsdtt3


@unittest.skipUnless(lsdtt3.spectral_available(), "requires an FFTW-enabled build")
class SpectralTest(unittest.TestCase):
    def raster(self):
        metadata = lsdtt3.RasterMeta.from_shape(64, 64, cell_size=10.0)
        rows, cols = np.indices((64, 64), dtype=np.float32)
        data = np.sin(2.0 * np.pi * cols / 16.0) + 0.2 * np.sin(
            2.0 * np.pi * rows / 32.0
        )
        return lsdtt3.Raster.from_numpy(data, np.ones_like(data, dtype=bool), metadata)

    def test_filters_preserve_grid_and_mask(self):
        raster = self.raster()
        outputs = [
            lsdtt3.lowpass_filter(raster, 100.0, n_threads=1),
            lsdtt3.highpass_filter(raster, 100.0, n_threads=1),
            lsdtt3.bandpass_filter(raster, 20.0, 100.0, n_threads=1),
            lsdtt3.wiener_filter(
                raster, signal_wavelength_long_metres=300.0,
                signal_wavelength_short_metres=30.0, n_threads=1,
            ),
        ]
        self.assertTrue(all(output.shape == raster.shape for output in outputs))
        self.assertTrue(all(output.valid_pixel_count == raster.valid_pixel_count for output in outputs))

    def test_full_analysis_and_directional_periodicity(self):
        raster = self.raster()
        result = lsdtt3.spectral_analysis(
            raster, n_background_iterations=2, seed=7, n_threads=1
        )

        self.assertEqual(result.padded_shape, (64, 64))
        self.assertEqual(result.psd_2d.shape, result.padded_shape)
        self.assertEqual(result.background_psd_2d.shape, result.padded_shape)
        self.assertEqual(result.frequency_x.shape, (64,))
        self.assertEqual(result.frequency_y.shape, (64,))
        self.assertEqual(result.frequency.shape, result.psd.shape)
        self.assertTrue(np.isfinite(result.rollover_frequency))
        self.assertEqual(
            set(result.binned),
            {"frequency", "wavelength", "psd", "background_psd", "normalised_psd", "beta", "r_squared"},
        )

        periodicity = lsdtt3.detect_directional_periodicity(
            raster, 16, 8, 25.0, 100.0, n_threads=1
        )
        self.assertEqual(
            set(periodicity), {"strength", "spacing_metres", "orientation_degrees"}
        )
        self.assertTrue(all(output.shape == raster.shape for output in periodicity.values()))

    def test_masked_surface_recovery(self):
        raster = self.raster()
        mask = np.zeros(raster.shape, dtype=np.float32)
        mask[:, :32] = 1.0
        mask_raster = lsdtt3.Raster.from_numpy(
            mask, np.ones(raster.shape, dtype=bool), raster.metadata
        )
        recovered = lsdtt3.masked_surface_recovery(
            raster, mask_raster, 100.0, n_threads=1
        )
        self.assertEqual(recovered.shape, raster.shape)
        self.assertEqual(recovered.valid_pixel_count, raster.valid_pixel_count)


if __name__ == "__main__":
    unittest.main()
