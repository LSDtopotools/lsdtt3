"""Small offline example used by installed-package documentation checks."""

import numpy as np
import lsdtt3

meta = lsdtt3.RasterMeta.from_shape(3, 3, cell_size=1.0)
values = np.arange(9, dtype=np.float32).reshape(3, 3)
raster = lsdtt3.Raster.from_numpy(values, np.ones((3, 3), dtype=bool), meta)
np.testing.assert_array_equal(lsdtt3.scale(raster, 2.0).to_numpy(), values * 2.0)
