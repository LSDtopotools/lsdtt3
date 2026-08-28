# Load, inspect, and save a DEM

This is the smallest useful raster workflow: read a DEM, condition it, derive
a display product, and explicitly write results that you want to keep.

```python
import numpy as np
import lsdtt3

dem = lsdtt3.Raster.read("dem.tif")
filled = lsdtt3.fill_pits(dem)
hillshade = lsdtt3.multidirectional_hillshade(filled)

values = np.where(hillshade.mask_numpy(), hillshade.to_numpy(), np.nan)
print(filled.shape, filled.metadata.cell_size)

filled.write("dem_filled.tif")
hillshade.write("dem_hillshade.tif")
```

`Raster.to_numpy()` returns the numerical values. Pair it with
`Raster.mask_numpy()` before plotting or analysis when nodata must become
`numpy.nan`.

For additional terrain products, see [Raster and terrain API](../python_api/raster.md).
