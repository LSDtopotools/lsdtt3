"""Runs a GDAL command-line program, through the approved launcher."""
import lsdtt_test_env

lsdtt_test_env.run(["gdalinfo", "-stats", "x.tif"], check=False)
