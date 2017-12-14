"""
read_map is used to read in maps using the osgeo - gdal module, which requires 
installation. This allows python to interpret GIS file format data such as 
raster or ascii, and store these as numpy arrays
"""

from osgeo import gdal
import sys
import numpy as np

def read_map(map_path):
    # 1. Open source file path
    try:
        src_map = gdal.Open(map_path)
    except RuntimeError as e:
        print("Unable to read map file")
        print(e)
        sys.exit(1)

    # 2. Convert map to array
    try:
        array_map = np.array(src_map.GetRasterBand(1).ReadAsArray())
    except RuntimeError as e:
        # Example, try GetRasterBand(10)
        print("Band (%i) not found' %band_num")
        print(e)
        sys.exit(1)

    return array_map
