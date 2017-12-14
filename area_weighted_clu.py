"""
This module is used to calculate the area weighted clumpiness error for a
given input map. It is developed for testing purposes.
"""

from read_map import read_map
import numpy as np
from clumpy_module import clumpiness_index


def area_weighted_clu_error(map1, map2, mask, luc, pas, act, luc_count):
    # Calculate the clumpiness.
    map1_clumpiness = clumpiness_index(map1, mask, luc)
    map2_clumpiness = clumpiness_index(map2, mask, luc)
    clu_error = [0] * luc
    for i in range(0, luc):
        clu_error[i] = abs(map1_clumpiness[i] - map2_clumpiness[i])
    # Extract the clumpiness error of the active classes.
    act_clu_error = [0] * (act)
    for i in range(0, act):
        act_clu_error[i] = clu_error[i + pas]
    # Calculate the area-weighted clumpiness error.
    AWCE = 0
    for i in range(0, act):
        AWCE = AWCE + act_clu_error[i] * luc_count[i + pas]
    AWCE = AWCE / (sum(luc_count))
    return  AWCE

# Module test.
"""
smap_path = ("C:\\Neighbourhood_sampling\\master\\Randstad_simpler\\"
             "Randstad\\Log\\Land_use\\sim_1000\\Land use map_2000.rst")
amap_path = ("C:\\Neighbourhood_sampling\\master\\Randstad_simpler\\"
             "Data\\lu2000.asc")
mask_path = ("C:\\Neighbourhood_sampling\\master\\Randstad_simpler\\"
             "Data\\region.asc")
# Read in the maps
smap = read_map(smap_path)
amap = read_map(amap_path)
mask = read_map(mask_path)
# Determine the map properties
map_dimensions = np.shape(amap)
rows = map_dimensions[0]
cols = map_dimensions[1]
luc = np.max(amap) + 1
pas = 1
fea = 2
act = luc - (pas + fea)

AWCE = area_weighted_clu_error(amap, smap, mask, rows, cols, luc, pas, act)
print AWCE
"""