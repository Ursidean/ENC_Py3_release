"""
Generates the contingency table for 2 input maps.
"""

import numpy as np


def contingency_table(map1, map2, mask, luc, rows, cols):
    cont_table = np.zeros(shape=(luc + 1, luc + 1))
    for i in range(0, rows):
        for j in range(0, cols):
            # Check to make sure the cell is not masked out
            if mask[i, j] != 0:
                # If not masked out, log the presence of the
                # cell between maps 1 & 2.
                x = map1[i, j]
                y = map2[i, j]
                cont_table[x, y] = cont_table[x, y] + 1
    # Determine the total number of each land-use class in map 1.
    for i in range(0, luc):
        for j in range(0, luc):
            x = cont_table[i, j]
            cont_table[i, luc] = cont_table[i, luc] + x
    # Determine the total number of each land-use class in map 2.
    for i in range(0, luc):
        for j in range(0, luc):
            x = cont_table[j, i]
            cont_table[luc, i] = cont_table[luc, i] + x
    # Determine the total number of cells in each map
    for i in range(0, luc):
        x = cont_table[i, luc]
        cont_table[luc, luc] = cont_table[luc, luc] + x
    return cont_table
