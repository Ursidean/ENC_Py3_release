"""
This module calculates the clumpiness metric (http://www.umass.edu/landeco/
research/fragstats/documents/Metrics/Contagion%20-%20Interspersion%20Metrics/
Metrics/C115%20-%20CLUMPY.htm) of each land-use class for a given input map.
For more information on this metric, download the free fragstats package:
www.umass.edu/landeco/research/fragstats/fragstats.html
"""


import numpy as np


def clumpiness_index(map1, mask, luc):
    # Determine the dimensions of the input map.
    mapshape = np.shape(map1)
    row = mapshape[0]
    column = mapshape[1]
    # Populate the adjacency matrix, used to calculate the clumpiness.
    adjacency_table = np.zeros(shape=(luc, luc))
    count_luc = [0]*luc
    for i in range(0, row):
        for j in range(0, column):
            if mask[i, j] > 0:
                central = map1[i, j]
                count_luc[central] = count_luc[central] + 1
                k = i - 1
                if k < 0:
                    pass
                elif mask[k, j] > 0:
                    link = map1[k, j]
                    adjacency_table[central, link] = (
                        adjacency_table[central, link] + 1
                    )
                l = i + 1
                if l >= row:
                    pass
                elif mask[l, j] > 0:
                    link = map1[l, j]
                    adjacency_table[central, link] = (
                        adjacency_table[central, link] + 1
                    )
                m = j - 1
                if m < 0:
                    pass
                elif mask[i, m] > 0:
                    link = map1[i, m]
                    adjacency_table[central, link] = (
                        adjacency_table[central, link] + 1
                    )
                n = j + 1
                if n >= column:
                    pass
                elif mask[i, n] > 0:
                    link = map1[i, n]
                    adjacency_table[central, link] = (
                        adjacency_table[central, link] + 1
                    )
    # Calculate the composition and proportion of different land-use classes
    # in the map
    prop = [0]*luc
    for i in range(0, luc):
        prop[i] = prop[i] + float(count_luc[i])/sum(count_luc)
    # Determine the background and edge segment arrays for calculating
    # clumpiness.
    background = [0]*luc
    edge_segment_matrix = np.zeros(shape=(row, column))
    for i in range(0, row):
        for j in range(0, column):
            if mask[i, j] == 0:
                pass
            else:
                if i-1 < 0:
                    edge_segment_matrix[i, j] = edge_segment_matrix[i, j] + 1
                else:
                    k = mask[i - 1, j]
                    if k == 0:
                        edge_segment_matrix[i, j] = (
                            edge_segment_matrix[i, j] + 1
                        )
                if i + 1 >= row:
                    edge_segment_matrix[i, j] = edge_segment_matrix[i, j] + 1
                else:
                    l = mask[i + 1, j]
                    if l == 0:
                        edge_segment_matrix[i, j] = (
                            edge_segment_matrix[i, j] + 1
                        )
                if j - 1 < 0:
                    edge_segment_matrix[i, j] = edge_segment_matrix[i, j] + 1
                else:
                    m = mask[i, j - 1]
                    if m == 0:
                        edge_segment_matrix[i, j] = (
                            edge_segment_matrix[i, j] + 1
                        )
                if j + 1 >= column:
                    edge_segment_matrix[i, j] = edge_segment_matrix[i, j] + 1
                else:
                    n = mask[i, j + 1]
                    if n == 0:
                        edge_segment_matrix[i, j] = (
                            edge_segment_matrix[i, j] + 1
                        )
    for i in range(0, row):
        for j in range(0, column):
            x = map1[i, j]
            if x != -9999:
                background[x] = background[x] + edge_segment_matrix[i, j]
    # Calculate the clumpiness metric.
    clump = np.zeros(shape=(luc, 8))
    for i in range(0, luc):
        clump[i, 1] = clump[i, 1] + background[i]
        clump[i, 2] = count_luc[i]
        clump[i, 3] = (clump[i, 2])**0.5
        clump[i, 4] = int(clump[i, 3])
        clump[i, 5] = clump[i, 2] - (clump[i, 4])**2
        if clump[i, 5] == 0:
            clump[i, 6] = 4*clump[i, 4]
        elif (clump[i, 4])**2 <= clump[i, 2] <= clump[i, 4]*(1 + clump[i, 4]):
            clump[i, 6] = 4*clump[i, 4] + 2
        else:
            clump[i, 6] = 4*clump[i, 4] + 4
        for j in range(0, luc):
            if i == j:
                clump[i, 0] = adjacency_table[i, j]
            clump[i, 1] = clump[i, 1] + adjacency_table[i, j]
    for i in range(0, luc):
        if (clump[i, 1] - clump[i, 6]) != 0:
            clump[i, 7] = clump[i, 0]/(clump[i, 1] - clump[i, 6])
    
    clumpy = [0]*luc
    for i in range(0, luc):
        if clump[i, 7] < prop[i] and prop[i] < 0.5:
            clumpy[i] = (clump[i, 7] - prop[i])/prop[i]
        else:
            clumpy[i] = (clump[i, 7] - prop[i])/(1 - prop[i])
            
    # Return a list of the clumpiness values per land-use class.
    return clumpy
