"""
This module calculates the enrichment factor for a given input map for transitions.
"""

import numpy as np


def ef(luc, max_d, cdl, cd, N, omap, amap, mask, row, col):
    # Initialise a dictionary to track the composition of the neighbourhood of
    # cells that transitioned between time steps
    transition_dictionary = {}
    for p in range(0, luc):
        for q in range(0, luc):
            for c in range(0, max_d):
                key = "c-" + str(p) + "|n-" + str(q) + "|d-" + str(c)
                value = [0]*(1 + N[c])
                transition_dictionary[key] = value
    # Initialise a dictionary to store the sum of enrichments for
    # neighbourhood analysis
    ef_sum = np.zeros(shape=(max_d, luc, luc))
    # Initialise lists to store the count of appeared land uses (apl) and the
    # count of each land use class.
    apl_count = [0]*luc
    luc_count = [0]*luc
    # Extract a count of the presence of different land-use classes in the
    # neighbourhood of cells that transitioned between time-slices.
    for i in range(0, row):
        for j in range(0, col):
            # Skip if masked out of region, or not included
            if mask[i,j] < 1:
                pass
            elif omap[i,j] > luc - 1:
                pass
            elif amap[i,j] > luc - 1:
                pass
            else:
                # Update the land-use count.
                luc_count[omap[i,j]] = luc_count[omap[i,j]] + 1
                # Evaluate whether a transition has occured.
                if omap[i,j] == amap[i,j]:
                    pass
                # If a transition has occured, conduct neighbourhood
                # composition analysis
                else:
                    # Initialise an array to store the count of different
                    # land-use classes in the neighbourhood of the cell of
                    # interest.
                    float_store_count = np.zeros(shape=(luc, cdl))
                    neighbourhood_size = [0]*max_d
                    for c in range(0, max_d + 1):
                        # Analyse the presence at the location of interest
                        if c == 0:
                            x = i
                            y = j
                            if mask[x, y] == 0:
                                pass
                            else:
                                float_store_count[omap[x, y],0] = (
                                    float_store_count[omap[x, y], 0] + 1
                                )
                        # Analyse the surrounding neighbourhood.
                        else:
                            x = i - c
                            if x < 0:
                                pass
                            else:
                                for y in range(j - c, j + c + 1):
                                    if y < 0:
                                        pass
                                    elif y >= col:
                                        pass
                                    elif mask[x, y] == 0:
                                        pass
                                    else:
                                        idx = ((x - i)**2 + (y - j)**2)**0.5
                                        for e in range(0, cdl):
                                            if idx == cd[e]:
                                                float_store_count[omap[x, y], e] = (
                                                    float_store_count[omap[x, y], e]
                                                    + 1
                                                )
                            x = i + c
                            if x >= row:
                                pass
                            else:
                                for y in range(j - c, j + c + 1):
                                    if y < 0:
                                        pass
                                    elif y >= col:
                                        pass
                                    elif mask[x, y] == 0:
                                        pass
                                    else:
                                        idx = ((x - i) ** 2 + (y - j) ** 2) ** 0.5
                                        for e in range(0, cdl):
                                            if idx == cd[e]:
                                                float_store_count[omap[x, y], e] = (
                                                    float_store_count[omap[x, y], e]
                                                    + 1
                                                )
                            y = j - c
                            if y < 0:
                                pass
                            else:
                                for x in range(i - c + 1, i + c):
                                    if x < 0:
                                        pass
                                    elif x >= row:
                                        pass
                                    elif mask[x, y] == 0:
                                        pass
                                    else:
                                        idx = ((x - i) ** 2 + (y - j) ** 2) ** 0.5
                                        for e in range(0, cdl):
                                            if idx == cd[e]:
                                                float_store_count[omap[x, y], e] = (
                                                    float_store_count[omap[x, y], e]
                                                    + 1
                                                )
                            y = j + c
                            if y >= col:
                                pass
                            else:
                                for x in range(i - c + 1, i + c):
                                    if x < 0:
                                        pass
                                    elif x >= row:
                                        pass
                                    elif mask[x, y] == 0:
                                        pass
                                    else:
                                        idx = ((x - i) ** 2 + (y - j) ** 2) ** 0.5
                                        for e in range(0, cdl):
                                            if idx == cd[e]:
                                                float_store_count[omap[x, y], e] = (
                                                    float_store_count[omap[x, y], e]
                                                    + 1
                                                )
                    # Initialise an array to aggregate the different distances
                    # into discrete unit rings
                    float_store_aggregated = np.zeros(shape=(luc, max_d))
                    neighbourhood_size = [0]*max_d
                    for p in range(0, luc):
                        for c in range(0, cdl):
                            if float_store_count[p, c] > 0:
                                for e in range(0, max_d):
                                    x = cd[c]
                                    if (e - 0.5) < x < (e + 0.5):
                                        float_store_aggregated[p, e] = (
                                            float_store_aggregated[p, e] +
                                            float_store_count[p, c]
                                        )
                                        neighbourhood_size[e] = (
                                            neighbourhood_size[e] +
                                            float_store_count[p, c]
                                        )
                    # Convert the float stored to a proportion for processing.
                    fsa_proportion = np.zeros(shape=(luc, max_d))
                    for p in range(0, luc):
                        for c in range(0, max_d):
                            if float_store_aggregated[p,c] == 0:
                                pass
                            else:
                                fsa_proportion[p,c] = (
                                    float_store_aggregated[p,c]/
                                    neighbourhood_size[c]
                                )
                    # Summate the values for calculation of the enrichment factor.
                    central = amap[i,j]
                    apl_count[central] = apl_count[central] + 1
                    for p in range(0, luc):
                        for c in range(0, max_d):
                            ef_sum[c, central, p] = (
                                ef_sum[c, central, p] + fsa_proportion[p, c]
                            )
                            key = (
                                "c-" + str(central) + "|n-" + str(p) + "|d-" +
                                str(c)
                            )
                            for a in range(0, N[c] + 1):
                                lower_bound = float(a)/float(N[c]) - (0.5/N[c])
                                upper_bound = float(a)/float(N[c]) + (0.5/N[c])
                                if lower_bound < fsa_proportion[p, c] < upper_bound:
                                    transition_dictionary[key][a] = (
                                        transition_dictionary[key][a] + 1
                                    )
    # Calculate the enrichment factor values per unit distance.
    enrichment_factors = np.zeros(shape=(max_d, luc, luc))
    for p in range(0, luc):
        for q in range(0, luc):
            for c in range(0, max_d):
                if apl_count[p] > 0 and luc_count[q] > 0:
                    x = ef_sum[c, p, q]/apl_count[p]
                    y = float(luc_count[q])/sum(luc_count[:])
                    enrichment_factors[c, p, q] = (x/y)
    return enrichment_factors
