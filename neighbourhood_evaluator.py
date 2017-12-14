""" This module, neighbourhood_evaluator, evaluates neighbourhood composition
for each land-use class in the data
"""
# The following modules are used.
import numpy as np


# Function: neighbourhood_evaluator
def neighbourhood_evaluator(luc, max_d, cdl, cd, N, omap, amap, mask, row, col):
    # Code block: Initialisation
    # Initialise a dictionary to log the neighbourhood land-use class (n)
    # presence for cells that transitioned to a particular class (c) at
    # each specific distance (d).
    transition_dictionary = {}
    for c in range(0, luc):
        for n in range(0, luc):
            for d in range(0, max_d):
                key = "c-" + str(c) + "|n-" + str(n) + "|d-" + str(d)
                value = [0]*(1 + N[d])
                transition_dictionary[key] = value
    # Initialise a set of dictionaries to log the neighbourhood land-use
    # class (n) presence for cells at each specific distance (d) for
    # three possible baselines:
    # 1. All cells considered;
    # 2. All cells except the cells of a particular class of interest (coi); &
    # 3. All cells except the newly allocated cells of a class of interest.
    all_cells_baseline = {}
    no_cells_ci_baseline = {}
    no_new_cells_ci_baseline = {}
    for n in range(0, luc):
        for d in range(0, max_d):
            key = "n-" + str(n) + "|d-" + str(d)
            value = [0]*(1 + N[d])
            all_cells_baseline[key] = value
            no_cells_ci_baseline[key] = value
            no_new_cells_ci_baseline[key] = value
    # Initialise an array to store the sum of enrichment
    EF_sum = np.zeros(shape=(max_d,luc,luc))
    # Initialise a list to store the count of appeared land-use classes (apl).
    apl_count = [0]*luc
    # Initialise a list to store the count of each land-use count.
    luc_count = [0]*luc

    # Code block: Map evaluation
    # Extract a count of the presence of different land-use classes
    # in neighbourhoods
    for i in range(0, row):
        for j in range(0, col):
            # Skip if masked out of region map.
            if mask[i, j] < 1:
                pass
            # Skip if land-use class not considered (avoids potential errors)
            elif omap[i, j] > (luc - 1):
                pass
            elif amap[i, j] > (luc - 1):
                pass
            else:
                # Add one to the land-use counter taken, from the original data.
                luc_count[omap[i,j]] = luc_count[omap[i,j]] + 1
                # Use an array to store the count of different land-use classes
                # in the neighbourhood of the cell of interest.
                float_store_count = np.zeros(shape=(luc,cdl))
                neighbourhood_size = [0]*max_d
                for d in range(0, (max_d + 1)):
                    # Analyse the presence at the location of interest.
                    if d == 0:
                        x = i
                        y = j
                        if mask[x, y] == 0:
                            pass
                        else:
                            float_store_count[omap[x, y], 0] = (
                                float_store_count[omap[x, y], 0] + 1
                            )
                    # Iterate around the neighbourhood of the cell of interest
                    else:
                        x = i - d
                        if x < 0:
                            pass
                        else:
                            for y in range((j - d), (j + d + 1)):
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
                                                float_store_count[omap[x, y], e] + 1
                                            )
                        x = i + d
                        if x >= row:
                            pass
                        else:
                            for y in range((j - d), (j + d + 1)):
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
                                                float_store_count[omap[x, y], e] + 1
                                            )
                        y = j - d
                        if y < 0:
                            pass
                        else:
                            for x in range((i - d + 1), (i + d)):
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
                                                float_store_count[omap[x, y], e] + 1
                                            )
                        y = j + d
                        if y >= col:
                            pass
                        else:
                            for x in range((i - d + 1), (i + d)):
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
                                                float_store_count[omap[x, y], e] + 1
                                            )
                # Initialise an array to aggregate the composition values into
                # discrete unit distances (1,2,3) by half distances
                float_store_aggregated = np.zeros(shape=(luc, max_d))
                neighbourhood_size = [0]*max_d
                for p in range(0, luc):
                    for c in range(0, cdl):
                        if float_store_count[p, c] > 0:
                            for e in range(0, max_d):
                                k = cd[c]
                                if (e - 0.5) < k < (e + 0.5):
                                    float_store_aggregated[p, e] = (
                                        float_store_aggregated[p, e] +
                                        float_store_count[p, c]
                                    )
                                    neighbourhood_size[e] = (
                                        neighbourhood_size[e] +
                                        float_store_count[p, c]
                                    )
                # Convert the aggregated values to proportional values
                fsa_proportion = np.zeros(shape=(luc, max_d))
                for p in range(0, luc):
                    for c in range(0, max_d):
                        if float_store_aggregated[p, c] == 0:
                            pass
                        else:
                            fsa_proportion[p, c] = (
                                float_store_aggregated[p, c] /
                                neighbourhood_size[c]
                            )
                # Store values into requisite baseline dictionary bins
                for p in range(0, luc):
                    for c in range(0, max_d):
                        key = "n-" + str(p) + "|d-" + str(c)
                        for a in range(0, (N[c] + 1)):
                            lower_bound = float(a)/float(N[c]) - (0.5/N[c])
                            upper_bound = float(a)/float(N[c]) + (0.5/N[c])
                            if lower_bound < fsa_proportion[p, c] < upper_bound:
                                all_cells_baseline[key][a] = (
                                    all_cells_baseline[key][a] + 1
                                )
                                if omap[i, j] != amap[i, j] == p:
                                    pass
                                else:
                                    no_new_cells_ci_baseline[key][a] = (
                                        no_new_cells_ci_baseline[key][a] + 1
                                    )
                                if omap[i, j] == p:
                                    pass
                                else:
                                    no_cells_ci_baseline[key][a] = (
                                        no_cells_ci_baseline[key][a] + 1
                                    )
                # Store values into the transition dictionary bins
                if omap[i, j] != amap[i, j]:
                    central = amap[i, j]
                    apl_count[central] = apl_count[central] + 1
                    for p in range(0, luc):
                        for c in range(0, max_d):
                            EF_sum[c, central, p] = (
                                EF_sum[c, central, p] + fsa_proportion[p, c]
                            )
                            key = (
                                "c-" + str(central) + "|n-" + str(p) + "|d-" + str(c)
                            )
                            for a in range(0, (N[c] + 1)):
                                lower_bound = float(a)/float(N[c]) - (0.5/N[c])
                                upper_bound = float(a)/float(N[c]) + (0.5/N[c])
                                if lower_bound < fsa_proportion[p, c] < upper_bound:
                                    transition_dictionary[key][a] = (
                                        transition_dictionary[key][a] + 1
                                    )
    # Process the enrichment factor calculation
    EF = np.zeros(shape=(max_d, luc, luc))
    for p in range(0, luc):
        for q in range(0, luc):
            for c in range(0, max_d):
                if apl_count[p] > 0 and luc_count[q] > 0:
                    x = EF_sum[c, p, q]/apl_count[p]
                    y=float(luc_count[q])/sum(luc_count[:])
                    EF[c, p, q] = x/y
# Return the values stored in the dictionaries.
    return (
        all_cells_baseline, no_new_cells_ci_baseline, no_cells_ci_baseline,
        transition_dictionary, EF
    )
