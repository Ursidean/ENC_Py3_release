"""
This module, considered distances, parses the neighbourhood of cells into
aggregated unit distances, based on the maximum distance specified.
"""


def considered_distances(max_distance):

    # The paritions of distances are tracked in a list.
    distances = []
    for i in range(1, max_distance):
        for j in range(0, max_distance):
            x = (i**2 + j**2)**0.5
            if x in distances:
                pass
            elif x >= max_distance:
                pass
            else:
                distances.append(x)
    # The maximum distance considered is added to the list of distances
    distances.append(max_distance)
    # The distances are sorted into ascending order
    distances.sort()
    # The minimum distance (0) is added to the list
    # (note: not included due to double-count error)
    distances = [0] + distances
    # The length of distances is determined, specifying the number of
    # distances to be considered.
    cdl = len(distances)
    return distances, cdl
