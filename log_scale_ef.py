"""
Log scale the enrichment factor value numpy array to the specified base for
the active classes. The output is arranged to be consistent with the
Metronamica neighbourhood rule format of from all land-use
classes (rows) ... to active land-use classes (cols).
"""

import math
import numpy as np


def log_scale_ef(enrichment_factors, base, luc, act, pas, max_d,):
    log_ef = np.zeros(shape=(max_d, luc, act))
    for c in range(0, max_d):
        for i in range(0, act):
            for j in range(0, luc):
                if enrichment_factors[c, i + pas, j] == 0:
                    log_ef[c, j, i] = -9999
                else:
                    log_ef[c, j, i] = math.log(
                        enrichment_factors[c, i + pas, j], base
                    )
    return log_ef
