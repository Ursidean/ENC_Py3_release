"""
This module conducts the Mann-Whitney U-test on the specified data.
It requires familiarity with how the test is conducted. The recommended reference
is 'Nonparametric statistics A step-by-step approach,' Corder and Foreman, 2014
"""


import numpy as np


def mwu_test(max_d, luc, transition_dictionary, baseline_dictionary, N):
    mwu_z_score = np.zeros(shape=(max_d, luc, luc))
    for p in range(0, luc):
        for q in range(0, luc):
            for c in range(0, max_d):
                key_baseline = "n-" + str(q) + "|d-" + str(c)
                key_transition = "c-" + str(p) + "|n-" + str(q) + "|d-" + str(c)
                # Determine the sample sizes
                n1 = sum(transition_dictionary[key_transition])
                n2 = sum(baseline_dictionary[key_baseline])
                n = n1 + n2
                # Create a composite list for counting the presence, then rank
                # the composite list.
                count_list = [0]*(N[c] + 1)
                rank_list = [0]*(N[c] + 1)
                for a in range(0, N[c] + 1):
                    count_list[a] = (
                        count_list[a] + baseline_dictionary[key_baseline][a] +
                        transition_dictionary[key_transition][a]
                    )
                store = 0.0
                for a in range(0, N[c]+1):
                    if a == 0:
                        store = store + count_list[a]
                        rank_list[a] = (store + 1)/2
                    else:
                        rank_list[a] = (2*store + count_list[a] + 1)/2
                        store = store + count_list[a]
                # Evaluate the rank sum to calculate the U score
                R1 = 0.0
                R2 = 0.0
                for a in range(0, N[c] + 1):
                    R1 = R1 + transition_dictionary[key_transition][a]*rank_list[a]
                    R2 = R2 + baseline_dictionary[key_baseline][a]*rank_list[a]
                U1 = (n1*n2) + (0.5*n1*(n1 + 1)) - R1
                U2 = (n1*n2) + (0.5*n2*(n2 + 1)) - R2
                U = min(U1, U2)
                # Normalise the U statistic to generate a z score.
                # Calculate the mean mu.
                mu = (n1*n2)/2
                # Calculate the standard deviation sigma
                t_sum = 0.0
                for a in range(0, N[c] + 1):
                    t_sum = t_sum + ((count_list[a]**3 - count_list[a])/(n*(n - 1)))
                sigma = ((((n1*n2)/12))*((n + 1) - t_sum))**0.5
                if sigma > 0:
                    z = (U - mu)/sigma
                    mwu_z_score[c, p, q] = z
    return mwu_z_score
