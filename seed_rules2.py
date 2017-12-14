"""
This module is used to seed the specified Metronamica model based on a set
of input meta-parameter values with neighbourhood rules.
"""

# Allows for multi-dimensional array handling.
import numpy as np
# Calculate the enrichment factor.
from enrichment_factor import ef
# Log scale the enrichment factor values (default is to base 10).
from log_scale_ef import log_scale_ef
# Generate a contingency table for the data.
from contingency_table import contingency_table
# Determines the size (# of cells) of square neighbourhoods.
from considered_distances import considered_distances
# Set neighbourhood rules based on a four-point structure.
from set_NR import set_lp_rule


def seed_rules2(omap, amap, mask, max_distance, luc_names, luc, act, pas,
                int_rules, theta_it, theta_cp, theta_ct, project_file):
    # Specify the bands levels (Change as neccessary).
    # The inertia band levels, for setting inertia values.
    high_inertia_band = 0.95
    mid_inertia_band = 0.90
    # The conversion band levels, for setting conversion values.
    high_conversion_band = 0.5
    mid_conversion_band = 0.1
    low_conversion_band = 0.025
    # The enrichment factor band levels, for setting the tail values.
    high_ef = 1.0
    mid_ef = 0.5
    low_ef = 0.0
    # Specify high, medium and low inertia values.
    # The default settings are a high inertia of 1000, med 500, low 250.
    high_inertia = 1000.0
    mid_inertia = 500.0
    low_inertia = 250.0

    # Analyse the input maps for evaluation purposes
    map_dimensions = np.shape(omap)
    rows = map_dimensions[0]
    cols = map_dimensions[1]            
    # Determine the distances that will be analysed using the module considered
    # distances.
    temp = considered_distances(max_distance)
    # Store the list of considered distances as a variable.
    cd = temp[0]
    # Store the total number of distances considered
    cdl = temp[1]
    # Determine the maximum neighbourhood size (unit) from considered distances
    N_all = [1, 8, 12, 16, 32, 28, 40, 40, 20]
    N = []
    for c in range(0, max_distance):
        N.append(N_all[c])
        
    # Generate the Enrichment Factor and contingency table.
    data_ef = ef(luc, max_distance, cdl, cd, N, omap, amap, mask, rows, cols)
    log_data_ef = log_scale_ef(data_ef, 10, luc, act, pas, max_distance)
    cont_table = contingency_table(omap, amap, mask, luc, rows, cols)
    
    # Determine the rates of inertia and conversion between different land-uses.
    ic_rates = np.zeros(shape=(luc, luc))
    for i in range(0, luc):
        for j in range(0, luc):
            if i == j:
                if cont_table[i, luc] > 0:
                    ic_rates[i, j] = cont_table[i, j] / cont_table[i, luc]
            else:
                conversions = abs(float(cont_table[j, j]) - float(cont_table[luc, j]))
                if conversions > 0:
                    ic_rates[i, j] = float(cont_table[i, j]) / float(conversions)
    # Calculate the different influence values.
    # Set influence values at distance 1 for self-influence tails.
    d1_high_it_value = high_inertia * theta_it
    d1_mid_it_value = mid_inertia * theta_it
    d1_low_it_value = low_inertia * theta_it
    # Set influence values at distance 2 for self-influence tails.
    d2_high_it_value = high_inertia * theta_it * 0.1
    d2_mid_it_value = mid_inertia * theta_it * 0.1
    d2_low_it_value = low_inertia * theta_it * 0.1
    # Set the conversion parameter values.
    high_conversion = high_inertia * theta_cp
    mid_conversion = mid_inertia * theta_cp
    low_conversion = low_inertia * theta_cp
    # Set influence values at distance 1 for interaction tails.
    d1_high_ct_value = high_inertia * theta_ct
    d1_mid_ct_value = mid_inertia * theta_ct
    d1_low_ct_value = low_inertia * theta_ct
    # Set influence values at distance 2 for interaction tails.
    d2_high_ct_value = high_inertia * theta_ct * 0.1
    d2_mid_ct_value = mid_inertia * theta_ct * 0.1
    d2_low_ct_value = low_inertia * theta_ct * 0.1
    # Initialise a dictionary to track the rules.
    rules = {}
    for i in range(0, act):
        for j in range(0, luc):
            key = "from " + luc_names[j] + " to " + luc_names[i + pas]
            rules[key] = [0, 0, 0, 5]

    # Set the values for inertia and conversion.
    for i in range(0, act):
        for j in range(0, luc):
            # Specify the neighbourhood rule key.
            key = (
                "from " + luc_names[j] + " to " + luc_names[i + pas]
            )
            # If a self-influence rule, set the inertia value.
            if i + pas == j:
                if (cont_table[i + pas, luc] >
                        cont_table[luc, i + pas]):
                    rules[key][0] = low_inertia
                else:
                    inertia_rate = ic_rates[j, i + pas]
                    if inertia_rate > high_inertia_band:
                        rules[key][0] = high_inertia
                    elif inertia_rate > mid_inertia_band:
                        rules[key][0] = mid_inertia
                    else:
                        rules[key][0] = low_inertia
            # If an interactive rule, set the conversion rule.
            else:
                conversion_rate = ic_rates[j, i + pas]
                if conversion_rate > high_conversion_band:
                    rules[key][0] = high_conversion
                elif conversion_rate > mid_conversion_band:
                    rules[key][0] = mid_conversion
                elif conversion_rate > low_conversion_band:
                    rules[key][0] = low_conversion
    # Set the values for self-influence and conversion tails.
    for i in range(0, act):
        for j in range(0, luc):
            # Specify the neighbourhood rule key.
            key = "from " + luc_names[j] + " to " + luc_names[i + pas]
            # If a self-influence rule, set the self-influence attraction values.
            if i + pas == j:
                if int_rules[j, i] == 1:
                    if log_data_ef[1, j, i] > high_ef:
                        rules[key][1] = d1_high_it_value
                        rules[key][2] = d2_high_it_value
                    elif log_data_ef[1, j, i] > mid_ef:
                        rules[key][1] = d1_mid_it_value
                        rules[key][2] = d2_mid_it_value
                    else:
                        rules[key][1] = d1_low_it_value
                        rules[key][2] = d2_low_it_value
            # If a conversion rule, set the interaction values.
            else:
                if int_rules[j, i] == 1:
                    if log_data_ef[1, j, i] > high_ef:
                        rules[key][1] = d1_high_ct_value
                        rules[key][2] = d2_high_ct_value
                    elif log_data_ef[1, j, i] > mid_ef:
                        rules[key][1] = d1_mid_ct_value
                        rules[key][2] = d2_mid_ct_value
                    elif log_data_ef[1, j, i] > low_ef:
                        rules[key][1] = d1_low_ct_value
                        rules[key][2] = d2_low_ct_value
    # Return the rules.
    return rules
