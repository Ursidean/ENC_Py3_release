"""
Stage 1 of the Empirical Neighbourhood Calibration method:
Conduct the significance test on the specified data and eliminate
interactions.
"""

# Modules
# Allows for multi-dimensional array handling.
import numpy as np
# Allows for maps to be stored as 2d arrays.
from read_map import read_map
# Determines the size (# of cells) of square neighbourhoods.
from considered_distances import considered_distances
# Evaluates the composition of different neighbourhoods.
from neighbourhood_evaluator import neighbourhood_evaluator
# Conducts the Mann-Whitney U-test for input neighbourhood composition
# dictionaries.
from MWU_test import mwu_test
# Conducts complex mathematical operations.
import math
# Generate the contingency table.
from contingency_table import contingency_table

# Specify the base path to the directory containing the empirical neighbourhood
# calibration tool-pack.
base_path = "C:\\Users\\charl\OneDrive\\Documents\\ENC_Py3_release\\"
# Set the case study
case_study = "Berlin"
# Set the paths to the directories and relevant data
data_path = base_path + "EU_data\\" + case_study + "\\"
output_path = base_path + "EU_output\\" + case_study + "\\"
# Specify the original map (data at time slice 1).
omap_path = data_path + case_study.lower() + "_1990.asc"
# Specify the actual map (data at time slice 2).
amap_path = data_path + case_study.lower() + "_2000.asc"
# Specify the masking map.
mask_path = data_path + case_study.lower() + "_mask.asc"
# Set the land-use class names.
luc_names = ["Natural areas", "Arable land", "Permanent crops", "Pastures",
             "Agricultural areas", "Residential", "Industry & commerce",
             "Recreation areas", "Forest", "Road & rail", "Seaports",
             "Airports", "Mine & dump sites", "Fresh water", "Marine water"]
# Set the land-use class parameters: number of land-use classes, passive,
# feature, and active.
luc = len(luc_names)
pas = 1
fea = 6
act = luc - (pas + fea)
# Specify the maximum neighbourhood size distance considered
max_distance = 5

# Set the minimum required percentage of conversions.
min_convo_rate = 0.025
# Set the minimum required Enrichment Factor value at distance 1.
min_EF_1 = 0.00
# Set the significance limit (recommended value is 1.96
z_limit = 1.96

# Read in the map for time slice 1.
omap = read_map(omap_path)
# Read in the map for time slice 2.
amap = read_map(amap_path)
# Read in the masking map.
mask = read_map(mask_path)
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

# Data analysis -----------------------------------------------------
# Generate a contingency table for the data.
cont_table = contingency_table(omap, amap, mask, luc, rows, cols)

# Evaluate the composition of the neighbourhoods of cells in the data.
dummy = neighbourhood_evaluator(
    luc, max_distance, cdl, cd, N, omap, amap, mask, rows, cols
)
# Store the requisite output into specific dictionaries
all_cells_baseline = dummy[0]
no_new_cells_ci_baseline = dummy[1]
no_cells_ci_baseline = dummy[2]
transition_dictionary = dummy[3]
ef = dummy[4]
# Log scale the enrichment factor values.
log_ef = np.zeros(shape=(max_distance, luc, luc))
for p in range(0, luc):
    for q in range(0, luc):
        for c in range(0, max_distance):
            # If the enrichment factor is not evaluated a value of 0 is given.
            # Hence, a did not error value of -9999 is used.
            if ef[c, p, q] == 0:
                log_ef[c, p, q] = -9999
            else:
                log_ef[c, p, q] = math.log(ef[c, p, q], 10)

# Conduct the MWU test using the mwu_test module.
z_scores = mwu_test(
    max_distance, luc, transition_dictionary, all_cells_baseline, N
)
# Conversion point elimination --------------------------------------
# Initialise a matrix to store the conversion points.
conversion_points = np.zeros(shape = (luc, act))
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
# Eliminate the conversions, set the included conversions at 1.
c = 0
for i in range(0, act):
    for j in range(0, luc):
        conversion_rate = ic_rates[j, i + pas]
        if i + pas == j:
            conversion_points[j, i] = 1
        elif conversion_rate > min_convo_rate:
            if abs(z_scores[c, i + pas, j]) > z_limit:
                conversion_points[j, i] = 1
# Interaction tail elimination --------------------------------------
# Initialise an array to track the distances evaluated as significant.
sig_distances = np.zeros(shape=(act, luc))
# Conduct the analysis
for c in range(1, max_distance):
    for i in range(0, act):
        for j in range(0, luc):
            if (
                abs(z_scores[c, i + pas, j]) > z_limit and
                log_ef[c, i + pas, j] > (min_EF_1 / 1)
            ):
                sig_distances[i, j] = sig_distances[i, j] + 1
            if(
                abs(z_scores[c, i + pas, j]) > z_limit and
                log_ef[c, i + pas, j] < (-1*min_EF_1 / 1)
            ):
                sig_distances[i, j] = sig_distances[i, j] - 1
# Initialise an array to track rules that are included.
int_rules = np.zeros(shape=(act, luc))
# Analyse the results
for i in range(0, act):
    for j in range(0, luc):
        if i + pas == j:
            int_rules[i, j] = 1
        elif sig_distances[i, j] == (max_distance - 1):
            int_rules[i, j] = 1

# Transpose the output for logging.
int_rules = np.transpose(int_rules)
# Save the output.
tails_file = output_path + "Rules\\tails.txt"
points_file = output_path + "Rules\\points.txt"
np.savetxt(tails_file, int_rules, fmt="%d")
np.savetxt(points_file, conversion_points, fmt="%d")

# Completed!
