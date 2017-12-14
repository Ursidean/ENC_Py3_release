"""
Stage 3 of the Empirical Neighbourhood Calibration method:
Conduct the coarse adjustment stage via structured sampling of meta-
parameter values.
"""

# Allows for multi-dimensional array handling.
import numpy as np
# Allows for maps to be stored as 2d arrays.
from read_map import read_map
# Determines the size (# of cells) of square neighbourhoods.
from considered_distances import considered_distances
# Calculate the enrichment factor.
from enrichment_factor import ef
# Log scale the enrichment factor values (default is to base 10).
from log_scale_ef import log_scale_ef
# Generate a contingency table for the data.
from contingency_table import contingency_table
# Set neighbourhood rules based on a four-point structure.
from set_NR import set_lp_rule
# Set the random number seed in the Metronamica file.
from set_rand import set_rand
# Run the Metronamica model.
from run_metro import run_metro
# Module for calculation of Fuzzy Kappa.
from fuzzy_kappa import fuzzy_kappa
# Module for calculation of Fuzzy Kappa Simulation.
from fuzzy_kappa import fks
# Calculate the absolute area-weighted clumpiness error (AAWCE).
from area_weighted_clu import area_weighted_clu_error
# Interact with csv files.
import csv

# Specify the base path to the directory containing the empirical neighbourhood
# calibration tool-pack.
base_path = "C:\\Users\\charl\OneDrive\\Documents\\ENC_Py3_release\\"
# Set the case study
case_study = "Berlin"
# Set the paths to the directories and relevant data
data_path = base_path + "EU_data\\" + case_study + "\\"
output_path = base_path + "EU_output\\" + case_study + "\\Coarse_cal_output\\"
# Specify the original map (data at time slice 1).
omap_path = data_path + case_study.lower() + "_1990.asc"
# Specify the actual map (data at time slice 2).
amap_path = data_path + case_study.lower() + "_2000.asc"
# Specify the masking map.
mask_path = data_path + case_study.lower() + "_mask.asc"
# Specify the fuzzy coefficients.
# Specify the fuzzy weights for the calculation of fuzzy Kappa.
fuzzy_coefficients = base_path + "EU_data\\coeff13.txt"
# Specify the fuzzy transition weights for the calculation of FKS.
fuzzy_trans_coefficients = base_path + "EU_data\\coefficients13.txt"

# Set the working directory, which contains the geoproject file.
working_directory = ("C:\\Geonamica\\Metronamica\\" + case_study)
# Set the project file path.
project_file = working_directory + "\\" + case_study + ".geoproj"
# Set the path to the command line version of Geonamica
geo_cmd = "C:\\Program Files (x86)\\Geonamica\\Metronamica\\GeonamicaCmd.exe"
# Set the path to the log file.
log_file = base_path + "LogSettings.xml"
# Set the path to the simulated output map
smap_path = (
    working_directory + "\\Log\\Land_use\\"
                        "Land use map_2000-Jan-01 00_00_00.rst"
)

# Set the land-use class names.
luc_names = ["Natural areas", "Arable land", "Permanent crops", "Pastures",
             "Agricultural areas", "Residential", "Industry & commerce",
             "Recreation areas", "Forest", "Road & rail", "Port area",
             "Airports", "Mine & dump sites", "Fresh water", "Marine water"]
# Set the land-use class parameters: number of land-use classes, passive,
# feature, and active.
luc = len(luc_names)
pas = 1
fea = 6
act = luc - (pas + fea)
# Specify the maximum neighbourhood size distance considered
max_distance = 5

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

# Count the presence of each land-use class in the actual map. This is
# used in the calculation of area-weighted average clumpiness across the
# active classes.
luc_count = [0] * luc
for i in range(0, rows):
    for j in range(0, cols):
        if mask[i, j] > 0:
            luc_count[amap[i, j]] = luc_count[amap[i, j]] + 1

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

# Specify the point levels file
point_levels_file_path = (base_path + "EU_output\\" + case_study +
                          "\\Rules\\point_levels.txt")
# Specify the tail levels file path
tail_levels_file_path = (base_path + "EU_output\\" + case_study +
                          "\\Rules\\tail_levels.txt")
# Load the specified files.
point_levels = np.loadtxt(point_levels_file_path)
tail_levels = np.loadtxt(tail_levels_file_path)

# Coarse adjustment settings ----------------------------------------
# Specify high, medium and low inertia values.
# The default settings are a high inertia of 1000, med 500, low 250.
high_inertia_point = 1000.0
mid_inertia_point = 500.0
low_inertia_point = 250.0

# Set the varied parameter. Must be one of theta_it, theta_cp or theta_ct
vp = "theta_ct"
# Set the fixed parameters.
if vp == "theta_it":
    # Values set by user
    theta_cp = 0.025
    theta_ct = 0.005
    # Set the range and interval size for the selected parameter.
    min_value = 0.000
    max_value = 0.100
    interval_size = 0.005
elif vp == "theta_cp":
    # Values set by user.
    theta_it = 0.050
    theta_ct = 0.005
    # Set the range and interval size for the selected parameter.
    min_value = 0.000
    max_value = 0.050
    interval_size = 0.0025
elif vp == "theta_ct":
    # Values set by user.
    theta_it = 0.050
    theta_cp = 0.025
    # Set the range and interval size for the selected parameter.
    min_value = 0.000
    max_value = 0.020
    interval_size = 0.001

# Coarse adjustment initialisation ----------------------------------
# Initialise a dictionary to track the rules.
rules = {}
for i in range(0, act):
    for j in range(0, luc):
        key = "from " + luc_names[j] + " to " + luc_names[i + pas]
        rules[key] = [0, 0, 0, 5]
# Set the base random number seed.
base_seed = 1000
# Set the number of simulation runs per iteration.
max_runs = 1
# Initialise a dictionary to store metric values.
coarse_metrics = {}
# Generate a list of values to test.
testing_range = []
testing_pts = max_value / interval_size
for i in range(0, int(testing_pts + 1)):
    temp = i * (max_value - min_value) / testing_pts
    testing_range.append(temp)
    coarse_metrics[temp] = []
# Coarse adjustment -------------------------------------------------
# Output to user.
print("Testing parameter: " + vp)
# Perform the iterative testing.
for x in range(0, len(testing_range)):
    if vp == "theta_it":
        # Set the theta_st value from the testing list.
        theta_it = testing_range[x]
        # Calculate all the relevant meta-parameters.
        # Set the tested parameter values for the self-influence tails.
        d1_high_it_value = high_inertia_point * theta_it
        d1_mid_it_value = mid_inertia_point * theta_it
        d1_low_it_value = low_inertia_point * theta_it
        # Set influence values at distance 2 for self-influence tails.
        d2_high_it_value = high_inertia_point * theta_it * 0.1
        d2_mid_it_value = mid_inertia_point * theta_it * 0.1
        d2_low_it_value = low_inertia_point * theta_it * 0.1
        # Set the conversion parameter values.
        high_conversion_point = high_inertia_point * theta_cp
        mid_conversion_point = mid_inertia_point * theta_cp
        low_conversion_point = low_inertia_point * theta_cp
        # Set influence values at distance 1 for interaction tails.
        d1_high_ct_value = high_inertia_point * theta_ct
        d1_mid_ct_value = mid_inertia_point * theta_ct
        d1_low_ct_value = low_inertia_point * theta_ct
        # Set influence values at distance 2 for interaction tails.
        d2_high_ct_value = high_inertia_point * theta_ct * 0.1
        d2_mid_ct_value = mid_inertia_point * theta_ct * 0.1
        d2_low_ct_value = low_inertia_point * theta_ct * 0.1
    elif vp == "theta_cp":
        theta_cp = testing_range[x]
        # Calculate all the relevant meta-parameters.
        # Set the tested parameter values for the self-influence tails.
        d1_high_it_value = high_inertia_point * theta_it
        d1_mid_it_value = mid_inertia_point * theta_it
        d1_low_it_value = low_inertia_point * theta_it
        # Set influence values at distance 2 for self-influence tails.
        d2_high_it_value = high_inertia_point * theta_it * 0.1
        d2_mid_it_value = mid_inertia_point * theta_it * 0.1
        d2_low_it_value = low_inertia_point * theta_it * 0.1
        # Set the conversion parameter values.
        high_conversion_point = high_inertia_point * theta_cp
        mid_conversion_point = mid_inertia_point * theta_cp
        low_conversion_point = low_inertia_point * theta_cp
        # Set influence values at distance 1 for interaction tails.
        d1_high_ct_value = high_inertia_point * theta_ct
        d1_mid_ct_value = mid_inertia_point * theta_ct
        d1_low_ct_value = low_inertia_point * theta_ct
        # Set influence values at distance 2 for interaction tails.
        d2_high_ct_value = high_inertia_point * theta_ct * 0.1
        d2_mid_ct_value = mid_inertia_point * theta_ct * 0.1
        d2_low_ct_value = low_inertia_point * theta_ct * 0.1
    elif vp == "theta_ct":
        theta_ct = testing_range[x]
        # Calculate all the relevant meta-parameters.
        # Set the tested parameter values for the self-influence tails.
        d1_high_it_value = high_inertia_point * theta_it
        d1_mid_it_value = mid_inertia_point * theta_it
        d1_low_it_value = low_inertia_point * theta_it
        # Set influence values at distance 2 for self-influence tails.
        d2_high_it_value = high_inertia_point * theta_it * 0.1
        d2_mid_it_value = mid_inertia_point * theta_it * 0.1
        d2_low_it_value = low_inertia_point * theta_it * 0.1
        # Set the conversion parameter values.
        high_conversion_point = high_inertia_point * theta_cp
        mid_conversion_point = mid_inertia_point * theta_cp
        low_conversion_point = low_inertia_point * theta_cp
        # Set influence values at distance 1 for interaction tails.
        d1_high_ct_value = high_inertia_point * theta_ct
        d1_mid_ct_value = mid_inertia_point * theta_ct
        d1_low_ct_value = low_inertia_point * theta_ct
        # Set influence values at distance 2 for interaction tails.
        d2_high_ct_value = high_inertia_point * theta_ct * 0.1
        d2_mid_ct_value = mid_inertia_point * theta_ct * 0.1
        d2_low_ct_value = low_inertia_point * theta_ct * 0.1
    # Provide user feedback.
    print("Parameter value: " + str(testing_range[x]))
    coarse_metrics_key = testing_range[x]
    # Set the values for inertia and conversion.
    for i in range(0, act):
        for j in range(0, luc):
            # Specify the neighbourhood rule key.
            key = "from " + luc_names[j] + " to " + luc_names[i + pas]
            # If an inertia rule, set the inertial point and tail.
            if i + pas == j:
                # Evaluate the points and set values.
                if point_levels[j, i] == 3:
                    rules[key][0] = high_inertia_point
                elif point_levels[j, i] == 2:
                    rules[key][0] = mid_inertia_point
                else:
                    rules[key][0] = low_inertia_point
                # Evaluate the tails and set values.
                if tail_levels[j, i] == 3:
                    rules[key][1] = d1_high_it_value
                    rules[key][2] = d2_high_it_value
                elif tail_levels[j, i] == 2:
                    rules[key][1] = d1_mid_it_value
                    rules[key][2] = d2_mid_it_value
                else:
                    rules[key][1] = d1_low_it_value
                    rules[key][2] = d2_low_it_value
            # If a conversion rule, set the conversion point and tail.
            else:
                # Evaluate the points and set values.
                if point_levels[j, i] == 3:
                    rules[key][0] = high_conversion_point
                elif point_levels[j, i] == 2:
                    rules[key][0] = mid_conversion_point
                elif point_levels[j, i] == 1:
                    rules[key][0] = low_conversion_point
                # Evaluate the tails and set values.
                if tail_levels[j, i] == 3:
                    rules[key][1] = d1_high_ct_value
                    rules[key][2] = d2_high_ct_value
                elif tail_levels[j, i] == 2:
                    rules[key][1] = d1_mid_ct_value
                    rules[key][2] = d2_mid_ct_value
                elif tail_levels[j, i] == 1:
                    rules[key][1] = d1_low_ct_value
                    rules[key][2] = d2_low_ct_value
    # Input the rules into the model.
    for i in range(0, luc):
        for j in range(0, act):
            key = "from " + luc_names[i] + " to " + luc_names[j + pas]
            fu_elem = j
            lu_elem = i
            y0 = rules[key][0]
            y1 = rules[key][1]
            y2 = rules[key][2]
            xe = rules[key][3]
            set_lp_rule(project_file, fu_elem, lu_elem, y0, y1, y2, xe)

    # Initialise a dictionary to store the iteration results.
    it_metrics = {}
    # Initialise the tracking of the Fuzzy Kappa values.
    key = "fk"
    it_metrics[key] = [0] * max_runs
    # Initialise the tracking of the Fuzzy Kappa Simulation values.
    key = "fks"
    it_metrics[key] = [0] * max_runs
    # Initialise the tracking of the Area-Weighted Clumpiness error.
    key = "awce"
    it_metrics[key] = [0] * max_runs
    # Generate the simulated output and record the results.
    # Reset the run count to zero.
    run_count = 0
    while run_count < max_runs:
        # Provide user feedback on the run being performed.
        print("Run: " + str(run_count + 1))
        # Generate the seed, input into the model, and run to generate output.
        rseed = base_seed + run_count
        set_rand(project_file, rseed)
        run_metro(project_file, log_file, working_directory,
                  geo_cmd)
        # Read in the simulated map.
        smap = read_map(smap_path)
        # Calculate metrics included for analysis.
        # Calculate Fuzzy Kappa.
        it_metrics["fk"][run_count] = fuzzy_kappa(amap_path, smap_path,
                                                  fuzzy_coefficients)
        # Calculate Fuzzy Kappa Simulation.
        it_metrics["fks"][run_count] = fks(omap_path, amap_path,
                                           smap_path,
                                           fuzzy_trans_coefficients)
        # Calculate AAWCE.
        it_metrics["awce"][run_count] = (
            area_weighted_clu_error(amap, smap, mask, luc, pas, act,
                                    luc_count)
        )
        # Add 1 to iterator to prevent infinite loop!
        run_count = run_count + 1
    # Log the output metrics in the dictionary.
    coarse_metrics[coarse_metrics_key].append(sum(it_metrics["fk"]) /
                                              len(it_metrics["fk"]))
    coarse_metrics[coarse_metrics_key].append(sum(it_metrics["fks"]) /
                                              len(it_metrics["fks"]))
    coarse_metrics[coarse_metrics_key].append(sum(it_metrics["awce"]) /
                                              len(it_metrics["awce"]))
# Logging of results ------------------------------------------------
# Specify the csv file.
if vp == "theta_it":
    metrics_output_file = output_path + "1. theta_it_coarse_cal_output.csv"
elif vp == "theta_cp":
    metrics_output_file = output_path + "2. theta_cp_coarse_cal_output.csv"
elif vp == "theta_ct":
    metrics_output_file = output_path + "3. theta_ct_coarse_cal_output.csv"
# Determine how many metrics to write to output.
store_len = 3
# Generate an empty list to store metric values.
store = [0]*(1 + store_len)
# Write to csv file.
with open (metrics_output_file, "w", newline='') as csv_file:
    writer = csv.writer(csv_file)
    values = [vp]
    for i in range(0, store_len):
        values.append("Metric_" + str(i + 1))
    writer.writerow(values)
    for x in range(0, len(testing_range)):
        store[0] = testing_range[x]
        coarse_metrics_key = testing_range[x]
        for i in range(0, store_len):
            store[i+1] = coarse_metrics[coarse_metrics_key][i]
        writer.writerow(store)
