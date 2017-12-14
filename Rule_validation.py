"""
This module validates the rules obtained at the conclusion of the calibration
method. Validation is performed with assessment of both the calibration &
independent validation data set and models.
"""


# Allows for multi-dimensional array handling.
import numpy as np
# Allows for maps to be stored as 2d arrays.
from read_map import read_map
# Set neighbourhood rules based on a four-point structure.
from set_NR import set_lp_rule
# Set the random number seed in the Metronamica file.
from set_rand import set_rand
# Run the Metronamica model to generate output.
from run_metro import run_metro
# Module for the calculation of Fuzzy Kappa.
from fuzzy_kappa import fuzzy_kappa
# Module for calculation of Fuzzy Kappa Simulation.
from fuzzy_kappa import fks
# Module for calculation of Absolute Area Weighted Avg. Clumpiness Error (AAWCE)
from area_weighted_clu import area_weighted_clu_error
# Interaction with csv file format.
import csv

# Initialisation.
# Specify the base path to the directory containing the empirical neighbourhood
# calibration tool-pack.
base_path = "C:\\Users\\charl\OneDrive\\Documents\\ENC_Py3_release\\"
# Set the case study.
case_study = "Berlin"
# Set the paths to the directories and relevant data.
data_path = base_path + "EU_data\\"
output_path = base_path + "EU_output\\"
# Specify the data map at time slice 1.
map1_path = data_path + case_study + "\\" + case_study.lower() + "_1990.asc"
# Specify the data map at time slice 2.
map2_path = data_path + case_study + "\\" + case_study.lower() + "_2000.asc"
# Specify the data map at time slice 3.
map3_path = data_path + case_study + "\\" + case_study.lower() + "_2006.asc"
# Specify the masking map.
mask_path = data_path + case_study + "\\" + case_study.lower() + "_mask.asc"
# Specify the fuzzy weights for the calculation of Fuzzy Kappa & FKS.
fuzzy_co = data_path + "coeff13.txt"
fuzzy_trans_co = data_path + "coefficients13.txt"
# Set the command line version of Geonamica.
geo_cmd = "C:\\Program Files (x86)\\Geonamica\\Metronamica\\GeonamicaCmd.exe"

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
# Load the rules file, reading the inputs from the specified .csv file.
final_rule_file = (output_path + case_study +
                   "\\Rules\\final_rules.csv")
# Store the rules to be analysed in a dictionary.
rules = {}
# Read inputs from csv file
with open(final_rule_file, 'r', newline='') as f:
    readCSV = csv.reader(f)
    next(f)  # This skips the header line
    for row in readCSV:
        i = row[0]
        j = row[1]
        key = "from " + i + " to " + j
        rules[key] = [0, 0, 0, 5]
        rules[key][0] = float(row[2])
        rules[key][1] = float(row[3])
        rules[key][2] = float(row[4])
        rules[key][3] = float(row[5])
# Specify the fine-tuning calibration method parameters, the base random seed,
# & the maximum number of simulation runs.
base_seed = 1000
max_runs = 5
# Provide user feedback about the case study being tested.
print("Testing case study: " + case_study)
# ----------------------------------------------------------------------------
# Evaluation of calibration data period.
print("Evaluating the calibration data.")
# Specify the map paths.
omap_path = map1_path
amap_path = map2_path
# Specify the maps for the calibration period.
omap = read_map(omap_path)
amap = read_map(amap_path)
mask = read_map(mask_path)
# Specify the log file.
log_file = base_path + "LogSettings.xml"
# Set the working directory, which contains the geoproject file.
working_directory = ("C:\\Geonamica\\Metronamica\\" + case_study)
# Set the project file path.
project_file = working_directory + "\\" + case_study + ".geoproj"
# Specify the simulated output map path.
smap_path = (
    working_directory + "\\Log\\Land_use\\"
    "Land use map_2000-Jan-01 00_00_00.rst"
)
# Set a dictionary to store metric values.
cal_metrics = {}
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
# Set an iterator to track the number of model evaluations performed.
run_count = 0
# Generate the simulated output and track the results.
while run_count < max_runs:
    # Set the random number seed.
    seed = base_seed + run_count
    set_rand(project_file, seed)
    # Run the model.
    run_metro(project_file, log_file, working_directory, geo_cmd)
    # Read in the simulated map.
    smap = read_map(smap_path)
    # Specify the AAWCE clumpiness key.
    key1 = "AAWCE-" + str(run_count)
    # Calculate & store the AAWCE.
    cal_metrics[key1] = area_weighted_clu_error(amap, smap, mask, luc, pas,
                                                act, luc_count)
    # Specify the Fuzzy Kappa key.
    key2 = "FK-" + str(run_count)
    # Calculate & store the Fuzzy Kappa.
    cal_metrics[key2] = fuzzy_kappa(amap_path, smap_path, fuzzy_co)
    # Specify the Fuzzy Kappa Simulation key.
    key3 = "FKS-" + str(run_count)
    # Calculate & store the Fuzzy Kappa Simulation.
    cal_metrics[key3] = fks(omap_path, amap_path, smap_path,
                            fuzzy_trans_co)
    # Add one to iterator (prevents an infinite loop!)
    run_count = run_count + 1
    # Provide user feedback.
    print("Iterations completed: " + str(run_count))
# ----------------------------------------------------------------------------
# Evaluation of validation data period.
print("Evaluating the validation data.")
# Specify the map paths.
omap_path = map2_path
amap_path = map3_path
# Specify the maps for the calibration period.
omap = read_map(omap_path)
amap = read_map(amap_path)
mask = read_map(mask_path)
# Specify the log file.
log_file = base_path + "LogSettings2006.xml"
# Set the working directory, which contains the geoproject file.
working_directory = ("C:\\Geonamica\\Metronamica\\" + case_study + "_2000")
# Set the project file path.
project_file = working_directory + "\\" + case_study + "_2000.geoproj"
# Specify the simulated output map path.
smap_path = (
    working_directory + "\\Log\\Land_use\\"
    "Land use map_2006-Jan-01 00_00_00.rst"
)
# Set a dictionary to store metric values.
val_metrics = {}
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
# Set an iterator to track the number of model evaluations performed.
run_count = 0
# Generate the simulated output and track the results.
while run_count < max_runs:
    # Set the random number seed.
    seed = base_seed + run_count
    set_rand(project_file, seed)
    # Run the model.
    run_metro(project_file, log_file, working_directory, geo_cmd)
    # Read in the simulated map.
    smap = read_map(smap_path)
    # Specify the AAWCE clumpiness key.
    key1 = "AAWCE-" + str(run_count)
    # Calculate & store the AAWCE.
    val_metrics[key1] = area_weighted_clu_error(amap, smap, mask, luc, pas,
                                                act, luc_count)
    # Specify the Fuzzy Kappa key.
    key2 = "FK-" + str(run_count)
    # Calculate & store the Fuzzy Kappa.
    val_metrics[key2] = fuzzy_kappa(amap_path, smap_path, fuzzy_co)
    # Specify the Fuzzy Kappa Simulation key.
    key3 = "FKS-" + str(run_count)
    # Calculate & store the Fuzzy Kappa Simulation.
    val_metrics[key3] = fks(omap_path, amap_path, smap_path,
                            fuzzy_trans_co)
    # Add one to iterator (prevents an infinite loop!)
    run_count = run_count + 1
    # Provide user feedback.
    print("Iterations completed: " + str(run_count))
# ----------------------------------------------------------------------------
# Write metrics to output
metric_output_file = (output_path + case_study + "\\" + case_study +
                      "_final_rule_metrics.csv")
store = [0] * 5
with open (metric_output_file, 'w', newline='') as csv_file:
    writer = csv.writer(csv_file)
    # Save and write a header line to the file
    values = ["Year", "Seed", "FK", "FKS", "AAWCE"]
    writer.writerow(values)
    # Write the calibration metrics.
    for i in range(0, max_runs):
        store[0] = "2000"
        store[1] = i
        key1 = "FK-" + str(i)
        key2 = "FKS-" + str(i)
        key3 = "AAWCE-" + str(i)
        store[2] = cal_metrics[key1]
        store[3] = cal_metrics[key2]
        store[4] = cal_metrics[key3]
        writer.writerow(store)
    # Write the validation metrics.
    for i in range(0, max_runs):
        store[0] = "2006"
        store[1] = i
        key1 = "FK-" + str(i)
        key2 = "FKS-" + str(i)
        key3 = "AAWCE-" + str(i)
        store[2] = val_metrics[key1]
        store[3] = val_metrics[key2]
        store[4] = val_metrics[key3]
        writer.writerow(store)
