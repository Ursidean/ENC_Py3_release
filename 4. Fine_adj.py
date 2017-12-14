"""
This component of the Empirical Neighbourhood Calibration method conducts the
fine calibration, setting neighbourhood rule parameters.
"""

# Allows for multi-dimensional array handling.
import numpy as np
# Allows for maps to be stored as 2d arrays.
from read_map import read_map
# Determines the size (# of cells) of square neighbourhoods.
from considered_distances import considered_distances
# Seed the model with a set of neighbourhood rules.
from seed_rules2 import seed_rules2
# Calculate the enrichment factor.
from enrichment_factor import ef
# Set neighbourhood rules based on a four-point structure.
from set_NR import set_lp_rule
# Set the random number seed in the Metronamica file.
from set_rand import set_rand
# Run the Metronamica model to generate output.
from run_metro import run_metro
# Module for calculation of Fuzzy Kappa.
from fuzzy_kappa import fuzzy_kappa
# Module for calculation of Fuzzy Kappa Simulation.
from fuzzy_kappa import fks
# Module for calculation of Absolute Area Weighted Avg. Clumpiness Error (AAWCE)
from area_weighted_clu import area_weighted_clu_error
# Module for calculating a weighted sum for 3 objectives.
from calc_weighted_sum import ws_3_metrics
# Import the time module to track the duration of the calibration method.
import time
# Import operator module for evaluating lists.
import operator
# Calculate the square-root of a number.
from math import sqrt
# Find the index of values in a numpy array.
from numpy import unravel_index
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
# Specify the original map (data at time slice 1).
omap_path = data_path + case_study + "\\" + case_study.lower() + "_1990.asc"
# Specify the actual map (data at time slice 2).
amap_path = data_path + case_study + "\\" + case_study.lower() + "_2000.asc"
# Specify the masking map.
mask_path = data_path + case_study + "\\" + case_study.lower() + "_mask.asc"
# Specify the fuzzy weights for the calculation of Fuzzy Kappa.
fuzzy_coefficients = data_path + "coeff13.txt"
# Specify the fuzzy transition weights for the calculation of FKS.
fuzzy_trans_coefficients = data_path + "coefficients13.txt"

# Set the land-use class names.
luc_names = ["Natural areas", "Arable land", "Permanent crops", "Pastures",
             "Agricultural areas", "Residential", "Industry & commerce",
             "Recreation areas", "Forest", "Road & rail", "Port area",
             "Airports", "Mine & dump sites", "Fresh water", "Marine water"]
# Set the socio-economic group for each land-use class. This is used to set the
# order of calibration.
socio_eco_group = ["Recreational", "Other", "Other", "Other", "Other",
                   "Residential", "Work", "Recreational", "Recreational",
                   "Other", "Work", "Work", "Work", "Other", "Other"]
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

# Determine the enrichment factor values for the data.
data_ef = ef(luc, max_distance, cdl, cd, N, omap, amap, mask, rows, cols)

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
# Load the points file.
points_file = output_path + case_study + "\\Rules\\points.txt"
points = np.loadtxt(points_file)
# Load the tails file.
tails_file = output_path + case_study + "\\Rules\\tails.txt"
tails = np.loadtxt(tails_file)

# Generate initial rule-set (coarse calibration ending) ----------------------
# Generate the rules in one of two ways. Rules can either be generated
# by seeding the model with the specified meta-parameters, or read in
# the rules from a csv file.
seed = True
if seed is True:
    theta_it = 0.050
    theta_cp = 0.025
    theta_ct = 0.005
    rules = seed_rules2(omap, amap, mask, max_distance, luc_names, luc, act, pas,
                        tails, theta_it, theta_cp, theta_ct, project_file)
else:
    # If false read the rules from an input file.
    initial_rule_file = (output_path + "\\" + case_study +
                         "\\Rules\\initial_rules.csv")
    # Initialise a dictionary for storing rule values.
    rules = {}
    for i in range(0, luc):
        for j in range(0, act):
            key = "from " + luc_names[i] + " to " + luc_names[j + pas]
            rules[key] = [0, 0, 0, 5]
    # Read inputs from csv file
    with open(initial_rule_file, 'r', newline='') as f:
        readCSV = csv.reader(f)
        next(f)  # This skips the header line
        for row in readCSV:
            i = row[0]
            j = row[1]
            key = "from " + i + " to " + j
            rules[key][0] = float(row[2])
            rules[key][1] = float(row[3])
            rules[key][2] = float(row[4])
            rules[key][3] = float(row[5])

# Input the rules into the model, analyse the input points for inclusion.
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

# Fine tuning ----------------------------------------------------------------
# Initialisation of variables.
# Set the Golden Ratio.
gr = (sqrt(5) + 1)/2
# Initialise a counter to track the number of model simulations performed.
iterations_counter = 0
# Specify the maximum number of model simulations.
max_iterations = 1000
# Set the base random number seed & max number of runs.
base_seed = 1000
max_runs = 1
# Specify the maximum and minimum bounding values for different neighbourhood
# components, and the golden section search tolerance:
# Inertia point
max_ip = 1000
min_ip = 250
gss_ip_tol = (max_ip - min_ip)/100
# The self-influence tail influence value at distance 1.
max_st = 100
min_st = 0
gss_st_tol = (max_st - min_st)/100
# The conversion point
max_cp = 100
min_cp = 0
gss_cp_tol = (max_cp - min_cp)/100
# The interaction tail at distance 1.
max_it = 100
min_it = 0
gss_it_tol = (max_it - min_it)/100
# Initialise two lists to track the rule adjusted by name.
rule_from_tracker = []
rule_to_tracker = []
# Initialise a list to track the point that is fine tuned,
# and the final value taken.
pt_tracker = []
pt_value = []
# Initialisation of metrics (One Kappa, One Clumpiness.
# Set the weighting values for weighted sum calculations.
w1 = 1/4
w2 = 1/4
w3 = 1/2
# Set the transformation metric ranges.
r1 = [0.900, 0.950]
r2 = [0.000, 0.100]
r3 = [0.000, 0.050]
# Set whether to maximise or minimise the objectives.
s1 = "maximise"
s2 = "maximise"
s3 = "minimise"
# Initialise a list to track the individual metrics averaged over the number of
# runs that are performed.
clu_log = []
kappa_log = []
kappa_sim_log =[]
# Initialise a set of temp variables for storing metric values.
clu_temp = [0] * max_runs
kappa_temp = [0] * max_runs
kappa_sim_temp = [0] * max_runs
# Initialise a list to track the weighted sum.
ws_log = []

# Determine the starting point metrics.
# Reset the run count.
run_count = 0
while run_count < max_runs:
    # Set the random seed.
    rseed = base_seed + run_count
    set_rand(project_file, rseed)
    # Run Metronamica to generate output.
    run_metro(project_file, log_file, working_directory, geo_cmd)
    # Read in the simulated map.
    smap = read_map(smap_path)
    # Calculate the AAWCE.
    clu_temp[run_count] = area_weighted_clu_error(amap, smap, mask, luc, pas,
                                                  act, luc_count)
    # Calculate the run Fuzzy Kappa.
    kappa_temp[run_count] = fuzzy_kappa(amap_path, smap_path,
                                        fuzzy_coefficients)
    # Calculate the run Fuzzy Kappa Sim.
    kappa_sim_temp[run_count] = fks(omap_path, amap_path, smap_path,
                                    fuzzy_trans_coefficients)
    # Add one to iterator (prevents an infinite loop!)
    run_count = run_count + 1
# Find the average over the number of runs performed.
clu_avg = sum(clu_temp) / len(clu_temp)
kappa_avg = sum(kappa_temp) / len(kappa_temp)
kappa_sim_avg = sum(kappa_sim_temp) / len(kappa_sim_temp)
# Now determine the weighted sum (ws).
ws_avg = ws_3_metrics(kappa_avg, kappa_sim_avg, clu_avg, w1, w2, w3,
                      r1, r2, r3, s1, s2, s3)

# Track the start of the calibration.
start = time.time()

# Fine tuning of the inertia points.
print("Fine tuning the inertial points")
# Determine the order of calibration of points.
ip_order = np.array([0]*act)
# Evaluate the socio-economic groups to determine the order for
# inertia points and self-influence rules.
# Allocate a numerical value.
for i in range(0, act):
    if socio_eco_group[i + pas] == "Residential":
        ip_order[i] = 4
    elif socio_eco_group[i + pas] == "Work":
        ip_order[i] = 3
    elif socio_eco_group[i + pas] == "Recreation":
        ip_order[i] = 2
    else:
        ip_order[i] = 1
# Rank and determine order.
temp = ip_order.argsort()
ranks = np.empty(len(ip_order), int)
ranks[temp] = np.arange(len(ip_order))
ip_order = ranks + 1
# Apply the Golden Section Search to all inertia points.
while sum(ip_order) > 0:
    if iterations_counter > max_iterations:
        # Set values to 0 to prevent further calibration.
        # Identify max value, corresponding to parameter that will be
        # tested.
        par_index, par_value = max(enumerate(ip_order), key=operator.itemgetter(1))
        # Set the corresponding value in the list to 0.
        ip_order[par_index] = 0
    else:
        # Identify max value, corresponding to parameter that will be
        # tested.
        par_index, par_value = max(enumerate(ip_order), key = operator.itemgetter(1))
        # Set the corresponding value in the list to 0.
        ip_order[par_index] = 0
        # Specify the function and land-use element index values for
        # appending the geoproject file.
        fu_elem = par_index
        lu_elem = par_index + 1
        rule_name = ("from " + luc_names[lu_elem] +
                     " to " + luc_names[fu_elem + pas])
        print("Calibrating rule: " + rule_name)
        # Track the rule being adjusted
        rule_from_tracker.append(lu_elem)
        rule_to_tracker.append(fu_elem)
        pt_tracker.append(0)
        # Find the influence values for the specified rule
        adjust_rule_values = rules[rule_name]
        # Fix the requisite parameter values.
        old_y0 = adjust_rule_values[0]
        y1 = adjust_rule_values[1]
        y2 = adjust_rule_values[2]
        xe = adjust_rule_values[3]
        # Apply the Golden Section Search to optimise the output.
        a = min_ip
        b = max_ip
        while (abs(a - b) > gss_ip_tol):
            # Calculate the interval values.
            c = b - (b - a) / gr
            d = a + (b - a) / gr
            # Set the neighbourhood rule to the lower interval values.
            set_lp_rule(project_file, fu_elem, lu_elem, c, y1, y2, xe)
            # Run the model, evaluate performance.
            kappa_temp_c = [0] * max_runs
            kappa_sim_temp_c = [0] * max_runs
            clu_temp_c = [0] * max_runs
            # Reset the run count.
            run_count = 0
            while run_count < max_runs:
                # Set the random seed
                rseed = base_seed + run_count
                set_rand(project_file, rseed)
                # Run Metronamica to generate output
                run_metro(project_file, log_file, working_directory,
                          geo_cmd)
                iterations_counter = iterations_counter + 1
                # Read in output
                smap = read_map(smap_path)
                # Determine objectives.
                # Calculate the run AAWCE.
                clu_temp_c[run_count] = area_weighted_clu_error(amap, smap, mask,
                                                                luc, pas, act,
                                                                luc_count)
                # Calculate the run Fuzzy Kappa.
                kappa_temp_c[run_count] = fuzzy_kappa(amap_path, smap_path,
                                                      fuzzy_coefficients)
                # Calculate the run Fuzzy Kappa Simulation.
                kappa_sim_temp_c[run_count] = fks(omap_path, amap_path,
                                                  smap_path,
                                                  fuzzy_trans_coefficients)
                # Add one to iterator to prevent infinite loop!
                run_count = run_count + 1
            # Find the average over the number of runs performed.
            clu_avg_c = sum(clu_temp_c) / len(clu_temp_c)
            kappa_avg_c = sum(kappa_temp_c) / len(kappa_temp_c)
            kappa_sim_avg_c = sum(kappa_sim_temp_c) / len(kappa_sim_temp_c)
            # Now determine the weighted sum.
            ws_c = ws_3_metrics(kappa_avg_c, kappa_sim_avg_c, clu_avg_c,
                                w1, w2, w3, r1, r2, r3, s1, s2, s3)
            # Now set the neighbourhood rule to the higher interval value.
            set_lp_rule(project_file, fu_elem, lu_elem, d, y1, y2, xe)
            # Run the model, evaluate performance.
            kappa_temp_d = [0] * max_runs
            kappa_sim_temp_d = [0] * max_runs
            clu_temp_d = [0] * max_runs
            # Reset the run count.
            run_count = 0
            while run_count < max_runs:
                # Set the random seed
                rseed = base_seed + run_count
                set_rand(project_file, rseed)
                # Run Metronamica to generate output
                run_metro(project_file, log_file, working_directory,
                          geo_cmd)
                iterations_counter = iterations_counter+ 1
                # Read in output
                smap = read_map(smap_path)
                # Determine objectives.
                # Calculate the run AAWCE.
                clu_temp_d[run_count] = area_weighted_clu_error(amap, smap, mask,
                                                                luc, pas, act,
                                                                luc_count)
                # Calculate the run Fuzzy Kappa.
                kappa_temp_d[run_count] = fuzzy_kappa(amap_path, smap_path,
                                                      fuzzy_coefficients)
                # Calculate the run Fuzzy Kappa Simulation.
                kappa_sim_temp_d[run_count] = fks(omap_path, amap_path,
                                                  smap_path,
                                                  fuzzy_trans_coefficients)
                # Add one to iterator to prevent infinite loop!
                run_count = run_count + 1
            # Find the average over the number of runs performed.
            clu_avg_d = sum(clu_temp_d) / len(clu_temp_d)
            kappa_avg_d = sum(kappa_temp_d) / len(kappa_temp_d)
            kappa_sim_avg_d = sum(kappa_sim_temp_d) / len(kappa_sim_temp_d)
            # Now determine the weighted sum.
            ws_d = ws_3_metrics(kappa_avg_d, kappa_sim_avg_d, clu_avg_d,
                                w1, w2, w3, r1, r2, r3, s1, s2, s3)
            # Add step here to break out if objectives are not changing.
            if ws_avg == ws_c and ws_avg == ws_d:
                break
            # Otherwise evaluate the result to close the bracket.
            # As we want to maximise, evaluation is as below:
            else:
                if ws_c > ws_d:
                    b = d
                else:
                    a = c
        # This section of code is run until the tolerance criteria are violated.
        # When this happens, the final evaluation is conducted.
        # Set the final value, input into the model.
        if ws_avg >= ws_c and ws_avg >= ws_d:
            # Reset the value if there was no improvement.
            rules[rule_name][0] = old_y0
            set_lp_rule(project_file, fu_elem, lu_elem, old_y0, y1, y2, xe)
            kappa_log.append(kappa_avg)
            kappa_sim_log.append(kappa_sim_avg)
            clu_log.append(clu_avg)
            ws_log.append(ws_avg)
            pt_value.append(old_y0)
        # If the weighted sum for c is better, take that value. Insert into
        # model and log output.
        elif ws_c > ws_d:
            # First track the rule.
            final_value = c
            set_lp_rule(project_file, fu_elem, lu_elem, final_value, y1, y2, xe)
            pt_value.append(final_value)
            rules[rule_name][0] = final_value
            # Now update the best set of objectives.
            kappa_avg = kappa_avg_c
            kappa_sim_avg = kappa_sim_avg_c
            clu_avg = clu_avg_c
            ws_avg = ws_c
            # Log the output.
            kappa_log.append(kappa_avg)
            kappa_sim_log.append(kappa_sim_avg)
            clu_log.append(clu_avg)
            ws_log.append(ws_avg)
        # If the weighted sum for d is better, take that value. Insert into
        # model and log output.
        else:
            # First track the rule.
            final_value = d
            set_lp_rule(project_file, fu_elem, lu_elem, final_value, y1, y2, xe)
            pt_value.append(final_value)
            rules[rule_name][0] = final_value
            # Now update the best set of objectives.
            kappa_avg = kappa_avg_d
            kappa_sim_avg = kappa_sim_avg_d
            clu_avg = clu_avg_d
            ws_avg = ws_d
            # Log the output.
            kappa_log.append(kappa_avg)
            kappa_sim_log.append(kappa_sim_avg)
            clu_log.append(clu_avg)
            ws_log.append(ws_avg)
# This section runs until all inertia points have been appended.
# -----------------------------------------------------------------------------
# Fine tuning of self-influence tails.
print("Fine tuning the inertial tails")
# Initialise a list for self-influence tail rule order tuning.
st_order = np.array([0]*act)
# Evaluate the socio-economic groups to determine the order for
# inertia points and self-influence rules.
# Allocate a numerical value.
for i in range(0, act):
    if socio_eco_group[i + pas] == "Residential":
        st_order[i] = 4
    elif socio_eco_group[i + pas] == "Work":
        st_order[i] = 3
    elif socio_eco_group[i + pas] == "Recreation":
        st_order[i] = 2
    else:
        st_order[i] = 1
# Rank and determine order.
temp = st_order.argsort()
ranks = np.empty(len(st_order), int)
ranks[temp] = np.arange(len(st_order))
st_order = ranks + 1
# Apply the Golden Section Search to self-influence tails.
while sum(st_order) > 0:
    if iterations_counter > max_iterations:
        # Set values to 0 to prevent further calibration.
        # Identify max value, corresponding to parameter that will be
        # tested.
        par_index, par_value = max(enumerate(st_order), key=operator.itemgetter(1))
        # Set the corresponding value in the list to 0.
        st_order[par_index] = 0
    else:
        # Identify max value, corresponding to parameter that will be
        # tested.
        par_index, par_value = max(enumerate(st_order), key=operator.itemgetter(1))
        # Set the corresponding value in the list to 0.
        st_order[par_index] = 0
        # Specify the function and land-use element index values for
        # appending the geoproject file.
        fu_elem = par_index
        lu_elem = par_index + 1
        rule_name = ("from " + luc_names[lu_elem] +
                     " to " + luc_names[fu_elem + pas])
        print("Calibrating rule: " + rule_name)
        # Track the rule being adjusted
        rule_from_tracker.append(lu_elem)
        rule_to_tracker.append(fu_elem)
        pt_tracker.append(1)
        # Find the influence values for the specified rule
        adjust_rule_values = rules[rule_name]
        # Fix the requisite parameter values.
        y0 = adjust_rule_values[0]
        old_y1 = adjust_rule_values[1]
        old_y2 = adjust_rule_values[2]
        xe = adjust_rule_values[3]
        # Apply the Golden Section Search to optimise the output.
        a = min_st
        b = max_st
        while (abs(a - b) > gss_st_tol):
            # Calculate the interval values.
            c = b - (b - a) / gr
            d = a + (b - a) / gr
            # Set the neighbourhood rule to the lower interval values.
            c2 = 0.1 * c
            set_lp_rule(project_file, fu_elem, lu_elem, y0, c, c2, xe)
            # Run the model, evaluate performance.
            kappa_temp_c = [0] * max_runs
            kappa_sim_temp_c = [0] * max_runs
            clu_temp_c = [0] * max_runs
            # Reset the run count.
            run_count = 0
            while run_count < max_runs:
                # Set the random seed
                rseed = base_seed + run_count
                set_rand(project_file, rseed)
                # Run Metronamica to generate output
                run_metro(project_file, log_file, working_directory,
                          geo_cmd)
                iterations_counter = iterations_counter + 1
                # Read in output
                smap = read_map(smap_path)
                # Determine objectives.
                # Calculate the run AAWCE.
                clu_temp_c[run_count] = area_weighted_clu_error(amap, smap, mask,
                                                                luc, pas, act,
                                                                luc_count)
                # Calculate the run Fuzzy Kappa.
                kappa_temp_c[run_count] = fuzzy_kappa(amap_path, smap_path,
                                                      fuzzy_coefficients)
                # Calculate the run Fuzzy Kappa Simulation.
                kappa_sim_temp_c[run_count] = fks(omap_path, amap_path,
                                                  smap_path,
                                                  fuzzy_trans_coefficients)
                # Add one to iterator to prevent infinite loop!
                run_count = run_count + 1
            # Find the average over the number of runs performed.
            clu_avg_c = sum(clu_temp_c) / len(clu_temp_c)
            kappa_avg_c = sum(kappa_temp_c) / len(kappa_temp_c)
            kappa_sim_avg_c = sum(kappa_sim_temp_c) / len(kappa_sim_temp_c)
            # Now determine the weighted sum.
            ws_c = ws_3_metrics(kappa_avg_c, kappa_sim_avg_c, clu_avg_c,
                                w1, w2, w3, r1, r2, r3, s1, s2, s3)
            # Now set the neighbourhood rule to the higher interval value.
            d2 = 0.1 * d
            set_lp_rule(project_file, fu_elem, lu_elem, y0, d, d2, xe)
            # Run the model, evaluate performance.
            kappa_temp_d = [0] * max_runs
            kappa_sim_temp_d =[0] * max_runs
            clu_temp_d = [0] * max_runs
            # Reset the run count.
            run_count = 0
            while run_count < max_runs:
                # Set the random seed
                rseed = base_seed + run_count
                set_rand(project_file, rseed)
                # Run Metronamica to generate output
                run_metro(project_file, log_file, working_directory,
                          geo_cmd)
                iterations_counter = iterations_counter + 1
                # Read in output
                smap = read_map(smap_path)
                # Determine objectives.
                # Calculate the run AAWCE.
                clu_temp_d[run_count] = area_weighted_clu_error(amap, smap, mask,
                                                                luc, pas, act,
                                                                luc_count)
                # Calculate the run Fuzzy Kappa.
                kappa_temp_d[run_count] = fuzzy_kappa(amap_path, smap_path,
                                                      fuzzy_coefficients)
                # Calculate the run Fuzzy Kappa.
                kappa_sim_temp_d[run_count] = fks(omap_path, amap_path,
                                                  smap_path,
                                                  fuzzy_trans_coefficients)
                # Add one to iterator to prevent infinite loop!
                run_count = run_count + 1
            # Find the average over the number of runs performed.
            clu_avg_d = sum(clu_temp_d) / len(clu_temp_d)
            kappa_avg_d = sum(kappa_temp_d) / len(kappa_temp_d)
            kappa_sim_avg_d = sum(kappa_sim_temp_d) / len(kappa_sim_temp_d)
            # Now determine the weighted sum.
            ws_d = ws_3_metrics(kappa_avg_d, kappa_sim_avg_d, clu_avg_d,
                                w1, w2, w3, r1, r2, r3, s1, s2, s3)
            # Add step here to break out if objectives are not changing.
            if ws_avg == ws_c and ws_avg == ws_d:
                break
            # Otherwise evaluate the result to close the bracket.
            # As we want to maximise, evaluation is as below:
            else:
                if ws_c > ws_d:
                    b = d
                else:
                    a = c
        # This section of code is run until the tolerance criteria are violated.
        # When this happens, the final evaluation is conducted.
        # Set the final value, input into the model.
        if ws_avg >= ws_c and ws_avg >= ws_d:
            rules[rule_name][1] = old_y1
            rules[rule_name][2] = old_y2
            set_lp_rule(project_file, fu_elem, lu_elem, y0, old_y1,
                        old_y2, xe)
            pt_value.append(old_y1)
            kappa_log.append(kappa_avg)
            kappa_sim_log.append(kappa_sim_avg)
            clu_log.append(clu_avg)
            ws_log.append(ws_avg)
        # If the weighted sum for c is better, take that value. Insert into
        # model and log output.
        elif ws_c > ws_d:
            # First track the rule.
            final_value = c
            final_value_2 = final_value * 0.1
            set_lp_rule(project_file, fu_elem, lu_elem, y0, final_value,
                        final_value_2, xe)
            pt_value.append(final_value)
            rules[rule_name][1] = final_value
            rules[rule_name][2] = final_value_2
            # Now update the best set of objectives.
            kappa_avg = kappa_avg_c
            kappa_sim_avg = kappa_sim_avg_c
            clu_avg = clu_avg_c
            ws_avg = ws_c
            # Log the output.
            kappa_log.append(kappa_avg)
            kappa_sim_log.append(kappa_sim_avg)
            clu_log.append(clu_avg)
            ws_log.append(ws_avg)
        # If the weighted sum for d is better, take that value. Insert into
        # model and log output.
        else:
            # First track the rule.
            final_value = d
            final_value_2 = final_value * 0.1
            set_lp_rule(project_file, fu_elem, lu_elem, y0, final_value,
                        final_value_2, xe)
            pt_value.append(final_value)
            rules[rule_name][1] = final_value
            rules[rule_name][2] = final_value_2
            # Now update the best set of objectives.
            kappa_avg = kappa_avg_d
            kappa_sim_avg = kappa_sim_avg_d
            clu_avg = clu_avg_d
            ws_avg = ws_d
            # Log the output.
            kappa_log.append(kappa_avg)
            kappa_sim_log.append(kappa_sim_avg)
            clu_log.append(clu_avg)
            ws_log.append(ws_avg)
# This section runs until all self-influence tails have been appended.
# -----------------------------------------------------------------------------
# Fine tuning of conversion points.
print("Fine tuning the conversion points")
temp_cp_order = np.zeros(shape=(luc,act))
# Evaluate the socio-economic groups to determine the order for
# conversion points. First, allocate a numerical value depending on
# class converted to.
for i in range(0, luc):
    for j in range(0, act):
        if i == j + pas:
            # Skip if an inertia point.
            pass
        elif points[i, j] == 0:
            # Skip if no conversion point in model
            pass
        else:
            # Evaluate the conversions, assign rank based on class that is
            # being converted to.
            if socio_eco_group[j + pas] == "Residential":
                temp_cp_order[i, j] = 4
            elif socio_eco_group[j + pas] == "Work":
                temp_cp_order[i, j] = 3
            elif socio_eco_group[j + pas] == "Recreation":
                temp_cp_order[i, j] = 2
            else:
                temp_cp_order[i, j] = 1
# Rank and determine order. Start by initialising an array to track order.
cp_order = np.zeros(shape=(luc, act))
# Iterate through and assign rank. Start at 1.
ranking = 1
iterator = 1
while iterator < 5:
    # Iterate through the ranks.
    for i in range(0, luc):
        for j in range(0, act):
            if temp_cp_order[i, j] == iterator:
                # If the equivalent value is found, assign current
                # ranking value.
                cp_order[i, j] = ranking
                # Add one to the ranking value.
                ranking = ranking + 1
    # Add one to the iterator.
    iterator = iterator + 1
# Apply the Golden-Section Search to optimise the output.
while np.sum(cp_order > 0):
    if iterations_counter > max_iterations:
        # Set values to 0 to prevent further calibration.
        # Identify max value, corresponding to parameter that will be
        # tested.
        par_index = unravel_index(cp_order.argmax(), cp_order.shape)
        # Set the corresponding value in the list to 0.
        cp_order[par_index[0], par_index[1]] = 0
    else:
        # Identify max value, corresponding to parameter that will be
        # tested.
        par_index = unravel_index(cp_order.argmax(), cp_order.shape)
        # Set the corresponding value in the list to 0.
        cp_order[par_index[0], par_index[1]] = 0
        # Specify the function and land-use element index values for
        # appending the geoproject file.
        lu_elem = par_index[0]
        fu_elem = par_index[1]
        rule_name = ("from " + luc_names[lu_elem] +
                     " to " + luc_names[fu_elem + pas])
        print("Calibrating rule: " + rule_name)
        # Track the rule being adjusted
        rule_from_tracker.append(lu_elem)
        rule_to_tracker.append(fu_elem)
        pt_tracker.append(0)
        # Find the influence values for the specified rule
        adjust_rule_values = rules[rule_name]
        # Fix the requisite parameter values.
        old_y0 = adjust_rule_values[0]
        y1 = adjust_rule_values[1]
        y2 = adjust_rule_values[2]
        xe = adjust_rule_values[3]
        # Apply the Golden Section Search to optimise the output.
        a = min_cp
        b = max_cp
        while (abs(a - b) > gss_cp_tol):
            # Calculate the interval values.
            c = b - (b - a) / gr
            d = a + (b - a) / gr
            # Set the neighbourhood rule to the lower interval values.
            set_lp_rule(project_file, fu_elem, lu_elem, c, y1, y2, xe)
            # Run the model, evaluate performance.
            kappa_temp_c = [0] * max_runs
            kappa_sim_temp_c = [0] * max_runs
            clu_temp_c = [0] * max_runs
            # Reset the run count.
            run_count = 0
            while run_count < max_runs:
                # Set the random seed
                rseed = base_seed + run_count
                set_rand(project_file, rseed)
                # Run Metronamica to generate output
                run_metro(project_file, log_file, working_directory, geo_cmd)
                iterations_counter = iterations_counter + 1
                # Read in output
                smap = read_map(smap_path)
                # Determine objectives.
                # Calculate the run AAWCE.
                clu_temp_c[run_count] = area_weighted_clu_error(amap, smap, mask,
                                                                luc, pas, act,
                                                                luc_count)
                # Calculate the run Fuzzy Kappa.
                kappa_temp_c[run_count] = fuzzy_kappa(amap_path, smap_path,
                                                      fuzzy_coefficients)
                # Calculate the run Fuzzy Kappa Simulation.
                kappa_sim_temp_c[run_count] = fks(omap_path, amap_path,
                                                  smap_path,
                                                  fuzzy_trans_coefficients)
                # Add one to iterator to prevent infinite loop!
                run_count = run_count + 1
            # Find the average over the number of runs performed.
            clu_avg_c = sum(clu_temp_c) / len(clu_temp_c)
            kappa_avg_c = sum(kappa_temp_c) / len(kappa_temp_c)
            kappa_sim_avg_c = sum(kappa_sim_temp_c) / len(kappa_sim_temp_c)
            # Now determine the weighted sum.
            ws_c = ws_3_metrics(kappa_avg_c, kappa_sim_avg_c, clu_avg_c,
                                w1, w2, w3, r1, r2, r3, s1, s2, s3)
            # Now set the neighbourhood rule to the higher interval value.
            set_lp_rule(project_file, fu_elem, lu_elem, d, y1, y2, xe)
            # Run the model, evaluate performance.
            kappa_temp_d = [0] * max_runs
            kappa_sim_temp_d = [0] * max_runs
            clu_temp_d = [0] * max_runs
            # Reset the run count.
            run_count = 0
            while run_count < max_runs:
                # Set the random seed
                rseed = base_seed + run_count
                set_rand(project_file, rseed)
                # Run Metronamica to generate output
                run_metro(project_file, log_file, working_directory, geo_cmd)
                iterations_counter = iterations_counter+ 1
                # Read in output
                smap = read_map(smap_path)
                # Determine objectives.
                # Calculate the run AAWCE.
                clu_temp_d[run_count] = area_weighted_clu_error(amap, smap, mask,
                                                                luc, pas, act,
                                                                luc_count)
                # Calculate the run Fuzzy Kappa.
                kappa_temp_d[run_count] = fuzzy_kappa(amap_path, smap_path,
                                                      fuzzy_coefficients)
                # Calculate the run Fuzzy Kappa Simulation.
                kappa_sim_temp_d[run_count] = fks(omap_path, amap_path,
                                                  smap_path,
                                                  fuzzy_trans_coefficients)
                # Add one to iterator to prevent infinite loop!
                run_count = run_count + 1
            # Find the average over the number of runs performed.
            clu_avg_d = sum(clu_temp_d) / len(clu_temp_d)
            kappa_avg_d = sum(kappa_temp_d) / len(kappa_temp_d)
            kappa_sim_avg_d = sum(kappa_sim_temp_d) / len(kappa_sim_temp_d)
            # Now determine the weighted sum.
            ws_d = ws_3_metrics(kappa_avg_d, kappa_sim_avg_d, clu_avg_d,
                                w1, w2, w3, r1, r2, r3, s1, s2, s3)
            # Add step here to break out if objectives are not changing.
            if ws_avg == ws_c and ws_avg == ws_d:
                break
            # Otherwise evaluate the result to close the bracket.
            # As we want to maximise, evaluation is as below:
            else:
                if ws_c > ws_d:
                    b = d
                else:
                    a = c
        # This section of code is run until the tolerance criteria are violated.
        # When this happens, the final evaluation is conducted.
        # Set the final value, input into the model.
        if ws_avg >= ws_c and ws_avg >= ws_d:
            # Reset the value if there was no improvement.
            rules[rule_name][0] = old_y0
            set_lp_rule(project_file, fu_elem, lu_elem, old_y0, y1, y2, xe)
            kappa_log.append(kappa_avg)
            kappa_sim_log.append(kappa_sim_avg)
            clu_log.append(clu_avg)
            ws_log.append(ws_avg)
            pt_value.append(old_y0)
        # If the weighted sum for c is better, take that value. Insert into
        # model and log output.
        elif ws_c > ws_d:
            # First track the rule.
            final_value = c
            set_lp_rule(project_file, fu_elem, lu_elem, final_value, y1, y2, xe)
            pt_value.append(final_value)
            rules[rule_name][0] = final_value
            # Now update the best set of objectives.
            kappa_avg = kappa_avg_c
            kappa_sim_avg = kappa_sim_avg_c
            clu_avg = clu_avg_c
            ws_avg = ws_c
            # Log the output.
            kappa_log.append(kappa_avg)
            kappa_sim_log.append(kappa_sim_avg)
            clu_log.append(clu_avg)
            ws_log.append(ws_avg)
        # If the weighted sum for d is better, take that value. Insert into
        # model and log output.
        else:
            # First track the rule.
            final_value = d
            set_lp_rule(project_file, fu_elem, lu_elem, final_value, y1, y2, xe)
            pt_value.append(final_value)
            rules[rule_name][0] = final_value
            # Now update the best set of objectives.
            kappa_avg = kappa_avg_d
            kappa_sim_avg = kappa_sim_avg_d
            clu_avg = clu_avg_d
            ws_avg = ws_d
            # Log the output.
            kappa_log.append(kappa_avg)
            kappa_sim_log.append(kappa_sim_avg)
            clu_log.append(clu_avg)
            ws_log.append(ws_avg)
# This section runs until all conversion points have been appended.
# -----------------------------------------------------------------------------
# Fine tuning of interaction tails.
print("Fine tuning the interaction tails")
# Evaluate the socio-economic groups to determine the order for
# interactive rules. First, allocate a numerical value depending on
# class converted to.
# Rank and determine order. Start by initialising an array to track order.
temp_it_order = np.zeros(shape=(luc, act))
for i in range(0, luc):
    for j in range(0, act):
        if i == j + pas:
            # Skip if an inertia point.
            pass
        elif tails[i, j] == 0:
            # Skip if no interactive rule.
            pass
        else:
            # Evaluate the conversions, assign rank based on class that is
            # being converted to.
            if socio_eco_group[j + pas] == "Residential":
                temp_it_order[i, j] = 4
            elif socio_eco_group[j + pas] == "Work":
                temp_it_order[i, j] = 3
            elif socio_eco_group[j + pas] == "Recreation":
                temp_it_order[i, j] = 2
            else:
                temp_it_order[i, j] = 1
# Rank and determine order. Start by initialising an array to track order.
it_order = np.zeros(shape=(luc, act))
# Iterate through and assign rank. Start at 1.
ranking = 1
iterator = 1
while iterator < 5:
    # Iterate through the ranks.
    for i in range(0, luc):
        for j in range(0, act):
            if temp_it_order[i, j] == iterator:
                # If the equivalent value is found, assign current
                # ranking value.
                it_order[i, j] = ranking
                # Add one to the ranking value.
                ranking = ranking + 1
    # Add one to the iterator.
    iterator = iterator + 1
# Apply the Golden Section Search to interactive tails.
while np.sum(it_order) > 0:
    if iterations_counter > max_iterations:
        # Set values to 0 to prevent further testing.
        # Identify max value, corresponding to parameter that will be
        # tested.
        par_index = unravel_index(it_order.argmax(), it_order.shape)
        # Set the corresponding value in the list to 0.
        it_order[par_index[0], par_index[1]] = 0
    else:
        # Identify max value, corresponding to parameter that will be
        # tested.
        par_index = unravel_index(it_order.argmax(), it_order.shape)
        # Set the corresponding value in the list to 0.
        it_order[par_index[0], par_index[1]] = 0
        # Specify the function and land-use element index values for
        # appending the geoproject file.
        lu_elem = par_index[0]
        fu_elem = par_index[1]
        # Specify the rule for tuning.
        rule_name = ("from " + luc_names[lu_elem] +
                     " to " + luc_names[fu_elem + pas])
        print("Calibrating rule: " + rule_name)
        # Track the rule being adjusted
        rule_from_tracker.append(lu_elem)
        rule_to_tracker.append(fu_elem)
        pt_tracker.append(1)
        # Find the influence values for the specified rule
        adjust_rule_values = rules[rule_name]
        # Fix the requisite parameter values.
        y0 = adjust_rule_values[0]
        old_y1 = adjust_rule_values[1]
        old_y2 = adjust_rule_values[2]
        xe = adjust_rule_values[3]
        # Apply the Golden Section Search to optimise the output.
        a = min_it
        b = max_it
        while (abs(a - b) > gss_it_tol):
            # Calculate the interval values.
            c = b - (b - a) / gr
            d = a + (b - a) / gr
            # Set the neighbourhood rule to the lower interval values.
            c2 = 0.1*c
            set_lp_rule(project_file, fu_elem, lu_elem, y0, c, c2, xe)
            # Run the model, evaluate performance.
            kappa_temp_c = [0] * max_runs
            kappa_sim_temp_c = [0] * max_runs
            clu_temp_c = [0] * max_runs
            # Reset the run count.
            run_count = 0
            while run_count < max_runs:
                # Set the random seed
                rseed = base_seed + run_count
                set_rand(project_file, rseed)
                # Run Metronamica to generate output
                run_metro(project_file, log_file, working_directory,
                          geo_cmd)
                iterations_counter = iterations_counter + 1
                # Read in output
                smap = read_map(smap_path)
                # Determine objectives.
                # Calculate the run AAWCE.
                clu_temp_c[run_count] = area_weighted_clu_error(amap, smap, mask,
                                                                luc, pas, act,
                                                                luc_count)
                # Calculate the run Fuzzy Kappa.
                kappa_temp_c[run_count] = fuzzy_kappa(amap_path, smap_path,
                                                      fuzzy_coefficients)
                # Calculate the run Fuzzy Kappa Simulation.
                kappa_sim_temp_c[run_count] = fks(omap_path, amap_path,
                                                  smap_path,
                                                  fuzzy_trans_coefficients)
                # Add one to iterator to prevent infinite loop!
                run_count = run_count + 1
            # Find the average over the number of runs performed.
            clu_avg_c = sum(clu_temp_c) / len(clu_temp_c)
            kappa_avg_c = sum(kappa_temp_c) / len(kappa_temp_c)
            kappa_sim_avg_c = sum(kappa_sim_temp_c) / len(kappa_sim_temp_c)
            # Now determine the weighted sum.
            ws_c = ws_3_metrics(kappa_avg_c, kappa_sim_avg_c, clu_avg_c,
                                w1, w2, w3, r1, r2, r3, s1, s2, s3)
            # Now set the neighbourhood rule to the higher interval value.
            d2 = 0.1 * d
            set_lp_rule(project_file, fu_elem, lu_elem, y0, d, d2, xe)
            # Run the model, evaluate performance.
            kappa_temp_d = [0] * max_runs
            kappa_sim_temp_d = [0] * max_runs
            clu_temp_d = [0] * max_runs
            # Reset the run count.
            run_count = 0
            while run_count < max_runs:
                # Set the random seed
                rseed = base_seed + run_count
                set_rand(project_file, rseed)
                # Run Metronamica to generate output
                run_metro(project_file, log_file, working_directory,
                          geo_cmd)
                iterations_counter = iterations_counter + 1
                # Read in output
                smap = read_map(smap_path)
                # Determine objectives.
                # Calculate the run AAWCE.
                clu_temp_d[run_count] = area_weighted_clu_error(amap, smap, mask,
                                                                luc, pas, act,
                                                                luc_count)
                # Calculate the run Fuzzy Kappa.
                kappa_temp_d[run_count] = fuzzy_kappa(amap_path, smap_path,
                                                      fuzzy_coefficients)
                # Calculate the run Fuzzy Kappa Simulation.
                kappa_sim_temp_d[run_count] = fks(omap_path, amap_path,
                                                  smap_path,
                                                  fuzzy_trans_coefficients)
                # Add one to iterator to prevent infinite loop!
                run_count = run_count + 1
            # Find the average over the number of runs performed.
            clu_avg_d = sum(clu_temp_d) / len(clu_temp_d)
            kappa_avg_d = sum(kappa_temp_d) / len(kappa_temp_d)
            kappa_sim_avg_d = sum(kappa_sim_temp_d) / len(kappa_sim_temp_d)
            # Now determine the weighted sum.
            ws_d = ws_3_metrics(kappa_avg_d, kappa_sim_avg_d, clu_avg_d,
                                w1, w2, w3, r1, r2, r3, s1, s2, s3)
            # Add step here to break out if objectives are not changing.
            if ws_avg == ws_c and ws_avg == ws_d:
                break
            # Otherwise evaluate the result to close the bracket.
            # As we want to maximise, evaluation is as below:
            else:
                if ws_c > ws_d:
                    b = d
                else:
                    a = c
        # This section of code is run until the tolerance criteria are violated.
        # When this happens, the final evaluation is conducted.
        # Set the final value, input into the model.
        if ws_avg >= ws_c and ws_avg >= ws_d:
            rules[rule_name][1] = old_y1
            rules[rule_name][2] = old_y2
            set_lp_rule(project_file, fu_elem, lu_elem, y0, old_y1, old_y2, xe)
            pt_value.append(old_y1)
            kappa_log.append(kappa_avg)
            kappa_sim_log.append(kappa_sim_avg)
            clu_log.append(clu_avg)
            ws_log.append(ws_avg)
        # If the weighted sum for c is better, take that value. Insert into
        # model and log output.
        elif ws_c > ws_d:
            # First track the rule.
            final_value = c
            final_value_2 = final_value * 0.1
            set_lp_rule(project_file, fu_elem, lu_elem, y0, final_value,
                        final_value_2, xe)
            pt_value.append(final_value)
            rules[rule_name][1] = final_value
            rules[rule_name][2] = final_value_2
            # Now update the best set of objectives.
            kappa_avg = kappa_avg_c
            kappa_sim_avg = kappa_sim_avg_c
            clu_avg = clu_avg_c
            ws_avg = ws_c
            # Log the output.
            kappa_log.append(kappa_avg)
            kappa_sim_log.append(kappa_sim_avg)
            clu_log.append(clu_avg)
            ws_log.append(ws_avg)
        # If the weighted sum for d is better, take that value. Insert into
        # model and log output.
        else:
            # First track the rule.
            final_value = d
            final_value_2 = final_value * 0.1
            set_lp_rule(project_file, fu_elem, lu_elem, y0, final_value,
                        final_value_2, xe)
            pt_value.append(final_value)
            rules[rule_name][1] = final_value
            rules[rule_name][2] = final_value_2
            # Now update the best set of objectives.
            kappa_avg = kappa_avg_d
            kappa_sim_avg = kappa_sim_avg_d
            clu_avg = clu_avg_d
            ws_avg = ws_d
            # Log the output.
            kappa_log.append(kappa_avg)
            kappa_sim_log.append(kappa_sim_avg)
            clu_log.append(clu_avg)
            ws_log.append(ws_avg)
# This section runs until all interaction tails have been appended.
# -----------------------------------------------------------------------------
# Record all the output information.
# Track the end time of the calibration.
end = time.time()
# Determine the duration of the calibration method.
duration = end - start
duration_h = duration/3600

# Record the duration of calibration and number of iterations performed.
output_duration_file = (output_path + case_study + "\\Fine_cal_output\\" 
                        "duration.txt")
f = open(output_duration_file, "w")
f.write("duration: " + str(duration_h) + "\n" +
        "iterations:" + str(iterations_counter))
f.close()

# Write the output rules.
output_rules_file = output_path + case_study + "\\Rules\\final_rules.csv"
store = [0]*6
with open(output_rules_file, "w", newline='') as csv_file:
    writer = csv.writer(csv_file)
    # Write a header line to the file
    values = ["from", "to", "y0", "y1", "y2", "xe"]
    writer.writerow(values)
    # Now write the neighbourhood rules in the form from ... to ...
    for i in range(0, luc):
        for j in range(0, act):
            key = "from " + luc_names[i] + " to " + luc_names[j + pas]
            store[0] = luc_names[i]
            store[1] = luc_names[j + pas]
            store[2] = rules[key][0]
            store[3] = rules[key][1]
            store[4] = rules[key][2]
            store[5] = rules[key][3]
            writer.writerow(store)

# Write the output log.
log_output_file = (output_path + "\\" + case_study + "\\Fine_cal_output\\" +
                   case_study + "_fine_tuning_output.csv")
store = [0]*8
with open(log_output_file, "w", newline='') as csv_file:
    writer = csv.writer(csv_file)
    # Save and write a header line to the file
    values = ["From", "To", "Distance", "Value", "Fuzzy Kappa", "FKS", "AAWCE", "WS"]
    writer.writerow(values)
    # Now write the log
    for i in range(0, len(pt_value)):
        store[0] = luc_names[rule_from_tracker[i]]
        store[1] = luc_names[rule_to_tracker[i] + pas]
        store[2] = pt_tracker[i]
        store[3] = pt_value[i]
        store[4] = kappa_log[i]
        store[5] = kappa_sim_log[i]
        store[6] = clu_log[i]
        store[7] = ws_log[i]
        writer.writerow(store)

# Indicate completion with a beep.
import winsound
Freq = 2500 # Set Frequency To 2500 Hertz
Dur = 1000 # Set Duration To 1000 ms == 1 second
winsound.Beep(Freq,Dur)

# Completed!
