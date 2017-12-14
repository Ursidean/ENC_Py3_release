"""
This module, run_metro, runs Metronamica geoproject files
"""

import os
import subprocess


def run_metro(project_file, log_file, working_directory, geo_cmd):
    # Move to the working directory.
    os.chdir(working_directory)
    # Reset the project file to the first year.
    p1 = subprocess.Popen(
        [geo_cmd, "--Reset", project_file],
        shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )
    for line in p1.stdout.readlines():
        print(line)
    retval = p1.wait()
    # Run the project file, logging the output map for analysis.
    p2 = subprocess.Popen(
        [geo_cmd, "--Run", "--LogSettings", log_file, project_file],
        shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )
    for line in p2.stdout.readlines():
        print(line)
    retval = p2.wait()
