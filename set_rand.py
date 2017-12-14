"""
Fix the random number seed in the Metronamica project file to a specified
value.
"""

import xml.etree.ElementTree as ET


def set_rand(project_path, rseed):

    source = open(project_path)
    tree = ET.parse(source)
    root = tree.getroot()

    spline_fixed = root[2][1][3][0][0][10][0][0]
    spline_value = root[2][1][3][0][0][10][0][1]
    
    fixed = 1

    spline_fixed.text = str(fixed)
    spline_value.text = str(rseed)
    tree.write(project_path)
