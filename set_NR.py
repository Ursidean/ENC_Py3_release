'''Adjust neighbourhood rules in the Metronamica project file.'''


#Modules
import numpy as np
import xml.etree.ElementTree as ET

#Function
def set_lp_rule(project_path,fu_elem,lu_elem,y0,y1,y2,xe):
    #y0 is the influence at a distance of 0
    #y1 is the influence at a distance of 1
    #y2 is the influence at a distance of 2
    #xe is the distance for an influence of 0
    source=open(project_path)
    tree=ET.parse(source)                                                 #Parse georpoject file as xml type data structure for modification
    root=tree.getroot()
    #Access specified neighbourhood rule
    spline=root[2][1][3][0][0][1][0][1][0][0][fu_elem][0][lu_elem][0]
    #Store array of values
    arrange=np.array([[0,y0],[1,y1],[2,y2],[xe,0]])
    #Find length of arrange
    size4=len(arrange)
    #Remove current neighbourhood rule values
    for point in spline.findall('point'):
        spline.remove(point)
    #Dictionary to store inputs for neighbourhood rules
    inputs={}                                                             
    for i in range(0,size4):
        key='new_value'+str(i)
        inputs[key]={'y':str(arrange[i,1]),'x':str(arrange[i,0])}
    #Insert points, commit to case study    
    for i in range(0,size4):                                              
        key='new_value'+str(i)
        child=ET.SubElement(spline, 'point', attrib=inputs[key]) 
    
    tree.write(project_path)
    source.close()