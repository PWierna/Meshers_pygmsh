import numpy as np
import meshio
import sys

original_stdout = sys.stdout #Original standard output


#·# Get & break down input data ----------------------------------------------------------
inputfile= str(input('\nEnter .msh file name/path (include .msh extension): '))
inputmsh = meshio.read(inputfile)
n_coords = inputmsh.points                    #Nodal coordinates



connect_h27  = inputmsh.cells_dict["hexahedron27"] + 1 #Conectivities for 27-noded-hexahedrons
nelem        = np.shape(connect_h27)[0]
npe          = np.shape(connect_h27)[1]

#Assemble connectivities matrix:
connectivities = np.zeros([nelem,npe+1])
connectivities[:,0] = (inputmsh.get_cell_data("gmsh:physical","hexahedron27")) + 1 #Layers type
    
#Reordering:
idx_reord = [2,3,0,1,6,7,4,5,10,11,8,9,18,19,16,17,14,15,12,13,24,23,20,22,21,25,26]
connectivities[:,1:] = connect_h27[:,idx_reord] 

#Check all elements have a material set:
if (all(connectivities[:,0]==0)==False) != True:
    print('Some elements were not assigned a material set. Check')
    


#·# WRITE OUTPUT (.m file)
print('\nWriting MATLAB .m file: Connectivities_and_Coordinates_3D.m , wait...\n')
with open('Connectivities_and_Coordinates_3D.m','w') as f:
    sys.stdout = f #Change the standard output to the file created

    np.set_printoptions( precision=20, threshold=sys.maxsize , suppress=1 , linewidth=1e3)
    
    print('\nMODEL.Conectivity = ...')
    print(connectivities)
    print(';\n')
    
    print('\nMODEL.Coordinates = ...')
    print(n_coords)
    print(';\n')
    #print("MODEL.Conectivity = [ (1:size(MODEL.Conectivity,1))' MODEL.Conectivity ];\n")
    #print("MODEL.Coordinates = [ (1:size(MODEL.Coordinates,1))' MODEL.Coordinates ];\n")

sys.stdout = original_stdout
print('\nDone writing Connectivities_and_Coordinates_3D.m\n')
