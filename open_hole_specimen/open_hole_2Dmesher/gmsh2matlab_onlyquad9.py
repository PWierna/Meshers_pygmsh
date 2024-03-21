###############################################################################
#  Python script to create
#



import numpy as np                                                        
import meshio                                     
import sys                                          
                                                    
original_stdout = sys.stdout #Get Original standard output                    
                                                    
                                   
#·# PARSE INPUT DATA ----------------------------------------------------------
inp_file = str(input('\nEnter .msh file name/path (include .msh extension): '))    
inp_msh  = meshio.read(inp_file)    
ncoords   = inp_msh.points               #Get Nodal coordinates matrix   
cells_con = inp_msh.cells_dict["quad9"]  #Get Conectivities    
npe = 9   #Nodes per element in hexahedron    

#Report number of elements found:
nelem  = np.shape(cells_con)[0] #Number of elements
print('\n'+18*'-'+' {0} 9-NODED QUAD ELEMENTS DETECTED '.format(nelem)+18*'-'+'\n')


#·# BUILD CONNECTIVITIES MATRIX -----------------------------------------------
resortidx = [0,4,1,5,2,6,3,7,8]  #nodes resorting indexes for 9-node quad
true_con  = cells_con + 1        #pyth indices start in 0 & n°nod must start in 1
connectivities = np.zeros([nelem , npe+1]) #Add aditional column (for mat type)
connectivities[:,0]  = 1  #1st Column = Material ID
connectivities[:,1:] = true_con[:,resortidx]   #connectivities (w/Reordering)


#·# WRITE OUTPUT (.txt file)
print('\nWriting .m file: Connectivities_and_Coordinates_2D.m , wait...\n')
with open('Connectivities_and_Coordinates_2D.m','w') as f:
    sys.stdout = f #Change the standard output to the file created

    np.set_printoptions( precision=20, threshold=sys.maxsize , suppress=1 , linewidth=1e3)
    
    print('\nMACRO_MODEL.Conectivity = ...')
    print(connectivities)
    print(';\n')
    
    print('\nMACRO_MODEL.Coordinates = ...')
    print(ncoords)
    print(';\n')
    
    # print("MACRO_MODEL.Conectivity = [ (1:size(MACRO_MODEL.Conectivity,1))' MACRO_MODEL.Conectivity ];\n")
    #print("MACRO_MODEL.Coordinates = [ (1:size(MACRO_MODEL.Coordinates,1))' MACRO_MODEL.Coordinates ];\n")

sys.stdout = original_stdout
print('\nDone writing Connectivities_and_Coordinates_2D.m\n')


