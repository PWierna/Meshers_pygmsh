###############################################################################
#  Script to create a structured two-dimensional mesh of an open-hole         #
# specimen with either linear or quadratic quadrilateral elements (via gmsh)  #
# in terms of given geometry and mesh parameters                              #
#                                                                             # 
#                                                        Author: Pablo Wierna #
#                                                             13-III-2023     #
###############################################################################

import gmsh
import numpy as np
import math

#·# Inputs --------------------------------------------------------------------
#·# Specimen Geometry parameters:
X0 = 0      #
Y0 = 0      #
Z0 = 0      #
total_width  = 500   #Specimen's grip (total) width 
hole_diam    = 250 #Hole diameter (Always located at the geom center of the specimen)
grip_length  = 250   #Specimen's total length
alpha_ratio  = 1.0  #Grip to Hole-Zone lengths ratio  
geom_type    = 'Whole' #Options: "Quarter","Half","Whole"

#·# Discretization parameters:
nelem_transv = 20#
nelem_diag   = 15#
nelem_long_holezone = 20
nelem_long_grip     = 10#
#------------------------------------------------------------------------------

     
def get_geometry_data( geometry_type ):
    
    #(I)-POINTS DEFINITION:
    npoints = 23
    pcoords = np.zeros([npoints,3])
    #idx_xy = np.array([0,1]) #dummy index for x and y coords
    
    #Center point: stays in [0,0,0]
    #Points 2 to 9 (on the hole):
    theta_vect = math.atan2(1.0,alpha_ratio) * np.array([0,1,0,-1,0,1,0,-1]) + math.pi/2 * np.array([0,0,1,2,2,2,3,4]) 
    pcoords[ np.arange(2,10)-1 , 0 ] = (hole_diam/2) * np.cos(theta_vect)
    pcoords[ np.arange(2,10)-1 , 1 ] = (hole_diam/2) * np.sin(theta_vect)
    
    #Points 10 to 17 (on the hole-zone boundary):
    radius_vect = (alpha_ratio*total_width/2) * np.array([1,0,0,0,1,0,0,0]) + (total_width*math.sqrt(1+alpha_ratio**2)/2) * np.array([0,1,0,1,0,1,0,1]) + (total_width/2) * np.array([0,0,1,0,0,0,1,0])
    pcoords[ np.arange(10,18)-1 , 0 ] = radius_vect * np.cos(theta_vect)
    pcoords[ np.arange(10,18)-1 , 1 ] = radius_vect * np.sin(theta_vect)
            
    #Points 18 to 23 (on the grip's boundaries):
    pcoords[ np.arange(18,24)-1 , 0 ] = (alpha_ratio*total_width/2 + grip_length) * np.array([1,1,-1,-1,-1,1])
    pcoords[ np.arange(18,24)-1 , 1 ] = (total_width/2) * np.array([0,1,1,0,-1,-1])
    
    #(II)-CIRCLE ARCS DEFINITION:
    ncirclearcs = 8
    #ca_ids      = np.zeros(ncirclearcs,dtype=int)
    ca_conect   = np.zeros([ncirclearcs,4],dtype=int) #Cols = [Startpt,Centerpt,Endpt,ndivisions] 
    #Connectivities:
    for i in range(1,8):
        ca_conect[i-1,np.array([0,1,2])] = np.array([i+1,1,i+2]) 
    ca_conect[7,np.array([0,1,2])] = np.array([8+1,1,2]) #8thCA
    #Divisions:
    ca_conect[ np.array([1,4,5,8])-1 , 3 ] = (nelem_transv/2) + 1  
    ca_conect[ np.array([2,3,6,7])-1 , 3 ] = (nelem_long_holezone/2) + 1   
    
    #(II)-LINES DEFINITION:
    nlines = 26
    ln_conect = np.zeros([nlines,3],dtype=int) #Cols = [Startpt,Endpt,ndivisions] 
    
    #Connectivities Lines 1 to 8:
    i = np.arange(1,8)
    ln_conect[ i-1 , 0 ] = i + 9
    ln_conect[ i-1 , 1 ] = i + 10
    i = 8
    ln_conect[ i-1 , 0 ] = i + 9
    ln_conect[ i-1 , 1 ] = 10
    
    #Connectivities Lines 9 to 16:
    i = np.arange(9,17)
    ln_conect[ i-1 , 0 ] = i - 7
    ln_conect[ i-1 , 1 ] = i + 1
    
    #Conectivities Lines 17 to 26:
    i = np.arange(17,27)
    ln_conect[ i-1 , 0 ] = np.array([10,14,18,19,13,20,21,22,17,23],dtype=int)
    ln_conect[ i-1 , 1 ] = np.array([18,21,19,11,20,21,22,15,23,18],dtype=int)
    
    #Lines Divisions:
    idx_transvln = np.array([1,4,5,8,19,22,23,26])  #Transversal lines
    idx_hzlongln = np.array([2,3,6,7])              #Hole-zone long lines
    idx_diagln   = np.arange(9,17)                  #Diagonal lines
    idx_glongln  = np.array([17,18,20,21,24,25])    #Grip-zone long lines
    ln_conect[ idx_transvln-1 , 2 ] = (nelem_transv/2) + 1
    ln_conect[ idx_hzlongln-1 , 2 ] = (nelem_long_holezone/2) + 1
    ln_conect[ idx_diagln-1 , 2 ]   = nelem_diag + 1
    ln_conect[ idx_glongln-1 , 2 ]  = nelem_long_grip + 1
    
        
    #(III)-CURVE LOOPS DEFINITIONS:
        
    #Curve-Loops conectivity: list where each element is a numpy array w/the conectivities of 
    #the respective curve loop. 2nd dim is type of the geom entity (0=line,1=circe_arc)
    #Cls 1 to 7:
    cl_conect = []
    for i in range(1,8):
        cl_conect.append({ 'geometry_types': np.array([  2 ,  1  , 1 ,   1   ]),
                           'entities_ids':   np.array([  i , i+8 , i , 8+i+1 ]),
                           'signs':          np.array([ -1 ,  1  , 1 ,  -1   ])
                         })
    #Cl 8:
    for i in range(8,9):
        cl_conect.append({ 'geometry_types': np.array([  2 ,  1  , 1 ,   1   ]),
                           'entities_ids':   np.array([  i , i+8 , i ,  8+1 ]),
                           'signs':          np.array([ -1 ,  1  , 1 ,  -1   ])
                         })
    #Cls 9 to 12:
    aux_conect = np.vstack(( np.array([-1,17,19,20]) , np.array([-4,21,22,-18]) , np.array([-5,18,23,24]) , np.array([-8,25,26,-17]) ))
    for i in range(0,4):        
        cl_conect.append({ 'geometry_types': np.ones(np.shape(aux_conect)[1],dtype='int'),
                           'entities_ids':   np.abs(aux_conect[i,:],dtype='int'),
                           'signs':          np.sign(aux_conect[i,:],dtype='int')
                         })
    #Additional CLs (13 to 16): 
    cl_conect.append({ 'geometry_types': np.array([ 1 , 1 , 1 ,  1 ,  2 ,  2 ]),
                       'entities_ids'  : np.array([ 9 , 1 , 2 , 11 ,  2 ,  1 ]),
                       'signs'         : np.array([ 1 , 1 , 1 , -1 , -1 , -1 ])})
     
    cl_conect.append({ 'geometry_types': np.array([  1 , 1 , 1 ,  1 ,  2 ,  2 ]),
                       'entities_ids'  : np.array([ 11 , 3 , 4 , 13 ,  4 ,  3 ]),
                       'signs'         : np.array([  1 , 1 , 1 , -1 , -1 , -1 ])})
    
    cl_conect.append({ 'geometry_types': np.array([  1 , 1 , 1 ,  1 ,  2 ,  2 ]),
                       'entities_ids'  : np.array([ 13 , 5 , 6 , 15 ,  6 ,  5 ]),
                       'signs'         : np.array([  1 , 1 , 1 , -1 , -1 , -1 ])})
    
    cl_conect.append({ 'geometry_types': np.array([  1 , 1 , 1 ,  1 ,  2 ,  2 ]),
                       'entities_ids'  : np.array([ 15 , 7 , 8 ,  9 ,  8 ,  7 ]),
                       'signs'         : np.array([  1 , 1 , 1 , -1 , -1 , -1 ])})
    
    #(IV)-RETRIEVE SURFACES TO BE ACTUALLY MESHED:
    if geometry_type == 'Quarter' or geom_type == 'quarter':
        sf_connect = np.array([1,2,9],dtype='int') 
    elif geometry_type == 'Half' or geom_type == 'half':
        sf_connect = np.array([1,2,7,8,9,12])
    elif geometry_type == 'Whole' or geom_type == 'whole':
        sf_connect = np.arange(1,13,dtype='int')
    elif geometry_type == 'Whole2' or geom_type == 'whole':
        sf_connect = np.arange(9,17,dtype='int')
    
    #(V)-BUILD OUTPUT DICT:
    output_geom_data = { "points" : pcoords ,
                         "circle_arcs" : ca_conect ,
                         "lines" : ln_conect ,
                         "curve_loops" : cl_conect ,
                         "surfaces" : sf_connect }
    
    return output_geom_data
        

    

#Build geometry data:    
geometrydata = get_geometry_data(geom_type)
        




geometrydata["points"][:,0] += X0
geometrydata["points"][:,1] += Y0
geometrydata["points"][:,2] += Z0


#·# Initialize Geometric Model and Mesh Algorithm 
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.model.add('OpenHole')


#-# Create entities:

#Points:
npoints = np.shape(geometrydata["points"])[0] #Get number of points
pt_ids  = np.zeros(npoints,dtype=int)         #Initialize IDs array
for pt in range(0,npoints):
    pt_ids[pt] = gmsh.model.geo.addPoint( geometrydata["points"][pt,0] , geometrydata["points"][pt,1] , geometrydata["points"][pt,2] )

#CircleArcs:
ncirclearcs = np.shape(geometrydata["circle_arcs"])[0] #Get number of circle arcs
ca_ids = np.zeros(ncirclearcs,dtype=int)               #Initialize IDs array
for ca in range(0,ncirclearcs):
    ca_ids[ca] = gmsh.model.geo.addCircleArc( geometrydata["circle_arcs"][ca,0] , geometrydata["circle_arcs"][ca,1] , geometrydata["circle_arcs"][ca,2])
    gmsh.model.geo.mesh.setTransfiniteCurve( ca_ids[ca] , geometrydata["circle_arcs"][ca,3] )

#Lines:
nlines = np.shape(geometrydata["lines"])[0] #Get number of lines
ln_ids = np.zeros(nlines,dtype=int)         #Initialize IDs array
for ln in range(0,nlines):
    ln_ids[ln] = gmsh.model.geo.addLine( geometrydata["lines"][ln,0] , geometrydata["lines"][ln,1] )
    gmsh.model.geo.mesh.setTransfiniteCurve( ln_ids[ln], geometrydata["lines"][ln,2] ) 
    
#CurveLoops:
ncloops = np.shape(geometrydata["curve_loops"])[0] #Get number of curve-loops
cl_ids  = np.zeros(ncloops,dtype=int)              #Initialize IDs array
for cl in range(0,ncloops): 
    
    #Initialize connectivities for the i-th curve loop:
    cl_connect_i = np.zeros( np.shape( geometrydata["curve_loops"][cl]['entities_ids'] ) , dtype='int' )
    
    #Get connectivities for the i-th curve loop:
    idx_lns = geometrydata["curve_loops"][cl]['geometry_types']==1 #Indexes for "line-type" curves
    idx_cas = geometrydata["curve_loops"][cl]['geometry_types']==2 #Indexes for "circle-arc-type" curves
    cl_connect_i[idx_lns] = ln_ids[ geometrydata["curve_loops"][cl]['entities_ids'][idx_lns]-1 ]
    cl_connect_i[idx_cas] = ca_ids[ geometrydata["curve_loops"][cl]['entities_ids'][idx_cas]-1 ]
    
    #Add i-th curve loop:
    cl_ids[cl] = gmsh.model.geo.addCurveLoop( cl_connect_i*geometrydata["curve_loops"][cl]['signs'] )


#Create surfaces to be actually meshed:
nsurfs = np.shape(geometrydata["surfaces"])[0] #Get number of surfaces for the mesh
sf_ids = np.zeros(nsurfs,dtype=int)            #Initialize IDs array
for sf in range(0,nsurfs): #try: in [2,3]:   
    sf_ids[sf] = gmsh.model.geo.addPlaneSurface( [ geometrydata["surfaces"][sf] ] ) 


#Assign surfaces a physical entity
gmsh.model.addPhysicalGroup( 2 , sf_ids )

#Synchronize model 
gmsh.model.geo.synchronize()

#Set transfinite surfaces & recombine
for sf in range(0,nsurfs):   
    gmsh.model.geo.mesh.setTransfiniteSurface(sf_ids[sf])
    gmsh.model.geo.mesh.setRecombine(2, sf_ids[sf])

#Synchronize model 
gmsh.model.geo.synchronize()

#Finalize meshing and run GUI
gmsh.option.setNumber("Mesh.RecombineAll", 2)
gmsh.model.mesh.generate(2)
gmsh.model.mesh.setOrder(2)
gmsh.option.setNumber('Mesh.SurfaceFaces', 1)
gmsh.option.setNumber('Mesh.Points', 1)
gmsh.write('open_hole2D.msh')
gmsh.fltk.run()
gmsh.finalize()



