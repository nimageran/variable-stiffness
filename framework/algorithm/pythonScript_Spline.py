# -*- coding: mbcs -*-
from __future__ import division 
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import math
import numpy
import sys
from abaqus import *
from abaqusConstants import *
import __main__
import regionToolset
import displayGroupMdbToolset as dgm
import xyPlot
import displayGroupOdbToolset as dgo
import re
import mesh
import os 
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
import fileinput
from scipy import interpolate

FittingMethod='Spline' # Choose method : 'Spline' or 'Linear'


jj=0
################################################
# Counting the number of files in the directory (nFiles)
################################################

path = "C:/optLaminatedComp/JOB"
dirs = os.listdir( path )

# This would print all the files and directories
nFiles=0 # Number of files existed in the directory path
for file in dirs:
        nFiles=nFiles+1
#print(nFiles)
################################################
# Making folder for the current Model: 
################################################
os.mkdir('C:\\optLaminatedComp\\JOB\\%d'%(jj+1+nFiles))

################################################
# Changing the default working directory to an arbitrary working directory :
################################################
os.chdir(r"C:\\optLaminatedComp\\JOB\\%d"%(jj+1+nFiles))

##################################### Regenerating  : 
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
meshTechnique=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
referenceRepresentation=OFF)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON, 
optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
meshTechnique=ON)

###############################################
#####            Creating Part            #####
#####(Units used in Abaqus [N] and [mm])  #####
###############################################

NS=134 #Number of Total narrow segments (must be even)
D=457.200012# L=D=18 in
L=D*1 # D [m] L/D=L/R*0.5 
thicknessOfEachPly=0.127



s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
sheetSize=1000.0)
g, v1, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.ConstructionLine(point1=(0.0, -500.0), point2=(0.0, 500.0))
s.FixedConstraint(entity=g[2])
s.Line(point1=(D/2, 0.0), point2=(D/2, L))
s.VerticalConstraint(entity=g[3], addUndoState=False)
p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=THREE_D, 
type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-1']
p.BaseShellRevolve(sketch=s, angle=360.0, flipRevolveDirection=OFF)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']


##########################################################
################################### Creating Assembly :
##########################################################
a1 = mdb.models['Model-1'].rootAssembly
a1.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Part-1']
a1.Instance(name='Part-1-1', part=p, dependent=OFF)

############################################################
##################################### Creating Step(Buckling):
############################################################
mdb.models['Model-1'].BuckleStep(name='Step-1', previous='Initial', numEigen=1, 
        eigensolver=LANCZOS, minEigen=0.0, blockSize=DEFAULT, 
        maxBlocks=DEFAULT)

############################################################
##################################### Material Definition :
############################################################
nCell=1
for l in range(nCell):
    mdb.models['Model-1'].Material(name='Material_%d'%(l+1))
    mdb.models['Model-1'].materials['Material_%d'%(l+1)].Elastic(type=ENGINEERING_CONSTANTS, table=((134000.0, 7710.0, 7710.0, 0.301, 
0.301, 0.396, 4310, 4310, 2760), ))

#######################################################################
####################          Partitioning           ##################
#######################################################################

## Defining Datum Planes for Partitioning:

# Creating ZX Datum Plane:
p = mdb.models['Model-1'].parts['Part-1']
p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.0)

# Creating line of Rotation for the rest of Datum Planes:
p = mdb.models['Model-1'].parts['Part-1']
p.DatumAxisByPrincipalAxis(principalAxis=YAXIS)

# Showing the Part Model:
p = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)

# Creating the rest of Datum Planes by rotating the First Plane:
NDP=(NS/2)-1 #Number of Datum Planes created by the Rotation of the XY Plane 
NDP=int(NDP)
for i in range(NDP):
	p = mdb.models['Model-1'].parts['Part-1']
	d1 = p.datums
	myAngle=(360/NS)*(i+1)
	p.DatumPlaneByRotation(plane=d1[2], axis=d1[3], angle=myAngle)
	
# Partitioning:	
NDPaLL=int((NS/2)+1) # Number of All Datum Planes
for i in range(NDPaLL): 
	if i==1: 
		continue  # Datum index number 3 is the Datum Axis and 

	p = mdb.models['Model-1'].parts['Part-1']
	d1 = p.datums
	f = p.faces
	pickedFaces=f.getByBoundingBox(-10000 ,-10000 ,-10000 ,10000 ,10000 ,10000)               

	p.PartitionFaceByDatumPlane(datumPlane=d1[i+2], faces=pickedFaces) 

	
##################################### Regenerating  : 
a = mdb.models['Model-1'].rootAssembly
a.regenerate()
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF, 
	interactions=ON, constraints=ON, connectors=ON, engineeringFeatures=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
	meshTechnique=OFF)
session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])	

#############################################################################################
#############################################################################################
##########      Creating MPC Constraint for the bottom and top circular Edges      ##########
#############################################################################################
#############################################################################################


####################################################################################
#### Creating a Set for the bottom edge to define an MPC constrait on that edge ####
#######################################################################3############

a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Part-1-1'].edges
# Following is the conversion between polar coordinates of 
# a point in the Top edges into Cartesian Cordinate:
# (Not that for the Top edge, z=L and the radius is unique for all 
#the edge segments)
rho=D/2
# But phi is different:
phi=[]
for i in range(NS):
    phi.append((360/(2*NS))+(i*(360/NS)))

# Therefore, x and y:
x=[]
for i in range(NS):
    x.append(rho * numpy.cos(numpy.deg2rad(phi[i])))

z=[]
for i in range(NS):
    z.append(rho * numpy.sin(numpy.deg2rad(phi[i])))


lvalues = []
for item in range(NS):
    lvalues.append(((x[item], 0, z[item]),))
aaa = {'coordinates': lvalues}

BOTTOMEDGE = e1.findAt(*aaa['coordinates'])
BOTTOMEDGE=a.Set(edges=BOTTOMEDGE, name='BOTTOMEDGE')


# Creating a Reference Point for the Center of the BOTTOM Circular Edge :
a = mdb.models['Model-1'].rootAssembly
myRefPoint=a.ReferencePoint(point=(0, 0, 0))

# Creating a Set for the above Reference Point :
RP_id=myRefPoint.id # # Creating a Set for Reference Point Needs the id of Reference Point, so we Create it with ((id)) syntax .
refpoint1= (mdb.models['Model-1'].rootAssembly.referencePoints[RP_id], ) # Creating a Set for Reference Point Needs a Tuple, So We Create a Tuple
a = mdb.models['Model-1'].rootAssembly
region1REFPOINTBOTTOM=a.Set(referencePoints=refpoint1, name='REFPOINTSET_BOTTOM')

# Creating MPC Constraint for the BOTTOM Circular Edge :
mdb.models['Model-1'].MultipointConstraint(name='Constraint-1', 
        controlPoint=region1REFPOINTBOTTOM, surface=BOTTOMEDGE, mpcType=BEAM_MPC, 
        userMode=DOF_MODE_MPC, userType=0, csys=None)

print('Bottom MPC Constraint Has Been Created')
###################################################################################
###### Creating a Set for the top edge to define an MPC constrait on that edge ####
#########################################################################3#########
a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Part-1-1'].edges
# Following is the conversion between polar coordinates of 
# a point in the Top edges into Cartesian Cordinate:
# (Not that for the Top edge, z=L and the radius is unique for all 
#the edge segments)
rho=D/2
# But phi is different:
phi=[]
for i in range(NS):
    phi.append((360/(2*NS))+(i*(360/NS)))

# Therefore, x and y:
x=[]
for i in range(NS):
    x.append(rho * numpy.cos(numpy.deg2rad(phi[i])))

z=[]
for i in range(NS):
    z.append(rho * numpy.sin(numpy.deg2rad(phi[i])))

lvalues = []

for item in range(NS):
    lvalues.append(((x[item], L, z[item]),))
aaa = {'coordinates': lvalues}

TOPEDGE = e1.findAt(*aaa['coordinates'])
TOPEDGE=a.Set(edges=TOPEDGE, name='TOPEDGE')


# Creating a Reference Point for the Center of the Top Circular Edge :
a = mdb.models['Model-1'].rootAssembly
myRefPoint=a.ReferencePoint(point=(0, L, 0))

# Creating a Set for the above Reference Point :
RP_id=myRefPoint.id # # Creating a Set for Reference Point Needs the id of Reference Point, so we Create it with ((id)) syntax .
refpoint1= (mdb.models['Model-1'].rootAssembly.referencePoints[RP_id], ) # Creating a Set for Reference Point Needs a Tuple, So We Create a Tuple
a = mdb.models['Model-1'].rootAssembly
region1REFPOINTTOP=a.Set(referencePoints=refpoint1, name='REFPOINTSET_TOP')

# Creating MPC Constraint for the Top Circular Edge :
mdb.models['Model-1'].MultipointConstraint(name='Constraint-2', 
        controlPoint=region1REFPOINTTOP, surface=TOPEDGE, mpcType=BEAM_MPC, 
        userMode=DOF_MODE_MPC, userType=0, csys=None)

print('TOP MPC Constraint Has Been Created')


## Turning the datume planes off:
session.viewports['Viewport: 1'].assemblyDisplay.geometryOptions.setValues(
        datumPlanes=OFF)
######################################################################
#######    Applying Moment on the Top Reference Point (Z=L)     ######
######################################################################
mdb.models['Model-1'].Moment(name='Load-1', createStepName='Step-1', 
        region=region1REFPOINTBOTTOM, cm3=-1, distributionType=UNIFORM, field='', 
        localCsys=None) # Moment load = 1 [kN.m] or 1000 [N.m]
mdb.models['Model-1'].Moment(name='Load-2', createStepName='Step-1', 
        region=region1REFPOINTTOP, cm3=1, distributionType=UNIFORM, field='', 
        localCsys=None) # Moment load = 1 [kN.m] or 1000 [N.m]

############################################################
################        Asigning Mesh Controls #############
############################################################

###############################################################
################          Meshing                ##############
###############################################################

############################ Seeding the Top and Bottom Edge:
if NS<60:
	nSeed=int(360/(NS*6))
else: 
	nSeed=1

a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Part-1-1'].edges

BOTTOMEDGE=mdb.models['Model-1'].rootAssembly.sets['BOTTOMEDGE'].edges
a.seedEdgeByNumber(edges=BOTTOMEDGE, number=nSeed, constraint=FINER)

TOPEDGE=mdb.models['Model-1'].rootAssembly.sets['TOPEDGE'].edges
a.seedEdgeByNumber(edges=TOPEDGE, number=nSeed, constraint=FINER)

##########################################################################
#################         Creating Boundary Conditions         ###########
##########################################################################
# Bottom Edge
a = mdb.models['Model-1'].rootAssembly
region = a.sets['REFPOINTSET_BOTTOM']
mdb.models['Model-1'].DisplacementBC(name='BC-REFPOINTSET_BOTTOM', createStepName='Step-1', 
        region=region, u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=UNSET, 
        amplitude=UNSET, buckleCase=PERTURBATION_AND_BUCKLING, fixed=OFF, 
        distributionType=UNIFORM, fieldName='', localCsys=None)
# Top Edge
a = mdb.models['Model-1'].rootAssembly
region = a.sets['REFPOINTSET_TOP']
mdb.models['Model-1'].DisplacementBC(name='BC-REFPOINTSET_TOP', createStepName='Step-1', 
        region=region, u1=0.0, u2=UNSET, u3=0.0, ur1=0.0, ur2=0.0, ur3=UNSET, 
        amplitude=UNSET, buckleCase=PERTURBATION_AND_BUCKLING, fixed=OFF, 
        distributionType=UNIFORM, fieldName='', localCsys=None)

##########################################################
###############           Layup           ################
##########################################################

###############################  Creating Set for the Segments:
# Similar to the case of Creating a Set for the Top edge, we have
#(The only difference is the z=L/2 (or other values except fo L or 0)):
rho=D/2
phi=[]
for i in range(NS):
    phi.append((360/(2*NS))+(i*(360/NS)))

x=[]
for i in range(NS):
    x.append(rho * numpy.cos(numpy.deg2rad(phi[i])))

z=[]
for i in range(NS):
    z.append(rho * numpy.sin(numpy.deg2rad(phi[i])))
	
for i in range(NS):
	p = mdb.models['Model-1'].parts['Part-1']
	f = p.faces
	faces = f.findAt(((x[i], L/2, z[i]), ))
	p.Set(faces=faces, name='narrowSegment_%d'%(i+1))

############################ Meshing:
a = mdb.models['Model-1'].rootAssembly
partInstances =(a.instances['Part-1-1'], )
a1 = mdb.models['Model-1'].rootAssembly
a1.regenerate()
a.generateMesh(regions=partInstances)
######################### S8R5
for i in range(NS):
        elemType1 = mesh.ElemType(elemCode=S8R5, elemLibrary=STANDARD)
        elemType2 = mesh.ElemType(elemCode=STRI65, elemLibrary=STANDARD)
        a = mdb.models['Model-1'].rootAssembly
        f1 = a.instances['Part-1-1'].faces
        faces1 = f1.findAt(((x[i], L/2, z[i]), ))
        pickedRegions =(faces1, )
        a.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
###################################  Layup:
######################## Reading the orientations from TEXT2 file:
#
f=open('C:\\optLaminatedComp\\proposed_T.txt')
T = []
for l in f:
	row = l.split()
	print(row)
	T=row
T = list(map(float, T))
print(T)
#
        
myOrientationsList = [ [] for e in range(NS) ]

## Linear 
m=7 # Number of variables
n=10 # Number of narrow bands in each wide segment

print('Total number of sections = ',2*(m+(n)*(m-1)))
theta1=[]

n=n+1
for k in range(0,n*(m-1)): 

        alpha_k=k
        i=math.floor(k/n)
        i=int(i)
        print(i)
        alpha_i=i*n
        alpha_ip1=(i+1)*n
        TT=T[i]+((alpha_k-alpha_i)/(alpha_ip1-alpha_i))*(T[i+1]-T[i])
        theta1.append(TT)
theta1.append(T[m-1])
###

## Spline

x_Points=numpy.linspace(0.0, 180.0, num=m)
x=numpy.linspace(0.0, 180.0, num=(m+(n-1)*(m-1)))
tck = interpolate.splrep(x_Points, T)
theta2=interpolate.splev(x, tck)
###
b ={
    'Linear'  : theta1,
    'Spline' : theta2
}
theta=b.get(FittingMethod, -1)
print(theta)
###
Flip_theta=theta[::-1] # Flipping the theta1
theta.extend(Flip_theta)
#print(theta1)
for i in range(NS):    
    myOrientationsList[i].extend([0,theta[i],90,-theta[i],-theta[i],90,theta[i],0])
print(myOrientationsList)

##################################

nLayers=8 # Number of layers in each segments
sectionLayer = [ [] for i in range(NS) ]
for i in range(NS): 
	for j in range(nLayers): 
                sectionLayer[i].extend([ section.SectionLayer(material='Material_1', thickness=0.127, 
                orientAngle=myOrientationsList[i][j], numIntPts=3, plyName='Ply-%d'%(j+1))] )
for i in range(NS): 
	for j in range(nLayers): 
                mdb.models['Model-1'].CompositeShellSection(name='narrowSegment_%d'%(i+1), 
                preIntegrate=OFF, idealization=NO_IDEALIZATION, symmetric=True, 
                thicknessType=UNIFORM, poissonDefinition=DEFAULT, 
                thicknessModulus=None, temperature=GRADIENT, useDensity=OFF, 
                integrationRule=SIMPSON, layup=(sectionLayer[i][0], sectionLayer[i][1], 
                sectionLayer[i][2], sectionLayer[i][3], sectionLayer[i][4], sectionLayer[i][5], 
                sectionLayer[i][6], sectionLayer[i][7], ))
                ###


p = mdb.models['Model-1'].parts['Part-1']
for i in range(NS):
        region = p.sets['narrowSegment_%d'%(i+1)]
        p = mdb.models['Model-1'].parts['Part-1']
        p.SectionAssignment(region=region, sectionName='narrowSegment_%d'%(i+1), offset=0.0, 
        offsetType=SINGLE_VALUE, offsetField='', 
        thicknessAssignment=FROM_SECTION)

#######################################################################################################
################################ Defining cylindrical coordinate systems for each narrow segments   ###
#######################################################################################################

##############################
p = mdb.models['Model-1'].parts['Part-1']
##################### Creating a CYLINDRICAL coordinate system at the origin (where we previously deinded refpint1):
p.DatumCsysByThreePoints(origin=refpoint1, point1=(0.0, 0.0, D/4), point2=(D/4, 0.0, 0.0), 
        name='myFirst_CSYS', coordSysType=CYLINDRICAL)
#################### Creating the CYLINDRICAL reference coordinate system based on the above coordinate system:

d = p.datums
Datume_Keys=d.keys() # Datume_Keys[-1] gives the key of last datum pint

p.DatumCsysByOffset(datumCoordSys=d[Datume_Keys[-1]], name='mySecond_CSYS', 
        coordSysType=CYLINDRICAL, vector=(D/2, 0.0, 0.0))
###############################
d = p.datums
Datume_Keys=d.keys() # Datume_Keys[-1] gives the key of last datum pint

p = mdb.models['Model-1'].parts['Part-1']
for i in range(NS):
        region = p.sets['narrowSegment_%d'%(i+1)]
        orientation = mdb.models['Model-1'].parts['Part-1'].datums[Datume_Keys[-1]]
        mdb.models['Model-1'].parts['Part-1'].MaterialOrientation(region=region, 
                orientationType=SYSTEM, axis=AXIS_2, localCsys=orientation, 
                fieldName='', additionalRotationType=ROTATION_ANGLE, 
                additionalRotationField='', angle=0.0)
##########################################################################################
##################################### Creating and submiting the current job     #########
##########################################################################################

jj=0
mdb.Job(name='Job-%d'%(jj+1+nFiles), model='Model-%d'%(jj+1), description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=False, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, 
        numGPUs=0)

#################################################################
############# Saving the completed job in arbitrary folder path :
#################################################################
mdb.saveAs(pathName='C:/optLaminatedComp/JOB/%d'%(jj+1+nFiles) + '/Job-%d'%(jj+1+nFiles))

###################################################
############# Submiting the current job :       ###
###################################################
mdb.jobs['Job-%d'%(jj+1+nFiles)].submit(consistencyChecking=OFF)

###########################################################
############# Waiting for completion of the current job  ##
###########################################################
mdb.jobs['Job-%d'%(jj+1+nFiles)].waitForCompletion()

#############################################################
###################   Result extraction : Buckling Loads  ###
#############################################################
import odbAccess

for i in range(1,101): # Number of desired eigenvalues are 100
        a1=odbAccess.openOdb('Job-%d.odb'%(jj+1+nFiles))
        a1.steps['Step-1'].frames[i].mode
        a1.steps['Step-1'].frames[i].description
        f1=a1.steps['Step-1'].frames[i].description
        First_Buckling_Load=float(f1[28:48])
        if numpy.sign(First_Buckling_Load) == -1:
                continue        
        else:
                modeNumber=i
                break

########################################################################
########################################################################
##r#####      Witing modeNumber for every simulations       ############
########         into the associated folder of that Job     ############
########################################################################
########################################################################
file_id1=open('C:/optLaminatedComp/JOB/%d'%(jj+1+nFiles) + '/modeNumber.txt','w')
file_id1.write(str(modeNumber))
file_id1.close()
########################################################################
########################################################################
########################################################################
file_id1=open('C:/optLaminatedComp/JOB/%d'%(jj+1+nFiles) + '/bucklingMOMENT.txt','w')
file_id1.write(str(First_Buckling_Load))
file_id1.close()

########################################################################
########################################################################
##r#####      Witing Ply Orientations for every simulations ############
########         into the associated folder of that Job     ############
########################################################################
########################################################################
file_id1=open('C:/optLaminatedComp/JOB/%d'%(jj+1+nFiles) + '/thisJobOrientations.txt','w')
file_id1.write(str(myOrientationsList))
file_id1.close()
########################################################################
########################################################################
##r#####                         Writing T                   ###########
########         into the associated folder of that Job      ###########
########################################################################
########################################################################
file_id1=open('C:/optLaminatedComp/JOB/%d'%(jj+1+nFiles) + '/T.txt','w')
file_id1.write(str(T))
file_id1.close()

