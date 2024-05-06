#
# Copyright ANSYS, Inc. Unauthorized use, distribution, or duplication is prohibited.
#
import os
import sys

def addSyCDLLs():
    if sys.platform.startswith("win"):
        for p in os.environ["PYTHON_DLL_PATH"].split(os.pathsep):
            try:
                os.add_dll_directory(p)
            except (FileNotFoundError, OSError):
                pass  # skip any paths that don't exist
                
addSyCDLLs()
                
from pyExt import SystemCouplingParticipant as sysc
import numpy as np



'''
Mesh representation:
  10 nodes
  4 quad elements
  Rectangle, aligned with XY-plane,
  stretching from the origin to point (1, 0.1, 0)
(0,0.1,0)                    (1,0.1,0)
         6----7----8----9----10
         |    |    |    |    |
         |    |    |    |    |
         1----2----3----4----5
(0,0,0)                       (1,0,0)
'''

''' Array of node coordinates in split array format. '''
nodeCoordsX = np.array([0, 0.25, 0.5, 0.75, 1, 0, 0.25, 0.5, 0.75, 1],
    dtype = np.float64)
nodeCoordsY = np.array([0, 0, 0, 0, 0, 0.1, 0.1, 0.1, 0.1, 0.1],
    dtype = np.float64)
nodeCoordsZ = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    dtype = np.float64)

''' Array of element node counts. '''
elemNodeCounts = np.array([4, 4, 4, 4])

''' Array of node ids for each element. '''
elemNodeIds = np.array(
   [1, 2, 7, 6,
    2, 3, 8, 7,
    3, 4, 9, 8,
    4, 5, 10, 9])

''' Array of temperature values (defined at nodes). '''
nodeTemperature = np.array([300.0] * 10)

'''
Array of heat rate values (defined at elements)
values are in Watts - total of 20 [W] generated
'''
elementHeatRate = np.array([5.0] * 4)

''' Function that returns surface mesh, given a region. '''
def getSurfaceMesh(regionName):
    if regionName == "FSI":
        nodes = sysc.NodeData(sysc.OutputVectorData(nodeCoordsX, nodeCoordsY, nodeCoordsZ))
        faces = sysc.FaceData(
            sysc.ElementNodeCountData(sysc.OutputIntegerData(elemNodeCounts)),
            sysc.ElementNodeConnectivityData(sysc.OutputIntegerData(elemNodeIds)))
        return sysc.SurfaceMesh(nodes, faces)
    errorMessage = "getSurfaceMesh error: "
    errorMessage += "Unknown region " + regionName
    raise RuntimeError(errorMessage)

'''
Function that returns scalar input variable values,
given a region and a variable.
'''
def getInputScalarData(regionName, variableName):
    if regionName == "FSI" and variableName == "Temperature":
        return sysc.InputScalarData(nodeTemperature)
    errorMessage = "getInputScalarData error: "
    errorMessage += "unknown variable " + variableName
    errorMessage += " on region " + regionName
    raise RuntimeError(errorMessage)

'''
Function that returns scalar output variable values,
given a region and a variable.
'''
def getOutputScalarData(regionName, variableName):
    if regionName == "FSI" and variableName == "Heat Rate":
        return sysc.OutputScalarData(elementHeatRate)
    errorMessage = "getOutputScalarData error: "
    errorMessage += "unknown variable " + variableName
    errorMessage += " on region " + regionName
    raise RuntimeError(errorMessage)

'''
This is the mock solver script.
It requires the following command line arguments to get
System Coupling host (string), port (integer),
and name (string) given to this participant:
--schost <host>
--scport <port>
--scname <name>

This program is going to run as part of the overall coupled
analysis. System Coupling is going to start this program
 automatically, so it will pass in the appropriate arguments as needed.
'''

# Parse and print input arguments (host, port, participant name).
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--schost', type=str, default="")
parser.add_argument('--scport', type=int, default=0)
parser.add_argument('--scname', type=str, default="")
parser.add_argument('--scsetup', default=False, action='store_true')
args = parser.parse_args()

print("Host: " + args.schost)
print("Port: " + str(args.scport))
print("Name: " + args.scname)

buildInfo = "Channel Flow Mock Solver, language: Python"

'''
Connection to System Coupling.

Create SystemCoupling object - this will establish the
connection to System Coupling.
'''
sc = sysc.SystemCoupling(args.schost, args.scport, args.scname, buildInfo)

''' Coupled analysis mode. '''

'''
Heavyweight data access and initialization.

Register access to mesh and variable data.
NOTE: The participant solver retains ownership
of mesh and variable field data arrays and is thus
responsible for allocating and deallocating memory
for those arrays. System Coupling only gets
read and/or write access to those arrays.
'''

'''
Register surface mesh access by providing pointer
to a function that returns SurfaceMesh object corresponding
to the requested region.
'''
sc.registerSurfaceMeshAccess(getSurfaceMesh)

'''
Register variable data access by providing pointer
to a function that provides access to variable data for
a region.
'''

''' Output scalar variable access. '''
sc.registerInputScalarDataAccess(getInputScalarData)

''' Input scalar variable access. '''
sc.registerOutputScalarDataAccess(getOutputScalarData)

''' Coupled analysis initialization. '''
sc.initializeAnalysis()

''' Coupled analysis loop. '''
while sc.doIteration():
    '''
    Inputs update.
    All input variables will be up-to-date after this call.
    '''
    sc.updateInputs()
    ''' 
    Solver iterations can be peformed here.
    Input variable values can be consumed.
    Output variable values must be generated.
    ...
    '''
    '''
    Outputs update.
    Provide convergence state to System Coupling,
    and notify that outputs can be consumed.
    '''
    sc.updateOutputs(sysc.Converged)

'''
Coupled analysis shutdown.
'''
sc.disconnect()
