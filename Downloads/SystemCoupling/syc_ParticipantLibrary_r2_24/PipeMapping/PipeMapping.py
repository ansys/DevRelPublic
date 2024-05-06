'''
Copyright ANSYS, Inc. Unauthorized use, distribution, or duplication is prohibited.
'''

'''
This is the pipe mapping script.

It generates two different meshes, each for a region representing a 
pipe geometry. It then performs non-conformal mesh mapping between 
these two regions.

The mesh refinement level, pipe offset, pipe radius, etc. can be 
controlled by the command line arguments. The resulting arrays can 
be written for further post-processing.
'''
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

from PipeMeshGenerator import getMeshProperties
from PipeMeshGenerator import getNodeCoords
from PipeMeshGenerator import getElemNodeCounts
from PipeMeshGenerator import getElemNodeIds
from PipeMeshGenerator import getSolutionData
from PipeMeshGenerator import initializeMesh

from pyExt import SystemCouplingParticipant as sysc
import time

def getSurfaceMesh(regionName):
    ''' Function that returns surface mesh, given a region. '''
    ncs = sysc.OutputVectorData(getNodeCoords(regionName))
    encs= sysc.OutputIntegerData(getElemNodeCounts(regionName))
    enis= sysc.OutputIntegerData(getElemNodeIds(regionName))
    nodes = sysc.NodeData(ncs)
    faces = sysc.FaceData(
        sysc.ElementNodeCountData(encs),
        sysc.ElementNodeConnectivityData(enis))
    return sysc.SurfaceMesh(nodes, faces)

def getInputScalarData(regionName, variableName):
    '''
    Function that returns scalar input variable values,
    given a region and a variable.
    '''
    return sysc.InputScalarData(getSolutionData(regionName, variableName))

def getOutputScalarData(regionName, variableName):
    '''
    Function that returns scalar output variable values,
    given a region and a variable.
    '''
    return sysc.OutputScalarData(getSolutionData(regionName, variableName))

def printMaxError(regionName):
    ''' Function to calculate and print max error on a given region. '''
    expectedVar = "linear1" if regionName == "quad" else "linear2"
    testVar = "linear1" if expectedVar == "linear2" else "linear2"
    expectedValues = getSolutionData(regionName, expectedVar)
    testValues = getSolutionData(regionName, testVar)
    numVals = len(testValues)
    maxError = 0.0
    for n in range(numVals):
        testValue = testValues[n]
        expectedValue = expectedValues[n]
        maxError = max(maxError, abs(expectedValue - testValue))
    print(f"Maximum error on {regionName}: {maxError:.4e}")

'''
Parse the command line arguments.
'''
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--quad', type=int, default=20, help='Refinement level for quad region.')
parser.add_argument('--tri', type=int, default=30, help='Refinement level for tri region.')
parser.add_argument('--offset', type=float, default=0.0, help='Axial offset for quad region.')
parser.add_argument('--qradius', type=float, default=0.05, help='Pipe radius for quad region.')
parser.add_argument('--tradius', type=float, default=0.05, help='Pipe radius for tri region.')
parser.add_argument('--write', type=int, default=1, help='Flag to turn on/off .dat file writing.')
'''
Add overlap layer for consistency 
# (it has no effect here, since this script can only run in serial).
'''
parser.add_argument('--overlap', type=int, default=0, help='Overlap layer (not applicable).')
args = parser.parse_known_args()[0]

meshProperties = getMeshProperties()
meshProperties.quadRefine = args.quad
meshProperties.triRefine = args.tri
meshProperties.pipeOffset = args.offset
meshProperties.quadPipeRadius = args.qradius
meshProperties.triPipeRadius = args.tradius

print('Setup:')
print(f'Quad refinement: {meshProperties.quadRefine}')
print(f'Tri refinement: {meshProperties.triRefine}')
print(f'Pipe offset: {meshProperties.pipeOffset}')
print(f'Quad radius: {meshProperties.quadPipeRadius}')
print(f'Tri radius: {meshProperties.triPipeRadius}')
print(f'Write: {args.write}')
print("==========")

'''
Create SystemCoupling object.
'''
sc = sysc.SystemCoupling()

''' Add variables, regions, and the interface. '''
variableLinear1 = sysc.Variable("linear1", sysc.Scalar, False, sysc.Node)
variableLinear2 = sysc.Variable("linear2", sysc.Scalar, False, sysc.Node)

quadRegion = sysc.Region("quad", sysc.Surface)
quadRegion.addOutputVariable(variableLinear1)
quadRegion.addInputVariable(variableLinear2)

triRegion = sysc.Region("tri", sysc.Surface)
triRegion.addOutputVariable(variableLinear2)
triRegion.addInputVariable(variableLinear1)

couplingInterface = sysc.CouplingInterface("interface")
couplingInterface.addSideOneRegion(quadRegion)
couplingInterface.addSideTwoRegion(triRegion)

sc.addCouplingInterface(couplingInterface)

''' Register data access functions. '''
sc.registerSurfaceMeshAccess(getSurfaceMesh)
sc.registerInputScalarDataAccess(getInputScalarData)
sc.registerOutputScalarDataAccess(getOutputScalarData)

''' Initialize mesh and solution data on both sides of the interface. '''
startTime = time.time()
initializeMesh()
print(f"Time to generate mesh: {time.time() - startTime:.4e} [s]")

''' Update inputs - this will perform mapping and bring inputs up-to-date. '''
startTime = time.time()
sc.updateInputs()
print(f"Time to update inputs: {time.time() - startTime:.4e} [s]")

''' Write results, if desired. '''
if bool(args.write):
    startTime = time.time()
    sc.writeResults(sysc.ResultsInfo("pipe"))
    print(f"Time to write results: {time.time() - startTime:.4e} [s]")

printMaxError("quad")
printMaxError("tri")
