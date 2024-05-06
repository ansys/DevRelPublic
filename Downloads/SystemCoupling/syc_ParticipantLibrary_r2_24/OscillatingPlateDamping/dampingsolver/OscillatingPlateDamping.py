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



"""
See Participant Library's Oscillating Plate Damping tutorial
for a detailed description of the problem.
"""

""" Define surface mesh. """

"""
Define node coordinates in a single array (compact form).
There are 24 nodes.
"""
nodeCoords = np.array(
    [
        1.00600004e001,
        1.00000000e000,
        4.00000006e-001,
        1.00600004e001,
        1.00000000e000,
        0.00000000e000,
        1.00000000e001,
        1.00000000e000,
        4.00000006e-001,
        1.00000000e001,
        1.00000000e000,
        0.00000000e000,
        1.00600004e001,
        0.00000000e000,
        4.00000006e-001,
        1.00600004e001,
        2.00000003e-001,
        4.00000006e-001,
        1.00600004e001,
        4.00000006e-001,
        4.00000006e-001,
        1.00600004e001,
        6.00000024e-001,
        4.00000006e-001,
        1.00600004e001,
        8.00000012e-001,
        4.00000006e-001,
        1.00600004e001,
        8.00000012e-001,
        0.00000000e000,
        1.00600004e001,
        6.00000024e-001,
        0.00000000e000,
        1.00600004e001,
        4.00000006e-001,
        0.00000000e000,
        1.00600004e001,
        2.00000003e-001,
        0.00000000e000,
        1.00600004e001,
        0.00000000e000,
        0.00000000e000,
        1.00000000e001,
        8.00000012e-001,
        4.00000006e-001,
        1.00000000e001,
        6.00000024e-001,
        4.00000006e-001,
        1.00000000e001,
        4.00000006e-001,
        4.00000006e-001,
        1.00000000e001,
        2.00000003e-001,
        4.00000006e-001,
        1.00000000e001,
        0.00000000e000,
        4.00000006e-001,
        1.00000000e001,
        0.00000000e000,
        0.00000000e000,
        1.00000000e001,
        8.00000012e-001,
        0.00000000e000,
        1.00000000e001,
        6.00000024e-001,
        0.00000000e000,
        1.00000000e001,
        4.00000006e-001,
        0.00000000e000,
        1.00000000e001,
        2.00000003e-001,
        0.00000000e000,
    ]
)

nodeCount = int(len(nodeCoords) / 3)

""" Define element node counts (11 quad elements). """
elementNodeCounts = np.array([4] * 11)

""" Define element node ids for each element (11 quad elements). """
elementNodeIds = np.array(
    [
        13,
        4,
        5,
        12,
        12,
        5,
        6,
        11,
        11,
        6,
        7,
        10,
        10,
        7,
        8,
        9,
        9,
        8,
        0,
        1,
        23,
        17,
        18,
        19,
        22,
        16,
        17,
        23,
        21,
        15,
        16,
        22,
        20,
        14,
        15,
        21,
        3,
        2,
        14,
        20,
        2,
        3,
        1,
        0,
    ]
)

""" Define array of nodal displacements (24 nodes, 3 components per node). """
displacement = np.array([0.0] * 3 * nodeCount)

""" Define array of nodal forces (24 nodes, 3 components per node). """
force = np.array([0.0] * 3 * nodeCount)

""" Keeps track of current time step. """
currStep = 0

""" Define function that returns surface mesh, given a region. """


def getSurfaceMesh(regionName):
    if regionName == "plate":
        nodes = sysc.NodeData(sysc.OutputVectorData(nodeCoords))
        faces = sysc.FaceData(
            sysc.ElementNodeCountData(sysc.OutputIntegerData(elementNodeCounts)),
            sysc.ElementNodeConnectivityData(sysc.OutputIntegerData(elementNodeIds)),
        )
        return sysc.SurfaceMesh(nodes, faces)

    raise RuntimeError("getSurfaceMesh")


"""
Define a function that returns vector input variable values,
given a region and a variable.
"""


def getInputVectorData(regionName, variableName):
    if regionName == "plate" and variableName == "Mesh Displacement":
        return sysc.InputVectorData(displacement)
    raise RuntimeError("getInputVectorData")


"""
Define a function that returns vector output variable values,
given a region and a variable.
"""


def getOutputVectorData(regionName, variableName):
    if regionName == "plate" and variableName == "Damping Force":
        return sysc.OutputVectorData(force)
    raise RuntimeError("getOutputVectorData")


def createRestartPoint():
    """
    Define a function that creates a restart point and returns
    a string that identifies that restart point.
    In this implementation, the file with up-to-date
    nodal coordinates and force values
    will be written and the name of that file is returned.
    Only those values are required to perform the restart
    for this problem.
    """
    restartFileName = "restart." + str(currStep) + ".txt"
    with open(restartFileName, "w") as f:
        for n in range(nodeCount):
            line = f"{nodeCoords[n * 3 + 0]:.8e}, "
            line += f"{nodeCoords[n * 3 + 1]:.8e}, "
            line += f"{nodeCoords[n * 3 + 2]:.8e}, "
            line += f"{force[n * 3 + 0]:.8e}, "
            line += f"{force[n * 3 + 1]:.8e}, "
            line += f"{force[n * 3 + 2]:.8e}\n"
            f.write(line)
    return restartFileName


def resumeState(restartPoint):
    global currStep
    """
    # Error checking is skipped here for simplicity.
    # Assume the file exists and is in the format
    # consistent with how it gets written
    # in the createRestartPoint() function.
    """

    """ Set current time step. """
    currStep = int(restartPoint.split(".")[1])

    """ Set node coordinates and forces from the restart file. """
    with open(restartPoint) as f:
        for n in range(nodeCount):
            line = f.readline()
            tokens = line.split(",")
            for component in range(3):
                nodeCoords[n * 3 + component] = float(tokens[component])
                force[n * 3 + component] = float(tokens[3 + component])


""" Parse and print input arguments. """
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--schost", type=str, default="")
parser.add_argument("--scport", type=int, default=0)
parser.add_argument("--scname", type=str, default="")
parser.add_argument("--scsetup", default=False, action="store_true")
parser.add_argument("--dampcoeff", type=float, default=0.0)
parser.add_argument("--screstart", type=str)
args = parser.parse_args()

dampCoeff = args.dampcoeff

print("Host: " + args.schost)
print("Port: " + str(args.scport))
print("Name: " + args.scname)
if args.screstart is not None:
    print("Restart: " + args.screstart)
print("Damping coefficient: " + str(dampCoeff))

buildInfo = "Oscillating Plate Damping Solver, language: Python"

""" Connect to System Coupling. """
sc = sysc.SystemCoupling(args.schost, args.scport, args.scname, buildInfo)

if args.scsetup:

    """Setup mode."""

    """ Create the surface region "plate". """
    region = sysc.Region("plate", sysc.Surface)

    """ Create the force variable. """
    forceVar = sysc.Variable("Damping Force", sysc.Vector, True, sysc.Node)

    """ Create the displacement variable. """
    displacementVar = sysc.Variable("Mesh Displacement", sysc.Vector, False, sysc.Node)

    """ Add force as output and displacement as input to plate. """
    region.addOutputVariable(forceVar)
    region.addInputVariable(displacementVar)

    """ Add region to the setup. """
    sc.addRegion(region)

    """
    Create setup info structure. Set analysis type to transient.
    Set restarts supported to true.
    """
    setupInfo = sysc.SetupInfo(sysc.Transient, True)

    """ Complete the setup. """
    sc.completeSetup(setupInfo)

else:
    """Register heavyweight data access functions."""
    sc.registerSurfaceMeshAccess(getSurfaceMesh)
    sc.registerInputVectorDataAccess(getInputVectorData)
    sc.registerOutputVectorDataAccess(getOutputVectorData)

    """ Register restart point creation funciton. """
    sc.registerRestartPointCreation(createRestartPoint)

    """ If doing a restart, resume the state from the given restart point. """
    if args.screstart:
        resumeState(args.screstart)

    """ Initialize coupled analysis. """
    sc.initializeAnalysis()

    """ Enter coupled analysis loop. """

    """ Coupling time step loop. """
    while sc.doTimeStep():
        """Increment time step."""
        currStep += 1
        """ Update nodal coordinates. """
        for n in range(len(nodeCoords)):
            nodeCoords[n] += displacement[n]
        """ Coupling iteration loop. """
        while sc.doIteration():
            """Inputs update."""
            sc.updateInputs()
            """
            Get time step size, which is used to caclulate
            nodal velocities.
            """
            ts = sc.getCurrentTimeStep()
            dt = ts.timeStepSize
            """ Calculate the damping force. """
            for n in range(len(force)):
                force[n] = -dampCoeff * displacement[n] / dt
            """
            Outputs update.
            Note: convergence status is Complete, since this is
            not an iterative solver.
            """
            sc.updateOutputs(sysc.Complete)

""" Coupled analyis shutdown. """
sc.disconnect()
