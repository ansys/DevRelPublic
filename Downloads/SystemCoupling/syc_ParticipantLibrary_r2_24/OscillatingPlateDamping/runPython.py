#
# Copyright ANSYS, Inc. Unauthorized use, distribution, or duplication is prohibited.
#

# Participant Library #
# Transient Oscillating Plate Damping Tutorial #

import os
import platform

def getDampingSolverExecutable():
    solverScript = 'OscillatingPlateDampingPython.bat'
    if platform.system() != 'Windows':
        solverScript = 'OscillatingPlateDampingPython.sh'
    return os.path.join('dampingsolver', solverScript)

# Setup

# Add the Mechanical solver as a coupling particpant.
mapdlName = AddParticipant(InputFile = os.path.join('mapdl', 'mapdl.scp'))

# Add the Damping Solver as a coupling participant
# via the direct setup approach.
dampingSolverName = AddParticipant(
    Executable = getDampingSolverExecutable(),
    WorkingDirectory = 'dampingsolver',
    AdditionalArguments = '--dampcoeff 1.0')

# Define the coupled analysis controls.
DatamodelRoot().SolutionControl.EndTime = '10.0 [s]'
DatamodelRoot().SolutionControl.TimeStepSize = '0.1 [s]'

# Add the coupling interface.
iName = AddInterface(
    SideOneParticipant = mapdlName, 
    SideOneRegions = ['FSIN_1'],
    SideTwoParticipant = dampingSolverName,
    SideTwoRegions = ['plate'])

# Add displacement transfer.
AddDataTransfer(
    Interface = iName,
    SideOneVariable = 'INCD',
    SideTwoVariable = 'Mesh Displacement',
    TargetSide = 'Two')

# Add force transfer.
AddDataTransfer(
    Interface = iName,
    SideOneVariable = 'FORC',
    SideTwoVariable = 'Damping Force',
    TargetSide = 'One')

# Define restart point controls - generate a
# restart point at the end of every 10th time step.
DatamodelRoot().OutputControl.Option = 'StepInterval'
DatamodelRoot().OutputControl.OutputFrequency = 10

# Solve the coupled analysis.
Solve()
