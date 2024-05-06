#
# Copyright ANSYS, Inc. Unauthorized use, distribution, or duplication is prohibited.
#

# Participant Library #
# Heat Transfer in Square Channel Air Flow Tutorial - Parallel #

import os
import platform

# Set the environment variable that points
# to the mock solver executable script.
def getMockSolverExecutable():
    solverScript = 'ChannelFlowMockSolverMPI.bat'
    if platform.system() != 'Windows':
        solverScript = 'ChannelFlowMockSolverMPI.sh'
    return os.path.join('mocksolver', solverScript)

# Add Fluent participant.
fluentName = AddParticipant(InputFile = os.path.join('fluent', 'fluent.scp'))

# Add mock solver participant.
mockSolverName = AddParticipant(InputFile = os.path.join('mocksolver', 'mocksolver.scp'))

# Set mock solver executable.
ec = DatamodelRoot().CouplingParticipant[mockSolverName].ExecutionControl
ec.Executable = getMockSolverExecutable()

# Define the coupling interface.
interfaceName = AddInterface(
    SideOneParticipant = fluentName,
    SideOneRegions = ['sc'],
    SideTwoParticipant = mockSolverName,
    SideTwoRegions = ['FSI'])

# Define heat flow data transfer.
AddDataTransfer(
    Interface = interfaceName,
    TargetSide = 'One',
    SideOneVariable = 'heatflow',
    SideTwoVariable = 'Heat Rate')

# Define temperature data transfer.
AddDataTransfer(
    Interface = interfaceName,
    TargetSide = 'Two',
    SideOneVariable = 'temperature',
    SideTwoVariable = 'Temperature')

# Set maximum number of coupling iterations.
DatamodelRoot().SolutionControl.MaximumIterations = 15

# Solution
Solve()
