'''
Copyright ANSYS, Inc. Unauthorized use, distribution, or duplication is prohibited.
'''

import math
import numpy as np

class PipeMeshProperties:
    ''' Define structure to hold mesh properties. '''
    def __init__(self):
        self.quadRefine = 20
        self.triRefine = 30
        self.pipeOffset = 0.0
        self.quadPipeRadius = 0.05
        self.triPipeRadius = 0.05
        self.pipeLength = 1.0

meshProperties = PipeMeshProperties()

''' Define data structures for mesh arrays. '''
nodeCoords = {}
elemNodeCounts = {}
elemNodeIds = {}

''' Define data structures to hold solution data. '''
solutionData = {}

nodeIdOffset = 0

def getMeshProperties():
    return meshProperties

def initializeQuad():
    '''
    Initializes mesh on "quad" region (quad mesh).

    Here is the schematic of the mesh with 4 cirumference edges and 2 axial
    edges.
 
       circumrefential direction
                 ^
                 |
                 |
                 |
                 x------> axial direction
  
    0--------4--------8
    |        |        |
    |        |        |
    |        |        |
    3--------7-------11
    |        |        |
    |        |        |
    |        |        |
    2--------6-------10
    |        |        |
    |        |        |
    |        |        |
    1--------5--------9
    |        |        |
    |        |        |
    |        |        |
    0--------4--------8
    '''

    global nodeCoords, elemNodeCounts, elemNodeIds
    global nodeIdOffset, solutionData, nodeIdOffset

    circumferenceEdges = meshProperties.quadRefine
    axialEdges = meshProperties.quadRefine * 3

    regionName = "quad"
    nodesPerFace = 4
    axialEdgeLen = meshProperties.pipeLength / axialEdges

    ''' Total number of nodes in this region (across all partitions). '''
    totalNumNodes = circumferenceEdges * (axialEdges + 1)

    elemAxialRange = list(range(axialEdges))
    elemCircRange = list(range(circumferenceEdges))

    ''' Calculate the number of elements and nodes in the current partition. '''
    numElems = len(elemCircRange) * len(elemAxialRange)
    numNodes = 0
    if numElems > 0:
        nodesAlongCirc = len(elemCircRange)
        nodesAlongAxial = (len(elemAxialRange) + 1)
        numNodes = nodesAlongCirc * nodesAlongAxial

    ''' Allocate and initialize mesh and solution data arrays. '''
    nodeCoords[regionName] = np.array([0.0] * numNodes * 3, dtype = np.float64)
    elemNodeCounts[regionName] = np.array([nodesPerFace] * numElems, dtype = np.int32)
    elemNodeIds[regionName] = np.array([0] * numElems * nodesPerFace, dtype = np.int32)
    solutionData[regionName] = {}
    solutionData[regionName]["linear1"] = np.array([0.0] * numNodes, dtype = np.float64)
    solutionData[regionName]["linear2"] = np.array([0.0] * numNodes, dtype = np.float64)

    ''' If there are no elements in the current partition, return. '''
    if numElems == 0:
        nodeIdOffset = totalNumNodes
        return

    '''
    Calculate axial node range
    (same as element axial range plus one row of nodes at the end).
    '''
    nodeAxialRange = list(elemAxialRange) + [elemAxialRange[-1] + 1]

    ''' Calculate circumferential node range. '''
    nodeCircRange = list(elemCircRange)

    ''' Fill node coordinates. '''
    currNodeCoords = nodeCoords[regionName]
    nodeIndex = 0
    for na in nodeAxialRange:
        for nc in nodeCircRange:
            angle = 2.0 * math.pi * nc / circumferenceEdges
            x = math.cos(angle) * meshProperties.quadPipeRadius
            y = math.sin(angle) * meshProperties.quadPipeRadius
            z = na * axialEdgeLen + meshProperties.pipeOffset
            currNodeCoords[nodeIndex * 3 + 0] = x
            currNodeCoords[nodeIndex * 3 + 1] = y
            currNodeCoords[nodeIndex * 3 + 2] = z
            nodeIndex += 1

    ''' Fill element to node connectivity. '''
    currElemNodeIds = elemNodeIds[regionName]
    elemIndex = 0
    for ea in elemAxialRange:
        for ec in elemCircRange:
            ''' Base node. '''
            baseNode = ea * circumferenceEdges + ec
            ''' Get next node along circumference. '''
            nextCirc = baseNode + 1 if ec + 1 < circumferenceEdges else ea * circumferenceEdges
            ''' Get next node along axial direction, relative to nextCirc node. '''
            nextAxialNextCirc = nextCirc + circumferenceEdges
            ''' Get next node along axial direction, relative to the baseNode. '''
            nextAxial = baseNode + circumferenceEdges
            ''' Fill element node ids array, taking node id offset into account. '''
            currElemNodeIds[elemIndex * nodesPerFace] = baseNode + nodeIdOffset
            currElemNodeIds[elemIndex * nodesPerFace + 1] = nextCirc + nodeIdOffset
            currElemNodeIds[elemIndex * nodesPerFace + 2] = nextAxialNextCirc + nodeIdOffset
            currElemNodeIds[elemIndex * nodesPerFace + 3] = nextAxial + nodeIdOffset
            elemIndex += 1

    ''' Fill solution data for linear1 variable with the linear profile. '''
    currSolution = solutionData[regionName]["linear1"]
    for v in range(numNodes):
        x = currNodeCoords[v * 3]
        y = currNodeCoords[v * 3 + 1]
        z = currNodeCoords[v * 3 + 2]
        currSolution[v] = 1.0 * x + 2.0 * y + 3.0 * z + 4.0

    nodeIdOffset = totalNumNodes

def initializeTri():
    '''
    Initializes mesh on "tri" region (tri mesh).

    Here is the schematic of the mesh with 4 cirumference edges and 2 axial
    edges.
 
       circumrefential direction
                 ^
                 |
                 |
                 |
                 x------> axial direction
 
    0                 0
    |  \     |     /  |
    |     \  |  /     |
    |        0        |
    |     /  |  \     |
    |  /     |     \  |
    3        |        3
    |  \     |     /  |
    |     \  |  /     |
    |        3        |
    |     /  |  \     |
    |  /     |     \  |
    2        |        2
    |  \     |     /  |
    |     \  |  /     |
    |        2        |
    |     /  |  \     |
    |  /     |     \  |
    1        |        1
    |  \     |     /  |
    |     \  |  /     |
    |        1        |
    |     /  |  \     |
    |  /     |     \  |
    0        |        0
       \     |     /  
          \  |  /     
             0        
    '''

    global nodeCoords, elemNodeCounts, elemNodeIds
    global nodeIdOffset, solutionData, nodeIdOffset

    circumferenceEdges = meshProperties.triRefine
    axialEdges = meshProperties.triRefine * 3

    regionName = "tri"
    nodesPerFace = 3
    axialEdgeLen = meshProperties.pipeLength / axialEdges

    ''' Total number of nodes in this region (across all partitions). '''
    totalNumNodes = circumferenceEdges * (axialEdges + 1)

    elemAxialRange = list(range(axialEdges))

    ''' There are two tri elements for every circumference edge. '''
    elemCircRange = list(range(circumferenceEdges * 2))

    nodeCircRange = list(range(circumferenceEdges))

    ''' Calculate the number of elements and nodes in the current partition. '''
    numElems = len(elemCircRange) * len(elemAxialRange)

    numNodes = 0
    '''
    Calculate axial node range
    (same as element axial range plus one row of nodes at the end).
    '''
    nodeAxialRange = []
    if numElems > 0:
        nodeAxialRange = list(elemAxialRange) + [elemAxialRange[-1] + 1]
        numNodes = len(nodeAxialRange) * len(nodeCircRange)

    ''' Allocate and initialize mesh and solution data arrays. '''
    nodeCoords[regionName] = np.array([0] * numNodes * 3, dtype = np.float64)
    elemNodeCounts[regionName] = np.array([nodesPerFace] * numElems, dtype = np.int32)
    elemNodeIds[regionName] = np.array([0] * numElems * nodesPerFace, dtype = np.int32)
    solutionData[regionName] = {}
    solutionData[regionName]["linear1"] = np.array([0.0] * numNodes, dtype = np.float64)
    solutionData[regionName]["linear2"] = np.array([0.0] * numNodes, dtype = np.float64)

    ''' If there are no elements in this partition, return. '''
    if numElems == 0:
        nodeIdOffset = totalNumNodes

    ''' Fill node coordinates. '''
    currNodeCoords = nodeCoords[regionName]
    nodeIndex = 0
    for na in nodeAxialRange:
        for nc in nodeCircRange:
            offset = (na % 2) * (-0.5)
            angle = 2 * math.pi * (nc + offset) / circumferenceEdges
            x = math.cos(angle) * meshProperties.triPipeRadius
            y = math.sin(angle) * meshProperties.triPipeRadius
            z = na * axialEdgeLen
            currNodeCoords[nodeIndex * 3] = x
            currNodeCoords[nodeIndex * 3 + 1] = y
            currNodeCoords[nodeIndex * 3 + 2] = z
            nodeIndex += 1

    ''' Fill element-to-node connectivity. '''
    currElemNodeIds = elemNodeIds[regionName]
    elemIndex = 0
    for ea in elemAxialRange:
        for ec in elemCircRange:
            ''' Base node. '''
            baseNode = ea * circumferenceEdges + int(ec / 2)
            ''' Next node along circumferential direction. '''
            nextCirc = baseNode + 1 if baseNode + 1 < (ea + 1) * circumferenceEdges else ea * circumferenceEdges
            ''' Next node along axial direction, relative to base node. '''
            nextAxial = baseNode + circumferenceEdges
            ''' Next node along circumferential direction, relative to nextCirc node. '''
            nextAxialNextCirc = nextCirc + circumferenceEdges
            '''
            Get element node ids. Connectivity is different, depending on whether
            axial row and circumferential row are odd or even.
            '''
            if ea % 2 == 0:
                if ec % 2 == 0:
                    a = baseNode
                    b = nextAxialNextCirc
                    c = nextAxial
                else:
                    a = baseNode
                    b = nextCirc
                    c = nextAxialNextCirc
            else:
                if ec % 2 == 0:
                    a = baseNode
                    b = nextCirc
                    c = nextAxial
                else:
                    a = nextCirc
                    b = nextAxialNextCirc
                    c = nextAxial
            currElemNodeIds[elemIndex * nodesPerFace + 0] = a + nodeIdOffset
            currElemNodeIds[elemIndex * nodesPerFace + 1] = b + nodeIdOffset
            currElemNodeIds[elemIndex * nodesPerFace + 2] = c + nodeIdOffset
            elemIndex += 1

    ''' Fill solution data for linear2 variable with linear profile. '''
    currSolution = solutionData[regionName]["linear2"]
    for v in range(numNodes):
        x = currNodeCoords[v * 3]
        y = currNodeCoords[v * 3 + 1]
        z = currNodeCoords[v * 3 + 2]
        currSolution[v] = 1.0 * x + 2.0 * y + 3.0 * z + 4.0

    nodeIdOffset = totalNumNodes

def initializeMesh():
    initializeQuad()
    initializeTri()

def getNodeCoords(regionName):
    return nodeCoords[regionName]

def getElemNodeCounts(regionName):
    return elemNodeCounts[regionName]

def getElemNodeIds(regionName):
    return elemNodeIds[regionName]

def getSolutionData(regionName, variableName):
    return solutionData[regionName][variableName]
