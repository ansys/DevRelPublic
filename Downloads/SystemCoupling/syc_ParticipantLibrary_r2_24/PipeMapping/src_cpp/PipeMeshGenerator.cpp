/*
 * Copyright ANSYS, Inc. Unauthorized use, distribution, or duplication is prohibited.
 */

#include "PipeMeshGenerator.hpp"

#include <cmath>
#include <map>

// Define useful constants.
static double pi = 3.14159265358979323846;

// Set mesh properties.
static PipeMeshProperties meshProperties;

static int nodeIdOffset = 0;

// Define data structures for mesh arrays.
static std::map<std::string, std::vector<double>> nodeCoords;
static std::map<std::string, std::vector<int>> elemNodeCounts;
static std::map<std::string, std::vector<int>> elemNodeIds;

// Define data structures to hold solution data.
static std::map<std::string, std::map<std::string, std::vector<double>>> solutionData;

// Define structure to hold partition bounds.
struct PartitionBounds {
  double axialStart;
  double axialEnd;
  double circStartAngle;
  double circEndAngle;

  PartitionBounds(
    double axialStart,
    double axialEnd,
    double circStartAngle,
    double circEndAngle) :
      axialStart(axialStart),
      axialEnd(axialEnd),
      circStartAngle(circStartAngle),
      circEndAngle(circEndAngle)
  {
  }
};

// Get reference to mesh properties.
PipeMeshProperties& getMeshProperties()
{
  return meshProperties;
}

/*
 * Recursive function to calculate partition bounds for a patch of a
 * pipe.
 * 
 * Arguments are as follows:
 *
 * myrank - current rank 
 * 
 * partitions - number of partitions to break the current patch into 
 *
 * axialStart - start Z-coordinate of a patch to partition
 *
 * axialEnd - end Z-coordinate of a patch to partition 
 *
 * circStart - start angle (in radians) of a patch to partition 
 *
 * circEnd - end angle (in radians) of a patch to partition
 *
 * For example, patch can be defined as [0.5, 1.0, 0, 2.0 * PI]. This
 * defines a patch of the pipe from 0.5 to 1.0 in axial direction and full
 * circle in circumferential direction.
 */
static PartitionBounds getPartitionBounds(
  int rank, int partitions,
  double axialStart, double axialEnd,
  double circStart, double circEnd,
  double pipeRadius)
{
  // Base case for the recursive function.
  if (partitions == 1) {
    return PartitionBounds(axialStart, axialEnd, circStart, circEnd);
  }

  double axialLen = axialEnd - axialStart;
  double circLen = (circEnd - circStart) * pipeRadius;

  // Determine if the current rank will be in upper (1) or lower (0) division.
  int splitBucket = rank < partitions / 2 ? 0 : 1;
  int newRank = rank;
  int newPartitions = partitions / 2;
  if (splitBucket == 1) {
    newRank = rank - partitions / 2;
    newPartitions += partitions % 2;
  }

  // Determine the split direction (0 is axial, 1 is circumferential).
  int splitDir = axialLen >= circLen ? 0 : 1;
  double newAxialStart = axialStart;
  double newAxialEnd = axialEnd;
  double newCircStart = circStart;
  double newCircEnd = circEnd;
  if (splitDir == 0) {
    if (splitBucket == 0) {
      newAxialEnd = axialStart + 0.5 * axialLen;
    }
    else {
      newAxialStart = axialStart + 0.5 * axialLen;
    }
  }
  else {
    if (splitBucket == 0) {
      newCircEnd = circStart + 0.5 * (circEnd - circStart);
    }
    else {
      newCircStart = circStart + 0.5 * (circEnd - circStart);
    }
  }

  return getPartitionBounds(
    newRank, newPartitions, newAxialStart, newAxialEnd, newCircStart, newCircEnd, pipeRadius);
}

// Function to add overlap layers when running in parallel.
static void addOverlapLayers(
  int circumferenceElems,
  int axialEdges,
  std::vector<int>& elemAxialRange,
  std::vector<int>& elemCircRange)
{
  // Add overlap layers along axial direction.
  if (elemAxialRange.size() > 0) {
    for (int overlap = 0; overlap < meshProperties.overlapLayers; ++overlap) {
      int firstAxialRow = elemAxialRange.front();
      int lastAxialRow = elemAxialRange.back();
      if (firstAxialRow > 0) {
        elemAxialRange.insert(elemAxialRange.begin(), firstAxialRow - 1);
      }
      if (lastAxialRow < axialEdges - 1) {
        elemAxialRange.push_back(lastAxialRow + 1);
      }
    }
  }

  // Add overlap layers along circumferential direction.
  if (elemCircRange.size() > 0) {
    std::vector<int> circElemFlags(circumferenceElems, 0);

    // Mark elements already in the elemCircRange.
    for (std::size_t e = 0; e < elemCircRange.size(); ++e) {
      circElemFlags[elemCircRange[e]] = 1;
    }

    // Mark elements in the overlap layer.
    for (int overlap = 0; overlap < meshProperties.overlapLayers; ++overlap) {
      // Add element at the beginning.
      int firstCircRow = elemCircRange.front();
      int prevCircRow = (firstCircRow + circumferenceElems - overlap) % circumferenceElems;
      if (prevCircRow >= 0) {
        circElemFlags[prevCircRow] = 1;
      }
      //Add element at the end.
      int lastCircRow = elemCircRange.back();
      int nextCircRow = (lastCircRow + overlap) % circumferenceElems;
      circElemFlags[nextCircRow] = 1;
    }

    // Clear and re-populate elemCircRange
    elemCircRange.clear();
    for (int e = 0; e < circumferenceElems; ++e) {
      if (circElemFlags[e] == 1) {
        elemCircRange.push_back(e);
      }
    }
  }
}

/*
 * Initializes mesh on "quad" region (quad mesh).
 *
 * Here is the schematic of the mesh with 4 cirumference edges and 2 axial
 * edges.
 *
 *    circumrefential direction
 *              ^
 *              |
 *              |
 *              |
 *              x------> axial direction
 *
 * 0--------4--------8
 * |        |        |
 * |        |        |
 * |        |        |
 * 3--------7-------11
 * |        |        |
 * |        |        |
 * |        |        |
 * 2--------6-------10
 * |        |        |
 * |        |        |
 * |        |        |
 * 1--------5--------9
 * |        |        |
 * |        |        |
 * |        |        |
 * 0--------4--------8
 *
 */
static void initializeQuad()
{
  int circumferenceEdges = meshProperties.quadRefine;
  int axialEdges = meshProperties.quadRefine * 3;

  int myrank = meshProperties.myrank;
  int numranks = meshProperties.numranks;

  std::string regionName("quad");
  int nodesPerFace = 4;
  double axialEdgeLen = meshProperties.pipeLength / axialEdges;
  double circEdgeAngle = 2.0 * pi / circumferenceEdges;

  // Total number of nodes in this region (across all partitions).
  int totalNumNodes = circumferenceEdges * (axialEdges + 1);

  // Get the partition bounds.
  PartitionBounds partBounds = getPartitionBounds(
    myrank, numranks, meshProperties.pipeOffset,
    meshProperties.pipeOffset + meshProperties.pipeLength,
    0.0, 2.0 * pi, meshProperties.quadPipeRadius);
	
  std::vector<int> elemAxialRange;
  for (int e = 0; e < axialEdges; ++e) {
    double midpointAxial =
      meshProperties.pipeOffset + 0.5 * axialEdgeLen + e * axialEdgeLen;
    if (midpointAxial >= partBounds.axialStart &&
        midpointAxial < partBounds.axialEnd) {
      elemAxialRange.push_back(e);
    }
  }

  std::vector<int> elemCircRange;
  for (int c = 0; c < circumferenceEdges; ++c) {
    double midpointCirc = 0.5 * circEdgeAngle + c * circEdgeAngle;
    if (midpointCirc >= partBounds.circStartAngle &&
        midpointCirc < partBounds.circEndAngle) {
      elemCircRange.push_back(c);
    }
  }

  // Add overlap layers
  addOverlapLayers(circumferenceEdges, axialEdges, elemAxialRange, elemCircRange);

  // Calculate the number of elements and nodes in the current partition.
  int numElems = static_cast<int>(elemCircRange.size() * elemAxialRange.size());
  int numNodes = 0;
  if (numElems > 0) {
    int nodesAlongCirc = static_cast<int>(elemCircRange.size());
    if (static_cast<int>(elemCircRange.size()) < circumferenceEdges) {
      nodesAlongCirc = static_cast<int>(elemCircRange.size() + 1);
    }
    int nodesAlongAxial = static_cast<int>(elemAxialRange.size() + 1);
    numNodes = nodesAlongCirc * nodesAlongAxial;
  }

  // Allocate and initialize mesh and solution data arrays.
  nodeCoords.insert(
    std::make_pair(regionName, std::vector<double>(numNodes * 3)));
  elemNodeCounts.insert(
    std::make_pair(regionName, std::vector<int>(numElems, nodesPerFace)));
  elemNodeIds.insert(
    std::make_pair(regionName, std::vector<int>(numElems * nodesPerFace)));
  solutionData.insert(
    std::make_pair(regionName, std::map<std::string, std::vector<double>>()));
  solutionData[regionName].insert(
    std::make_pair("linear1", std::vector<double>(numNodes)));
  solutionData[regionName].insert(
    std::make_pair("linear2", std::vector<double>(numNodes)));

  // If there are no elements in the current partition, return.
  if (numElems == 0) {
    nodeIdOffset = totalNumNodes;
    return;
  }

  // Calculate axial node range
  // (same as element axial range plus one row of nodes at the end).
  std::vector<int> nodeAxialRange(elemAxialRange.begin(), elemAxialRange.end());
  nodeAxialRange.push_back(elemAxialRange.back() + 1);

  // Calculate circumferential node range.
  std::vector<int> nodeCircFlags(circumferenceEdges, 0);
  for (int e : elemCircRange) {
    nodeCircFlags[e] = 1;
    nodeCircFlags[(e + 1) % circumferenceEdges] = 1;
  }
  std::vector<int> nodeCircRange;
  for (std::size_t n = 0; n < nodeCircFlags.size(); ++n) {
    if (nodeCircFlags[n] == 1) {
      nodeCircRange.push_back(static_cast<int>(n));
    }
  }

  // Fill node coordinates.
  std::vector<double>& currNodeCoords = nodeCoords.at(regionName);
  int nodeIndex = 0;
  for (const int& na : nodeAxialRange) {
    for (const int& nc : nodeCircRange) {
      double angle = 2.0 * pi *
                     static_cast<double>(nc) /
                     static_cast<double>(circumferenceEdges);
      double x = std::cos(angle) * meshProperties.quadPipeRadius;
      double y = std::sin(angle) * meshProperties.quadPipeRadius;
      double z = na * axialEdgeLen + meshProperties.pipeOffset;
      currNodeCoords[nodeIndex * 3 + 0] = x;
      currNodeCoords[nodeIndex * 3 + 1] = y;
      currNodeCoords[nodeIndex * 3 + 2] = z;
      ++nodeIndex;
    }
  }

  // Fill element to node connectivity.
  std::vector<int>& currElemNodeIds = elemNodeIds.at(regionName);
  int elemIndex = 0;
  for (const int& ea : elemAxialRange) {
    for (const int& ec : elemCircRange) {
      // Base node.
      int baseNode = ea * circumferenceEdges + ec;
      // Get next node along circumference.
      int nextCirc = baseNode + 1;
      if (ec + 1 >= circumferenceEdges) {
        nextCirc = ea * circumferenceEdges;
      }
      // Get next node along axial direction, relative to nextCirc node.
      int nextAxialNextCirc = nextCirc + circumferenceEdges;
      // Get next node along axial direction, relative to the baseNode.
      int nextAxial = baseNode + circumferenceEdges;
      // Fill element node ids array, taking node id offset into account.
      currElemNodeIds[elemIndex * nodesPerFace] = baseNode + nodeIdOffset;
      currElemNodeIds[elemIndex * nodesPerFace + 1] = nextCirc + nodeIdOffset;
      currElemNodeIds[elemIndex * nodesPerFace + 2] = nextAxialNextCirc + nodeIdOffset;
      currElemNodeIds[elemIndex * nodesPerFace + 3] = nextAxial + nodeIdOffset;
      ++elemIndex;
    }
  }

  // Fill solution data for linear1 variable with the linear profile.
  std::vector<double>& currSolution = solutionData.at(regionName).at("linear1");
  for (int v = 0; v < numNodes; ++v) {
    double x = currNodeCoords[v * 3];
    double y = currNodeCoords[v * 3 + 1];
    double z = currNodeCoords[v * 3 + 2];
    currSolution[v] = 1.0 * x + 2.0 * y + 3.0 * z + 4.0;
  }

  nodeIdOffset = totalNumNodes;
}

/* 
 * Helper function to be used by initializeTri() function -
 * get range of nodes along circumference (local to this partition). 
 */
static std::vector<int> getNodeCircRange(
  int axialRow, const std::vector<int>& elemCircRange, int circumferenceEdges)
{
  std::vector<int> circNodeFlags(circumferenceEdges, 0);
  for (int e : elemCircRange) {
    if (!(axialRow % 2 == 1 && e % 2 == 1)) {
      circNodeFlags[(e / 2) % circumferenceEdges] = 1;
    }
    if (!(axialRow % 2 == 0 && e % 2 == 0)) {
      circNodeFlags[(e / 2 + 1) % circumferenceEdges] = 1;
    }
  }

  std::vector<int> nodeCircRange;
  for (int n = 0; n < circumferenceEdges; ++n) {
    if (circNodeFlags[n] == 1) {
      nodeCircRange.push_back(n);
    }
  }

  return nodeCircRange;
}

/*
 * Initializes mesh on "tri" region (tri mesh).
 *
 * Here is the schematic of the mesh with 4 cirumference edges and 2 axial
 * edges.
 *
 *    circumrefential direction
 *              ^
 *              |
 *              |
 *              |
 *              x------> axial direction
 *
 * 0                 0
 * |  \     |     /  |
 * |     \  |  /     |
 * |        0        |
 * |     /  |  \     |
 * |  /     |     \  |
 * 3        |        3
 * |  \     |     /  |
 * |     \  |  /     |
 * |        3        |
 * |     /  |  \     |
 * |  /     |     \  |
 * 2        |        2
 * |  \     |     /  |
 * |     \  |  /     |
 * |        2        |
 * |     /  |  \     |
 * |  /     |     \  |
 * 1        |        1
 * |  \     |     /  |
 * |     \  |  /     |
 * |        1        |
 * |     /  |  \     |
 * |  /     |     \  |
 * 0        |        0
 *    \     |     /  
 *       \  |  /     
 *          0        
 *
 */
static void initializeTri()
{
  int circumferenceEdges = meshProperties.triRefine;
  int axialEdges = meshProperties.triRefine * 3;

  int myrank = meshProperties.myrank;
  int numranks = meshProperties.numranks;

  std::string regionName("tri");
  int nodesPerFace = 3;
  double axialEdgeLen = meshProperties.pipeLength / axialEdges;
  double circEdgeAngle = 2.0 * pi / circumferenceEdges;

  // Total number of nodes in this region (across all partitions).
  int totalNumNodes = circumferenceEdges * (axialEdges + 1);

  // Get partition bounds.
  PartitionBounds partBounds = getPartitionBounds(
    myrank, numranks, 0.0, meshProperties.pipeLength,
    0.0, 2.0 * pi, meshProperties.triPipeRadius);

  std::vector<int> elemAxialRange;
  for (int e = 0; e < axialEdges; ++e) {
    double midpointAxial = 0.5 * axialEdgeLen + e * axialEdgeLen;
    if (midpointAxial >= partBounds.axialStart &&
        midpointAxial < partBounds.axialEnd) {
      elemAxialRange.push_back(e);
    }
  }

  std::vector<int> elemCircRange;
  for (int c = 0; c < circumferenceEdges; ++c) {
    // There are two tri elements for every circumference edge.
    double midpointCirc1 = static_cast<double>(c) * circEdgeAngle;
    double midpointCirc2 = (static_cast<double>(c) + 0.5) * circEdgeAngle;
    if (midpointCirc1 >= partBounds.circStartAngle &&
        midpointCirc1 < partBounds.circEndAngle) {
      elemCircRange.push_back(c * 2);
    }
    if (midpointCirc2 >= partBounds.circStartAngle &&
        midpointCirc2 < partBounds.circEndAngle) {
      elemCircRange.push_back(c * 2 + 1);
    }
  }

  // Add overlap layers
  addOverlapLayers(circumferenceEdges * 2, axialEdges, elemAxialRange, elemCircRange);

  // Calculate the number of elements and nodes in the current partition.
  int numElems = static_cast<int>(elemCircRange.size() * elemAxialRange.size());
  int numNodes = 0;
  // Calculate axial node range
  // (same as element axial range plus one row of nodes at the end).
  std::vector<int> nodeAxialRange;
  if (numElems > 0) {
    nodeAxialRange.assign(elemAxialRange.begin(), elemAxialRange.end());
    nodeAxialRange.push_back(elemAxialRange.back() + 1);
    // Number of nodes along circumference can vary, depending on axial row.
    // Loop over axial rows, use helper function to get number of nodes.
    for (const int& an : nodeAxialRange) {
      std::vector<int> nodeCircRange = getNodeCircRange(
        an, elemCircRange, circumferenceEdges);
      numNodes += static_cast<int>(nodeCircRange.size());
    }
  }

  // Allocate and initialize mesh and solution data arrays.
  nodeCoords.insert(
    std::make_pair(regionName, std::vector<double>(numNodes * 3)));
  elemNodeCounts.insert(
    std::make_pair(regionName, std::vector<int>(numElems, nodesPerFace)));
  elemNodeIds.insert(
    std::make_pair(regionName, std::vector<int>(numElems * nodesPerFace)));
  solutionData.insert(
    std::make_pair(regionName, std::map<std::string, std::vector<double>>()));
  solutionData[regionName].insert(
    std::make_pair("linear1", std::vector<double>(numNodes)));
  solutionData[regionName].insert(
    std::make_pair("linear2", std::vector<double>(numNodes)));

  // If there are no elements in this partition, return.
  if (numElems == 0) {
    nodeIdOffset = totalNumNodes;
    return;
  }

  // Fill node coordinates.
  std::vector<double>& currNodeCoords = nodeCoords.at(regionName);
  int nodeIndex = 0;
  for (const int& na : nodeAxialRange) {
    std::vector<int> nodeCircRange = getNodeCircRange(
      na, elemCircRange, circumferenceEdges);
    for (const int& nc : nodeCircRange) {
      double offset = (na % 2) * (-0.5);
      double angle = 2.0 * pi * (nc + offset) / circumferenceEdges;
      double x = std::cos(angle) * meshProperties.triPipeRadius;
      double y = std::sin(angle) * meshProperties.triPipeRadius;
      double z = na * axialEdgeLen;
      currNodeCoords[nodeIndex * 3] = x;
      currNodeCoords[nodeIndex * 3 + 1] = y;
      currNodeCoords[nodeIndex * 3 + 2] = z;
      ++nodeIndex;
    }
  }

  // Fill element-to-node connectivity.
  std::vector<int>& currElemNodeIds = elemNodeIds.at(regionName);
  int elemIndex = 0;
  for (const int& ea : elemAxialRange) {
    for (const int& ec : elemCircRange) {
      // Base node.
      int baseNode = ea * circumferenceEdges + ec / 2;
      // Next node along circumferential direction.
      int nextCirc = baseNode + 1;
      if (baseNode + 1 >= (ea + 1) * circumferenceEdges) {
        nextCirc = ea * circumferenceEdges;
      }
      // Next node along axial direction, relative to base node.
      int nextAxial = baseNode + circumferenceEdges;
      // Next node along circumferential direction, relative to nextCirc node.
      int nextAxialNextCirc = nextCirc + circumferenceEdges;
      // Get element node ids. Connectivity is different, depending on whether
      // axial row and circumferential row are odd or even.
      int a, b, c;
      if (ea % 2 == 0) {
        if (ec % 2 == 0) {
          a = baseNode;
          b = nextAxialNextCirc;
          c = nextAxial;
        }
        else {
          a = baseNode;
          b = nextCirc;
          c = nextAxialNextCirc;
        }
      }
      else {
        if (ec % 2 == 0) {
          a = baseNode;
          b = nextCirc;
          c = nextAxial;
        }
        else {
          a = nextCirc;
          b = nextAxialNextCirc;
          c = nextAxial;
        }
      }
      currElemNodeIds[elemIndex * nodesPerFace + 0] = a + nodeIdOffset;
      currElemNodeIds[elemIndex * nodesPerFace + 1] = b + nodeIdOffset;
      currElemNodeIds[elemIndex * nodesPerFace + 2] = c + nodeIdOffset;
      ++elemIndex;
    }
  }

  // Fill solution data for linear2 variable with linear profile.
  std::vector<double>& currSolution = solutionData.at(regionName).at("linear2");
  for (int v = 0; v < numNodes; ++v) {
    double x = currNodeCoords[v * 3];
    double y = currNodeCoords[v * 3 + 1];
    double z = currNodeCoords[v * 3 + 2];
    currSolution[v] = 1.0 * x + 2.0 * y + 3.0 * z + 4.0;
  }

  nodeIdOffset = totalNumNodes;
}

void initializeMesh()
{
  initializeQuad();
  initializeTri();
}

const std::vector<double>& getNodeCoords(const std::string& regionName)
{
  return nodeCoords.at(regionName);
}

const std::vector<int>& getElemNodeCounts(const std::string& regionName)
{
  return elemNodeCounts.at(regionName);
}

const std::vector<int>& getElemNodeIds(const std::string& regionName)
{
  return elemNodeIds.at(regionName);
}

std::vector<double>& getSolutionData(
  const std::string& regionName, const std::string& variableName)
{
  return solutionData.at(regionName).at(variableName);
}
