/*
 * Copyright ANSYS, Inc. Unauthorized use, distribution, or duplication is prohibited.
 */

#include "PipeMeshGenerator.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Define structure to hold partition bounds. */
typedef struct PartitionBounds {
  double axialStart;
  double axialEnd;
  double circStartAngle;
  double circEndAngle;
} PartitionBounds;

static PipeMeshProperties meshProperties = { 
  .quadRefine = 20, 
  .triRefine = 30, 
  .pipeOffset = 0.0, 
  .quadPipeRadius = 0.05, 
  .triPipeRadius = 0.05, 
  .overlapLayers = 0, 
  .pipeLength = 1.0,
  .myrank = 0, 
  .numranks = 1 };

static double pi = 3.14159265358979323846;

/* Define data structures for mesh arrays. */
static double* nodeCoordsQuad = NULL;
static double* nodeCoordsTri = NULL;
static int* elemNodeCountsQuad = NULL;
static int* elemNodeCountsTri = NULL;
static int* elemNodeIdsQuad = NULL;
static int* elemNodeIdsTri = NULL;

static size_t numNodesQuad = 0;
static size_t numNodesTri = 0;
static size_t numElemsQuad = 0;
static size_t numElemsTri = 0;

static int nodeIdOffset = 0;

/* Define data structures to hold solution data. */
static double* linear1Quad = NULL;
static double* linear2Quad = NULL;
static double* linear1Tri = NULL;
static double* linear2Tri = NULL;

PipeMeshProperties* getMeshProperties()
{
  return &meshProperties;
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
  /* Base case for the recursive function. */
  if (partitions == 1) {
    PartitionBounds bounds;
    bounds.axialStart = axialStart;
    bounds.axialEnd = axialEnd;
    bounds.circStartAngle = circStart;
    bounds.circEndAngle = circEnd;
    return bounds;
  }

  double axialLen = axialEnd - axialStart;
  double circLen = (circEnd - circStart) * pipeRadius;

  /* Determine if the current rank will be in upper (1) or lower (0) division. */
  int splitBucket = rank < partitions / 2 ? 0 : 1;
  int newRank = rank;
  int newPartitions = partitions / 2;
  if (splitBucket == 1) {
    newRank = rank - partitions / 2;
    newPartitions += partitions % 2;
  }

  /* Determine the split direction (0 is axial, 1 is circumferential). */
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

/* Convenience structure to hold an array of ints (for ranges). */
typedef struct Array {
  int* data;
  size_t size;
} Array;

static Array getArray(size_t size, int initValue)
{
  Array myArr;
  myArr.size = size;
  myArr.data = (int*) malloc (size * sizeof(int));
  for (size_t i = 0; i < size; ++i) {
    myArr.data[i] = initValue;
  }
  return myArr;
}

static void clearArray(Array* array)
{
  if (array->size > 0) {
    free(array->data);
    array->size = 0;
  }
}

static void addBeginning(Array* array, int value)
{
  int* newData = (int*) malloc ((array->size + 1) * sizeof(int));
  newData[0] = value;
  for (size_t i = 0; i < array->size; ++i) {
    newData[i + 1] = array->data[i];
  }
  size_t arrayNewSize = array->size + 1;
  clearArray(array);
  array->size = arrayNewSize;
  array->data = newData;
}

static void addEnd(Array* array, int value)
{
  int* newData = (int*) malloc ((array->size + 1) * sizeof(int));
  for (size_t i = 0; i < array->size; ++i) {
    newData[i] = array->data[i];
  }
  newData[array->size] = value;
  size_t arrayNewSize = array->size + 1;
  clearArray(array);
  array->size = arrayNewSize;
  array->data = newData;
}

/* Function to add overlap layers when running in parallel. */
static void addOverlapLayers(
  int circumferenceElems,
  int axialEdges,
  Array* elemAxialRange,
  Array* elemCircRange)
{
  /* Add overlap layers along axial direction. */
  if (elemAxialRange->size > 0) {
    for (int overlap = 0; overlap < meshProperties.overlapLayers; ++overlap) {
      int firstAxialRow = elemAxialRange->data[0];
      int lastAxialRow = elemAxialRange->data[elemAxialRange->size - 1];
      if (firstAxialRow > 0) {
        addBeginning(elemAxialRange, firstAxialRow - 1);
      }
      if (lastAxialRow < axialEdges - 1) {
        addEnd(elemAxialRange, lastAxialRow + 1);
      }
    }
  }

  /* Add overlap layers along circumferential direction. */
  if (elemCircRange->size > 0) {
    Array circElemFlags = getArray((size_t)circumferenceElems, 0);

    /* Mark elements already in the elemCircRange. */
    for (int e = 0; e < elemCircRange->size; ++e) {
      circElemFlags.data[elemCircRange->data[e]] = 1;
    }

    /* Mark elements in the overlap layer. */
    for (int overlap = 0; overlap < meshProperties.overlapLayers; ++overlap) {
      /* Add element at the beginning. */
      int firstCircRow = elemCircRange->data[0];
      int prevCircRow = (firstCircRow + circumferenceElems - overlap) % circumferenceElems;
      if (prevCircRow >= 0) {
        circElemFlags.data[prevCircRow] = 1;
      }
      /* Add element at the end. */
      int lastCircRow = elemCircRange->data[elemCircRange->size - 1];
      int nextCircRow = (lastCircRow + overlap) % circumferenceElems;
      circElemFlags.data[nextCircRow] = 1;
    }

    /* Clear and re-populate elemCircRange */
    clearArray(elemCircRange);
    for (int e = 0; e < circumferenceElems; ++e) {
      if (circElemFlags.data[e] == 1) {
        addEnd(elemCircRange, e);
      }
    }

    clearArray(&circElemFlags);
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

  int nodesPerFace = 4;
  double axialEdgeLen = meshProperties.pipeLength / (double) axialEdges;
  double circEdgeAngle = 2.0 * pi / (double) circumferenceEdges;

  /* Total number of nodes in this region (across all partitions). */
  int totalNumNodes = circumferenceEdges * (axialEdges + 1);

  /* Get the partition bounds. */
  PartitionBounds partBounds = getPartitionBounds(
    meshProperties.myrank, meshProperties.numranks, meshProperties.pipeOffset, 
    meshProperties.pipeOffset + meshProperties.pipeLength, 0.0, 2.0 * pi, 
    meshProperties.quadPipeRadius);

  Array elemAxialRange = getArray(0, 0);
  for (int e = 0; e < axialEdges; ++e) {
    double midpointAxial = meshProperties.pipeOffset + 0.5 * axialEdgeLen + e * axialEdgeLen;
    if (midpointAxial >= partBounds.axialStart &&
        midpointAxial < partBounds.axialEnd) {
      addEnd(&elemAxialRange, e);
    }
  }

  Array elemCircRange = getArray(0, 0);
  for (int c = 0; c < circumferenceEdges; ++c) {
    double midpointCirc = 0.5 * circEdgeAngle + c * circEdgeAngle;
    if (midpointCirc >= partBounds.circStartAngle &&
        midpointCirc < partBounds.circEndAngle) {
      addEnd(&elemCircRange, c);
    }
  }

  /* Add overlap layers. */
  addOverlapLayers(circumferenceEdges, axialEdges, &elemAxialRange, &elemCircRange);

  /* Calculate the number of elements and nodes in the current partition. */
  numElemsQuad = elemCircRange.size * elemAxialRange.size;
  numNodesQuad = 0;
  if (numElemsQuad > 0) {
    int nodesAlongCirc = (int)elemCircRange.size;
    if (elemCircRange.size < circumferenceEdges) {
      nodesAlongCirc = (int)elemCircRange.size + 1;
    }
    int nodesAlongAxial = (int)elemAxialRange.size + 1;
    numNodesQuad = (size_t) (nodesAlongCirc * nodesAlongAxial);
  }

  /* Allocate and initialize mesh and solution data arrays. */
  nodeCoordsQuad = (double*) malloc (numNodesQuad * 3 * sizeof(double));
  elemNodeCountsQuad = (int*) malloc (numElemsQuad * sizeof(int));
  for (size_t e = 0; e < numElemsQuad; e++) {
    elemNodeCountsQuad[e] = nodesPerFace;
  }
  elemNodeIdsQuad = (int*) malloc (numElemsQuad * (size_t)nodesPerFace * sizeof(int));

  linear1Quad = (double*) malloc (numNodesQuad * sizeof(double));
  linear2Quad = (double*) malloc (numNodesQuad * sizeof(double));

  /* If there are no elements in the current partition, return. */
  if (numElemsQuad > 0) {

    /* Calculate axial node range
      (same as element axial range plus one row of nodes at the end). */
    Array nodeAxialRange = getArray(elemAxialRange.size, 0);
    for (size_t i = 0; i < elemAxialRange.size; ++i) {
      nodeAxialRange.data[i] = elemAxialRange.data[i];
    }
    addEnd(&nodeAxialRange, elemAxialRange.data[elemAxialRange.size - 1] + 1);

    /* Calculate circumferential node range. */
    Array nodeCircFlags = getArray((size_t)circumferenceEdges, 0);
    for (size_t i = 0; i < elemCircRange.size; ++i) {
      int ec = elemCircRange.data[i];
      nodeCircFlags.data[ec] = 1;
      nodeCircFlags.data[(ec + 1) % circumferenceEdges] = 1;
    }

    Array nodeCircRange = getArray(0, 0);
    for (int n = 0; n < nodeCircFlags.size; n++) {
      if (nodeCircFlags.data[n] == 1) {
        addEnd(&nodeCircRange, n);
      }
    }

    /* Fill node coordinates. */
    int nodeIndex = 0;
    for (size_t nai = 0; nai < nodeAxialRange.size; ++nai) {
      int na = nodeAxialRange.data[nai];
      for (size_t nci = 0; nci < nodeCircRange.size; ++nci) {
        int nc = nodeCircRange.data[nci];
        double angle = 2.0 * pi *
                      (double)(nc) /
                      (double)(circumferenceEdges);
        double x = cos(angle) * meshProperties.quadPipeRadius;
        double y = sin(angle) * meshProperties.quadPipeRadius;
        double z = na * axialEdgeLen + meshProperties.pipeOffset;
        nodeCoordsQuad[nodeIndex * 3 + 0] = x;
        nodeCoordsQuad[nodeIndex * 3 + 1] = y;
        nodeCoordsQuad[nodeIndex * 3 + 2] = z;
        ++nodeIndex;
      }
    }

    /* Fill element to node connectivity. */
    int elemIndex = 0;
    for (size_t eai = 0; eai < elemAxialRange.size; ++eai) {
      int ea = elemAxialRange.data[eai];
      for (size_t eci = 0; eci < elemCircRange.size; ++eci) {
        int ec = elemCircRange.data[eci];
        /* Base node. */
        int baseNode = ea * circumferenceEdges + ec;
        /* Get next node along circumference. */
        int nextCirc = baseNode + 1;
        if (ec + 1 >= circumferenceEdges) {
          nextCirc = ea * circumferenceEdges;
        }
        /* Get next node along axial direction, relative to nextCirc node. */
        int nextAxialNextCirc = nextCirc + circumferenceEdges;
        /* Get next node along axial direction, relative to the baseNode. */
        int nextAxial = baseNode + circumferenceEdges;
        /* Fill element node ids array, taking node id offset into account. */
        elemNodeIdsQuad[elemIndex * nodesPerFace] = baseNode + nodeIdOffset;
        elemNodeIdsQuad[elemIndex * nodesPerFace + 1] = nextCirc + nodeIdOffset;
        elemNodeIdsQuad[elemIndex * nodesPerFace + 2] = nextAxialNextCirc + nodeIdOffset;
        elemNodeIdsQuad[elemIndex * nodesPerFace + 3] = nextAxial + nodeIdOffset;
        ++elemIndex;
      }
    }

    clearArray(&nodeAxialRange);
    clearArray(&nodeCircRange);
    clearArray(&nodeCircFlags);

    /* Fill solution data for linear1 variable with the linear profile. */
    for (int v = 0; v < numNodesQuad; ++v) {
      double x = nodeCoordsQuad[v * 3];
      double y = nodeCoordsQuad[v * 3 + 1];
      double z = nodeCoordsQuad[v * 3 + 2];
      linear1Quad[v] = 1.0 * x + 2.0 * y + 3.0 * z + 4.0;
    }

  }

  clearArray(&elemAxialRange);
  clearArray(&elemCircRange);

  nodeIdOffset = totalNumNodes;
}

/* 
 * Helper function to be used by initializeTri() function -
 * get range of nodes along circumference (local to this partition). 
 */
static Array getNodeCircRange(
  int axialRow, 
  const Array* elemCircRange, 
  int circumferenceEdges)
{
  Array circNodeFlags = getArray((size_t)circumferenceEdges, 0);

  for (size_t eci = 0; eci < elemCircRange->size; ++eci) {
    int e = elemCircRange->data[eci];
    if (!(axialRow % 2 == 1 && e % 2 == 1)) {
      circNodeFlags.data[(e / 2) % circumferenceEdges] = 1;
    }
    if (!(axialRow % 2 == 0 && e % 2 == 0)) {
      circNodeFlags.data[(e / 2 + 1) % circumferenceEdges] = 1;
    }
  }

  Array nodeCircRange = getArray(0, 0);
  for (int n = 0; n < circumferenceEdges; ++n) {
    if (circNodeFlags.data[n] == 1) {
      addEnd(&nodeCircRange, n);
    }
  }

  clearArray(&circNodeFlags);

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

  int nodesPerFace = 3;
  double axialEdgeLen = meshProperties.pipeLength / axialEdges;
  double circEdgeAngle = 2.0 * pi / circumferenceEdges;

  /* Total number of nodes in this region (across all partitions). */
  int totalNumNodes = circumferenceEdges * (axialEdges + 1);

  /* Get partition bounds. */
  PartitionBounds partBounds = getPartitionBounds(
    meshProperties.myrank, meshProperties.numranks, 0.0,
    meshProperties.pipeLength, 0.0, 2.0 * pi,
    meshProperties.triPipeRadius);

  Array elemAxialRange = getArray(0, 0);
  for (int e = 0; e < axialEdges; ++e) {
    double midpointAxial = 0.5 * axialEdgeLen + e * axialEdgeLen;
    if (midpointAxial >= partBounds.axialStart &&
        midpointAxial < partBounds.axialEnd) {
      addEnd(&elemAxialRange, e);
    }
  }

  Array elemCircRange = getArray(0, 0);
  for (int c = 0; c < circumferenceEdges; ++c) {
    /* There are two tri elements for every circumference edge. */
    double midpointCirc1 = (double)(c) * circEdgeAngle;
    double midpointCirc2 = ((double)(c) + 0.5) * circEdgeAngle;
    if (midpointCirc1 >= partBounds.circStartAngle &&
        midpointCirc1 < partBounds.circEndAngle) {
      addEnd(&elemCircRange, c * 2);
    }
    if (midpointCirc2 >= partBounds.circStartAngle &&
        midpointCirc2 < partBounds.circEndAngle) {
      addEnd(&elemCircRange, c * 2 + 1);
    }
  }

  /* Add overlap layers. */
  addOverlapLayers(circumferenceEdges * 2, axialEdges, &elemAxialRange, &elemCircRange);

  /* Calculate the number of elements and nodes in the current partition. */
  numElemsTri = elemCircRange.size * elemAxialRange.size;
  numNodesTri = 0;
  /* Calculate axial node range. */
  /* (same as element axial range plus one row of nodes at the end). */
  Array nodeAxialRange = getArray(0, 0);
  if (numElemsTri > 0) {
    for (size_t n = 0; n < elemAxialRange.size; ++n) {
      addEnd(&nodeAxialRange, elemAxialRange.data[n]);
    }
    addEnd(&nodeAxialRange, elemAxialRange.data[elemAxialRange.size - 1] + 1);

    /* Number of nodes along circumference can vary, depending on axial row.
      Loop over axial rows, use helper function to get number of nodes. */
    for (int nai = 0; nai < nodeAxialRange.size; ++nai) {
      int na = nodeAxialRange.data[nai];
      Array nodeCircRange = getNodeCircRange(na, &elemCircRange, circumferenceEdges);
      numNodesTri += nodeCircRange.size;
      clearArray(&nodeCircRange);
    }
  }

  /* Allocate and initialize mesh and solution data arrays. */
  nodeCoordsTri = (double*) malloc (numNodesTri * 3 * sizeof(double));
  elemNodeCountsTri = (int*) malloc (numElemsTri * sizeof(int));
  for (size_t e = 0; e < numElemsTri; ++e) {
    elemNodeCountsTri[e] = nodesPerFace;
  }
  elemNodeIdsTri = (int*) malloc (numElemsTri * (size_t)nodesPerFace * sizeof(int));

  linear1Tri = (double*) malloc (numNodesTri * sizeof(double));
  linear2Tri = (double*) malloc (numNodesTri * sizeof(double));

  /* If there are no elements in this partition, return. */
  if (numElemsTri == 0) {
    nodeIdOffset = totalNumNodes;
    return;
  }

  /* Fill node coordinates. */
  int nodeIndex = 0;
  for (size_t nai = 0; nai < nodeAxialRange.size; ++nai) {
    int na = nodeAxialRange.data[nai];
    Array nodeCircRange = getNodeCircRange(na, &elemCircRange, circumferenceEdges);
    for (size_t nci = 0; nci < nodeCircRange.size; ++nci) {
      int nc = nodeCircRange.data[nci];
      double offset = (na % 2) * (-0.5);
      double angle = 2.0 * pi * (nc + offset) / circumferenceEdges;
      double x = cos(angle) * meshProperties.triPipeRadius;
      double y = sin(angle) * meshProperties.triPipeRadius;
      double z = na * axialEdgeLen;
      nodeCoordsTri[nodeIndex * 3] = x;
      nodeCoordsTri[nodeIndex * 3 + 1] = y;
      nodeCoordsTri[nodeIndex * 3 + 2] = z;
      ++nodeIndex;
    }
    clearArray(&nodeCircRange);
  }

  /* Fill element-to-node connectivity. */
  int elemIndex = 0;
  for (size_t eai = 0; eai < elemAxialRange.size; ++eai) {
    int ea = elemAxialRange.data[eai];
    for (size_t eci = 0; eci < elemCircRange.size; ++eci) {
      int ec = elemCircRange.data[eci];
      /* Base node. */
      int baseNode = ea * circumferenceEdges + ec / 2;
      /* Next node along circumferential direction. */
      int nextCirc = baseNode + 1;
      if (baseNode + 1 >= (ea + 1) * circumferenceEdges) {
        nextCirc = ea * circumferenceEdges;
      }
      /* Next node along axial direction, relative to base node. */
      int nextAxial = baseNode + circumferenceEdges;
      /* Next node along circumferential direction, relative to nextCirc node. */
      int nextAxialNextCirc = nextCirc + circumferenceEdges;
      /* Get element node ids. Connectivity is different, depending on whether
        axial row and circumferential row are odd or even. */
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
      elemNodeIdsTri[elemIndex * nodesPerFace + 0] = a + nodeIdOffset;
      elemNodeIdsTri[elemIndex * nodesPerFace + 1] = b + nodeIdOffset;
      elemNodeIdsTri[elemIndex * nodesPerFace + 2] = c + nodeIdOffset;
      ++elemIndex;
    }
  }

  clearArray(&elemAxialRange);
  clearArray(&elemCircRange);
  clearArray(&nodeAxialRange);

  /* Fill solution data for linear2 variable with linear profile. */
  for (size_t v = 0; v < numNodesTri; ++v) {
    double x = nodeCoordsTri[v * 3];
    double y = nodeCoordsTri[v * 3 + 1];
    double z = nodeCoordsTri[v * 3 + 2];
    linear2Tri[v] = 1.0 * x + 2.0 * y + 3.0 * z + 4.0;
  }

  nodeIdOffset = totalNumNodes;
}

void initializeMesh()
{
  initializeQuad();
  initializeTri();
}

void clearMesh()
{
  if (numElemsQuad > 0) {
    free(nodeCoordsQuad);
    free(elemNodeCountsQuad);
    free(elemNodeIdsQuad);
    free(linear1Quad);
    free(linear2Quad);
  }

  if (numElemsTri > 0) {
    free(nodeCoordsTri);
    free(elemNodeCountsTri);
    free(elemNodeIdsTri);
    free(linear1Tri);
    free(linear2Tri);
  }
}

const double* getNodeCoords(const char* regionName)
{
  if (strcmp(regionName, "quad") == 0) {
    return nodeCoordsQuad;
  }
  else if (strcmp(regionName, "tri") == 0) {
    return nodeCoordsTri;
  }
  else {
    return NULL;
  }
}

size_t getNumNodes(const char* regionName)
{
  if (strcmp(regionName, "quad") == 0) {
    return numNodesQuad;
  }
  else if (strcmp(regionName, "tri") == 0) {
    return numNodesTri;
  }
  else {
    return 0;
  }
}

const int* getElemNodeCounts(const char* regionName)
{
  if (strcmp(regionName, "quad") == 0) {
    return elemNodeCountsQuad;
  }
  else if (strcmp(regionName, "tri") == 0) {
    return elemNodeCountsTri;
  }
  else {
    return NULL;
  }
}

size_t getNumElems(const char* regionName)
{
  if (strcmp(regionName, "quad") == 0) {
    return numElemsQuad;
  }
  else if (strcmp(regionName, "tri") == 0) {
    return numElemsTri;
  }
  else {
    return 0;
  }
}

const int* getElemNodeIds(const char* regionName)
{
  if (strcmp(regionName, "quad") == 0) {
    return elemNodeIdsQuad;
  }
  else if (strcmp(regionName, "tri") == 0) {
    return elemNodeIdsTri;
  }
  else {
    return NULL;
  }
}

size_t getNumElemNodeIds(const char* regionName)
{
  if (strcmp(regionName, "quad") == 0) {
    size_t nodesPerElem = 4;
    return numElemsQuad * nodesPerElem;
  }
  else if (strcmp(regionName, "tri") == 0) {
    size_t nodesPerElem = 3;
    return numElemsTri * nodesPerElem;
  }
  else {
    return 0;
  }
}

double* getSolutionData(const char* regionName, const char* variableName)
{
  if (strcmp(regionName, "quad") == 0) {
    if (strcmp(variableName, "linear1") == 0) {
      return linear1Quad;
    }
    else if (strcmp(variableName, "linear2") == 0) {
      return linear2Quad;
    }
  }
  else if (strcmp(regionName, "tri") == 0) {
    if (strcmp(variableName, "linear1") == 0) {
      return linear1Tri;
    }
    else if (strcmp(variableName, "linear2") == 0) {
      return linear2Tri;
    }
  }
  return NULL;
}
