/*
 * Copyright ANSYS, Inc. Unauthorized use, distribution, or duplication is prohibited.
 */

#include "PipeMeshGenerator.h"

#include "SystemCouplingParticipant/syscSystemCoupling.h"

#include "mpi.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(_MSC_VER)
/* Disable Microsoft compiler warnings about standard C functions. */
#pragma warning (disable : 4996)
#endif

/* Define useful constants. */
#define STRING_MAX_SIZE 4096

/* Function that returns surface mesh, given a region. */
SyscSurfaceMesh getSurfaceMesh(const char* regionName)
{
  /* node coords */
  SyscOutputVectorData nodeCoordsData = syscGetOutputVectorDataCompactDouble(
    getNodeCoords(regionName), getNumNodes(regionName));

  /* elem node counts */
  SyscOutputIntegerData elemNodeCountsData = syscGetOutputIntegerDataInt32(
    getElemNodeCounts(regionName), getNumElems(regionName));

  /* elem node ids */
  SyscOutputIntegerData elemNodeIdsData = syscGetOutputIntegerDataInt32(
    getElemNodeIds(regionName), getNumElemNodeIds(regionName));

  /* nodes */
  SyscNodeData nodes = syscGetNodeDataC(nodeCoordsData);

  /* faces */
  SyscFaceData faces = syscGetFaceDataCN(
    syscGetElementNodeCountData(elemNodeCountsData), 
    syscGetElementNodeConnectivityData(elemNodeIdsData));

  return syscGetSurfaceMeshNF(nodes, faces);
}

/*    
 * Function that returns scalar input variable values,
 * given a region and a variable.
 */
SyscInputScalarData getInputScalarData(
  const char* regionName,
  const char* variableName)
{
  return syscGetInputScalarDataDouble(
    getSolutionData(regionName, variableName), 
    getNumNodes(regionName));
}

/*
 * Function that returns scalar output variable values,
 * given a region and a variable.
 */
SyscOutputScalarData getOutputScalarData(
  const char* regionName,
  const char* variableName)
{
  return syscGetOutputScalarDataDouble(
    getSolutionData(regionName, variableName), 
    getNumNodes(regionName));
}

/* Function to calculate and print max error on a given region. */
void printMaxError(const char* regionName)
{
  int isQuadRegion = strcmp(regionName, "quad") == 0;
  const double* testValues;
  const double* expectedValues;
  size_t numVals = getNumNodes(regionName);
  if (isQuadRegion) {
    testValues = getSolutionData("quad", "linear2");
    expectedValues = getSolutionData("quad", "linear1");
  }
  else {
    testValues = getSolutionData("tri", "linear1");
    expectedValues = getSolutionData("tri", "linear2");
  }

  double maxError = 0.0, totalMaxError;
  for (size_t n = 0; n < numVals; ++n) {
    double testValue = testValues[n];
    double expectedValue = expectedValues[n];
    double currDiff = fabs(expectedValue - testValue);
    maxError = currDiff > maxError ? currDiff : maxError;
  }

  MPI_Reduce(&maxError, &totalMaxError, 1, 
    MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (getMeshProperties()->myrank == 0) {
    printf("Maximum error on %s: %.4e\n", regionName, totalMaxError);
  }
}

/* Function to write out mesh and solution data arrays to text files. */
void writeArrays(const char* regionName)
{
  const double* currNodeCoords = getNodeCoords(regionName);
  const int* currElemNodeCounts = getElemNodeCounts(regionName);
  const int* currElemNodeIds = getElemNodeIds(regionName);
  size_t numNodes = getNumNodes(regionName);
  size_t numElems = getNumElems(regionName);

  char rankStr[STRING_MAX_SIZE];
  sprintf(rankStr, "%d", getMeshProperties()->myrank);

  char fileName[STRING_MAX_SIZE];
  fileName[0] = '\0';
  strcat(fileName, regionName);
  strcat(fileName, ".");
  strcat(fileName, rankStr);
  strcat(fileName, ".dat");

  FILE* fileHandle = fopen(fileName, "w");

  /* Write number of nodes, elements, and solution data arrays. */
  size_t numScalarNodeArrays = 2;  /* "linear1" and "linear2" */
  fprintf(fileHandle, "%zu, %zu, %zu\n", numNodes, numElems, 
    numScalarNodeArrays);

  /* Write node coordinates. */
  for (size_t n = 0; n < numNodes; ++n) {
    fprintf(fileHandle, "%.8e, %.8e, %.8e\n",
      currNodeCoords[n * 3 + 0],
      currNodeCoords[n * 3 + 1],
      currNodeCoords[n * 3 + 2]);
  }

  /* Write element-to-node connectivity. */
  size_t elemNodeIndex = 0;
  for (size_t e = 0; e < numElems; ++e) {
    /* Write the first node id for the current element. */
    fprintf(fileHandle, "%d", currElemNodeIds[elemNodeIndex]);
    elemNodeIndex++;
    /* Write subsequent node ids for the current element. */
    for (int en = 1; en < currElemNodeCounts[e]; ++en) {
      fprintf(fileHandle, ", %d", currElemNodeIds[elemNodeIndex]);
      elemNodeIndex++;
    }
    fprintf(fileHandle, "\n");
  }

  /* Write solution data arrays. */
  fprintf(fileHandle, "linear1\n");
  const double* currData = getSolutionData(regionName, "linear1");
  for (size_t n = 0; n < numNodes; n++) {
    fprintf(fileHandle, "%.8e\n", currData[n]);
  }

  fprintf(fileHandle, "linear2\n");
  currData = getSolutionData(regionName, "linear2");
  for (size_t n = 0; n < numNodes; n++) {
    fprintf(fileHandle, "%.8e\n", currData[n]);
  }
}

void processRetCode(SyscError ret)
{
  if (ret.retcode != 0) {
    printf("ERROR: %s\n", ret.message);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
}

/*
 * This is the main function for the pipe mapping program.
 *
 * It generates two different meshes, each for a region representing a pipe
 * geometry. It then performs non-conformal mesh mapping between these two
 * regions.
 *
 * The mesh refinement level, pipe offset, pipe radius, etc. can be controlled
 * by the command line arguments. The resulting arrays can be written for
 * further post-processing.
 *
 */
int main(int argc, char* argv[])
{
  /* Parse command line arguments. */
  int write = 1;

  PipeMeshProperties* meshProperties = getMeshProperties();

  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "--quad") == 0 && i + 1 < argc) {
      meshProperties->quadRefine = atoi(argv[i + 1]);
    }
    else if (strcmp(argv[i], "--tri") == 0 && i + 1 < argc) {
      meshProperties->triRefine = atoi(argv[i + 1]);
    }
    else if (strcmp(argv[i], "--offset") == 0 && i + 1 < argc) {
      meshProperties->pipeOffset = (double) atof(argv[i + 1]);
    }
    else if (strcmp(argv[i], "--qradius") == 0 && i + 1 < argc) {
      meshProperties->quadPipeRadius = (double) atof(argv[i + 1]);
    }
    else if (strcmp(argv[i], "--tradius") == 0 && i + 1 < argc) {
      meshProperties->triPipeRadius = (double) atof(argv[i + 1]);
    }
    else if (strcmp(argv[i], "--write") == 0 && i + 1 < argc) {
      write = atoi(argv[i + 1]);
    }
    else if (strcmp(argv[i], "--overlap") == 0 && i + 1 < argc) {
      meshProperties->overlapLayers = atoi(argv[i + 1]);
    }
  }

  /* Initialize MPI. */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &meshProperties->myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &meshProperties->numranks);

  if (meshProperties->myrank == 0) {
    printf("Setup:\n");
    printf("Quad refinement: %d\n", meshProperties->quadRefine);
    printf("Tri refinement: %d\n", meshProperties->triRefine);
    printf("Pipe offset: %f\n", meshProperties->pipeOffset);
    printf("Quad radius: %f\n", meshProperties->quadPipeRadius);
    printf("Tri radius: %f\n", meshProperties->triPipeRadius);
    printf("Write: %d\n", write);
    printf("Overlap layers: %d\n", meshProperties->overlapLayers);
    printf("Processes: %d\n", meshProperties->numranks);
    printf("==========\n");
  }

  SyscError ret;

  /* Add variables, regions, and the interface. */
  SyscVariable variableLinear1 = syscGetVariableTE(
    "linear1", SyscScalar, 0, SyscNode);
  SyscVariable variableLinear2 = syscGetVariableTE(
    "linear2", SyscScalar, 0, SyscNode);
 
  SyscRegion quadRegion = syscGetRegionT("quad", SyscSurface);

  ret = syscAddOutputVariable(quadRegion, variableLinear1);
  processRetCode(ret);

  ret = syscAddInputVariable(quadRegion, variableLinear2);
  processRetCode(ret);

  SyscRegion triRegion = syscGetRegionT("tri", SyscSurface);

  ret = syscAddOutputVariable(triRegion, variableLinear2);
  processRetCode(ret);

  ret = syscAddInputVariable(triRegion, variableLinear1);
  processRetCode(ret);

  SyscCouplingInterface couplingInterface = syscGetCouplingInterface(
    "interface");

  ret = syscAddSideOneRegion(couplingInterface, quadRegion);
  processRetCode(ret);

  ret = syscAddSideTwoRegion(couplingInterface, triRegion);
  processRetCode(ret);

  ret = syscAddCouplingInterface(couplingInterface);
  processRetCode(ret);

  /* Register data access functions. */
  ret = syscRegisterSurfMeshAccess(&getSurfaceMesh);
  processRetCode(ret);

  ret = syscRegisterInputScalarDataAccess(&getInputScalarData);
  processRetCode(ret);

  ret = syscRegisterOutputScalarDataAccess(&getOutputScalarData);
  processRetCode(ret);

  double startTime, localTime, totalTime;

  /* Initialize mesh and solution data on both sides of the interface. */
  startTime = MPI_Wtime();
  initializeMesh();
  localTime = MPI_Wtime() - startTime;
  MPI_Reduce(&localTime, &totalTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (meshProperties->myrank == 0) {
    printf("Time to generate mesh: %.4e [s]\n", totalTime);
  }

  /* Update inputs - this will perform mapping and bring inputs up-to-date. */
  startTime = MPI_Wtime();
  ret = syscUpdateInputs();
  localTime = MPI_Wtime() - startTime;
  processRetCode(ret);
  MPI_Reduce(&localTime, &totalTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (meshProperties->myrank == 0) {
    printf("Time to update inputs: %.4e [s]\n", totalTime);
  }

  /* Write results, if desired. */
  if (write != 0) {
    startTime = MPI_Wtime();
    writeArrays("quad");
    writeArrays("tri");
    ret = syscWriteResults(syscGetResultsInfo("pipe"));
    localTime = MPI_Wtime() - startTime;
    MPI_Reduce(&localTime, &totalTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (meshProperties->myrank == 0) {
      printf("Time to write results: %.4e [s]\n", totalTime);
    }
  }

  printMaxError("quad");
  printMaxError("tri");

  clearMesh();

  MPI_Finalize();
  return EXIT_SUCCESS;
}
