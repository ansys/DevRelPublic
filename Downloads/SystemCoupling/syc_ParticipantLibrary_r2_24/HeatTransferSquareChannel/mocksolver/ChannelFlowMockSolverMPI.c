/*
 * Copyright ANSYS, Inc. Unauthorized use, distribution, or duplication is prohibited.
 */

#include "SystemCouplingParticipant/syscSystemCoupling.h"

#include "mpi.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(_MSC_VER)
/* Disable Microsoft compiler warnings about standard C functions. */
#pragma warning (disable : 4996)
#endif

#define STRING_MAX_SIZE 4096

/*
* Mesh representation:
*   10 nodes
*   4 quad elements
*   Rectangle, aligned with XY-plane,
*   stretching from the origin to point (1, 0.1, 0)
*
* (0,0.1,0)                    (1,0.1,0)
*          6----7----8----9----10
*          | A  | B  | C  | D  |
*          |    |    |    |    |
*          1----2----3----4----5
* (0,0,0)                       (1,0,0)
*
* In parallel, the above mesh will be partitioned as follows:
*
* numranks : 1,  elements on rank 0 : A, B, C, D
* numranks : 2,  elements on rank 0 : A, B
*                elements on rank 1 : C, D
* numranks : 3,  elements on rank 0 : A
*                elements on rank 1 : B
*                elements on rank 2 : C, D
* numranks : 4+, elements on rank 0 : A
*                elements on rank 1 : B
*                elements on rank 2 : C
*                elements on rank 3 : D
*                elements on rank 4+: no elements
*
*/

/* 
*  Declare arrays used for this problem.
*  The contents are populated in the initializeArray()
*  function, depending on parallel settings.
*/ 

size_t numNodes = 0;
size_t numElems = 0;

/* Array of node coordinates in split array format. */
double* nodeCoordsX = NULL;
double* nodeCoordsY = NULL;
double* nodeCoordsZ = NULL;

/* Array of element node counts. */
uint16_t* elementNodeCounts = NULL;

/* Array of node ids for each element. */
uint64_t* elementNodeIds = NULL;

/* Array of temperature values (defined at nodes). */
double* nodeTemperature = NULL;

/*
 * Array of heat rate values (defined at elements)
 * values are in Watts - total of 20 [W] generated
 */
double* elementHeatRate = NULL;

/*
*  Function that partitions the mesh based
*  on the number of parallel processes, and
*  initializes the arrays.
*/
void initializeArrays()
{
  /* Find out total number of ranks and local rank. */
  int myrank = 0;
  int numranks = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &numranks);
  
  /* Partition the mesh based on number of parallel processes. */
  if (numranks == 1) {

    /* All four faces are present on rank 0. */

    numNodes = 10;
    numElems = 4;

    nodeCoordsX = (double*) malloc (numNodes * sizeof(double));
    nodeCoordsY = (double*) malloc (numNodes * sizeof(double));
    nodeCoordsZ = (double*) malloc (numNodes * sizeof(double));
    elementNodeCounts = (uint16_t*) malloc (numElems * sizeof(uint16_t));
    elementNodeIds = (uint64_t*) malloc (numElems * 4 * sizeof(uint64_t));

    nodeCoordsX[0] = 0.00; nodeCoordsY[0] = 0.0;
    nodeCoordsX[1] = 0.25; nodeCoordsY[1] = 0.0;
    nodeCoordsX[2] = 0.50; nodeCoordsY[2] = 0.0;
    nodeCoordsX[3] = 0.75; nodeCoordsY[3] = 0.0;
    nodeCoordsX[4] = 1.00; nodeCoordsY[4] = 0.0;
    nodeCoordsX[5] = 0.00; nodeCoordsY[5] = 0.1;
    nodeCoordsX[6] = 0.25; nodeCoordsY[6] = 0.1;
    nodeCoordsX[7] = 0.50; nodeCoordsY[7] = 0.1;
    nodeCoordsX[8] = 0.75; nodeCoordsY[8] = 0.1;
    nodeCoordsX[9] = 1.00; nodeCoordsY[9] = 0.1;

    elementNodeIds[0] = 1;
    elementNodeIds[1] = 2;
    elementNodeIds[2] = 7;
    elementNodeIds[3] = 6;

    elementNodeIds[4] = 2;
    elementNodeIds[5] = 3;
    elementNodeIds[6] = 8;
    elementNodeIds[7] = 7;

    elementNodeIds[8] = 3;
    elementNodeIds[9] = 4;
    elementNodeIds[10] = 9;
    elementNodeIds[11] = 8;

    elementNodeIds[12] = 4;
    elementNodeIds[13] = 5;
    elementNodeIds[14] = 10;
    elementNodeIds[15] = 9;
  }
  else if (numranks == 2) {

    if (myrank == 0) {

      /* Faces A and B are on rank 0. */

      numNodes = 6;
      numElems = 2;

      nodeCoordsX = (double*) malloc (numNodes * sizeof(double));
      nodeCoordsY = (double*) malloc (numNodes * sizeof(double));
      nodeCoordsZ = (double*) malloc (numNodes * sizeof(double));
      elementNodeCounts = (uint16_t*) malloc (numElems * sizeof(uint16_t));
      elementNodeIds = (uint64_t*) malloc (numElems * 4 * sizeof(uint64_t));

      nodeCoordsX[0] = 0.00; nodeCoordsY[0] = 0.0;
      nodeCoordsX[1] = 0.25; nodeCoordsY[1] = 0.0;
      nodeCoordsX[2] = 0.50; nodeCoordsY[2] = 0.0;
      nodeCoordsX[3] = 0.00; nodeCoordsY[3] = 0.1;
      nodeCoordsX[4] = 0.25; nodeCoordsY[4] = 0.1;
      nodeCoordsX[5] = 0.50; nodeCoordsY[5] = 0.1;

      elementNodeIds[0] = 1;
      elementNodeIds[1] = 2;
      elementNodeIds[2] = 7;
      elementNodeIds[3] = 6;

      elementNodeIds[4] = 2;
      elementNodeIds[5] = 3;
      elementNodeIds[6] = 8;
      elementNodeIds[7] = 7;

    }
    else {

      /* Faces C and D are on rank 1. */

      numNodes = 6;
      numElems = 2;

      nodeCoordsX = (double*) malloc (numNodes * sizeof(double));
      nodeCoordsY = (double*) malloc (numNodes * sizeof(double));
      nodeCoordsZ = (double*) malloc (numNodes * sizeof(double));
      elementNodeCounts = (uint16_t*) malloc (numElems * sizeof(uint16_t));
      elementNodeIds = (uint64_t*) malloc (numElems * 4 * sizeof(uint64_t));

      nodeCoordsX[0] = 0.50; nodeCoordsY[0] = 0.0;
      nodeCoordsX[1] = 0.75; nodeCoordsY[1] = 0.0;
      nodeCoordsX[2] = 1.00; nodeCoordsY[2] = 0.0;
      nodeCoordsX[3] = 0.50; nodeCoordsY[3] = 0.1;
      nodeCoordsX[4] = 0.75; nodeCoordsY[4] = 0.1;
      nodeCoordsX[5] = 1.00; nodeCoordsY[5] = 0.1;

      elementNodeIds[0] = 3;
      elementNodeIds[1] = 4;
      elementNodeIds[2] = 9;
      elementNodeIds[3] = 8;

      elementNodeIds[4] = 4;
      elementNodeIds[5] = 5;
      elementNodeIds[6] = 10;
      elementNodeIds[7] = 9;

    }

  }
  else if (numranks == 3) {

    if (myrank == 0) {

      /* Face A is on rank 0. */

      numNodes = 4;
      numElems = 1;

      nodeCoordsX = (double*) malloc (numNodes * sizeof(double));
      nodeCoordsY = (double*) malloc (numNodes * sizeof(double));
      nodeCoordsZ = (double*) malloc (numNodes * sizeof(double));
      elementNodeCounts = (uint16_t*) malloc (numElems * sizeof(uint16_t));
      elementNodeIds = (uint64_t*) malloc (numElems * 4 * sizeof(uint64_t));

      nodeCoordsX[0] = 0.00; nodeCoordsY[0] = 0.0;
      nodeCoordsX[1] = 0.25; nodeCoordsY[1] = 0.0;
      nodeCoordsX[2] = 0.00; nodeCoordsY[2] = 0.1;
      nodeCoordsX[3] = 0.25; nodeCoordsY[3] = 0.1;

      elementNodeIds[0] = 1;
      elementNodeIds[1] = 2;
      elementNodeIds[2] = 7;
      elementNodeIds[3] = 6;

    }
    else if (myrank == 1) {

      /* Face B is on rank 1. */

      numNodes = 4;
      numElems = 1;

      nodeCoordsX = (double*) malloc (numNodes * sizeof(double));
      nodeCoordsY = (double*) malloc (numNodes * sizeof(double));
      nodeCoordsZ = (double*) malloc (numNodes * sizeof(double));
      elementNodeCounts = (uint16_t*) malloc (numElems * sizeof(uint16_t));
      elementNodeIds = (uint64_t*) malloc (numElems * 4 * sizeof(uint64_t));

      nodeCoordsX[0] = 0.25; nodeCoordsY[0] = 0.0;
      nodeCoordsX[1] = 0.50; nodeCoordsY[1] = 0.0;
      nodeCoordsX[2] = 0.25; nodeCoordsY[2] = 0.1;
      nodeCoordsX[3] = 0.50; nodeCoordsY[3] = 0.1;

      elementNodeIds[0] = 2;
      elementNodeIds[1] = 3;
      elementNodeIds[2] = 8;
      elementNodeIds[3] = 7;

    }
    else {

      /* Faces C and D are on rank 2. */

      numNodes = 6;
      numElems = 2;

      nodeCoordsX = (double*) malloc (numNodes * sizeof(double));
      nodeCoordsY = (double*) malloc (numNodes * sizeof(double));
      nodeCoordsZ = (double*) malloc (numNodes * sizeof(double));
      elementNodeCounts = (uint16_t*) malloc (numElems * sizeof(uint16_t));
      elementNodeIds = (uint64_t*) malloc (numElems * 4 * sizeof(uint64_t));

      nodeCoordsX[0] = 0.50; nodeCoordsY[0] = 0.0;
      nodeCoordsX[1] = 0.75; nodeCoordsY[1] = 0.0;
      nodeCoordsX[2] = 1.00; nodeCoordsY[2] = 0.0;
      nodeCoordsX[3] = 0.50; nodeCoordsY[3] = 0.1;
      nodeCoordsX[4] = 0.75; nodeCoordsY[4] = 0.1;
      nodeCoordsX[5] = 1.00; nodeCoordsY[5] = 0.1;

      elementNodeIds[0] = 3;
      elementNodeIds[1] = 4;
      elementNodeIds[2] = 9;
      elementNodeIds[3] = 8;

      elementNodeIds[4] = 4;
      elementNodeIds[5] = 5;
      elementNodeIds[6] = 10;
      elementNodeIds[7] = 9;

    }

  }
  else {

    /* Four or more processes.
     * Ranks 0, 1, 2, 3 get faces A, B, C, D, respectively.
     * Any higher ranks (if exist) have to mesh.
     */

    if (myrank == 0) {

      /* Face A is on rank 0. */

      numNodes = 4;
      numElems = 1;

      nodeCoordsX = (double*) malloc (numNodes * sizeof(double));
      nodeCoordsY = (double*) malloc (numNodes * sizeof(double));
      nodeCoordsZ = (double*) malloc (numNodes * sizeof(double));
      elementNodeCounts = (uint16_t*) malloc (numElems * sizeof(uint16_t));
      elementNodeIds = (uint64_t*) malloc (numElems * 4 * sizeof(uint64_t));

      nodeCoordsX[0] = 0.00; nodeCoordsY[0] = 0.0;
      nodeCoordsX[1] = 0.25; nodeCoordsY[1] = 0.0;
      nodeCoordsX[2] = 0.00; nodeCoordsY[2] = 0.1;
      nodeCoordsX[3] = 0.25; nodeCoordsY[3] = 0.1;

      elementNodeIds[0] = 1;
      elementNodeIds[1] = 2;
      elementNodeIds[2] = 7;
      elementNodeIds[3] = 6;

    }
    else if (myrank == 1) {

      /* Face B is on rank 1. */

      numNodes = 4;
      numElems = 1;

      nodeCoordsX = (double*) malloc (numNodes * sizeof(double));
      nodeCoordsY = (double*) malloc (numNodes * sizeof(double));
      nodeCoordsZ = (double*) malloc (numNodes * sizeof(double));
      elementNodeCounts = (uint16_t*) malloc (numElems * sizeof(uint16_t));
      elementNodeIds = (uint64_t*) malloc (numElems * 4 * sizeof(uint64_t));

      nodeCoordsX[0] = 0.25; nodeCoordsY[0] = 0.0;
      nodeCoordsX[1] = 0.50; nodeCoordsY[1] = 0.0;
      nodeCoordsX[2] = 0.25; nodeCoordsY[2] = 0.1;
      nodeCoordsX[3] = 0.50; nodeCoordsY[3] = 0.1;

      elementNodeIds[0] = 2;
      elementNodeIds[1] = 3;
      elementNodeIds[2] = 8;
      elementNodeIds[3] = 7;

    }
    else if (myrank == 2) {

      /* Face C is on rank 2. */

      numNodes = 4;
      numElems = 1;

      nodeCoordsX = (double*) malloc (numNodes * sizeof(double));
      nodeCoordsY = (double*) malloc (numNodes * sizeof(double));
      nodeCoordsZ = (double*) malloc (numNodes * sizeof(double));
      elementNodeCounts = (uint16_t*) malloc (numElems * sizeof(uint16_t));
      elementNodeIds = (uint64_t*) malloc (numElems * 4 * sizeof(uint64_t));

      nodeCoordsX[0] = 0.50; nodeCoordsY[0] = 0.0;
      nodeCoordsX[1] = 0.75; nodeCoordsY[1] = 0.0;
      nodeCoordsX[2] = 0.50; nodeCoordsY[2] = 0.1;
      nodeCoordsX[3] = 0.75; nodeCoordsY[3] = 0.1;

      elementNodeIds[0] = 3;
      elementNodeIds[1] = 4;
      elementNodeIds[2] = 9;
      elementNodeIds[3] = 8;

    }
    else if (myrank == 3) {

      /* Face D is on rank 3. */

      numNodes = 4;
      numElems = 1;

      nodeCoordsX = (double*) malloc (numNodes * sizeof(double));
      nodeCoordsY = (double*) malloc (numNodes * sizeof(double));
      nodeCoordsZ = (double*) malloc (numNodes * sizeof(double));
      elementNodeCounts = (uint16_t*) malloc (numElems * sizeof(uint16_t));
      elementNodeIds = (uint64_t*) malloc (numElems * 4 * sizeof(uint64_t));

      nodeCoordsX[0] = 0.75; nodeCoordsY[0] = 0.0;
      nodeCoordsX[1] = 1.00; nodeCoordsY[1] = 0.0;
      nodeCoordsX[2] = 0.75; nodeCoordsY[2] = 0.1;
      nodeCoordsX[3] = 1.00; nodeCoordsY[3] = 0.1;

      elementNodeIds[0] = 4;
      elementNodeIds[1] = 5;
      elementNodeIds[2] = 10;
      elementNodeIds[3] = 9;

    }
    else {
      /* All arrays remain empty for ranks 4 and higher. */
    }

  }

  /* Resize element heat rate to however many elements there are,
   * initialize all entries to 5.0 [W]. */
  if (numElems > 0) {
    elementHeatRate = (double*) malloc (numElems * sizeof(double));
    for (size_t e = 0; e < numElems; e += 1) {
      elementNodeCounts[e] = 4;
      elementHeatRate[e] = 5.0;
    }
  }

  /* Resize node temperature to however many nodes there are,
   * initialize all entries to 300.0 [K]. */
  if (numNodes > 0) {
    nodeTemperature = (double*) malloc (numNodes * sizeof(double));
    for (size_t n = 0; n < numNodes; n += 1) {
      nodeCoordsZ[n] = 0.0;
      nodeTemperature[n] = 300.0;
    }
  }
}

/* Function that returns surface mesh, given a region. */
SyscSurfaceMesh getSurfaceMesh(
    const char* regionName)
{
  if (strcmp(regionName, "FSI") == 0) {

    /* node coords */
    SyscOutputVectorData nodeCoordsData = syscGetOutputVectorDataSplitDouble(
        &nodeCoordsX[0],
        &nodeCoordsY[0],
        &nodeCoordsZ[0],
        numNodes);

    /* elem node counts */
    SyscOutputIntegerData elemNodeCountsData = syscGetOutputIntegerDataUInt16(
        &elementNodeCounts[0],
        numElems);

    /* elem node ids */
    SyscOutputIntegerData elemNodeIdsData = syscGetOutputIntegerDataUInt64(
        &elementNodeIds[0],
        numElems * 4);

    /* nodes */
    SyscNodeData nodes = syscGetNodeDataC(nodeCoordsData);

    /* faces */
    SyscFaceData faces = syscGetFaceDataCN(
      syscGetElementNodeCountData(elemNodeCountsData),
      syscGetElementNodeConnectivityData(elemNodeIdsData));

    return syscGetSurfaceMeshNF(nodes, faces);
  }

  syscFatalError("getSurfaceMesh: incorrect region name");
  return syscGetSurfaceMesh();
}

/*
 * Function that returns scalar input variable values,
 * given a region and a variable.
 */
SyscInputScalarData getInputScalarData(
    const char* regionName,
    const char* variableName)
{
  SyscInputScalarData isd = syscGetInputScalarData();

  if (strcmp(regionName, "FSI") == 0) {
    if (strcmp(variableName, "Temperature") == 0) {
      isd = syscGetInputScalarDataDouble(nodeTemperature, numNodes);
    }
    else {
      syscFatalError("getInputScalarVariableData: incorrect variable name");
    }
  }
  else {
    syscFatalError("getInputScalarVariableData: incorrect region name");
  }

  return isd;
}

/* Function that returns scalar output variable values,
 * given a region and a variable.
 */
SyscOutputScalarData getOutputScalarData(
    const char* regionName,
    const char* variableName)
{
  if (strcmp(regionName, "FSI") == 0) {
    if (strcmp(variableName, "Heat Rate") == 0) {
      return syscGetOutputScalarDataDouble(elementHeatRate, numElems);
    }
  }

  syscFatalError("getOutputScalarVariableData: incorrect region or variable name");
  return syscGetOutputScalarData();
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
 * This is the main function for mock solver program.
 * It takes the following command line arguments to get
 * System Coupling host (string), port (integer),
 * and name (string) given to this participant:
 * --schost <host>
 * --scport <port>
 * --scname <name>
 *
 * This program is going to run as part of the overall coupled
 * analysis. System Coupling is going to start this program
 * automatically, so it will pass in the appropriate arguments as needed.
 */
int main(int argc, char* argv[])
{
  /* Parse and print input arguments (host, port, participant name). */
  char scHost[STRING_MAX_SIZE];
  scHost[0] = '\0';
  unsigned short scPort = 0;
  char scName[STRING_MAX_SIZE];
  scName[0] = '\0';
  int i = 1;
  for (i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "--schost") == 0 && i + 1 < argc) {
      strcpy(scHost, argv[i + 1]);
    }
    else if (strcmp(argv[i], "--scport") == 0 && i + 1 < argc) {
      scPort = (unsigned short) atoi(argv[i + 1]);
    }
    else if (strcmp(argv[i], "--scname") == 0 && i + 1 < argc) {
      strcpy(scName, argv[i + 1]);
    }
  }

  const char* buildInfo = "Channel Flow Mock Solver MPI, language: C";
  SyscError ret;

  printf("Host: %s\n", scHost);
  printf("Port: %d\n", scPort);
  printf("Name: %s\n", scName);

  /* Initialize MPI. */
  int myrank = 0;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  /* Connection to System Coupling. */
  ret = syscConnectParallel(scHost, scPort, scName, MPI_COMM_WORLD, buildInfo);
  processRetCode(ret);

  /*
  * Heavyweight data access and initialization.
  *
  * Register access to mesh and variable data.
  * NOTE: The participant solver retains ownership
  * of mesh and variable field data arrays and is thus
  * responsible for allocating and deallocating memory
  * for those arrays. System Coupling only gets
  * read and/or write access to those arrays.
  */

  /* Register surface mesh access by providing pointer
  * to a function that provides mesh information for a given region.
  */
  ret = syscRegisterSurfMeshAccess(&getSurfaceMesh);
  processRetCode(ret);

  /* Register variable data access by providing pointer
  * to a function that provides access to variable data for
  * a region.
  */

  /* Input scalar variable access. */
  ret = syscRegisterInputScalarDataAccess(&getInputScalarData);
  processRetCode(ret);

  /* Output scalar variable access. */
  ret = syscRegisterOutputScalarDataAccess(&getOutputScalarData);
  processRetCode(ret);

  /* Initialize the arrays. */
  initializeArrays();

  /* Coupled analysis initialization. */
  ret = syscInitializeAnalysis();
  processRetCode(ret);

  /* Coupled analysis loop. */
  while (syscDoIteration() == 1) {

    /*
    * Inputs update.
    * All input variables will be up-to-date after this call.
    */
    ret = syscUpdateInputs();
    processRetCode(ret);

    /* Solver iterations can be performed here.
    * Input variable values can be consumed.
    * Output variable values must be generated.
    * ...
    */

    /*
    * Outputs update.
    * Provide convergence state to System Coupling,
    * and notify that outputs can be consumed.
    */
    ret = syscUpdateOutputs(SyscConverged);
    processRetCode(ret);

  }

  /* deallocate arrays */
  if (numNodes > 0) {
    free(nodeCoordsX);
    free(nodeCoordsY);
    free(nodeCoordsZ);
    free(nodeTemperature);
  }

  if (numElems > 0) {
    free(elementNodeCounts);
    free(elementNodeIds);
    free(elementHeatRate);
  }

  /* Coupled analysis shutdown. */
  ret = syscDisconnect();
  processRetCode(ret);

  MPI_Finalize();
  return EXIT_SUCCESS;
}
