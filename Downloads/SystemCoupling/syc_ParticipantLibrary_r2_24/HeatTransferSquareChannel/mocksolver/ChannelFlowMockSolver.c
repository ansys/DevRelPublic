/*
 * Copyright ANSYS, Inc. Unauthorized use, distribution, or duplication is prohibited.
 */

#include "SystemCouplingParticipant/syscSystemCoupling.h"

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
*          |    |    |    |    |
*          |    |    |    |    |
*          1----2----3----4----5
* (0,0,0)                       (1,0,0)
*/

/* Array of node coordinates in split array format. */
double nodeCoordsX[] = {0, 0.25, 0.5, 0.75, 1, 0, 0.25, 0.5, 0.75, 1};
double nodeCoordsY[] = {0, 0, 0, 0, 0, 0.1, 0.1, 0.1, 0.1, 0.1};
double nodeCoordsZ[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

/* Array of element node counts. */
const uint64_t elementNodeCounts[] = { 4, 4, 4, 4 };

/* Array of node ids for each element. */
const uint64_t elementNodeIds[] = { 1, 2, 7, 6,
                                        2, 3, 8, 7,
                                        3, 4, 9, 8,
                                        4, 5, 10, 9
                                      };

/* Array of temperature values (defined at nodes). */
double nodeTemperature[] = { 300, 300, 300, 300, 300,
                             300, 300, 300, 300, 300
                           };

/*
 * Array of heat rate values (defined at elements)
 * values are in Watts - total of 20 [W] generated
 */
const double elementHeatRate[] = { 5, 5, 5, 5 };

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
        10);

    /* elem node counts */
    SyscOutputIntegerData elemNodeCountsData = syscGetOutputIntegerDataUInt64(
        &elementNodeCounts[0],
        4);

    /* elem node ids */
    SyscOutputIntegerData elemNodeIdsData = syscGetOutputIntegerDataUInt64(
        &elementNodeIds[0],
        16);

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
 * Function that returns scalar input variable access,
 * given a region and a variable.
 */
SyscInputScalarData getInputScalarData(
    const char* regionName,
    const char* variableName)
{
  SyscInputScalarData isd = syscGetInputScalarData();

  if (strcmp(regionName, "FSI") == 0) {
    if (strcmp(variableName, "Temperature") == 0) {
      isd = syscGetInputScalarDataDouble(nodeTemperature, 10);
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
      return syscGetOutputScalarDataDouble(elementHeatRate, 4);
    }
  }

  syscFatalError("getOutputScalarVariableData: incorrect region or variable name");
  return syscGetOutputScalarData();
}

void processRetCode(SyscError ret)
{
  if (ret.retcode != 0) {
    printf("ERROR: %s\n", ret.message);
    exit(EXIT_FAILURE);
  }
}

/*
 * This is the main function for mock solver program.
 * It requires the following command line arguments to get
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

  const char* buildInfo = "Channel Flow Mock Solver, language: C";
  SyscError ret;

  printf("Host: %s\n", scHost);
  printf("Port: %d\n", scPort);
  printf("Name: %s\n", scName);

  /* Connection to System Coupling. */
  ret = syscConnect(scHost, scPort, scName, buildInfo);
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

  /* Coupled analysis shutdown. */
  ret = syscDisconnect();
  processRetCode(ret);

  return EXIT_SUCCESS;
}
