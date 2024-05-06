/*
 * Copyright ANSYS, Inc. Unauthorized use, distribution, or duplication is prohibited.
 */

#include "SystemCouplingParticipant/syscSystemCoupling.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*
 * See Participant Library's Oscillating Plate Damping tutorial
 * for a detailed description of the problem.
 */

#if defined(_MSC_VER)
/* Disable Microsoft compiler warnings about standard C functions. */
#pragma warning (disable : 4996)
#endif

#define STRING_MAX_SIZE 4096

/* Define surface mesh */

size_t nodeCount = 24;
size_t elemCount = 11;

/* Define node coordinates in a single array (compact form)
   There are 24 nodes. */
double nodeCoords[] = {
1.00600004e+001, 1.00000000e+000, 4.00000006e-001,
1.00600004e+001, 1.00000000e+000, 0.00000000e+000,
1.00000000e+001, 1.00000000e+000, 4.00000006e-001,
1.00000000e+001, 1.00000000e+000, 0.00000000e+000,
1.00600004e+001, 0.00000000e+000, 4.00000006e-001,
1.00600004e+001, 2.00000003e-001, 4.00000006e-001,
1.00600004e+001, 4.00000006e-001, 4.00000006e-001,
1.00600004e+001, 6.00000024e-001, 4.00000006e-001,
1.00600004e+001, 8.00000012e-001, 4.00000006e-001,
1.00600004e+001, 8.00000012e-001, 0.00000000e+000,
1.00600004e+001, 6.00000024e-001, 0.00000000e+000,
1.00600004e+001, 4.00000006e-001, 0.00000000e+000,
1.00600004e+001, 2.00000003e-001, 0.00000000e+000,
1.00600004e+001, 0.00000000e+000, 0.00000000e+000,
1.00000000e+001, 8.00000012e-001, 4.00000006e-001,
1.00000000e+001, 6.00000024e-001, 4.00000006e-001,
1.00000000e+001, 4.00000006e-001, 4.00000006e-001,
1.00000000e+001, 2.00000003e-001, 4.00000006e-001,
1.00000000e+001, 0.00000000e+000, 4.00000006e-001,
1.00000000e+001, 0.00000000e+000, 0.00000000e+000,
1.00000000e+001, 8.00000012e-001, 0.00000000e+000,
1.00000000e+001, 6.00000024e-001, 0.00000000e+000,
1.00000000e+001, 4.00000006e-001, 0.00000000e+000,
1.00000000e+001, 2.00000003e-001, 0.00000000e+000
};

/* Define element node counts (11 quad elements). */
int elementNodeCounts[] = {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

/* Define element node ids for each element (11 quad elements). */
int elementNodeIds[] = {
13, 4, 5, 12,
12, 5, 6, 11,
11, 6, 7, 10,
10, 7, 8, 9,
9, 8, 0, 1,
23, 17, 18, 19,
22, 16, 17, 23,
21, 15, 16, 22,
20, 14, 15, 21,
3, 2, 14, 20,
2, 3, 1, 0
};

/* Define array of nodal displacements
   (24 nodes, 3 components per node). */
double displacement[] = {
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

/* Define array of nodal forces
   (24 nodes, 3 components per node). */
double force[] = {
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

/* Keeps track of current time step. */
int currStep = 0;

char restartPoint[STRING_MAX_SIZE];

/* Define function that returns surface mesh, given a region. */
SyscSurfaceMesh getSurfaceMesh(
    const char* regionName)
{
  if (strcmp(regionName, "plate") == 0) {

    /* Node coordinates. */
    SyscOutputVectorData nodeCoordsData = syscGetOutputVectorDataCompactDouble(
        &nodeCoords[0],
        nodeCount);

    /* Element node counts. */
    SyscOutputIntegerData elemNodeCountsData = syscGetOutputIntegerDataInt32(
        &elementNodeCounts[0],
        elemCount);

    /* Element node ids. */
    SyscOutputIntegerData elemNodeIdsData = syscGetOutputIntegerDataInt32(
        &elementNodeIds[0],
        elemCount * 4);

    /* Nodes */
    SyscNodeData nodes = syscGetNodeDataC(nodeCoordsData);

    /* Faces */
    SyscFaceData faces = syscGetFaceDataCN(
      syscGetElementNodeCountData(elemNodeCountsData),
      syscGetElementNodeConnectivityData(elemNodeIdsData));

    return syscGetSurfaceMeshNF(nodes, faces);

  }

  syscFatalError("getSurfaceMesh: incorrect region name");
  return syscGetSurfaceMesh();
}

/*
 * Define a function that returns vector input variable values,
 * given a region and a variable.
 */
SyscInputVectorData getInputVectorData(
    const char* regionName,
    const char* variableName)
{
  if (strcmp(regionName, "plate") == 0) {
    if (strcmp(variableName, "Mesh Displacement") == 0) {
      return syscGetInputVectorDataCompactDouble(displacement, nodeCount);
    }
  }

  syscFatalError("getInputVectorData: incorrect region or variable name");
  return syscGetInputVectorData();
}

/*
 * Define a function that returns vector output variable values,
 * given a region and a variable.
 */
SyscOutputVectorData getOutputVectorData(
    const char* regionName,
    const char* variableName)
{
  if (strcmp(regionName, "plate") == 0) {
    if (strcmp(variableName, "Damping Force") == 0) {
      return syscGetOutputVectorDataCompactDouble(force, nodeCount);
    }
  }

  syscFatalError("getOutputScalarData: incorrect region or variable name");
  return syscGetOutputVectorData();
}

/*
 * Define a function that creates a restart point and returns
 * a string that identifies that restart point.
 * In this implementation, the file with up-to-date
 * nodal coordinates and force values
 * will be written and the name of that file is returned.
 * Only those values are required to perform the restart
 * for this problem.
 */
const char* createRestartPoint()
{
  char stepStr[STRING_MAX_SIZE];
  sprintf(stepStr, "%d", currStep);
  restartPoint[0] = '\0';
  strcat(restartPoint, "restart.");
  strcat(restartPoint, stepStr);
  strcat(restartPoint, ".txt");

  FILE* restartFileHandle = fopen(restartPoint, "w");

  for (size_t n = 0; n < nodeCount; n += 1) {
    fprintf(restartFileHandle, "%.8e, ", nodeCoords[n * 3 + 0]);
    fprintf(restartFileHandle, "%.8e, ", nodeCoords[n * 3 + 1]);
    fprintf(restartFileHandle, "%.8e, ", nodeCoords[n * 3 + 2]);
    fprintf(restartFileHandle, "%.8e, ", force[n * 3 + 0]);
    fprintf(restartFileHandle, "%.8e, ", force[n * 3 + 1]);
    fprintf(restartFileHandle, "%.8e\n", force[n * 3 + 2]);
  }

  fclose(restartFileHandle);
  return restartPoint;
}

void resumeState()
{
  /*
   * Error checking is skipped here for simplicity.
   * Assume the file exists and is in the format
   * consistent with how it gets written
   * in the createRestartPoint() function.
   */
  char restartPointCopy[STRING_MAX_SIZE];
  strcpy(restartPointCopy, restartPoint);

  // Set current time step.
  strtok(restartPointCopy, ".");
  currStep = atoi(strtok(NULL, "."));

  // Set node coordinates and forces from the restart file.

  FILE* restartFileHandle = fopen(restartPoint, "r");

  for (size_t n = 0; n < nodeCount; n += 1) {
    char line[STRING_MAX_SIZE];
    fgets(line, STRING_MAX_SIZE, restartFileHandle);
    nodeCoords[n * 3 + 0] = (double) atof(strtok(line, ","));
    nodeCoords[n * 3 + 1] = (double) atof(strtok(NULL, ","));
    nodeCoords[n * 3 + 2] = (double) atof(strtok(NULL, ","));
    force[n * 3 + 0] = (double) atof(strtok(NULL, ","));
    force[n * 3 + 1] = (double) atof(strtok(NULL, ","));
    force[n * 3 + 2] = (double) atof(strtok(NULL, ","));
  }

  fclose(restartFileHandle);
}

void processRetCode(SyscError ret)
{
  if (ret.retcode != 0) {
    printf("ERROR: %s\n", ret.message);
    exit(EXIT_FAILURE);
  }
}

/*
 * This is the main function for oscillating plate damping solver program.
 *
 * It requires the following command line arguments to get
 * System Coupling host (string), port (integer),
 * name (string), and damping coefficient (double) given to this participant:
 * --schost <host>
 * --scport <port>
 * --scname <name>
 * --dampcoeff <damping coefficient>
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
  int isSetupMode = 0;
  restartPoint[0] = '\0';
  double dampCoeff = 0.0;
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
    else if (strcmp(argv[i], "--dampcoeff") == 0 && i + 1 < argc) {
      dampCoeff = (double) atof(argv[i + 1]);
    }
    else if (strcmp(argv[i], "--screstart") == 0 && i + 1 < argc) {
      strcpy(restartPoint, argv[i + 1]);
    }
    else if (strcmp(argv[i], "--scsetup") == 0) {
      isSetupMode = 1;
    }
  }

  const char* buildInfo = "Oscillating Plate Damping Solver, language: C";
  SyscError ret;

  printf("Host: %s\n", scHost);
  printf("Port: %d\n", scPort);
  printf("Name: %s\n", scName);
  if (strlen(restartPoint) != 0) {
    printf("Restart: %s\n", restartPoint);
  }
  printf("Damping coefficient: %f\n", dampCoeff);

  /* Connection to System Coupling. */
  ret = syscConnect(scHost, scPort, scName, buildInfo);
  processRetCode(ret);

  if (isSetupMode == 1) {
    /* Setup mode. */

    /* Create surface region "plate". */
    SyscRegion region = syscGetRegionT("plate", SyscSurface);

    /* Create the force variable. */
    SyscVariable forceVar = syscGetVariableTE("Damping Force", SyscVector, 1, SyscNode);

    /* Create the displacement variable. */
    SyscVariable displacementVar = syscGetVariableTE("Mesh Displacement", SyscVector, 0, SyscNode);

    /* Add force as output and displacement as input to plate. */
    ret = syscAddOutputVariable(region, forceVar);
    processRetCode(ret);    
    ret = syscAddInputVariable(region, displacementVar);
    processRetCode(ret);

    /* Add region to the setup. */
    ret = syscAddRegion(region);
    processRetCode(ret);

    /* Create setup info structure. Set analysis type to transient. 
    Set restarts supported to true. */
    int restartsSupported = 1;
    SyscSetupInfo setupInfo = syscGetSetupInfoAR(SyscTransient, restartsSupported);

    /* Complete the setup. */
    ret = syscCompleteSetup(setupInfo);
    processRetCode(ret);
  }
  else {

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

    /* Input vector variable access. */
    ret = syscRegisterInputVectorDataAccess(&getInputVectorData);
    processRetCode(ret);

    /* Output vector variable access. */
    ret = syscRegisterOutputVectorDataAccess(&getOutputVectorData);
    processRetCode(ret);

    /* Register restart point creation function. */
    ret = syscRegisterRestartPointCreation(&createRestartPoint);
    processRetCode(ret);

    /* If doing a restart, resume the state from the given restart point. */
    if (strlen(restartPoint) != 0) {
      resumeState();
    }

    /* Coupled analysis initialization. */
    ret = syscInitializeAnalysis();
    processRetCode(ret);

    /* Coupling analysis loop. */

    /* Coupling time step loop. */
    while (syscDoTimeStep() == 1) {

      /* Increment time step. */
      currStep += 1;

      /* Update nodal coordinates. */
      size_t n = 0;
      for (n = 0; n < nodeCount * 3; ++n) {
        nodeCoords[n] += displacement[n];
      }

      /* Coupling iteration loop. */
      while (syscDoIteration() == 1) {

        /*
        * Inputs update.
        * All input variables will be up-to-date after this call.
        */
        ret = syscUpdateInputs();
        processRetCode(ret);

        /* Get time step size, which is used to caclulate
        * nodal velocities.
        */
        SyscTimeStep timeStep = syscGetCurrentTimeStep();
        double dt = timeStep.timeStepSize;

        /* Calculate the damping force. */
        for (n = 0; n < nodeCount * 3; n = n + 1) {
          force[n] = - dampCoeff * displacement[n] / dt;
        }

        /*
        * Outputs update.
        * Note: convergence status is Complete, since this is
        * not an iterative solver.
        */
        ret = syscUpdateOutputs(SyscComplete);
        processRetCode(ret);

      }

    }

  }

  /* Coupled analysis shutdown. */
  ret = syscDisconnect();
  processRetCode(ret);

  return EXIT_SUCCESS;
}
