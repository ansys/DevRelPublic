/*
 * Copyright ANSYS, Inc. Unauthorized use, distribution, or duplication is prohibited.
 */

#include "SystemCouplingParticipant/SystemCoupling.hpp"

#include <mpi.h>

#include <iostream>
#include <string>
#include <vector>

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

// Declare arrays used for this problem.
// The contents are populated in the initializeArray()
// function, depending on parallel settings.

// Array of node coordinates in split array format.
std::vector<double> nodeCoordsX;
std::vector<double> nodeCoordsY;
std::vector<double> nodeCoordsZ;

// Array of element node counts.
std::vector<uint16_t> elementNodeCounts;

// Array of node ids for each element.
std::vector<uint64_t> elementNodeIds;

// Array of temperature values (defined at nodes).
std::vector<double> nodeTemperature;

// Array of heat rate values (defined at elements)
// values are in Watts - total of 20 [W] generated
std::vector<double> elementHeatRate;

// Function that partitions the mesh based
// on the number of parallel processes, and
// initializes the arrays.
void initializeArrays()
{
  // Find out total number of ranks and local rank.
  int myrank = 0;
  int numranks = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &numranks);

  std::size_t numNodes = 0;

  // Partition the mesh based on number of parallel processes.
  if (numranks == 1) {
    // All four faces are present on rank 0.

    numNodes = 10;
    nodeCoordsX.resize(numNodes, 0);
    nodeCoordsY.resize(numNodes, 0);
    nodeCoordsZ.resize(numNodes, 0);
    elementNodeCounts.resize(4, 4);
    elementNodeIds.resize(16, 0);

    nodeCoordsX[0] = 0.00;
    nodeCoordsY[0] = 0.0;
    nodeCoordsX[1] = 0.25;
    nodeCoordsY[1] = 0.0;
    nodeCoordsX[2] = 0.50;
    nodeCoordsY[2] = 0.0;
    nodeCoordsX[3] = 0.75;
    nodeCoordsY[3] = 0.0;
    nodeCoordsX[4] = 1.00;
    nodeCoordsY[4] = 0.0;
    nodeCoordsX[5] = 0.00;
    nodeCoordsY[5] = 0.1;
    nodeCoordsX[6] = 0.25;
    nodeCoordsY[6] = 0.1;
    nodeCoordsX[7] = 0.50;
    nodeCoordsY[7] = 0.1;
    nodeCoordsX[8] = 0.75;
    nodeCoordsY[8] = 0.1;
    nodeCoordsX[9] = 1.00;
    nodeCoordsY[9] = 0.1;

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
      // Faces A and B are on rank 0.

      numNodes = 6;
      nodeCoordsX.resize(numNodes, 0);
      nodeCoordsY.resize(numNodes, 0);
      nodeCoordsZ.resize(numNodes, 0);
      elementNodeCounts.resize(2, 4);
      elementNodeIds.resize(8, 0);

      nodeCoordsX[0] = 0.00;
      nodeCoordsY[0] = 0.0;
      nodeCoordsX[1] = 0.25;
      nodeCoordsY[1] = 0.0;
      nodeCoordsX[2] = 0.50;
      nodeCoordsY[2] = 0.0;
      nodeCoordsX[3] = 0.00;
      nodeCoordsY[3] = 0.1;
      nodeCoordsX[4] = 0.25;
      nodeCoordsY[4] = 0.1;
      nodeCoordsX[5] = 0.50;
      nodeCoordsY[5] = 0.1;

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
      // Faces C and D are on rank 1.

      numNodes = 6;
      nodeCoordsX.resize(numNodes, 0);
      nodeCoordsY.resize(numNodes, 0);
      nodeCoordsZ.resize(numNodes, 0);
      elementNodeCounts.resize(2, 4);
      elementNodeIds.resize(8, 0);

      nodeCoordsX[0] = 0.50;
      nodeCoordsY[0] = 0.0;
      nodeCoordsX[1] = 0.75;
      nodeCoordsY[1] = 0.0;
      nodeCoordsX[2] = 1.00;
      nodeCoordsY[2] = 0.0;
      nodeCoordsX[3] = 0.50;
      nodeCoordsY[3] = 0.1;
      nodeCoordsX[4] = 0.75;
      nodeCoordsY[4] = 0.1;
      nodeCoordsX[5] = 1.00;
      nodeCoordsY[5] = 0.1;

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
      // Face A is on rank 0.

      numNodes = 4;
      nodeCoordsX.resize(numNodes, 0);
      nodeCoordsY.resize(numNodes, 0);
      nodeCoordsZ.resize(numNodes, 0);
      elementNodeCounts.resize(1, 4);
      elementNodeIds.resize(4, 0);

      nodeCoordsX[0] = 0.00;
      nodeCoordsY[0] = 0.0;
      nodeCoordsX[1] = 0.25;
      nodeCoordsY[1] = 0.0;
      nodeCoordsX[2] = 0.00;
      nodeCoordsY[2] = 0.1;
      nodeCoordsX[3] = 0.25;
      nodeCoordsY[3] = 0.1;

      elementNodeIds[0] = 1;
      elementNodeIds[1] = 2;
      elementNodeIds[2] = 7;
      elementNodeIds[3] = 6;
    }
    else if (myrank == 1) {
      // Face B is on rank 1.

      numNodes = 4;
      nodeCoordsX.resize(numNodes, 0);
      nodeCoordsY.resize(numNodes, 0);
      nodeCoordsZ.resize(numNodes, 0);
      elementNodeCounts.resize(1, 4);
      elementNodeIds.resize(4, 0);

      nodeCoordsX[0] = 0.25;
      nodeCoordsY[0] = 0.0;
      nodeCoordsX[1] = 0.50;
      nodeCoordsY[1] = 0.0;
      nodeCoordsX[2] = 0.25;
      nodeCoordsY[2] = 0.1;
      nodeCoordsX[3] = 0.50;
      nodeCoordsY[3] = 0.1;

      elementNodeIds[0] = 2;
      elementNodeIds[1] = 3;
      elementNodeIds[2] = 8;
      elementNodeIds[3] = 7;
    }
    else {
      // Faces C and D are on rank 2.

      numNodes = 6;
      nodeCoordsX.resize(numNodes, 0);
      nodeCoordsY.resize(numNodes, 0);
      nodeCoordsZ.resize(numNodes, 0);
      elementNodeCounts.resize(2, 4);
      elementNodeIds.resize(8, 0);

      nodeCoordsX[0] = 0.50;
      nodeCoordsY[0] = 0.0;
      nodeCoordsX[1] = 0.75;
      nodeCoordsY[1] = 0.0;
      nodeCoordsX[2] = 1.00;
      nodeCoordsY[2] = 0.0;
      nodeCoordsX[3] = 0.50;
      nodeCoordsY[3] = 0.1;
      nodeCoordsX[4] = 0.75;
      nodeCoordsY[4] = 0.1;
      nodeCoordsX[5] = 1.00;
      nodeCoordsY[5] = 0.1;

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
    // Four or more processes.
    // Ranks 0, 1, 2, 3 get faces A, B, C, D, respectively.
    // Any higher ranks (if exist) have to mesh.

    if (myrank == 0) {
      // Face A is on rank 0.

      numNodes = 4;
      nodeCoordsX.resize(numNodes, 0);
      nodeCoordsY.resize(numNodes, 0);
      nodeCoordsZ.resize(numNodes, 0);
      elementNodeCounts.resize(1, 4);
      elementNodeIds.resize(4, 0);

      nodeCoordsX[0] = 0.00;
      nodeCoordsY[0] = 0.0;
      nodeCoordsX[1] = 0.25;
      nodeCoordsY[1] = 0.0;
      nodeCoordsX[2] = 0.00;
      nodeCoordsY[2] = 0.1;
      nodeCoordsX[3] = 0.25;
      nodeCoordsY[3] = 0.1;

      elementNodeIds[0] = 1;
      elementNodeIds[1] = 2;
      elementNodeIds[2] = 7;
      elementNodeIds[3] = 6;
    }
    else if (myrank == 1) {
      // Face B is on rank 1.

      numNodes = 4;
      nodeCoordsX.resize(numNodes, 0);
      nodeCoordsY.resize(numNodes, 0);
      nodeCoordsZ.resize(numNodes, 0);
      elementNodeCounts.resize(1, 4);
      elementNodeIds.resize(4, 0);

      nodeCoordsX[0] = 0.25;
      nodeCoordsY[0] = 0.0;
      nodeCoordsX[1] = 0.50;
      nodeCoordsY[1] = 0.0;
      nodeCoordsX[2] = 0.25;
      nodeCoordsY[2] = 0.1;
      nodeCoordsX[3] = 0.50;
      nodeCoordsY[3] = 0.1;

      elementNodeIds[0] = 2;
      elementNodeIds[1] = 3;
      elementNodeIds[2] = 8;
      elementNodeIds[3] = 7;
    }
    else if (myrank == 2) {
      // Face C is on rank 2.

      numNodes = 4;
      nodeCoordsX.resize(numNodes, 0);
      nodeCoordsY.resize(numNodes, 0);
      nodeCoordsZ.resize(numNodes, 0);
      elementNodeCounts.resize(1, 4);
      elementNodeIds.resize(4, 0);

      nodeCoordsX[0] = 0.50;
      nodeCoordsY[0] = 0.0;
      nodeCoordsX[1] = 0.75;
      nodeCoordsY[1] = 0.0;
      nodeCoordsX[2] = 0.50;
      nodeCoordsY[2] = 0.1;
      nodeCoordsX[3] = 0.75;
      nodeCoordsY[3] = 0.1;

      elementNodeIds[0] = 3;
      elementNodeIds[1] = 4;
      elementNodeIds[2] = 9;
      elementNodeIds[3] = 8;
    }
    else if (myrank == 3) {
      // Face D is on rank 3.

      numNodes = 4;
      nodeCoordsX.resize(numNodes, 0);
      nodeCoordsY.resize(numNodes, 0);
      nodeCoordsZ.resize(numNodes, 0);
      elementNodeCounts.resize(1, 4);
      elementNodeIds.resize(4, 0);

      nodeCoordsX[0] = 0.75;
      nodeCoordsY[0] = 0.0;
      nodeCoordsX[1] = 1.00;
      nodeCoordsY[1] = 0.0;
      nodeCoordsX[2] = 0.75;
      nodeCoordsY[2] = 0.1;
      nodeCoordsX[3] = 1.00;
      nodeCoordsY[3] = 0.1;

      elementNodeIds[0] = 4;
      elementNodeIds[1] = 5;
      elementNodeIds[2] = 10;
      elementNodeIds[3] = 9;
    }
    else {
      // All arrays remain empty for ranks 4 and higher.
    }
  }

  // Resize element heat rate to however many elements there are,
  // initialize all entries to 5.0 [W].
  elementHeatRate.resize(elementNodeCounts.size(), 5.0);

  // Resize node temperature to however many nodes there are,
  // initialize all entries to 300.0 [K].
  nodeTemperature.resize(numNodes, 300.0);
}

// Function that returns surface mesh, given a region.
sysc::SurfaceMesh getSurfaceMesh(
  const std::string& regionName)
{
  if (regionName == "FSI") {
    sysc::NodeData nodes(sysc::OutputVectorData(nodeCoordsX, nodeCoordsY, nodeCoordsZ));
    sysc::ElementNodeCountData faceNodeCounts(elementNodeCounts);
    sysc::ElementNodeConnectivityData faceNodeIds(elementNodeIds);
    sysc::FaceData faces(faceNodeCounts, faceNodeIds);
    sysc::SurfaceMesh mesh(nodes, faces);
    return mesh;
  }

  std::string errorMessage = "getSurfaceMesh error: ";
  errorMessage += "Unknown region " + regionName;
  throw std::runtime_error(errorMessage);
}

// Function that returns scalar input variable values,
// given a region and a variable.
sysc::InputScalarData getInputScalarData(
  const std::string& regionName,
  const std::string& variableName)
{
  if (regionName == "FSI" && variableName == "Temperature") {
    return sysc::InputScalarData(nodeTemperature);
  }

  std::string errorMessage = "getInputScalarData error: ";
  errorMessage += "unknown variable " + variableName;
  errorMessage += " on region " + regionName;
  throw std::runtime_error(errorMessage);
}

// Function that returns scalar output variable values,
// given a region and a variable.
sysc::OutputScalarData getOutputScalarData(
  const std::string& regionName,
  const std::string& variableName)
{
  if (regionName == "FSI" && variableName == "Heat Rate") {
    return sysc::OutputScalarData(elementHeatRate);
  }

  std::string errorMessage = "getOutputScalarData error: ";
  errorMessage += "unknown variable " + variableName;
  errorMessage += " on region " + regionName;
  throw std::runtime_error(errorMessage);
}

// This is the parallel implementation of the mock
// solver program.
// It takes the following command line arguments to get
// System Coupling host (string), port (integer),
// and name (string) given to this participant:
// --schost <host>
// --scport <port>
// --scname <name>
//
// This program is going to run as part of the overall coupled
// analysis. System Coupling is going to start this program
// automatically, so it will pass in the appropriate arguments as needed.
int main(int argc, char* argv[])
{
  // Parse and print input arguments (host, port, participant name).
  std::string scHost, scName;
  unsigned short int scPort = 0;
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "--schost" && i + 1 < argc) {
      scHost = std::string(argv[i + 1]);
    }
    else if (std::string(argv[i]) == "--scport" && i + 1 < argc) {
      scPort = (unsigned short int)std::stoi(argv[i + 1]);
    }
    else if (std::string(argv[i]) == "--scname" && i + 1 < argc) {
      scName = std::string(argv[i + 1]);
    }
  }

  std::cout << "Host: " << scHost << std::endl;
  std::cout << "Port: " << scPort << std::endl;
  std::cout << "Name: " << scName << std::endl;
  std::string buildInfo = "Channel Flow Mock Solver MPI, language: C++";

  int myrank = 0;

  try {
    // Initialize MPI.

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // Connection to System Coupling.

    // Create SystemCoupling object - this will establish the
    // connection to System Coupling in parallel.
    sysc::SystemCoupling sc(scHost, scPort, scName, MPI_COMM_WORLD, buildInfo);

    // Heavyweight data access and initialization.

    // Register access to mesh and variable data.
    // NOTE: The participant solver retains ownership
    // of mesh and variable field data arrays and is thus
    // responsible for allocating and deallocating memory
    // for those arrays. System Coupling only gets
    // read and/or write access to those arrays.

    // Register surface mesh access by providing pointer
    // to a function that returns SurfaceMesh object corresponding
    // to the requested region.
    sc.registerSurfaceMeshAccess(&getSurfaceMesh);

    // Register variable data access by providing pointer
    // to a function that provides access to variable data for
    // a region.

    // Output scalar variable access.
    sc.registerOutputScalarDataAccess(&getOutputScalarData);

    // Input scalar variable access.
    sc.registerInputScalarDataAccess(&getInputScalarData);

    // Initialize the arrays.
    initializeArrays();

    // Coupled analysis initialization.
    sc.initializeAnalysis();

    // Coupled analysis loop.
    while (sc.doIteration()) {
      // Inputs update.
      // All input variables will be up-to-date after this call.
      sc.updateInputs();

      // Solver iterations can be peformed here.
      // Input variable values can be consumed.
      // Output variable values must be generated.
      // ...

      // Outputs update.
      // Provide convergence state to System Coupling,
      // and notify that outputs can be consumed.
      sc.updateOutputs(sysc::Converged);
    }

    // Coupled analysis shutdown.
    sc.disconnect();
  }
  catch (const std::exception& exc) {
    std::cout << myrank << " ERROR: " << exc.what() << std::endl;
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  MPI_Finalize();
  return EXIT_SUCCESS;
}
