/*
 * Copyright ANSYS, Inc. Unauthorized use, distribution, or duplication is prohibited.
 */

#include "SystemCouplingParticipant/SystemCoupling.hpp"

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
*          |    |    |    |    |
*          |    |    |    |    |
*          1----2----3----4----5
* (0,0,0)                       (1,0,0)
*/

// Array of node coordinates in split array format.
std::vector<double> nodeCoordsX{0, 0.25, 0.5, 0.75, 1, 0, 0.25, 0.5, 0.75, 1};
std::vector<double> nodeCoordsY{0, 0, 0, 0, 0, 0.1, 0.1, 0.1, 0.1, 0.1};
std::vector<double> nodeCoordsZ{0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

// Array of element node counts.
std::vector<std::uint16_t> elementNodeCounts{4, 4, 4, 4};

// Array of node ids for each element.
std::vector<std::uint64_t> elementNodeIds{
  1, 2, 7, 6,
  2, 3, 8, 7,
  3, 4, 9, 8,
  4, 5, 10, 9};

// Array of temperature values (defined at nodes).
std::vector<double> nodeTemperature{
  300, 300, 300, 300, 300, 300, 300, 300, 300, 300};

// Array of heat rate values (defined at elements)
// values are in Watts - total of 20 [W] generated
const std::vector<double> elementHeatRate{5, 5, 5, 5};

// Function that returns surface mesh, given a region.
sysc::SurfaceMesh getSurfaceMesh(
  const std::string& regionName)
{
  if (regionName == "FSI") {
    return sysc::SurfaceMesh(
      sysc::NodeData(sysc::OutputVectorData(nodeCoordsX, nodeCoordsY, nodeCoordsZ)),
      sysc::FaceData(
        sysc::ElementNodeCountData(elementNodeCounts),
        sysc::ElementNodeConnectivityData(elementNodeIds)));
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

// This is the main function for mock solver program.
// It requires the following command line arguments to get
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
  std::string buildInfo = "Channel Flow Mock Solver, language: C++";

  try {
    // Connection to System Coupling.

    // Create SystemCoupling object - this will establish the
    // connection to System Coupling.
    sysc::SystemCoupling sc(scHost, scPort, scName, buildInfo);

    // Coupled analysis mode.

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
    std::cout << "ERROR: " << exc.what() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
