/*
 * Copyright ANSYS, Inc. Unauthorized use, distribution, or duplication is prohibited.
 */

#include "SystemCouplingParticipant/SystemCoupling.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

/*
 * See Participant Library's Oscillating Plate Damping tutorial
 * for a detailed description of the problem.
 */

// Define surface mesh.

/*
 * Define node coordinates in a single array (compact form).
 * There are 24 nodes.
 */
std::vector<double> nodeCoords{
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
  1.00000000e+001, 2.00000003e-001, 0.00000000e+000};

// Define element node counts (11 quad elements).
const std::vector<int> elementNodeCounts{4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

// Define element node ids for each element (11 quad elements).
const std::vector<int> elementNodeIds{
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
  2, 3, 1, 0};

// Define array of nodal displacements (24 nodes, 3 components per node).
std::vector<double> displacement(24 * 3, 0.0);

// Define array of nodal forces (24 nodes, 3 components per node).
std::vector<double> force(24 * 3, 0.0);

// Keeps track of current time step.
int currStep = 0;

// Keeps track of current restart point.
std::string restartPoint;

// Define function that returns surface mesh, given a region.
sysc::SurfaceMesh getSurfaceMesh(
  const std::string& regionName)
{
  if (regionName == "plate") {
    return sysc::SurfaceMesh(
      sysc::NodeData(sysc::OutputVectorData(nodeCoords)),
      sysc::FaceData(
        sysc::ElementNodeCountData(sysc::OutputScalarData(elementNodeCounts)),
        sysc::ElementNodeConnectivityData(sysc::OutputScalarData(elementNodeIds))));
  }

  std::string errorMessage = "getSurfaceMesh error: ";
  errorMessage += "Unknown region " + regionName;
  throw std::runtime_error(errorMessage);
}

/*
 * Define a function that returns vector input variable values,
 * given a region and a variable.
 */
sysc::InputVectorData getInputVectorData(
  const std::string& regionName,
  const std::string& variableName)
{
  if (regionName == "plate" && variableName == "Mesh Displacement") {
    return sysc::InputVectorData(displacement);
  }

  std::string errorMessage = "getInputVectorData error: ";
  errorMessage += "unknown variable " + variableName;
  errorMessage += " on region " + regionName;
  throw std::runtime_error(errorMessage);
}

/*
 * Define a function that returns vector output variable values,
 * given a region and a variable.
 */
sysc::OutputVectorData getOutputVectorData(
  const std::string& regionName,
  const std::string& variableName)
{
  if (regionName == "plate" && variableName == "Damping Force") {
    return sysc::OutputVectorData(force);
  }

  std::string errorMessage = "getOutputVectorData error: ";
  errorMessage += "unknown variable " + variableName;
  errorMessage += " on region " + regionName;
  throw std::runtime_error(errorMessage);
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
std::string createRestartPoint()
{
  restartPoint = "restart." + std::to_string(currStep) + ".txt";
  std::ofstream restartFileHandle(restartPoint.c_str());
  restartFileHandle << std::setprecision(8) << std::scientific;
  const int dimensions = 3;
  const int nodeCount = static_cast<int>(nodeCoords.size() / dimensions);
  for (int n = 0; n < nodeCount; ++n) {
    restartFileHandle
      << nodeCoords[n * 3 + 0] << ", "
      << nodeCoords[n * 3 + 1] << ", "
      << nodeCoords[n * 3 + 2] << ", "
      << force[n * 3 + 0] << ", "
      << force[n * 3 + 1] << ", "
      << force[n * 3 + 2] << std::endl;
  }

  restartFileHandle.close();
  return restartPoint;
}

void resumeState()
{
  // Error checking is skipped here for simplicity.
  // Assume the file exists and is in the format
  // consistent with how it gets written
  // in the createRestartPoint() function.

  // Set current time step.
  std::string token;
  std::stringstream restartPointStream(restartPoint);
  std::getline(restartPointStream, token, '.');
  std::getline(restartPointStream, token, '.');
  currStep = std::stoi(token);

  // Set node coordinates and forces from the restart file.

  std::ifstream restartFile(restartPoint.c_str());

  const int dimensions = 3;
  const int nodeCount = static_cast<int>(nodeCoords.size() / dimensions);
  for (int n = 0; n < nodeCount; ++n) {
    std::string line;
    std::getline(restartFile, line);
    std::stringstream linestream(line);

    // read-in node coordinates
    for (int component = 0; component < 3; ++component) {
      std::getline(linestream, token, ',');
      nodeCoords[n * 3 + component] = std::stod(token);
    }
    // read-in forces
    for (int component = 0; component < 3; ++component) {
      std::getline(linestream, token, ',');
      force[n * 3 + component] = std::stod(token);
    }
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
  // Parse and print input arguments.
  std::string scHost, scName;
  unsigned short int scPort = 0;
  double dampCoeff = 0.0;
  bool isSetupMode = false;
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
    else if (std::string(argv[i]) == "--dampcoeff" && i + 1 < argc) {
      dampCoeff = std::stod(argv[i + 1]);
    }
    else if (std::string(argv[i]) == "--screstart" && i + 1 < argc) {
      restartPoint = std::string(argv[i + 1]);
    }
    else if (std::string(argv[i]) == "--scsetup") {
      isSetupMode = true;
    }
  }

  std::cout << "Host: " << scHost << std::endl;
  std::cout << "Port: " << scPort << std::endl;
  std::cout << "Name: " << scName << std::endl;
  if (!restartPoint.empty()) {
    std::cout << "Restart: " << restartPoint << std::endl;
  }
  std::cout << "Damping coefficient: " << dampCoeff << std::endl;

  std::string buildInfo = "Oscillating Plate Damping Solver, language: C++";

  try {
    // Connect to System Coupling.
    sysc::SystemCoupling sc(scHost, scPort, scName, buildInfo);

    if (isSetupMode) {
      // Setup mode.

      // Create the surface region "plate".
      sysc::Region region("plate", sysc::Surface);

      // Create the force variable.
      sysc::Variable forceVar("Damping Force", sysc::Vector, true, sysc::Node);

      // Create the displacement variable.
      sysc::Variable displacementVar("Mesh Displacement", sysc::Vector, false, sysc::Node);

      // Add force as output and displacement as input to plate.
      region.addOutputVariable(forceVar);
      region.addInputVariable(displacementVar);

      // Add region to the setup.
      sc.addRegion(region);

      // Create setup info structure. Set analysis type to transient.
      // Set restarts supported to true.
      sysc::SetupInfo setupInfo(sysc::Transient, true);

      // Complete the setup.
      sc.completeSetup(setupInfo);
    }
    else {
      // Register heavyweight data access functions.
      sc.registerSurfaceMeshAccess(&getSurfaceMesh);
      sc.registerOutputVectorDataAccess(&getOutputVectorData);
      sc.registerInputVectorDataAccess(&getInputVectorData);

      // Register restart point creation function.
      sc.registerRestartPointCreation(&createRestartPoint);

      // If doing a restart, resume the state from the given restart point.
      if (!restartPoint.empty()) {
        resumeState();
      }

      // Initialize coupled analysis.
      sc.initializeAnalysis();

      // Enter coupled analysis loop.

      // Coupling time step loop.
      while (sc.doTimeStep()) {
        // Increment time step.
        ++currStep;

        // Update nodal coordinates.
        for (std::size_t n = 0; n < displacement.size(); ++n) {
          nodeCoords[n] += displacement[n];
        }

        // Coupling iteration loop.
        while (sc.doIteration()) {
          // Inputs update.
          sc.updateInputs();

          /*
         * Get time step size, which is used to caclulate
         * nodal velocities.
         */
          sysc::TimeStep timeStep = sc.getCurrentTimeStep();
          double dt = timeStep.timeStepSize;

          // Calculate the damping force.
          for (std::size_t n = 0; n < force.size(); ++n) {
            force[n] = -dampCoeff * displacement[n] / dt;
          }

          /*
         * Outputs update.
         * Note: convergence status is Complete, since this is
         * not an iterative solver.
         */
          sc.updateOutputs(sysc::Complete);
        }
      }
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
