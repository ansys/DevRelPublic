/*
 * Copyright ANSYS, Inc. Unauthorized use, distribution, or duplication is prohibited.
 */

#include "PipeMeshGenerator.hpp"

#include "SystemCouplingParticipant/SystemCoupling.hpp"

#include <mpi.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// Define parallel properties.
int myrank = 0;
int numranks = 1;

/* Function that returns surface mesh, given a region. */
sysc::SurfaceMesh getSurfaceMesh(
  const std::string& regionName)
{
  return sysc::SurfaceMesh(
    sysc::NodeData(sysc::OutputVectorData(getNodeCoords(regionName))),
    sysc::FaceData(
      sysc::ElementNodeCountData(sysc::OutputScalarData(getElemNodeCounts(regionName))),
      sysc::ElementNodeConnectivityData(sysc::OutputScalarData(getElemNodeIds(regionName)))));
}

/*    
 * Function that returns scalar input variable values,
 * given a region and a variable.
 */
sysc::InputScalarData getInputScalarData(
  const std::string& regionName,
  const std::string& variableName)
{
  return sysc::InputScalarData(getSolutionData(regionName, variableName));
}

/*
 * Function that returns scalar output variable values,
 * given a region and a variable.
 */
sysc::OutputScalarData getOutputScalarData(
  const std::string& regionName,
  const std::string& variableName)
{
  return sysc::OutputScalarData(getSolutionData(regionName, variableName));
}

/* Function to calculate and print max error on a given region. */
void printMaxError(const std::string& regionName)
{
  std::string expectedVar = regionName == "quad" ? "linear1" : "linear2";
  std::string testVar = expectedVar == "linear2" ? "linear1" : "linear2";
  const std::vector<double>& expectedValues = getSolutionData(regionName, expectedVar);
  const std::vector<double>& testValues = getSolutionData(regionName, testVar);
  std::size_t numVals = testValues.size();
  double maxError = 0.0, totalMaxError;
  for (std::size_t n = 0; n < numVals; ++n) {
    const double& testValue = testValues[n];
    const double& expectedValue = expectedValues[n];
    maxError = std::fmax(maxError, std::abs(expectedValue - testValue));
  }

  MPI_Reduce(&maxError, &totalMaxError, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (myrank == 0) {
    std::cout << "Maximum error on " << regionName << ": " << totalMaxError << "\n";
  }
}

/* Function to write out mesh and solution data arrays to text files. */
void writeArrays(const std::string& regionName)
{
  const std::vector<double>& currNodeCoords = getNodeCoords(regionName);
  const std::vector<int>& currElemNodeCounts = getElemNodeCounts(regionName);
  const std::vector<int>& currElemNodeIds = getElemNodeIds(regionName);

  std::string fileName = regionName + "." + std::to_string(myrank) + ".dat";
  std::ofstream fileHandle(fileName.c_str());
  fileHandle << std::setprecision(8) << std::scientific;

  // Write number of nodes, elements, and solution data arrays.
  int numNodes = static_cast<int>(currNodeCoords.size() / 3);
  int numElems = static_cast<int>(currElemNodeCounts.size());
  int numScalarNodeArrays = 2;  // "linear1" and "linear2"
  fileHandle << numNodes << ", "
             << numElems << ", "
             << numScalarNodeArrays << "\n";

  // Write node coordinates.
  for (int n = 0; n < numNodes; ++n) {
    fileHandle << currNodeCoords[n * 3 + 0] << ", ";
    fileHandle << currNodeCoords[n * 3 + 1] << ", ";
    fileHandle << currNodeCoords[n * 3 + 2] << "\n";
  }

  // Write element-to-node connectivity.
  std::size_t elemNodeIndex = 0;
  for (int e = 0; e < numElems; ++e) {
    // Write the first node id for the current element.
    fileHandle << currElemNodeIds[elemNodeIndex];
    ++elemNodeIndex;
    // Write subsequent node ids for the current element.
    for (int en = 1; en < currElemNodeCounts[e]; ++en) {
      fileHandle << ", " << currElemNodeIds[elemNodeIndex];
      ++elemNodeIndex;
    }
    fileHandle << "\n";
  }

  // Write solution data arrays.
  std::vector<std::string> variableNames(2);
  variableNames[0] = "linear1";
  variableNames[1] = "linear2";
  for (const std::string& variableName : variableNames) {
    const std::vector<double>& currData = getSolutionData(regionName, variableName);
    fileHandle << variableName << "\n";
    for (int n = 0; n < numNodes; ++n) {
      fileHandle << currData[n] << "\n";
    }
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
  try {
    // Parse command line arguments.
    bool write = true;

    PipeMeshProperties& meshProperties = getMeshProperties();

    for (int i = 1; i < argc; ++i) {
      if (std::string(argv[i]) == "--quad" && i + 1 < argc) {
        meshProperties.quadRefine = std::stoi(argv[i + 1]);
      }
      else if (std::string(argv[i]) == "--tri" && i + 1 < argc) {
        meshProperties.triRefine = std::stoi(argv[i + 1]);
      }
      else if (std::string(argv[i]) == "--offset" && i + 1 < argc) {
        meshProperties.pipeOffset = std::stod(argv[i + 1]);
      }
      else if (std::string(argv[i]) == "--qradius" && i + 1 < argc) {
        meshProperties.quadPipeRadius = std::stod(argv[i + 1]);
      }
      else if (std::string(argv[i]) == "--tradius" && i + 1 < argc) {
        meshProperties.triPipeRadius = std::stod(argv[i + 1]);
      }
      else if (std::string(argv[i]) == "--write" && i + 1 < argc) {
        write = static_cast<bool>(std::stoi(argv[i + 1]));
      }
      else if (std::string(argv[i]) == "--overlap" && i + 1 < argc) {
        meshProperties.overlapLayers = std::stoi(argv[i + 1]);
      }
    }

    // Initialize MPI.
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numranks);

    meshProperties.myrank = myrank;
    meshProperties.numranks = numranks;

    if (myrank == 0) {
      std::cout << "Setup:\n";
      std::cout << "Quad refinement: " << meshProperties.quadRefine << '\n';
      std::cout << "Tri refinement: " << meshProperties.triRefine << '\n';
      std::cout << "Pipe offset: " << meshProperties.pipeOffset << '\n';
      std::cout << "Quad radius: " << meshProperties.quadPipeRadius << '\n';
      std::cout << "Tri radius: " << meshProperties.triPipeRadius << '\n';
      std::cout << "Write: " << write << '\n';
      std::cout << "Overlap layers: " << meshProperties.overlapLayers << '\n';
      std::cout << "Processes: " << numranks << '\n';
      std::cout << "==========\n";
    }

    // Create SystemCoupling object.
    sysc::SystemCoupling sc;

    // Add variables, regions, and the interface.
    sysc::Variable variableLinear1("linear1", sysc::Scalar, false, sysc::Node);
    sysc::Variable variableLinear2("linear2", sysc::Scalar, false, sysc::Node);

    sysc::Region quadRegion("quad", sysc::Surface);
    quadRegion.addOutputVariable(variableLinear1);
    quadRegion.addInputVariable(variableLinear2);

    sysc::Region triRegion("tri", sysc::Surface);
    triRegion.addOutputVariable(variableLinear2);
    triRegion.addInputVariable(variableLinear1);

    sysc::CouplingInterface couplingInterface("interface");
    couplingInterface.addSideOneRegion(quadRegion);
    couplingInterface.addSideTwoRegion(triRegion);

    sc.addCouplingInterface(couplingInterface);

    // Register data access functions.
    sc.registerSurfaceMeshAccess(&getSurfaceMesh);
    sc.registerInputScalarDataAccess(&getInputScalarData);
    sc.registerOutputScalarDataAccess(&getOutputScalarData);

    double startTime, localTime, totalTime;
    std::cout << std::setprecision(4) << std::scientific;

    // Initialize mesh and solution data on both sides of the interface.
    startTime = MPI_Wtime();
    initializeMesh();
    localTime = MPI_Wtime() - startTime;
    MPI_Reduce(&localTime, &totalTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (myrank == 0) {
      std::cout << "Time to generate mesh: " << totalTime << " [s]\n";
    }

    // Update inputs - this will perform mapping and bring inputs up-to-date.
    startTime = MPI_Wtime();
    sc.updateInputs();
    localTime = MPI_Wtime() - startTime;
    MPI_Reduce(&localTime, &totalTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (myrank == 0) {
      std::cout << "Time to update inputs: " << totalTime << " [s]\n";
    }

    // Write results, if desired.
    if (write) {
      startTime = MPI_Wtime();
      writeArrays("quad");
      writeArrays("tri");
      sc.writeResults(sysc::ResultsInfo("pipe"));
      localTime = MPI_Wtime() - startTime;
      MPI_Reduce(&localTime, &totalTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      if (myrank == 0) {
        std::cout << "Time to write results: " << totalTime << " [s]\n";
      }
    }

    printMaxError("quad");
    printMaxError("tri");
  }
  catch (const std::exception& exc) {
    std::cout << "ERROR: " << exc.what() << std::endl;
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  MPI_Finalize();
  return EXIT_SUCCESS;
}
