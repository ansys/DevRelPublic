/*
 * Copyright ANSYS, Inc. Unauthorized use, distribution, or duplication is prohibited.
 */

#pragma once

#include <string>
#include <vector>

// Define structure to hold mesh properties.
struct PipeMeshProperties {
  int quadRefine;
  int triRefine;
  double pipeOffset;
  double quadPipeRadius;
  double triPipeRadius;
  int overlapLayers;
  double pipeLength;
  int myrank;
  int numranks;

  PipeMeshProperties() :
      quadRefine(20),
      triRefine(30),
      pipeOffset(0.0),
      quadPipeRadius(0.05),
      triPipeRadius(0.05),
      overlapLayers(0),
      pipeLength(1.0),
      myrank(0),
      numranks(1)
  {
  }
};

PipeMeshProperties& getMeshProperties();

void initializeMesh();

const std::vector<double>& getNodeCoords(const std::string& regionName);

const std::vector<int>& getElemNodeCounts(const std::string& regionName);

const std::vector<int>& getElemNodeIds(const std::string& regionName);

std::vector<double>& getSolutionData(
  const std::string& regionName, const std::string& variableName);
