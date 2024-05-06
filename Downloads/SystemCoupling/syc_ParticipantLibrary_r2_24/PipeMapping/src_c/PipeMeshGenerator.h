/*
 * Copyright ANSYS, Inc. Unauthorized use, distribution, or duplication is prohibited.
 */

#pragma once

#include <stddef.h>

typedef struct PipeMeshProperties {
  int quadRefine;
  int triRefine;
  double pipeOffset;
  double quadPipeRadius;
  double triPipeRadius;
  int overlapLayers;
  double pipeLength;
  int myrank;
  int numranks;
} PipeMeshProperties;

PipeMeshProperties* getMeshProperties();

void initializeMesh();

void clearMesh();

const double* getNodeCoords(const char* regionName);
size_t getNumNodes(const char* regionName);

const int* getElemNodeCounts(const char* regionName);
size_t getNumElems(const char* regionName);

const int* getElemNodeIds(const char* regionName);
size_t getNumElemNodeIds(const char* regionName);

double* getSolutionData(const char* regionName, const char* variableName);
