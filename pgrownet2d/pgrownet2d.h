#ifndef PGROWNET2D_H
#define PGROWNET2D_H

#include "../graph/growingnetwork.h"
#include "../growingnetwork3d/spatialvertex.h"
#include <cmath>
#include <cstdlib>
#include <cfloat>
#include <ctime>
#include <omp.h>
#include <iostream>

#define X(node) node->position[0]
#define Y(node) node->position[1]
#define DISTANCE_SQUARED_2D(a, b)  ((X(b) - X(a)) * (X(b) - X(a)) + (Y(b) - Y(a)) * (Y(b) - Y(a)))
#define DISTANCE_2D(a, b) sqrt(DISTANCE_SQUARED_2D(a, b))

using namespace std;

class PositiveChargeGrowingNetwork2D : public GrowingNetwork{

public:

	long int DIM;
	float radius;
	float alpha; //electron-electron force constant
	float beta; //electron-cloud force constant
	float gamma;
	float tolerance;
	long int maxItr;

	PositiveChargeGrowingNetwork2D(long int n, long int m, float gamma, float tolerance, long int maxItr);
	SpatialVertex* getNode(long int i);
	void grow(long int n);
	float* randomLocation();
	float linearDistance(Vertex* a, Vertex* b);
	float* sumForces(SpatialVertex* node);
	SpatialVertex** findMNearestNeighbors(SpatialVertex* start);
	float calculatePotential();
	void equalize();
	void gradientDescent(float gamma, float tolerance, long int maxItr);
	~PositiveChargeGrowingNetwork2D();

};

#endif
