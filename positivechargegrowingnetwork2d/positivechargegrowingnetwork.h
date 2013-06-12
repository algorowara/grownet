#ifndef POSITIVECHARGEGROWINGNETWORK2D_H
#define POSITIVECHARGEGROWINGNETWORK2D_H

#include "../graph/growingnetwork.h"
#include "../growingnetwork3d/spatialvertex.h"
#include <cmath>
#include <cstdlib>
#include <cfloat>
#include <ctime>
#include <omp.h>
#include <iostream>

#define DIM(2)
#define X(node) node->position[0]
#define Y(node) node->position[1]
#define DISTANCE_SQUARED(a, b)  ((X(b) - X(a)) * (X(b) - X(a)) + (Y(b) - Y(a)) * (Y(b) - Y(a)))
#define DISTANCE(a, b) sqrt(DISTANCE_SQUARED(a, b))

using namespace std;

class PositiveChargeGrowingNetwork2D : public GrowingNetwork{

public:

	vector<SpatialVertex*> nodes;
	double radius;

	PositiveChargeGrowingNetwork2D(long int n, long int m);
	void grow(long int n);
	double* randomLocation();
	double linearDistance(Vertex* a, Vertex* b);
	double* sumForces(SpatialVertex* node);
	SpatialVertex** findMNearestNeighbors(SpatialVertex* start);
	double calculatePotential();
	void equalize();
	void gradientDescent();

};

#endif
