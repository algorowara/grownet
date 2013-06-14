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

#define DIM 2
#define X(node) node->position[0]
#define Y(node) node->position[1]
#define DISTANCE_SQUARED_2D(a, b)  ((X(b) - X(a)) * (X(b) - X(a)) + (Y(b) - Y(a)) * (Y(b) - Y(a)))
#define DISTANCE_2D(a, b) sqrt(DISTANCE_SQUARED_2D(a, b))

using namespace std;

class PositiveChargeGrowingNetwork2D : public GrowingNetwork{

public:

	vector<SpatialVertex*> nodes;
	double radius;
	double alpha; //electron-electron force constant
	double beta; //electron-cloud force constant
	double gamma;
	double tolerance;
	long int maxItr;

	PositiveChargeGrowingNetwork2D(long int n, long int m, double gamma, double tolerance, long int maxItr);
	void grow(long int n);
	double* randomLocation();
	double linearDistance(Vertex* a, Vertex* b);
	double* sumForces(SpatialVertex* node);
	SpatialVertex** findMNearestNeighbors(SpatialVertex* start);
	double calculatePotential();
	void equalize();
	void gradientDescent(double gamma, double tolerance, long int maxItr);

};

#endif
