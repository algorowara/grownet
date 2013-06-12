#ifndef GROWINGNETWORK3D_H
#define GROWINGNETWORK3D_H

#include "../graph/growingnetwork.h"
#include "../growingnetwork3d/spatialvertex.h"
#include <cmath>
#include <cstdlib>
#include <cfloat>
#include <ctime>
#include <omp.h>
#include <iostream>

/**
 * for specifications of the coordinate system,
 * please see doc.odt
 */
#define DIM (3)
#define X(node) node->position[0]
#define Y(node) node->position[1]
#define Z(node) node->position[2]
#define DISTANCE_SQUARED(a, b) ((X(b) - X(a)) * (X(b) - X(a)) + (Y(b) - Y(a)) * (Y(b) - Y(a)) + (Z(b) - Z(a)) * (Z(b) - Z(a)))
#define DISTANCE(a, b) sqrt(DISTANCE_SQUARED(a, b))

using namespace std;

class GrowingNetwork3D : public GrowingNetwork{
	
public:

	vector<SpatialVertex*> nodes;
	double radius;

	GrowingNetwork3D(long int n, long int m);
	void grow(long int n);
	double* randomLocation();
	double linearDistance(Vertex* a, Vertex* b);
	double* sumForces(SpatialVertex* node);
	SpatialVertex** findMNearestNeighbors(SpatialVertex* start);
	void normalizeRadius(SpatialVertex* node);
	double calculatePotential();
	void equalize();
	void gradientDescent(double gamma, double tolerance, long int maxItr);
	
};

#endif
