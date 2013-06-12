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

#define DEFAULT_GAMMA 1.0
#define DEFAULT_TOLERANCE 0.0
#define DEFAULT_ITR 36

using namespace std;

class GrowingNetwork3D : public GrowingNetwork{
	
public:

	vector<SpatialVertex*> nodes;
	double radius;
	double baseGam;
	double baseTol;
	long int baseItr;

	GrowingNetwork3D(long int n, long int m, double gam = DEFAULT_GAMMA, double tol = DEFAULT_TOLERANCE, long int itr = DEFAULT_ITR);
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
