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
#define X(node) node->position[0]
#define Y(node) node->position[1]
#define Z(node) node->position[2]

#define DEFAULT_GAMMA 1.0
#define DEFAULT_TOLERANCE 0.01
#define DEFAULT_ITR 36

using namespace std;

class GrowingNetwork3D : public GrowingNetwork{
	
public:

	const static long int DIM = 3;
	double radius;
	double baseGam;
	double baseTol;
	long int baseItr;

	GrowingNetwork3D(long int n, long int m, double gam = DEFAULT_GAMMA, double tol = DEFAULT_TOLERANCE, long int itr = DEFAULT_ITR);
	SpatialVertex* getNode(long int i);
	void grow(long int n);
	double* randomLocation();
	double distanceSquared(SpatialVertex* a, SpatialVertex* b);
	double distance(SpatialVertex* a, SpatialVertex* b);
	double linearDistance(Vertex* a, Vertex* b);
	double* sumForces(SpatialVertex* node);
	SpatialVertex** findMNearestNeighbors(SpatialVertex* start);
	void normalizeRadius(SpatialVertex* node);
	double calculatePotential();
	void equalize();
	void gradientDescent(double gamma, double baseTolerance, long int maxItr);
	static double calculateMinimumPotential(long int n);
	static double calculateMinimumPotentialDifference(long int init_n, long int final_n);
	~GrowingNetwork3D();
	
};

#endif
