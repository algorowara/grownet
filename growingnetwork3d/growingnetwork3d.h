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

#define DEFAULT_RADIUS 1.0
#define DEFAULT_GAMMA 1.0
#define DEFAULT_TOLERANCE 1.0
#define DEFAULT_ITR 36

#define GUIDED_N 16

using namespace std;

class GrowingNetwork3D : public GrowingNetwork{
	
public:

	long int DIM;
	float radius;
	float baseGam;
	float baseTol;
	long int baseItr;

	GrowingNetwork3D();
	GrowingNetwork3D(long int n, long int m, float gam = DEFAULT_GAMMA, float tol = DEFAULT_TOLERANCE, long int itr = DEFAULT_ITR);
	SpatialVertex* getNode(long int i) const;
	void grow(long int n);
	virtual float* randomLocation() const;
	virtual float distanceSquared(SpatialVertex* a, SpatialVertex* b) const;
	virtual float distance(SpatialVertex* a, SpatialVertex* b) const;
	float linearDistance(Vertex* a, Vertex* b);
	virtual float* sumForces(SpatialVertex* node) const;
	SpatialVertex** findMNearestNeighbors(SpatialVertex* start) const;
	virtual void normalizeRadius(SpatialVertex* node);
	float calculatePotential() const;
	void equalize();
	void gradientDescent(float gamma, float baseTolerance, long int maxItr);
	virtual float calculateMinimumPotential(long int n, long int d) const;
	float calculateMinimumPotentialDifference(long int init_n, long int final_n, long int d) const;
	~GrowingNetwork3D();
	
};

#endif
