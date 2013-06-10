#ifndef GROWINGNETWORK3D_H
#define GROWINGNETWORK3D_H

#include "../graph/growingnetwork.h"

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
	double* calculateAngularDisplacement(SpatialVertex* a, SpatialVertex* b);
	double* sumForces(SpatialVertex* node);
	SpatialVertex** findMNearestNeighbors(SpatialVertex* start);
	void equalize();
	void normalizeRaidus(SpatialVertex* node);
	
};

#endif
