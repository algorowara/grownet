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
#define DIST_SQUARED(d) (d[0] * d[0] + d[1] * d[1] + d[2] * d[2])

using namespace std;

class GrowingNetwork3D : public GrowingNetwork{
	
public:

	vector<SpatialVertex*> nodes;
	double radius;

	GrowingNetwork3D(long int n, long int m);
	void grow(long int n);
	double* randomLocation();
	double* calculateAngularDisplacement(SpatialVertex* a, SpatialVertex* b);
	double* calculateForce(double* disp);
	SpatialVertex** findMNearestNeighbors(SpatialVertex* start);
	void equalize();
	double linearDistance(SpatialVertex* a, SpatialVertex* b);
	
};

#endif
