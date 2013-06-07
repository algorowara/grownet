#ifndef GROWINGNETWORK3D_H
#define GROWINGNETWORK3D_H

#include "../graph/growingnetwork.h"

/**
 * for specifications of the coordinate system,
 * please see doc.odt
 */
#define DIM (2)
#define THETA(node) (node->position[0])
#define PHI(node) (node->position[1])
#define DIST_SQUARED(disp) (disp[0] * disp[0] + disp[1] * disp[1])

using namespace std;

class GrowingNetwork3D : public GrowingNetwork{
	
public:

	vector<SpatialVertex*> nodes;
	double radius;

	GrowingNetwork3D(long int n, long int m);
	void grow(long int n);
	double* randomLocation();
	double* calculateDisplacement(SpatialVertex* a, SpatialVertex* b);
	double* calculateForce(double* disp);
	SpatialVertex** findMNearestNeighbors(SpatialVertex* start);
	void equalize();
	double edgeLinearDistance(SpatialVertex* a, SpatialVertex* b);
	
};

#endif
