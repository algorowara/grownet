#ifndef GROWINGNETWORK3D_H
#define GROWINGNETWORK3D_H

#include "../graph/growingnetwork.h"

/**
 * though the nodes inhabit three-dimensional space,
 * they do so on the surface of a sphere,
 * making the radial component of their sphereical coordinates
 * trivial/common/irrelevant
 */
#define DIM (2)

using namespace std;

class GrowingNetwork3D : public GrowingNetwork{
	
public:

	vector<SpatialVertex*> neighbors;

	GrowingNetwork3D(long int n, long int m);
	void grow(long int n);
	double* randomLocation();
	double* displacement(SpatialVertex* a, SpatialVertex* b);
	double* force(double* disp);
	SpatialVertex** findMNearestNeighbors(SpatialVertex* start);
	void equalize();
	
};

#endif
