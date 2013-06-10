#ifndef GROWINGNETWORK2D_H
#define GROWINGNETWORK2D_H

#include "../graph/growingnetwork.h"
#include <omp.h>

using namespace std;

class GrowingNetwork2D : public GrowingNetwork {

public:

	GrowingNetwork2D(long int n, long int m);
	void grow(long int n);
	long int edgeArcDistance(Vertex* a, Vertex* b);
	double* edgeAgeVsArcDistance();
	double linearDistance(Vertex* a, Vertex* b);
	double* edgeArcDistanceDistribution();

};

#endif
