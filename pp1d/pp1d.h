#ifndef PP1D_H
#define PP1D_H

#include "../graph/growingnetwork.h"
#include "../growingnetwork3d/spatialvertex.h"
#include <cmath>
#include <cstdlib>
#include <cfloat>
#include <ctime>
#include <omp.h>
#include <iostream>

#define X(node) node->position[0]
#define DISTANCE_1D(a,b) abs(X(a) - X(b))

using namespace std;

class PPGrowingNetwork1D : public GrowingNetwork{

public:

	long int DIM;
	float radius;
	float alpha; //electron electron force constant
	float beta; //electron cloud force constant
	float gamma; //step size for gradient descent
	float tolerance; //for grad des
	long int maxItr; //for grad des

	PPGrowingNetwork1D(long int n, long int m, float gamma, float tolerance, long int maxItr);
	SpatialVertex* getNode(long int i);
	void grow(long int n);
	float* randomLocation();
	float linearDistance(Vertex* a, Vertex* b);
	float* sumForces(SpatialVertex* node);
	SpatialVertex** findMNearestNeighbors(SpatialVertex* start);
	float calculatePotential();
	void equalize();
	void gradientDescent(float gamma, float tolerance, long int maxItr);
	~PPGrowingNetwork1D();

};

#endif	
