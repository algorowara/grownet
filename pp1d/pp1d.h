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
	double radius;
	double alpha; //electron electron force constant
	double beta; //electron cloud force constant
	double gamma; //step size for gradient descent
	double tolerance; //for grad des
	long int maxItr; //for grad des

	PPGrowingNetwork1D(long int n, long int m, double gamma, double tolerance, long int maxItr);
	SpatialVertex* getNode(long int i);
	void grow(long int n);
	double* randomLocation();
	double linearDistance(Vertex* a, Vertex* b);
	double* sumForces(SpatialVertex* node);
	SpatialVertex** findMNearestNeighbors(SpatialVertex* start);
	double calculatePotential();
	void equalize();
	void gradientDescent(double gamma, double tolerance, long int maxItr);
	~PPGrowingNetwork1D();

};

#endif	
