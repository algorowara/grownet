#ifndef NSPHERE_H
#define NSPHERE_H

#include "../growingnetwork3d/growingnetwork3d.h"

#define NSPHERE_DEFAULT_GAMMA 1.0
#define NSPHERE_DEFAULT_TOLERANCE 0.01
#define NSPHERE_DEFAULT_ITERATIONS 36

using namespace std;

class NSphere : public GrowingNetwork {

	public:
	
		const long int DIM;
		double radius;
		double baseGam;
		double baseTol;
		double baseItr;
		double iterationWeights;
		
		NSphere(long int n, long int m, long int d, double baseGam = NSPHERE_DEFAULT_GAMMA, double baseTol = NSPHERE_DEFAULT_TOLERANCE, long int baseItr = NSPHERE_DEFAULT_ITERATIONS);
		SpatialVertex* getNode(long int i);
		void grow(long int n);
		double* randomLocation();
		double linearDistance(Vertex* a, Vertex* b);
		double linearDistance(double* a, double* b);
		SpatialVertex** findMNearestNeighbors(SpatialVertex* start);
		double* sumForces(SpatialVertex* node);
		void equalize();
		void gradientDescent(double gamma, double tolerance, long int maxItr);
		void normalizeRadius(SpatialVertex* node);
	
};



#endif
