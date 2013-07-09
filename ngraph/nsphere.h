#ifndef NSPHERE_H
#define NSPHERE_H

#include "../growingnetwork3d/growingnetwork3d.h"

#define NSPHERE_DEFAULT_RADIUS 1.0
#define NSPHERE_DEFAULT_GAMMA 1.0
#define NSPHERE_DEFAULT_TOLERANCE 0.001
#define NSPHERE_DEFAULT_ITERATIONS 36
#define NSPHERE_DEFAULT_THRESHOLD 100
#define NSPHERE_DEFAULT_PERIOD 1

using namespace std;

class NSphere : public GrowingNetwork {

	public:
	
		const long int DIM;
		double radius;
		double baseGam;
		double baseTol;
		double baseItr;
		long int equalizationThreshold;	// the number of nodes below which GradientDescent is called once per node added
		long int equalizationPeriod;	// the number of nodes added per call to GradientDescent
		double iterationWeights;	// count of weighted iterations over this object's lifetime
		
		NSphere(long int n, long int m, long int d, double baseGam = NSPHERE_DEFAULT_GAMMA, double baseTol = NSPHERE_DEFAULT_TOLERANCE, long int baseItr = NSPHERE_DEFAULT_ITERATIONS, long int threshold = NSPHERE_DEFAULT_THRESHOLD, long int period = NSPHERE_DEFAULT_PERIOD);
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
