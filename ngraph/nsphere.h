#ifndef NSPHERE_H
#define NSPHERE_H

#include "ngraph.h"

#define NSPHERE_DEFAULT_RADIUS 1.0
#define NSPHERE_DEFAULT_GAMMA 1.0
#define NSPHERE_DEFAULT_TOLERANCE 0.01
#define NSPHERE_DEFAULT_ITERATIONS 100
#define NSPHERE_DEFAULT_THRESHOLD 100
#define NSPHERE_DEFAULT_PERIOD 1

using namespace std;

class NSphere : public NGraph {

	public:
		
		NSphere(long int n, long int m, long int d, double r = NSPHERE_DEFAULT_RADIUS, double baseGam = NSPHERE_DEFAULT_GAMMA, double baseTol = NSPHERE_DEFAULT_TOLERANCE, long int baseItr = NSPHERE_DEFAULT_ITERATIONS, long int threshold = NSPHERE_DEFAULT_THRESHOLD, long int period = NSPHERE_DEFAULT_PERIOD);
		void grow(long int n);
		double* randomLocation();
		double linearDistance(Vertex* a, Vertex* b);
		double linearDistance(double* a, double* b);
		SpatialVertex** findMNearestNeighbors(SpatialVertex* start);
		double* sumForces(SpatialVertex* node);
		void equalize();
		void gradientDescent(double gamma, double tolerance, long int maxItr);
		void normalizeRadius(SpatialVertex* node);
		static void exportObject(const NSphere* ns, const char* filename);
		static NSphere* importObject(const char* filename);
	
};



#endif
