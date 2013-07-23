#ifndef NSPHERE_H
#define NSPHERE_H

#include "ngraph.h"

#define NSPHERE_DEFAULT_RADIUS 1.0
#define NSPHERE_DEFAULT_GAMMA 6.0
#define NSPHERE_DEFAULT_TOLERANCE 0.01
#define NSPHERE_DEFAULT_ITERATIONS 10000
#define NSPHERE_DEFAULT_THRESHOLD 100
#define NSPHERE_DEFAULT_PERIOD 1

#define NSPHERE_LINEAR_DISTANCE(nodePos, otherPos, dist) \
	dist = 0; \
	for(long int itr_var = 0; itr_var < DIM+1; itr_var++){ \
		dist += (nodePos[itr_var] - otherPos[itr_var]) * (nodePos[itr_var] - otherPos[itr_var]); \
	} \
	dist = sqrt(dist); \

using namespace std;

class NSphere : public NGraph {

	public:
	
		long int forceExp;
		
		NSphere(long int n, long int m, int d, float r = NSPHERE_DEFAULT_RADIUS, float baseGam = NSPHERE_DEFAULT_GAMMA, float baseTol = NSPHERE_DEFAULT_TOLERANCE, long int baseItr = NSPHERE_DEFAULT_ITERATIONS, long int threshold = NSPHERE_DEFAULT_THRESHOLD, long int period = NSPHERE_DEFAULT_PERIOD);
		void grow(long int n);
		float* randomLocation();
		float linearDistance(Vertex* a, Vertex* b);
		float linearDistance(float* a, float* b);
		SpatialVertex** findMNearestNeighbors(SpatialVertex* start);
		float* sumForces(SpatialVertex* node);
		void equalize();
		void gradientDescent(float gamma, float tolerance, long int maxItr);
		void normalizeRadius(SpatialVertex* node);
		static void exportObject(const NSphere* ns, const char* filename);
		static NSphere* importObject(const char* filename);
	
};



#endif
