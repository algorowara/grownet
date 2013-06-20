#ifndef NSPHERE_H
#define NSPHERE_H

#include "../growingnetwork3d/growingnetwork3d.h"

#define GROWTYPE_RANDOM 0
#define GROWTYPE_BUD 1

using namespace std;

class NSphere : public GrowingNetwork3D {

	public:
	
		long int DIM;
		long int growtype;
		
		NSphere(long int DIM, long int n, long int m, long int growtype = GROWTYPE_RANDOM, double baseGam = 1.0, double baseTol = 0.1, long int baseItr = 100);
		double* randomLocation();
		double distanceSquared(SpatialVertex* a, SpatialVertex* b);
		double distance(SpatialVertex* a, SpatialVertex* b);
		double* sumForces(SpatialVertex* node);
		void normalizeRadius(SpatialVertex* node);
		SpatialVertex* bud(SpatialVertex* source, double dist);
		static calculateMinimumPotential(long int n, long int d);
	
};



#endif
