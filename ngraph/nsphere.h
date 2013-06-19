#ifndef NSPHERE_H
#define NSPHERE_H

#include "../growingnetwork3d/growingnetwork3d.h"

#define GROWTYPE_RANDOM 0
#define GROWTYPE_BUD 1

class NSphere : public GrowingNetwork3D {

	public:
	
		long int DIM;
		long int growtype;
		
		NSphere(long int dim, long int n, long int m, long int growtype = GROWTYPE_RANDOM, double baseGam = 1.0, double baseTol = 0.01, long int baseItr = 100);
		double* randomLocation();
		double distanceSquared(SpatialVertex* a, SpatialVertex* b);
		double distance(SpatialVertex* a, SpatialVertex* b);
		double* sumForces(SpatialVertex* node);
		void normalizeRadius(SpatialVertex* node);
		SpatialVertex* bud(SpatialVertex* source, double dist);
	
};



#endif
