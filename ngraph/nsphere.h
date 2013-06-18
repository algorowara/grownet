#ifndef NSPHERE_H
#define NSPHERE_H

#include "../growingnetwork3d/growingnetwork3d.h"

class NSphere : public GrowingNetwork3D {

	public:
	
		final long int DIM;
		
		double* randomLocation();
		double distanceSquared(SpatialVertex* a, SpatialVertex* b);
		double distance(SpatialVertex* a, SpatialVertex* b);
	
}



#endif
