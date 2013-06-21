#ifndef NBALL_H
#define NBALL_H

#include "../graph/growingnetwork.h"

class NBall : public GrowingNetwork {
	
	public:
	
		const long int DIM;
		double radius;
		double baseAlpha;	// node-node force constant
		double baseBeta;	// attractive cloud force constant
		double baseGamma;	// size of timestep per iteration of GradientDescent
		double baseTolerance;
		long int baseItr;
	
}

#endif
