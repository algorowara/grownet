#ifndef NBALL_H
#define NBALL_H

#include "../graph/growingnetwork.h"

#define DEFAULT_RADIUS 1.0
#define DEFAULT_GAMMA 1.0
#define DEFAULT_ALPHA 0.1
#define DEFAULT_BETA 0.1
#define DEFAULT_TOLERANCE 1.0
#define DEFAULT_ITERATIONS 36

class NBall : public GrowingNetwork {
	
	public:
	
		const long int DIM;
		double radius;
		double baseAlpha;	// node-node force constant
		double baseBeta;	// attractive cloud force constant
		double baseGamma;	// size of timestep per iteration of GradientDescent
		double baseTolerance;
		long int baseItr;
		
		NBall(long int d, long int n, long int m, double r = DEFAULT_RADIUS, double a = DEFAULT_ALPHA, double b = DEFAULT_BETA, double g = DEFAULT_GAMMA, double t = DEFAULT_TOLERANCE, long int i = DEFAULT_ITERATIONS);
		
};

#endif
