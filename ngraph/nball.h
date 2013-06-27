#ifndef NBALL_H
#define NBALL_H

#include "../graph/growingnetwork.h"
#include "../growingnetwork3d/spatialvertex.h"

#define NUM_NODES_ERR 13

#define DEFAULT_RADIUS 1.0
#define DEFAULT_GAMMA 0.1
#define DEFAULT_ALPHA 0.1
#define DEFAULT_TOLERANCE 1.0
#define DEFAULT_ITERATIONS 1000

#define GUIDED_N 16

class NBall : public GrowingNetwork {
	
	public:
	
		const long int DIM;	// the dimension of the NBall; immutable after initialization
		double radius;	// the radius of the NBall in N-dimensional space
		double alpha;	// node-node force constant
		double beta;	// attractive cloud force constant; always equal to the node-node force constant multiplied by the number of nodes
		double baseGamma;	// size of timestep per iteration of GradientDescent
		double baseTolerance;	// tolerance of the GradientDescent method
		long int baseItr;	// maximum number of iterations allowed for the GradientDescent method
		
		NBall(long int n, long int m, long int d, double r = DEFAULT_RADIUS, double a = DEFAULT_ALPHA, double g = DEFAULT_GAMMA, double t = DEFAULT_TOLERANCE, long int i = DEFAULT_ITERATIONS);
		SpatialVertex* getNode(long int i);
		void grow(long int n);
		double* randomLocation();
		double linearDistance(Vertex* a, Vertex* b);
		SpatialVertex** findMNearestNeighbors(SpatialVertex* start);
		double* sumForces(SpatialVertex* node);
		double calculatePotential();
		double calculateMinimumPotential();
		void equalize();
		void gradientDescent(double gamma, double tolerance, long int maxItr);
		~NBall();
		
};

#endif
