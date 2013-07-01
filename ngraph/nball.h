#ifndef NBALL_H
#define NBALL_H

#include "../graph/growingnetwork.h"
#include "../growingnetwork3d/spatialvertex.h"

#define NUM_NODES_ERR 13

#define NBALL_DEFAULT_RADIUS 1.0
#define NBALL_DEFAULT_GAMMA 0.25
#define NBALL_DEFAULT_ALPHA 1.0
#define NBALL_DEFAULT_TOLERANCE 0.01
#define NBALL_DEFAULT_ITERATIONS 100

class NBall : public GrowingNetwork {
	
	public:
	
		const long int DIM;	// the dimension of the NBall; immutable after initialization
		double radius;	// the radius of the NBall in N-dimensional space
		double alpha;	// node-node force constant
		double beta;	// attractive cloud force constant; always equal to the node-node force constant multiplied by the number of nodes
		double baseGamma;	// size of timestep per iteration of GradientDescent
		double baseTolerance;	// tolerance of the GradientDescent method
		long int baseItr;	// maximum number of iterations allowed for the GradientDescent method
		double iterationWeights;
		
		NBall(long int n, long int m, long int d, double r = NBALL_DEFAULT_RADIUS, double a = NBALL_DEFAULT_ALPHA, double g = NBALL_DEFAULT_GAMMA, double t = NBALL_DEFAULT_TOLERANCE, long int i = NBALL_DEFAULT_ITERATIONS);
		NBall(const NBall* obj);
		SpatialVertex* getNode(long int i) const;
		void grow(long int n);
		double* randomLocation();
		double linearDistance(Vertex* a, Vertex* b);
		double linearDistance(double* a, double* b);
		SpatialVertex** findMNearestNeighbors(SpatialVertex* start);
		double* sumForces(SpatialVertex* node);
		double calculatePotential();
		double calculateMinimumPotential();
		void equalize();
		void gradientDescent(double gamma, double tolerance, long int maxItr);
		~NBall();
		
};

#endif
