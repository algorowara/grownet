#ifndef NBALL_H
#define NBALL_H

#include "ngraph.h"

#define NUM_NODES_ERR 13

#define NBALL_DEFAULT_RADIUS 1.0
#define NBALL_DEFAULT_GAMMA 4.0
#define NBALL_DEFAULT_ALPHA 1.0
#define NBALL_DEFAULT_TOLERANCE 0.01
#define NBALL_DEFAULT_ITERATIONS 200
#define NBALL_DEFAULT_THRESHOLD 100
#define NBALL_DEFAULT_PERIOD 1

class NBall : public NGraph {
	
	public:
	
		float alpha;	// node-node force constant
		float beta;	// attractive cloud force constant; always equal to the node-node force constant multiplied by the number of nodes
		
		NBall(long int n, long int m, long int d, float r = NBALL_DEFAULT_RADIUS, float g = NBALL_DEFAULT_GAMMA, float t = NBALL_DEFAULT_TOLERANCE, long int i = NBALL_DEFAULT_ITERATIONS, long int et = NBALL_DEFAULT_THRESHOLD, long int ep = NBALL_DEFAULT_PERIOD);
		NBall(const NBall* obj);
		void grow(long int n);
		float* randomLocation();
		float linearDistance(Vertex* a, Vertex* b);
		float linearDistance(float* a, float* b);
		SpatialVertex** findMNearestNeighbors(SpatialVertex* start);
		float* sumForces(SpatialVertex* node);
		void equalize();
		void gradientDescent(float gamma, float tolerance, long int maxItr);
		static void exportObject(const NBall* nb, const char* filename);
		static NBall* importObject(const char* filename);
		~NBall();
		
};

#endif
