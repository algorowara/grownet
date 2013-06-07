#ifndef SPATIALVERTEX_H
#define SPATIALVERTEX_H

#include "../graph/vertex.h"

class SpatialVertex : public Vertex{
	
public:

	long int dimension;
	double* position;

	SpatialVertex(long int dimension, double* position, long int startTime);
	~SpatialVertex();
	
};

#endif
