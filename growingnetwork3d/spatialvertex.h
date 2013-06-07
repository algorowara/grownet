#ifndef SPATIALVERTEX_H
#define SPATIALVERTEX_H

#include "../graph/vertex.h"

class SpatialVertex : public Vertex{
	
public:

	SpatialVertex(long int dimension, double* position, long int startTime);
	long int getDimension();
	double* getPosition();
	~SpatialVertex();

protected:

	long int dimension;
	double* position;
	
};

#endif
