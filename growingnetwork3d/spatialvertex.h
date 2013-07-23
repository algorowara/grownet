#ifndef SPATIALVERTEX_H
#define SPATIALVERTEX_H

#include "../graph/vertex.h"

class SpatialVertex : public Vertex{
	
public:

	long int dimension;
	float* position;

	SpatialVertex(long int dimension, float* position, long int startTime);
	SpatialVertex* getNeighbor(long int i) const;
	float radialDistance() const;
	~SpatialVertex();
	
};

#endif
