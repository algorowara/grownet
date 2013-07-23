#include "spatialvertex.h"
#include <cmath>

/**
 * constructor which specifies
 * dimension (the size of the position vector)
 * the position vector itself (allocated elsewhere)
 * the creation time
 */
SpatialVertex::SpatialVertex(long int dimension, float* position, long int startTime){
	
	this->dimension = dimension;
	this->position = position;
	this->startTime = startTime;
	
}

/**
 * method to retrieve a neighbor at a given index
 */
SpatialVertex* SpatialVertex::getNeighbor(long int i) const{
	
	return (SpatialVertex*)(neighbors.at(i));
	
}

/**
 * method to determine the scalar distance between this position and the origin (0, 0...0)
 */
float SpatialVertex::radialDistance() const{
	
	float rsq = 0;
	
	for(long int i = 0; i < dimension; i++){
		
		rsq += position[i] * position[i];
		
	}
	
	return sqrt(rsq);
	
}

/**
 * destructor which deallocates the position vector, which is this object's responsibility
 */
SpatialVertex::~SpatialVertex(){
	
	delete[] this->position;
	
}
