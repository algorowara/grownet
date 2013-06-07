#include "spatialvertex.h"

SpatialVertex::SpatialVertex(long int dimension, double* position, long int startTime){
	
	this->dimension = dimension;
	this->position = position;
	this->startTime = startTime;
	
}

SpatialVertex::~SpatialVertex(){
	
	delete this->position;
	
}
