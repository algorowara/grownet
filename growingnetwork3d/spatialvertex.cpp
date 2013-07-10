#include "spatialvertex.h"
#include <cmath>

SpatialVertex::SpatialVertex(long int dimension, float* position, long int startTime){
	
	this->dimension = dimension;
	this->position = position;
	this->startTime = startTime;
	
}

SpatialVertex* SpatialVertex::getNeighbor(long int i){
	
	return (SpatialVertex*)(neighbors.at(i));
	
}

float SpatialVertex::radialDistance(){
	
	float rsq = 0;
	
	for(long int i = 0; i < dimension; i++){
		
		rsq += pow(position[i], 2.0);
		
	}
	
	return sqrt(rsq);
	
}

SpatialVertex::~SpatialVertex(){
	
	delete[] this->position;
	
}
