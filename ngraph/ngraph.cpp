#include "ngraph.h"
#include "../growingnetwork3d/spatialvertex.h"
#include <fstream>
#include <cmath>

NGraph::NGraph(long int d) : DIM(d){
	
	
	
}

SpatialVertex* NGraph::getNode(const long int i) const{
	
	return (SpatialVertex*)nodes.at(i);
	
}

/**
 * method of calculating the linear distance between two nodes in a space
 * behavior is undefined and dangerous if both nodes are not SpatialVertices
 * if one node is of a higher dimension than the other, the lower dimension will be the space considered
 */
double NGraph::linearDistance(Vertex* a, Vertex* b){
	
	double sum = 0;
	
	for(long int i = 0; i < ((SpatialVertex*)a)->dimension && i < ((SpatialVertex*)b)->dimension; i++){
		
		sum += pow(((SpatialVertex*)a)->position[i] - ((SpatialVertex*)b)->position[i], 2.0);
		
	}
	
	return sqrt(sum);
	
}
