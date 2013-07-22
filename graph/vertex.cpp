#include "vertex.h"
#include <iostream>
#include <list>
#include <vector>
#include <climits>

using namespace std;

/**
 * default constructor
 * does not initialize the time field
 * sets distanceFromInitial to LONG_MAX, meaning "not yet memoized"
 */
Vertex::Vertex(){

	distanceFromInitial = LONG_MAX;
	
}

/**
 * constructor for a Vertex in a GrowingNetwork
 * sets the start time, usually to the GrowingNetwork's current time
 * also sets distanceFromInitial to the signal value of "not yet memoized"
 */
Vertex::Vertex(long int time){
	
	distanceFromInitial = LONG_MAX;
	startTime = time;
	
}

/**
 * adds a neighbor to this Vertex by placing it in its list of neighbors
 * and placing this Vertex in the neighbor's list of neighbors
 * if and only if
 *   a) the two Vertices are not already neighbors and
 *   b) the two Vertices are different (no edges to the self allowed)
 */
void Vertex::addNeighbor(Vertex* neighbor){
	
	if(!(this->hasNeighbor(neighbor)) && this != neighbor){
		
		this->neighbors.push_back(neighbor);
		neighbor->neighbors.push_back(this);
		
	}
	
}

/**
 * method to calculate the clustering coefficient of this Vertex
 * defined as the number of edges between two neighbors of this Vertex
 * normalized over the number of possible edges between such neighbors
 */
float Vertex::clusteringCoefficient(){
	
	long int q = 0;
	long int k = neighbors.size();
	
	for(int a = 0; a < k; a++){	// for every neighboring node
		
		for(int b = 0; b < getNeighbor(a)->neighbors.size(); b++){	// examine each of that neighbor's neighbors
			
			for(int c = 0; c < k; c++){	// search through each neighboring node
				
					if(getNeighbor(a)->getNeighbor(b) == getNeighbor(c)){	// if the neighbor's neighbor is also a neighbor
					
						q++;	// increment the number of known neighbor-to-neighbor edges
						
					}
				
			}
			
		}
	
	}
	
	q /= 2;	// divide q by two, as it will have counted each edge twice
	
	return ((float)q / ( ( k * (k-1) ) / 2) );	// return the number of neighbor-to-neighbor edges over the number of possible such edges
	
}

/**
 * accessor method to get the starting time
 * used to emphasize that while Vertex::startTime is a public field,
 * it should not be written to by other entities
 */
long int Vertex::getStartTime(){

	return startTime;
	
}

/**
 * method to check if this Vertex and another Vertex are neighbors
 * checks this Vertex's neighbor vector only
 */
bool Vertex::hasNeighbor(Vertex* neighbor){
	
	for(int i = 0; i < neighbors.size(); i++){
		
		if(getNeighbor(i) == neighbor){
		
			return true;
			
		}
		
	}
	
	return false;
	
}

/**
 * accesses Vertex::neighbors to return a given neighbor
 * without requiring outside access to Vertex::neighbors
 */
Vertex* Vertex::getNeighbor(long int i){
	
	return neighbors.at(i);
	
}
