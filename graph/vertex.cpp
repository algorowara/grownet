#include "vertex.h"
#include <iostream>
#include <list>
#include <vector>
#include <climits>

using namespace std;

Vertex::Vertex(){

	distanceFromInitial = LONG_MAX;
	
}

Vertex::Vertex(long int time){
	
	distanceFromInitial = LONG_MAX;
	startTime = time;
	
}
	
void Vertex::addNeighbor(Vertex* neighbor){

	neighbors.push_back(neighbor);
	neighbor->neighbors.push_back(this);
	
}

double Vertex::clusteringCoefficient(){
	
	long int q = 0;
	long int k = neighbors.size();
	
	for(int a = 0; a < k; a++){	// for every neighboring node
		
		for(int b = 0; b < neighbors.at(a)->neighbors.size(); b++){	// examine each of that neighbor's neighbors
			
			for(int c = 0; c < k; c++){	// search through each neighboring node
				
					if(neighbors.at(a)->neighbors.at(b) == neighbors.at(c)){	// if the neighbor's neighbor is also a neighbor
					
						q++;	// increment the number of known neighbor-to-neighbor edges
						break;	// just in case a neighbor is included twice
						
					}
				
			}
			
		}
	
	}
	
	q /= 2;	// divide q by two, as it will have counted each edge twice
	
	return (2 * (double)q / ( k * (k-1) ) );	// return the number of neighbor-to-neighbor edges over the number of possible such edges
	
}

long int Vertex::getStartTime(){

	return startTime;
	
}

bool Vertex::hasNeighbor(Vertex* neighbor){
	
	for(int i = 0; i < neighbors.size(); i++){
		
		if(neighbors.at(i) == neighbor){
		
			return true;
			
		}
		
	}
	
	return false;
	
}
