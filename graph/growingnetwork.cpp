#include "growingnetwork.h"

/*
 * as the time field of a GrowingNetwork is important and not to be tampered with,
 * there is a getter method here to return a copy of it
 */
long int GrowingNetwork::getTime(){

	return time;

}

/*
 * to better keep track of the time field, there is a method to increment time
 */
void GrowingNetwork::tick(){

	time++;

}

/**
 * method to determine the age of an edge between two nodes
 * where the age of each node is defined as the difference between the current time and its start time
 * and the age of the edge is the minimum of the two nodes' ages, as all new edges are brought into existence by one of their nodes
 */
long int GrowingNetwork::edgeAge(Vertex* a, Vertex* b){
	
	if(!(a->hasNeighbor(b))){
	
		return -1;
		
	}
	
	return min(getTime() - a->getStartTime(), getTime() - b->getStartTime());
	
}

/**
 * generate an array relating edge age (the index) to average edge betweenness (the content)
 * where edge betweenness is defined as the proportion of shortest paths of which the edge is a part
 * as there are as many ages as nodes, the resulting dynamically allocated array will be of size N
 */
double* GrowingNetwork::edgeAgeVsBetweenness(){

	long int betweenness[N];	// this array is local and temporary to minimize the effects of floating-point error
								// the index signifies an age
								// the element signifies the number of shortest paths which use that edge

	double* dbet = new double[N];
	
	for(int i = 0; i < N; i++){
		
		betweenness[i] = 0;
		
	}
	
	for(int i = 0; i < N; i++){	// for every node, construct the shortest paths from that node
	
		memoize(nodes.at(i));
		
		for(int j = 0; j < N; j++){	// for every node, use its shortest path from the root
			
			for(int k = 1; k < nodes.at(j)->pathFromInitial.size(); k++){	// for every node in the shortest path except the first
																			// note that starting at k = 1 ensures that node i is not included
																			
				Vertex* a = nodes.at(j)->pathFromInitial.at(k-1);	// use the previous node
				Vertex* b = nodes.at(j)->pathFromInitial.at(k);	// and the current node
				
				betweenness[edgeAge(a, b)]++;	// to determine the age of the edge, and increment the count accordingly
			
			}
			
		}
		
		clean(nodes.at(i));
		
	}
	
	for(int i = 0; i < N; i++){	// average and normalize each value
		
		dbet[i] = betweenness[i]/(double)m;	// average over the m edges of each age
		dbet[i] /= N * (N-1);	// normalize over the number of shortest paths
								// which will each be counted twice:
								// once from node A to node B, and once from node B to node A
		
	}
	
	return dbet;
	
}

/**
 * returns a dynamically allocated array of doubles relating
 * edge age (the indices of the array) to
 * average linear distance for that edge (the elements of the array)
 */
double* GrowingNetwork::edgeAgeVsLinearDistance(){
	
	double* lin = new double[N];
	
	for(int i = 0; i < N; i++){	// iterating over all nodes in the GrowingNetwork
		
		for(int j = 0; j < K(i); i++){	// for each neighbor (and the corresponding edge)
		
			Vertex* a = nodes.at(i);
			Vertex* b = nodes.at(i)->neighbors.at(j);
			lin[edgeAge(a, b)] += edgeLinearDistance(a, b);	// add the length of the edge to the sum
			
		}
		
	}
	
	for(int i = 0; i < N; i++){	// for each sum of lengths
		
		lin[i] /= 2 * m;	// average over the number of edges of age i (m edges)
							// and the number of time each edge is counted (twice)
		
	}
	
	return lin;
	
}
