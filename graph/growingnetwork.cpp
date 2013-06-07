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

/*
 * generate an array relating edge age (the index) to average edge connectivity (the content)
 * where edge connectivity is defined as the proportion of shortest paths of which the edge is a part
 * as there are as many ages as nodes, the resulting dynamically allocated array will be of size N
 */
double* GrowingNetwork::edgeAgeVsConnectivity(){

	long int connectivity[N];	// this array is local and temporary to minimize the effects of floating-point error
											// the index signifies an age
											// the element signifies the number of shortest paths which use that edge

	double* dconn = new double[N];
	
	for(int i = 0; i < N; i++){
		
		connectivity[i] = 0;
		
	}
	
	for(int i = 0; i < N; i++){	// for every node, construct the shortest paths from that node
	
		memoize(nodes.at(i));
		
		for(int j = 0; j < N; j++){	// for every node, use its shortest path from the root
			
			for(int k = 1; k < nodes.at(j)->pathFromInitial.size(); k++){	// for every node in the shortest path except the first
																			// note that starting at k = 1 ensures that node i is not included
																			
				Vertex* a = nodes.at(j)->pathFromInitial.at(k-1);	// use the previous node
				Vertex* b = nodes.at(j)->pathFromInitial.at(k);	// and the current node
				
				connectivity[edgeAge(a, b)]++;	// to determine the age of the edge, and increment the count accordingly
			
			}
			
		}
		
		clean(nodes.at(i));
		
	}
	
	for(int i = 0; i < N; i++){	// average and normalize each value
		
		dconn[i] = connectivity[i]/(double)m;	// average over the m edges of each age
		dconn[i] /= N * (N-1);	// normalize over the number of shortest paths
														// which will each be counted twice:
														// once from node A to node B, and once from node B to node A
		
	}
	
	return dconn;
	
}
