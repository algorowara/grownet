#include "growingnetwork.h"
#include <iostream>
#include <cstring>

/*
 * as the time field of a GrowingNetwork is important and not to be tampered with,
 * there is a getter method here to return a copy of it
 */
long int GrowingNetwork::getTime() const{

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
 * where the age of each node is defined as the difference between the current time and its start time minus one
 * as the current time has not "happened" (had a node added at that timestep) yet
 * and the age of the edge is the minimum of the two nodes' ages, as all new edges are brought into existence by one of their nodes
 */
long int GrowingNetwork::edgeAge(Vertex* a, Vertex* b){
	
	if(!(a->hasNeighbor(b))){
	
		return -1;
		
	}
	
	return min(getTime() - a->getStartTime() - 1, getTime() - b->getStartTime() - 1);
	
}

/**
 * generate an array relating edge age (the index) to average edge betweenness (the content)
 * where edge betweenness is defined as the proportion of shortest paths of which the edge is a part
 * as there are as many ages as nodes, the resulting dynamically allocated array will be of size N
 */
float* GrowingNetwork::edgeAgeVsBetweenness(){

	long int betweenness[N];	// this array is local and temporary to minimize the effects of floating-point error
								// the index signifies an age
								// the element signifies the number of shortest paths which use that edge
	long int num_edges[N];	// local array to hold the number of edges (the element) of each age (the index)
	float* dbet = new float[N];
	memset(betweenness, 0, N * sizeof(long int));
	memset(&num_edges, 0, N * sizeof(long int));
	
	for(long int i = 0; i < N; i++){	// populate the num_edges array
		
		for(long int j = 0; j < N; j++){	// for each pair of nodes
			
			if(edgeAge(getNode(i), getNode(j)) >= 0){	// if an edge exists
				
				num_edges[edgeAge(getNode(i), getNode(j))]++;	// record it, which will float-count each edge: once as i->j and once as j->i
				
			}
			
		}
		
	}
	
	#pragma omp parallel shared(betweenness)
	{
	
		#pragma omp for schedule(guided)
		for(long int i = 0; i < N; i++){	// for every node, construct the shortest paths from that node
		
			Graph* dup = new Graph();	// create a duplicate network which can be memoized
			
			for(long int j = 0; j < N; j++){	// with N duplicate nodes
				
				dup->addNode(new Vertex(this->getNode(j)->getStartTime()));	// with the same start time as the original
				
			}
			
			for(long int j = 0; j < N; j++){	// for each duplicate node
				
				for(long int k = 0; k < K(j); k++){	// for each neighbor
					
					long int index = this->indexOf(this->getNode(j)->getNeighbor(k));	// find the index of that neighbor
					dup->getNode(j)->addNeighbor(dup->getNode(index));	// link the appropriate duplicates
					
				}
				
			}
		
			memoize(dup->getNode(i));	// with the duplicate constructed, memoize its nodes
			
			for(long int j = 0; j < N; j++){	// for every node, use its shortest path from the root
				
				for(int k = 1; k < dup->getNode(j)->pathFromInitial.size(); k++){	// for every node in the shortest path except the first
																					// note that starting at k = 1 ensures that node i is not included
																				
					long int a = dup->indexOf(dup->getNode(j)->pathFromInitial.at(k-1));	// use the previous node's index
					long int b = dup->indexOf(dup->getNode(j)->pathFromInitial.at(k));	// and the current node's index
					
					#pragma omp atomic
					betweenness[this->edgeAge(this->getNode(a), this->getNode(b))]++;	// to determine the age of the edge, and increment the count accordingly
																						// because the duplicate cannot be made a GrowingNetwork (which is a virtual class)
																						// the ages must be referenced from this GrowingNetwork
				}
				
			}
			
			delete dup;	// destroy the duplicate GrowingNetwork once the need for it is over
			
		}
		
	}
	
	for(int i = 0; i < N; i++){	// average and normalize each value
		
		if(num_edges[i] != 0){	// if there exist edges of this age
		
			dbet[i] = betweenness[i] / num_edges[i];	// average over the (float-counted) number of edges
			
		}
		
		else{	// if no edges of this age have been observed
			
			dbet[i] = 0;	// wipe the data for this age
			
		}
		
		dbet[i] /= (N * (N-1));	// normalize over the number of shortest paths
								// which will each be counted twice:
								// once from node A to node B, and once from node B to node A
	}
	
	return dbet;
	
}

/**
 * returns a dynamically allocated array of floats relating
 * edge age (the indices of the array) to
 * average linear distance for that edge (the elements of the array)
 */
float* GrowingNetwork::edgeAgeVsLinearDistance(){
	
	float* lin = new float[N];
	long int num_edges[N];
	memset(lin, 0, sizeof(float) * N);
	memset(num_edges, 0, N * sizeof(long int));
	
	for(int i = 0; i < N; i++){	// iterating over all nodes in the GrowingNetwork
		
		for(int j = 0; j < K(i); j++){	// for each neighbor (and the corresponding edge)
		
			Vertex* a = getNode(i);
			Vertex* b = getNode(i)->getNeighbor(j);
			lin[edgeAge(a, b)] += linearDistance(a, b);	// add the length of the edge to the sum
			num_edges[edgeAge(a, b)]++;	// increment the count of edges of this age seen
			
		}
		
	}
	
	for(int i = 0; i < N; i++){	// for each sum of lengths
		
		if(num_edges[i] != 0){	// if some edges were counted for this age
			
			lin[i] /= num_edges[i];	// average the length over all the times an edge of this age was counted
									// which will naturally include float-counting, once for each node
			
		}
		
		else{	// otherwise, if no edges of this age were seen
			
			lin[i] = 0;	// make sure there is no misleading information
			
		}
		
	}
	
	return lin;
	
}

/**
 * method to output a dynamically allocated array of
 * node age (the index and independent variable)
 * vs node degree (the content at index and dependent variable)
 */
float* GrowingNetwork::nodeAgeVsDegree(){
	
	float* degree = new float[N];	// allocate a new array with as many indices as there are ages
										// since each node has a different age, there are N different ages
	
	for(long int i = 0; i < N; i++){	// for each node
		
		long int age = (this->getTime()-1) - (getNode(i)->getStartTime());	// calculate the time since its creation
																			// noting that the current time has not "happened" yet
		
		degree[age] = K(i);	// note the degree of that node, and store it at the appropriate index
		
	}
	
	return degree;
	
}
