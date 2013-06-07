#include "../graph/vertex.h"
#include "../graph/growingnetwork.h"
#include "growingnetwork2d.h"
#include <vector>
#include <iostream>
#include <ctime>
#include <cmath>
#include <stdlib.h>

using namespace std;

GrowingNetwork2D::GrowingNetwork2D(long int n, long int m){
	
	srand(std::time(NULL));

	if(n <= m){	// a growing network cannot be instantiated with as many or fewer nodes than its parameter m
	
		cerr<<"The number of nodes in a growing network must initially be greater than the parameter m.\n";
	
	}
	
	else if(m%2 != 0){	// a two-dimensional growing network must have an even parameter m
		
	
		cerr<<"The parameter m must be even.\n";
		
	}
	
	else{
		
		this->time = 0;	// set the time to zero
		nodes.reserve(m);	// prepare the vector containing nodes to contain the first m elements
		this->m = m;	// remember the value of the parameter m
	
		for(int i = 0; i < this->m+1; i++){	// insert the first m+1 nodes, which will form a clique
		
			Vertex* newNode = new Vertex(time);
			insertNode(newNode, i);	// add a new node to the list
			tick();	// increment the current time
			n--;	// note that one less node needs to be added
		
			for(int j = 0; j < i; j++){	// link all previously created nodes to this new node
			
				nodes.at(i)->addNeighbor(nodes.at(j));	// link node i with node j (0 <= j < i)
				
			}
			
		}
		
		grow(n);
		
	}
	
}

void GrowingNetwork2D::grow(long int n){

		nodes.reserve(N + n);	// prepare the vector to receive n additional elements
	
	while(n > 0){	// while there remain nodes to be added
	
		int pos = rand()%N;	// select a random interval between nodes
										// note that while this will never displace the zeroth node,
										// it is still possible for new nodes to grow inside its interval
	
		Vertex* newNode = new Vertex(time);	// create a new node
		tick();	// increment the current time
		n--;	// decrement the number of nodes not yet added
	
		for(int i = 0; i < m/2; i++){	// for the m/2 nodes on either side of the interval
			
			newNode->addNeighbor(nodes.at(pos-i));	// connect the new node to its neighbor i nodes behind
			newNode->addNeighbor(nodes.at((pos+i+1)%N));	// and the other neighbor i+1 nodes ahead
		
		}
	
		insertNode(newNode, pos+1);	// place the new node in the middle of the interval, between pos and pos+1 (which is displaced)
	
	}
	
}

/*
 * method to determine the arc distance of an edge between two given nodes
 * where the arc distance is defined as the number of intervals between adjacent nodes on the circle
 * between the two nodes that the edge connects
 */
long int GrowingNetwork2D::edgeArcDistance(Vertex* a, Vertex* b){
	
	// find the minimum of distance in ascending order vs. descending order of index
	return min(abs(indexOf(b) - indexOf(a)), N - abs(indexOf(b) - indexOf(a)));
	
}

double* GrowingNetwork2D::edgeAgeVsArcDistance(){
	
	long int dist[N];	// temporary array to minimize the later effects of floating-point rounding
									// where the index is the age of the edge
									// and the distance is the sum of the distances
									
	double* ddist = new double[N];
	
	for(int i = 0; i < N; i++){
		
		dist[i] = 0;
		
	}
	
	for(int i = 0; i < N; i++){
		
		for(int j = 0; j < nodes.at(i)->neighbors.size(); j++){
			
			long int age = edgeAge(nodes.at(i), nodes.at(i)->neighbors.at(j));
			long int distance = edgeArcDistance(nodes.at(i), nodes.at(i)->neighbors.at(j));
			dist[age] += distance;
			
		}
		
	}
	
	for(int i = 0; i < N; i++){
	
		ddist[i] = dist[i]/(double)(2 * m);	// average the sum of double-counted distances over the m edges of the same age
		
	}
	
	return ddist;
	
}

double GrowingNetwork2D::edgeLinearDistance(Vertex* a, Vertex* b){

	double theta = 2 * M_PI * edgeArcDistance(a, b)/N;
	double r = N/(2 * M_PI);
	
	return 2 * r * sin(theta/2);
	
}

/*
 * method to return the proportion of edges of a given arc distance, for all possible arc distances
 * where the resulting array is N/2 + 1 long, as N/2 is the maximum possible arc distance
 */
double* GrowingNetwork2D::edgeArcDistanceDistribution(){
	
	long int num[N/2];	// local array relating the arc distance of an edge (the index)
						// to the number of edges with that arc distance (the content)
	double* prop = new double[N/2 + 1];	// dynamically allocated array not used until the end in order to minimize floating-point rounding error
	
	for(int i = 0; i < N/2 + 1; i++){
	
		num[i] = 0;
		
	}
	
	for(int i = 0; i < N; i++){	// for every node in the graph
		
		for(int j = 0; j < nodes.at(i)->neighbors.size(); j++){	// for all of its neighbors
			
			long int arc = edgeArcDistance(nodes.at(i), nodes.at(i)->neighbors.at(j));	// determine the arc distance between them
			num[arc]++;	// increment the count of edges with the appropriate length
			
		}
		
	}
	
	for(int i = 0; i < N/2 + 1; i++){	// convert all count values into a proportion of the total number of edges
		
		prop[i] = ((double)num[i]) / (2 * ((m * (m-1))/2 + (N - m) * m));	// recall that for every node added (N), there are m edges, each of which will be double-counted
																			// except for the first m edges, which form a clique with (m * (m-1))/2 edges
		
	}
	
	return prop;
	
}
