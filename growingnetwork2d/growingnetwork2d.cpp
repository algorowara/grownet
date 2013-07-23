#include "../growingnetwork2d/growingnetwork2d.h"
#include <vector>
#include <iostream>
#include <ctime>
#include <cmath>
#include <stdlib.h>

using namespace std;

/**
 * constructor for a GrowingNetwork2D
 * in which n nodes are placed on a circle
 * and each newly formed node (after the initial m nodes) create m new links
 */
GrowingNetwork2D::GrowingNetwork2D(long int n, long int m){
	
	static bool randSeeded = false;	// static variable to track this process's RNG seed
	
	if(!randSeeded){	// if the random number generator has not yet been seeded
		
		srand(std::time(NULL));	// do so using the current time as a seed
		randSeeded = true;	// make a note of this
		
	}

	if(n <= m){	// a growing network cannot be instantiated with as many or fewer nodes than its parameter m
	
		cerr<<"The number of nodes in a growing network must initially be greater than the parameter m.\n";
	
	}
	
	else{	// if there are at least m+1 nodes
		
		this->time = 0;	// set the time to zero
		this->m = m;	// store the value of the parameter m
		nodes.reserve(m+1);	// prepare the nodes vector to hold m+1 elements
		
		for(int i = 0; i < this->m+1; i++){	// insert the first m+1 nodes, which will form a clique
		
			Vertex* newNode = new Vertex(time);	// initialize a new Vertex at the current time
			insertNode(newNode, i);	// add a new node to the list
			tick();	// increment the current time
		
			for(int j = 0; j < i; j++){	// link all previously created nodes to this new node
			
				getNode(i)->addNeighbor(getNode(j));	// link node i with node j (0 <= j < i)
				
			}
			
		}
		
		grow(n - (m+1));	// grow the remaining n nodes, less the m+1 initial nodes
		
	}
	
}

/**
 * method to grow an arbitrary number of nodes using the following procedure:
 *   1) randomly select an inter-node interval
 *   2) place a new node in that interval, between the two pre-existing nodes on either side
 *   3) connect the new node to its m nearest neighbors
 */
void GrowingNetwork2D::grow(long int n){

	nodes.reserve(N + n);	// prepare the vector to receive n additional elements
	
	while(n > 0){	// while there remain nodes to be added
	
		long int pos = rand()%N;	// select a random interval between nodes
									// note that while this will never displace the 0th node,
									// it is still possible for new nodes to grow inside its interval
									// if the interval following the (N-1)th node is selected
	
		Vertex* newNode = new Vertex(time);	// create a new node
		tick();	// increment the current time
		n--;	// decrement the number of nodes not yet added
	
		for(int i = 0; i < m/2; i++){	// for the m/2 nodes on either side of the interval
			
			newNode->addNeighbor(getNode((pos-i)%N));	// connect the new node to its neighbor i nodes behind
			newNode->addNeighbor(getNode((pos+i+1)%N));	// and the other neighbor i+1 nodes ahead
		
		}
		
		if(m%2 > 0){	// if m is odd, randomly select the nearest left or right unconnected neighbor to be the mth connected neighbor
			
			long int oddIndex = (rand()%2 == 0) ? (pos - m/2) : (pos + m/2 + 1);	// randomly select an additional node beyond those connected
																					// where pos - m/2 and pos + m/2 + 1 (mod N) are the indices of the nearest left and right unconnected neighbors
			newNode->addNeighbor(getNode((pos+oddIndex)%N));	// add an edge between the nodes
			
		}
	
		insertNode(newNode, pos+1);	// place the new node in the middle of the interval, between pos and pos+1 (which is displaced, if it exists)
	
	}
	
	return;
	
}

/**
 * method to determine the arc distance of an edge between two given nodes
 * where the arc distance is defined as the number of intervals between adjacent nodes on the circle
 * between the two nodes that the edge connects
 */
long int GrowingNetwork2D::edgeArcDistance(Vertex* a, Vertex* b){
	
	// find the minimum of distance in ascending order vs. descending order of index
	long int oneway = abs(indexOf(b) - indexOf(a));
	long int another = N - abs(indexOf(b) - indexOf(a));
	return min(oneway, another);
	
}

/**
 * method to fill an array representing average arc distance (the element)
 * vs. the age of an edge (the index)
 * returns a dynamically allocated array of size N, as N-1 is the maximum edge age
 */
float* GrowingNetwork2D::edgeAgeVsArcDistance(){
	
	long int dist[N];	// temporary array to minimize the later effects of floating-point rounding
						// where the index is the age of the edge
						// and the distance is the sum of the distances
									
	float* ddist = new float[N];	// dynamically allocated array to store the age/arc distance relation
	
	for(int i = 0; i < N; i++){	// for each sum of distances
		
		dist[i] = 0;	// set it to zero to start with
		
	}
	
	for(int i = 0; i < N; i++){	// for each node
		
		for(int j = 0; j < K(i); j++){	// for each node's neighbors
			
			long int age = edgeAge(getNode(i), getNode(i)->getNeighbor(j));	// determine the age of the edge between them
			long int distance = edgeArcDistance(getNode(i), getNode(i)->getNeighbor(j));	// determine the distance covered by the edge between them
			dist[age] += distance;	// add the distance to the sum of distances for edges of a given age
			
		}
		
	}
	
	for(int i = 0; i < N; i++){	// for each possible age, 0:N-1 inclusive
	
		ddist[i] = dist[i]/(float)(2 * m);	// average the sum of float-counted distances over the m edges of the same age
		
	}
	
	return ddist;	// return that average
	
}

/**
 * method to calculate the linear distance between two Vertices, connected or not
 * where the network is a circle of circumference N
 */
float GrowingNetwork2D::linearDistance(Vertex* a, Vertex* b){

	float theta = 2 * M_PI * edgeArcDistance(a, b)/N;
	float r = N/(2 * M_PI);
	
	return 2 * r * sin(theta/2);
	
}

/**
 * method to return the proportion of edges of a given arc distance, for all possible arc distances
 * where the resulting array is N/2 + 1 long, as N/2 is the maximum possible arc distance
 */
float* GrowingNetwork2D::edgeArcDistanceDistribution(){
	
	long int num[N/2];	// local array relating the arc distance of an edge (the index)
						// to the number of edges with that arc distance (the content)
	float* prop = new float[N/2 + 1];	// dynamically allocated array not used until the end in order to minimize floating-point rounding error
	
	for(int i = 0; i < N/2 + 1; i++){
	
		num[i] = 0;
		
	}
	
	for(int i = 0; i < N; i++){	// for every node in the graph
		
		for(int j = 0; j < K(i); j++){	// for all of its neighbors
			
			long int arc = edgeArcDistance(getNode(i), getNode(i)->getNeighbor(j));	// determine the arc distance between them
			num[arc]++;	// increment the count of edges with the appropriate length
			
		}
		
	}
	
	for(int i = 0; i < N/2 + 1; i++){	// convert all count values into a proportion of the total number of edges
		
		prop[i] = ((float)num[i]) / (2 * ((m * (m-1))/2 + (N - m) * m));	// recall that for every node added (N), there are m edges, each of which will be float-counted
																			// except for the first m edges, which form a clique with (m * (m-1))/2 edges
		
	}
	
	return prop;
	
}
