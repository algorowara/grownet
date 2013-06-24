#include "nball.h"

NBall::NBall(long int d, long int n, long int m, double r, double a, double g, double t, long int i) : DIM(d) {
	
	static bool randSeeded = false;	// static variable to check if the pseudorandom number generator has been seeded yet
	
	if(!randSeeded){	// if it has not
		
		srand(time(NULL));	// do so
		randSeeded = true;	// and reflect this by setting the variable to true
		
	}
	
	if(n < m+1){
		
		cerr<<"No growing network may be initialized with fewer than m+1 nodes; otherwise, no node may have degree m."<<endl;
		return NULL;
		
	}
	
	// initialize the non-constant variables
	this->m = m;
	this->radius = r;
	this->alpha = a;
	this->beta = a * N;
	this->baseGamma = g;
	this->baseTolerance = t;
	this->baseItr = i;
	
	for(long int i = 0; i < m+1; i++){	// add the first m+1 nodes, which will form a clique, each node having m links
		
		SpatialVertex* newNode = new SpatialVertex(DIM, randomLocation(), getTime());	// generate a new node at a random location
		
		for(long int j = 0; j < i; j++){	// for all nodes created prior to this one
			
			newNode->addNeighbor(getNode(j));	// link the new node and the old node
			
		}
		
		addNode(newNode);	// add the new node to the graph
		beta = alpha * N;	// update the attractive cloud force constant
		equalize();	// ensure that the nodes are spaced appropriately
		tick();	// move forward in time
		
	}
	
	grow(n - (m+1));	// grow the remaining nodes normally
	
}

/**
 * accessor method to retrieve a node, automatically casting it to a SpatialVertex*
 */
SpatialVertex* NBall::getNode(long int i){
	
	return nodes.at(i);
	
}

/**
 * method to grow n additional nodes in the current graph, one by one
 */
void NBall::grow(long int n){
	
	while(n > 0){	// for every one of n nodes that has been requested
		
		SpatialVertex* newNode = new SpatialVertex(DIM, randomLocation(), getTime());	// generate a new node at a random location
		SpatialVertex** nearNeighbors = findMNearestNeighbors(newNode);	// find its m nearest neighbors
		
		for(long int i = 0; i < m; i++){	// for each of those nearest neighbors
			
			newNode->addNeighbor(nearNeighbors[i]);	// link it to the new node
			
		}
		
		addNode(newNode);	// add the new node to the graph
		equalize();
		
		delete[] nearNeighbors;
	
		beta = alpha * N;	// update the attractive cloud force constant
		tick();
		n--;
		
	}
	
}

/**
 * method to generate a random position vector within the ball
 * generates a point from a uniform distribution over an n-cube
 * then rejects it and tries again if it lies outside the ball
 */
double* NBall::randomLocation(){
	
	double* pos = new double[DIM];
	double r_squared = 0;	// local sum to track the square of the radial length of the position vector
	
	for(long int i = 0; i < DIM; i++){	// for every dimension
		
		pos[i] = ((rand() / ((double)RAND_MAX)) - 0.5) * (2 * radius);	// generate a random number from 0 to RAND_MAX
																		// scale it to the range (0, 1) by dividing by RAND_MAX
																		// shift it to the range (-0.5, 0.5) by subtracting 0.5
																		// scale it to the range (-radius, radius) by multiplying by 2 * radius
		
		r_squared += pos[i] * pos[i];
		
	}
	
	if(r_squared > radius * radius){	// if the position lies outside the NBall
		
		delete pos;	// throw it out
		return randomLocation();	// try again
		
	}
	
	else{	// otherwise, return this position
		
		return pos;
		
	}
	
}

/**
 * method to determine the distance in Cartesian coordinates 
 * (rather than arc length) between two nodes
 * unsafe for use with classes without a position field
 */
double NBall::linearDistance(Vertex* a, Vertex* b){
	
	double d_squared = 0;
	
	for(long int i = 0; i < DIM; i++){
		
		d_squared += pow(((SpatialVertex*)a_loc)->position[i] - ((SpatialVertex*)b_loc)->position[i], 2);
		
	}
	
	return sqrt(d_squared);
	
}

/**
 * method to find the m nearest neighbors to a given SpatialVertex*,
 * where m is not given as a parameter but instead is taken from the field in NBall
 * allocates memory for an array of m pointers which should be subsequently deallocated
 */
SpatialVertex** NBall::findMNearestNeighbors(SpatialVertex* start){

	SpatialVertex** near = new SpatialVertex*[m];
	double dist[m];	//local record of the distance from the start of the m nearest neighbors

	for(int i = 0; i < m; i++){

		dist[i] = DBL_MAX;
	}

	for(int i = 0; i < N; i++){

		if(getNode(i) == start){	//do not consider the distance between the starting node and itself

			continue;

		}

		double d = linearDistance(start, getNode(i));

		for(int j = 0; j < m; j++){	//iterate through all distance squared records

			if(d < dist[j]){	//if this node is closer than the nearest neighbor j

				for (int k = m-1; k > j; k--){	//shift everything from j to m one position back

					dist[k] = dist[k-1];
					near[k] = near[k-1];

				}	
			
 				dist[j] = d;	// put the data for this new nearest neighbor 
				near[j] = getNode(i)	// in the spot previously occupied by neighbor j
				break;	// stop looking for others

			}

		}

	}

	return near;
	
}
