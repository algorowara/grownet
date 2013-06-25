#include "nball.h"
#include "../growingnetwork3d/spatialvertex.h"
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cfloat>
#include <cstring>
#include <iostream>

NBall::NBall(long int n, long int m, long int d, double r, double a, double g, double t, long int i) : DIM(d) {
	
	static bool randSeeded = false;	// static variable to check if the pseudorandom number generator has been seeded yet
	
	if(!randSeeded){	// if it has not
		
		srand(std::time(NULL));	// do so
		randSeeded = true;	// and reflect this by setting the variable to true
		
	}
	
	if(n < m+1){
		
		cerr<<"No growing network may be initialized with fewer than m+1 nodes; otherwise, no node could have degree m."<<endl;
		throw NUM_NODES_ERR;
		
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
		beta = (alpha * N)/pow(radius, (double)DIM);	// update the attractive cloud force constant
		equalize();	// ensure that the nodes are spaced appropriately
		tick();	// move forward in time
		
	}
	
	grow(n - (m+1));	// grow the remaining nodes normally
	
}

/**
 * accessor method to retrieve a node, automatically casting it to a SpatialVertex*
 */
SpatialVertex* NBall::getNode(long int i){
	
	return (SpatialVertex*)nodes.at(i);
	
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
	
		beta = (alpha * N)/pow(radius, (double)DIM);	// update the attractive cloud force constant
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
	
		delete[] pos;	// throw it out
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
		
		d_squared += pow(((SpatialVertex*)a)->position[i] - ((SpatialVertex*)b)->position[i], 2.0);
		
	}
	
	return sqrt(d_squared);
	
}

/**
 * method to determine the distance from the center to a given node
 * relies on the linearDistance method
 */
double NBall::radialDistance(SpatialVertex* node){
	
	static double* centerPosition = new double[DIM];	
	static SpatialVertex* centralNode = new SpatialVertex(DIM, centerPosition, -1);
	
	return linearDistance(node, centralNode);
	
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
				near[j] = getNode(i);	// in the spot previously occupied by neighbor j
				break;	// stop looking for others

			}

		}

	}

	return near;
	
}

/**
 * method to calculate the electrostatic force on a given node (1/r^(DIM-1))
 * using the inter-node repulsive forces and the attractive cloud force
 * allocates and returns an array which should be subsequently deallocated
 */
double* NBall::sumForces(SpatialVertex* node){
	
	double* force = new double[DIM];
	SpatialVertex* other;	// placeholder for a pointer to some other node in two-body interactions
	double mag, dist;	// local fields to hold the force magnitude and distance of a two-body interaction
	
	memset(force, 0, DIM * sizeof(double));	// ensure that the components of force are set to zero
	
	for(long int i = 0; i < N; i++){	// for every node in the graph
		
		other = getNode(i);
		
		if(other == node){	// except this one
			
			continue;	// the self-force is not considered here
			
		}
		
		dist = linearDistance(node, other);	// remember the distance
		mag = alpha / pow(dist, (double)DIM-1);	// and the magnitude of the electrostatic repulsive force
		
		for(long int j = 0; j < DIM; j++){	// for every dimension
			
			force[j] += mag * (node->position[j] - other->position[j])/dist;	// multiply the magnitude by the appropriate component of the unit displacement vector
																				// where each component is the difference in position normalized by the total distance
			
		}
		
	}
	
	dist = radialDistance(node);	// calculate the attractive cloud force
	mag = -beta * dist;	// since the force is attractive, it will be negative with respect to the radially outwards vector
	
	for(long int i = 0; i < DIM; i++){
		
		force[i] += mag * ((node->position[i]) / dist);	// here, the unit displacement vector is the unit position vector
		
	}
	
	return force;
	
}

/**
 * method to calculate the potential of the NBall
 * with the convention that potential from repulsive forces is negative
 * and potential from the attractive cloud is positive
 */
double NBall::calculatePotential(){
	
	double pot = 0;
	
	#pragma omp parallel shared(pot)
	{
		
		double localSum = 0;
		
		#pragma omp for schedule(guided)
		for(long int i = 0; i < N; i++){	// for all nodes
			
			for(long int j = i+1; j < N; j++){	// for all two-body interactions and all N * (N-1) / 2 possible order-independent pairs
				
				if(DIM != 2){	// if the dimension is not equal to two, integrating the electrostatic force over distance is simple
					
					localSum += -1 * alpha / pow(linearDistance(getNode(i), getNode(j)), (double)DIM-2);
					
				}
				
				else if(DIM == 2){	// if the dimension is equal to two, the potential is logarithmic, because the force goes as 1/r
					
					localSum += -1 * alpha * log(linearDistance(getNode(i), getNode(j)));
					
				}
				
			}
			
			localSum += (beta / 2.0) * pow(radialDistance(getNode(i)), 2.0);	// calculate the potential due to the attractive cloud
			
		}
		
		#pragma omp atomic
		pot += localSum;
		
	}
	
	return pot;
	
}

/**
 * method which properly tunes the relevant parameters before passing them to gradientDescent
 * gamma scales with 1/(N^(1 + 1/D)) because the number of nodes, and therefore the force, will increase as N,
 * and the average separation between nodes will decrease with N^(1/D), as the D-dimensional volume per unit decreases with N
 */
void NBall::equalize(){
	
	gradientDescent(baseGamma / (N * pow(N, 1.0/DIM)), baseTolerance, baseItr);
	
}

void NBall::gradientDescent(double gamma, double tolerance, long int maxItr){
	
	double* netForce[N];	// a record of the net force on each node
	double previousPotential = DBL_MAX;	// record of the previous potential; set to an impossibly large value to ensure that at least one iteration occurs
	double toleratedPotential;	// if the potential is above this value, continue iterating
	double temp_minimum_potential = DBL_MAX;
	
	memset(netForce, 0, N * sizeof(double*));
	
	if(N > GUIDED_N){	// if the graph size is large enough that an estimate of the minimum possible potential is reasonably accurate
		
		toleratedPotential = -DBL_MAX;	/* TODO: find a general solution for the minimum possible potential */
		
	}
	
	else{	// otherwise, if there are few enough nodes that the minimum potential is uncertain
		
		toleratedPotential = -DBL_MAX;	// go through the maximum number of iterations anyway
		
	}
	
	while(previousPotential > toleratedPotential && maxItr > 0){	// while there remains some excess energy above tolerance
																	// and the hard limit of iterations has not been passed

		#pragma omp parallel shared(netForce)
		{
		
			#pragma omp for schedule(guided)
			for(long int i = 0; i < N; i++){	// for each node
				
				netForce[i] = this->sumForces(getNode(i));	// store the current net force on each node
				
			}
			
			
		}
		
		#pragma omp parallel shared(netForce)
		{
			
			#pragma omp for schedule(guided)
			for(long int i = 0; i < N; i++){	// for every node
				
				for(long int j = 0; j < DIM; j++){	// for every dimension
					
					getNode(i)->position[j] += gamma * netForce[i][j];	// displace the node in that dimension
																		// according to the force on it and the timestep size
					
				}
				
				delete[] netForce[i];	// deallocate the now-outdated force vector
				
			}
			
		}
		
		previousPotential = calculatePotential();	// update the record of the potential
		
		if(previousPotential < temp_minimum_potential){
			
			temp_minimum_potential = previousPotential;
			
		}
		
		maxItr--;	// move one iteration closer to stopping regardless of potential
		
	}
	
	//cout<<N<<" "<<temp_minimum_potential<<endl;
	
	return;	// return statement placed here for clarity only
	
}

NBall::~NBall(){
	
}
