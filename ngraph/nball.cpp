#include "nball.h"
#include "../growingnetwork3d/spatialvertex.h"
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cfloat>
#include <cstring>
#include <iostream>

NBall::NBall(long int n, long int m, long int d, double r, double a, double g, double t, long int i, long int et, long int ep) : NGraph(d) {
	
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
	this->baseGam = g;
	this->baseTol = t;
	this->baseItr = i;
	this->equalizationThreshold = et;
	this->equalizationPeriod = ep;
	this->iterationWeights = 0;
	
	for(long int i = 0; i < m+1; i++){	// add the first m+1 nodes, which will form a clique, each node having m links
		
		SpatialVertex* newNode = new SpatialVertex(DIM, randomLocation(), getTime());	// generate a new node at a random location
		
		for(long int j = 0; j < i; j++){	// for all nodes created prior to this one
			
			newNode->addNeighbor(getNode(j));	// link the new node and the old node
			
		}
		
		addNode(newNode);	// add the new node to the graph
		beta = (alpha * N)/pow(radius, (double)DIM);	// update the attractive cloud force constant
		equalize();	
		tick();	// move forward in time
		
	}
	
	grow(n - (m+1));	// grow the remaining nodes normally
	
}

NBall::NBall(const NBall* obj) : NGraph(obj->DIM){
	
	this->m = obj->m;
	this->radius = obj->radius;
	this->alpha = obj->alpha;
	this->beta = obj->beta;
	this->baseGam = obj->baseGam;
	this->baseTol = obj->baseTol;
	this->baseItr = obj->baseItr;
	this->iterationWeights = obj->iterationWeights;
	this->equalizationThreshold = obj->equalizationThreshold;
	this->equalizationPeriod = obj->equalizationPeriod;
	
	for(long int i = 0; i < obj->N; i++){	// for all nodes in the original NBall
	
		double* copyPos = new double[DIM];	// create a duplicate position vector of the appropriate dimension
		memcpy(copyPos, obj->getNode(i)->position, DIM * sizeof(double));	// and populate it
		SpatialVertex* copyNode = new SpatialVertex(DIM, copyPos, obj->getNode(i)->getStartTime());	// create a new node with identical position and start time
		this->addNode(copyNode);	// add it to the duplicate NBall
		
	}
	
	for(long int i = 0; i < this->N; i++){	// for all nodes in the duplicate NBall
		
		SpatialVertex* duplicate = this->getNode(i);	// for this duplicate node
		SpatialVertex* original = obj->getNode(i);	// find the corresponding original node
		
		for(long int j = 0; j < original->neighbors.size(); j++){	// for all of the original node's neighbors
			
			long int index = obj->indexOf(original->getNeighbor(j));	// find their index in the original NBall
			duplicate->addNeighbor(this->getNode(index));	// connect this duplicate node to the duplicate node of the appropriate index
			
		}
		
	}
	
}

/**
 * accessor method to retrieve a node, automatically casting it to a SpatialVertex*
 */
SpatialVertex* NBall::getNode(long int i) const{
	
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
		
		if(N < this->equalizationThreshold || N%(this->equalizationPeriod) == 0){
			
			equalize();	// ensure that the nodes are spaced appropriately
		
		}
		
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
	
	return linearDistance(((SpatialVertex*)a)->position, ((SpatialVertex*)b)->position);
	
}

/**
 * method to determine the scalar distance between two position vectors
 * assumed to be of the same dimension as the NBall; unsafe/undefined otherwise
 */
double NBall::linearDistance(double* a, double* b){
	
	double d_squared = 0;
	
	for(long int i = 0; i < DIM; i++){
		
		d_squared += pow(a[i] - b[i], 2.0);
		
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
		
		if(dist == 0){	// if these nodes are right on top of each other
			
			continue;	// skip to avoid the divide-by-zero error
			
		}
		
		mag = alpha / pow(dist, (double)DIM-1);	// and the magnitude of the electrostatic repulsive force
		
		for(long int j = 0; j < DIM; j++){	// for every dimension
			
			force[j] += mag * (node->position[j] - other->position[j])/dist;	// multiply the magnitude by the appropriate component of the unit displacement vector
																				// where each component is the difference in position normalized by the total distance
			
		}
		
	}
	
	dist = node->radialDistance();	// calculate the attractive cloud force
	mag = -beta * dist;	// since the force is attractive, it will be negative with respect to the radially outwards vector
	
	if(dist > 0){	// if this node is not at the origin, in which case finding the unit vector will result in a nan

		for(long int i = 0; i < DIM; i++){
			
			force[i] += mag * ((node->position[i]) / dist);	// here, the unit displacement vector is the unit position vector
			
		}
		
	}
	
	return force;
	
}

/**
 * method which properly tunes the relevant parameters before passing them to gradientDescent
 * gamma scales with 1/(N^(1 + 1/D)) because the number of nodes, and therefore the force, will increase as N,
 * and the average separation between nodes will decrease with N^(1/D), as the D-dimensional volume per unit decreases with N
 */
void NBall::equalize(){
	
	gradientDescent(baseGam /pow(N, 1 + 1.0/DIM), baseTol, baseItr);
	
}

void NBall::gradientDescent(const double gamma, const double tolerance, long int maxItr){
	
	double* netForce[N];	// a record of the net force on each node
	double prevMaxDisp = DBL_MAX;	// record of the maximum scalar displacement of a node on the previous iteration
									// initially set to an impossibly large value to ensure at least one iteration occurs
	double tolMaxDisp = (radius * pow(N, -1.0/DIM) * tolerance) * baseGam;	// the maximum tolerated displacement
																			// if the maximum displacement is above this value, continue iterating
																			// scales with gamma to provide equal treatment across timestep sizes

	memset(netForce, 0, N * sizeof(double*));
	
	while(prevMaxDisp > tolMaxDisp && maxItr > 0){	// while there remains some excess energy above tolerance
													// and the hard limit of iterations has not been passed

		double maxDisp = 0;

		#pragma omp parallel shared(netForce, maxDisp)
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
			
				double oldPos[DIM];
				
				for(long int j = 0; j < DIM; j++){	// for every dimension
				
					oldPos[j] = getNode(i)->position[j];	// store the old position
					getNode(i)->position[j] += gamma * netForce[i][j];	// displace the node in that dimension
																		// according to the force on it and the timestep size
					
				}
				
				double disp = linearDistance(oldPos, getNode(i)->position);	// calculate the displacement
																				// between the old and current positions
				
				#pragma omp critial (maximumDisplacement)	// encase this comparison in a critical region
				{
					
					if(disp > maxDisp){	// if the displacement of this particular node is larger than the largest seen this iteration
						
						maxDisp = disp;	// declare this displacement the largest of the iteration
						
					}
					
				}
				
				delete[] netForce[i];	// deallocate the now-outdated force vector
				
			}
			
		}
		
		prevMaxDisp = maxDisp;	// record the largest displacement encountered here	
		maxItr--;	// move one iteration closer to stopping regardless of potential
		
	}
	
	this->iterationWeights += (baseItr - maxItr) * N * N;
	
	return;	// return statement placed here for clarity only
	
}

NBall::~NBall(){
	
}
