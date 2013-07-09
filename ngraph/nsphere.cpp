#include "nsphere.h"
#include <cstring>

using namespace std;

NSphere::NSphere(long int n, long int m, long int d, double baseGam, double baseTol, long int baseItr, long int threshold, long int period) : DIM(d){
	
	static bool randSeeded = false;
	
	if(!randSeeded){	// if the random number generator has not yet been seeded
		
		srand(std::time(NULL));	// do so with the current time as the seed
		randSeeded = true;	// make a note for future initializations
		
	}
	
	this->time = 0;
	this->radius = NSPHERE_DEFAULT_RADIUS;
	this->m = m;
	this->baseGam = baseGam;
	this->baseTol = baseTol;
	this->baseItr = baseItr;
	this->equalizationThreshold = threshold;
	this->equalizationPeriod = period;
	this->iterationWeights = 0;
	
	for(long int i = 0; i < m+1; i++){	// for the first m+1 nodes, which form a clique
	
		SpatialVertex* newNode = new SpatialVertex(DIM+1, randomLocation(), getTime());	// generate a new node in DIM dimensions
		addNode(newNode);
		
		for(long int j = 0; j < nodes.size()-1; j++){	// link it to all previously created nodes
			
			newNode->addNeighbor(getNode(j));
			
		}
		
		
		equalize();	// equalize the distribution of nodes
		tick();
		n--;
		
	}
	
	grow(n);	// grow the remaining nodes normally
	
}

/**
 * method to retrieve a node and cast it to SpatialVertex*
 */
SpatialVertex* NSphere::getNode(long int i){
	
	return (SpatialVertex*)(nodes.at(i));
	
}

/**
 * method to randomly add n new nodes to the NSphere, one by one
 * equalizing the system after each addition
 */
void NSphere::grow(long int n){
	
	while(n > 0){	// while there remain nodes to add
		
		SpatialVertex* newNode = new SpatialVertex(DIM+1, randomLocation(), getTime());
		SpatialVertex** nearNeighbors = findMNearestNeighbors(newNode);
		
		for(long int i = 0; i < m; i++){	// for each of the new node's nearest neighbors
			
			newNode->addNeighbor(nearNeighbors[i]);	// link the new node and its neighbor
			
		}
		
		addNode(newNode);	// add the new node to the NSphere
		
		if(N < this->equalizationThreshold || N%(this->equalizationPeriod) == 0){
			
			equalize();	// equalize the sphere periodically
			
		}
		
		delete[] nearNeighbors;
		
		tick();
		n--;
		
	}
	
}

/**
 * uses the Box-Muller transform to randomly generate a position on the surface of an NSphere
 * dynamically allocates a double array of length DIM
 */
double* NSphere::randomLocation(){
	
	double* position = new double[DIM+1];
	double actual_radius = 0;
	double ratio;
	
	for(long int i = 0; i < DIM+1; i++){	// for every dimension
		
		double u = ((double)rand())/RAND_MAX, v = ((double)rand())/RAND_MAX;	// pick two values from a uniform distribution (0, 1)
		position[i] = sqrt(-2 * log(u)) * cos(2 * M_PI * v);	// set position component in dimension i to the result of the transform, distributed on N(0, 1)
		actual_radius += position[i] * position[i];	// track the actual radial distance of the position
		
	}
	
	actual_radius = sqrt(actual_radius);	// remember to take the square root of the sum of squares
	ratio = radius/actual_radius;	// to convert the vector to one of the appropriate radius, divide out the actual radius and multiply by the intended radius
	
	for(long int i = 0; i < DIM+1; i++){	// for every dimension
		
		position[i] *= ratio;	// apply the correction
		
	}
	
	
	return position;
	
}

/**
 * calculates the linear distance between two SpatialVertices
 * assumes that both SpatialVertices are of the dimension specified by the NSphere (DIM+1)
 * otherwise behavior is undefined and possibly fatal
 */
double NSphere::linearDistance(Vertex* a, Vertex* b){
	
	return linearDistance(((SpatialVertex*)a)->position, ((SpatialVertex*)b)->position);
	
}

/**
 * method to calculate the linear distance between two position vectors
 * assumes that both arrays are of size DIM+1, such that they are of the dimension of the NSphere's space
 * otherwise behavior is undefined and possibly fatal
 */
double NSphere::linearDistance(double* a, double* b){
	
	double sumOfSquares = 0;
	
	for(long int i = 0; i < DIM+1; i++){
		
		sumOfSquares += pow(a[i] - b[i], 2.0);
		
	}
	
	return sqrt(sumOfSquares);
	
}

/**
 * method to find and return an array of the m spatially nearest neighbors to a given starting node (excluding itself)
 * returns a dynamically allocated array which should be deleted after use
 */
SpatialVertex** NSphere::findMNearestNeighbors(SpatialVertex* start){
	
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
 * calculates the sum of the forces on a particular SpatialVertex
 * currently uses a force law of 1/r^(DIM-1)
 * returns a dynamically allocated array (force vector) which must be deleted after use
 */
double* NSphere::sumForces(SpatialVertex* node){
	
	double* force = new double[DIM+1];	// allocate a new force vector of a size appropriate to the space
	SpatialVertex* other;	// local placeholder for the other node in two-body interactions
	double magnitude, dist;	// local fields to hold the magnitude of the force and distance between two nodes
	
	for(long int i = 0; i < DIM+1; i++){	// first, ensure that all components of the force vector are zero
		
		force[i] = 0;
		
	}
	
	for(long int i = 0; i < N; i++){	// for all nodes on the NSphere
		
		if(node == getNode(i)){	// except the node for which force is being calculated
			
			continue;	// seriously, skip the infinite self-force
			
		}
		
		other = getNode(i);
		dist = linearDistance(node, other);
		magnitude = 1.0 / pow(dist, (double)DIM);	// calculate the magnitude of the force as 1/r^DIM
		
		for(long int j = 0; j < DIM+1; j++){	// for each dimension, add the component of force from this interaction
			
			double unitComponent = (node->position[j] - other->position[j]) / dist;	// the component of the unit
																						// vector is r(j)/|r|
			force[j] += unitComponent * magnitude;	// the component of the force in the direction of dimension j
													// is the value of the component of the unit vector
													// multiplied by the magnitude of the force vector
													// unit vector * magnitude = vector
			
		}
		
	}
	
	return force;
	
}

void NSphere::equalize(){
	
	gradientDescent(baseGam/pow(N, 1 + 1.0/DIM), baseTol, baseItr);
	
}

void NSphere::gradientDescent(double gamma, double tolerance, long int maxItr){

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
				
				for(long int j = 0; j < DIM+1; j++){	// for every dimension
				
					oldPos[j] = getNode(i)->position[j];	// store the old position
					getNode(i)->position[j] += gamma * netForce[i][j];	// displace the node in that dimension
																		// according to the force on it and the timestep size
					
				}
				
				normalizeRadius(getNode(i));	// normalize the radius
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

/**
 * method to normalize the radial distance of a given node to the radius of the NSphere
 */
void NSphere::normalizeRadius(SpatialVertex* node){
	
	double actual_radius_squared = 0;
	double ratio;	// the ratio of the desired radius to the actual radius
	
	for(long int i = 0; i < DIM+1; i++){	// for each dimension
		
		actual_radius_squared += node->position[i] * node->position[i];	// sum the square of the appropriate component
		
	}
	
	ratio = radius / sqrt(actual_radius_squared);
	
	for(long int i = 0; i < DIM+1; i++){	// for each dimension
		
		node->position[i] *= ratio;	// multiply the appropriate component by the ratio of radii
									// |vector * scalar/|vector|| = |vector/|vector| * scalar| = |unit vector * scalar| = scalar
		
	}
	
}
