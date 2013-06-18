#include "../growingnetwork3d/growingnetwork3d.h"
#include <iostream>

using namespace std;

GrowingNetwork3D::GrowingNetwork3D(long int n, long int m, double gam, double tol, long int itr){
	
	static bool randSeeded = false;
	
	if(!randSeeded){	// if the random number generator has not yet been seeded
		
		srand(std::time(NULL));	// do so with the current time as the seed
		randSeeded = true;	// make a note for future initializations
		
	}
	
	this->radius = 1;
	this->m = m;
	baseGam = gam;
	baseTol = tol;
	baseItr = itr;
	
	for(long int i = 0; i < m+1; i++){	// for the first m+1 nodes, which form a clique
	
		SpatialVertex* newNode = new SpatialVertex(DIM, randomLocation(), getTime());	// generate a new node in DIM dimensions
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

SpatialVertex* GrowingNetwork3D::getNode(long int i){
	
	return (SpatialVertex*)nodes.at(i);
	
}

void GrowingNetwork3D::grow(long int n){
	
	while(n > 0){
		
		SpatialVertex* newNode = new SpatialVertex(DIM, randomLocation(), getTime());
		SpatialVertex** nearNeighbors = findMNearestNeighbors(newNode);
		
		for(int i = 0; i < m; i++){
			
			newNode->addNeighbor(nearNeighbors[i]);
			
		}
		
		addNode(newNode);		
		equalize();
		
		delete[] nearNeighbors;
		tick();		
		n--;
		
	}
	
}

/**
 * randomly select a point on the surface of a sphere
 * using theta = 2 * pi * u, phi = acos(2 * v - 1), where
 * u and v are taken uniformly from the interval (0, 1)
 * see doc.odt for specification of the coordinate system
 */
double* GrowingNetwork3D::randomLocation(){
	
	double* position = new double[DIM];
	
	double theta = 2 * M_PI * (((double)rand())/RAND_MAX);
	double phi = acos((2 * ((double)rand())/RAND_MAX - 1));
	
	position[0] = radius * sin(phi) * cos(theta);
	position[1] = radius * sin(phi) * sin(theta);
	position[2] = radius * cos(phi);
	
	return position;
	
}

/**
 * method to calculate the linear distance between two nodes
 * if the nodes are not of a type to have a position in space
 * return 0, as they are unknown
 */
double GrowingNetwork3D::linearDistance(Vertex* a, Vertex* b){
	
	if(sizeof(a) == sizeof(SpatialVertex) && sizeof(b) == sizeof(SpatialVertex)){
	
		SpatialVertex* a_loc = (SpatialVertex*)a;
		SpatialVertex* b_loc = (SpatialVertex*)b;
	
		return DISTANCE(a_loc, b_loc);
	
	}
	
	return 0;
	
}

/**
 * method to calculate the repulsive force between two nodes based on their displacement
 * currently falls off as 1/r^2
 * dynamically allocates an array which should be deleted by the caller after use
 */
double* GrowingNetwork3D::sumForces(SpatialVertex* node){
	
	double* force = new double[DIM];	// allocate a new array for the force vector
	SpatialVertex* other;	// local placeholder for any other node in two-body interactions
	double magnitude, distance;	// local placeholders for the magnitude of a force and the distance between two nodes
	
	for(int i = 0; i < DIM; i++){	// ensure all components are set to zero first
		
		force[i] = 0;
		
	}
	
	for(int i = 0; i < N; i++){	// for every node in the graph
		
		if(getNode(i) == node){	// except this one
			
			continue;	// skip this one
			
		}
		
		other = getNode(i);
		magnitude = 1.0 / DISTANCE_SQUARED(node, other);	// calculate the unitless magnitude of the repulsive force
		distance = DISTANCE(node, other);
		
		force[0] += magnitude * (X(node) - X(other))/distance;	// multiply the magnitude
		force[1] += magnitude * (Y(node) - Y(other))/distance;	// by the unit vector
		force[2] += magnitude * (Z(node) - Z(other))/distance;	// of the displacement
		
	}
	
	return force;
	
}

/**
 * method to find the m nearest neighbors to some starting node
 */
SpatialVertex** GrowingNetwork3D::findMNearestNeighbors(SpatialVertex* start){
	 
	SpatialVertex** near = new SpatialVertex*[m];
	double dsquare[m];	// local record of the distance-squared of the m nearest neighbors
	
	for(int i = 0; i < m; i++){
	
		dsquare[i] = DBL_MAX;
		
	}
	
	
	for(int i = 0; i < N; i++){	// for all nodes in the graph
		
		if(getNode(i) == start){	// do not consider the distance between the starting node and itself
		
			continue;
			
		}
	
		double square = DISTANCE_SQUARED(start, getNode(i));	// calculate the square of the distance
																// as minimizing distance squared is the same
																// as minimizing distance (as a signless scalar)
		
		for(int j = 0; j < m; j++){	// iterate through all distance-squared records
			
			if(square < dsquare[j]){	// if this node is closer than the nearest neighbor j
				
				for(int k = m-1; k > j; k--){	// shift everything from j to m one position back
					
					dsquare[k] = dsquare[k-1];
					near[k] = near[k-1];
					
				}
				
				dsquare[j] = square;	// put the data for this new nearest neighbor
				near[j] = getNode(i);	// in the spot previously occupied by neighbor j
				
				break;	// stop iterating through records
				
			}
			
		}
		 
	}
	
	return near;
	 
}

/**
 * method to ensure that the radial distance of a node from the origin
 * is that specified in the radius field
 * without disturbing the unit direction
 */
void GrowingNetwork3D::normalizeRadius(SpatialVertex* node){
	
	// find the ratio of the ideal radius to the current radial distance
	double ratio = radius / sqrt(X(node) * X(node) + Y(node) * Y(node) + Z(node) * Z(node));
	
	for(long int i = 0; i < DIM; i++){	// for each dimension
		
		node->position[i] *= ratio;	// multiply it by that ratio
		
	}
	
}

/**
 * method to calculate the potential of all the nodes in the network
 * where the potential of any node pair is defined as 1/r
 * where r is the magnitude of the displacement between the two nodes
 */
double GrowingNetwork3D::calculatePotential(){
	
	double potential = 0;
	
	#pragma omp parallel shared(potential)
	{
		
		double localSum = 0;	// sum of potential energies local to this thread
		SpatialVertex* a;
		SpatialVertex* b;
		
		#pragma omp for schedule(guided)
		
		for(int i = 0; i < N; i++){	// for every node
			
			for(int j = i+1; j < N; j++){	// for every node which has not yet been visited by the outer loop
				
				a = getNode(i);
				b = getNode(j);				
				localSum += 1.0/DISTANCE(a, b);	// add their potential energy to the local sum
								
			}
			
		}
		
		#pragma omp atomic
		potential += localSum;	// add the local sums together
		
	}
	
	return potential;
	
}

void GrowingNetwork3D::equalize(){
	
	// externally provided solutions to the Thompson problem show Energy proportional to N^2
	// force and average separation should both be proportional to 1/N, so gamma should be proportional to 1/N^2
	
	gradientDescent(baseGam/(N * N), baseTol, baseItr);
	
}

/**
 * uses a gradient descent algorithm where the potential is 1/r
 * where gamma is the conversion factor from force to movement; delta(position) = gamma * netForce
 * where baseTolerance is the maximum ratio of the difference between the potential and the minimum, and the minimum
 * and maxItr is the maximum number of iterations allowed before the function exits, regardless of tolerance
 */
void GrowingNetwork3D::gradientDescent(double gamma, double baseTolerance, long int maxItr){

	if(gamma < 0){	// check for negative gamma, which would turn gradient descent into gradient ascent
	
		gamma *= -1;	// fix if necessary, assume the negative sign is in error
		
	}
	
	if(baseTolerance < 0){	// check for negative tolerance, which would demand an energy below the minimum
		
		baseTolerance *= -1;	// fix if necessary, assume the negative sign is in error
		
	}
	
	if(maxItr < 0){	// check for negative iterations, which is just silly
		
		maxItr *= -1;	// fix the darn thing
		
	}
	
	double* netForce[N];	// local array to store the net forces on each node
	double previousPotential = DBL_MAX;	// local field to store the previous known potential
	double minPot = calculateMinimumPotential();	// storage of the minimum potential for this N
	double tolerance = minPot * baseTolerance;	// actual value of tolerance
	
	while(previousPotential - minPot > tolerance && maxItr > 0){

		#pragma omp parallel shared(netForce)
		{
			
			#pragma omp for schedule(guided)
			for(int i = 0; i < N; i++){	// for every node
				
				netForce[i] = sumForces(getNode(i));	// store the net force on the node
														// thread safety is not an issue here
														// because the values of i are divided among threads
														// and no two will ever write to the same index
				
			}
			
		}
		
		#pragma omp parallel shared(netForce)
		{
			
			#pragma omp for schedule(guided)
			for(int i = 0; i < N; i++){	// for every node
				
				for(int j = 0; j < DIM; j++){	// for every dimension
					
					getNode(i)->position[j] += gamma * netForce[i][j];	// displace the node by gamma * netForce
					
				}
				
				normalizeRadius(getNode(i));	// return the node to the surface of the sphere
				delete[] netForce[i];	// delete the now-outdated force vector for this node
				
			}
			
		}
	
		previousPotential = calculatePotential();	// update the record of the latest potential
		maxItr--;	// move one iteration closer to ending regardless of the potential energy
	
	}
	
}

/**
 * return the potential energy as a function of N
 * as calculated from the equation given in the design document
 */
double GrowingNetwork3D::calculateMinimumPotential(){
	
	return 0.063594969640041382 * sqrt(N) + -0.55213813866389005 * pow(N, 1.5) + 0.49998893897252450 * (N * N);
	
}

GrowingNetwork3D::~GrowingNetwork3D(){
	
}
