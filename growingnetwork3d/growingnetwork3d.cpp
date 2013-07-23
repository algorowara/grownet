#include "../growingnetwork3d/growingnetwork3d.h"

using namespace std;

/**
 * empty constructor which exists solely for derived classes to use without conflict
 */
GrowingNetwork3D::GrowingNetwork3D(){
	
}

/**
 * constructor method specifying:
 *   n, the number of nodes to be grown in the network
 *   m, the number of new edges per newly created node
 *   gam, the base timestep size (the actual timestep size will scale as a function of N)
 *   tol, the base tolerance in terms of a fraction of the estimated minimum potential difference between N-1 and N
 *   itr, the maximum number of iterations to be used in gradientDescent before giving up
 */
GrowingNetwork3D::GrowingNetwork3D(long int n, long int m, float gam, float tol, long int itr) : DIM(3){
	
	this->time = 0;
	this->radius = DEFAULT_RADIUS;
	this->m = m;
	baseGam = gam;
	baseTol = tol;
	baseItr = itr;
	
	for(long int i = 0; i < m+1; i++){	// for the first m+1 nodes, which form a clique
	
		SpatialVertex* newNode = new SpatialVertex(DIM, randomLocation(), getTime());	// generate a new node in DIM dimensions
		addNode(newNode);	// place it in the network
		
		for(long int j = 0; j < nodes.size()-1; j++){	// link it to all previously created nodes
			
			newNode->addNeighbor(getNode(j));
			
		}
		
		
		equalize();	// equalize the distribution of nodes
		tick();	// advance the current time
		
	}
	
	grow(n  - (m+1));	// grow the remaining nodes, less these m+1 initial nodes, normally
	
}

/**
 * method to retrieve a node and automatically cast it to a SpatialVertex
 */
SpatialVertex* GrowingNetwork3D::getNode(long int i) const{
	
	return (SpatialVertex*)nodes.at(i);
	
}

/**
 * method to grow a given number of nodes and equalize them using gradient descent
 */
void GrowingNetwork3D::grow(long int n){
	
	while(n > 0){	// while there remain nodes to be added
		
		SpatialVertex* newNode = new SpatialVertex(DIM, randomLocation(), getTime());	// create a new node at a random location
		SpatialVertex** nearNeighbors = findMNearestNeighbors(newNode);	// find its m nearest neighbors
		
		for(int i = 0; i < m; i++){	// for each such neighbor
			
			newNode->addNeighbor(nearNeighbors[i]);	// create an edge between it and the new node
			
		}
		
		addNode(newNode);	// add the new node to the network so that it can participate in inter-node interactions
		equalize();	// perform gradient descent on the nodes to space them evenly
		
		delete[] nearNeighbors;	// deallocate the array of nearest neighbors
		tick();	// advance the current time of the network
		n--;	// recognize that one less node needs to be added
		
	}
	
	return;
	
}

/**
 * randomly select a point on the surface of a sphere
 * using theta = 2 * pi * u, phi = acos(2 * v - 1), where
 * u and v are taken uniformly from the interval (0, 1)
 * see doc.odt for specification of the coordinate system
 */
float* GrowingNetwork3D::randomLocation() const{
	
	float* position = new float[DIM];	// allocate a position vector
	
	float theta = 2 * M_PI * (((float)rand())/RAND_MAX);	// randomly select a theta
	float phi = acos((2 * ((float)rand())/RAND_MAX - 1));	// and a phi
	
	position[0] = radius * sin(phi) * cos(theta);
	position[1] = radius * sin(phi) * sin(theta);
	position[2] = radius * cos(phi);
	
	return position;
	
}

/**
 * method to quickly calculate the square of the distance between two SpatialVertices in three dimensions
 * used either to compare distance for ordering or for inverse-square forces
 * for uses which actually require distance, use the distance method
 * if a and b are not both three-dimensional, this method's behavior is undefined
 */
float GrowingNetwork3D::distanceSquared(SpatialVertex* a, SpatialVertex* b) const{
	
	return ((X(a) - X(b)) * (X(a) - X(b)) + (Y(a) - Y(b)) * (Y(a) - Y(b)) + (Z(a) - Z(b)) * (Z(a) - Z(b)));
	
}

/**
 * method to calculate the distance (as a scalar) between two SpatialVertices in three dimensions
 * if a and b are not both three-dimensional, this method's behavior is undefined
 * relies on the distanceSquared method
 */
float GrowingNetwork3D::distance(SpatialVertex* a, SpatialVertex* b) const{
	
	return sqrt(distanceSquared(a, b));
	
}

/**
 * method to calculate the linear distance between two nodes
 * if the nodes are not of a type to have a position in space
 * return 0, as they are unknown
 * this method exists solely to satisfy the requirements of base classes
 * for all other uses, use the distance method
 */
float GrowingNetwork3D::linearDistance(Vertex* a, Vertex* b){

	SpatialVertex* a_loc = (SpatialVertex*)a;
	SpatialVertex* b_loc = (SpatialVertex*)b;
		
	return distance(a_loc, b_loc);
	
}

/**
 * method to calculate the repulsive force between two nodes based on their displacement
 * currently falls off as 1/r^2
 * dynamically allocates an array which should be deleted by the caller after use
 */
float* GrowingNetwork3D::sumForces(SpatialVertex* node) const{
	
	float* force = new float[DIM];	// allocate a new array for the force vector
	SpatialVertex* other;	// local placeholder for any other node in two-body interactions
	float magnitude, dist;	// local placeholders for the magnitude of a force and the distance between two nodes
	
	for(long int i = 0; i < DIM; i++){	// ensure all components are set to zero first
		
		force[i] = 0;
		
	}
	
	for(int i = 0; i < N; i++){	// for every node in the graph
		
		if(getNode(i) == node){	// except this one
			
			continue;	// skip this one
			
		}
		
		other = getNode(i);
		magnitude = 1.0 / distanceSquared(node, other);	// calculate the unitless magnitude of the repulsive force
		dist = distance(node, other);
		
		force[0] += magnitude * (X(node) - X(other))/dist;	// multiply the magnitude
		force[1] += magnitude * (Y(node) - Y(other))/dist;	// by the unit vector
		force[2] += magnitude * (Z(node) - Z(other))/dist;	// of the displacement
		
	}
	
	return force;
	
}

/**
 * method to find the m nearest neighbors to some starting node
 */
SpatialVertex** GrowingNetwork3D::findMNearestNeighbors(SpatialVertex* start) const{
	 
	SpatialVertex** near = new SpatialVertex*[m];	// dynamically allocated array to hold the m nearest neighbors
	float dsquare[m];	// local record of the distance-squared of the m nearest neighbors
	
	for(int i = 0; i < m; i++){
	
		dsquare[i] = FLT_MAX;
		
	}
	
	
	for(int i = 0; i < N; i++){	// for all nodes in the graph
		
		if(getNode(i) == start){	// do not consider the distance between the starting node and itself
		
			continue;
			
		}
	
		float square = distanceSquared(start, getNode(i));	// calculate the square of the distance
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
	float ratio = radius / sqrt(X(node) * X(node) + Y(node) * Y(node) + Z(node) * Z(node));
	
	for(long int i = 0; i < DIM; i++){	// for each dimension
		
		node->position[i] *= ratio;	// multiply it by that ratio
		
	}
	
}

/**
 * method to calculate the potential of all the nodes in the network
 * where the potential of any node pair is defined as 1/r
 * where r is the magnitude of the displacement between the two nodes
 */
float GrowingNetwork3D::calculatePotential() const{
	
	float potential = 0;
	
	#pragma omp parallel shared(potential)
	{
		
		float localSum = 0;	// sum of potential energies local to this thread
		SpatialVertex* a;
		SpatialVertex* b;
		
		#pragma omp for schedule(guided)
		
		for(int i = 0; i < N; i++){	// for every node
			
			for(int j = i+1; j < N; j++){	// for every node which has not yet been visited by the outer loop
				
				a = getNode(i);
				b = getNode(j);				
				localSum += 1.0/distance(a, b);	// add their potential energy to the local sum
								
			}
			
		}
		
		#pragma omp atomic
		potential += localSum;	// add the local sums together
		
	}
	
	return potential;
	
}

/**
 * call the gradient descent algorithm
 * with a timestep which scales with N^-3/2
 */
void GrowingNetwork3D::equalize(){
	
	gradientDescent(baseGam/(N * sqrt(N)), baseTol, baseItr);
	
}

/**
 * uses a gradient descent algorithm where the potential is 1/r
 * where gamma is the conversion factor from force to movement; delta(position) = gamma * netForce
 * where baseTolerance is the maximum ratio of the difference between the potential and the minimum, and the minimum
 * and maxItr is the maximum number of iterations allowed before the function exits, regardless of tolerance
 */
void GrowingNetwork3D::gradientDescent(float gamma, float baseTolerance, long int maxItr){
	
	float* netForce[N];	// local array to store the net forces on each node
	float previousPotential = FLT_MAX;	// local field to store the previous known potential; set to an arbitrary maximum to ensure that at least one iteration occurs
	float toleratedPotential = 0;	// if N is less than the guided N, the algorithm will perform the maximum number of iterations anyway
	
	if(N > GUIDED_N){	// otherwise, if N is greater, the tolerated potential is equal to the minimum potential plus some fraction of the estimated potential difference due to the addition of the last node
		toleratedPotential = calculateMinimumPotential(N, DIM) + baseTolerance * calculateMinimumPotentialDifference(N-1, N, DIM);
	}
	
	while(previousPotential > toleratedPotential && maxItr > 0){	// until either the potential is sufficiently low or the maximum number of iterations have been performed

		#pragma omp parallel shared(netForce)
		{
			
			#pragma omp for schedule(guided)
			for(int i = 0; i < N; i++){	// for every node
				
				netForce[i] = this->sumForces(getNode(i));	// store the net force on the node
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
 * return the potential energy as a function of 
 * as calculated from the equation given in the design document
 */
float GrowingNetwork3D::calculateMinimumPotential(long int n, long int d) const{
	
	return 0.063594969640041382 * sqrt(n) + -0.55213813866389005 * pow(n, 1.5) + 0.49998893897252450 * (n * n);
	
}

/**
 * return the difference between the minimum potentials at two different values of n
 */
float GrowingNetwork3D::calculateMinimumPotentialDifference(long int init_n, long int final_n, long int d) const{
	
	return calculateMinimumPotential(final_n, d) - calculateMinimumPotential(init_n, d);
	
}

/**
 * default destructor which does nothing beyond calling its base classes' destructors
 */
GrowingNetwork3D::~GrowingNetwork3D(){
	
}
