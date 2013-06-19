using namespace std;

/**
 * uses the Box-Muller transform to randomly generate a position on the surface of an NSphere
 * dynamically allocates a double array of length DIM
 */
double* NSphere::randomLocation(){
	
	double* position = new double[DIM];
	double actual_radius = 0;
	double ratio;
	
	for(long int i = 0; i < DIM; i++){	// for every dimension
		
		double u = ((double)rand())/RAND_MAX, v = ((double)rand())/RAND_MAX;	// pick two values from a uniform distribution (0, 1)
		position[i] = sqrt(-2 * ln(u)) * cos(2 * M_PI * v);	// set position component in dimension i to the result of the transform, distributed on N(0, 1)
		actual_radius += position[i] * position[i];	// track the actual radial distance of the position
		
	}
	
	actual_radius = sqrt(actual_radius);	// remember to take the square root of the sum of squares
	ratio = radius/actual_radius;	// to convert the vector to one of the appropriate radius, divide out the actual radius and multiply by the intended radius
	
	for(long int i = 0; i < DIM; i++){	// for every dimension
		
		position[i] *= ratio;	// apply the correction
		
	}
	
	
	return position;
	
}

/**
 * calculates the square of the linear distance between two SpatialVertices
 * assumes that both SpatialVertices are of the dimension specified by the NSphere
 * otherwise behavior is undefined and possibly fatal
 */
double NSphere::distanceSquared(SpatialVertex* a, SpatialVertex* b){
	
	double square_sum = 0;
	
	for(long int i = 0; i < DIM; i++){	// for all components of the two position vectors
		
		square_sum += (a->position[i] - b->position[i]) * (a->position[i] - b->position[i]);	// add the square
																								// of the difference
		
	}
	
	return square_sum;	// return the resulting sum
	
}

/**
 * calculates the linear distnace between two SpatialVertices
 * assumes that both SpatialVertices are of the dimension specified by the NSphere
 * otherwise behavior is undefined and possibly fatal
 */
double NSphere::distance(SpatialVertex* a, SpatialVertex* b){
	
	return sqrt(distanceSquared(a, b));
	
}

/**
 * calculates the sum of the forces on a particular SpatialVertex
 * currently uses a force law of 1/r^(DIM-1)
 * returns a dynamically allocated array (force vector) which must be deleted after use
 */
double* sumForces(SpatialVertex* node){
	
	double* force = new double[DIM];	// allocate a new force vector of a size appropriate to the space
	SpatialVertex* other;	// local placeholder for the other node in two-body interactions
	double magnitude, dist;	// local fields to hold the magnitude of the force and distance between two nodes
	
	for(long int i = 0; i < DIM; i++){	// first, ensure that all components of the force vector are zero
		
		force[i] = 0;
		
	}
	
	for(long int i = 0; i < N; i++){	// for all nodes on the NSphere
		
		if(node == getNode(i)){	// except the node for which force is being calculated
			
			continue;	// seriously, skip the infinite self-force
			
		}
		
		other = getNode(i);
		dist = distance(node, other);
		magnitude = 1.0 / pow(dist, DIM-1);	// calculate the magnitude of the force as 1/r^(DIM-1)
		
		for(long int j = 0; j < DIM; j++){	// for each dimension, add the component of force from this interaction
			
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

/**
 * method to normalize the radial distance of a given node to the radius of the NSphere
 */
void NSphere::normalizeRadius(SpatialVertex* node){
	
	double actual_radius_squared = 0;
	double ratio;	// the ratio of the desired radius to the actual radius
	
	for(long int i = 0; i < DIM; i++){	// for each dimension
		
		actual_radius_squared += node->position[i] * node->position[i];	// sum the square of the appropriate component
		
	}
	
	ratio = radius / sqrt(radius_squared);
	
	for(long int i = 0; i < DIM; i++){	// for each dimension
		
		node->position[i] *= ratio;	// multiply the appropriate component by the ratio of radii
									// |vector * scalar/|vector|| = |vector/|vector| * scalar| = |unit vector * scalar| = scalar
		
	}
	
}
