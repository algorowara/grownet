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
	
	for(long int i = 0; i < DIM; i++){
		
		square_sum += (a->position[i] - b->position[i]) * (a->position[i] - b->position[i]);
		
	}
	
	return square_sum;
	
}

/**
 * calculates the linear distnace between two SpatialVertices
 * assumes that both SpatialVertices are of the dimension specified by the NSphere
 * otherwise behavior is undefined and possibly fatal
 */
double NSphere::distance(SpatialVertex* a, SpatialVertex* b){
	
	return sqrt(distanceSquared(a, b));
	
}
