#include "nsphere.h"
#include <cstring>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <cfloat>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <fstream>

using namespace std;

NSphere::NSphere(long int n, long int m, int d, float r, float baseGam, float baseTol, long int baseItr, long int threshold, long int period) : NGraph(d){
	
	static bool randSeeded = false;
	
	if(!randSeeded){	// if the random number generator has not yet been seeded
		
		srand(std::time(NULL));	// do so with the current time as the seed
		randSeeded = true;	// make a note for future initializations
		
	}
	
	this->time = 0;
	this->radius = r;
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
		
		tick();
		
	}
	
	equalize();
	grow(n - (m+1));	// grow the remaining nodes normally
	
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
		
		if(m > DIM){	// if there are enough nearest neighbors
						// set the position of this node to their center
			
			for(long int i = 0; i < DIM+1; i++){	// for each dimension
				
				newNode->position[i] = 0;
				
				for(long int j = 0; j < m; j++){	// slightly pre-equalize the new node
					
					newNode->position[i] += (nearNeighbors[j]->position[i])/m;	// by placing it at the average position of its nearest neighbors
					
				}
					
			}
			
		}
		
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
 * dynamically allocates a float array of length DIM
 */
float* NSphere::randomLocation(){
	
	float* position = new float[DIM+1];
	float actual_radius = 0;
	float ratio;
	
	for(long int i = 0; i < DIM+1; i++){	// for every dimension
		
		float u = ((float)rand())/RAND_MAX, v = ((float)rand())/RAND_MAX;	// pick two values from a uniform distribution (0, 1)
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
float NSphere::linearDistance(Vertex* a, Vertex* b){
	
	return linearDistance(((SpatialVertex*)a)->position, ((SpatialVertex*)b)->position);
	
}

/**
 * method to calculate the linear distance between two position vectors
 * assumes that both arrays are of size DIM+1, such that they are of the dimension of the NSphere's space
 * otherwise behavior is undefined and possibly fatal
 */
float NSphere::linearDistance(float* a, float* b){
	
	float sumOfSquares = 0;
	
	for(long int i = 0; i < DIM+1; i++){
		
		sumOfSquares += (a[i] - b[i]) * (a[i] - b[i]);
		
	}
	
	return sqrt(sumOfSquares);
	
}

/**
 * method to find and return an array of the m spatially nearest neighbors to a given starting node (excluding itself)
 * returns a dynamically allocated array which should be deleted after use
 */
SpatialVertex** NSphere::findMNearestNeighbors(SpatialVertex* start){
	
	SpatialVertex** near = new SpatialVertex*[m];
	float dist[m];	//local record of the distance from the start of the m nearest neighbors

	for(int i = 0; i < m; i++){

		dist[i] = FLT_MAX;
	}

	for(int i = 0; i < N; i++){

		if(getNode(i) == start){	//do not consider the distance between the starting node and itself

			continue;

		}

		float d = linearDistance(start, getNode(i));

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
float* NSphere::sumForces(SpatialVertex* node){
	
	float* force = new float[DIM+1];	// allocate a new force vector of a size appropriate to the space
	SpatialVertex* other;	// local placeholder for the other node in two-body interactions
	float magnitude, dist;	// local fields to hold the magnitude of the force and distance between two nodes
	
	for(long int i = 0; i < DIM+1; i++){	// first, ensure that all components of the force vector are zero
		
		force[i] = 0;
		
	}
	
	for(long int i = 0; i < N; i++){	// for all nodes on the NSphere
		
		if(node == getNode(i)){	// except the node for which force is being calculated
			
			continue;	// seriously, skip the infinite self-force
			
		}
		
		other = getNode(i);
		dist = linearDistance(node, other);
		magnitude = 1.0 / pow(dist, (float)DIM);	// calculate the magnitude of the force as 1/r^DIM
		
		for(long int j = 0; j < DIM+1; j++){	// for each dimension, add the component of force from this interaction
			
			float unitComponent = (node->position[j] - other->position[j]) / dist;	// the component of the unit
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
	
	gradientDescent(baseGam / pow(N, 1 + 1.0/DIM), baseTol, baseItr);
	
}

void NSphere::gradientDescent(float gamma, float tolerance, long int maxItr){

	float* netForce[N];	// a record of the net force on each node
	float prevMaxDisp = FLT_MAX;	// record of the maximum scalar displacement of a node on the previous iteration
									// initially set to an impossibly large value to ensure at least one iteration occurs
	float tolMaxDisp = (radius * pow(N, -1.0/DIM) * tolerance) * baseGam;	// the maximum tolerated displacement
																				// if the maximum displacement is above this value, continue iterating
																				// scales with gamma to provide equal treatment across timestep sizes

	memset(netForce, 0, N * sizeof(float*));
	
	while(prevMaxDisp > tolMaxDisp && maxItr > 0){	// while there remains some excess energy above tolerance
													// and the hard limit of iterations has not been passed

		float maxDisp = 0;

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
			
				float oldPos[DIM];
				
				for(long int j = 0; j < DIM+1; j++){	// for every dimension
				
					oldPos[j] = getNode(i)->position[j];	// store the old position
					getNode(i)->position[j] += gamma * netForce[i][j];	// displace the node in that dimension
																		// according to the force on it and the timestep size
					
				}
				
				normalizeRadius(getNode(i));	// normalize the radius
				float disp = linearDistance(oldPos, getNode(i)->position);	// calculate the displacement
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
	
	float actual_radius_squared = 0;
	float ratio;	// the ratio of the desired radius to the actual radius
	
	for(long int i = 0; i < DIM+1; i++){	// for each dimension
		
		actual_radius_squared += node->position[i] * node->position[i];	// sum the square of the appropriate component
		
	}
	
	ratio = radius / sqrt(actual_radius_squared);
	
	for(long int i = 0; i < DIM+1; i++){	// for each dimension
		
		node->position[i] *= ratio;	// multiply the appropriate component by the ratio of radii
									// |vector * scalar/|vector|| = |vector/|vector| * scalar| = |unit vector * scalar| = scalar
		
	}
	
}

void NSphere::exportObject(const NSphere* ns, const char* filename){
	
	ofstream outfile(filename, ios::out | ios::trunc);	// open the output file
	
	// output all of the parameters of the model
	outfile<<"# dimension = "<<ns->DIM<<endl;
	outfile<<"# radius = "<<ns->radius<<endl;
	outfile<<"# base gamma = "<<ns->baseGam<<endl;
	outfile<<"# base tolerance = "<<ns->baseTol<<endl;
	outfile<<"# base iterations = "<<ns->baseItr<<endl;
	outfile<<"# equalization threshold = "<<ns->equalizationThreshold<<endl;
	outfile<<"# equalization period = "<<ns->equalizationPeriod<<endl;
	outfile<<"# iteration weights = "<<ns->iterationWeights<<endl;
	outfile<<"# time = "<<ns->getTime()<<endl;
	outfile<<"# m = "<<ns->m<<endl;
	outfile<<"# N = "<<ns->N<<endl;
	
	outfile<<endl;
	
	for(long int i = 0; i < ns->N; i++){	// for each node in the graph
		
		outfile<<i<<" ";	// output the index of the node, now its only identifier
		
		outfile<<ns->getNode(i)->dimension<<" ";	// output the dimension of the node so that is position can be accurately read
		
		outfile<<ns->getNode(i)->getStartTime()<<" ";	// output the starting time of this node
		
		for(long int j = 0; j < ns->getNode(i)->dimension; j++){	// for each dimension of the space
			
			outfile<<ns->getNode(i)->position[j]<<" ";	// output the position
			
		}
		
		for(long int j = 0; j < ns->K(i); j++){	// for each neighbor of this node
			
			long int index = ns->indexOf(ns->getNode(i)->getNeighbor(j));	// identify it by index
			
			outfile<<index<<" ";
			
		}
		
		outfile<<endl;	// signify the end of this node's information with a new line
		
	}
	
	outfile.close();	// close the output file
	
}

NSphere* NSphere::importObject(const char* filename){
	
	ifstream infile(filename, ios::in);
	
	NSphere* ns;
	long int bufsize = 2048;
	long int m, time, dim, iterations, threshold, period, size;
	float radius, gamma, tolerance, weights;
	char* line = new char[bufsize];	// create a character buffer larger than could be reasonably used
	bool** adjacency;	// create an adjacency matrix to record edges
	
	// retrieve the parameters of the model as they were output
	infile.getline(line, bufsize);
	sscanf(line, "# dimension = %li", &dim);
	memset(line, 0, bufsize * sizeof(char));
	
	infile.getline(line, bufsize);
	sscanf(line, "# radius = %f", &radius);
	memset(line, 0, bufsize * sizeof(char));
	
	infile.getline(line, bufsize);
	sscanf(line, "# base gamma = %f", &gamma);
	memset(line, 0, bufsize * sizeof(char));
	
	infile.getline(line, bufsize);
	sscanf(line, "# base tolerance = %f", &tolerance);
	memset(line, 0, bufsize * sizeof(char));
	
	infile.getline(line, bufsize);
	sscanf(line, "# base iterations = %li", &iterations);
	memset(line, 0, bufsize * sizeof(char));
	
	infile.getline(line, bufsize);
	sscanf(line, "# equalization threshold = %li", &threshold);
	memset(line, 0, bufsize * sizeof(char));
	
	infile.getline(line, bufsize);
	sscanf(line, "# equalization period = %li", &period);
	memset(line, 0, bufsize * sizeof(char));
	
	infile.getline(line, bufsize);
	sscanf(line, "# iteration weights = %f", &weights);
	memset(line, 0, bufsize * sizeof(char));
	
	infile.getline(line, bufsize);
	sscanf(line, "# time = %li", &time);
	memset(line, 0, bufsize * sizeof(char));
	
	infile.getline(line, bufsize);
	sscanf(line, "# m = %li", &m);
	memset(line, 0, bufsize * sizeof(char));
	
	infile.getline(line, bufsize);
	sscanf(line, "# N = %li", &size);
	memset(line, 0, bufsize * sizeof(char));
	
	// create a new NSphere (with m+1 nodes so as not to anger the constructor)
	ns = new NSphere(m+1, m, dim, radius, gamma, tolerance, iterations, threshold, period);
	ns->time = time;
	
	while(ns->N > 0){	// remove any nodes in the NSphere
		
		ns->nodes.pop_back();
		
	}
	
	adjacency = new bool*[size];	// create an array of N adjacency arrays
	
	while(!infile.eof()){	// while there remain lines in the file
	
		long int index, dimension, starttime;
		float* position;
		SpatialVertex* node;
		
		infile.getline(line, bufsize);	// retrieve them from the file
		
		vector<string> tokens;	// create a vector of string tokens
		istringstream linestream(line);	// and execute the appropriate steps
		copy(istream_iterator<string>(linestream), istream_iterator<string>(), back_inserter<vector <string> >(tokens));	// to read each whitespace-delimited string as a separate token
		
		if(tokens.size() == 0){	// if this is a line of whitespace, or otherwise contains no information
			
			memset(line, 0, bufsize * sizeof(char));	// delete the information
			continue;	// do not attempt to parse it
			
		}
		
		index = atol(tokens.at(0).c_str());	// record the index to properly populate the adjacency matrix
		dimension = atol(tokens.at(1).c_str());	// record the dimension to read the appropriate number of position components
		starttime = atol(tokens.at(2).c_str());	// record the starttime to use as a parameter to the new node's constructor
		
		position = new float[dimension];	// initialize the new position array
		adjacency[index] = new bool[size];	// initialize this row/column of the adjacency matrix
		memset(adjacency[index], false, size * sizeof(bool));
		
		for(long int i = 3; i < 3 + dimension; i++){	// for each component of the node's position
			
			position[i - 3] = atof(tokens.at(i).c_str());	// the first component is the fourth value given
														// since the tokens vector is zero-indexed, this means that the position starts from index 3
			
		}
		
		for(long int i = 3 + dimension; i < tokens.size(); i++){	// for all remaining tokens
			
			long int neighborindex = atol(tokens.at(i).c_str());	// record the index given by the token
			adjacency[index][neighborindex] = true;	// record that this node (index) is linked to some other node (neighborindex)
			
		}
		
		node = new SpatialVertex(dimension, position, starttime);
		ns->addNode(node);
		
		memset(line, 0, bufsize * sizeof(char));
		
	}
	
	for(long int i = 0; i < size; i++){	// for each adjacency list in the adjacency matrix
		
		for(long int j = 0; j < size; j++){	// for each element in that list
			
			if(adjacency[i][j]){	// if that element is true (i.e. nodes i and j are linked)
				
				ns->getNode(i)->addNeighbor(ns->getNode(j));	// create an edge (if one does not exist aleady)
				
			}
			
		}
		
	}
	
	infile.close();
	
	for(long int i = 0; i < size; i++){
		
		delete adjacency[i];
		
	}
	
	delete adjacency;
	delete line;
	
	return ns;
	
}
