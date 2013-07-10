#include "../pp1d/pp1d.h"
#include <vector>
#include <iostream>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <assert.h>

using namespace std;

PPGrowingNetwork1D::PPGrowingNetwork1D(long int n, long int m, float gamma, float tolerance, long int maxItr) : DIM(1){

	static bool randSeeded = false;

	if(!randSeeded){

		srand(std::time(NULL));
		randSeeded = true;
	}

	this->gamma = gamma;
	this->tolerance = tolerance;
	this->maxItr = maxItr;
	this->m = m;
	radius = 1;
	alpha = 1;
	if(N == 0){
		beta = alpha;
	}
	else {
		beta = alpha*N;
	}

	for(long int i = 0; i < m+1; i++){ //for the first m+1 nodes, which form a clique

		SpatialVertex* newNode = new SpatialVertex(DIM, randomLocation(), getTime()); //generate a new node in DIM dimensions (1 here)

		addNode(newNode);

		beta = alpha*N; //update beta so positive and negative  charges are balanced

		for(long int j = 0; j < nodes.size()-1; j++){ //link to the others in this clique

			newNode->addNeighbor(getNode(j));

		}

		equalize(); //minimize the potential energy

		tick();
		n--;

	}

	grow(n);  //grow the remaining nodes normall

}

SpatialVertex* PPGrowingNetwork1D::getNode(long int i){

		return (SpatialVertex*)nodes.at(i);

}

void PPGrowingNetwork1D::grow(long int n){

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
		beta = alpha*N;

	}

}

/**
 * randomly select a point on the line using p_radius = radius * v, where v is taken uniformly from the interval [0,1]
*/
float* PPGrowingNetwork1D::randomLocation(){

	float* position = new float[DIM]; 
	
	float p_radius = radius * (((float)rand())/RAND_MAX);
	long int p_sign = (rand() % 2) - 1;
		if(p_sign == 0){
			p_sign = 1;
		}

	position[0] = p_radius * p_sign;

	return position;

}

/**
 * method to find linear distance, I don't use it but it is needed in here to compile properly
*/
float PPGrowingNetwork1D::linearDistance(Vertex* a, Vertex* b){
	if(sizeof(a) == sizeof(SpatialVertex) && sizeof(b) == sizeof(SpatialVertex)){

		SpatialVertex* a_loc = (SpatialVertex*)a;
		SpatialVertex* b_loc = (SpatialVertex*)b;

		return DISTANCE_1D(a_loc, b_loc);

	}
	
	return 0;

}

/**
 * calculate the repulsive force between two nodes based on their displacement. It dynamically allocated an array which should be deleted by caller after use
*/
float* PPGrowingNetwork1D::sumForces(SpatialVertex* node){

	float* force = new float[DIM]; //allocate a float for the force vector
	SpatialVertex* other; //local placeholder for any other node in 2 body interactions
	float magnitude, distance, pmagnitude, pdistance; //local placeholders for magnitude of force and location of nodes for electron - electron interactions and electron cloud interactions

	force[0] = 0; //force is initially 0

	for(int i = 0; i < N; i++){  //for each node in the graph

		if(getNode(i) == node){  //except this one

			continue;

		}

		other = getNode(i);
		magnitude = alpha; //unitless magnitude of electron electron force, it is constant in 1D
		distance = DISTANCE_1D(node, other);
		force[0] += magnitude * (X(node) - X(other))/distance; //multiply the magnitude by the unit vector of displacement

	}
	pdistance = X(node); //the distance from the origin
	pmagnitude = - beta * pdistance; //unitless magnitude of attractive electron-cloud interaction

	force[0] += pmagnitude; //we don't even need to use a unit vector, it points radially outwards

	return force;

}

/**
 * method to find the m nearest neighbors to some starting node
 */
SpatialVertex** PPGrowingNetwork1D::findMNearestNeighbors(SpatialVertex* start){

	SpatialVertex** near = new SpatialVertex*[m];
	float dnormal[m]; //local record of the distance of the m nearest neighbors

	for(int i = 0; i < m; i++){

		dnormal[i] = DBL_MAX; //everything is initially far away

	}

	for(int i = 0; i < N; i++){

		if(getNode(i) == start){  //do not consider the starting node...

			continue;

		}

		float dist = DISTANCE_1D(start, getNode(i));

		for(int j = 0; j < m; j++){  //iterate through all distance records
			if(dist < dnormal[j]){  //if this node is closer than the nearest neighbor j

				for (int k = m-1; k > j; k--){  //shift everything from j to m one position back

					dnormal[k] = dnormal[k-1];
					near[k] = near[k-1];

				}

				dnormal[j] = dist;  //put the data for this new nearest neighbor in the spot previously occupied by neighbor j
				near[j] = getNode(i);

				break;  //stop looking for others

			}	

		}

	}

	return near;

}

/**
 * method to calculate the potential of all the nodes in the network, where the potential of any node pair is defined as alpha * r and the potential between the clouse of positive charge and any node is beta/2*r^2
*/
float PPGrowingNetwork1D::calculatePotential(){

	float potential = 0;

	#pragma omp parallel shared(potential)
	{

		float localSum = 0;  //sum of potential energies local to this thread
		SpatialVertex* a;
		SpatialVertex* b;

		#pragma omp for schedule(guided)

		for(int i = 0; i < N; i++){  //for each node

		
			a = getNode(i);

			for(int j = i+1; j < N; j++){  //for each pair of nodes

				b = getNode(j);

				if(a != b){  //if they are not the same
					localSum += -alpha*abs(DISTANCE_1D(a ,b));  //add the potential energy for this pair

				}

			}

			localSum += (beta/2)*(X(a) * X(a));  //add the potential energy due to the cloud of positive charge

		}

		#pragma omp atomic

		potential += localSum;  //add the local sums together

	}

	return potential;

}

/**
 * method to call whichever equalization algorith is in use at the time, currently gradientDescent()
 */
void PPGrowingNetwork1D::equalize(){

	gradientDescent(gamma/(N*N), tolerance*N, maxItr);

}

void PPGrowingNetwork1D::gradientDescent(float gamma, float tolerance, long int maxItr){

	float* netForce[N];	//local array to store forces on each node
	float previousPotential = DBL_MAX;  //record of the last potential for minimization

	while(abs(previousPotential - calculatePotential()) > tolerance && maxItr > 0){

		
		previousPotential = calculatePotential();

		#pragma omp parallel shared(netForce)
		{

			#pragma omp for schedule(guided)
			for(int i = 0; i < N; i++){  //for each node

				netForce[i] = sumForces(getNode(i));

			}

		}

		#pragma omp parallel shared(netForce)
		{

			#pragma omp for schedule(guided)
			for(int i = 0; i < N; i++){  //for every node

					getNode(i)->position[0] += gamma * netForce[i][0];  //displace the node by gamma * netForce

				delete[] netForce[i]; //clear netForce to avoid memory leak

			}

		}
		
		maxItr--;

	}

}

PPGrowingNetwork1D::~PPGrowingNetwork1D(){

}



