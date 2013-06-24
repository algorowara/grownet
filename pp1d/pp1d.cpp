#include "../pp1d/pp1d.h"
#include <vector>
#include <iostream>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <assert.h>

using namespace std;

PPGrowingNetwork1D::PPGrowingNetwork1D(long int n, long int m, double gamma, double tolerance, long int maxItr){

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
	alpha = .01;
	if(N == 0){
		beta = alpha;
	}
	else {
		beta = alpha*N;
	}

	for(long int i = 0; i < m+1; i++){ //for the first m+1 nodes, which form a clique

		SpatialVertex* newNode = new SpatialVertex(DIM, randomLocation(), getTime()); //generate a new node in DIM dimensions (1 here)

		addNode(newNode);

		for(long int j = 0; j < nodes.size()-1; j++){ //link to the others in this clique

			newNode->addNeighbor(getNode(j));

		}

		equalize(); //minimize the potential energy

		tick();
		n--;

	}

	grow(n) //grow the remaining nodes normall

}

SpatialVertex* PPGrowingNetwork1D::getNode(long int i){

		return (SpatialVertex*)nodes.at(i);

}

void PPGrowingNetwork1D::grow(long int n){

	while(n > 0){

		SpatialVertex* newNode = new SpatialVertex(DIM, randomLocation(); getTime());
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
double* PPGrowingNetwork1D::randomLocation(){

	double* position = new double[DIM];
	
	double p_radius = radius*radius * (((double)rand())/RAND_MAX);
	long int p_sign = (rand() % 2) - 1;

	position[0] = sqrt(p_radius) * p_sign;

	return position;

}

/**
 * calculate the repulsive force between two nodes based on their displacement. It dynamically allocated an array which should be deleted by caller after use
*/
double* PPGrowingNetwork1D::sumForces(SpatialVertex* node){

	double* force = new double[DIM]; //allocate an array for the force vector
	SpatialVertex* other; //local placeholder for any other node in 2 body interactions
	double magnitude, distance, pmagnitude, pdistance //local placeholders for magnitude of force and location of nodes for electron - electron interactions and electron cloud interactions

	for(int i = 0; i < DIM; i++){ //set the forces initially to 0

		force[i] = 0;

	}

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
	double dnormal[m]; //local record of the distance of the m nearest neighbors

	for(int i = 0; i < m; i++){

		dnormal[i] = DBL_MAX; //everything is initially far away

	}

	for(int i = 0; i < N; i++){

		if(getNode(i) == start){  //do not consider the starting node...

			continue;

		}

		double dist = DISTANCE_1D(start, getNode(i));

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
double PPGrowingNetwork1D::calculatePotential(){

	double potential = 0;

	#pragma omp parallel shared(potential)
	{

		double localSum = 0;  //sum of potential energies local to this thread
		SpatialVertex* a;
		SpatialVertex* b;

		#pragma omp for schedule(guided)

		for(int i = 0; i < N; i++){  //for each node

			for(int j = i+1; j < N; j++){  //for each pair of nodes

				a = getNode(i);
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

	gradientDescent(gamma, tolerance, maxItr);

}

void PPGrowingNetwork1D::gradientDescent(double gamma, double tolerance, long int maxItr){

	double* netForce[N];	//local array to store forces on each node
	double previousPotential = DBL_MAX;  //record of the last potential for minimization

	while(abs(previousPotential - calculatePotential()) > tolerance && maxItr > 0){

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

				for(int j = 0; j < DIM; j++){  //for every dimension (only one here)

					getNode(i)->position[j] += gamma * netForce[i][j];  //displace the node by gamma * netForce

				}

				delete[] netForce[i]; //clear netForce to avoid memory leak

			}

		}

		previousPotential = calculatePotential();

		maxItr--;

	}

}

PPGrowingNetwork1D::~PPGrowingNetwork1D(){

}



