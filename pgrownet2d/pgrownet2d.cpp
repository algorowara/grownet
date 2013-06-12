#include "../pgrownet2d/pgrownet2d.h"
#include <vector>
#include <iostream>
#include <ctime>
#include <cmath>
#include <stdlib.h>

using namespace std;

PositiveChargeGrowingNetwork2D::PositiveChargeGrowingNetwork2D(long int n, long int m){

	srand(std::time(NULL));

	radius = 1;

	for(long int i = 0; i < m+1; i++){	//for the first m+1 nodes, which form a clique

		SpatialVertex* newNode = new SpatialVertex(DIM, randomLocation(), getTime()); //generate a new node in DIM dimensions

		nodes.push_back(newNode);

		for(long int j = 0; j < nodes.size()-1; j++){	//link it to all previously created nodes

			newNode->addNeighbor(nodes.at(j));

		}

		
		equalize();	//equalize the distribution of nodes
		tick();
		n--;

	}

	grow(n);	//grow the remaining nodes normally

}

void PositiveChargeGrowingNetwork2D::grow(long int n){

	while(n > 0){

		SpatialVertex* newNode = new SpatialVertex(DIM, randomLocation(), getTime());
		SpatialVertex** neighbors = findMNearestNeighbors(newNode);

		for(int i = 0; i < m; i++){

			newNode->addNeighbor(neighbors[i]);

		}

		nodes.push_back(newNode);
		equalize();

		delete neighbors;
		tick();
		n--;

	}

}

/**
 * randomly select a point within the circle using theta = 2 * pi * * u and p_radius = radius * v, where u and v are taken uniformly   * from the interval (0, 1)
 */
double* PositiveChargeGrowingNetwork2D::randomLocation(){

	double* position = new double[DIM];

	double theta = 2 * M_PI * (((double)rand())/RAND_MAX);
	double p_radius = radius * (((double)rand())/RAND_MAX); 

	position[0] = p_radius * cos(theta);
	position[1] = p_radius * sin(theta);

	return position;

}

/**
 * method to calculate the linear distance between two nodes if the * nodes are not of a type to have a position in space return 0, as * they are unknown
 */ 
double PositiveChargeGrowingNetwork2D::linearDistance(Vertex* a, Vertex* b){

	if(sizeof(a) == sizeof(SpatialVertex) && sizeof(b) == sizeof(SpatialVertex)){

		SpatialVertex* a_loc = (SpatialVertex*)a;
		SpatialVertex* b_loc = (SpatialVertex*)b;

		return DISTANCE(a_loc, b_loc);

	}

	return 0;

}

/**
 * method to calculate the repulzive force between two nodes based  * on their displacement, here it falls off as 1/r
 * dynamically allocates an array which should be deleted by the
 * caller after use
 */
double* PositiveChargeGrowingNetwork2D::sumForces(SpatialVertex* node){

	double* force = new double[DIM];	//allocate a new array for the force vector
	SpatialVertex* other;	//local placeholder for any other node in two-body interactions
	double magnitude, distance, pmagnitude, pdistance;	//local placeholders for the magnitude of a force (electron-electron and electron-positive) and the distance between two nodes

	for(int i = 0; i < DIM; i++){	//ensure all components are set to zero first

		force[i] = 0;

	}

	for(int i = 0; i < N; i++){	//for every node in the graph

		if(nodes.at(i) == node){	//except this one

			continue;	//skip this one

		}

		other = nodes.at(i);
		magnitude = alpha / DISTANCE(node, other);	//calculate the unitless magnitude of the repulsive electron-electron force
		distance = DISTANCE(node, other);

		force[0] += magnitude * (X(node) - X(other))/distance; 
		force[1] += magnitude * (Y(node) - Y(other))/distance;  //multiply the magnitude by the unit vector of the distancement

	}
	pdistance = sqrt((X(node) * X(node)) + (Y(node) * Y(node)));	//this is the distance from the origin to the node (i.e. radius of node)
	pmagnitude = - beta / pdistance;	//calcualte the unitless magnitude of the attractive electron-positive charge force
	
	force[0] += pmagnitude * (X(node) / pdistance);
	force[1] += pmagnitude * (Y(node) / pdistance);	//multiple the magnitude by a vector pointin radially outwards

	return force;

}

/**
 * method to find the m nearest neighbors to some starting node
 */
SpatialVertex** PositiveChargeGrowingNetwork2D::findMNearestNeighbors(SpatialVertex* start){

	SpatialVertex** near = new SpatialVertex*[m];
	double dsquare[m];	//local record of the distance-squared of the m nearest neighbors

	for(int i = 0; i < m; i++){

		dsquare[i] = DBL_MAX; //everything is intitially far away

	}

	for(int i = 0; i < N; i++){

		if(nodes.at(i) == start){	//do not consider the distance between the starting node and itself

			continue;

		}

		double square = DISTANCE_SQUARED(start, nodes.at(i));

		for(int j = 0; j < m; j++){	//iterate through all distance squared records

			if(square < dsquare[j]){	//if this node is closer than the nearest neighbor j

				for int(k = m-1; k > j; k++){	//shift everything from j to m one position back

					dsquare[k] = dsquare[k-1];
					near[k] = near[k-1];

				}	
			
 				dsquare[j] = square;	//put the data for this new nearest neighbor in the spot previously occupied by neighbor j
				near[j] = nodes.at(i);

			}

		}

	}

	return near;

}

/** 
 * method to calculate the potential of all the nodes in the
 * network, where the potential of any node pair is defined as 
 * log(r-r') and the potential between the cloud of positive charge * and any node is log(r) where r is the position of the node and 
 * r' is the position of any other node
 */
double PositiveChargeGrowingNetwork2D::calculatePotential(){

	double potential = 0;

	#pragma omp parallel shared(potential)
	{

		double localSum = 0;	//sum of potential energies local to this thread
		SpatialVertex* a;
		Spatialvertex* b;

		#pragma omp for schedule(guided)

		for(int i = 0; i < N; i++){	//for every node

			for(int j = 0; j < N; j++){	//for every pair of nodes

				a = nodes.at(i);
				b = nodes.at(j);

				if(a != b){	//if they are not the same

					localSum += -alpha*log(DISTANCE(a, b); //add their potential energy to the local sum

				}

			}

			localSum += beta*log(sqrt((X(i) * X(i)) + (Y(i) * Y(i)))); //add the potential energy due to the cloud of positive charge

		}	

		#pragma omp atomic
		potential += localSum;	// add the local sums together

	}

	return potential;

}

/**
 * method to call whichever equalization algorithm is in use at the * time; currently called gradientDescent()
 */
void Growingnetwork3D::equalize(){

	gradientDescent();

}
