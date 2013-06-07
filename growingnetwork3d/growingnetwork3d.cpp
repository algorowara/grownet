#include "spatialvertex.h"
#include "growingnetwork3d.h"
#include <cmath>
#include <cstdlib>

using namespace std;

GrowingNetwork3D::GrowingNetwork3D(long int n, long int m){
	
	for(int i = 0; i < m+1; i++){	// for the first m+1 nodes, which form a clique
		
		SpatialVertex* newNode = new SpatialVertex(DIM, randomLocation(), getTime());	// generate a new node in DIM dimensions
		
		for(int j = 0; j < nodes.size()-1; j++){	// link it to all previously created nodes
			
			newNode->addNeighbor(nodes.at(j));
			
		}
		
		equalize();	// equalize the distribution of nodes
		tick();
		n--;
		
	}
	
	grow(n);	// grow the remaining nodes normally
	
}

void GrowingNetwork3D::grow(long int n){
	
	while(n > 0){
		
		SpatialVertex* newNode = new SpatialVertex(DIM, randomLocation(), getTime());
		vector<SpatialVertex*>* neighbors = findMNearestNeighbors(newNode);
		
		for(int i = 0; i < m; i++){
			
			newNode->addNeighbor(neighbors->at(i));
			
		}
		
		addNode(newNode);
		equalize();
		
		delete neighbors;
		tick();		
		n--;
		
	}
	
}

/**
 * randomly select a point on the surface of a sphere
 * using theta = 2 * pi * u, phi = acos(2 * v - 1), where
 * u and v are taken uniformly from the interval (0, 1)
 */
double* GrowingNetwork3D::randomLocation(){
	
	double* position = new double[DIM];
	
	position[0] = 2 * M_PI * (((double)rand())/RAND_MAX);
	position[1] = acos(2 * (((double)rand())/RAND_MAX) - 1);
	
	return position;
	
}

/**
 * calculates the linear displacement in the theta- and phi-directions
 * between two vertices (from a to b)
 * assumes a unit sphere, so a radius of one is implied
 * dynamically allocates an array which should be deleted after use
 */
double* GrowingNetwork3D::displacement(Vertex* a, Vertex* b){
	
	double* displacement = new double[DIM]
	
	double delsig = 
	
	return displacement;
	
}
