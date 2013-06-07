#include "spatialvertex.h"
#include "growingnetwork3d.h"
#include <cmath>
#include <cstdlib>
#include <cfloat>

#define DIST_SQUARED(displacement) (displacement[0] * displacement[0] + displacement[1] * displacement[1])

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
	

	
	return displacement;
	
}

/**
 * method to calculate the repulsive force between two nodes based on their displacement
 * currently falls off as 1/r^2
 */
double* GrowingNetwork3D::force(double* disp){
	
	double* force = new double[DIM];
	double magnitude = 1.0 / DIST_SQUARED(disp);	// calculate the magnitude as 1/r^2, units are irrelevant
	
	force[0] = magnitude * disp[0] / sqrt(DIST_SQUARED(disp));	// multiply the magnitude
	force[1] = magnitude * displacement[1] / sqrt(DIST_SQUARED(disp));	// by the unit vector
	
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
	
	
	for(int i = 0; i < N; i++){
		
		if(nodes.at(i) == start){	// do not consider the distance between the starting node and itself
		
			continue;
			
		}
	
		double* disp = displacement(start, nodes.at(i));
		double square = DIST_SQUARED(disp);
		
		for(int j = 0; j < m; j++){	// iterate through all distance-squared records
			
			if(square < dsquare[j]){	// if this node is closer than the nearest neighbor j
				
				for(int k = m-1; k > j; k++){	// shift everything from j to m one position back
					
					dsquare[k] = dsquare[k-1];
					near[k] = near[k-1];
					
				}
				
				dsquare[j] = square;	// put the data for this new nearest neighbor
				near[j] = nodes.at(i);	// in the spot previously occupied by neighbor j
				
			}
			
		}
		
		delete disp;
		 
	}
	
	return near;
	 
}
