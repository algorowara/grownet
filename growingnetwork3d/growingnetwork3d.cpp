#include "spatialvertex.h"
#include "growingnetwork3d.h"
#include <cmath>
#include <cstdlib>
#include <cfloat>

using namespace std;

GrowingNetwork3D::GrowingNetwork3D(long int n, long int m){
	
	radius = 1;
	
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
		SpatialVertex** neighbors = findMNearestNeighbors(newNode);
		
		for(int i = 0; i < m; i++){
			
			newNode->addNeighbor(neighbors[i]);
			
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
 * see doc.odt for specification of the coordinate system
 */
double* GrowingNetwork3D::randomLocation(){
	
	double* position = new double[DIM];
	
	double theta = 2 * M_PI * (((double)rand())/RAND_MAX);
	double phi = acos(2 * (((double)rand())/RAND_MAX - 1);
	
	position[0] = radius * sin(phi) * cos(theta);
	position[1] = radius * sin(phi) * sin(theta);
	position[2] = radius * cos(phi);
	
	return position;
	
}

/**
 * method to calculate the repulsive force between two nodes based on their displacement
 * currently falls off as 1/r^2
 * dynamically allocates an array which should be deleted by the caller after use
 */
double* GrowingNetwork3D::calculateForce(double* disp){
	
	double* force = new double[DIM];
	double magnitude = 1.0 / DIST_SQUARED(disp);	// calculate the magnitude as 1/r^2, units are irrelevant
	
	for(long int i = 0; i < DIM; i++){

		force[i] = magnitude * disp[i] / sqrt(DIST_SQUARED(disp));	// multiply the magnitude
																	// by the unit vector
																	
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
	
	
	for(int i = 0; i < N; i++){
		
		if(nodes.at(i) == start){	// do not consider the distance between the starting node and itself
		
			continue;
			
		}
	
		double* disp = calculateAngularDisplacement(start, nodes.at(i));
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

/**
 * method to calculate the length of the straight line between two nodes
 * note that the result is linear with radius;
 * if the radius is dependent on N, results for different N
 * will not be directly comparable
 */
double GrowingNetwork3D::linearDistance(SpatialVertex* a, SpatialVertex* b){
	
	return ((X(b) - X(a)) * (X(b) - X(a)) + (Y(b) - Y(a)) * (Y(b) - Y(a)) + (Z(b) - Z(a)) * (Z(b) - Z(a)));
	
}
