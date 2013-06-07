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
	
	position[0] = 2 * M_PI * (((double)rand())/RAND_MAX);
	position[1] = acos(2 * (((double)rand())/RAND_MAX) - 1);
	
	return position;
	
}

/**
 * calculates the linear displacement in the theta- and phi-directions
 * between two vertices (from a to b)
 * dynamically allocates an array which should be deleted after use
 */
double* GrowingNetwork3D::calculateDisplacement(SpatialVertex* a, SpatialVertex* b){
	
	//TODO: implement this
	return NULL;
	
}

/**
 * method to calculate the repulsive force between two nodes based on their displacement
 * currently falls off as 1/r^2
 */
double* GrowingNetwork3D::calculateForce(double* disp){
	
	double* force = new double[DIM];
	double magnitude = 1.0 / DIST_SQUARED(disp);	// calculate the magnitude as 1/r^2, units are irrelevant
	
	force[0] = magnitude * disp[0] / sqrt(DIST_SQUARED(disp));	// multiply the magnitude
	force[1] = magnitude * disp[1] / sqrt(DIST_SQUARED(disp));	// by the unit vector
	
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
	
		double* disp = calculateDisplacement(start, nodes.at(i));
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
 * method to calculate the linear distance covered by an edge
 * for some edge that is a straight line
 * if no edge exists between the specified nodes, returns 0
 * note that the result is linear with radius;
 * if the radius is dependent on N, results for different N
 * will not be directly comparable
 */
double GrowingNetwork3D::edgeLinearDistance(SpatialVertex* a, SpatialVertex* b){
	
	if(!a->hasNeighbor(b)){	// if the specified pair of nodes does is not linked
	
		return 0;	// no linear distance is covered
		
	}
	
	// convert the spherical coordinates of both nodes to Cartesian coordinates
	// and store the Cartesian displacement
	double aCart[3];
	double bCart[3];
	double disCart[3];
	
	aCart[0] = radius * sin(PHI(a)) * cos(THETA(a));
	aCart[1] = radius * sin(PHI(a)) * sin(THETA(a));
	aCart[2] = radius * cos(PHI(a));
	
	bCart[0] = radius * sin(PHI(b)) * cos(THETA(b));
	bCart[1] = radius * sin(PHI(b)) * sin(THETA(b));
	bCart[2] = radius * cos(PHI(a));
	
	for(int i = 0; i < 3; i++){
		
		disCart[i] = bCart[i] - aCart[i];
		
	}
	
	// return the square root of the sum of the squares of the cartesian components
	return sqrt(disCart[0] * disCart[0] + disCart[1] * disCart[1]  + disCart[2] * disCart[2]);
	
}
