#include "spatialvertex.h"
#include "growingnetwork3d.h"

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
