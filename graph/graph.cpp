#include "vertex.h"
#include "graph.h"
#include <vector>
#include <climits>

using namespace std;

/*
 * method to add a node to the graph at the end of the current set of nodes
 */
void Graph::addNode(Vertex* node){

	nodes.push_back(node);
	
}

/**
 * method to calculate the mean length of the shortest path between any two given nodes in the graph
 */
double Graph::averagePathLength(){

	long double sum = 0;
	long int num = 0;
	
	for(int i = 0; i < N; i++){
	
		memoize(nodes.at(i));
		
		for(int j = 0; j < N; j++){
		
			sum += nodes.at(j)->distanceFromInitial;
			num++;
			
		}
		
		clean(nodes.at(i));
		
	}
	
	return sum/num;
	
}

/**
 * method to insert a node at a specific position in the graph's set of nodes
 */
void Graph::insertNode(Vertex* node, long int position){
	
	vector<Vertex*>::iterator itr = nodes.begin();
	itr += position;
	
	nodes.insert(itr, node);

}

/**
 * method to return a count of the number of nodes with exactly k edges
 */
long int Graph::nodesWithDegree(long int k){
	
	long int num = 0;
	
	for(int i = 0; i < N; i++){
	
		if(nodes.at(i)->neighbors.size() == k){
		
			num++;
			
		}
		
	}	
	
	return num;
	
}

/*
 * method to return the proportion of nodes of a given degree, for all possible degrees
 * where the resulting array is N long, as N-1 is the maximum possible degree
 */
double* Graph::degreeDistribution(){
	
	long int dist[N];	// this array is local and temporary to minimize the effects of floating-point rounding
									// where the index is the degree k
									// and the content is the number of nodes with degree k
									
	double* ddist = new double[N];
	
	for(int i = 0; i < N; i++){
	
		dist[i] = 0;
		
	}

	for(int i = 0; i < N; i++){
		
		long int index = nodes.at(i)->neighbors.size();
		
		dist[nodes.at(i)->neighbors.size()]++;
		
	}
	
	for(int i = 0; i < N; i++){
	
		ddist[i] = dist[i]/(double)N;
		
	}
	
	return ddist;
	
}

double Graph::averageDegree(){
	
	long int sum = 0;
	
	for(int i = 0; i < N; i++){
		
			sum += nodes.at(i)->neighbors.size();
		
	}
	
	return ((double)sum)/N;
	
}

double Graph::averageClusteringCoefficient(){
	
	double sum = 0;
	
	for(int i = 0; i < N; i++){
	
		sum += nodes.at(i)->clusteringCoefficient();
		
	}
	
	return sum/N;
	
}

/*
 * method to find the index of a given node in the vector<Vertex*> nodes
 * returns -1 if the given node is not in the data structure
 */
long int Graph::indexOf(Vertex* node){
	
	for(int i = 0; i < N; i++){
	
		if(node == nodes.at(i)){
		
			return i;
			
		}
		
	}
	
	return -1;
	
}

void Graph::memoize(Vertex* root){

	vector<Vertex*>* path = new vector<Vertex*>;
	root->pathFromInitial = *path;
	memoize(root, *path, 0);
	
}

/*
 * memoize all vertices with distance notations
 */
void Graph::memoize(Vertex* root, vector<Vertex*> path, long int distance){
	
	if(root->distanceFromInitial > distance){
		
		root->distanceFromInitial = distance;
		path.push_back(root);
		root->pathFromInitial = path;
		
		for(long int i = 0; i < root->neighbors.size(); i++){
		
			memoize(root->neighbors.at(i), path, distance + 1);
		
		}
		
	}
		
}

/*
 * remove all distance notations from this connected graph
 * and all records of paths
 */
void Graph::clean(Vertex* root){
	
	if(root->distanceFromInitial != LONG_MAX){
	
		root->distanceFromInitial = LONG_MAX;
		
		while(root->pathFromInitial.size() > 0){
		
			root->pathFromInitial.pop_back();
			
		}
		
		for(long int i = 0; i < root->neighbors.size(); i++){
		
			clean(root->neighbors.at(i));
			
		}
		
	}
	
}