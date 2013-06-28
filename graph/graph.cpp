#include "vertex.h"
#include "graph.h"
#include <vector>
#include <climits>

using namespace std;

/**
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
	
	#pragma omp parallel shared(sum, num)
	{
		
		#pragma omp for schedule(guided)
		for(long int i = 0; i < N; i++){	// for each node
		
			Graph* dup = new Graph();	// create a duplicate graph
			
			for(long int j = 0; j < N; j++){	// with N duplicate nodes
				
				dup->addNode(new Vertex());
				
			}
			
			for(long int j = 0; j < N; j++){	// for each duplicate node
				
				for(long int k = 0; k < K(j); k++){	// for each neighbor
					
					long int index = this->indexOf(this->getNode(j)->getNeighbor(k));	// find the index of that neighbor
					dup->getNode(j)->addNeighbor(dup->getNode(index));	// link the appropriate duplicates
					
				}
				
			}
		
			memoize(dup->getNode(i));	// memoize the duplicate graph, starting from node i
			
			for(long int j = 0; j < N; j++){	// for each duplicate node
			
				if(j == i){	// if the duplicate node is the same as the starting node
					
					continue;	// skip it; the shortest path of zero is irrelevant
					
				}
				
				#pragma omp atomic
				sum += dup->getNode(j)->distanceFromInitial;	// add the length of the shortest path
								
				#pragma omp atomic
				num++;	// note that one additional node has been counted
				
			}
			
			delete dup;	// destroy the duplicate
			
		}
		
	}
	
	return sum/num;	// return the sum of the shortest path lengths divided by the number of such paths
	
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
	
		if(getNode(i)->neighbors.size() == k){
		
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
		
		long int index = getNode(i)->neighbors.size();
		
		dist[getNode(i)->neighbors.size()]++;
		
	}
	
	for(int i = 0; i < N; i++){
	
		ddist[i] = dist[i]/(double)N;
		
	}
	
	return ddist;
	
}

double Graph::averageDegree(){
	
	long int sum = 0;
	
	for(int i = 0; i < N; i++){
		
			sum += getNode(i)->neighbors.size();
		
	}
	
	return ((double)sum)/N;
	
}

double Graph::averageClusteringCoefficient(){
	
	double sum = 0;
	
	#pragma omp parallel shared(sum)
	{
	
		#pragma omp for schedule(guided)
		for(int i = 0; i < N; i++){
		
			double coef = getNode(i)->clusteringCoefficient();
			
			#pragma omp atomic
			sum += coef;
			
		}
		
	}
	
	return sum/N;
	
}

/*
 * method to find the index of a given node in the vector<Vertex*> nodes
 * returns -1 if the given node is not in the data structure
 */
long int Graph::indexOf(Vertex* node){
	
	for(int i = 0; i < N; i++){
	
		if(node == getNode(i)){
		
			return i;
			
		}
		
	}
	
	return -1;
	
}

Vertex* Graph::getNode(long int i){
	
	return nodes.at(i);
	
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
		
			memoize(root->getNeighbor(i), path, distance + 1);
		
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
		
			clean(root->getNeighbor(i));
			
		}
		
	}
	
}

Graph::~Graph(){
	
	for(long int i = 0; i < N; i++){
		
		delete getNode(i);
		
	}
	
}
