#include "vertex.h"
#include "graph.h"
#include <vector>
#include <climits>
#include <cstring>
#include <list>
#include <queue>
#include <stack>

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
 * method to find the betweenness centrality of a given node. Returns the unitless value of betweenness
 * which is defined as the sum over node pairs of the ratio of shortest paths passing through this node
 * to all such shortest paths between pairs
 */
double* Graph::nodeBetweenness(){

	double* betweennessCentrality = new double[N];	//a dynamically allocated array giving the betweennesscentrality of each node
	memset(betweennessCentrality, 0, N*sizeof(double));	//start with 0 everywhere
	

	for(long int i = 0; i < N; i++){	//for all of the nodes (we're finding betweeneness centrality of node i), it's s in the pseudocode
	
		stack<Vertex*>* stackS = new stack<Vertex*>;	//in the pseudocode this is S

		list<Vertex*>* vertices = new list<Vertex*>[N];	//in the pseudocode this is P[w], and it actually needs to be an array of lists indexed by w

		long int* sigma = new long int[N];	//in the pseudocode this is sigma
		long int* dist = new long int[N];	//in the pseudocode this is d
	
		for(long int j = 0; j < N; j++){	//for the other nodes (we're comparing to node j), it's t in the pseudocode

			if (j == i){	//if this is the node we're on, set it to one

				sigma[j] = 1;
				dist[j] = 0;
	
			}
			else {	//otherwise sigma = 0

				sigma[j] = 0;
				dist[j] = -1;

			}	

		}

		queue<Vertex*>* myQ = new queue<Vertex*>;	//this is Q in the pseudocode
		(*myQ).push(getNode(i));	//put node s in the Q						
		
		while(!(*myQ).empty()){

			Vertex* vertex = (*myQ).front();	//retrieve the first member in the Q (this is v in pseudocode)
			(*myQ).pop();	//and remove it from the Q
			
			(*stackS).push(vertex);	//add vertex to S

			for(long int k = 0; k < vertex->neighbors.size(); k++){	//for the neighbors of vertex

				Vertex* neighbor = vertex->getNeighbor(k);	//this is w in the pseudocode

				if(dist[indexOf(neighbor)] < 0){	//if w was found for the first time

					(*myQ).push(neighbor);	//add w to the Q

					dist[indexOf(neighbor)] = dist[indexOf(vertex)] + 1;	//d[w] = d[v] + 1

				}

				if(dist[indexOf(neighbor)] == dist[indexOf(vertex)] + 1){	//shortest path to w via v?

					sigma[indexOf(neighbor)] += sigma[indexOf(vertex)];	//add this node to the current tally
	
					vertices[indexOf(neighbor)].push_back(vertex);	//add v to the list P[w]					
				}

			}

		}

		long int* delta = new long int[N];	//initialize the delta array

		while(!(*stackS).empty()){	//while S is not empty

			Vertex* top = (*stackS).top();	//this is w again in pseudocode

			(*stackS).pop();	//remove w from the stack

			for(long int l = 0; l < vertices[indexOf(top)].size(); l++){	//for each element in the list P[w]	
				Vertex* vert = vertices[indexOf(top)].back();	//this is v here
				vertices[indexOf(top)].pop_back();	//remove this element (v) from the list				

				delta[indexOf(vert)] += (sigma[indexOf(vert)]/sigma[indexOf(top)]) * (1+delta[indexOf(top)]);	//delta[v] = delta[v] + sigma[v]/sigma[w] * (1+delta[w])
		
			}

			if(top != getNode(i)){	//if w != s

				betweennessCentrality[indexOf(top)] += delta[indexOf(top)];	//Cb[w] = Cb[w] + delta[w]

			}

		}

	}	

	return betweennessCentrality;

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
