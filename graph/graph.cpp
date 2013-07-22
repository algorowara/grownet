#include "vertex.h"
#include "graph.h"
#include <vector>
#include <climits>
#include <cstring>
#include <cmath>
#include <list>
#include <queue>
#include <stack>
#include <algorithm>
#include <iostream>

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
float Graph::averagePathLength(){

	long int sum = 0;	// the sum of all shortest path distances
	long int num = 0;	// the number of shortest paths recorded
	
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
		
			dup->memoize(dup->getNode(i));	// memoize the duplicate graph, starting from node i
			
			for(long int j = 0; j < N; j++){	// for each duplicate node
			
				if(j == i){	// if the duplicate node is the same as the starting node
					
					continue;	// skip it; the shortest path of zero is irrelevant
					
				}
				
				#pragma omp critical (summing)
				{
					
					sum += dup->getNode(j)->distanceFromInitial;	// add the length of the shortest path
					num++;	// note that one additional node has been counted
					
				}
				
			}
			
			delete dup;	// destroy the duplicate
			
		}
		
	}
	
	return ((float)sum)/num;	// return the sum of the shortest path lengths divided by the number of such paths
	
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
float* Graph::nodeBetweenness(){

	float* betweennessCentrality = new float[N];	//a dynamically allocated array giving the betweennesscentrality of each node
	memset(betweennessCentrality, 0, N*sizeof(float));	//start with 0 everywhere
	
	#pragma omp parallel shared(betweennessCentrality)
	{

	#pragma omp for schedule(guided)
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
		myQ->push(getNode(i));	//put node s in the Q						
		
		while(!myQ->empty()){

			Vertex* vertex = myQ->front();	//retrieve the first member in the Q (this is v in pseudocode)
			myQ->pop();	//and remove it from the Q
			
			stackS->push(vertex);	//add vertex to S

			for(long int k = 0; k < vertex->neighbors.size(); k++){	//for the neighbors of vertex

				Vertex* neighbor = vertex->getNeighbor(k);	//this is w in the pseudocode

				if(dist[indexOf(neighbor)] < 0){	//if w was found for the first time

					myQ->push(neighbor);	//add w to the Q

					dist[indexOf(neighbor)] = dist[indexOf(vertex)] + 1;	//d[w] = d[v] + 1

				}

				if(dist[indexOf(neighbor)] == (dist[indexOf(vertex)] + 1)){	//shortest path to w via v?

					sigma[indexOf(neighbor)] += sigma[indexOf(vertex)];	//add this node to the current tally
	
					vertices[indexOf(neighbor)].push_back(vertex);	//add v to the list P[w]					
				}

			}

		}

		long int* delta = new long int[N];	//initialize the delta array
		memset(delta, 0, N*sizeof(long int));	//make sure it's empty

		while(!stackS->empty()){	//while S is not empty

			Vertex* top = stackS->top();	//this is w again in pseudocode

			stackS->pop();	//remove w from the stack

			for(long int l = 0; l < vertices[indexOf(top)].size(); l++){	//for each element in the list P[w]	
				Vertex* vert = vertices[indexOf(top)].back();	//this is v here
				vertices[indexOf(top)].pop_back();	//remove this element (v) from the list				

				delta[indexOf(vert)] += (sigma[indexOf(vert)]/sigma[indexOf(top)]) * (1+delta[indexOf(top)]);	//delta[v] = delta[v] + sigma[v]/sigma[w] * (1+delta[w])
		
			}

			if(top != getNode(i)){	//if w != s

				#pragma omp atomic
				betweennessCentrality[indexOf(top)] += delta[indexOf(top)];	//Cb[w] = Cb[w] + delta[w]

			}

		}

		//clear memory to avoid memleak
		delete stackS;
		delete[] vertices; 
		delete[] sigma;
		delete[] dist;
		delete myQ;
		delete[] delta;
		
	}

	}	
	//now we just need to normalize these results
	for(long int q = 0; q < N; q++){

		betweennessCentrality[q] = betweennessCentrality[q] / ((N-1)*(N-2)/2);

	}

	return betweennessCentrality;

}

/**
 * method to return a count of the number of nodes with exactly k edges
 */
long int Graph::nodesWithDegree(long int k) const{
	
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
float* Graph::degreeDistribution() const{
	
	long int dist[N];	// this array is local and temporary to minimize the effects of floating-point rounding
						// where the index is the degree k
						// and the content is the number of nodes with degree k
									
	float* ddist = new float[N];	// dynamically allocated array to store the degree distribution
									// normalized over N
	memset(dist, (long int)0, N * sizeof(long int));

	for(int i = 0; i < N; i++){	// for each node in the graph
		
		long int index = getNode(i)->neighbors.size();	// record the number of neighbors it has
		
		dist[getNode(i)->neighbors.size()]++;	// note that there is one more node with that number of neighbors
		
	}
	
	for(int i = 0; i < N; i++){	// for each possible degree (maximum degree is N-1)
	
		ddist[i] = dist[i]/(float)N;	// the proportion of nodes with that degree
										// is the number of nodes with that degree divided by N
		
	}
	
	return ddist;
	
}

/**
 * method to find the mean degree of a node in the graph
 */
float Graph::averageDegree() const{
	
	long int sum = 0;	// total number of "degrees"
	
	for(int i = 0; i < N; i++){	// for each node
		
			sum += K(i);	// add the degree of that node to the sum
		
	}
	
	return ((float)sum)/N;	// average over all nodes
	
}

/**
 * multithreaded method to calculate the transitivity ratio of the graph
 * obtained by calculating a weighted average of the local clustering coefficient
 * where each node is weighted by k(k-1)/2, the number of possible triangles
 */
float Graph::transitivity() const{
	
	float sum = 0;	// the weighted sum of local clustering coefficients, equal to the number of triangles
	long int weightsum = 0;	// the sum of the weights, equal to the number of possible triangles
	
	#pragma omp parallel shared(sum, weightsum)
	{
	
		#pragma omp for schedule(guided)
		for(int i = 0; i < N; i++){	// for each node
		
			float coef = getNode(i)->clusteringCoefficient();	// calculate the clustering coefficient
			long int weight = (K(i) * (K(i)-1))/2;	// and the weight
			
			#pragma omp critical (addition)
			{
				
				sum += coef * weight;	// add the weighted contribution (number of triangles) to the sum
				weightsum += weight;	// and the weight (number of possible triangles) to the sum)
				
			}
			
		}
		
	}
	
	return sum/weightsum;	// return the weighted average; the weighted sum over the weight
	
}

/**
 * multithreaded method to calculate the mean clustering coefficient of a graph
 * weighing the contributions of each node equally
 */
float Graph::clustering() const{
	
	float sum = 0;	// sum of all local clustering coefficients
	
	#pragma omp parallel shared(sum)
	{
	
		#pragma omp for schedule(guided)
		for(int i = 0; i < N; i++){	// for each node
		
			float coef = getNode(i)->clusteringCoefficient();	// calculate the clustering coefficient
			
			#pragma omp atomic
			sum += coef;	// and add it to the sum
			
		}
		
	}
	
	return sum/N;	// return the simple mean
	
}

/**
 * method to find the index of a given node in the vector<Vertex*> nodes
 * returns -1 if the given node is not in the data structure
 */
long int Graph::indexOf(Vertex* node) const{
	
	for(int i = 0; i < N; i++){	// for each node
	
		if(node == getNode(i)){	// if this is the target node
		
			return i;	// return its index
			
		}
		
	}
	
	return -1;	// if the target node is not found, return a clearly invalid index
	
}

/**
 * method to retrieve a node by index without requiring external access to the Graph::nodes
 * will be overridden in derivde classes to cast nodes to the appropriate derived class of Vertex
 */
Vertex* Graph::getNode(long int i) const{
	
	return nodes.at(i);
	
}

/**
 * method to use a breadth-first queue-based search to memoize all nodes
 * with the minimum distance from some root node and the corresponding shortest path
 */
void Graph::memoize(Vertex* root){

	root->distanceFromInitial = 0;	// the root is always a distance of zero from itself
	root->pathFromInitial.push_back(root);	// the root is always the first node in any shortest path from the root
	
	vector<Vertex*> set = this->nodes;	// create a temporary container for all the nodes
	
	while(set.size() > 0){	// while there are elements in the set
	
		long int index = 0;	// the index of the node which is closest to the root (possibly the root itself)
		long int minimumDistance = LONG_MAX;	// the shortest known distance from the root
		
		for(long int i = 0; i < set.size(); i++){	// for all nodes
			
			if(set.at(i)->distanceFromInitial < minimumDistance){	// find the one with the shortest distance
				
				index = i;	// and remember it
				minimumDistance = set.at(i)->distanceFromInitial;	// and remember the shortest distance for future comparisons
				
			}
			
		}

		Vertex* node = set.at(index);	// take this element which is closest to the root
		set.erase(set.begin() + index);	// and remove it from the set
		
		if(node->distanceFromInitial == LONG_MAX){	// if, somehow, the least distance is "infinite"
			
			break;	// stop running, as somehow an unconnected component has been encountered
			
		}
		
		for(long int i = 0; i < node->neighbors.size(); i++){	// for each neighbor of this closest node
			
			long int alt = node->distanceFromInitial + 1;	// there is a path from the node which is one additional edge longer
			
			if(alt < node->getNeighbor(i)->distanceFromInitial){	// if this path is shorter than the one previously found
				
				node->getNeighbor(i)->distanceFromInitial = alt;	// use it 
				node->getNeighbor(i)->pathFromInitial = node->pathFromInitial;	// and remember the path
				node->getNeighbor(i)->pathFromInitial.push_back(node->getNeighbor(i));	// from the previous node to this one
				
			}
			
		}
		
	}
	
	return;	// return statement placed here for clarity only
	
}

/**
 * remove all distance notations from this connected graph
 * and all records of paths
 */
void Graph::clean(){
	
	for(long int i = 0; i < N; i++){	// for each node
		
		getNode(i)->distanceFromInitial = LONG_MAX;	// set the distanceFromInitial to the value which signals "not yet memoized"
		
		while(getNode(i)->pathFromInitial.size() > 0){	// remove all elements of the pathFromInitial vector
			
			getNode(i)->pathFromInitial.pop_back();
			
		}
		
	}
	
	return;	// return statement placed here for clarity only
	
}

/**
 * method to compare two nodes to check which one has a greater distance from some root node
 */
bool Graph::compareDistancesFromInitial(const Vertex* a, const Vertex* b){
	
	return a->distanceFromInitial < b->distanceFromInitial;
	
}

/**
 * destructor which also deletes all nodes
 */
Graph::~Graph(){
	
	for(long int i = 0; i < N; i++){	// for each node
		
		delete getNode(i);	// call its destructor
		
	}
	
}
