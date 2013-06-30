#ifndef GRAPH_H
#define GRAPH_H

#include "vertex.h"
#include <vector>

using namespace std;

#define N nodes.size()
#define K(i) getNode(i)->neighbors.size()

class Graph{

public:

	vector<Vertex*> nodes;
	
	void addNode(Vertex* node);
	double averagePathLength();
	long int nodesWithDegree(long int k);
	double* degreeDistribution();
	double averageDegree();
	double averageClusteringCoefficient();
	double* nodeBetweenness();
	void insertNode(Vertex* node, long int position);
	long int indexOf(Vertex* node);
	Vertex* getNode(long int i);	
	static void memoize(Vertex* root);
	static void memoize(Vertex* root, vector<Vertex*> path, long int distance);
	static void clean(Vertex* root);
	~Graph();
	
};

#endif
