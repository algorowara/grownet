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
	double weightedClusteringCoefficient();
	double unweightedClusteringCoefficient();
	double* nodeBetweenness();
	void insertNode(Vertex* node, long int position);
	long int indexOf(Vertex* node) const;
	Vertex* getNode(long int i) const;	
	void memoize(Vertex* root);
	void clean();
	static bool compareDistancesFromInitial(const Vertex* a, const Vertex* b);
	~Graph();
	
};

#endif
