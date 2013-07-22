#ifndef VERTEX_H
#define VERTEX_H

#include <vector>
#include <list>

using namespace std;

class Vertex{

public:
	
	long int distanceFromInitial;
	vector<Vertex*> pathFromInitial;
	vector<Vertex*> neighbors;
	long int startTime;
		
	Vertex();
	Vertex(long int time);
	
	void addNeighbor(Vertex* neighbor);
	float clusteringCoefficient();
	long int getStartTime();
	bool hasNeighbor(Vertex* neighbor);
	Vertex* getNeighbor(long int i);
	
};

#endif
