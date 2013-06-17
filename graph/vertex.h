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
		
	Vertex();
	Vertex(long int time);
	
	void addNeighbor(Vertex* neighbor);
	double clusteringCoefficient();
	long int getStartTime();
	bool hasNeighbor(Vertex* neighbor);
	Vertex* getNeighbor(long int i);
	
protected:

	long int startTime;
	
};

#endif
