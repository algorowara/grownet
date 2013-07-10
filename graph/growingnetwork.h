#ifndef GROWINGNETWORK_H
#define GROWINGNETWORK_H

#include "graph.h"

class GrowingNetwork : public Graph {

public:

	long int m;

	long int getTime() const;	
	virtual void grow(long int n) = 0;
	long int edgeAge(Vertex* a, Vertex* b);
	double* edgeAgeVsBetweenness();
	virtual double linearDistance(Vertex* a, Vertex* b) = 0;
	double* edgeAgeVsLinearDistance();
	double* nodeAgeVsDegree();
	
protected:

	long int time;
	
	void tick();

};

#endif
