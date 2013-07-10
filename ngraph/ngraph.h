#ifndef NGRAPH_H
#define NGRAPH_H

#include "../graph/growingnetwork.h"
#include "../growingnetwork3d/spatialvertex.h"

#define NGRAPH_ID 0
#define NBALL_ID 1
#define NSPHERE_ID 2

using namespace std;

class NGraph : public GrowingNetwork {
	
	public:

		const long int DIM;
		double radius;
		double baseGam;
		double baseTol;
		double baseItr;
		long int equalizationThreshold;	// the number of nodes below which GradientDescent is called once per node added
		long int equalizationPeriod;	// the number of nodes added per call to GradientDescent
		double iterationWeights;	// count of weighted iterations over this object's lifetime
		
		NGraph(long int d);
		SpatialVertex* getNode(const long int i) const;
		virtual double* randomLocation() = 0;
		virtual double linearDistance(Vertex* a, Vertex* b) = 0;
		virtual double linearDistance(double* a, double* b) = 0;
	
};

#endif
