#ifndef NGRAPH_H
#define NGRAPH_H

#include "../graph/growingnetwork.h"
#include "../growingnetwork3d/spatialvertex.h"

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
		static void exportObject(const NGraph* g, const char* filename);
		static NGraph* importObject(const char* filename);
		virtual void equalize() = 0;
	
};

#endif
