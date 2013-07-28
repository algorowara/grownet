#ifndef NGRAPH_H
#define NGRAPH_H

#include "../graph/growingnetwork.h"
#include "../growingnetwork3d/spatialvertex.h"

#define NGRAPH_ID 0
#define NBALL_ID 1
#define NSPHERE_ID 2

/**
 * macro to quickly retrieve a node and cast it to a SpatialVertex*
 * should be used only in the object itself; not suitable for being called by other objects or functions
 */
#define GET_NODE(i) ((SpatialVertex*)(nodes.at(i)))

/**
 * macro to quickly calculate the result of base^exponent
 * where exponent is a non-negative integer (can be zero)
 */
#define POSITIVE_INTEGER_POWER(base, exponent, output){ \
  output = 1; \
  for(long int itr_var = 0; itr_var < exponent; itr_var++){ \
    output *= base; \
  } \
}

using namespace std;

class NGraph : public GrowingNetwork {
	
	public:

		const int DIM;	// the spatial dimension of the graph (not of any embedding space)
		float radius;	// the radius of the graph (if applicable); has not been for radius != 1 on the NBall
		float baseGam;	// the base (non-scaled) timestep of any equalization algorithm
		float baseTol;	// the base (non-scaled) tolerance for inaccuracy of any equalization algorithm
		long int baseItr;	// the base (non-scaled) maximum number of iterations for any equalization algorithm
		long int equalizationThreshold;	// the number of nodes below which GradientDescent is called once per node added
		long int equalizationPeriod;	// the number of nodes added per call to GradientDescent
		float iterationWeights;	// count of weighted iterations over this object's lifetime
		
		NGraph(long int d);
		SpatialVertex* getNode(const long int i) const;
		virtual float* randomLocation() = 0;
		virtual float linearDistance(Vertex* a, Vertex* b) = 0;
		virtual float linearDistance(float* a, float* b) = 0;
	
};

#endif
