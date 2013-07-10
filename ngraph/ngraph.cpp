#include "ngraph.h"
#include "../growingnetwork3d/spatialvertex.h"
#include <fstream>
#include <cmath>

NGraph::NGraph(long int d) : DIM(d){
	
	
	
}

SpatialVertex* NGraph::getNode(const long int i) const{
	
	return (SpatialVertex*)nodes.at(i);
	
}
