#include "ngraph.h"
#include "../growingnetwork3d/spatialvertex.h"
#include <fstream>

NGraph::NGraph(long int d) : DIM(d){
	
	
	
}

SpatialVertex* NGraph::getNode(const long int i) const{
	
	return (SpatialVertex*)nodes.at(i);
	
}

void exportObject(const NGraph* g, const char* filename){
	
	ofstream outfile(filename, ios::out | ios::trunc);	// open the output file
	
	// output all of the parameters of the model
	outfile<<"# dimension ="<<g->DIM<<endl;
	outfile<<"# radius = "<<g->radius<<endl;
	outfile<<"# base gamma = "<<g->baseGam<<endl;
	outfile<<"# base tolerance = "<<g->baseTol<<endl;
	outfile<<"# base iterations = "<<g->baseItr<<endl;
	outfile<<"# equalization threshold = "<<g->equalizationThreshold<<endl;
	outfile<<"# equalization period = "<<g->equalizationPeriod<<endl;
	outfile<<"# iteration weights = "<<g->iterationWeights<<endl;
	outfile<<"# time = "<<g->getTime()<<endl;
	outfile<<"# m = "<<g->m<<endl;
	
	for(long int i = 0; i < g->N; i++){	// for each node in the graph
		
		outfile<<i<<" ";	// output the index of the node, now its only identifier
		
		outfile<<g->getNode(i)->dimension<<" ";	// output the dimension of the node so that is position can be accurately read
		
		for(long int j = 0; j < g->getNode(i)->dimension; j++){	// for each dimension of the space
			
			outfile<<g->getNode(i)->position[j]<<" ";	// output the position
			
		}
		
		for(long int j = 0; j < g->K(i); j++){	// for each neighbor of this node
			
			long int index = g->indexOf(g->getNode(i)->getNeighbor(j));	// identify it by index
			
			outfile<<index<<" ";
			
		}
		
		outfile<<endl;	// signify the end of this node's information with a new line
		
	}
	
	outfile.close();	// close the output file
	
}
