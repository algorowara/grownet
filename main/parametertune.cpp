#include "../growingnetwork3d/growingnetwork3d.h"
#include <omp.h>
#include <cfloat>

using namespace std;

double utility(GrowingNetwork3D* net);

int main(){
	
	GrowingNetwork3D* net;
	long int n = 100, m = 3, num = 16;
	double minGam, minTol;
	long int minItr;
	double avgUtil = 0;
	double maxUtil = DBL_MIN;
	
	#pragma omp parallel shared(minGam, minTol, minItr, maxUtil) private(net, avgUtil)
	{
		
		for(double baseGam = 0.2; baseGam <= 3.0; baseGam += 0.2){
			
			for(double baseTol = 0.0; baseTol <= 2.0; baseTol += 0.1){
				
				#pragma omp for schedule(dynamic, 1)
				for(long int baseItr = 10; baseItr <= 200; baseItr += 10){
					
					avgUtil = 0.0;
					
					for(long int i = 0; i < num; i++){
						
						net = new GrowingNetwork3D(n, m, baseGam, baseTol, baseItr);
						avgUtil += utility(net)/num;
						delete net;
						
					}
					
					#pragma omp critical
					{
						
						if(avgUtil > maxUtil){
							
							minGam = baseGam;
							minTol = baseTol;
							minItr = baseItr;
							maxUtil = avgUtil;
							
						}
						
						cout<<"Completed test of "<<baseGam<<", "<<baseTol<<", "<<baseItr<<endl;
						
					}
					
				}
				
			}
			
		}
	
	}
	
	cout<<"Maximizing Gamma: "<<minGam<<endl;
	cout<<"Maximizing Tolerance: "<<minTol<<endl;
	cout<<"Maximizing Iterations: "<<minItr<<endl;
	cout<<"Maximum Utility: "<<maxUtil<<endl;
	
}

double utility(GrowingNetwork3D* net){
	
	double excess = net->calculatePotential() - net->calculateMinimumPotential();
	long int iterations = net->baseItr;
	
	return -1 * (excess + ((double)iterations)/net->N);
	
}
