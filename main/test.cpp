#include "../ngraph/nball.h"
#include "../growingnetwork3d/growingnetwork3d.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <sstream>
#include <sys/time.h>
#include <cfloat>

using namespace std;

double calculateVariance(double* data, long int len);
double* findNormalizedNearestNeighborDistances(NBall* net);

int main(){
	
	long int n = 1000, m = 3;
	GrowingNetwork3D* net = new GrowingNetwork3D(n, m);
	
	for(long int i = 0; i < n; i++){
		
		for(long int j = i; j < n; j++){
			
			cout<<"("<<i<<"->"<<j<<"): "<<net->distance(net->getNode(i), net->getNode(j))<<endl;
			
		}
		
	}
	
}

double calculateVariance(double* data, long int len){
	
	double mean = 0;
	double var = 0;
	
	for(long int i = 0; i < len; i++){
		
		mean += data[i];
		
	}
	
	mean /= len;
	
	for(long int i = 0; i < len; i++){
		
		var += pow(mean - data[i], 2.0);
		
	}
	
	var /= (len-1);
	
	return var;
	
	
	
}

double* findNormalizedNearestNeighborDistances(NBall* net){
	
	double* distances = new double[net->N];
	
	#pragma omp parallel shared(distances)
	{
		
		#pragma omp parallel for schedule(guided)
		for(long int i = 0; i < net->N; i++){
			
			double min = DBL_MAX;
			
			for(long int j = 0; j < net->N; j++){
				
				if(net->getNode(i) == net->getNode(j)){
					
					continue;
					
				}
				
				if(net->linearDistance(net->getNode(i), net->getNode(j)) < min){
					
					min = net->linearDistance(net->getNode(i), net->getNode(j));
					
				}
				
			}
			
			distances[i] = min;
					
		}
			
	}
	
	for(long int i = 0; i < net->N; i++){
		
		distances[i] /= ((net->radius) * pow(net->N, -1.0/net->DIM));	
		
	}
	
	return distances;
	
}
