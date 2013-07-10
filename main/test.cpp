#include "../ngraph/ngraph.h"
#include "../ngraph/nball.h"
#include "../ngraph/nsphere.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <sstream>
#include <sys/time.h>
#include <cfloat>
#include <algorithm>
#include <iterator>

using namespace std;

float calculateVariance(float* data, long int len);
float* findNormalizedNearestNeighborDistances(NBall* net);

int main(){
	
	long int n = 1000, m = 3, d = 2;
	NBall* net = new NBall(n, m, d);
	
}

float calculateVariance(float* data, long int len){
	
	float mean = 0;
	float var = 0;
	
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

float* findNormalizedNearestNeighborDistances(NBall* net){
	
	float* distances = new float[net->N];
	
	#pragma omp parallel shared(distances)
	{
		
		#pragma omp parallel for schedule(guided)
		for(long int i = 0; i < net->N; i++){
			
			float min = DBL_MAX;
			
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
