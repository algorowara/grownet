#include "../ngraph/nball.h"
#include "../growingnetwork2d/growingnetwork2d.h"
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
	
	long int nend = 10000, m = 3, d = 2;
	long int s = 16, nstep = 500;
	double data[nend/nstep][s];
	double mean[nend/nstep];
	double variance[nend/nstep];
	
	for(long int i = 0; i < s; i++){
		
		NBall* net = NULL;
		
		for(long int n = 0; n < nend; n += nstep){
		
			if(net == NULL){
				
				net = new NBall(nstep, m, d);
				net->equalizationPeriod = 10;
				
			}
			
			else{
				
				net->grow(nstep);
				
			}
			
			data[n/nstep][i] = net->averagePathLength();
		
		}
		
	}
	
	for(long int n = 0; n < nend; n += nstep){
		
		for(long int j = 0; j < s; j++){
			
			mean[n/nstep] += data[n/nstep][j]/s;
			
		}
		
		variance[n/nstep] = calculateVariance(data[n/nstep], s);
		
	}
	
	for(long int n = 0; n < nend; n += nstep){
		
		cout<<(n + nstep)<<" "<<mean[n/nstep]<<" "<<variance[n/nstep]<<endl;
		
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
