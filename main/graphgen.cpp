#include "../growingnetwork2d/growingnetwork2d.h"
#include <iostream>
#include <omp.h>

int main(){	// average degree as a function of node age
	
	GrowingNetwork2D* net;
	long int n = 1000000, m = 2;
	double average[n];
	int num_reps = 100;
	
	for(int i = 0; i < n; i++){
		
		average[i] = 0;
		
	}
	
	#pragma omp parallel shared(average) private(net, i)
	{
		
		#pragma omp for schedule(guided)
	
		for(int i = 0; i < num_reps; i++){
			
			cout<<"Sample "<<i<<" being produced by thread "<<omp_get_thread_num()<<endl;
			net = new GrowingNetwork2D(n, m);
			
			for(int j = 0; j < n; j++){
			
				#pragma omp atomic
				average[n - net->nodes.at(j)->getStartTime()] += net->K(j);
				
			}
			
			delete net;
			
		}
		
	}
	
	for(int i = 0; i < n; i+= 1000){
		
		average[i] /= num_reps;
		cout<<i<<" "<<average[i]<<endl;
		
	}
	
}
