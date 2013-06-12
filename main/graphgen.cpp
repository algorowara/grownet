#include "../growingnetwork2d/growingnetwork2d.h"
#include "../growingnetwork3d/growingnetwork3d.h"
#include <iostream>
#include <cfloat>
#include <omp.h>

int main(){	// potential as a function of N
	
	long int num_reps = 1000;
	GrowingNetwork3D* net;
	
	for(long int n = 2; n < 100; n++){
		
		double max = DBL_MIN;
		double avg = 0;
		double min = DBL_MAX;
		
		for(long int i = 0; i < num_reps; i++){
			
			net = new GrowingNetwork3D(n, 0);
			double potential = net->calculatePotential();
			
			if(potential < min){
			
				min = potential;
				
			}
			
			if(potential > max){
				
				max = potential;
				
			}
			
			
			
			avg += potential;
			delete net;
		
		}
		
		
		avg /= (num_reps * sqrt(n));
		cout<<n<<" "<<min<<endl;
		
	}
	
}
