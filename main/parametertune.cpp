#include "../growingnetwork3d/growingnetwork3d.h"
#include <omp.h>
#include <cfloat>

using namespace std;

double utility(GrowingNetwork3D* net);

int main(){
	
	long int sample_size = 25;
	long int n = 100, m = 3;
	double base_gamma_step = 0.01, base_tol_step = 0.01;
	long int base_itr_step = 1;
	double gamma = 2.8, gamma_step = base_gamma_step, tolerance = 0.6, tol_step = base_tol_step;
	long int iterations = 40, itr_step = base_itr_step;
	long int gen = 0, genmax = 400;
	GrowingNetwork3D* net = new GrowingNetwork3D(n, m, gamma, tolerance, iterations);

	double prev_util = utility(net);
	double best_gamma, best_tolerance;
	long int best_iterations;

	while(gen < genmax){
		
		best_gamma = gamma;
		best_tolerance = tolerance;
		best_iterations = iterations;
		
		#pragma omp parallel shared(best_gamma, best_tolerance, best_iterations, prev_util) private(net)
		{
		
			#pragma omp for schedule(dynamic)
			for(long int i = 0; i < 8; i++){
				
				double test_gamma = abs(gamma + ((i & 1<<0) ? -1:1) * gamma_step);
				double test_tolerance = abs(tolerance + ((i & 1<<1) ? -1:1) * tol_step);
				long int test_iterations = abs(iterations + ((i & 1<<2) ? -1:1) * itr_step);
				double test_util = 0;
				
				for(long int j = 0; j < sample_size; j++){
					
					net = new GrowingNetwork3D(n, m, test_gamma, test_tolerance, test_iterations);
					test_util += utility(net)/sample_size;
					delete net;
					
				}
				
				#pragma omp critical (comparison)
				{
					
					if(test_util > prev_util){
						
						prev_util = test_util;
						best_gamma = test_gamma;
						best_tolerance = test_tolerance;
						best_iterations = test_iterations;
						
					}
				
				}
				
			}
			
		}
		
		if(best_gamma == gamma && best_tolerance == tolerance && best_iterations == iterations){
			
			gamma_step += base_gamma_step;
			tol_step += base_tol_step;
			itr_step += base_itr_step;
			
		}
		
		else{
			
			gamma = best_gamma;
			tolerance = best_tolerance;
			iterations = best_iterations;
			
			if(gamma_step > base_gamma_step && tol_step > base_tol_step && itr_step > base_itr_step){
				
				gamma_step = base_gamma_step;
				tol_step = base_tol_step;
				itr_step = base_itr_step;
				
			}
			
		}
		
		gen++;
		//cout<<"End generation "<<gen<<": "<<gamma<<", "<<tolerance<<", "<<iterations<<endl;
		
	}
	
	cout<<"Found after "<<gen<<" generations:"<<endl;
	cout<<"\tGamma = "<<gamma<<endl;
	cout<<"\tTolerance = "<<tolerance<<endl;
	cout<<"\tIterations = "<<iterations<<endl;
	cout<<"\twith "<<prev_util<<" utility"<<endl;
	
}

double utility(GrowingNetwork3D* net){
	
	return net->calculateMinimumPotential(net->N, 3) - net->calculatePotential();
	
}
