#include "../growingnetwork3d/growingnetwork3d.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(){
	
	double tol_min = 0.0, tol_max = 5.0, tol_step = 0.05;
	double gam_min = 0.0, gam_max = 100, gam_step = 1.0;
	long int sample = 4;
	long int n = 100, m = 3;
	
	for(double tol = tol_min + tol_step; tol <= tol_max; tol += tol_step){
		
		for(double gam = gam_min + gam_step; gam <= gam_max; gam += gam_step){
			
			double avg = 0.0;
			
			#pragma omp parallel shared(avg)
			{
				
				#pragma omp for schedule(dynamic)
				for(long int i = 0; i < sample; i++){
					
					GrowingNetwork3D* net = new GrowingNetwork3D(n, m, gam, tol);
					
					#pragma omp critical
					{
						
						avg += ((net->calculatePotential() - net->calculateMinimumPotential())/net->calculateMinimumPotential())/sample;
						
					}
					
				}
				
			}
			
			cout<<tol<<" "<<gam<<" "<<avg<<endl;
			
		}
		
	}
	
}
