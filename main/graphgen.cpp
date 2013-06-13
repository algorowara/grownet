#include "../growingnetwork2d/growingnetwork2d.h"
#include "../growingnetwork3d/growingnetwork3d.h"
#include <iostream>
#include <cfloat>
#include <omp.h>

int main(){	// potential as a function of N

	long int n = 100, m = 3;
	GrowingNetwork3D* net;
	double minPot = DBL_MAX;
	double minGam, minTol;
	long int minItr;
	
	for(double gamma = 0.1; gamma <= 2.0; gamma += 0.1){
		
		for(double tolerance = 0.0; tolerance <= 1.0; tolerance += 0.05){
			
			for(long int maxItr = 0; maxItr <= 100; maxItr += 5){
				
				net = new GrowingNetwork3D(n, m, gamma, tolerance, maxItr);
				double pot = net->calculatePotential();
				
				if(pot < minPot){
					
					minPot = pot;
					minGam = gamma;
					minTol = tolerance;
					minItr = maxItr;
					
				}
				
				delete net;
				
			}
			
		}
		
	}
	
	cout<<"Minimum Potential for n = 16: "<<minPot<<endl;
	cout<<"\t at gamma = "<<minGam<<endl;
	cout<<"\t at tolerance = "<<minTol<<endl;
	cout<<"\t with "<<minItr<<" iterations"<<endl;
	
}
